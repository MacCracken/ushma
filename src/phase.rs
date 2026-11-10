//! Phase transitions — Clausius-Clapeyron, latent heat, phase diagram lookup.
//!
//! All SI units: kelvins, pascals, joules, kilograms.

use std::borrow::Cow;

use serde::{Deserialize, Serialize};

use crate::error::{Result, UshmaError};
use crate::state::{ATM, GAS_CONSTANT};
use crate::validate::{finite_ge, finite_gt, invalid_param};

/// Thermodynamic phase of matter.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum Phase {
    Solid,
    Liquid,
    Gas,
    SupercriticalFluid,
}

/// Phase transition reference data for a substance.
///
/// All values at standard conditions unless noted. Sources: NIST Chemistry WebBook.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SubstanceData {
    /// Substance name.
    pub name: Cow<'static, str>,
    /// Triple point temperature (K).
    pub triple_point_t: f64,
    /// Triple point pressure (Pa).
    pub triple_point_p: f64,
    /// Critical point temperature (K).
    pub critical_t: f64,
    /// Critical point pressure (Pa).
    pub critical_p: f64,
    /// Critical density (kg/m³).
    pub critical_density: f64,
    /// Latent heat of fusion at 1 atm (J/kg).
    pub latent_heat_fusion: f64,
    /// Latent heat of vaporization at 1 atm (J/kg).
    pub latent_heat_vaporization: f64,
    /// Normal melting point at 1 atm (K).
    pub melting_point: f64,
    /// Normal boiling point at 1 atm (K).
    pub boiling_point: f64,
    /// Molar mass (kg/mol).
    pub molar_mass: f64,
}

impl SubstanceData {
    /// `true` when the triple-point pressure exceeds 1 atm (e.g. CO₂): at 1 atm
    /// the solid sublimes at [`boiling_point`](Self::boiling_point) and no
    /// liquid exists.
    #[inline]
    #[must_use]
    pub fn sublimes_at_atmospheric(&self) -> bool {
        self.triple_point_p > ATM
    }

    /// Vapor (liquid–gas) pressure at `temperature` (Pa), valid from the triple
    /// point to the critical point.
    ///
    /// Piecewise integrated Clausius-Clapeyron (ln P linear in 1/T) through the
    /// triple point, the normal boiling point (when it lies between them) and
    /// the critical point, so the curve meets all three exactly. Typical
    /// accuracy is a few percent (water: within 5 % of IAPWS-IF97; CO₂: within
    /// 1 %). For water, prefer `steam::saturation_pressure` (exact IF97)
    /// when the `steam` feature is enabled.
    #[inline]
    pub fn vapor_pressure(&self, temperature: f64) -> Result<f64> {
        if !finite_ge(temperature, self.triple_point_t) || temperature > self.critical_t {
            return Err(UshmaError::UndefinedPhaseTransition {
                temperature,
                pressure: f64::NAN,
            });
        }
        let tp = (self.triple_point_t, self.triple_point_p);
        let crit = (self.critical_t, self.critical_p);
        let nbp = (self.boiling_point, ATM);
        let p = if self.boiling_point > self.triple_point_t && self.boiling_point < self.critical_t
        {
            if temperature <= self.boiling_point {
                log_interp_inv_t(temperature, tp, nbp)
            } else {
                log_interp_inv_t(temperature, nbp, crit)
            }
        } else {
            log_interp_inv_t(temperature, tp, crit)
        };
        Ok(p)
    }

    /// Sublimation (solid–gas) pressure below the triple point (Pa).
    ///
    /// For substances that sublime at 1 atm the curve runs through the triple
    /// point and (`boiling_point`, 1 atm). Otherwise it is integrated
    /// Clausius-Clapeyron from the triple point with L_sub ≈ L_fus + L_vap.
    pub fn sublimation_pressure(&self, temperature: f64) -> Result<f64> {
        if !finite_gt(temperature, 0.0) || temperature > self.triple_point_t {
            return Err(UshmaError::UndefinedPhaseTransition {
                temperature,
                pressure: f64::NAN,
            });
        }
        let tp = (self.triple_point_t, self.triple_point_p);
        if self.sublimes_at_atmospheric() && self.boiling_point < self.triple_point_t {
            return Ok(log_interp_inv_t(temperature, tp, (self.boiling_point, ATM)));
        }
        let l_sub_molar =
            (self.latent_heat_fusion + self.latent_heat_vaporization) * self.molar_mass;
        clausius_clapeyron_pressure(
            self.triple_point_p,
            self.triple_point_t,
            temperature,
            l_sub_molar,
        )
    }

    /// Specific latent heat of sublimation (J/kg).
    ///
    /// For substances that sublime at 1 atm this is the Clausius-Clapeyron slope
    /// of the sublimation curve through the triple point and the 1 atm
    /// sublimation point (CO₂: ≈ 594 kJ/kg vs ≈ 571 kJ/kg measured). Otherwise
    /// L_fus + L_vap.
    #[must_use]
    pub fn latent_heat_sublimation(&self) -> f64 {
        let (t1, p1) = (self.boiling_point, ATM);
        let (t2, p2) = (self.triple_point_t, self.triple_point_p);
        let d_inv_t = 1.0 / t1 - 1.0 / t2;
        if self.sublimes_at_atmospheric()
            && t1 < t2
            && finite_gt(d_inv_t, 0.0)
            && self.molar_mass > 0.0
        {
            GAS_CONSTANT * (p2 / p1).ln() / d_inv_t / self.molar_mass
        } else {
            self.latent_heat_fusion + self.latent_heat_vaporization
        }
    }

    /// Determine the thermodynamic phase at a given temperature and pressure.
    ///
    /// Uses simplified phase boundaries:
    /// - Supercritical: T > T_c and P > P_c
    /// - Below the triple-point temperature: sublimation curve (solid above it, gas below)
    /// - Between triple and critical point: [`vapor_pressure`](Self::vapor_pressure)
    ///   curve (liquid above it — including compressed liquid at P > P_c — gas below)
    /// - Melting curve approximated as vertical at the melting point
    /// - T > T_c and P ≤ P_c: gas
    ///
    /// On a liquid–gas or solid–gas boundary, returns the condensed phase.
    #[tracing::instrument(level = "debug", skip(self), fields(substance = %self.name))]
    pub fn phase_at(&self, temperature: f64, pressure: f64) -> Result<Phase> {
        if !finite_gt(temperature, 0.0) {
            return Err(UshmaError::InvalidTemperature {
                kelvin: temperature,
            });
        }
        if !finite_gt(pressure, 0.0) {
            return Err(UshmaError::InvalidPressure { pascals: pressure });
        }

        // Supercritical region
        if temperature > self.critical_t && pressure > self.critical_p {
            return Ok(Phase::SupercriticalFluid);
        }

        // Below the triple point: solid or vapor, split by the sublimation curve
        if temperature < self.triple_point_t {
            let p_sub = self.sublimation_pressure(temperature)?;
            return Ok(if pressure >= p_sub {
                Phase::Solid
            } else {
                Phase::Gas
            });
        }

        // Between triple point and critical point: check vaporization curve
        if temperature <= self.critical_t {
            let p_sat = self.vapor_pressure(temperature)?;
            if pressure >= p_sat {
                // Solid-liquid boundary: vertical approximation at the melting point
                return Ok(if temperature <= self.melting_point {
                    Phase::Solid
                } else {
                    Phase::Liquid
                });
            }
            return Ok(Phase::Gas);
        }

        // T > T_c but P <= P_c: gas-like behavior
        Ok(Phase::Gas)
    }
}

/// Integrated Clausius-Clapeyron between two points: ln P linear in 1/T.
#[inline]
fn log_interp_inv_t(t: f64, (t1, p1): (f64, f64), (t2, p2): (f64, f64)) -> f64 {
    let span = 1.0 / t2 - 1.0 / t1;
    if span.abs() < f64::EPSILON {
        return p1;
    }
    let frac = (1.0 / t - 1.0 / t1) / span;
    p1 * (frac * (p2 / p1).ln()).exp()
}

// --- Substance data (NIST Chemistry WebBook) ---

/// Water (H₂O).
pub const WATER_PHASE: SubstanceData = SubstanceData {
    name: Cow::Borrowed("Water"),
    triple_point_t: 273.16,
    triple_point_p: 611.657,
    critical_t: 647.096,
    critical_p: 22_064_000.0,
    critical_density: 322.0,
    latent_heat_fusion: 334_000.0,
    latent_heat_vaporization: 2_260_000.0,
    melting_point: 273.15,
    boiling_point: 373.15,
    molar_mass: 0.018_015,
};

/// Nitrogen (N₂).
pub const NITROGEN: SubstanceData = SubstanceData {
    name: Cow::Borrowed("Nitrogen"),
    triple_point_t: 63.15,
    triple_point_p: 12_520.0,
    critical_t: 126.19,
    critical_p: 3_390_000.0,
    critical_density: 313.3,
    latent_heat_fusion: 25_700.0,
    latent_heat_vaporization: 199_000.0,
    melting_point: 63.15,
    boiling_point: 77.36,
    molar_mass: 0.028_014,
};

/// Oxygen (O₂).
pub const OXYGEN: SubstanceData = SubstanceData {
    name: Cow::Borrowed("Oxygen"),
    triple_point_t: 54.36,
    triple_point_p: 146.3,
    critical_t: 154.58,
    critical_p: 5_043_000.0,
    critical_density: 436.1,
    latent_heat_fusion: 13_900.0,
    latent_heat_vaporization: 213_000.0,
    melting_point: 54.36,
    boiling_point: 90.19,
    molar_mass: 0.031_998,
};

/// Carbon dioxide (CO₂).
pub const CARBON_DIOXIDE: SubstanceData = SubstanceData {
    name: Cow::Borrowed("Carbon dioxide"),
    triple_point_t: 216.55,
    triple_point_p: 518_500.0,
    critical_t: 304.13,
    critical_p: 7_375_000.0,
    critical_density: 467.6,
    latent_heat_fusion: 196_000.0,
    latent_heat_vaporization: 234_000.0,
    // CO₂ sublimes at 1 atm; melting point at 5.18 atm
    melting_point: 216.55,
    boiling_point: 194.65, // sublimation point at 1 atm
    molar_mass: 0.044_009,
};

/// Ammonia (NH₃).
pub const AMMONIA: SubstanceData = SubstanceData {
    name: Cow::Borrowed("Ammonia"),
    triple_point_t: 195.40,
    triple_point_p: 6_060.0,
    critical_t: 405.56,
    critical_p: 11_280_000.0,
    critical_density: 225.0,
    latent_heat_fusion: 332_000.0,
    latent_heat_vaporization: 1_370_000.0,
    melting_point: 195.42,
    boiling_point: 239.82,
    molar_mass: 0.017_031,
};

/// Methane (CH₄).
pub const METHANE: SubstanceData = SubstanceData {
    name: Cow::Borrowed("Methane"),
    triple_point_t: 90.69,
    triple_point_p: 11_696.0,
    critical_t: 190.56,
    critical_p: 4_599_000.0,
    critical_density: 162.7,
    latent_heat_fusion: 58_700.0,
    latent_heat_vaporization: 510_000.0,
    melting_point: 90.69,
    boiling_point: 111.66,
    molar_mass: 0.016_043,
};

/// Ethanol (C₂H₅OH).
pub const ETHANOL: SubstanceData = SubstanceData {
    name: Cow::Borrowed("Ethanol"),
    triple_point_t: 159.0,
    triple_point_p: 0.43,
    critical_t: 513.9,
    critical_p: 6_148_000.0,
    critical_density: 276.0,
    latent_heat_fusion: 108_000.0,
    latent_heat_vaporization: 841_000.0,
    melting_point: 159.05,
    boiling_point: 351.44,
    molar_mass: 0.046_068,
};

/// Mercury (Hg).
pub const MERCURY: SubstanceData = SubstanceData {
    name: Cow::Borrowed("Mercury"),
    triple_point_t: 234.32,
    triple_point_p: 0.000_17,
    critical_t: 1750.0,
    critical_p: 167_000_000.0,
    critical_density: 5500.0,
    latent_heat_fusion: 11_300.0,
    latent_heat_vaporization: 295_000.0,
    melting_point: 234.32,
    boiling_point: 629.88,
    molar_mass: 0.200_592,
};

/// All built-in substances for iteration.
pub const ALL_SUBSTANCES: &[&SubstanceData] = &[
    &WATER_PHASE,
    &NITROGEN,
    &OXYGEN,
    &CARBON_DIOXIDE,
    &AMMONIA,
    &METHANE,
    &ETHANOL,
    &MERCURY,
];

/// Heat required to melt a given mass: Q = m·L_fusion (J).
///
/// - `substance`: substance with latent heat data
/// - `mass`: mass in kg
#[inline]
#[must_use]
pub fn heat_of_fusion(substance: &SubstanceData, mass: f64) -> f64 {
    mass * substance.latent_heat_fusion
}

/// Heat required to vaporize a given mass: Q = m·L_vap (J).
///
/// - `substance`: substance with latent heat data
/// - `mass`: mass in kg
#[inline]
#[must_use]
pub fn heat_of_vaporization(substance: &SubstanceData, mass: f64) -> f64 {
    mass * substance.latent_heat_vaporization
}

/// Total heat to move a substance from t_start to t_end at 1 atm,
/// including sensible heating and latent heat at phase boundaries.
///
/// This is the "ice to steam" calculation: sum of sensible heat segments
/// plus latent heats at each phase boundary crossed. Substances that sublime
/// at 1 atm (see [`SubstanceData::sublimes_at_atmospheric`]) get a single
/// sublimation event at `boiling_point` instead of fusion + vaporization.
///
/// - `mass`: mass in kg
/// - `specific_heat`: average c_p across all phases (J/(kg·K)) — simplified assumption
/// - `t_start`: starting temperature (K)
/// - `t_end`: ending temperature (K)
/// - `substance`: substance with melting/boiling points and latent heats
#[tracing::instrument(level = "debug", skip(substance), fields(substance = %substance.name))]
pub fn heat_for_phase_change(
    mass: f64,
    specific_heat: f64,
    t_start: f64,
    t_end: f64,
    substance: &SubstanceData,
) -> Result<f64> {
    if !finite_gt(mass, 0.0) {
        return Err(invalid_param(format_args!(
            "mass {mass} kg must be positive"
        )));
    }
    if !finite_gt(specific_heat, 0.0) {
        return Err(UshmaError::InvalidSpecificHeat {
            value: specific_heat,
        });
    }
    if !finite_gt(t_start, 0.0) {
        return Err(UshmaError::InvalidTemperature { kelvin: t_start });
    }
    if !finite_gt(t_end, 0.0) {
        return Err(UshmaError::InvalidTemperature { kelvin: t_end });
    }

    let (t_lo, t_hi, sign) = if t_end >= t_start {
        (t_start, t_end, 1.0)
    } else {
        (t_end, t_start, -1.0)
    };

    let mut q = 0.0;

    // Sensible heat for entire range
    q += mass * specific_heat * (t_hi - t_lo);

    if substance.sublimes_at_atmospheric() {
        // No liquid at 1 atm (e.g. CO₂): one solid → gas transition at the
        // sublimation point.
        if t_lo < substance.boiling_point && t_hi >= substance.boiling_point {
            q += mass * substance.latent_heat_sublimation();
        }
    } else {
        // Latent heat of fusion if crossing melting point
        if t_lo < substance.melting_point && t_hi >= substance.melting_point {
            q += mass * substance.latent_heat_fusion;
        }

        // Latent heat of vaporization if crossing boiling point
        if t_lo < substance.boiling_point && t_hi >= substance.boiling_point {
            q += mass * substance.latent_heat_vaporization;
        }
    }

    Ok(sign * q)
}

/// Clausius-Clapeyron slope: dP/dT = L/(T·Δv) (Pa/K).
///
/// Exact form for the slope of a phase boundary in P-T space.
///
/// - `latent_heat`: specific latent heat L (J/kg)
/// - `temperature`: phase transition temperature T (K)
/// - `delta_v`: specific volume change Δv = v_gas - v_liquid (m³/kg)
#[inline]
pub fn clausius_clapeyron_slope(latent_heat: f64, temperature: f64, delta_v: f64) -> Result<f64> {
    if !finite_gt(temperature, 0.0) {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    if !finite_ge(delta_v.abs(), 1e-30) {
        return Err(UshmaError::DivisionByZero {
            context: "Δv cannot be zero for Clausius-Clapeyron".into(),
        });
    }
    if !finite_gt(latent_heat, 0.0) {
        return Err(invalid_param(format_args!(
            "latent heat {latent_heat} J/kg must be positive"
        )));
    }
    Ok(latent_heat / (temperature * delta_v))
}

/// Clausius-Clapeyron integrated form: vapor pressure at a target temperature.
///
/// ln(P₂/P₁) = -(L/R)·(1/T₂ - 1/T₁)
///
/// Assumes ideal gas vapor and negligible liquid volume.
///
/// - `p_ref`: reference pressure P₁ (Pa)
/// - `t_ref`: reference temperature T₁ (K)
/// - `t_target`: target temperature T₂ (K)
/// - `latent_heat_molar`: molar latent heat L (J/mol)
#[inline]
pub fn clausius_clapeyron_pressure(
    p_ref: f64,
    t_ref: f64,
    t_target: f64,
    latent_heat_molar: f64,
) -> Result<f64> {
    if !finite_gt(p_ref, 0.0) {
        return Err(UshmaError::InvalidPressure { pascals: p_ref });
    }
    if !finite_gt(t_ref, 0.0) {
        return Err(UshmaError::InvalidTemperature { kelvin: t_ref });
    }
    if !finite_gt(t_target, 0.0) {
        return Err(UshmaError::InvalidTemperature { kelvin: t_target });
    }
    if !finite_gt(latent_heat_molar, 0.0) {
        return Err(invalid_param(format_args!(
            "molar latent heat {latent_heat_molar} J/mol must be positive"
        )));
    }
    let exponent = -(latent_heat_molar / GAS_CONSTANT) * (1.0 / t_target - 1.0 / t_ref);
    let p = p_ref * exponent.exp();
    if !p.is_finite() {
        return Err(invalid_param(format_args!(
            "Clausius-Clapeyron pressure overflows (exponent {exponent})"
        )));
    }
    Ok(p)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_phase_enum_traits() {
        let p = Phase::Liquid;
        assert_eq!(p, Phase::Liquid);
        assert_ne!(p, Phase::Gas);
        let cloned = p;
        assert_eq!(cloned, Phase::Liquid);
        let debug = format!("{p:?}");
        assert!(debug.contains("Liquid"));
    }

    #[test]
    fn test_all_substances_count() {
        assert_eq!(ALL_SUBSTANCES.len(), 8);
    }

    // --- NIST reference value checks ---

    #[test]
    fn test_water_triple_point() {
        assert!((WATER_PHASE.triple_point_t - 273.16).abs() < 0.01);
        assert!((WATER_PHASE.triple_point_p - 611.657).abs() < 1.0);
    }

    #[test]
    fn test_water_critical_point() {
        assert!((WATER_PHASE.critical_t - 647.096).abs() < 0.01);
        assert!((WATER_PHASE.critical_p - 22_064_000.0).abs() < 1000.0);
    }

    #[test]
    fn test_water_latent_heats() {
        // Fusion: 334 kJ/kg (NIST)
        assert!((WATER_PHASE.latent_heat_fusion - 334_000.0).abs() < 1000.0);
        // Vaporization: 2260 kJ/kg (NIST)
        assert!((WATER_PHASE.latent_heat_vaporization - 2_260_000.0).abs() < 10_000.0);
    }

    #[test]
    fn test_nitrogen_boiling_point() {
        // N₂ boils at 77.36 K (NIST)
        assert!((NITROGEN.boiling_point - 77.36).abs() < 0.1);
    }

    #[test]
    fn test_co2_critical_point() {
        // CO₂ critical: 304.13 K, 7.375 MPa (NIST)
        assert!((CARBON_DIOXIDE.critical_t - 304.13).abs() < 0.1);
        assert!((CARBON_DIOXIDE.critical_p - 7_375_000.0).abs() < 10_000.0);
    }

    #[test]
    fn test_ammonia_boiling_point() {
        assert!((AMMONIA.boiling_point - 239.82).abs() < 0.1);
    }

    #[test]
    fn test_mercury_melting_point() {
        assert!((MERCURY.melting_point - 234.32).abs() < 0.1);
    }

    #[test]
    fn test_substance_serde_roundtrip() {
        let json = serde_json::to_string(&WATER_PHASE).unwrap();
        let back: SubstanceData = serde_json::from_str(&json).unwrap();
        assert_eq!(back.name, "Water");
        assert!((back.critical_t - 647.096).abs() < 0.001);
    }

    #[test]
    fn test_phase_serde_roundtrip() {
        let json = serde_json::to_string(&Phase::SupercriticalFluid).unwrap();
        let back: Phase = serde_json::from_str(&json).unwrap();
        assert_eq!(back, Phase::SupercriticalFluid);
    }

    // --- Latent heat tests ---

    #[test]
    fn test_heat_of_fusion_water() {
        // 1 kg water: 334 kJ
        let q = heat_of_fusion(&WATER_PHASE, 1.0);
        assert!((q - 334_000.0).abs() < 1000.0);
    }

    #[test]
    fn test_heat_of_vaporization_water() {
        // 1 kg water: 2260 kJ
        let q = heat_of_vaporization(&WATER_PHASE, 1.0);
        assert!((q - 2_260_000.0).abs() < 10_000.0);
    }

    #[test]
    fn test_heat_of_fusion_scales_with_mass() {
        let q1 = heat_of_fusion(&WATER_PHASE, 1.0);
        let q2 = heat_of_fusion(&WATER_PHASE, 2.0);
        assert!((q2 - 2.0 * q1).abs() < 1e-10);
    }

    #[test]
    fn test_heat_for_phase_change_ice_to_steam() {
        // 1 kg water from 263 K (-10°C) to 383 K (110°C)
        // Simplified: use average c_p ≈ 2500 J/(kg·K) across all phases
        // Sensible: 2500 * 120 = 300,000 J
        // Fusion: 334,000 J (crosses 273.15 K)
        // Vaporization: 2,260,000 J (crosses 373.15 K)
        // Total ≈ 2,894,000 J
        let q = heat_for_phase_change(1.0, 2500.0, 263.0, 383.0, &WATER_PHASE).unwrap();
        assert!(q > 2_500_000.0); // must include latent heats
        assert!((q - 2_894_000.0).abs() < 10_000.0);
    }

    #[test]
    fn test_heat_for_phase_change_no_transition() {
        // 1 kg water from 280 K to 350 K — no phase boundary crossed
        let q = heat_for_phase_change(1.0, 4186.0, 280.0, 350.0, &WATER_PHASE).unwrap();
        // Pure sensible: 4186 * 70 = 293,020 J
        assert!((q - 293_020.0).abs() < 1.0);
    }

    #[test]
    fn test_heat_for_phase_change_cooling() {
        // Cooling: 383 K → 263 K should be negative of heating
        let q_heat = heat_for_phase_change(1.0, 2500.0, 263.0, 383.0, &WATER_PHASE).unwrap();
        let q_cool = heat_for_phase_change(1.0, 2500.0, 383.0, 263.0, &WATER_PHASE).unwrap();
        assert!((q_heat + q_cool).abs() < 1e-10);
    }

    #[test]
    fn test_heat_for_phase_change_invalid() {
        assert!(heat_for_phase_change(0.0, 4186.0, 280.0, 350.0, &WATER_PHASE).is_err());
        assert!(heat_for_phase_change(1.0, 4186.0, 0.0, 350.0, &WATER_PHASE).is_err());
        assert!(heat_for_phase_change(1.0, 4186.0, 280.0, 0.0, &WATER_PHASE).is_err());
        // Negative or zero specific heat
        assert!(heat_for_phase_change(1.0, 0.0, 280.0, 350.0, &WATER_PHASE).is_err());
        assert!(heat_for_phase_change(1.0, -100.0, 280.0, 350.0, &WATER_PHASE).is_err());
    }

    // --- Clausius-Clapeyron tests ---

    #[test]
    fn test_clausius_clapeyron_slope_water() {
        // Water at 100°C (373.15 K), 1 atm:
        // L_vap = 2260 kJ/kg, Δv ≈ v_g - v_f ≈ 1.673 - 0.001 ≈ 1.672 m³/kg
        // dP/dT = 2_260_000 / (373.15 * 1.672) ≈ 3622 Pa/K
        // Literature: ~3570 Pa/K (the approximation is within 2%)
        let dp_dt = clausius_clapeyron_slope(2_260_000.0, 373.15, 1.672).unwrap();
        assert!((dp_dt - 3570.0).abs() / 3570.0 < 0.02);
    }

    #[test]
    fn test_clausius_clapeyron_slope_invalid() {
        assert!(clausius_clapeyron_slope(2_260_000.0, 0.0, 1.672).is_err());
        assert!(clausius_clapeyron_slope(2_260_000.0, 373.15, 0.0).is_err());
        assert!(clausius_clapeyron_slope(-100.0, 373.15, 1.672).is_err());
    }

    #[test]
    fn test_clausius_clapeyron_pressure_water() {
        // Water: P_ref = 101325 Pa at T_ref = 373.15 K
        // Molar L_vap = 2260 kJ/kg * 0.018015 kg/mol ≈ 40714 J/mol
        // Predict P at 383.15 K (110°C)
        let l_molar = 2_260_000.0 * 0.018_015;
        let p = clausius_clapeyron_pressure(101_325.0, 373.15, 383.15, l_molar).unwrap();
        // At 110°C, P_sat ≈ 143.4 kPa (NIST). CC approximation should be within 5%.
        assert!((p - 143_400.0).abs() / 143_400.0 < 0.05);
    }

    #[test]
    fn test_clausius_clapeyron_pressure_roundtrip() {
        let l_molar = 2_260_000.0 * 0.018_015;
        // Go from 373.15K to 393.15K
        let p2 = clausius_clapeyron_pressure(101_325.0, 373.15, 393.15, l_molar).unwrap();
        // Come back from 393.15K to 373.15K
        let p_back = clausius_clapeyron_pressure(p2, 393.15, 373.15, l_molar).unwrap();
        assert!((p_back - 101_325.0).abs() / 101_325.0 < 1e-10);
    }

    #[test]
    fn test_clausius_clapeyron_pressure_invalid() {
        assert!(clausius_clapeyron_pressure(0.0, 373.15, 383.15, 40_000.0).is_err());
        assert!(clausius_clapeyron_pressure(101_325.0, 0.0, 383.15, 40_000.0).is_err());
        assert!(clausius_clapeyron_pressure(101_325.0, 373.15, 0.0, 40_000.0).is_err());
        assert!(clausius_clapeyron_pressure(101_325.0, 373.15, 383.15, -1.0).is_err());
    }

    // --- Phase lookup tests ---

    #[test]
    fn test_water_liquid_at_room_temp() {
        let phase = WATER_PHASE.phase_at(300.0, 101_325.0).unwrap();
        assert_eq!(phase, Phase::Liquid);
    }

    #[test]
    fn test_water_gas_above_boiling() {
        let phase = WATER_PHASE.phase_at(400.0, 101_325.0).unwrap();
        assert_eq!(phase, Phase::Gas);
    }

    #[test]
    fn test_water_solid_below_freezing() {
        let phase = WATER_PHASE.phase_at(200.0, 101_325.0).unwrap();
        assert_eq!(phase, Phase::Solid);
    }

    #[test]
    fn test_water_supercritical() {
        // T > 647 K, P > 22 MPa
        let phase = WATER_PHASE.phase_at(700.0, 25_000_000.0).unwrap();
        assert_eq!(phase, Phase::SupercriticalFluid);
    }

    #[test]
    fn test_water_high_temp_low_pressure_is_gas() {
        // T > T_c but P < P_c → gas
        let phase = WATER_PHASE.phase_at(700.0, 101_325.0).unwrap();
        assert_eq!(phase, Phase::Gas);
    }

    #[test]
    fn test_co2_supercritical() {
        // CO₂: T_c=304.13K, P_c=7.375 MPa
        let phase = CARBON_DIOXIDE.phase_at(310.0, 8_000_000.0).unwrap();
        assert_eq!(phase, Phase::SupercriticalFluid);
    }

    #[test]
    fn test_nitrogen_gas_at_room_temp() {
        // N₂ boils at 77 K — room temperature is well above
        let phase = NITROGEN.phase_at(300.0, 101_325.0).unwrap();
        assert_eq!(phase, Phase::Gas);
    }

    #[test]
    fn test_nitrogen_liquid_at_cryogenic() {
        // N₂ at 70 K under moderate pressure should be liquid
        let phase = NITROGEN.phase_at(70.0, 101_325.0).unwrap();
        assert_eq!(phase, Phase::Liquid);
    }

    #[test]
    fn test_phase_at_invalid_inputs() {
        assert!(WATER_PHASE.phase_at(0.0, 101_325.0).is_err());
        assert!(WATER_PHASE.phase_at(300.0, 0.0).is_err());
        assert!(WATER_PHASE.phase_at(-10.0, 101_325.0).is_err());
    }

    #[test]
    fn test_water_high_pressure_liquid() {
        // Water at 350 K under very high pressure → liquid
        let phase = WATER_PHASE.phase_at(350.0, 10_000_000.0).unwrap();
        assert_eq!(phase, Phase::Liquid);
    }

    #[test]
    fn test_all_substances_have_valid_data() {
        for substance in ALL_SUBSTANCES {
            assert!(
                substance.triple_point_t > 0.0,
                "{} triple T",
                substance.name
            );
            assert!(
                substance.triple_point_p >= 0.0,
                "{} triple P",
                substance.name
            );
            assert!(
                substance.critical_t > substance.triple_point_t,
                "{} critical > triple",
                substance.name
            );
            assert!(substance.critical_p > 0.0, "{} critical P", substance.name);
            assert!(
                substance.latent_heat_fusion > 0.0,
                "{} L_fusion",
                substance.name
            );
            assert!(
                substance.latent_heat_vaporization > 0.0,
                "{} L_vap",
                substance.name
            );
            assert!(substance.molar_mass > 0.0, "{} molar_mass", substance.name);
        }
    }

    // --- Regression tests from the 1.4.0 audit ---

    #[test]
    fn test_vapor_pressure_meets_anchor_points() {
        let w = &WATER_PHASE;
        assert!((w.vapor_pressure(373.15).unwrap() / crate::state::ATM - 1.0).abs() < 1e-12);
        assert!((w.vapor_pressure(w.critical_t).unwrap() / w.critical_p - 1.0).abs() < 1e-12);
        assert!(
            (w.vapor_pressure(w.triple_point_t).unwrap() / w.triple_point_p - 1.0).abs() < 1e-12
        );
        assert!(w.vapor_pressure(200.0).is_err());
        assert!(w.vapor_pressure(700.0).is_err());
        assert!(w.vapor_pressure(f64::NAN).is_err());
    }

    #[test]
    fn test_water_vapor_pressure_never_exceeds_critical() {
        // The single-reference curve used to give 26 MPa at T_c (> P_c = 22.064 MPa).
        let w = &WATER_PHASE;
        let p = w.vapor_pressure(640.0).unwrap();
        assert!((p / 20.266e6 - 1.0).abs() < 0.02, "{p}"); // IF97: 20.27 MPa
    }

    #[test]
    fn test_co2_vapor_pressure() {
        // NIST: 1.785 MPa at 250 K, 6.713 MPa at 300 K.
        let c = &CARBON_DIOXIDE;
        assert!((c.vapor_pressure(250.0).unwrap() / 1.785e6 - 1.0).abs() < 0.02);
        assert!((c.vapor_pressure(300.0).unwrap() / 6.713e6 - 1.0).abs() < 0.02);
    }

    #[test]
    fn test_compressed_liquid_above_critical_pressure() {
        // T < T_c, P > P_c → compressed liquid (was Gas).
        assert_eq!(WATER_PHASE.phase_at(640.0, 23e6).unwrap(), Phase::Liquid);
    }

    #[test]
    fn test_co2_phase_regions() {
        let c = &CARBON_DIOXIDE;
        assert_eq!(c.phase_at(300.0, 3e6).unwrap(), Phase::Gas); // was Liquid
        assert_eq!(c.phase_at(250.0, 1e6).unwrap(), Phase::Gas); // was Liquid
        assert_eq!(c.phase_at(250.0, 2.5e6).unwrap(), Phase::Liquid);
        // Dry ice sublimes at 194.65 K at 1 atm: 200 K is vapor (was Solid).
        assert_eq!(c.phase_at(200.0, 101_325.0).unwrap(), Phase::Gas);
        assert_eq!(c.phase_at(180.0, 101_325.0).unwrap(), Phase::Solid);
        assert_eq!(c.phase_at(320.0, 8e6).unwrap(), Phase::SupercriticalFluid);
    }

    #[test]
    fn test_co2_sublimation_curve() {
        let c = &CARBON_DIOXIDE;
        assert!((c.sublimation_pressure(194.65).unwrap() / crate::state::ATM - 1.0).abs() < 1e-12);
        // Measured: ≈ 160.6 kPa at 200 K.
        assert!((c.sublimation_pressure(200.0).unwrap() / 160.6e3 - 1.0).abs() < 0.05);
        assert!(c.sublimation_pressure(250.0).is_err());
    }

    #[test]
    fn test_water_below_triple_point_low_pressure_is_vapor() {
        // Ice sublimation pressure at 250 K is ≈ 76 Pa.
        assert_eq!(WATER_PHASE.phase_at(250.0, 10.0).unwrap(), Phase::Gas);
        assert_eq!(
            WATER_PHASE.phase_at(250.0, 101_325.0).unwrap(),
            Phase::Solid
        );
        let p_sub = WATER_PHASE.sublimation_pressure(250.0).unwrap();
        assert!(p_sub > 60.0 && p_sub < 100.0, "{p_sub}");
    }

    #[test]
    fn test_co2_heat_for_phase_change_uses_sublimation() {
        let c = &CARBON_DIOXIDE;
        assert!(c.sublimes_at_atmospheric());
        assert!(!WATER_PHASE.sublimes_at_atmospheric());
        // 200 → 220 K at 1 atm is all vapor: sensible heat only (fusion was added).
        let q = heat_for_phase_change(1.0, 800.0, 200.0, 220.0, c).unwrap();
        assert!((q - 800.0 * 20.0).abs() < 1e-9);
        // 190 → 300 K crosses the sublimation point once.
        let q = heat_for_phase_change(1.0, 800.0, 190.0, 300.0, c).unwrap();
        let l_sub = c.latent_heat_sublimation();
        assert!((l_sub - 571_000.0).abs() < 40_000.0, "{l_sub}");
        assert!((q - (800.0 * 110.0 + l_sub)).abs() < 1e-6);
    }

    #[test]
    fn test_clausius_clapeyron_pressure_overflow_is_error() {
        assert!(clausius_clapeyron_pressure(1e5, 1e-10, 300.0, 4e4).is_err());
    }

    #[test]
    fn test_phase_functions_reject_nan() {
        assert!(WATER_PHASE.phase_at(f64::NAN, 1e5).is_err());
        assert!(WATER_PHASE.phase_at(300.0, f64::NAN).is_err());
        assert!(clausius_clapeyron_pressure(f64::NAN, 373.15, 380.0, 4e4).is_err());
        assert!(clausius_clapeyron_slope(f64::NAN, 373.15, 1.67).is_err());
        assert!(heat_for_phase_change(1.0, 4186.0, f64::NAN, 350.0, &WATER_PHASE).is_err());
    }
}
