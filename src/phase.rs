//! Phase transitions — Clausius-Clapeyron, latent heat, phase diagram lookup.
//!
//! All SI units: kelvins, pascals, joules, kilograms.

use std::borrow::Cow;

use serde::{Deserialize, Serialize};

use crate::error::{Result, UshmaError};
use crate::state::GAS_CONSTANT;

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

/// Clausius-Clapeyron slope: dP/dT = L/(T·Δv) (Pa/K).
///
/// Exact form for the slope of a phase boundary in P-T space.
///
/// - `latent_heat`: specific latent heat L (J/kg)
/// - `temperature`: phase transition temperature T (K)
/// - `delta_v`: specific volume change Δv = v_gas - v_liquid (m³/kg)
pub fn clausius_clapeyron_slope(latent_heat: f64, temperature: f64, delta_v: f64) -> Result<f64> {
    if temperature <= 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    if delta_v.abs() < 1e-30 {
        return Err(UshmaError::DivisionByZero {
            context: "Δv cannot be zero for Clausius-Clapeyron".into(),
        });
    }
    if latent_heat <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("latent heat {latent_heat} J/kg must be positive"),
        });
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
pub fn clausius_clapeyron_pressure(
    p_ref: f64,
    t_ref: f64,
    t_target: f64,
    latent_heat_molar: f64,
) -> Result<f64> {
    if p_ref <= 0.0 {
        return Err(UshmaError::InvalidPressure { pascals: p_ref });
    }
    if t_ref <= 0.0 {
        return Err(UshmaError::InvalidTemperature { kelvin: t_ref });
    }
    if t_target <= 0.0 {
        return Err(UshmaError::InvalidTemperature { kelvin: t_target });
    }
    if latent_heat_molar <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("molar latent heat {latent_heat_molar} J/mol must be positive"),
        });
    }
    let exponent = -(latent_heat_molar / GAS_CONSTANT) * (1.0 / t_target - 1.0 / t_ref);
    Ok(p_ref * exponent.exp())
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
        }
    }
}
