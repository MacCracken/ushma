//! Equations of state — ideal gas, van der Waals, phase diagrams.
//!
//! All SI units: pascals, cubic meters, kelvins, moles.

use std::borrow::Cow;

use serde::{Deserialize, Serialize};

use crate::error::{Result, UshmaError};

/// Universal gas constant R (J/(mol⋅K)).
pub const GAS_CONSTANT: f64 = 8.314_462_618;

/// Standard atmospheric pressure (Pa).
pub const ATM: f64 = 101_325.0;

/// Standard temperature (K).
pub const STANDARD_TEMP: f64 = 273.15;

/// Ideal gas law: PV = nRT.
///
/// Returns pressure (Pa).
pub fn ideal_gas_pressure(moles: f64, temperature: f64, volume: f64) -> Result<f64> {
    if temperature < 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    if volume <= 0.0 {
        return Err(UshmaError::InvalidVolume {
            cubic_meters: volume,
        });
    }
    Ok(moles * GAS_CONSTANT * temperature / volume)
}

/// Ideal gas volume: V = nRT/P (m³).
pub fn ideal_gas_volume(moles: f64, temperature: f64, pressure: f64) -> Result<f64> {
    if temperature < 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    if pressure <= 0.0 {
        return Err(UshmaError::InvalidPressure { pascals: pressure });
    }
    Ok(moles * GAS_CONSTANT * temperature / pressure)
}

/// Ideal gas temperature: T = PV/(nR) (K).
pub fn ideal_gas_temperature(pressure: f64, volume: f64, moles: f64) -> Result<f64> {
    if moles.abs() < 1e-30 {
        return Err(UshmaError::DivisionByZero {
            context: "moles cannot be zero".into(),
        });
    }
    Ok(pressure * volume / (moles * GAS_CONSTANT))
}

/// Van der Waals equation: [P + a(n/V)²][V - nb] = nRT.
///
/// Returns pressure (Pa).
/// - `a`: attraction parameter (Pa⋅m⁶/mol²)
/// - `b`: volume exclusion parameter (m³/mol)
pub fn van_der_waals_pressure(
    moles: f64,
    temperature: f64,
    volume: f64,
    a: f64,
    b: f64,
) -> Result<f64> {
    if temperature < 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    let v_eff = volume - moles * b;
    if v_eff <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("effective volume {v_eff} m³ is non-positive (V - nb)"),
        });
    }
    let density = moles / volume;
    Ok(moles * GAS_CONSTANT * temperature / v_eff - a * density * density)
}

/// Work done by an ideal gas in isothermal expansion.
///
/// W = nRT⋅ln(V₂/V₁) (joules)
pub fn isothermal_work(moles: f64, temperature: f64, v1: f64, v2: f64) -> Result<f64> {
    if temperature < 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    if v1 <= 0.0 {
        return Err(UshmaError::InvalidVolume { cubic_meters: v1 });
    }
    if v2 <= 0.0 {
        return Err(UshmaError::InvalidVolume { cubic_meters: v2 });
    }
    Ok(moles * GAS_CONSTANT * temperature * (v2 / v1).ln())
}

/// Work done in isobaric (constant pressure) process.
///
/// W = P(V₂ - V₁) (joules)
#[inline]
#[must_use]
pub fn isobaric_work(pressure: f64, v1: f64, v2: f64) -> f64 {
    pressure * (v2 - v1)
}

/// Adiabatic process: PV^γ = constant.
///
/// Returns final temperature for adiabatic expansion/compression.
/// T₂ = T₁(V₁/V₂)^(γ-1)
pub fn adiabatic_temperature(t1: f64, v1: f64, v2: f64, gamma: f64) -> Result<f64> {
    if gamma <= 1.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("heat capacity ratio γ={gamma} must be > 1"),
        });
    }
    if v1 <= 0.0 {
        return Err(UshmaError::InvalidVolume { cubic_meters: v1 });
    }
    if v2 <= 0.0 {
        return Err(UshmaError::InvalidVolume { cubic_meters: v2 });
    }
    Ok(t1 * (v1 / v2).powf(gamma - 1.0))
}

/// Compressibility factor Z = PV/(nRT).
///
/// Z = 1 for ideal gas. Deviations indicate intermolecular forces.
pub fn compressibility_factor(
    pressure: f64,
    volume: f64,
    moles: f64,
    temperature: f64,
) -> Result<f64> {
    let denom = moles * GAS_CONSTANT * temperature;
    if denom.abs() < 1e-30 {
        return Err(UshmaError::DivisionByZero {
            context: "nRT cannot be zero for compressibility factor".into(),
        });
    }
    Ok(pressure * volume / denom)
}

// --- Real gas data ---

/// Critical properties and acentric factor for a gas.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GasData {
    /// Gas name.
    pub name: Cow<'static, str>,
    /// Critical temperature (K).
    pub critical_t: f64,
    /// Critical pressure (Pa).
    pub critical_p: f64,
    /// Acentric factor ω (dimensionless).
    pub acentric_factor: f64,
}

// NIST Chemistry WebBook values
pub const GAS_NITROGEN: GasData = GasData {
    name: Cow::Borrowed("Nitrogen"),
    critical_t: 126.19,
    critical_p: 3_390_000.0,
    acentric_factor: 0.037,
};
pub const GAS_OXYGEN: GasData = GasData {
    name: Cow::Borrowed("Oxygen"),
    critical_t: 154.58,
    critical_p: 5_043_000.0,
    acentric_factor: 0.022,
};
pub const GAS_CO2: GasData = GasData {
    name: Cow::Borrowed("Carbon dioxide"),
    critical_t: 304.13,
    critical_p: 7_375_000.0,
    acentric_factor: 0.224,
};
pub const GAS_METHANE: GasData = GasData {
    name: Cow::Borrowed("Methane"),
    critical_t: 190.56,
    critical_p: 4_599_000.0,
    acentric_factor: 0.011,
};
pub const GAS_WATER: GasData = GasData {
    name: Cow::Borrowed("Water"),
    critical_t: 647.096,
    critical_p: 22_064_000.0,
    acentric_factor: 0.344,
};
pub const GAS_AMMONIA: GasData = GasData {
    name: Cow::Borrowed("Ammonia"),
    critical_t: 405.56,
    critical_p: 11_280_000.0,
    acentric_factor: 0.253,
};
pub const GAS_ETHANE: GasData = GasData {
    name: Cow::Borrowed("Ethane"),
    critical_t: 305.32,
    critical_p: 4_872_000.0,
    acentric_factor: 0.099,
};
pub const GAS_PROPANE: GasData = GasData {
    name: Cow::Borrowed("Propane"),
    critical_t: 369.83,
    critical_p: 4_248_000.0,
    acentric_factor: 0.152,
};

/// All built-in gas data for iteration.
pub const ALL_GASES: &[&GasData] = &[
    &GAS_NITROGEN,
    &GAS_OXYGEN,
    &GAS_CO2,
    &GAS_METHANE,
    &GAS_WATER,
    &GAS_AMMONIA,
    &GAS_ETHANE,
    &GAS_PROPANE,
];

// --- Reduced coordinates ---

/// Reduced temperature: Tr = T/Tc (dimensionless).
#[inline]
#[must_use]
pub fn reduced_temperature(t: f64, t_critical: f64) -> f64 {
    t / t_critical
}

/// Reduced pressure: Pr = P/Pc (dimensionless).
#[inline]
#[must_use]
pub fn reduced_pressure(p: f64, p_critical: f64) -> f64 {
    p / p_critical
}

// --- Redlich-Kwong equation of state ---

/// Redlich-Kwong parameters from critical properties.
///
/// Returns (a, b) where:
/// - a = 0.42748 · R² · Tc^2.5 / Pc
/// - b = 0.08664 · R · Tc / Pc
#[must_use]
pub fn redlich_kwong_params(tc: f64, pc: f64) -> (f64, f64) {
    let a = 0.42748 * GAS_CONSTANT * GAS_CONSTANT * tc.powf(2.5) / pc;
    let b = 0.08664 * GAS_CONSTANT * tc / pc;
    (a, b)
}

/// Redlich-Kwong pressure: P = RT/(Vm - b) - a/(T^0.5 · Vm · (Vm + b)).
///
/// - `temperature`: T (K)
/// - `molar_volume`: Vm (m³/mol)
/// - `tc`: critical temperature (K)
/// - `pc`: critical pressure (Pa)
pub fn redlich_kwong_pressure(
    temperature: f64,
    molar_volume: f64,
    tc: f64,
    pc: f64,
) -> Result<f64> {
    if temperature <= 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    let (a, b) = redlich_kwong_params(tc, pc);
    let vm_b = molar_volume - b;
    if vm_b <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("molar volume {molar_volume} too small for RK (Vm must > b={b:.6})"),
        });
    }
    let repulsive = GAS_CONSTANT * temperature / vm_b;
    let attractive = a / (temperature.sqrt() * molar_volume * (molar_volume + b));
    Ok(repulsive - attractive)
}

// --- Peng-Robinson equation of state ---

/// Peng-Robinson parameters from critical properties and acentric factor.
///
/// Returns (a, b, kappa) where:
/// - a = 0.45724 · R² · Tc² / Pc
/// - b = 0.07780 · R · Tc / Pc
/// - κ = 0.37464 + 1.54226ω - 0.26992ω²
#[must_use]
pub fn peng_robinson_params(tc: f64, pc: f64, omega: f64) -> (f64, f64, f64) {
    let a = 0.45724 * GAS_CONSTANT * GAS_CONSTANT * tc * tc / pc;
    let b = 0.07780 * GAS_CONSTANT * tc / pc;
    let kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
    (a, b, kappa)
}

/// Peng-Robinson pressure.
///
/// P = RT/(Vm - b) - a·α(T) / (Vm² + 2bVm - b²)
/// where α(T) = [1 + κ(1 - √(T/Tc))]²
pub fn peng_robinson_pressure(
    temperature: f64,
    molar_volume: f64,
    tc: f64,
    pc: f64,
    omega: f64,
) -> Result<f64> {
    if temperature <= 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    let (a, b, kappa) = peng_robinson_params(tc, pc, omega);
    let vm_b = molar_volume - b;
    if vm_b <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("molar volume {molar_volume} too small for PR (Vm must > b={b:.6})"),
        });
    }
    let alpha_sqrt = 1.0 + kappa * (1.0 - (temperature / tc).sqrt());
    let alpha = alpha_sqrt * alpha_sqrt;
    let repulsive = GAS_CONSTANT * temperature / vm_b;
    let denom_attr = molar_volume * molar_volume + 2.0 * b * molar_volume - b * b;
    if denom_attr.abs() < 1e-30 {
        return Err(UshmaError::DivisionByZero {
            context: "PR attractive term denominator is zero".into(),
        });
    }
    Ok(repulsive - a * alpha / denom_attr)
}

// --- Virial equation of state ---

/// Pitzer correlation for second virial coefficient B (m³/mol).
///
/// B·Pc/(R·Tc) = B⁰ + ω·B¹
/// - B⁰ = 0.083 - 0.422/Tr^1.6
/// - B¹ = 0.139 - 0.172/Tr^4.2
#[must_use]
pub fn pitzer_second_virial(temperature: f64, tc: f64, pc: f64, omega: f64) -> f64 {
    let tr = temperature / tc;
    let b0 = 0.083 - 0.422 / tr.powf(1.6);
    let b1 = 0.139 - 0.172 / tr.powf(4.2);
    (b0 + omega * b1) * GAS_CONSTANT * tc / pc
}

/// Virial EOS pressure (2nd order): P = RT/Vm · (1 + B/Vm).
pub fn virial_pressure_2nd(temperature: f64, molar_volume: f64, b_coeff: f64) -> Result<f64> {
    if temperature <= 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    if molar_volume <= 0.0 {
        return Err(UshmaError::InvalidVolume {
            cubic_meters: molar_volume,
        });
    }
    let z = 1.0 + b_coeff / molar_volume;
    Ok(z * GAS_CONSTANT * temperature / molar_volume)
}

/// Virial EOS pressure (3rd order): P = RT/Vm · (1 + B/Vm + C/Vm²).
pub fn virial_pressure_3rd(
    temperature: f64,
    molar_volume: f64,
    b_coeff: f64,
    c_coeff: f64,
) -> Result<f64> {
    if temperature <= 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    if molar_volume <= 0.0 {
        return Err(UshmaError::InvalidVolume {
            cubic_meters: molar_volume,
        });
    }
    let vm2 = molar_volume * molar_volume;
    let z = 1.0 + b_coeff / molar_volume + c_coeff / vm2;
    Ok(z * GAS_CONSTANT * temperature / molar_volume)
}

// --- Generalized compressibility ---

/// Pitzer compressibility correlation: Z⁰ + ω·Z¹.
///
/// Simplified correlation valid for Tr > 0.8.
/// Z⁰ and Z¹ from the Lee-Kesler tables (simplified fit).
pub fn compressibility_pitzer(pr: f64, tr: f64, omega: f64) -> Result<f64> {
    if tr <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("reduced temperature {tr} must be positive"),
        });
    }
    // Simplified Lee-Kesler: Z⁰ = 1 + B⁰·Pr/Tr, Z¹ = B¹·Pr/Tr
    let b0 = 0.083 - 0.422 / tr.powf(1.6);
    let b1 = 0.139 - 0.172 / tr.powf(4.2);
    let z0 = 1.0 + b0 * pr / tr;
    let z1 = b1 * pr / tr;
    Ok(z0 + omega * z1)
}

// --- Mixture rules ---

/// Kay's rule for pseudo-critical temperature: Tc_mix = Σ yi·Tc_i (K).
pub fn kays_rule_tc(mole_fractions: &[f64], tc_values: &[f64]) -> Result<f64> {
    validate_mixture(mole_fractions, tc_values)?;
    Ok(mole_fractions
        .iter()
        .zip(tc_values.iter())
        .map(|(y, tc)| y * tc)
        .sum())
}

/// Kay's rule for pseudo-critical pressure: Pc_mix = Σ yi·Pc_i (Pa).
pub fn kays_rule_pc(mole_fractions: &[f64], pc_values: &[f64]) -> Result<f64> {
    validate_mixture(mole_fractions, pc_values)?;
    Ok(mole_fractions
        .iter()
        .zip(pc_values.iter())
        .map(|(y, pc)| y * pc)
        .sum())
}

/// Kay's rule for pseudo-acentric factor: ω_mix = Σ yi·ωi.
pub fn kays_rule_omega(mole_fractions: &[f64], omega_values: &[f64]) -> Result<f64> {
    validate_mixture(mole_fractions, omega_values)?;
    Ok(mole_fractions
        .iter()
        .zip(omega_values.iter())
        .map(|(y, w)| y * w)
        .sum())
}

fn validate_mixture(mole_fractions: &[f64], values: &[f64]) -> Result<()> {
    if mole_fractions.len() != values.len() {
        return Err(UshmaError::InvalidParameter {
            reason: format!(
                "mole_fractions length {} != values length {}",
                mole_fractions.len(),
                values.len()
            ),
        });
    }
    let sum: f64 = mole_fractions.iter().sum();
    if (sum - 1.0).abs() > 1e-10 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("mole fractions must sum to 1.0, got {sum}"),
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ideal_gas_stp() {
        // 1 mol at STP: V ≈ 0.02241 m³
        let v = ideal_gas_volume(1.0, STANDARD_TEMP, ATM).unwrap();
        assert!((v - 0.02241).abs() < 0.001);
    }

    #[test]
    fn test_ideal_gas_pressure() {
        let p = ideal_gas_pressure(1.0, STANDARD_TEMP, 0.02241).unwrap();
        assert!((p - ATM).abs() / ATM < 0.01);
    }

    #[test]
    fn test_ideal_gas_roundtrip() {
        let t = 300.0;
        let n = 2.0;
        let v = ideal_gas_volume(n, t, ATM).unwrap();
        let p = ideal_gas_pressure(n, t, v).unwrap();
        assert!((p - ATM).abs() / ATM < 1e-6);
    }

    #[test]
    fn test_van_der_waals() {
        // CO₂: a = 0.3658 Pa⋅m⁶/mol², b = 4.286e-5 m³/mol
        let p = van_der_waals_pressure(1.0, STANDARD_TEMP, 0.02241, 0.3658, 4.286e-5).unwrap();
        // Should be close to but not exactly ATM
        assert!((p - ATM).abs() / ATM < 0.05);
    }

    #[test]
    fn test_isothermal_work_expansion() {
        // Expansion: V₂ > V₁ → positive work
        let w = isothermal_work(1.0, 300.0, 0.01, 0.02).unwrap();
        assert!(w > 0.0);
    }

    #[test]
    fn test_isothermal_work_compression() {
        // Compression: V₂ < V₁ → negative work
        let w = isothermal_work(1.0, 300.0, 0.02, 0.01).unwrap();
        assert!(w < 0.0);
    }

    #[test]
    fn test_isobaric_work() {
        let w = isobaric_work(ATM, 0.01, 0.02);
        assert!((w - ATM * 0.01).abs() < 1.0);
    }

    #[test]
    fn test_adiabatic_temperature() {
        // Compression → temperature rises
        let t2 = adiabatic_temperature(300.0, 0.02, 0.01, 1.4).unwrap();
        assert!(t2 > 300.0);
    }

    #[test]
    fn test_compressibility_ideal() {
        let z = compressibility_factor(ATM, 0.02241, 1.0, STANDARD_TEMP).unwrap();
        assert!((z - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_negative_temperature() {
        assert!(ideal_gas_pressure(1.0, -10.0, 1.0).is_err());
    }

    #[test]
    fn test_zero_volume() {
        assert!(ideal_gas_pressure(1.0, 300.0, 0.0).is_err());
    }

    #[test]
    fn test_isothermal_work_negative_temp() {
        assert!(isothermal_work(1.0, -10.0, 0.01, 0.02).is_err());
    }

    #[test]
    fn test_isothermal_work_zero_temp() {
        // T=0 is valid (W=0 at absolute zero)
        let w = isothermal_work(1.0, 0.0, 0.01, 0.02).unwrap();
        assert!(w.abs() < 1e-30);
    }

    #[test]
    fn test_adiabatic_invalid_gamma() {
        // γ must be > 1
        assert!(adiabatic_temperature(300.0, 0.02, 0.01, 1.0).is_err());
        assert!(adiabatic_temperature(300.0, 0.02, 0.01, 0.5).is_err());
    }

    #[test]
    fn test_adiabatic_roundtrip() {
        // Compress then expand back — temperature should return
        let t2 = adiabatic_temperature(300.0, 0.02, 0.01, 1.4).unwrap();
        let t3 = adiabatic_temperature(t2, 0.01, 0.02, 1.4).unwrap();
        assert!((t3 - 300.0).abs() < 1e-10);
    }

    #[test]
    fn test_ideal_gas_temperature_zero_moles() {
        assert!(ideal_gas_temperature(ATM, 0.02241, 0.0).is_err());
    }

    #[test]
    fn test_van_der_waals_negative_temp() {
        assert!(van_der_waals_pressure(1.0, -10.0, 0.02241, 0.3658, 4.286e-5).is_err());
    }

    #[test]
    fn test_van_der_waals_volume_too_small() {
        // Volume < n*b → effective volume negative
        assert!(van_der_waals_pressure(1.0, 300.0, 1e-6, 0.3658, 4.286e-5).is_err());
    }

    #[test]
    fn test_compressibility_factor_zero_temp() {
        assert!(compressibility_factor(ATM, 0.02241, 1.0, 0.0).is_err());
    }

    #[test]
    fn test_ideal_gas_zero_temp() {
        // T=0 → P=0, which is physically valid
        let p = ideal_gas_pressure(1.0, 0.0, 0.02241).unwrap();
        assert!(p.abs() < 1e-30);
    }

    // --- Gas data ---

    #[test]
    fn test_all_gases_count() {
        assert_eq!(ALL_GASES.len(), 8);
    }

    #[test]
    fn test_gas_data_valid() {
        for g in ALL_GASES {
            assert!(g.critical_t > 0.0, "{} Tc", g.name);
            assert!(g.critical_p > 0.0, "{} Pc", g.name);
            assert!(g.acentric_factor >= 0.0, "{} ω", g.name);
        }
    }

    #[test]
    fn test_gas_data_serde_roundtrip() {
        let json = serde_json::to_string(&GAS_CO2).unwrap();
        let back: GasData = serde_json::from_str(&json).unwrap();
        assert_eq!(back.name, "Carbon dioxide");
        assert!((back.critical_t - 304.13).abs() < 0.01);
    }

    #[test]
    fn test_reduced_coordinates() {
        let tr = reduced_temperature(300.0, 304.13);
        assert!((tr - 300.0 / 304.13).abs() < 1e-10);
        let pr = reduced_pressure(5_000_000.0, 7_375_000.0);
        assert!((pr - 5_000_000.0 / 7_375_000.0).abs() < 1e-10);
    }

    // --- Redlich-Kwong ---

    #[test]
    fn test_rk_params() {
        let (a, b) = redlich_kwong_params(304.13, 7_375_000.0);
        assert!(a > 0.0);
        assert!(b > 0.0);
    }

    #[test]
    fn test_rk_pressure_co2_stp() {
        // CO₂ at STP: 1 mol in 0.02241 m³
        let p = redlich_kwong_pressure(STANDARD_TEMP, 0.02241, 304.13, 7_375_000.0).unwrap();
        // Should be close to ATM, closer than vdW
        assert!((p - ATM).abs() / ATM < 0.05);
    }

    #[test]
    fn test_rk_invalid() {
        assert!(redlich_kwong_pressure(0.0, 0.02241, 304.13, 7_375_000.0).is_err());
        assert!(redlich_kwong_pressure(300.0, 1e-7, 304.13, 7_375_000.0).is_err()); // Vm < b
    }

    // --- Peng-Robinson ---

    #[test]
    fn test_pr_params() {
        let (a, b, kappa) = peng_robinson_params(304.13, 7_375_000.0, 0.224);
        assert!(a > 0.0);
        assert!(b > 0.0);
        assert!(kappa > 0.0);
    }

    #[test]
    fn test_pr_pressure_co2_stp() {
        let p = peng_robinson_pressure(STANDARD_TEMP, 0.02241, 304.13, 7_375_000.0, 0.224).unwrap();
        assert!((p - ATM).abs() / ATM < 0.05);
    }

    #[test]
    fn test_pr_vs_rk_differ() {
        // At moderate molar volume, PR and RK give different pressures
        let vm = 0.005; // compressed but above b for both EOS
        let p_rk = redlich_kwong_pressure(400.0, vm, 304.13, 7_375_000.0).unwrap();
        let p_pr = peng_robinson_pressure(400.0, vm, 304.13, 7_375_000.0, 0.224).unwrap();
        // Both positive, but different
        assert!(p_rk > 0.0);
        assert!(p_pr > 0.0);
        assert!((p_rk - p_pr).abs() > 100.0);
    }

    #[test]
    fn test_pr_invalid() {
        assert!(peng_robinson_pressure(0.0, 0.02241, 304.13, 7_375_000.0, 0.224).is_err());
    }

    // --- Virial ---

    #[test]
    fn test_pitzer_second_virial() {
        let b = pitzer_second_virial(300.0, 304.13, 7_375_000.0, 0.224);
        // B should be negative at T near Tc (attractive forces dominate)
        assert!(b < 0.0);
    }

    #[test]
    fn test_virial_ideal_at_large_vm() {
        // Large Vm → B/Vm → 0 → P ≈ RT/Vm (ideal gas)
        let vm = 1.0; // very large molar volume
        let b = pitzer_second_virial(300.0, 304.13, 7_375_000.0, 0.224);
        let p_virial = virial_pressure_2nd(300.0, vm, b).unwrap();
        let p_ideal = GAS_CONSTANT * 300.0 / vm;
        assert!((p_virial - p_ideal).abs() / p_ideal < 0.01);
    }

    #[test]
    fn test_virial_3rd_order() {
        let p = virial_pressure_3rd(300.0, 0.02241, -1e-4, 1e-8).unwrap();
        assert!(p > 0.0);
    }

    #[test]
    fn test_virial_invalid() {
        assert!(virial_pressure_2nd(0.0, 0.02241, -1e-4).is_err());
        assert!(virial_pressure_2nd(300.0, 0.0, -1e-4).is_err());
    }

    // --- Compressibility + mixtures ---

    #[test]
    fn test_compressibility_pitzer_ideal() {
        // At low Pr, Z → 1
        let z = compressibility_pitzer(0.01, 2.0, 0.0).unwrap();
        assert!((z - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_compressibility_pitzer_invalid() {
        assert!(compressibility_pitzer(1.0, 0.0, 0.1).is_err());
    }

    #[test]
    fn test_kays_rule_pure() {
        // Pure component → identity
        let tc = kays_rule_tc(&[1.0], &[304.13]).unwrap();
        assert!((tc - 304.13).abs() < 1e-10);
    }

    #[test]
    fn test_kays_rule_binary() {
        // 50/50 N₂/O₂
        let tc = kays_rule_tc(&[0.5, 0.5], &[126.19, 154.58]).unwrap();
        assert!((tc - 140.385).abs() < 0.01);
    }

    #[test]
    fn test_kays_rule_invalid() {
        // Fractions don't sum to 1
        assert!(kays_rule_tc(&[0.3, 0.3], &[126.19, 154.58]).is_err());
        // Length mismatch
        assert!(kays_rule_tc(&[0.5, 0.5], &[126.19]).is_err());
    }
}
