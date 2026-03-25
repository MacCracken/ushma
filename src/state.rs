//! Equations of state — ideal gas, van der Waals, phase diagrams.
//!
//! All SI units: pascals, cubic meters, kelvins, moles.

use serde::{Deserialize, Serialize};

use crate::error::{UshmaError, Result};

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
}
