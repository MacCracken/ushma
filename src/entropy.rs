//! Entropy, free energy, thermodynamic potentials, equilibrium.
//!
//! The second law: total entropy of an isolated system never decreases.

use crate::error::{Result, UshmaError};
use crate::state::GAS_CONSTANT;

/// Entropy change for ideal gas: ΔS = nCv⋅ln(T₂/T₁) + nR⋅ln(V₂/V₁).
pub fn ideal_gas_entropy_change(
    moles: f64,
    cv: f64,
    t1: f64,
    t2: f64,
    v1: f64,
    v2: f64,
) -> Result<f64> {
    if t1 <= 0.0 {
        return Err(UshmaError::InvalidTemperature { kelvin: t1 });
    }
    if t2 <= 0.0 {
        return Err(UshmaError::InvalidTemperature { kelvin: t2 });
    }
    if v1 <= 0.0 {
        return Err(UshmaError::InvalidVolume { cubic_meters: v1 });
    }
    if v2 <= 0.0 {
        return Err(UshmaError::InvalidVolume { cubic_meters: v2 });
    }
    Ok(moles * (cv * (t2 / t1).ln() + GAS_CONSTANT * (v2 / v1).ln()))
}

/// Entropy change for isothermal process: ΔS = nR⋅ln(V₂/V₁).
pub fn isothermal_entropy_change(moles: f64, v1: f64, v2: f64) -> Result<f64> {
    if v1 <= 0.0 {
        return Err(UshmaError::InvalidVolume { cubic_meters: v1 });
    }
    if v2 <= 0.0 {
        return Err(UshmaError::InvalidVolume { cubic_meters: v2 });
    }
    Ok(moles * GAS_CONSTANT * (v2 / v1).ln())
}

/// Entropy change for heat transfer: ΔS = Q/T.
pub fn heat_transfer_entropy(heat: f64, temperature: f64) -> Result<f64> {
    if temperature <= 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    Ok(heat / temperature)
}

/// Carnot efficiency: η = 1 - T_cold/T_hot.
#[tracing::instrument(level = "debug")]
pub fn carnot_efficiency(t_hot: f64, t_cold: f64) -> Result<f64> {
    if t_hot <= 0.0 {
        return Err(UshmaError::InvalidTemperature { kelvin: t_hot });
    }
    if t_cold <= 0.0 {
        return Err(UshmaError::InvalidTemperature { kelvin: t_cold });
    }
    if t_cold >= t_hot {
        return Err(UshmaError::InvalidParameter {
            reason: format!("T_cold ({t_cold} K) must be less than T_hot ({t_hot} K)"),
        });
    }
    Ok(1.0 - t_cold / t_hot)
}

/// Carnot COP for refrigeration: COP = T_cold/(T_hot - T_cold).
pub fn carnot_cop_refrigeration(t_hot: f64, t_cold: f64) -> Result<f64> {
    if t_hot <= 0.0 {
        return Err(UshmaError::InvalidTemperature { kelvin: t_hot });
    }
    if t_cold <= 0.0 {
        return Err(UshmaError::InvalidTemperature { kelvin: t_cold });
    }
    let dt = t_hot - t_cold;
    if dt <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: "T_hot must exceed T_cold".into(),
        });
    }
    Ok(t_cold / dt)
}

/// Helmholtz free energy: A = U - TS (joules).
#[inline]
#[must_use]
pub fn helmholtz(internal_energy: f64, temperature: f64, entropy: f64) -> f64 {
    internal_energy - temperature * entropy
}

/// Gibbs free energy: G = H - TS (joules).
#[inline]
#[must_use]
pub fn gibbs(enthalpy: f64, temperature: f64, entropy: f64) -> f64 {
    enthalpy - temperature * entropy
}

/// Clausius inequality check: ΔS_total ≥ 0 for irreversible processes.
#[inline]
#[must_use]
pub fn is_spontaneous(delta_s_total: f64) -> bool {
    delta_s_total >= 0.0
}

/// Entropy of mixing for ideal gases: ΔS_mix = -nR Σ xᵢ ln(xᵢ).
#[tracing::instrument(level = "debug")]
pub fn entropy_of_mixing(moles: f64, mole_fractions: &[f64]) -> Result<f64> {
    for &x in mole_fractions {
        if !(0.0..=1.0).contains(&x) {
            return Err(UshmaError::InvalidParameter {
                reason: format!("mole fraction {x} must be in [0, 1]"),
            });
        }
    }
    let frac_sum: f64 = mole_fractions.iter().sum();
    if (frac_sum - 1.0).abs() > 1e-10 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("mole fractions must sum to 1.0, got {frac_sum}"),
        });
    }
    // x=0 contributes nothing: lim x→0+ of x·ln(x) = 0
    let sum: f64 = mole_fractions
        .iter()
        .filter(|&&x| x > 0.0)
        .map(|&x| x * x.ln())
        .sum();
    Ok(-moles * GAS_CONSTANT * sum)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_carnot_efficiency() {
        let eta = carnot_efficiency(500.0, 300.0).unwrap();
        assert!((eta - 0.4).abs() < 1e-10);
    }

    #[test]
    fn test_carnot_cop() {
        let cop = carnot_cop_refrigeration(300.0, 250.0).unwrap();
        assert!((cop - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_isothermal_entropy_expansion() {
        // Expansion → positive entropy change
        let ds = isothermal_entropy_change(1.0, 0.01, 0.02).unwrap();
        assert!(ds > 0.0);
    }

    #[test]
    fn test_heat_transfer_entropy() {
        let ds = heat_transfer_entropy(1000.0, 300.0).unwrap();
        assert!((ds - 1000.0 / 300.0).abs() < 1e-10);
    }

    #[test]
    fn test_helmholtz() {
        let a = helmholtz(1000.0, 300.0, 5.0);
        assert!((a - (1000.0 - 1500.0)).abs() < 1e-10);
    }

    #[test]
    fn test_gibbs() {
        let g = gibbs(2000.0, 300.0, 5.0);
        assert!((g - (2000.0 - 1500.0)).abs() < 1e-10);
    }

    #[test]
    fn test_spontaneous() {
        assert!(is_spontaneous(0.1));
        assert!(is_spontaneous(0.0));
        assert!(!is_spontaneous(-0.1));
    }

    #[test]
    fn test_entropy_of_mixing() {
        // Equal mixture of two gases
        let ds = entropy_of_mixing(1.0, &[0.5, 0.5]).unwrap();
        assert!(ds > 0.0); // mixing always increases entropy
    }

    #[test]
    fn test_ideal_gas_entropy_change() {
        // Isothermal expansion (T constant, V doubles)
        let ds = ideal_gas_entropy_change(1.0, 20.8, 300.0, 300.0, 0.01, 0.02).unwrap();
        assert!(ds > 0.0);
    }

    #[test]
    fn test_carnot_invalid_temps() {
        assert!(carnot_efficiency(300.0, 300.0).is_err());
        assert!(carnot_efficiency(200.0, 300.0).is_err());
    }

    #[test]
    fn test_carnot_cop_invalid() {
        assert!(carnot_cop_refrigeration(300.0, 300.0).is_err());
        assert!(carnot_cop_refrigeration(200.0, 300.0).is_err());
        assert!(carnot_cop_refrigeration(0.0, 250.0).is_err());
    }

    #[test]
    fn test_heat_transfer_entropy_zero_temp() {
        assert!(heat_transfer_entropy(1000.0, 0.0).is_err());
    }

    #[test]
    fn test_heat_transfer_entropy_negative_heat() {
        // Heat rejected: Q < 0 → ΔS < 0
        let ds = heat_transfer_entropy(-1000.0, 300.0).unwrap();
        assert!(ds < 0.0);
    }

    #[test]
    fn test_entropy_of_mixing_bad_fractions() {
        // Fractions don't sum to 1
        assert!(entropy_of_mixing(1.0, &[0.3, 0.3]).is_err());
    }

    #[test]
    fn test_entropy_of_mixing_zero_fraction() {
        // x=0 is valid (absent component contributes nothing)
        let ds = entropy_of_mixing(1.0, &[0.0, 1.0]).unwrap();
        // Pure substance: only x=1 contributes, ln(1)=0, so ΔS=0
        assert!(ds.abs() < 1e-10);
    }

    #[test]
    fn test_entropy_of_mixing_pure_substance() {
        // Single component: x=1, ln(1)=0, ΔS=0
        let ds = entropy_of_mixing(1.0, &[1.0]).unwrap();
        assert!(ds.abs() < 1e-10);
    }

    #[test]
    fn test_ideal_gas_entropy_zero_temp() {
        assert!(ideal_gas_entropy_change(1.0, 20.8, 0.0, 300.0, 0.01, 0.02).is_err());
        assert!(ideal_gas_entropy_change(1.0, 20.8, 300.0, 0.0, 0.01, 0.02).is_err());
    }

    #[test]
    fn test_ideal_gas_entropy_zero_volume() {
        assert!(ideal_gas_entropy_change(1.0, 20.8, 300.0, 300.0, 0.0, 0.02).is_err());
        assert!(ideal_gas_entropy_change(1.0, 20.8, 300.0, 300.0, 0.01, 0.0).is_err());
    }

    #[test]
    fn test_isothermal_entropy_compression() {
        // Compression: V₂ < V₁ → negative entropy change
        let ds = isothermal_entropy_change(1.0, 0.02, 0.01).unwrap();
        assert!(ds < 0.0);
    }

    #[test]
    fn test_is_spontaneous_boundary() {
        // Reversible process: ΔS = 0 is still "spontaneous" (not forbidden)
        assert!(is_spontaneous(0.0));
    }

    #[test]
    fn test_entropy_of_mixing_three_components() {
        let ds = entropy_of_mixing(1.0, &[0.2, 0.3, 0.5]).unwrap();
        assert!(ds > 0.0);
    }
}
