//! Heat transfer — conduction, convection, radiation.
//!
//! All SI units: watts, meters, kelvins, seconds.

use serde::{Deserialize, Serialize};

use crate::error::{UshmaError, Result};

/// Stefan-Boltzmann constant σ (W/(m²⋅K⁴)).
pub const STEFAN_BOLTZMANN: f64 = 5.670_374_419e-8;

/// Boltzmann constant k_B (J/K).
pub const BOLTZMANN_K: f64 = 1.380_649e-23;

/// Fourier's law: heat flux through conduction.
///
/// q = -k⋅A⋅(T_hot - T_cold)/L (watts)
///
/// - `conductivity`: thermal conductivity k (W/(m⋅K))
/// - `area`: cross-sectional area A (m²)
/// - `t_hot`, `t_cold`: temperatures (K)
/// - `thickness`: material thickness L (m)
pub fn conduction(
    conductivity: f64,
    area: f64,
    t_hot: f64,
    t_cold: f64,
    thickness: f64,
) -> Result<f64> {
    if conductivity < 0.0 {
        return Err(UshmaError::InvalidConductivity {
            value: conductivity,
        });
    }
    if thickness.abs() < 1e-30 {
        return Err(UshmaError::DivisionByZero {
            context: "thickness cannot be zero for conduction".into(),
        });
    }
    Ok(conductivity * area * (t_hot - t_cold) / thickness)
}

/// Newton's law of cooling: convective heat transfer.
///
/// q = h⋅A⋅(T_surface - T_fluid) (watts)
///
/// - `h`: convective heat transfer coefficient (W/(m²⋅K))
/// - `area`: surface area (m²)
/// - `t_surface`, `t_fluid`: temperatures (K)
#[inline]
#[must_use]
pub fn convection(h: f64, area: f64, t_surface: f64, t_fluid: f64) -> f64 {
    h * area * (t_surface - t_fluid)
}

/// Stefan-Boltzmann law: radiative heat transfer.
///
/// q = ε⋅σ⋅A⋅(T⁴ - T_surr⁴) (watts)
///
/// - `emissivity`: surface emissivity ε (0-1)
/// - `area`: surface area (m²)
/// - `t_surface`, `t_surrounding`: temperatures (K)
pub fn radiation(
    emissivity: f64,
    area: f64,
    t_surface: f64,
    t_surrounding: f64,
) -> Result<f64> {
    if t_surface < 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: t_surface,
        });
    }
    if t_surrounding < 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: t_surrounding,
        });
    }
    let t_s4 = t_surface.powi(4);
    let t_r4 = t_surrounding.powi(4);
    Ok(emissivity * STEFAN_BOLTZMANN * area * (t_s4 - t_r4))
}

/// Thermal resistance for conduction: R = L/(k⋅A) (K/W).
pub fn thermal_resistance_conduction(
    conductivity: f64,
    area: f64,
    thickness: f64,
) -> Result<f64> {
    let denom = conductivity * area;
    if denom.abs() < 1e-30 {
        return Err(UshmaError::DivisionByZero {
            context: "k⋅A cannot be zero for thermal resistance".into(),
        });
    }
    Ok(thickness / denom)
}

/// Thermal resistance for convection: R = 1/(h⋅A) (K/W).
pub fn thermal_resistance_convection(h: f64, area: f64) -> Result<f64> {
    let denom = h * area;
    if denom.abs() < 1e-30 {
        return Err(UshmaError::DivisionByZero {
            context: "h⋅A cannot be zero for thermal resistance".into(),
        });
    }
    Ok(1.0 / denom)
}

/// Series thermal resistance: R_total = R₁ + R₂ + ... (K/W).
#[must_use]
pub fn thermal_resistance_series(resistances: &[f64]) -> f64 {
    resistances.iter().sum()
}

/// Parallel thermal resistance: 1/R_total = 1/R₁ + 1/R₂ + ... (K/W).
pub fn thermal_resistance_parallel(resistances: &[f64]) -> Result<f64> {
    let sum: f64 = resistances.iter().map(|r| 1.0 / r).sum();
    if sum.abs() < 1e-30 {
        return Err(UshmaError::DivisionByZero {
            context: "parallel resistance sum is zero".into(),
        });
    }
    Ok(1.0 / sum)
}

/// Heat stored in a body: Q = mcΔT (joules).
#[inline]
#[must_use]
pub fn heat_stored(mass: f64, specific_heat: f64, delta_t: f64) -> f64 {
    mass * specific_heat * delta_t
}

/// Thermal diffusivity: α = k/(ρ⋅c_p) (m²/s).
///
/// How quickly heat diffuses through a material.
pub fn thermal_diffusivity(
    conductivity: f64,
    density: f64,
    specific_heat: f64,
) -> Result<f64> {
    let denom = density * specific_heat;
    if denom.abs() < 1e-30 {
        return Err(UshmaError::DivisionByZero {
            context: "ρ⋅c_p cannot be zero for thermal diffusivity".into(),
        });
    }
    Ok(conductivity / denom)
}

/// Biot number: Bi = hL/k (dimensionless).
///
/// Ratio of conduction resistance to convection resistance.
/// Bi < 0.1 → lumped capacitance model valid.
pub fn biot_number(h: f64, characteristic_length: f64, conductivity: f64) -> Result<f64> {
    if conductivity.abs() < 1e-30 {
        return Err(UshmaError::DivisionByZero {
            context: "conductivity cannot be zero for Biot number".into(),
        });
    }
    Ok(h * characteristic_length / conductivity)
}

/// Lumped capacitance cooling: T(t) = T_env + (T_0 - T_env)⋅e^(-t/τ).
///
/// Valid when Bi < 0.1.
/// - `t_initial`: starting temperature (K)
/// - `t_environment`: ambient temperature (K)
/// - `time`: elapsed time (s)
/// - `time_constant`: τ = ρVc/(hA) (s)
#[inline]
#[must_use]
pub fn lumped_capacitance(
    t_initial: f64,
    t_environment: f64,
    time: f64,
    time_constant: f64,
) -> f64 {
    t_environment + (t_initial - t_environment) * (-time / time_constant).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_conduction_positive_flux() {
        let q = conduction(401.0, 0.01, 373.15, 293.15, 0.1).unwrap();
        // Copper, 1cm² area, 80K difference, 10cm thick
        assert!(q > 0.0);
    }

    #[test]
    fn test_conduction_zero_gradient() {
        let q = conduction(401.0, 0.01, 300.0, 300.0, 0.1).unwrap();
        assert!(q.abs() < 1e-10);
    }

    #[test]
    fn test_convection() {
        let q = convection(25.0, 1.0, 373.15, 293.15);
        // h=25, A=1m², ΔT=80K → q=2000W
        assert!((q - 2000.0).abs() < 1.0);
    }

    #[test]
    fn test_radiation() {
        let q = radiation(1.0, 1.0, 373.15, 293.15).unwrap();
        // Blackbody, 1m², 100°C to 20°C
        assert!(q > 0.0);
    }

    #[test]
    fn test_radiation_negative_temp() {
        assert!(radiation(1.0, 1.0, -10.0, 293.15).is_err());
    }

    #[test]
    fn test_thermal_resistance_series() {
        let r = thermal_resistance_series(&[1.0, 2.0, 3.0]);
        assert!((r - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_thermal_resistance_parallel() {
        let r = thermal_resistance_parallel(&[2.0, 2.0]).unwrap();
        assert!((r - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_heat_stored() {
        // 1 kg water, 4186 J/(kg⋅K), 10K rise → 41860 J
        let q = heat_stored(1.0, 4186.0, 10.0);
        assert!((q - 41860.0).abs() < 1.0);
    }

    #[test]
    fn test_thermal_diffusivity() {
        // Copper: k=401, ρ=8960, c=385 → α ≈ 1.16e-4 m²/s
        let alpha = thermal_diffusivity(401.0, 8960.0, 385.0).unwrap();
        assert!((alpha - 1.16e-4).abs() / 1.16e-4 < 0.01);
    }

    #[test]
    fn test_biot_number() {
        let bi = biot_number(25.0, 0.01, 401.0).unwrap();
        // h=25, L=1cm, k=401 → Bi ≈ 0.0006 (lumped valid)
        assert!(bi < 0.1);
    }

    #[test]
    fn test_lumped_capacitance() {
        let t0 = 373.15;
        let t_env = 293.15;
        let tau = 100.0;
        // At t=0, temperature equals initial
        let t = lumped_capacitance(t0, t_env, 0.0, tau);
        assert!((t - t0).abs() < 1e-10);
        // At t→∞, temperature approaches environment
        let t_inf = lumped_capacitance(t0, t_env, 10000.0, tau);
        assert!((t_inf - t_env).abs() < 1.0);
    }

    #[test]
    fn test_invalid_conductivity() {
        assert!(conduction(-1.0, 1.0, 373.0, 293.0, 0.1).is_err());
    }
}
