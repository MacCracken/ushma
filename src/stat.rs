//! Statistical thermodynamics — Boltzmann, partition functions, speed distributions, solid Cv.
//!
//! Bridges macroscopic thermodynamics to molecular physics.
//! All SI units: joules, kelvins, kilograms, meters, seconds.

use crate::error::{Result, UshmaError};
use crate::state::GAS_CONSTANT;
use crate::transfer::BOLTZMANN_K;

/// Planck constant h (J·s).
pub const PLANCK: f64 = 6.626_070_15e-34;

/// Avogadro number (1/mol).
pub const AVOGADRO: f64 = 6.022_140_76e23;

// --- Boltzmann distribution ---

/// Boltzmann factor: exp(-E/(kT)).
///
/// Proportional to the probability of a state with energy E at temperature T.
pub fn boltzmann_probability(energy: f64, temperature: f64) -> Result<f64> {
    if temperature <= 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    Ok((-energy / (BOLTZMANN_K * temperature)).exp())
}

/// Population ratio of two energy levels: N₂/N₁ = exp(-(E₂-E₁)/(kT)).
pub fn boltzmann_ratio(e1: f64, e2: f64, temperature: f64) -> Result<f64> {
    if temperature <= 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    Ok((-(e2 - e1) / (BOLTZMANN_K * temperature)).exp())
}

/// Entropy from microstates: S = k·ln(W) (J/K).
pub fn entropy_from_microstates(num_microstates: f64) -> Result<f64> {
    if num_microstates <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("microstates {num_microstates} must be positive"),
        });
    }
    Ok(BOLTZMANN_K * num_microstates.ln())
}

// --- Equipartition ---

/// Equipartition energy per molecule: E = (f/2)·kT (J).
#[inline]
#[must_use]
pub fn equipartition_energy(degrees_of_freedom: f64, temperature: f64) -> f64 {
    0.5 * degrees_of_freedom * BOLTZMANN_K * temperature
}

/// Equipartition molar heat capacity: Cv = (f/2)·R (J/(mol·K)).
#[inline]
#[must_use]
pub fn equipartition_cv(degrees_of_freedom: f64) -> f64 {
    0.5 * degrees_of_freedom * GAS_CONSTANT
}

// --- Partition functions ---

/// Translational partition function: q_trans = V·(2πmkT/h²)^(3/2).
///
/// - `mass`: molecular mass (kg)
/// - `volume`: container volume (m³)
/// - `temperature`: T (K)
pub fn partition_translational(mass: f64, volume: f64, temperature: f64) -> Result<f64> {
    if mass <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("mass {mass} kg must be positive"),
        });
    }
    if volume <= 0.0 {
        return Err(UshmaError::InvalidVolume {
            cubic_meters: volume,
        });
    }
    if temperature <= 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    let lambda_sq =
        PLANCK * PLANCK / (2.0 * std::f64::consts::PI * mass * BOLTZMANN_K * temperature);
    Ok(volume / lambda_sq.powf(1.5))
}

/// Rotational partition function (linear molecule, high-T limit).
///
/// q_rot = T / (σ · Θ_rot) where Θ_rot = h²/(8π²Ik).
///
/// - `moment_of_inertia`: I (kg·m²)
/// - `temperature`: T (K)
/// - `symmetry_number`: σ (1 for heteronuclear, 2 for homonuclear)
pub fn partition_rotational(
    moment_of_inertia: f64,
    temperature: f64,
    symmetry_number: u32,
) -> Result<f64> {
    if moment_of_inertia <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("moment of inertia {moment_of_inertia} must be positive"),
        });
    }
    if temperature <= 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    if symmetry_number == 0 {
        return Err(UshmaError::InvalidParameter {
            reason: "symmetry number must be >= 1".into(),
        });
    }
    let theta_rot = PLANCK * PLANCK
        / (8.0 * std::f64::consts::PI * std::f64::consts::PI * moment_of_inertia * BOLTZMANN_K);
    Ok(temperature / (symmetry_number as f64 * theta_rot))
}

/// Vibrational partition function: q_vib = 1/(1 - exp(-hν/(kT))).
///
/// - `frequency`: vibrational frequency ν (Hz)
/// - `temperature`: T (K)
pub fn partition_vibrational(frequency: f64, temperature: f64) -> Result<f64> {
    if frequency <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("frequency {frequency} Hz must be positive"),
        });
    }
    if temperature <= 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    let x = PLANCK * frequency / (BOLTZMANN_K * temperature);
    let denom = 1.0 - (-x).exp();
    if denom.abs() < 1e-30 {
        return Err(UshmaError::DivisionByZero {
            context: "vibrational partition function denominator is zero".into(),
        });
    }
    Ok(1.0 / denom)
}

// --- Maxwell-Boltzmann speed distribution ---

/// Maxwell-Boltzmann speed probability density: f(v) = 4π·(m/(2πkT))^(3/2)·v²·exp(-mv²/(2kT)).
///
/// - `speed`: v (m/s)
/// - `mass`: molecular mass (kg)
/// - `temperature`: T (K)
pub fn maxwell_boltzmann_speed_pdf(speed: f64, mass: f64, temperature: f64) -> Result<f64> {
    if speed < 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("speed {speed} m/s must be non-negative"),
        });
    }
    if mass <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("mass {mass} kg must be positive"),
        });
    }
    if temperature <= 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    let a = mass / (2.0 * std::f64::consts::PI * BOLTZMANN_K * temperature);
    let coeff = 4.0 * std::f64::consts::PI * a.powf(1.5);
    let exponent = -mass * speed * speed / (2.0 * BOLTZMANN_K * temperature);
    Ok(coeff * speed * speed * exponent.exp())
}

/// Most probable speed: v_p = √(2kT/m) (m/s).
pub fn most_probable_speed(mass: f64, temperature: f64) -> Result<f64> {
    validate_mass_temp(mass, temperature)?;
    Ok((2.0 * BOLTZMANN_K * temperature / mass).sqrt())
}

/// Mean speed: v_avg = √(8kT/(πm)) (m/s).
pub fn mean_speed(mass: f64, temperature: f64) -> Result<f64> {
    validate_mass_temp(mass, temperature)?;
    Ok((8.0 * BOLTZMANN_K * temperature / (std::f64::consts::PI * mass)).sqrt())
}

/// Root-mean-square speed: v_rms = √(3kT/m) (m/s).
pub fn rms_speed(mass: f64, temperature: f64) -> Result<f64> {
    validate_mass_temp(mass, temperature)?;
    Ok((3.0 * BOLTZMANN_K * temperature / mass).sqrt())
}

fn validate_mass_temp(mass: f64, temperature: f64) -> Result<()> {
    if mass <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("mass {mass} kg must be positive"),
        });
    }
    if temperature <= 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    Ok(())
}

// --- Solid heat capacity models ---

/// Einstein model for solid heat capacity: Cv = 3R·(Θ_E/T)²·exp(Θ_E/T)/(exp(Θ_E/T)-1)² (J/(mol·K)).
///
/// - `einstein_temp`: Θ_E (K) — Einstein characteristic temperature
/// - `temperature`: T (K)
pub fn einstein_cv(einstein_temp: f64, temperature: f64) -> Result<f64> {
    if einstein_temp <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("Einstein temperature {einstein_temp} K must be positive"),
        });
    }
    if temperature <= 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    let x = einstein_temp / temperature;
    // For very large x (low T), exp(x) overflows — Cv → 0
    if x > 500.0 {
        return Ok(0.0);
    }
    let ex = x.exp();
    let denom = (ex - 1.0) * (ex - 1.0);
    if denom.abs() < 1e-30 {
        return Ok(3.0 * GAS_CONSTANT); // high-T limit
    }
    Ok(3.0 * GAS_CONSTANT * x * x * ex / denom)
}

/// Debye model for solid heat capacity (J/(mol·K)).
///
/// Cv = 9R·(T/Θ_D)³·∫₀^(Θ_D/T) x⁴eˣ/(eˣ-1)² dx
///
/// Uses [`hisab::calc::integral_simpson`] (200 intervals).
pub fn debye_cv(debye_temp: f64, temperature: f64) -> Result<f64> {
    if debye_temp <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("Debye temperature {debye_temp} K must be positive"),
        });
    }
    if temperature <= 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    let td_t = debye_temp / temperature;

    // High-T limit: Cv → 3R
    if td_t < 0.01 {
        return Ok(3.0 * GAS_CONSTANT);
    }

    // Numerical integration via hisab Simpson's rule
    let integrand = |x: f64| -> f64 {
        if x < 1e-30 {
            return 0.0; // limit as x→0
        }
        let ex = x.exp();
        let denom = (ex - 1.0) * (ex - 1.0);
        if denom < 1e-300 {
            return 0.0;
        }
        x * x * x * x * ex / denom
    };

    let integral = hisab::calc::integral_simpson(integrand, 0.0, td_t, 200).unwrap_or(0.0);

    let ratio = temperature / debye_temp;
    Ok(9.0 * GAS_CONSTANT * ratio * ratio * ratio * integral)
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- Constants ---

    #[test]
    fn test_constants() {
        let h = PLANCK;
        let na = AVOGADRO;
        assert!(h > 6.6e-34 && h < 6.7e-34);
        assert!(na > 6.02e23 && na < 6.03e23);
        // k_B * N_A ≈ R
        assert!((BOLTZMANN_K * na - GAS_CONSTANT).abs() / GAS_CONSTANT < 1e-6);
    }

    // --- Boltzmann ---

    #[test]
    fn test_boltzmann_ground_state() {
        // E=0 → P = 1
        let p = boltzmann_probability(0.0, 300.0).unwrap();
        assert!((p - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_boltzmann_ratio_equal_energies() {
        let r = boltzmann_ratio(1.0, 1.0, 300.0).unwrap();
        assert!((r - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_boltzmann_ratio_higher_energy_less_populated() {
        let r = boltzmann_ratio(0.0, 1e-20, 300.0).unwrap();
        assert!(r < 1.0);
    }

    #[test]
    fn test_boltzmann_high_temp_equal() {
        // Very high T → ratio → 1
        let r = boltzmann_ratio(0.0, 1e-20, 1e10).unwrap();
        assert!((r - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_boltzmann_invalid() {
        assert!(boltzmann_probability(0.0, 0.0).is_err());
        assert!(boltzmann_ratio(0.0, 1.0, 0.0).is_err());
    }

    #[test]
    fn test_entropy_from_microstates() {
        // S = k·ln(1) = 0
        let s = entropy_from_microstates(1.0).unwrap();
        assert!(s.abs() < 1e-30);
        // S increases with W
        let s2 = entropy_from_microstates(100.0).unwrap();
        assert!(s2 > 0.0);
    }

    #[test]
    fn test_entropy_microstates_invalid() {
        assert!(entropy_from_microstates(0.0).is_err());
        assert!(entropy_from_microstates(-1.0).is_err());
    }

    // --- Equipartition ---

    #[test]
    fn test_equipartition_monatomic() {
        // f=3 → Cv = 3R/2 = 12.47 J/(mol·K)
        let cv = equipartition_cv(3.0);
        assert!((cv - 1.5 * GAS_CONSTANT).abs() < 1e-10);
    }

    #[test]
    fn test_equipartition_diatomic() {
        // f=5 → Cv = 5R/2 = 20.79 J/(mol·K)
        let cv = equipartition_cv(5.0);
        assert!((cv - 2.5 * GAS_CONSTANT).abs() < 1e-10);
    }

    #[test]
    fn test_equipartition_energy_scales() {
        let e1 = equipartition_energy(3.0, 300.0);
        let e2 = equipartition_energy(3.0, 600.0);
        assert!((e2 - 2.0 * e1).abs() < 1e-30);
    }

    // --- Partition functions ---

    #[test]
    fn test_partition_translational_positive() {
        // N₂ molecule: mass = 28.014e-3 / 6.022e23 ≈ 4.65e-26 kg
        let m = 28.014e-3 / AVOGADRO;
        let q = partition_translational(m, 0.02241, 300.0).unwrap();
        assert!(q > 1e25); // very large for macroscopic volume
    }

    #[test]
    fn test_partition_rotational() {
        // N₂: I ≈ 1.4e-46 kg·m², σ=2
        let q = partition_rotational(1.4e-46, 300.0, 2).unwrap();
        assert!(q > 10.0 && q < 200.0);
    }

    #[test]
    fn test_partition_vibrational_low_t() {
        // At very low T, q_vib → 1 (ground state only)
        let q = partition_vibrational(1e13, 10.0).unwrap();
        assert!((q - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_partition_vibrational_high_t() {
        // At high T, q_vib → kT/(hν) >> 1
        let q = partition_vibrational(1e12, 3000.0).unwrap();
        assert!(q > 10.0);
    }

    #[test]
    fn test_partition_invalid() {
        assert!(partition_translational(0.0, 1.0, 300.0).is_err());
        assert!(partition_translational(1e-26, 0.0, 300.0).is_err());
        assert!(partition_rotational(1e-46, 300.0, 0).is_err());
        assert!(partition_vibrational(0.0, 300.0).is_err());
    }

    // --- Maxwell-Boltzmann ---

    #[test]
    fn test_mb_speed_ordering() {
        // v_p < v_avg < v_rms
        let m = 28.014e-3 / AVOGADRO;
        let vp = most_probable_speed(m, 300.0).unwrap();
        let va = mean_speed(m, 300.0).unwrap();
        let vr = rms_speed(m, 300.0).unwrap();
        assert!(vp < va);
        assert!(va < vr);
    }

    #[test]
    fn test_n2_rms_speed() {
        // N₂ at 300 K: v_rms ≈ 517 m/s
        let m = 28.014e-3 / AVOGADRO;
        let v = rms_speed(m, 300.0).unwrap();
        assert!((v - 517.0).abs() < 5.0, "N₂ v_rms={v:.0} m/s");
    }

    #[test]
    fn test_mb_pdf_zero_at_zero() {
        let m = 28.014e-3 / AVOGADRO;
        let f = maxwell_boltzmann_speed_pdf(0.0, m, 300.0).unwrap();
        assert!(f.abs() < 1e-30);
    }

    #[test]
    fn test_mb_pdf_positive_at_peak() {
        let m = 28.014e-3 / AVOGADRO;
        let vp = most_probable_speed(m, 300.0).unwrap();
        let f = maxwell_boltzmann_speed_pdf(vp, m, 300.0).unwrap();
        assert!(f > 0.0);
    }

    #[test]
    fn test_mb_speed_invalid() {
        assert!(most_probable_speed(0.0, 300.0).is_err());
        assert!(most_probable_speed(1e-26, 0.0).is_err());
        assert!(maxwell_boltzmann_speed_pdf(-1.0, 1e-26, 300.0).is_err());
    }

    // --- Einstein + Debye ---

    #[test]
    fn test_einstein_high_t_dulong_petit() {
        // T >> Θ_E → Cv → 3R ≈ 24.94 J/(mol·K)
        let cv = einstein_cv(200.0, 5000.0).unwrap();
        assert!((cv - 3.0 * GAS_CONSTANT).abs() < 0.1);
    }

    #[test]
    fn test_einstein_low_t_approaches_zero() {
        let cv = einstein_cv(1000.0, 10.0).unwrap();
        assert!(cv < 0.1, "Einstein Cv={cv} should be near 0 at low T");
    }

    #[test]
    fn test_debye_high_t_dulong_petit() {
        let cv = debye_cv(300.0, 3000.0).unwrap();
        assert!(
            (cv - 3.0 * GAS_CONSTANT).abs() < 0.5,
            "Debye Cv={cv}, expected ~{}",
            3.0 * GAS_CONSTANT
        );
    }

    #[test]
    fn test_debye_low_t_approaches_zero() {
        let cv = debye_cv(300.0, 5.0).unwrap();
        assert!(cv < 1.0, "Debye Cv={cv} should be near 0 at low T");
    }

    #[test]
    fn test_debye_vs_einstein_low_t() {
        // At low T, Debye follows T³ law (more accurate) while Einstein drops exponentially.
        // Debye gives HIGHER Cv than Einstein at low T because Einstein underestimates.
        let cv_e = einstein_cv(300.0, 50.0).unwrap();
        let cv_d = debye_cv(300.0, 50.0).unwrap();
        assert!(
            cv_d > cv_e,
            "Debye {cv_d} should exceed Einstein {cv_e} at low T (T³ vs exponential)"
        );
    }

    #[test]
    fn test_einstein_invalid() {
        assert!(einstein_cv(0.0, 300.0).is_err());
        assert!(einstein_cv(200.0, 0.0).is_err());
    }

    #[test]
    fn test_debye_invalid() {
        assert!(debye_cv(0.0, 300.0).is_err());
        assert!(debye_cv(300.0, 0.0).is_err());
    }
}
