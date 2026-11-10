//! Statistical thermodynamics — Boltzmann, partition functions, speed distributions, solid Cv.
//!
//! Bridges macroscopic thermodynamics to molecular physics.
//! All SI units: joules, kelvins, kilograms, meters, seconds.

use crate::error::{Result, UshmaError};
use crate::state::GAS_CONSTANT;
use crate::transfer::BOLTZMANN_K;
use crate::validate::{finite_ge, finite_gt, invalid_param};

/// Planck constant h (J·s).
pub const PLANCK: f64 = 6.626_070_15e-34;

/// Avogadro number (1/mol).
pub const AVOGADRO: f64 = 6.022_140_76e23;

// --- Boltzmann distribution ---

/// Boltzmann factor: exp(-E/(kT)).
///
/// Proportional to the probability of a state with energy E at temperature T.
#[inline]
pub fn boltzmann_probability(energy: f64, temperature: f64) -> Result<f64> {
    if !finite_gt(temperature, 0.0) {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    Ok((-energy / (BOLTZMANN_K * temperature)).exp())
}

/// Population ratio of two energy levels: N₂/N₁ = exp(-(E₂-E₁)/(kT)).
#[inline]
pub fn boltzmann_ratio(e1: f64, e2: f64, temperature: f64) -> Result<f64> {
    if !finite_gt(temperature, 0.0) {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    Ok((-(e2 - e1) / (BOLTZMANN_K * temperature)).exp())
}

/// Entropy from microstates: S = k·ln(W) (J/K).
#[inline]
pub fn entropy_from_microstates(num_microstates: f64) -> Result<f64> {
    if !finite_gt(num_microstates, 0.0) {
        return Err(invalid_param(format_args!(
            "microstates {num_microstates} must be positive"
        )));
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
#[inline]
pub fn partition_translational(mass: f64, volume: f64, temperature: f64) -> Result<f64> {
    if !finite_gt(mass, 0.0) {
        return Err(invalid_param(format_args!(
            "mass {mass} kg must be positive"
        )));
    }
    if !finite_gt(volume, 0.0) {
        return Err(UshmaError::InvalidVolume {
            cubic_meters: volume,
        });
    }
    if !finite_gt(temperature, 0.0) {
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
#[inline]
pub fn partition_rotational(
    moment_of_inertia: f64,
    temperature: f64,
    symmetry_number: u32,
) -> Result<f64> {
    if !finite_gt(moment_of_inertia, 0.0) {
        return Err(invalid_param(format_args!(
            "moment of inertia {moment_of_inertia} must be positive"
        )));
    }
    if !finite_gt(temperature, 0.0) {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    if symmetry_number == 0 {
        return Err(invalid_param(format_args!("symmetry number must be >= 1")));
    }
    let theta_rot = PLANCK * PLANCK
        / (8.0 * std::f64::consts::PI * std::f64::consts::PI * moment_of_inertia * BOLTZMANN_K);
    Ok(temperature / (symmetry_number as f64 * theta_rot))
}

/// Vibrational partition function: q_vib = 1/(1 - exp(-hν/(kT))).
///
/// Energies are measured from the vibrational ground state (zero-point energy
/// excluded). Computed with `exp_m1`, so it stays accurate when hν ≪ kT
/// (q_vib → kT/(hν)).
///
/// - `frequency`: vibrational frequency ν (Hz)
/// - `temperature`: T (K)
#[inline]
pub fn partition_vibrational(frequency: f64, temperature: f64) -> Result<f64> {
    if !finite_gt(frequency, 0.0) {
        return Err(invalid_param(format_args!(
            "frequency {frequency} Hz must be positive"
        )));
    }
    if !finite_gt(temperature, 0.0) {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    let x = PLANCK * frequency / (BOLTZMANN_K * temperature);
    // 1 − e^(−x) without cancellation for small x.
    let denom = -(-x).exp_m1();
    if !finite_gt(denom, 0.0) {
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
#[inline]
pub fn maxwell_boltzmann_speed_pdf(speed: f64, mass: f64, temperature: f64) -> Result<f64> {
    if !finite_ge(speed, 0.0) {
        return Err(invalid_param(format_args!(
            "speed {speed} m/s must be non-negative"
        )));
    }
    if !finite_gt(mass, 0.0) {
        return Err(invalid_param(format_args!(
            "mass {mass} kg must be positive"
        )));
    }
    if !finite_gt(temperature, 0.0) {
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
#[inline]
pub fn most_probable_speed(mass: f64, temperature: f64) -> Result<f64> {
    validate_mass_temp(mass, temperature)?;
    Ok((2.0 * BOLTZMANN_K * temperature / mass).sqrt())
}

/// Mean speed: v_avg = √(8kT/(πm)) (m/s).
#[inline]
pub fn mean_speed(mass: f64, temperature: f64) -> Result<f64> {
    validate_mass_temp(mass, temperature)?;
    Ok((8.0 * BOLTZMANN_K * temperature / (std::f64::consts::PI * mass)).sqrt())
}

/// Root-mean-square speed: v_rms = √(3kT/m) (m/s).
#[inline]
pub fn rms_speed(mass: f64, temperature: f64) -> Result<f64> {
    validate_mass_temp(mass, temperature)?;
    Ok((3.0 * BOLTZMANN_K * temperature / mass).sqrt())
}

fn validate_mass_temp(mass: f64, temperature: f64) -> Result<()> {
    if !finite_gt(mass, 0.0) {
        return Err(invalid_param(format_args!(
            "mass {mass} kg must be positive"
        )));
    }
    if !finite_gt(temperature, 0.0) {
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
#[inline]
pub fn einstein_cv(einstein_temp: f64, temperature: f64) -> Result<f64> {
    if !finite_gt(einstein_temp, 0.0) {
        return Err(invalid_param(format_args!(
            "Einstein temperature {einstein_temp} K must be positive"
        )));
    }
    if !finite_gt(temperature, 0.0) {
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

/// Θ_D/T above which the Debye integral is replaced by its value at ∞
/// (the tail beyond x = 50 is below 1e-15 of the total).
const DEBYE_X_MAX: f64 = 50.0;

/// Debye model for solid heat capacity (J/(mol·K)).
///
/// Cv = 9R·(T/Θ_D)³·∫₀^(Θ_D/T) x⁴eˣ/(eˣ-1)² dx
///
/// The integrand is evaluated as x⁴e⁻ˣ/expm1(−x)², which cannot overflow, and
/// integrated with [`hisab::calc::integral_simpson`] (200 intervals, relative
/// error ≈ 1e-10). For Θ_D/T > 50 the exact low-temperature T³ law
/// Cv = (12π⁴/5)·R·(T/Θ_D)³ is used.
pub fn debye_cv(debye_temp: f64, temperature: f64) -> Result<f64> {
    if !finite_gt(debye_temp, 0.0) {
        return Err(invalid_param(format_args!(
            "Debye temperature {debye_temp} K must be positive"
        )));
    }
    if !finite_gt(temperature, 0.0) {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    let td_t = debye_temp / temperature;

    // High-T limit: Cv → 3R
    if td_t < 0.01 {
        return Ok(3.0 * GAS_CONSTANT);
    }

    let ratio = temperature / debye_temp;
    let ratio3 = ratio * ratio * ratio;

    // Low-T limit: ∫₀^∞ x⁴eˣ/(eˣ−1)² dx = 4π⁴/15
    if td_t > DEBYE_X_MAX {
        let pi2 = std::f64::consts::PI * std::f64::consts::PI;
        return Ok(9.0 * GAS_CONSTANT * ratio3 * 4.0 * pi2 * pi2 / 15.0);
    }

    let integrand = |x: f64| -> f64 {
        if x < 1e-8 {
            return x * x; // x⁴eˣ/(eˣ−1)² → x² as x → 0
        }
        // e^(−x) = 1 + expm1(−x): one transcendental per point.
        let em1 = (-x).exp_m1();
        x * x * x * x * (1.0 + em1) / (em1 * em1)
    };

    let integral = hisab::calc::integral_simpson(integrand, 0.0, td_t, 200)
        .map_err(|e| invalid_param(format_args!("Debye integral failed: {e}")))?;

    Ok(9.0 * GAS_CONSTANT * ratio3 * integral)
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

    // --- Regression tests from the 1.4.0 audit ---

    #[test]
    fn test_debye_cv_low_temperature_t3_law() {
        // exp(Θ/T) used to overflow (inf/inf = NaN) for Θ/T > 709.
        let pi2 = std::f64::consts::PI * std::f64::consts::PI;
        for (theta, t) in [(2230.0, 3.0), (300.0, 0.42), (428.0, 1.0)] {
            let cv = debye_cv(theta, t).unwrap();
            let r: f64 = t / theta;
            let t3 = 12.0 * pi2 * pi2 / 5.0 * GAS_CONSTANT * r.powi(3);
            assert!(cv.is_finite());
            assert!(
                (cv / t3 - 1.0).abs() < 1e-12,
                "Θ={theta}, T={t}: {cv} vs {t3}"
            );
        }
    }

    #[test]
    fn test_debye_cv_reference_values() {
        // scipy.integrate.quad reference values of 9R(T/Θ)³∫₀^(Θ/T) …
        for (xd, cv_ref) in [
            (1.0, 23.739_423_794_400_19),
            (5.0, 9.195_001_381_844_16),
            (10.0, 1.891_232_686_118_677),
            (30.0, 0.071_991_488_298_527),
        ] {
            let cv = debye_cv(xd * 100.0, 100.0).unwrap();
            assert!((cv / cv_ref - 1.0).abs() < 1e-9, "Θ/T={xd}: {cv}");
        }
    }

    #[test]
    fn test_debye_cv_continuous_at_cutoff() {
        let below = debye_cv(49.999_999, 1.0).unwrap();
        let above = debye_cv(50.000_001, 1.0).unwrap();
        assert!((below / above - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_partition_vibrational_small_x() {
        // hν ≪ kT: q_vib → kT/(hν). The 1 − e^(−x) form lost all precision here.
        let t = 1e6;
        let nu = 1.0;
        let x = PLANCK * nu / (BOLTZMANN_K * t);
        let q = partition_vibrational(nu, t).unwrap();
        assert!((q * x - 1.0).abs() < 1e-9, "{q}");
    }

    #[test]
    fn test_stat_rejects_nan() {
        assert!(boltzmann_probability(1e-21, f64::NAN).is_err());
        assert!(debye_cv(f64::NAN, 300.0).is_err());
        assert!(einstein_cv(300.0, f64::NAN).is_err());
        assert!(partition_vibrational(f64::NAN, 300.0).is_err());
        assert!(maxwell_boltzmann_speed_pdf(f64::NAN, 4.65e-26, 300.0).is_err());
    }
}
