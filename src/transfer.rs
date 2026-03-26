//! Heat transfer — conduction, convection, radiation.
//!
//! All SI units: watts, meters, kelvins, seconds.

use crate::error::{Result, UshmaError};

/// Stefan-Boltzmann constant σ (W/(m²⋅K⁴)).
pub const STEFAN_BOLTZMANN: f64 = 5.670_374_419e-8;

/// Boltzmann constant k_B (J/K).
pub const BOLTZMANN_K: f64 = 1.380_649e-23;

/// Fourier's law: heat flux through conduction.
///
/// q = k⋅A⋅(T_hot - T_cold)/L (watts)
///
/// - `conductivity`: thermal conductivity k (W/(m⋅K))
/// - `area`: cross-sectional area A (m²)
/// - `t_hot`, `t_cold`: temperatures (K)
/// - `thickness`: material thickness L (m)
#[tracing::instrument(level = "debug")]
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
#[tracing::instrument(level = "debug")]
pub fn radiation(emissivity: f64, area: f64, t_surface: f64, t_surrounding: f64) -> Result<f64> {
    if !(0.0..=1.0).contains(&emissivity) {
        return Err(UshmaError::InvalidParameter {
            reason: format!("emissivity {emissivity} must be in [0, 1]"),
        });
    }
    if t_surface < 0.0 {
        return Err(UshmaError::InvalidTemperature { kelvin: t_surface });
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
pub fn thermal_resistance_conduction(conductivity: f64, area: f64, thickness: f64) -> Result<f64> {
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
#[inline]
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
pub fn thermal_diffusivity(conductivity: f64, density: f64, specific_heat: f64) -> Result<f64> {
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

// --- Dimensionless numbers ---

/// Reynolds number: Re = VL/ν (dimensionless).
///
/// Ratio of inertial to viscous forces. Re > ~4000 = turbulent.
/// - `velocity`: flow velocity (m/s)
/// - `length`: characteristic length (m)
/// - `kinematic_viscosity`: ν (m²/s)
pub fn reynolds_number(velocity: f64, length: f64, kinematic_viscosity: f64) -> Result<f64> {
    if kinematic_viscosity <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("kinematic viscosity {kinematic_viscosity} must be positive"),
        });
    }
    Ok(velocity * length / kinematic_viscosity)
}

/// Prandtl number: Pr = ν/α (dimensionless).
///
/// Ratio of momentum diffusivity to thermal diffusivity.
/// - `kinematic_viscosity`: ν (m²/s)
/// - `thermal_diffusivity`: α (m²/s)
pub fn prandtl_number(kinematic_viscosity: f64, thermal_diffusivity: f64) -> Result<f64> {
    if thermal_diffusivity <= 0.0 {
        return Err(UshmaError::DivisionByZero {
            context: "thermal diffusivity must be positive for Prandtl number".into(),
        });
    }
    Ok(kinematic_viscosity / thermal_diffusivity)
}

/// Nusselt number (definition): Nu = hL/k (dimensionless).
///
/// Ratio of convective to conductive heat transfer.
pub fn nusselt_number(h: f64, length: f64, conductivity: f64) -> Result<f64> {
    if conductivity.abs() < 1e-30 {
        return Err(UshmaError::DivisionByZero {
            context: "conductivity cannot be zero for Nusselt number".into(),
        });
    }
    Ok(h * length / conductivity)
}

/// Dittus-Boelter correlation: Nu = 0.023·Re^0.8·Pr^0.4.
///
/// For turbulent flow in smooth circular tubes (Re > 10,000, 0.6 < Pr < 160).
/// Heating case (fluid being heated, n=0.4).
pub fn nusselt_dittus_boelter(re: f64, pr: f64) -> Result<f64> {
    if re <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("Reynolds number {re} must be positive"),
        });
    }
    if pr <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("Prandtl number {pr} must be positive"),
        });
    }
    Ok(0.023 * re.powf(0.8) * pr.powf(0.4))
}

/// Churchill-Chu correlation for natural convection on a vertical plate.
///
/// Nu = {0.825 + 0.387·Ra^(1/6) / [1 + (0.492/Pr)^(9/16)]^(8/27)}²
///
/// Valid for all Ra (laminar and turbulent).
/// - `ra`: Rayleigh number (Ra = Gr·Pr)
/// - `pr`: Prandtl number
pub fn nusselt_natural_vertical(ra: f64, pr: f64) -> Result<f64> {
    if ra < 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("Rayleigh number {ra} must be non-negative"),
        });
    }
    if pr <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("Prandtl number {pr} must be positive"),
        });
    }
    let psi = (1.0 + (0.492 / pr).powf(9.0 / 16.0)).powf(-8.0 / 27.0);
    let nu = (0.825 + 0.387 * ra.powf(1.0 / 6.0) * psi).powi(2);
    Ok(nu)
}

// --- Boundary layer ---

/// Laminar boundary layer thickness: δ = 5x/√Re_x (Blasius solution, m).
///
/// Valid for laminar flow (Re_x < ~5×10⁵).
pub fn boundary_layer_thickness(x: f64, re_x: f64) -> Result<f64> {
    if x <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("distance x={x} m must be positive"),
        });
    }
    if re_x <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("local Reynolds number {re_x} must be positive"),
        });
    }
    Ok(5.0 * x / re_x.sqrt())
}

/// Thermal boundary layer thickness: δ_t = δ/Pr^(1/3) (m).
///
/// Ratio of velocity to thermal boundary layer.
pub fn thermal_boundary_layer(x: f64, re_x: f64, pr: f64) -> Result<f64> {
    if pr <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("Prandtl number {pr} must be positive"),
        });
    }
    let delta = boundary_layer_thickness(x, re_x)?;
    Ok(delta / pr.powf(1.0 / 3.0))
}

// --- Fin heat transfer ---

/// Fin parameter: m = √(hP/(kA_c)) (1/m).
///
/// - `h`: convective coefficient (W/(m²·K))
/// - `perimeter`: fin cross-section perimeter (m)
/// - `k`: fin thermal conductivity (W/(m·K))
/// - `cross_area`: fin cross-section area (m²)
pub fn fin_parameter(h: f64, perimeter: f64, k: f64, cross_area: f64) -> Result<f64> {
    let denom = k * cross_area;
    if denom <= 0.0 {
        return Err(UshmaError::DivisionByZero {
            context: "k·A_c must be positive for fin parameter".into(),
        });
    }
    if h < 0.0 || perimeter <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("h={h}, perimeter={perimeter} must be positive"),
        });
    }
    Ok((h * perimeter / denom).sqrt())
}

/// Heat transfer from a rectangular fin with insulated tip (W).
///
/// Q = √(hPkA_c)·(T_base - T_fluid)·tanh(mL)
#[tracing::instrument(level = "debug")]
pub fn fin_rectangular_heat(
    h: f64,
    perimeter: f64,
    k: f64,
    cross_area: f64,
    length: f64,
    t_base: f64,
    t_fluid: f64,
) -> Result<f64> {
    let m = fin_parameter(h, perimeter, k, cross_area)?;
    if length <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("fin length {length} m must be positive"),
        });
    }
    let coeff = (h * perimeter * k * cross_area).sqrt();
    Ok(coeff * (t_base - t_fluid) * (m * length).tanh())
}

/// Efficiency of a rectangular fin with insulated tip: η = tanh(mL)/(mL).
///
/// η = 1.0 for zero length, approaches 0 for very long fins.
pub fn fin_efficiency_rectangular(m: f64, length: f64) -> Result<f64> {
    let ml = m * length;
    if ml.abs() < 1e-30 {
        return Ok(1.0); // degenerate: zero-length fin is 100% efficient
    }
    Ok(ml.tanh() / ml)
}

/// Fin effectiveness: ε = Q_fin / Q_no_fin (dimensionless).
///
/// Ratio of fin heat transfer to what the base area would transfer without the fin.
/// ε > 1 means the fin improves heat transfer (which it should).
pub fn fin_effectiveness(
    q_fin: f64,
    h: f64,
    cross_area: f64,
    t_base: f64,
    t_fluid: f64,
) -> Result<f64> {
    let q_no_fin = h * cross_area * (t_base - t_fluid);
    if q_no_fin.abs() < 1e-30 {
        return Err(UshmaError::DivisionByZero {
            context: "no temperature difference for fin effectiveness".into(),
        });
    }
    Ok(q_fin / q_no_fin)
}

// --- Heat exchangers ---

/// Log mean temperature difference for parallel-flow heat exchanger (K).
///
/// LMTD = (ΔT₁ - ΔT₂) / ln(ΔT₁/ΔT₂)
/// where ΔT₁ = T_h_in - T_c_in, ΔT₂ = T_h_out - T_c_out.
#[tracing::instrument(level = "debug")]
pub fn lmtd_parallel(t_h_in: f64, t_h_out: f64, t_c_in: f64, t_c_out: f64) -> Result<f64> {
    let dt1 = t_h_in - t_c_in;
    let dt2 = t_h_out - t_c_out;
    lmtd_from_deltas(dt1, dt2)
}

/// Log mean temperature difference for counter-flow heat exchanger (K).
///
/// LMTD = (ΔT₁ - ΔT₂) / ln(ΔT₁/ΔT₂)
/// where ΔT₁ = T_h_in - T_c_out, ΔT₂ = T_h_out - T_c_in.
#[tracing::instrument(level = "debug")]
pub fn lmtd_counter(t_h_in: f64, t_h_out: f64, t_c_in: f64, t_c_out: f64) -> Result<f64> {
    let dt1 = t_h_in - t_c_out;
    let dt2 = t_h_out - t_c_in;
    lmtd_from_deltas(dt1, dt2)
}

fn lmtd_from_deltas(dt1: f64, dt2: f64) -> Result<f64> {
    if dt1 <= 0.0 || dt2 <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("temperature differences must be positive: ΔT₁={dt1}, ΔT₂={dt2}"),
        });
    }
    let ratio = dt1 / dt2;
    if (ratio - 1.0).abs() < 1e-10 {
        // ΔT₁ ≈ ΔT₂: use arithmetic mean to avoid 0/0
        return Ok(0.5 * (dt1 + dt2));
    }
    Ok((dt1 - dt2) / ratio.ln())
}

/// Heat transfer rate using LMTD method: Q = U·A·LMTD (W).
///
/// - `u`: overall heat transfer coefficient (W/(m²·K))
/// - `area`: heat exchanger area (m²)
/// - `lmtd`: log mean temperature difference (K)
#[inline]
#[must_use]
pub fn heat_exchanger_lmtd(u: f64, area: f64, lmtd: f64) -> f64 {
    u * area * lmtd
}

/// Number of transfer units: NTU = UA/C_min (dimensionless).
pub fn ntu(u: f64, area: f64, c_min: f64) -> Result<f64> {
    if c_min <= 0.0 {
        return Err(UshmaError::DivisionByZero {
            context: "C_min must be positive for NTU".into(),
        });
    }
    Ok(u * area / c_min)
}

/// Effectiveness of a parallel-flow heat exchanger.
///
/// ε = (1 - exp(-NTU·(1 + Cr))) / (1 + Cr)
/// Special case Cr = 0: ε = 1 - exp(-NTU).
pub fn effectiveness_parallel(ntu_val: f64, c_ratio: f64) -> Result<f64> {
    if ntu_val < 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("NTU {ntu_val} must be non-negative"),
        });
    }
    if !(0.0..=1.0).contains(&c_ratio) {
        return Err(UshmaError::InvalidParameter {
            reason: format!("capacity ratio Cr={c_ratio} must be in [0, 1]"),
        });
    }
    if c_ratio < 1e-30 {
        return Ok(1.0 - (-ntu_val).exp());
    }
    Ok((1.0 - (-(1.0 + c_ratio) * ntu_val).exp()) / (1.0 + c_ratio))
}

/// Effectiveness of a counter-flow heat exchanger.
///
/// ε = (1 - exp(-NTU·(1 - Cr))) / (1 - Cr·exp(-NTU·(1 - Cr)))
/// Special case Cr = 0: ε = 1 - exp(-NTU).
/// Special case Cr = 1: ε = NTU / (1 + NTU).
pub fn effectiveness_counter(ntu_val: f64, c_ratio: f64) -> Result<f64> {
    if ntu_val < 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("NTU {ntu_val} must be non-negative"),
        });
    }
    if !(0.0..=1.0).contains(&c_ratio) {
        return Err(UshmaError::InvalidParameter {
            reason: format!("capacity ratio Cr={c_ratio} must be in [0, 1]"),
        });
    }
    if c_ratio < 1e-30 {
        return Ok(1.0 - (-ntu_val).exp());
    }
    if (c_ratio - 1.0).abs() < 1e-10 {
        return Ok(ntu_val / (1.0 + ntu_val));
    }
    let exp_term = (-(1.0 - c_ratio) * ntu_val).exp();
    Ok((1.0 - exp_term) / (1.0 - c_ratio * exp_term))
}

/// Heat transfer rate using ε-NTU method: Q = ε·C_min·(T_h_in - T_c_in) (W).
#[inline]
#[must_use]
pub fn heat_exchanger_ntu(effectiveness: f64, c_min: f64, t_h_in: f64, t_c_in: f64) -> f64 {
    effectiveness * c_min * (t_h_in - t_c_in)
}

// --- Radiation view factors ---

/// View factor between two identical, directly opposed, parallel rectangles.
///
/// F₁₂ for two aligned parallel rectangles of width W and height H separated by distance D.
/// Uses the Hottel crossed-string formula.
#[tracing::instrument(level = "debug")]
pub fn view_factor_parallel_plates(width: f64, height: f64, distance: f64) -> Result<f64> {
    if width <= 0.0 || height <= 0.0 || distance <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("dimensions must be positive: W={width}, H={height}, D={distance}"),
        });
    }
    let x = width / distance;
    let y = height / distance;
    let x2 = x * x;
    let y2 = y * y;

    let term1 = ((1.0 + x2) * (1.0 + y2) / (1.0 + x2 + y2)).sqrt().ln();
    let term2 = x * (1.0 + y2).sqrt() * (x / (1.0 + y2).sqrt()).atan();
    let term3 = y * (1.0 + x2).sqrt() * (y / (1.0 + x2).sqrt()).atan();
    let term4 = x * x.atan();
    let term5 = y * y.atan();

    let f12 = (2.0 / (std::f64::consts::PI * x * y)) * (term1 + term2 + term3 - term4 - term5);
    Ok(f12)
}

/// View factor between two perpendicular rectangles sharing a common edge.
///
/// F₁₂ for rectangle 1 (W × H₁) perpendicular to rectangle 2 (W × H₂)
/// sharing edge of length W.
pub fn view_factor_perpendicular_plates(width: f64, h1: f64, h2: f64) -> Result<f64> {
    if width <= 0.0 || h1 <= 0.0 || h2 <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("dimensions must be positive: W={width}, H1={h1}, H2={h2}"),
        });
    }
    let h = h1 / width;
    let w = h2 / width;
    let h2v = h * h;
    let w2 = w * w;

    let a = (1.0 + h2v) * (1.0 + w2);
    let b = 1.0 + h2v + w2;

    let term1 = w * (1.0 + h2v).sqrt().atan() + h * (1.0 + w2).sqrt().atan();

    // Simplified Hottel formula
    let f12 = (1.0 / (std::f64::consts::PI * h)) * (term1 - b.sqrt().atan() + 0.25 * (a / b).ln());
    Ok(f12.clamp(0.0, 1.0))
}

/// View factor between two coaxial parallel disks.
///
/// F₁₂ for disk 1 (radius r₁) seeing disk 2 (radius r₂) at distance d.
pub fn view_factor_coaxial_disks(r1: f64, r2: f64, distance: f64) -> Result<f64> {
    if r1 <= 0.0 || r2 <= 0.0 || distance <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("dimensions must be positive: r1={r1}, r2={r2}, d={distance}"),
        });
    }
    let r1_d = r1 / distance;
    let r2_d = r2 / distance;
    let s = 1.0 + (1.0 + r2_d * r2_d) / (r1_d * r1_d);
    let f12 = 0.5 * (s - (s * s - 4.0 * (r2_d / r1_d).powi(2)).sqrt());
    Ok(f12)
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

    #[test]
    fn test_radiation_invalid_emissivity() {
        assert!(radiation(-0.1, 1.0, 373.0, 293.0).is_err());
        assert!(radiation(1.1, 1.0, 373.0, 293.0).is_err());
        // Boundary values should work
        assert!(radiation(0.0, 1.0, 373.0, 293.0).is_ok());
        assert!(radiation(1.0, 1.0, 373.0, 293.0).is_ok());
    }

    #[test]
    fn test_conduction_zero_thickness() {
        assert!(conduction(401.0, 1.0, 373.0, 293.0, 0.0).is_err());
    }

    #[test]
    fn test_thermal_resistance_convection_zero() {
        assert!(thermal_resistance_convection(0.0, 1.0).is_err());
        assert!(thermal_resistance_convection(25.0, 0.0).is_err());
    }

    #[test]
    fn test_thermal_diffusivity_zero_density() {
        assert!(thermal_diffusivity(401.0, 0.0, 385.0).is_err());
    }

    #[test]
    fn test_biot_number_zero_conductivity() {
        assert!(biot_number(25.0, 0.01, 0.0).is_err());
    }

    #[test]
    fn test_conduction_nan_propagation() {
        let result = conduction(f64::NAN, 1.0, 373.0, 293.0, 0.1);
        // NaN conductivity should be caught by validation (NaN < 0.0 is false)
        // so it passes through — result contains NaN
        assert!(result.unwrap().is_nan());
    }

    #[test]
    fn test_radiation_nan_emissivity() {
        // NaN is not in [0, 1] so should be rejected
        assert!(radiation(f64::NAN, 1.0, 373.0, 293.0).is_err());
    }

    #[test]
    fn test_heat_stored_negative_delta() {
        // Cooling: negative ΔT → negative heat stored (heat lost)
        let q = heat_stored(1.0, 4186.0, -10.0);
        assert!(q < 0.0);
    }

    #[test]
    fn test_thermal_resistance_parallel_single() {
        let r = thermal_resistance_parallel(&[5.0]).unwrap();
        assert!((r - 5.0).abs() < 1e-10);
    }

    // --- Dimensionless number tests ---

    #[test]
    fn test_reynolds_number() {
        // Water at 1 m/s, 0.1 m pipe, ν = 1e-6 m²/s → Re = 100,000
        let re = reynolds_number(1.0, 0.1, 1e-6).unwrap();
        assert!((re - 100_000.0).abs() < 1.0);
    }

    #[test]
    fn test_reynolds_invalid() {
        assert!(reynolds_number(1.0, 0.1, 0.0).is_err());
        assert!(reynolds_number(1.0, 0.1, -1e-6).is_err());
    }

    #[test]
    fn test_prandtl_number() {
        // Air at 300K: ν ≈ 1.6e-5, α ≈ 2.2e-5 → Pr ≈ 0.73
        let pr = prandtl_number(1.6e-5, 2.2e-5).unwrap();
        assert!((pr - 0.727).abs() < 0.01);
    }

    #[test]
    fn test_nusselt_number() {
        let nu = nusselt_number(25.0, 0.1, 0.6).unwrap();
        assert!((nu - 25.0 * 0.1 / 0.6).abs() < 1e-10);
    }

    #[test]
    fn test_dittus_boelter() {
        // Re=50000, Pr=0.7 → Nu = 0.023 * 50000^0.8 * 0.7^0.4 ≈ 119.7
        let nu = nusselt_dittus_boelter(50_000.0, 0.7).unwrap();
        assert!(nu > 100.0 && nu < 150.0);
    }

    #[test]
    fn test_dittus_boelter_invalid() {
        assert!(nusselt_dittus_boelter(0.0, 0.7).is_err());
        assert!(nusselt_dittus_boelter(50_000.0, 0.0).is_err());
    }

    #[test]
    fn test_nusselt_natural_vertical() {
        // Ra=1e9, Pr=0.7 (air, turbulent natural convection)
        let nu = nusselt_natural_vertical(1e9, 0.7).unwrap();
        assert!(nu > 50.0); // should be significant
    }

    #[test]
    fn test_nusselt_natural_vertical_zero_ra() {
        // Ra=0 → pure conduction, Nu should be small
        let nu = nusselt_natural_vertical(0.0, 0.7).unwrap();
        assert!((nu - 0.825_f64.powi(2)).abs() < 0.01);
    }

    // --- Boundary layer tests ---

    #[test]
    fn test_boundary_layer_thickness() {
        // x=0.5 m, Re_x=100,000 → δ = 5*0.5/√100000 ≈ 0.0079 m
        let delta = boundary_layer_thickness(0.5, 100_000.0).unwrap();
        assert!((delta - 0.00791).abs() < 0.001);
    }

    #[test]
    fn test_thermal_boundary_layer() {
        // Pr=0.7 (air) → δ_t > δ (thermal BL thicker for Pr < 1)
        let delta = boundary_layer_thickness(0.5, 100_000.0).unwrap();
        let delta_t = thermal_boundary_layer(0.5, 100_000.0, 0.7).unwrap();
        assert!(delta_t > delta);
    }

    #[test]
    fn test_boundary_layer_invalid() {
        assert!(boundary_layer_thickness(0.0, 100_000.0).is_err());
        assert!(boundary_layer_thickness(0.5, 0.0).is_err());
        assert!(thermal_boundary_layer(0.5, 100_000.0, 0.0).is_err());
    }

    // --- Fin tests ---

    #[test]
    fn test_fin_parameter() {
        // Aluminum fin: h=25, P=0.02 m (thin rectangular), k=237, A_c=1e-4 m²
        let m = fin_parameter(25.0, 0.02, 237.0, 1e-4).unwrap();
        // m = √(25*0.02 / (237*1e-4)) = √(0.5/0.0237) = √21.1 ≈ 4.59
        assert!((m - 4.59).abs() < 0.1);
    }

    #[test]
    fn test_fin_rectangular_heat() {
        let q = fin_rectangular_heat(25.0, 0.02, 237.0, 1e-4, 0.05, 373.15, 293.15).unwrap();
        assert!(q > 0.0); // heat flows from hot base to cooler fluid
    }

    #[test]
    fn test_fin_efficiency_short_fin() {
        // Very short fin → η ≈ 1.0
        let eta = fin_efficiency_rectangular(5.0, 0.001).unwrap();
        assert!((eta - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_fin_efficiency_long_fin() {
        // Long fin → η < 1.0
        let eta = fin_efficiency_rectangular(5.0, 0.5).unwrap();
        assert!(eta < 1.0);
        assert!(eta > 0.0);
    }

    #[test]
    fn test_fin_efficiency_zero_length() {
        let eta = fin_efficiency_rectangular(5.0, 0.0).unwrap();
        assert!((eta - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_fin_effectiveness_greater_than_one() {
        let q_fin = fin_rectangular_heat(25.0, 0.02, 237.0, 1e-4, 0.05, 373.15, 293.15).unwrap();
        let eff = fin_effectiveness(q_fin, 25.0, 1e-4, 373.15, 293.15).unwrap();
        assert!(eff > 1.0, "Fin effectiveness {eff} should exceed 1.0");
    }

    #[test]
    fn test_fin_parameter_invalid() {
        assert!(fin_parameter(25.0, 0.02, 0.0, 1e-4).is_err()); // k=0
        assert!(fin_parameter(25.0, 0.0, 237.0, 1e-4).is_err()); // P=0
    }

    // --- LMTD tests ---

    #[test]
    fn test_lmtd_parallel() {
        // Hot: 150→100, Cold: 30→80 → ΔT1=120, ΔT2=20
        let lmtd = lmtd_parallel(150.0, 100.0, 30.0, 80.0).unwrap();
        let expected = (120.0 - 20.0) / (120.0_f64 / 20.0).ln();
        assert!((lmtd - expected).abs() < 0.1);
    }

    #[test]
    fn test_lmtd_counter() {
        // Hot: 150→100, Cold: 30→80 → ΔT1=150-80=70, ΔT2=100-30=70
        // Equal ΔTs → LMTD = arithmetic mean = 70
        let lmtd = lmtd_counter(150.0, 100.0, 30.0, 80.0).unwrap();
        assert!((lmtd - 70.0).abs() < 0.1);
    }

    #[test]
    fn test_lmtd_counter_gte_parallel() {
        let lmtd_p = lmtd_parallel(200.0, 120.0, 50.0, 90.0).unwrap();
        let lmtd_c = lmtd_counter(200.0, 120.0, 50.0, 90.0).unwrap();
        assert!(lmtd_c >= lmtd_p);
    }

    #[test]
    fn test_lmtd_invalid() {
        // Cold side hotter than hot side at inlet → negative ΔT
        assert!(lmtd_parallel(50.0, 40.0, 60.0, 30.0).is_err());
    }

    #[test]
    fn test_heat_exchanger_lmtd_basic() {
        let q = heat_exchanger_lmtd(500.0, 2.0, 50.0);
        assert!((q - 50_000.0).abs() < 1e-10);
    }

    // --- ε-NTU tests ---

    #[test]
    fn test_ntu_basic() {
        let n = ntu(500.0, 2.0, 1000.0).unwrap();
        assert!((n - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_effectiveness_parallel_increases_with_ntu() {
        let e1 = effectiveness_parallel(1.0, 0.5).unwrap();
        let e2 = effectiveness_parallel(3.0, 0.5).unwrap();
        assert!(e2 > e1);
        assert!(e2 <= 1.0);
    }

    #[test]
    fn test_effectiveness_counter_gte_parallel() {
        let ep = effectiveness_parallel(2.0, 0.5).unwrap();
        let ec = effectiveness_counter(2.0, 0.5).unwrap();
        assert!(ec >= ep);
    }

    #[test]
    fn test_effectiveness_counter_balanced() {
        // Cr=1: ε = NTU/(1+NTU)
        let e = effectiveness_counter(2.0, 1.0).unwrap();
        assert!((e - 2.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_effectiveness_cr_zero() {
        // Cr=0 (condenser): ε = 1 - exp(-NTU)
        let ep = effectiveness_parallel(2.0, 0.0).unwrap();
        let ec = effectiveness_counter(2.0, 0.0).unwrap();
        let expected = 1.0 - (-2.0_f64).exp();
        assert!((ep - expected).abs() < 1e-10);
        assert!((ec - expected).abs() < 1e-10);
    }

    #[test]
    fn test_effectiveness_invalid() {
        assert!(effectiveness_parallel(-1.0, 0.5).is_err());
        assert!(effectiveness_parallel(2.0, 1.5).is_err());
        assert!(effectiveness_counter(2.0, -0.1).is_err());
    }

    #[test]
    fn test_heat_exchanger_ntu_basic() {
        let q = heat_exchanger_ntu(0.8, 1000.0, 200.0, 50.0);
        assert!((q - 0.8 * 1000.0 * 150.0).abs() < 1e-10);
    }

    // --- View factor tests ---

    #[test]
    fn test_view_factor_parallel_plates_close() {
        // Large plates close together → F approaches 1
        let f = view_factor_parallel_plates(10.0, 10.0, 0.1).unwrap();
        assert!(f > 0.9, "F={f} should be close to 1 for large close plates");
    }

    #[test]
    fn test_view_factor_parallel_plates_far() {
        // Small plates far apart → F approaches 0
        let f = view_factor_parallel_plates(0.1, 0.1, 10.0).unwrap();
        assert!(f < 0.01);
    }

    #[test]
    fn test_view_factor_parallel_plates_range() {
        let f = view_factor_parallel_plates(1.0, 1.0, 1.0).unwrap();
        assert!(f > 0.0 && f < 1.0);
    }

    #[test]
    fn test_view_factor_parallel_plates_invalid() {
        assert!(view_factor_parallel_plates(0.0, 1.0, 1.0).is_err());
        assert!(view_factor_parallel_plates(1.0, 0.0, 1.0).is_err());
        assert!(view_factor_parallel_plates(1.0, 1.0, 0.0).is_err());
    }

    #[test]
    fn test_view_factor_perpendicular_range() {
        let f = view_factor_perpendicular_plates(1.0, 1.0, 1.0).unwrap();
        assert!(f > 0.0 && f < 1.0);
    }

    #[test]
    fn test_view_factor_perpendicular_invalid() {
        assert!(view_factor_perpendicular_plates(0.0, 1.0, 1.0).is_err());
    }

    #[test]
    fn test_view_factor_coaxial_disks_equal() {
        // Two equal disks close together → F approaches 1
        let f = view_factor_coaxial_disks(1.0, 1.0, 0.1).unwrap();
        assert!(f > 0.9);
    }

    #[test]
    fn test_view_factor_coaxial_disks_far() {
        // Far apart → F approaches 0
        let f = view_factor_coaxial_disks(0.1, 0.1, 10.0).unwrap();
        assert!(f < 0.01);
    }

    #[test]
    fn test_view_factor_coaxial_disks_range() {
        let f = view_factor_coaxial_disks(0.5, 0.5, 1.0).unwrap();
        assert!(f > 0.0 && f < 1.0);
    }

    #[test]
    fn test_view_factor_invalid() {
        assert!(view_factor_coaxial_disks(0.0, 0.5, 1.0).is_err());
        assert!(view_factor_coaxial_disks(0.5, 0.0, 1.0).is_err());
        assert!(view_factor_coaxial_disks(0.5, 0.5, 0.0).is_err());
    }
}
