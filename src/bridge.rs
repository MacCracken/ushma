//! Cross-crate bridges — convert primitive values from other AGNOS science crates
//! into ushma thermodynamics parameters and vice versa.
//!
//! Always available — takes primitive values (f64), no science crate deps.
//!
//! # Architecture
//!
//! ```text
//! bijli   (electromagnetism)  ──┐
//! kimiya  (chemistry)           ┤
//! dravya  (material science)    ┼──> bridge ──> ushma thermo parameters
//! badal   (weather/atmosphere)  ┤
//! pravash (fluid dynamics)     ┘
//! ```

// ── Bijli bridges (electromagnetism) ───────────────────────────────────────

/// Convert Joule heating power (W) to temperature rise rate (K/s).
///
/// dT/dt = P / (m × c_p)
#[must_use]
#[inline]
pub fn joule_heating_to_temperature_rate(
    power_w: f64,
    mass_kg: f64,
    specific_heat_j_per_kg_k: f64,
) -> f64 {
    let denom = mass_kg * specific_heat_j_per_kg_k;
    if denom <= 0.0 {
        return 0.0;
    }
    power_w / denom
}

/// Convert EM absorption rate (W/m³) to volumetric heat source for
/// thermal simulation.
///
/// Returns the same value (units are already compatible), but validates
/// and clamps to non-negative.
#[must_use]
#[inline]
pub fn em_absorption_to_heat_source(absorption_w_m3: f64) -> f64 {
    absorption_w_m3.max(0.0)
}

// ── Kimiya bridges (chemistry) ─────────────────────────────────────────────

/// Convert reaction enthalpy (J/mol) and molar rate (mol/s) to heat
/// release rate (W).
///
/// Q̇ = ΔH × ṅ. Exothermic reactions have negative ΔH, so we negate
/// to get positive heat release.
#[must_use]
#[inline]
pub fn reaction_heat_rate(enthalpy_j_per_mol: f64, molar_rate_mol_per_s: f64) -> f64 {
    -enthalpy_j_per_mol * molar_rate_mol_per_s
}

/// Convert equilibrium constant K and standard Gibbs energy (J/mol)
/// to the temperature at which the reaction achieves that K.
///
/// From Van't Hoff: T = -ΔG° / (R × ln(K))
/// Returns temperature in Kelvin, or 0 if inputs are invalid.
#[must_use]
#[inline]
pub fn equilibrium_temperature(gibbs_energy_j_per_mol: f64, equilibrium_constant: f64) -> f64 {
    const R: f64 = 8.314; // J/(mol·K)
    if equilibrium_constant <= 0.0 {
        return 0.0;
    }
    let ln_k = equilibrium_constant.ln();
    if ln_k.abs() < 1e-15 {
        return 0.0;
    }
    let t = -gibbs_energy_j_per_mol / (R * ln_k);
    t.max(0.0)
}

// ── Dravya bridges (material science) ──────────────────────────────────────

/// Convert temperature (K) to thermal expansion strain.
///
/// ε = α × (T - T_ref), where T_ref is typically 293.15 K (20°C).
#[must_use]
#[inline]
pub fn temperature_to_thermal_strain(
    temperature_k: f64,
    expansion_coeff_per_k: f64,
    reference_temp_k: f64,
) -> f64 {
    expansion_coeff_per_k * (temperature_k - reference_temp_k)
}

/// Convert thermal gradient (K/m) and Young's modulus (Pa) and expansion
/// coefficient (1/K) to thermal stress (Pa).
///
/// σ = E × α × ΔT/Δx × L (for a constrained bar of length L).
/// Simplified: σ = E × α × gradient × characteristic_length.
#[must_use]
#[inline]
pub fn thermal_gradient_to_stress(
    gradient_k_per_m: f64,
    youngs_modulus_pa: f64,
    expansion_coeff_per_k: f64,
    length_m: f64,
) -> f64 {
    youngs_modulus_pa * expansion_coeff_per_k * gradient_k_per_m * length_m
}

// ── Badal bridges (weather/atmosphere) ──────────────────────────────────────

/// Convert altitude (m) and lapse rate (K/m) to atmospheric temperature (K).
///
/// T(h) = T_surface - Γ × h, where Γ is the lapse rate.
/// Standard tropospheric lapse rate ≈ 0.0065 K/m.
#[must_use]
#[inline]
pub fn altitude_to_temperature(
    surface_temperature_k: f64,
    altitude_m: f64,
    lapse_rate_k_per_m: f64,
) -> f64 {
    (surface_temperature_k - lapse_rate_k_per_m * altitude_m).max(0.0)
}

/// Convert humidity ratio (kg water / kg dry air) and dry-bulb
/// temperature (K) to approximate wet-bulb temperature (K).
///
/// Stull (2011) empirical formula for wet-bulb from RH and T.
/// `relative_humidity_fraction`: 0.0–1.0.
#[must_use]
pub fn humidity_to_wet_bulb(dry_bulb_k: f64, relative_humidity_fraction: f64) -> f64 {
    let t_c = dry_bulb_k - 273.15;
    let rh = (relative_humidity_fraction * 100.0).clamp(0.0, 100.0);
    // Stull (2011) formula
    let tw_c = t_c * (0.151_977 * (rh + 8.313_659).sqrt()).atan() + (t_c + rh).atan()
        - (rh - 1.676_331).atan()
        + 0.003_918_38 * rh.powf(1.5) * (0.023_101 * rh).atan()
        - 4.686_035;
    tw_c + 273.15
}

// ── Pravash bridges (fluid dynamics) ───────────────────────────────────────

/// Convert fluid velocity (m/s) to convective heat transfer coefficient
/// (W/(m²·K)) using a simplified turbulent flat plate correlation.
///
/// h ≈ C × V^0.8 for turbulent forced convection over a plate.
/// `fluid_conductivity_w_per_m_k`: thermal conductivity of the fluid.
/// `characteristic_length_m`: plate length.
/// `fluid_viscosity_pa_s`: dynamic viscosity.
#[must_use]
pub fn velocity_to_convection_coefficient(
    velocity_ms: f64,
    fluid_density_kg_m3: f64,
    fluid_viscosity_pa_s: f64,
    fluid_conductivity_w_per_m_k: f64,
    specific_heat_j_per_kg_k: f64,
    characteristic_length_m: f64,
) -> f64 {
    if fluid_viscosity_pa_s <= 0.0 || characteristic_length_m <= 0.0 {
        return 0.0;
    }
    let re =
        fluid_density_kg_m3 * velocity_ms.abs() * characteristic_length_m / fluid_viscosity_pa_s;
    let pr = specific_heat_j_per_kg_k * fluid_viscosity_pa_s / fluid_conductivity_w_per_m_k;
    // Dittus-Boelter: Nu = 0.023 × Re^0.8 × Pr^0.4 (heating)
    let nu = 0.023 * re.powf(0.8) * pr.powf(0.4);
    nu * fluid_conductivity_w_per_m_k / characteristic_length_m
}

/// Convert turbulent kinetic energy (m²/s²) to eddy thermal diffusivity (m²/s).
///
/// α_t = ν_t / Pr_t, using C_μ^(1/4) × sqrt(k) × L model.
#[must_use]
#[inline]
pub fn tke_to_eddy_thermal_diffusivity(tke: f64, length_scale_m: f64) -> f64 {
    if tke <= 0.0 {
        return 0.0;
    }
    let c_mu_quarter = 0.5477;
    let pr_t = 0.85;
    c_mu_quarter * tke.sqrt() * length_scale_m / pr_t
}

// ── Sharira bridges (physiology/biomechanics) ──────────────────────────────

/// Convert muscle mechanical work rate (W) to metabolic heat generation (W).
///
/// Muscles are ~25% efficient; ~75% of metabolic energy becomes heat.
/// P_heat = P_mechanical / efficiency - P_mechanical
/// = P_mechanical × (1/η - 1)
#[must_use]
#[inline]
pub fn muscle_work_to_heat(mechanical_power_w: f64, efficiency: f64) -> f64 {
    if efficiency <= 0.0 || efficiency >= 1.0 {
        return mechanical_power_w.abs() * 3.0; // default ~25% efficiency
    }
    (mechanical_power_w.abs() / efficiency - mechanical_power_w.abs()).max(0.0)
}

/// Estimate basal metabolic rate (W) from body mass (kg).
///
/// Kleiber's law: BMR = 3.5 × m^0.75 (watts).
/// Applies to mammals and birds across 7 orders of magnitude.
#[must_use]
#[inline]
pub fn basal_metabolic_rate(mass_kg: f64) -> f64 {
    if mass_kg <= 0.0 {
        return 0.0;
    }
    3.5 * mass_kg.powf(0.75)
}

/// Estimate body surface area (m²) from mass (kg) and height (m).
///
/// Du Bois formula: BSA = 0.007184 × m^0.425 × h_cm^0.725
/// (height is converted from metres to centimetres internally).
#[must_use]
pub fn body_surface_area(mass_kg: f64, height_m: f64) -> f64 {
    if mass_kg <= 0.0 || height_m <= 0.0 {
        return 0.0;
    }
    let height_cm = height_m * 100.0;
    0.007184 * mass_kg.powf(0.425) * height_cm.powf(0.725)
}

/// Convert skin surface area (m²) and temperatures to radiative heat loss (W).
///
/// Uses Stefan-Boltzmann: Q = ε × σ × A × (T_skin⁴ - T_env⁴)
/// Skin emissivity ≈ 0.98 for all skin types.
#[must_use]
pub fn skin_radiation_loss(
    surface_area_m2: f64,
    skin_temperature_k: f64,
    environment_temperature_k: f64,
) -> f64 {
    const STEFAN_BOLTZMANN: f64 = 5.670_374_419e-8;
    const SKIN_EMISSIVITY: f64 = 0.98;
    if surface_area_m2 <= 0.0 || skin_temperature_k <= 0.0 || environment_temperature_k <= 0.0 {
        return 0.0;
    }
    SKIN_EMISSIVITY
        * STEFAN_BOLTZMANN
        * surface_area_m2
        * (skin_temperature_k.powi(4) - environment_temperature_k.powi(4))
}

/// Convert muscle activation level (0-1) and max force (N) to estimated
/// heat generation rate (W) from muscle contraction.
///
/// Heat = force × shortening_velocity × (1 - efficiency).
/// Simplified: at moderate activation, heat ≈ 0.15 × activation × max_force (W).
#[must_use]
#[inline]
pub fn muscle_activation_heat(activation: f64, max_force_n: f64) -> f64 {
    let a = activation.clamp(0.0, 1.0);
    // Empirical: ~0.15 W per Newton of max force at full activation
    0.15 * a * max_force_n.abs()
}

/// Estimate evaporative cooling rate (W) from sweat rate (kg/s).
///
/// Q_evap = ṁ × L_v, where L_v ≈ 2,430,000 J/kg (latent heat of sweat).
#[must_use]
#[inline]
pub fn evaporative_cooling(sweat_rate_kg_per_s: f64) -> f64 {
    const LATENT_HEAT_SWEAT: f64 = 2_430_000.0; // J/kg
    sweat_rate_kg_per_s.max(0.0) * LATENT_HEAT_SWEAT
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Bijli ──────────────────────────────────────────────────────────

    #[test]
    fn joule_heating_basic() {
        // 100W into 1kg water (c_p=4186) → 0.0239 K/s
        let rate = joule_heating_to_temperature_rate(100.0, 1.0, 4186.0);
        assert!((rate - 0.0239).abs() < 0.001);
    }

    #[test]
    fn joule_heating_zero_mass() {
        assert_eq!(joule_heating_to_temperature_rate(100.0, 0.0, 4186.0), 0.0);
    }

    #[test]
    fn em_absorption_clamps() {
        assert_eq!(em_absorption_to_heat_source(-100.0), 0.0);
        assert_eq!(em_absorption_to_heat_source(500.0), 500.0);
    }

    // ── Kimiya ─────────────────────────────────────────────────────────

    #[test]
    fn reaction_heat_exothermic() {
        // ΔH = -100 kJ/mol, 0.1 mol/s → 10 kW heat release
        let q = reaction_heat_rate(-100_000.0, 0.1);
        assert!((q - 10_000.0).abs() < 0.1);
    }

    #[test]
    fn equilibrium_temp_basic() {
        // ΔG = -10000 J/mol, K = 100 → T = 10000 / (8.314 × ln(100)) ≈ 261 K
        let t = equilibrium_temperature(-10_000.0, 100.0);
        assert!((t - 261.3).abs() < 1.0);
    }

    #[test]
    fn equilibrium_temp_invalid_k() {
        assert_eq!(equilibrium_temperature(-10_000.0, 0.0), 0.0);
    }

    // ── Dravya ─────────────────────────────────────────────────────────

    #[test]
    fn thermal_strain_steel() {
        // Steel α ≈ 12e-6, T = 393K, ref = 293K → ε = 0.0012
        let e = temperature_to_thermal_strain(393.15, 12e-6, 293.15);
        assert!((e - 0.0012).abs() < 1e-6);
    }

    #[test]
    fn thermal_stress_basic() {
        // gradient = 10 K/m, E = 200 GPa, α = 12e-6, L = 1m
        // σ = 200e9 × 12e-6 × 10 × 1 = 24,000,000 Pa = 24 MPa
        let s = thermal_gradient_to_stress(10.0, 200e9, 12e-6, 1.0);
        assert!((s - 24_000_000.0).abs() < 1.0);
    }

    // ── Badal ──────────────────────────────────────────────────────────

    #[test]
    fn altitude_temperature_sea_level() {
        let t = altitude_to_temperature(288.15, 0.0, 0.0065);
        assert!((t - 288.15).abs() < 0.01);
    }

    #[test]
    fn altitude_temperature_1km() {
        let t = altitude_to_temperature(288.15, 1000.0, 0.0065);
        assert!((t - 281.65).abs() < 0.01);
    }

    #[test]
    fn wet_bulb_dry_conditions() {
        // At low humidity, wet-bulb should be lower than dry-bulb
        let tw = humidity_to_wet_bulb(300.0, 0.2);
        assert!(tw < 300.0);
    }

    #[test]
    fn wet_bulb_saturated() {
        // At 100% humidity, wet-bulb ≈ dry-bulb
        let tw = humidity_to_wet_bulb(300.0, 1.0);
        assert!((tw - 300.0).abs() < 2.0);
    }

    // ── Pravash ────────────────────────────────────────────────────────

    #[test]
    fn convection_coefficient_positive() {
        // Water at 1 m/s
        let h = velocity_to_convection_coefficient(1.0, 1000.0, 0.001, 0.6, 4186.0, 0.1);
        assert!(h > 0.0);
    }

    #[test]
    fn convection_coefficient_zero_velocity() {
        let h = velocity_to_convection_coefficient(0.0, 1000.0, 0.001, 0.6, 4186.0, 0.1);
        // Re=0 → Nu≈0 → h≈0
        assert!((0.0..1.0).contains(&h));
    }

    #[test]
    fn convection_coefficient_zero_viscosity() {
        assert_eq!(
            velocity_to_convection_coefficient(1.0, 1000.0, 0.0, 0.6, 4186.0, 0.1),
            0.0
        );
    }

    #[test]
    fn eddy_diffusivity_positive() {
        let alpha = tke_to_eddy_thermal_diffusivity(1.0, 0.1);
        assert!(alpha > 0.0);
    }

    #[test]
    fn eddy_diffusivity_zero_tke() {
        assert_eq!(tke_to_eddy_thermal_diffusivity(0.0, 0.1), 0.0);
    }

    // ── Sharira ────────────────────────────────────────────────────────

    #[test]
    fn muscle_heat_basic() {
        // 100W mechanical at 25% efficiency → 300W heat
        let heat = muscle_work_to_heat(100.0, 0.25);
        assert!((heat - 300.0).abs() < 0.1);
    }

    #[test]
    fn muscle_heat_default_efficiency() {
        let heat = muscle_work_to_heat(100.0, 0.0);
        assert!((heat - 300.0).abs() < 0.1); // defaults to 25%
    }

    #[test]
    fn bmr_human() {
        // 70kg human: BMR ≈ 3.5 × 70^0.75 ≈ 84.7W
        let bmr = basal_metabolic_rate(70.0);
        assert!((bmr - 84.7).abs() < 1.0);
    }

    #[test]
    fn bmr_zero_mass() {
        assert_eq!(basal_metabolic_rate(0.0), 0.0);
    }

    #[test]
    fn bsa_human() {
        // 70kg, 1.75m → BSA ≈ 1.85 m² (Du Bois)
        let bsa = body_surface_area(70.0, 1.75);
        assert!((bsa - 1.85).abs() < 0.1);
    }

    #[test]
    fn bsa_zero_inputs() {
        assert_eq!(body_surface_area(0.0, 1.75), 0.0);
        assert_eq!(body_surface_area(70.0, 0.0), 0.0);
    }

    #[test]
    fn skin_radiation_positive() {
        // 1.85 m², skin 310K, env 293K → should be positive (losing heat)
        let q = skin_radiation_loss(1.85, 310.0, 293.0);
        assert!(q > 0.0);
        // Should be ~100-200W range for a human
        assert!(q > 50.0 && q < 300.0, "radiation loss = {q}W");
    }

    #[test]
    fn skin_radiation_zero_area() {
        assert_eq!(skin_radiation_loss(0.0, 310.0, 293.0), 0.0);
    }

    #[test]
    fn muscle_activation_heat_basic() {
        let h = muscle_activation_heat(0.5, 1000.0);
        assert!((h - 75.0).abs() < 0.1); // 0.15 × 0.5 × 1000
    }

    #[test]
    fn evaporative_cooling_basic() {
        // 0.5 L/hr = 0.5/3600 kg/s ≈ 0.000139 kg/s → ~337W
        let q = evaporative_cooling(0.5 / 3600.0);
        assert!((q - 337.5).abs() < 1.0);
    }
}
