//! Chemical thermodynamics — formation data, Hess's law, equilibrium, flame temperature.
//!
//! All SI units: joules, moles, kelvins. Standard state: 298.15 K, 1 atm.

use std::borrow::Cow;

use serde::{Deserialize, Serialize};

use crate::error::{Result, UshmaError};
use crate::state::GAS_CONSTANT;

/// Standard reference temperature (K).
pub const STANDARD_T: f64 = 298.15;

/// Thermodynamic formation data for a chemical species at standard state (298.15 K, 1 atm).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Species {
    /// Species name / formula.
    pub name: Cow<'static, str>,
    /// Standard enthalpy of formation ΔH°f (J/mol).
    pub delta_hf: f64,
    /// Standard Gibbs energy of formation ΔG°f (J/mol).
    pub delta_gf: f64,
    /// Standard molar heat capacity C°p (J/(mol·K)).
    pub cp: f64,
}

// --- Standard formation data (NIST Chemistry WebBook, 298.15 K) ---

/// Hydrogen gas (H₂) — reference element.
pub const H2: Species = Species {
    name: Cow::Borrowed("H2"),
    delta_hf: 0.0,
    delta_gf: 0.0,
    cp: 28.84,
};
/// Oxygen gas (O₂) — reference element.
pub const O2: Species = Species {
    name: Cow::Borrowed("O2"),
    delta_hf: 0.0,
    delta_gf: 0.0,
    cp: 29.38,
};
/// Nitrogen gas (N₂) — reference element.
pub const N2: Species = Species {
    name: Cow::Borrowed("N2"),
    delta_hf: 0.0,
    delta_gf: 0.0,
    cp: 29.12,
};
/// Water vapor (H₂O, g).
pub const H2O_GAS: Species = Species {
    name: Cow::Borrowed("H2O(g)"),
    delta_hf: -241_820.0,
    delta_gf: -228_570.0,
    cp: 33.58,
};
/// Liquid water (H₂O, l).
pub const H2O_LIQUID: Species = Species {
    name: Cow::Borrowed("H2O(l)"),
    delta_hf: -285_830.0,
    delta_gf: -237_130.0,
    cp: 75.29,
};
/// Carbon dioxide (CO₂).
pub const CO2: Species = Species {
    name: Cow::Borrowed("CO2"),
    delta_hf: -393_510.0,
    delta_gf: -394_360.0,
    cp: 37.11,
};
/// Carbon monoxide (CO).
pub const CO: Species = Species {
    name: Cow::Borrowed("CO"),
    delta_hf: -110_530.0,
    delta_gf: -137_170.0,
    cp: 29.14,
};
/// Methane (CH₄).
pub const CH4: Species = Species {
    name: Cow::Borrowed("CH4"),
    delta_hf: -74_810.0,
    delta_gf: -50_720.0,
    cp: 35.31,
};
/// Ethane (C₂H₆).
pub const C2H6: Species = Species {
    name: Cow::Borrowed("C2H6"),
    delta_hf: -84_680.0,
    delta_gf: -32_820.0,
    cp: 52.63,
};
/// Propane (C₃H₈).
pub const C3H8: Species = Species {
    name: Cow::Borrowed("C3H8"),
    delta_hf: -103_850.0,
    delta_gf: -23_490.0,
    cp: 73.51,
};
/// Ammonia (NH₃).
pub const NH3: Species = Species {
    name: Cow::Borrowed("NH3"),
    delta_hf: -46_110.0,
    delta_gf: -16_450.0,
    cp: 35.06,
};
/// Nitric oxide (NO).
pub const NO: Species = Species {
    name: Cow::Borrowed("NO"),
    delta_hf: 90_250.0,
    delta_gf: 86_550.0,
    cp: 29.84,
};

/// All built-in species for iteration.
pub const ALL_SPECIES: &[&Species] = &[
    &H2,
    &O2,
    &N2,
    &H2O_GAS,
    &H2O_LIQUID,
    &CO2,
    &CO,
    &CH4,
    &C2H6,
    &C3H8,
    &NH3,
    &NO,
];

// --- Reaction thermodynamics ---

/// Standard reaction enthalpy via Hess's law: ΔH°rxn = Σ n·ΔH°f(products) - Σ n·ΔH°f(reactants) (J/mol).
///
/// Each entry is (stoichiometric coefficient, species reference).
#[must_use]
pub fn reaction_enthalpy(products: &[(f64, &Species)], reactants: &[(f64, &Species)]) -> f64 {
    let h_prod: f64 = products.iter().map(|(n, s)| n * s.delta_hf).sum();
    let h_react: f64 = reactants.iter().map(|(n, s)| n * s.delta_hf).sum();
    h_prod - h_react
}

/// Standard Gibbs energy of reaction: ΔG°rxn = Σ n·ΔG°f(products) - Σ n·ΔG°f(reactants) (J/mol).
#[must_use]
pub fn reaction_gibbs(products: &[(f64, &Species)], reactants: &[(f64, &Species)]) -> f64 {
    let g_prod: f64 = products.iter().map(|(n, s)| n * s.delta_gf).sum();
    let g_react: f64 = reactants.iter().map(|(n, s)| n * s.delta_gf).sum();
    g_prod - g_react
}

/// Gibbs energy from enthalpy and entropy: ΔG = ΔH - TΔS (J/mol).
#[inline]
#[must_use]
pub fn gibbs_from_enthalpy_entropy(delta_h: f64, delta_s: f64, temperature: f64) -> f64 {
    delta_h - temperature * delta_s
}

/// Equilibrium constant from standard Gibbs energy: K = exp(-ΔG°/(RT)).
#[tracing::instrument(level = "debug")]
pub fn equilibrium_constant(delta_g: f64, temperature: f64) -> Result<f64> {
    if temperature <= 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    Ok((-delta_g / (GAS_CONSTANT * temperature)).exp())
}

/// Van't Hoff equation: K at a new temperature from a reference K.
///
/// ln(K₂/K₁) = -(ΔH°/R)·(1/T₂ - 1/T₁)
pub fn vant_hoff_k(k_ref: f64, delta_h: f64, t_ref: f64, t_target: f64) -> Result<f64> {
    if k_ref <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("reference K={k_ref} must be positive"),
        });
    }
    if t_ref <= 0.0 {
        return Err(UshmaError::InvalidTemperature { kelvin: t_ref });
    }
    if t_target <= 0.0 {
        return Err(UshmaError::InvalidTemperature { kelvin: t_target });
    }
    let exponent = -(delta_h / GAS_CONSTANT) * (1.0 / t_target - 1.0 / t_ref);
    Ok(k_ref * exponent.exp())
}

/// Adiabatic flame temperature (K).
///
/// Finds T where the heat released by combustion equals the heat absorbed
/// by products: -ΔH°rxn + Σ n_r·Cp·(T_in - 298.15) = Σ n_p·Cp·(T_flame - 298.15)
///
/// Uses constant-Cp assumption (Cp at 298.15 K). Direct analytical solution.
#[tracing::instrument(level = "debug", skip(reactants, products))]
pub fn adiabatic_flame_temperature(
    reactants: &[(f64, &Species)],
    products: &[(f64, &Species)],
    t_initial: f64,
) -> Result<f64> {
    if t_initial <= 0.0 {
        return Err(UshmaError::InvalidTemperature { kelvin: t_initial });
    }

    let delta_h = reaction_enthalpy(products, reactants);
    let cp_products: f64 = products.iter().map(|(n, s)| n * s.cp).sum();

    if cp_products <= 0.0 {
        return Err(UshmaError::DivisionByZero {
            context: "product heat capacity sum must be positive".into(),
        });
    }

    // Sensible heat of reactants above standard state
    let q_reactants: f64 = reactants
        .iter()
        .map(|(n, s)| n * s.cp * (t_initial - STANDARD_T))
        .sum();

    // Energy balance: cp_products * (T_flame - 298.15) = -delta_h + q_reactants
    // T_flame = 298.15 + (-delta_h + q_reactants) / cp_products
    let t_flame = STANDARD_T + (-delta_h + q_reactants) / cp_products;

    if t_flame <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("computed flame temperature {t_flame:.1} K is non-positive"),
        });
    }

    Ok(t_flame)
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- Species data ---

    #[test]
    fn test_all_species_count() {
        assert_eq!(ALL_SPECIES.len(), 12);
    }

    #[test]
    fn test_elements_zero_formation() {
        assert_eq!(H2.delta_hf, 0.0);
        assert_eq!(O2.delta_hf, 0.0);
        assert_eq!(N2.delta_hf, 0.0);
        assert_eq!(H2.delta_gf, 0.0);
    }

    #[test]
    fn test_species_positive_cp() {
        for s in ALL_SPECIES {
            assert!(s.cp > 0.0, "{} Cp must be positive", s.name);
        }
    }

    #[test]
    fn test_water_formation_enthalpy() {
        // H₂O(l): -285.83 kJ/mol (NIST)
        assert!((H2O_LIQUID.delta_hf - (-285_830.0)).abs() < 10.0);
    }

    #[test]
    fn test_species_serde_roundtrip() {
        let json = serde_json::to_string(&CH4).unwrap();
        let back: Species = serde_json::from_str(&json).unwrap();
        assert_eq!(back.name, "CH4");
        assert!((back.delta_hf - (-74_810.0)).abs() < 1.0);
    }

    // --- Hess's law ---

    #[test]
    fn test_methane_combustion_enthalpy() {
        // CH₄ + 2O₂ → CO₂ + 2H₂O(g)
        // ΔH = (-393510 + 2*(-241820)) - (-74810 + 0) = -802340 J/mol
        let dh = reaction_enthalpy(&[(1.0, &CO2), (2.0, &H2O_GAS)], &[(1.0, &CH4), (2.0, &O2)]);
        assert!((dh - (-802_340.0)).abs() < 100.0);
    }

    #[test]
    fn test_hydrogen_combustion_enthalpy() {
        // H₂ + ½O₂ → H₂O(l): ΔH = -285830 J/mol
        let dh = reaction_enthalpy(&[(1.0, &H2O_LIQUID)], &[(1.0, &H2), (0.5, &O2)]);
        assert!((dh - (-285_830.0)).abs() < 10.0);
    }

    #[test]
    fn test_reaction_enthalpy_reverse() {
        // Reverse reaction has opposite ΔH
        let dh_fwd = reaction_enthalpy(&[(1.0, &CO2), (2.0, &H2O_GAS)], &[(1.0, &CH4), (2.0, &O2)]);
        let dh_rev = reaction_enthalpy(&[(1.0, &CH4), (2.0, &O2)], &[(1.0, &CO2), (2.0, &H2O_GAS)]);
        assert!((dh_fwd + dh_rev).abs() < 1e-10);
    }

    // --- Gibbs + equilibrium ---

    #[test]
    fn test_methane_combustion_gibbs() {
        let dg = reaction_gibbs(&[(1.0, &CO2), (2.0, &H2O_GAS)], &[(1.0, &CH4), (2.0, &O2)]);
        // Should be very negative (highly spontaneous combustion)
        assert!(dg < -700_000.0);
    }

    #[test]
    fn test_equilibrium_constant_spontaneous() {
        // Negative ΔG → K > 1
        let k = equilibrium_constant(-50_000.0, 298.15).unwrap();
        assert!(k > 1.0);
    }

    #[test]
    fn test_equilibrium_constant_nonspontaneous() {
        // Positive ΔG → K < 1
        let k = equilibrium_constant(50_000.0, 298.15).unwrap();
        assert!(k < 1.0);
    }

    #[test]
    fn test_equilibrium_constant_invalid() {
        assert!(equilibrium_constant(-50_000.0, 0.0).is_err());
    }

    #[test]
    fn test_gibbs_from_enthalpy_entropy() {
        let dg = gibbs_from_enthalpy_entropy(-100_000.0, -200.0, 298.15);
        // ΔG = -100000 - 298.15*(-200) = -100000 + 59630 = -40370
        assert!((dg - (-40_370.0)).abs() < 1.0);
    }

    // --- Van't Hoff ---

    #[test]
    fn test_vant_hoff_roundtrip() {
        let k1 = 100.0;
        let dh = -50_000.0;
        let k2 = vant_hoff_k(k1, dh, 298.15, 400.0).unwrap();
        let k_back = vant_hoff_k(k2, dh, 400.0, 298.15).unwrap();
        assert!((k_back - k1).abs() / k1 < 1e-10);
    }

    #[test]
    fn test_vant_hoff_exothermic() {
        // Exothermic (ΔH < 0): K decreases with temperature
        let k1 = 100.0;
        let k2 = vant_hoff_k(k1, -50_000.0, 298.15, 400.0).unwrap();
        assert!(k2 < k1);
    }

    #[test]
    fn test_vant_hoff_endothermic() {
        // Endothermic (ΔH > 0): K increases with temperature
        let k1 = 0.01;
        let k2 = vant_hoff_k(k1, 50_000.0, 298.15, 400.0).unwrap();
        assert!(k2 > k1);
    }

    #[test]
    fn test_vant_hoff_invalid() {
        assert!(vant_hoff_k(0.0, -50_000.0, 298.15, 400.0).is_err());
        assert!(vant_hoff_k(100.0, -50_000.0, 0.0, 400.0).is_err());
        assert!(vant_hoff_k(100.0, -50_000.0, 298.15, 0.0).is_err());
    }

    // --- Adiabatic flame temperature ---

    #[test]
    fn test_methane_air_flame_temperature() {
        // CH₄ + 2O₂ + 7.52N₂ → CO₂ + 2H₂O(g) + 7.52N₂
        // Constant-Cp model overestimates vs real (Cp increases with T).
        // Real ≈ 2230 K; constant-Cp at 298 K gives ~2780 K.
        let t_flame = adiabatic_flame_temperature(
            &[(1.0, &CH4), (2.0, &O2), (7.52, &N2)],
            &[(1.0, &CO2), (2.0, &H2O_GAS), (7.52, &N2)],
            298.15,
        )
        .unwrap();
        assert!(
            t_flame > 2500.0 && t_flame < 3200.0,
            "CH₄ flame T={t_flame:.0} K"
        );
    }

    #[test]
    fn test_hydrogen_flame_temperature() {
        // Constant-Cp model: real ≈ 2480 K, constant-Cp at 298 K gives ~3036 K.
        let t_flame = adiabatic_flame_temperature(
            &[(1.0, &H2), (0.5, &O2), (1.88, &N2)],
            &[(1.0, &H2O_GAS), (1.88, &N2)],
            298.15,
        )
        .unwrap();
        assert!(
            t_flame > 2700.0 && t_flame < 3500.0,
            "H₂ flame T={t_flame:.0} K"
        );
    }

    #[test]
    fn test_flame_temperature_invalid() {
        assert!(adiabatic_flame_temperature(&[(1.0, &CH4)], &[(1.0, &CO2)], 0.0,).is_err());
    }
}
