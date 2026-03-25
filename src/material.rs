//! Thermal material properties — conductivity, specific heat, density.
//!
//! Reference data for common materials at approximately 300 K.

use std::borrow::Cow;

use serde::{Deserialize, Serialize};

/// Thermal properties of a material.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThermalMaterial {
    /// Material name.
    pub name: Cow<'static, str>,
    /// Thermal conductivity k (W/(m⋅K)).
    pub conductivity: f64,
    /// Specific heat capacity c_p (J/(kg⋅K)).
    pub specific_heat: f64,
    /// Density ρ (kg/m³).
    pub density: f64,
    /// Melting point (K). 0.0 if not applicable.
    pub melting_point: f64,
    /// Boiling point (K). 0.0 if not applicable.
    pub boiling_point: f64,
}

impl ThermalMaterial {
    /// Thermal diffusivity α = k/(ρ⋅c_p) (m²/s).
    #[inline]
    #[must_use]
    pub fn diffusivity(&self) -> f64 {
        self.conductivity / (self.density * self.specific_heat)
    }

    /// Volumetric heat capacity ρ⋅c_p (J/(m³⋅K)).
    #[inline]
    #[must_use]
    pub fn volumetric_heat_capacity(&self) -> f64 {
        self.density * self.specific_heat
    }
}

// --- Common materials at ~300 K ---

pub const COPPER: ThermalMaterial = ThermalMaterial {
    name: Cow::Borrowed("Copper"),
    conductivity: 401.0,
    specific_heat: 385.0,
    density: 8960.0,
    melting_point: 1358.0,
    boiling_point: 2835.0,
};

pub const ALUMINUM: ThermalMaterial = ThermalMaterial {
    name: Cow::Borrowed("Aluminum"),
    conductivity: 237.0,
    specific_heat: 897.0,
    density: 2700.0,
    melting_point: 933.5,
    boiling_point: 2743.0,
};

pub const IRON: ThermalMaterial = ThermalMaterial {
    name: Cow::Borrowed("Iron"),
    conductivity: 80.2,
    specific_heat: 449.0,
    density: 7874.0,
    melting_point: 1811.0,
    boiling_point: 3134.0,
};

pub const STEEL: ThermalMaterial = ThermalMaterial {
    name: Cow::Borrowed("Steel (1% carbon)"),
    conductivity: 43.0,
    specific_heat: 490.0,
    density: 7850.0,
    melting_point: 1700.0,
    boiling_point: 3000.0,
};

pub const GLASS: ThermalMaterial = ThermalMaterial {
    name: Cow::Borrowed("Glass (soda-lime)"),
    conductivity: 1.0,
    specific_heat: 840.0,
    density: 2500.0,
    melting_point: 1000.0,
    boiling_point: 0.0,
};

pub const WATER: ThermalMaterial = ThermalMaterial {
    name: Cow::Borrowed("Water"),
    conductivity: 0.606,
    specific_heat: 4186.0,
    density: 997.0,
    melting_point: 273.15,
    boiling_point: 373.15,
};

pub const AIR: ThermalMaterial = ThermalMaterial {
    name: Cow::Borrowed("Air"),
    conductivity: 0.026,
    specific_heat: 1005.0,
    density: 1.225,
    melting_point: 0.0,
    boiling_point: 0.0,
};

pub const WOOD_OAK: ThermalMaterial = ThermalMaterial {
    name: Cow::Borrowed("Wood (oak)"),
    conductivity: 0.17,
    specific_heat: 2380.0,
    density: 600.0,
    melting_point: 0.0,
    boiling_point: 0.0,
};

pub const CONCRETE: ThermalMaterial = ThermalMaterial {
    name: Cow::Borrowed("Concrete"),
    conductivity: 1.7,
    specific_heat: 880.0,
    density: 2300.0,
    melting_point: 0.0,
    boiling_point: 0.0,
};

pub const DIAMOND: ThermalMaterial = ThermalMaterial {
    name: Cow::Borrowed("Diamond"),
    conductivity: 2200.0,
    specific_heat: 509.0,
    density: 3510.0,
    melting_point: 3820.0,
    boiling_point: 5100.0,
};

/// All built-in materials for iteration.
pub const ALL_MATERIALS: &[&ThermalMaterial] = &[
    &COPPER, &ALUMINUM, &IRON, &STEEL, &GLASS, &WATER, &AIR, &WOOD_OAK, &CONCRETE, &DIAMOND,
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_copper_diffusivity() {
        let alpha = COPPER.diffusivity();
        assert!((alpha - 1.16e-4).abs() / 1.16e-4 < 0.01);
    }

    #[test]
    fn test_water_specific_heat() {
        assert!((WATER.specific_heat - 4186.0).abs() < 1.0);
    }

    #[test]
    fn test_diamond_highest_conductivity() {
        for mat in ALL_MATERIALS {
            assert!(DIAMOND.conductivity >= mat.conductivity);
        }
    }

    #[test]
    fn test_volumetric_heat_capacity() {
        let vc = WATER.volumetric_heat_capacity();
        // Water: ~4.17 MJ/(m³⋅K)
        assert!((vc - 4_173_342.0).abs() / 4_173_342.0 < 0.01);
    }

    #[test]
    fn test_all_materials_count() {
        assert_eq!(ALL_MATERIALS.len(), 10);
    }

    #[test]
    fn test_material_serde_roundtrip() {
        let json = serde_json::to_string(&COPPER).unwrap();
        let back: ThermalMaterial = serde_json::from_str(&json).unwrap();
        assert_eq!(back.name, "Copper");
        assert!((back.conductivity - 401.0).abs() < 1e-10);
    }

    #[test]
    fn test_air_low_conductivity() {
        let k = AIR.conductivity;
        assert!(k < 0.1);
    }

    #[test]
    fn test_water_phase_transitions() {
        assert!((WATER.melting_point - 273.15).abs() < 0.01);
        assert!((WATER.boiling_point - 373.15).abs() < 0.01);
    }
}
