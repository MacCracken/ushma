//! Integration tests for ushma — cross-module thermodynamics.

use ushma::entropy;
use ushma::material;
use ushma::state;
use ushma::transfer;

#[test]
fn heat_conduction_through_copper_wall() {
    // 10cm copper wall, 1m² area, 100°C to 20°C
    let q = transfer::conduction(
        material::COPPER.conductivity,
        1.0,    // 1 m²
        373.15, // 100°C
        293.15, // 20°C
        0.1,    // 10 cm
    )
    .unwrap();

    // q = 401 * 1.0 * 80 / 0.1 = 320,800 W
    assert!((q - 320_800.0).abs() < 1.0);
}

#[test]
fn carnot_engine_efficiency() {
    // Steam engine: T_hot=500K, T_cold=300K → η=40%
    let eta = entropy::carnot_efficiency(500.0, 300.0).unwrap();
    assert!((eta - 0.4).abs() < 1e-10);

    // Real engines are always less efficient than Carnot
    let real_efficiency = 0.35;
    assert!(real_efficiency < eta);
}

#[test]
fn ideal_gas_law_roundtrip() {
    let n = 2.0;
    let t = 350.0;

    // Calculate volume at atmospheric pressure
    let v = state::ideal_gas_volume(n, t, state::ATM).unwrap();

    // Recover pressure from volume
    let p = state::ideal_gas_pressure(n, t, v).unwrap();
    assert!((p - state::ATM).abs() / state::ATM < 1e-10);

    // Recover temperature from P and V
    let t_back = state::ideal_gas_temperature(p, v, n).unwrap();
    assert!((t_back - t).abs() < 1e-10);
}

#[test]
fn thermal_equilibrium_lumped_capacitance() {
    let t_initial = 500.0; // hot object
    let t_env = 293.15; // room temperature
    let tau = 60.0; // 1 minute time constant

    // Temperature approaches environment over time
    let t_1min = transfer::lumped_capacitance(t_initial, t_env, 60.0, tau);
    let t_5min = transfer::lumped_capacitance(t_initial, t_env, 300.0, tau);
    let t_30min = transfer::lumped_capacitance(t_initial, t_env, 1800.0, tau);

    // Temperature must decrease monotonically toward t_env
    assert!(t_1min < t_initial);
    assert!(t_5min < t_1min);
    assert!(t_30min < t_5min);

    // After many time constants, should be very close to environment
    assert!((t_30min - t_env).abs() < 0.01);
}

#[test]
fn energy_conservation_adiabatic_process() {
    let t1 = 300.0;
    let v1 = 0.02;
    let v2 = 0.01;
    let gamma = 1.4; // diatomic gas

    // Adiabatic compression raises temperature
    let t2 = state::adiabatic_temperature(t1, v1, v2, gamma).unwrap();
    assert!(t2 > t1);

    // Expanding back should recover original temperature
    let t_back = state::adiabatic_temperature(t2, v2, v1, gamma).unwrap();
    assert!((t_back - t1).abs() / t1 < 1e-10);
}

#[test]
fn entropy_increase_irreversible_heat_transfer() {
    let q = 1000.0; // 1 kJ transferred
    let t_hot = 500.0;
    let t_cold = 300.0;

    // Hot body loses entropy
    let ds_hot = entropy::heat_transfer_entropy(-q, t_hot).unwrap();
    // Cold body gains entropy
    let ds_cold = entropy::heat_transfer_entropy(q, t_cold).unwrap();

    // Total entropy change must be positive (irreversible)
    let ds_total = ds_hot + ds_cold;
    assert!(ds_total > 0.0);
    assert!(entropy::is_spontaneous(ds_total));
}

#[test]
fn material_thermal_diffusivity_matches_manual() {
    // Copper: α = k/(ρ⋅c_p) = 401/(8960*385)
    let alpha_manual =
        material::COPPER.conductivity / (material::COPPER.density * material::COPPER.specific_heat);
    let alpha_method = material::COPPER.diffusivity();
    assert!((alpha_manual - alpha_method).abs() < 1e-15);

    // Cross-check with transfer module function
    let alpha_fn = transfer::thermal_diffusivity(
        material::COPPER.conductivity,
        material::COPPER.density,
        material::COPPER.specific_heat,
    )
    .unwrap();
    assert!((alpha_fn - alpha_method).abs() < 1e-15);
}

#[test]
fn van_der_waals_vs_ideal_gas() {
    let n = 1.0;
    let t = state::STANDARD_TEMP;
    let v = 0.02241;

    let p_ideal = state::ideal_gas_pressure(n, t, v).unwrap();

    // CO₂ van der Waals parameters
    let a = 0.3658;
    let b = 4.286e-5;
    let p_vdw = state::van_der_waals_pressure(n, t, v, a, b).unwrap();

    // Van der Waals pressure should differ from ideal
    // (intermolecular attractions reduce pressure)
    assert!(p_vdw < p_ideal);

    // But at standard conditions the difference should be small
    let relative_diff = (p_ideal - p_vdw).abs() / p_ideal;
    assert!(relative_diff < 0.05); // within 5%
}
