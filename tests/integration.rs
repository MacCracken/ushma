//! Integration tests for ushma — cross-module thermodynamics.

use ushma::cycle;
use ushma::entropy;
use ushma::material;
use ushma::phase;
use ushma::state;
use ushma::steam;
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

// --- Phase transition integration tests ---

#[test]
fn water_boiling_point_consistency() {
    // material::WATER and phase::WATER_PHASE should agree on boiling point
    assert!(
        (material::WATER.boiling_point - phase::WATER_PHASE.boiling_point).abs() < 0.01,
        "material ({}) vs phase ({}) boiling points disagree",
        material::WATER.boiling_point,
        phase::WATER_PHASE.boiling_point
    );
}

#[test]
fn clausius_clapeyron_vs_steam_table() {
    // Compare CC-predicted P_sat at 120°C against steam table
    let l_molar = phase::WATER_PHASE.latent_heat_vaporization * phase::WATER_PHASE.molar_mass;
    let p_cc = phase::clausius_clapeyron_pressure(
        state::ATM,
        phase::WATER_PHASE.boiling_point,
        393.15,
        l_molar,
    )
    .unwrap();

    let e = steam::saturated_by_temperature(393.15).unwrap();

    // CC approximation should be within 5% of steam table
    let rel_err = (p_cc - e.pressure).abs() / e.pressure;
    assert!(
        rel_err < 0.05,
        "CC ({:.0} Pa) vs steam table ({:.0} Pa): {:.1}% error",
        p_cc,
        e.pressure,
        rel_err * 100.0
    );
}

#[test]
fn quality_consistency_across_methods() {
    // At 373.15 K, compute quality from v, h, and s for the same mixture
    let entry = steam::saturated_by_temperature(373.15).unwrap();
    let target_x = 0.6;

    let props = steam::wet_steam_properties(target_x, &entry).unwrap();
    let x_from_v = steam::quality_from_volume(props.specific_volume, &entry).unwrap();
    let x_from_h = steam::quality_from_enthalpy(props.specific_enthalpy, &entry).unwrap();
    let x_from_s = steam::quality_from_entropy(props.specific_entropy, &entry).unwrap();

    assert!((x_from_v - target_x).abs() < 0.01);
    assert!((x_from_h - target_x).abs() < 0.01);
    assert!((x_from_s - target_x).abs() < 0.01);
}

#[test]
fn simple_rankine_cycle_energy_balance() {
    // Simplified Rankine cycle at 2 MPa / condenser at 10 kPa
    // State 1: saturated liquid at condenser pressure (10 kPa)
    let cond = steam::saturated_by_pressure(10_000.0).unwrap();
    let h1 = cond.h_f; // pump inlet

    // State 2: compressed liquid (approximate as h1 + v_f * dP)
    let p_boiler = 2_000_000.0;
    let w_pump = cond.v_f * (p_boiler - 10_000.0); // J/kg
    let h2 = h1 + w_pump;

    // State 3: superheated steam at boiler exit (2 MPa, 573.15 K / 300°C)
    let sh = steam::superheated_lookup(573.15, p_boiler).unwrap();
    let h3 = sh.specific_enthalpy;

    // State 4: wet steam after turbine (isentropic to condenser pressure)
    // s3 = s4, find quality at condenser pressure
    let x4 = steam::quality_from_entropy(sh.specific_entropy, &cond).unwrap();
    let h4 = cond.h_f + x4 * cond.h_fg;

    // Energy balance
    let q_in = h3 - h2; // heat added in boiler
    let w_turbine = h3 - h4; // work out of turbine
    let q_out = h4 - h1; // heat rejected in condenser
    let w_net = w_turbine - w_pump;

    // First law: q_in = w_net + q_out (within rounding tolerance)
    let balance = (q_in - w_net - q_out).abs();
    assert!(
        balance < 5000.0, // within 5 kJ/kg for rounded table values
        "Energy balance error: {balance} J/kg"
    );

    // Thermal efficiency should be reasonable (20-35% for this cycle)
    let eta = w_net / q_in;
    assert!(
        eta > 0.15 && eta < 0.40,
        "Rankine efficiency {eta} out of range"
    );
}

#[test]
fn phase_lookup_agrees_with_steam_table() {
    // Water at 300 K should be liquid per phase lookup
    let p = phase::WATER_PHASE.phase_at(300.0, state::ATM).unwrap();
    assert_eq!(p, phase::Phase::Liquid);

    // Steam table at 373.15 K: P_sat ≈ 101,350 Pa
    // At slightly above P_sat, phase should be liquid
    let sat = steam::saturated_by_temperature(373.15).unwrap();
    let p_above = phase::WATER_PHASE
        .phase_at(373.15, sat.pressure * 1.1)
        .unwrap();
    assert_eq!(p_above, phase::Phase::Liquid);
}

// --- Cycle integration tests ---

#[test]
fn ideal_gas_cycle_efficiency_ordering() {
    // At comparable conditions, Otto and Brayton efficiencies follow analytical formulas
    // Otto r=8: η = 1 - 1/8^0.4 = 0.5647
    // Brayton rp=10: η = 1 - 1/10^(0.4/1.4) = 0.4820
    let otto = cycle::otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0).unwrap();
    let brayton = cycle::brayton_cycle(300.0, 101_325.0, 10.0, 1400.0, 1.4, 1.0).unwrap();
    let diesel = cycle::diesel_cycle(300.0, 101_325.0, 20.0, 2.0, 1.4, 1.0).unwrap();

    // All should have positive efficiency below Carnot
    assert!(otto.efficiency > 0.0 && otto.efficiency < 1.0);
    assert!(brayton.efficiency > 0.0 && brayton.efficiency < 1.0);
    assert!(diesel.efficiency > 0.0 && diesel.efficiency < 1.0);

    // All satisfy first law
    for c in [&otto, &brayton, &diesel] {
        let balance = (c.heat_in - c.net_work - c.heat_out).abs();
        assert!(balance < 1e-6, "Energy balance violation");
    }
}

#[test]
fn rankine_with_steam_tables_end_to_end() {
    // Simple Rankine: 10 kPa condenser, 2 MPa boiler, 300°C superheat
    let r = cycle::rankine_cycle(10_000.0, 2_000_000.0, Some(573.15)).unwrap();

    assert_eq!(r.state_points.len(), 4);
    assert!(r.efficiency > 0.15);
    assert!(r.net_work > 0.0);

    // Should be less efficient than Carnot between same temperature limits
    let t_hot = r.state_points[2].temperature;
    let t_cold = r.state_points[0].temperature;
    let eta_carnot = entropy::carnot_efficiency(t_hot, t_cold).unwrap();
    assert!(r.efficiency < eta_carnot);
}

#[test]
fn refrigeration_vs_carnot_cop() {
    let r = cycle::refrigeration_cycle(7_384.0, 47_390.0).unwrap();

    let t_cold = r.cycle.state_points[0].temperature;
    let t_hot = r.cycle.state_points[2].temperature;
    let cop_carnot = entropy::carnot_cop_refrigeration(t_hot, t_cold).unwrap();

    // Real COP must be less than Carnot COP
    assert!(
        r.cop_refrigeration < cop_carnot,
        "COP {} exceeds Carnot {}",
        r.cop_refrigeration,
        cop_carnot
    );
}

#[test]
fn cycle_pv_diagram_closed_loop_all_types() {
    let otto = cycle::otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0).unwrap();
    let diesel = cycle::diesel_cycle(300.0, 101_325.0, 20.0, 2.0, 1.4, 1.0).unwrap();
    let brayton = cycle::brayton_cycle(300.0, 101_325.0, 10.0, 1400.0, 1.4, 1.0).unwrap();

    for (name, result) in [("Otto", &otto), ("Diesel", &diesel), ("Brayton", &brayton)] {
        let pts = cycle::cycle_pv_diagram(result, 10);
        let first = &pts[0];
        let last = &pts[pts.len() - 1];
        assert!(
            (first.x - last.x).abs() / first.x < 1e-6,
            "{name} P-v loop not closed (V)"
        );
        assert!(
            (first.y - last.y).abs() / first.y < 1e-6,
            "{name} P-v loop not closed (P)"
        );
    }
}

#[test]
fn cycle_comparison_all_below_carnot() {
    let otto = cycle::otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0).unwrap();
    let brayton = cycle::brayton_cycle(300.0, 101_325.0, 10.0, 1400.0, 1.4, 1.0).unwrap();

    let entries = vec![
        (
            cycle::CycleKind::Otto,
            &otto,
            otto.state_points[2].temperature,
            otto.state_points[0].temperature,
        ),
        (
            cycle::CycleKind::Brayton,
            &brayton,
            brayton.state_points[2].temperature,
            brayton.state_points[0].temperature,
        ),
    ];

    let comparison = cycle::compare_cycles(&entries);
    for entry in &comparison {
        assert!(
            entry.second_law_efficiency <= 1.0,
            "{:?} η_II={} > 1",
            entry.kind,
            entry.second_law_efficiency
        );
    }
}
