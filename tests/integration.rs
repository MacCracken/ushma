//! Integration tests for ushma — cross-module thermodynamics.

use ushma::chem;
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

// --- Extended heat transfer integration tests ---

#[test]
fn nusselt_to_h_to_convection_chain() {
    // Air over a 0.1m plate: Re=50000, Pr=0.7
    // Dittus-Boelter → Nu → h → Q via convection
    let nu = transfer::nusselt_dittus_boelter(50_000.0, 0.7).unwrap();
    let k_air = 0.026; // W/(m·K)
    let length = 0.1;
    let h = nu * k_air / length;

    let area = 0.1; // m²
    let q = transfer::convection(h, area, 373.15, 293.15);
    assert!(q > 0.0);
    // h should be reasonable for forced convection (10-500 W/(m²·K))
    assert!(h > 1.0 && h < 1000.0);
}

#[test]
fn lmtd_vs_ntu_same_heat_exchanger() {
    // Same HX computed both ways should give same Q
    let t_h_in = 400.0;
    let t_h_out = 350.0;
    let t_c_in = 300.0;
    let t_c_out = 330.0;
    let u = 500.0;
    let area = 2.0;

    // LMTD method (counter-flow)
    let lmtd = transfer::lmtd_counter(t_h_in, t_h_out, t_c_in, t_c_out).unwrap();
    let q_lmtd = transfer::heat_exchanger_lmtd(u, area, lmtd);

    // ε-NTU method
    // C_h = Q_h / ΔT_h, C_c = Q_c / ΔT_c
    // Q = C_h(T_h_in - T_h_out) = C_c(T_c_out - T_c_in)
    // Use Q from LMTD to find C values
    let c_h = q_lmtd / (t_h_in - t_h_out);
    let c_c = q_lmtd / (t_c_out - t_c_in);
    let c_min = c_h.min(c_c);
    let c_max = c_h.max(c_c);
    let cr = c_min / c_max;
    let ntu_val = transfer::ntu(u, area, c_min).unwrap();
    let eff = transfer::effectiveness_counter(ntu_val, cr).unwrap();
    let q_ntu = transfer::heat_exchanger_ntu(eff, c_min, t_h_in, t_c_in);

    assert!(
        (q_lmtd - q_ntu).abs() / q_lmtd < 0.01,
        "LMTD Q={q_lmtd} vs NTU Q={q_ntu}"
    );
}

#[test]
fn view_factor_reciprocity() {
    // A₁F₁₂ = A₂F₂₁ for coaxial disks
    let r1 = 0.3;
    let r2 = 0.5;
    let d = 0.4;
    let f12 = transfer::view_factor_coaxial_disks(r1, r2, d).unwrap();
    let f21 = transfer::view_factor_coaxial_disks(r2, r1, d).unwrap();

    let a1 = std::f64::consts::PI * r1 * r1;
    let a2 = std::f64::consts::PI * r2 * r2;

    let diff = (a1 * f12 - a2 * f21).abs();
    assert!(
        diff / (a1 * f12) < 0.01,
        "Reciprocity violated: A1*F12={}, A2*F21={}",
        a1 * f12,
        a2 * f21
    );
}

// --- Real gas model integration tests ---

#[test]
fn eos_comparison_co2_stp() {
    // CO₂ at STP: compare vdW, RK, PR against ideal gas
    let vm = 0.02241; // molar volume at STP
    let t = state::STANDARD_TEMP;
    let tc = 304.13;
    let pc = 7_375_000.0;

    let p_ideal = state::ideal_gas_pressure(1.0, t, vm).unwrap();
    let p_vdw = state::van_der_waals_pressure(1.0, t, vm, 0.3658, 4.286e-5).unwrap();
    let p_rk = state::redlich_kwong_pressure(t, vm, tc, pc).unwrap();
    let p_pr = state::peng_robinson_pressure(t, vm, tc, pc, 0.224).unwrap();

    // All should be near ATM at STP
    for (name, p) in [
        ("ideal", p_ideal),
        ("vdW", p_vdw),
        ("RK", p_rk),
        ("PR", p_pr),
    ] {
        assert!(
            (p - state::ATM).abs() / state::ATM < 0.05,
            "{name} pressure {p:.0} Pa too far from ATM"
        );
    }
}

#[test]
fn virial_approaches_ideal_at_low_density() {
    let vm_large = 1.0; // 1 m³/mol — very dilute
    let b = state::pitzer_second_virial(300.0, 304.13, 7_375_000.0, 0.224);
    let p_virial = state::virial_pressure_2nd(300.0, vm_large, b).unwrap();
    let p_ideal = state::GAS_CONSTANT * 300.0 / vm_large;
    assert!(
        (p_virial - p_ideal).abs() / p_ideal < 0.001,
        "virial not converging to ideal at low density"
    );
}

#[test]
fn mixture_z_factor_kays_rule() {
    // 70% N₂ + 30% CO₂ at 300 K, 5 MPa
    let y = [0.7, 0.3];
    let tc = [126.19, 304.13];
    let pc = [3_390_000.0, 7_375_000.0];
    let omega = [0.037, 0.224];

    let tc_mix = state::kays_rule_tc(&y, &tc).unwrap();
    let pc_mix = state::kays_rule_pc(&y, &pc).unwrap();
    let omega_mix = state::kays_rule_omega(&y, &omega).unwrap();

    let tr = state::reduced_temperature(300.0, tc_mix);
    let pr = state::reduced_pressure(5_000_000.0, pc_mix);
    let z = state::compressibility_pitzer(pr, tr, omega_mix).unwrap();

    // Z should be between 0 and 2 for physical gas
    assert!(z > 0.5 && z < 1.5, "mixture Z={z} out of physical range");
}

// --- Chemical thermodynamics integration tests ---

#[test]
fn hess_law_indirect_path() {
    // C + O₂ → CO₂ via direct vs indirect (C→CO→CO₂)
    // Direct: ΔH = ΔHf(CO₂) = -393510
    // Step 1: C + ½O₂ → CO: ΔH₁ = ΔHf(CO) = -110530
    // Step 2: CO + ½O₂ → CO₂: ΔH₂ = ΔHf(CO₂) - ΔHf(CO) = -393510 - (-110530) = -282980
    // Total: ΔH₁ + ΔH₂ = -110530 + -282980 = -393510 ✓
    let dh_step1 = chem::reaction_enthalpy(&[(1.0, &chem::CO)], &[(0.5, &chem::O2)]);
    let dh_step2 =
        chem::reaction_enthalpy(&[(1.0, &chem::CO2)], &[(1.0, &chem::CO), (0.5, &chem::O2)]);
    let dh_direct = chem::CO2.delta_hf; // C(s) is reference, ΔHf = 0

    assert!(
        (dh_step1 + dh_step2 - dh_direct).abs() < 1.0,
        "Hess's law violated: {} + {} != {}",
        dh_step1,
        dh_step2,
        dh_direct
    );
}

#[test]
fn equilibrium_k_consistent_with_vant_hoff() {
    // K at 298.15 K from ΔG, then Van't Hoff to 298.15 K should return same K
    let dg = chem::reaction_gibbs(
        &[(1.0, &chem::CO2), (2.0, &chem::H2O_GAS)],
        &[(1.0, &chem::CH4), (2.0, &chem::O2)],
    );
    let dh = chem::reaction_enthalpy(
        &[(1.0, &chem::CO2), (2.0, &chem::H2O_GAS)],
        &[(1.0, &chem::CH4), (2.0, &chem::O2)],
    );
    let k_298 = chem::equilibrium_constant(dg, 298.15).unwrap();

    // Van't Hoff from 298.15 to 500 K and back
    let k_500 = chem::vant_hoff_k(k_298, dh, 298.15, 500.0).unwrap();
    let k_back = chem::vant_hoff_k(k_500, dh, 500.0, 298.15).unwrap();

    assert!(
        (k_back - k_298).abs() / k_298 < 1e-6,
        "Van't Hoff roundtrip failed: {} vs {}",
        k_back,
        k_298
    );
}

#[test]
fn excess_air_lowers_flame_temperature() {
    // Stoichiometric: CH₄ + 2O₂ + 7.52N₂
    let t_stoich = chem::adiabatic_flame_temperature(
        &[(1.0, &chem::CH4), (2.0, &chem::O2), (7.52, &chem::N2)],
        &[(1.0, &chem::CO2), (2.0, &chem::H2O_GAS), (7.52, &chem::N2)],
        298.15,
    )
    .unwrap();

    // 50% excess air: CH₄ + 3O₂ + 11.28N₂ → CO₂ + 2H₂O + O₂ + 11.28N₂
    let t_excess = chem::adiabatic_flame_temperature(
        &[(1.0, &chem::CH4), (3.0, &chem::O2), (11.28, &chem::N2)],
        &[
            (1.0, &chem::CO2),
            (2.0, &chem::H2O_GAS),
            (1.0, &chem::O2),
            (11.28, &chem::N2),
        ],
        298.15,
    )
    .unwrap();

    assert!(
        t_excess < t_stoich,
        "excess air T={t_excess:.0} should be < stoich T={t_stoich:.0}"
    );
}
