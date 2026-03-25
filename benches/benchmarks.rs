//! Criterion benchmarks for ushma thermodynamics.

use criterion::{Criterion, black_box, criterion_group, criterion_main};

use ushma::cycle;
use ushma::entropy;
use ushma::material;
use ushma::numerical;
use ushma::phase;
use ushma::state;
use ushma::steam;
use ushma::transfer;

// --- Transfer benchmarks ---

fn bench_conduction(c: &mut Criterion) {
    c.bench_function("transfer/conduction", |b| {
        b.iter(|| {
            transfer::conduction(
                black_box(401.0),
                black_box(0.01),
                black_box(373.15),
                black_box(293.15),
                black_box(0.1),
            )
        });
    });
}

fn bench_convection(c: &mut Criterion) {
    c.bench_function("transfer/convection", |b| {
        b.iter(|| {
            transfer::convection(
                black_box(25.0),
                black_box(1.0),
                black_box(373.15),
                black_box(293.15),
            )
        });
    });
}

fn bench_radiation(c: &mut Criterion) {
    c.bench_function("transfer/radiation", |b| {
        b.iter(|| {
            transfer::radiation(
                black_box(0.9),
                black_box(1.0),
                black_box(473.15),
                black_box(293.15),
            )
        });
    });
}

fn bench_thermal_resistance_series(c: &mut Criterion) {
    let resistances = vec![0.5, 1.0, 2.0, 0.3, 1.5];
    c.bench_function("transfer/resistance_series", |b| {
        b.iter(|| transfer::thermal_resistance_series(black_box(&resistances)));
    });
}

fn bench_thermal_resistance_parallel(c: &mut Criterion) {
    let resistances = vec![2.0, 3.0, 6.0];
    c.bench_function("transfer/resistance_parallel", |b| {
        b.iter(|| transfer::thermal_resistance_parallel(black_box(&resistances)));
    });
}

fn bench_heat_stored(c: &mut Criterion) {
    c.bench_function("transfer/heat_stored", |b| {
        b.iter(|| transfer::heat_stored(black_box(1.0), black_box(4186.0), black_box(10.0)));
    });
}

fn bench_thermal_diffusivity(c: &mut Criterion) {
    c.bench_function("transfer/diffusivity", |b| {
        b.iter(|| {
            transfer::thermal_diffusivity(black_box(401.0), black_box(8960.0), black_box(385.0))
        });
    });
}

fn bench_biot_number(c: &mut Criterion) {
    c.bench_function("transfer/biot_number", |b| {
        b.iter(|| transfer::biot_number(black_box(25.0), black_box(0.01), black_box(401.0)));
    });
}

fn bench_lumped_capacitance(c: &mut Criterion) {
    c.bench_function("transfer/lumped_capacitance", |b| {
        b.iter(|| {
            transfer::lumped_capacitance(
                black_box(373.15),
                black_box(293.15),
                black_box(60.0),
                black_box(100.0),
            )
        });
    });
}

fn bench_reynolds_number(c: &mut Criterion) {
    c.bench_function("transfer/reynolds", |b| {
        b.iter(|| transfer::reynolds_number(black_box(1.0), black_box(0.1), black_box(1e-6)));
    });
}

fn bench_nusselt_dittus_boelter(c: &mut Criterion) {
    c.bench_function("transfer/nusselt_db", |b| {
        b.iter(|| transfer::nusselt_dittus_boelter(black_box(50_000.0), black_box(0.7)));
    });
}

fn bench_fin_rectangular_heat(c: &mut Criterion) {
    c.bench_function("transfer/fin_rect_heat", |b| {
        b.iter(|| {
            transfer::fin_rectangular_heat(
                black_box(25.0),
                black_box(0.02),
                black_box(237.0),
                black_box(1e-4),
                black_box(0.05),
                black_box(373.15),
                black_box(293.15),
            )
        });
    });
}

fn bench_lmtd_counter(c: &mut Criterion) {
    c.bench_function("transfer/lmtd_counter", |b| {
        b.iter(|| {
            transfer::lmtd_counter(
                black_box(400.0),
                black_box(350.0),
                black_box(300.0),
                black_box(330.0),
            )
        });
    });
}

fn bench_effectiveness_counter(c: &mut Criterion) {
    c.bench_function("transfer/eff_counter", |b| {
        b.iter(|| transfer::effectiveness_counter(black_box(2.0), black_box(0.5)));
    });
}

fn bench_view_factor_disks(c: &mut Criterion) {
    c.bench_function("transfer/view_factor_disks", |b| {
        b.iter(|| {
            transfer::view_factor_coaxial_disks(black_box(0.5), black_box(0.5), black_box(1.0))
        });
    });
}

// --- State benchmarks ---

fn bench_ideal_gas_pressure(c: &mut Criterion) {
    c.bench_function("state/ideal_gas_pressure", |b| {
        b.iter(|| state::ideal_gas_pressure(black_box(1.0), black_box(300.0), black_box(0.025)));
    });
}

fn bench_ideal_gas_volume(c: &mut Criterion) {
    c.bench_function("state/ideal_gas_volume", |b| {
        b.iter(|| state::ideal_gas_volume(black_box(1.0), black_box(300.0), black_box(101325.0)));
    });
}

fn bench_van_der_waals(c: &mut Criterion) {
    c.bench_function("state/van_der_waals", |b| {
        b.iter(|| {
            state::van_der_waals_pressure(
                black_box(1.0),
                black_box(300.0),
                black_box(0.025),
                black_box(0.3658),
                black_box(4.286e-5),
            )
        });
    });
}

fn bench_isothermal_work(c: &mut Criterion) {
    c.bench_function("state/isothermal_work", |b| {
        b.iter(|| {
            state::isothermal_work(
                black_box(1.0),
                black_box(300.0),
                black_box(0.01),
                black_box(0.02),
            )
        });
    });
}

fn bench_adiabatic_temperature(c: &mut Criterion) {
    c.bench_function("state/adiabatic_temperature", |b| {
        b.iter(|| {
            state::adiabatic_temperature(
                black_box(300.0),
                black_box(0.02),
                black_box(0.01),
                black_box(1.4),
            )
        });
    });
}

fn bench_compressibility_factor(c: &mut Criterion) {
    c.bench_function("state/compressibility_factor", |b| {
        b.iter(|| {
            state::compressibility_factor(
                black_box(101325.0),
                black_box(0.02241),
                black_box(1.0),
                black_box(273.15),
            )
        });
    });
}

// --- Entropy benchmarks ---

fn bench_carnot_efficiency(c: &mut Criterion) {
    c.bench_function("entropy/carnot_efficiency", |b| {
        b.iter(|| entropy::carnot_efficiency(black_box(500.0), black_box(300.0)));
    });
}

fn bench_carnot_cop(c: &mut Criterion) {
    c.bench_function("entropy/carnot_cop", |b| {
        b.iter(|| entropy::carnot_cop_refrigeration(black_box(300.0), black_box(250.0)));
    });
}

fn bench_helmholtz(c: &mut Criterion) {
    c.bench_function("entropy/helmholtz", |b| {
        b.iter(|| entropy::helmholtz(black_box(1000.0), black_box(300.0), black_box(5.0)));
    });
}

fn bench_gibbs(c: &mut Criterion) {
    c.bench_function("entropy/gibbs", |b| {
        b.iter(|| entropy::gibbs(black_box(2000.0), black_box(300.0), black_box(5.0)));
    });
}

fn bench_entropy_of_mixing(c: &mut Criterion) {
    let fractions = vec![0.3, 0.3, 0.4];
    c.bench_function("entropy/mixing", |b| {
        b.iter(|| entropy::entropy_of_mixing(black_box(1.0), black_box(&fractions)));
    });
}

fn bench_ideal_gas_entropy_change(c: &mut Criterion) {
    c.bench_function("entropy/ideal_gas_ds", |b| {
        b.iter(|| {
            entropy::ideal_gas_entropy_change(
                black_box(1.0),
                black_box(20.8),
                black_box(300.0),
                black_box(400.0),
                black_box(0.01),
                black_box(0.02),
            )
        });
    });
}

// --- Material benchmarks ---

fn bench_material_diffusivity(c: &mut Criterion) {
    c.bench_function("material/diffusivity", |b| {
        b.iter(|| material::COPPER.diffusivity());
    });
}

fn bench_material_volumetric_heat_capacity(c: &mut Criterion) {
    c.bench_function("material/volumetric_cp", |b| {
        b.iter(|| material::WATER.volumetric_heat_capacity());
    });
}

// --- Phase benchmarks ---

fn bench_clausius_clapeyron_slope(c: &mut Criterion) {
    c.bench_function("phase/clausius_clapeyron_slope", |b| {
        b.iter(|| {
            phase::clausius_clapeyron_slope(
                black_box(2_260_000.0),
                black_box(373.15),
                black_box(1.672),
            )
        });
    });
}

fn bench_clausius_clapeyron_pressure(c: &mut Criterion) {
    c.bench_function("phase/clausius_clapeyron_pressure", |b| {
        b.iter(|| {
            phase::clausius_clapeyron_pressure(
                black_box(101_325.0),
                black_box(373.15),
                black_box(383.15),
                black_box(40_714.0),
            )
        });
    });
}

fn bench_heat_of_fusion(c: &mut Criterion) {
    c.bench_function("phase/heat_of_fusion", |b| {
        b.iter(|| phase::heat_of_fusion(black_box(&phase::WATER_PHASE), black_box(1.0)));
    });
}

fn bench_heat_of_vaporization(c: &mut Criterion) {
    c.bench_function("phase/heat_of_vaporization", |b| {
        b.iter(|| phase::heat_of_vaporization(black_box(&phase::WATER_PHASE), black_box(1.0)));
    });
}

fn bench_heat_for_phase_change(c: &mut Criterion) {
    c.bench_function("phase/heat_for_phase_change", |b| {
        b.iter(|| {
            phase::heat_for_phase_change(
                black_box(1.0),
                black_box(2500.0),
                black_box(263.0),
                black_box(383.0),
                black_box(&phase::WATER_PHASE),
            )
        });
    });
}

fn bench_saturated_by_temp(c: &mut Criterion) {
    c.bench_function("steam/saturated_by_temp", |b| {
        b.iter(|| steam::saturated_by_temperature(black_box(373.15)));
    });
}

fn bench_saturated_by_pressure(c: &mut Criterion) {
    c.bench_function("steam/saturated_by_pressure", |b| {
        b.iter(|| steam::saturated_by_pressure(black_box(101_325.0)));
    });
}

fn bench_quality_from_enthalpy(c: &mut Criterion) {
    let entry = steam::saturated_by_temperature(373.15).unwrap();
    let h_mid = entry.h_f + 0.5 * entry.h_fg;
    c.bench_function("steam/quality_from_enthalpy", |b| {
        b.iter(|| steam::quality_from_enthalpy(black_box(h_mid), black_box(&entry)));
    });
}

fn bench_wet_steam_properties(c: &mut Criterion) {
    let entry = steam::saturated_by_temperature(373.15).unwrap();
    c.bench_function("steam/wet_steam_properties", |b| {
        b.iter(|| steam::wet_steam_properties(black_box(0.75), black_box(&entry)));
    });
}

fn bench_otto_cycle(c: &mut Criterion) {
    c.bench_function("cycle/otto", |b| {
        b.iter(|| {
            cycle::otto_cycle(
                black_box(300.0),
                black_box(101_325.0),
                black_box(8.0),
                black_box(50_000.0),
                black_box(1.4),
                black_box(1.0),
            )
        });
    });
}

fn bench_diesel_cycle(c: &mut Criterion) {
    c.bench_function("cycle/diesel", |b| {
        b.iter(|| {
            cycle::diesel_cycle(
                black_box(300.0),
                black_box(101_325.0),
                black_box(20.0),
                black_box(2.0),
                black_box(1.4),
                black_box(1.0),
            )
        });
    });
}

fn bench_brayton_cycle(c: &mut Criterion) {
    c.bench_function("cycle/brayton", |b| {
        b.iter(|| {
            cycle::brayton_cycle(
                black_box(300.0),
                black_box(101_325.0),
                black_box(10.0),
                black_box(1400.0),
                black_box(1.4),
                black_box(1.0),
            )
        });
    });
}

fn bench_rankine_cycle(c: &mut Criterion) {
    c.bench_function("cycle/rankine", |b| {
        b.iter(|| {
            cycle::rankine_cycle(
                black_box(10_000.0),
                black_box(2_000_000.0),
                black_box(Some(573.15)),
            )
        });
    });
}

fn bench_refrigeration_cycle(c: &mut Criterion) {
    c.bench_function("cycle/refrigeration", |b| {
        b.iter(|| cycle::refrigeration_cycle(black_box(7_384.0), black_box(47_390.0)));
    });
}

fn bench_superheated_lookup(c: &mut Criterion) {
    c.bench_function("steam/superheated_lookup", |b| {
        b.iter(|| steam::superheated_lookup(black_box(573.15), black_box(500_000.0)));
    });
}

fn bench_phase_lookup(c: &mut Criterion) {
    c.bench_function("phase/phase_lookup", |b| {
        b.iter(|| phase::WATER_PHASE.phase_at(black_box(300.0), black_box(101_325.0)));
    });
}

criterion_group!(
    benches,
    // transfer
    bench_conduction,
    bench_convection,
    bench_radiation,
    bench_thermal_resistance_series,
    bench_thermal_resistance_parallel,
    bench_heat_stored,
    bench_thermal_diffusivity,
    bench_biot_number,
    bench_lumped_capacitance,
    bench_reynolds_number,
    bench_nusselt_dittus_boelter,
    bench_fin_rectangular_heat,
    bench_lmtd_counter,
    bench_effectiveness_counter,
    bench_view_factor_disks,
    // state
    bench_ideal_gas_pressure,
    bench_ideal_gas_volume,
    bench_van_der_waals,
    bench_isothermal_work,
    bench_adiabatic_temperature,
    bench_compressibility_factor,
    // entropy
    bench_carnot_efficiency,
    bench_carnot_cop,
    bench_helmholtz,
    bench_gibbs,
    bench_entropy_of_mixing,
    bench_ideal_gas_entropy_change,
    // material
    bench_material_diffusivity,
    bench_material_volumetric_heat_capacity,
    // phase
    bench_clausius_clapeyron_slope,
    bench_clausius_clapeyron_pressure,
    bench_phase_lookup,
    bench_heat_of_fusion,
    bench_heat_of_vaporization,
    bench_heat_for_phase_change,
    // cycle
    bench_otto_cycle,
    bench_diesel_cycle,
    bench_brayton_cycle,
    bench_rankine_cycle,
    bench_refrigeration_cycle,
    // steam
    bench_saturated_by_temp,
    bench_saturated_by_pressure,
    bench_quality_from_enthalpy,
    bench_wet_steam_properties,
    bench_superheated_lookup,
    // numerical
    bench_explicit_1d_step,
    bench_crank_nicolson_1d_step,
    bench_gauss_seidel_2d,
    bench_thermal_network,
);

fn bench_explicit_1d_step(c: &mut Criterion) {
    let mut g = numerical::ThermalGrid1D::new(
        100,
        1.0,
        1e-4,
        300.0,
        numerical::BoundaryCondition::Fixed(400.0),
        numerical::BoundaryCondition::Fixed(300.0),
    )
    .unwrap();
    let dt = 0.4 * g.dx * g.dx / g.alpha;
    c.bench_function("numerical/explicit_1d", |b| {
        b.iter(|| g.step_explicit(black_box(dt)));
    });
}

fn bench_crank_nicolson_1d_step(c: &mut Criterion) {
    let mut g = numerical::ThermalGrid1D::new(
        100,
        1.0,
        1e-4,
        300.0,
        numerical::BoundaryCondition::Fixed(400.0),
        numerical::BoundaryCondition::Fixed(300.0),
    )
    .unwrap();
    c.bench_function("numerical/crank_nicolson_1d", |b| {
        b.iter(|| g.step_crank_nicolson(black_box(1.0)));
    });
}

fn bench_gauss_seidel_2d(c: &mut Criterion) {
    let mut g = numerical::ThermalGrid2D::new(20, 20, 1.0, 1.0, 300.0).unwrap();
    g.set_boundary(
        numerical::Side::Left,
        numerical::BoundaryCondition::Fixed(400.0),
    );
    g.set_boundary(
        numerical::Side::Right,
        numerical::BoundaryCondition::Fixed(300.0),
    );
    c.bench_function("numerical/gauss_seidel_2d_20x20", |b| {
        b.iter(|| {
            let mut g2 = g.clone();
            g2.solve_steady_state(black_box(1e-4), black_box(10_000))
        });
    });
}

fn bench_thermal_network(c: &mut Criterion) {
    let mut net = numerical::ThermalNetwork::new(5);
    net.add_resistance(0, 1, 1.0).unwrap();
    net.add_resistance(1, 2, 1.0).unwrap();
    net.add_resistance(2, 3, 1.0).unwrap();
    net.add_resistance(3, 4, 1.0).unwrap();
    net.set_fixed_temperature(0, 400.0).unwrap();
    net.set_fixed_temperature(4, 300.0).unwrap();
    c.bench_function("numerical/thermal_network_5", |b| {
        b.iter(|| net.solve());
    });
}

criterion_main!(benches);
