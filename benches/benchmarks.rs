//! Criterion benchmarks for ushma thermodynamics.

use criterion::{Criterion, black_box, criterion_group, criterion_main};

use ushma::entropy;
use ushma::material;
use ushma::state;
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
        b.iter(|| {
            transfer::heat_stored(black_box(1.0), black_box(4186.0), black_box(10.0))
        });
    });
}

fn bench_thermal_diffusivity(c: &mut Criterion) {
    c.bench_function("transfer/diffusivity", |b| {
        b.iter(|| {
            transfer::thermal_diffusivity(
                black_box(401.0),
                black_box(8960.0),
                black_box(385.0),
            )
        });
    });
}

fn bench_biot_number(c: &mut Criterion) {
    c.bench_function("transfer/biot_number", |b| {
        b.iter(|| {
            transfer::biot_number(black_box(25.0), black_box(0.01), black_box(401.0))
        });
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

// --- State benchmarks ---

fn bench_ideal_gas_pressure(c: &mut Criterion) {
    c.bench_function("state/ideal_gas_pressure", |b| {
        b.iter(|| {
            state::ideal_gas_pressure(black_box(1.0), black_box(300.0), black_box(0.025))
        });
    });
}

fn bench_ideal_gas_volume(c: &mut Criterion) {
    c.bench_function("state/ideal_gas_volume", |b| {
        b.iter(|| {
            state::ideal_gas_volume(black_box(1.0), black_box(300.0), black_box(101325.0))
        });
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
        b.iter(|| {
            entropy::helmholtz(black_box(1000.0), black_box(300.0), black_box(5.0))
        });
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
);

criterion_main!(benches);
