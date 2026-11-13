//! `bridge.tsv` (M7a): every public function of Rust 1.4.0's `bridge`.
//!
//! Each conversion seeded over its physical ranges, then each argument of its base cases (the 1.4.0 tests' values
//! among them) at every edge value, and each fallback's switch either side: m·c_p at 0 and its sign, K = 1 ± 12 ulps
//! across equilibrium_temperature's |ln K| < 1e-15 test, the altitude where T reaches 0, η at 1 ± 2 ulps, the
//! activation at 0 and 1, and the wet-bulb clamp near saturation over a grid of temperatures. Classes: B for + − × ÷ √
//! (with powi(4)), C where ln, powf or atan enters (`equilibrium_temperature`, `humidity_to_wet_bulb`,
//! `velocity_to_convection_coefficient`, `basal_metabolic_rate`, `body_surface_area`).

use ushma::bridge as br;

use crate::grids::Rng;
use crate::sweep::{around, row, tuples};
use crate::tsv::{Class, Module, f};

/// One row per argument tuple, of a function giving one f64.
fn rows(m: &mut Module, name: &str, class: Class, args: &[Vec<f64>], g: &dyn Fn(&[f64]) -> f64) {
    for a in args {
        row(m, name, class, a, vec![f(g(a))]);
    }
}

fn bijli(m: &mut Module) {
    // joule_heating_to_temperature_rate(P, m, c_p), with m·c_p at 0, negative, under- and overflowing
    let mut v = tuples(
        &[vec![100.0, 1.0, 4186.0], vec![100.0, 1.0, 385.0]],
        300,
        700,
        &|r| {
            vec![
                r.signed_log_uniform(-3.0, 6.0),
                r.log_uniform(-3.0, 3.0),
                r.log_uniform(2.0, 4.0),
            ]
        },
    );
    v.extend([
        vec![100.0, -1.0, -4186.0],
        vec![100.0, -1.0, 4186.0],
        vec![100.0, 1e-200, 1e-200],
        vec![100.0, 1e300, 1e300],
        vec![f64::INFINITY, 1.0, 1.0],
        vec![100.0, f64::from_bits(1), 1.0],
    ]);
    rows(m, "bridge::joule_heating_to_temperature_rate", Class::B, &v, &|a| {
        br::joule_heating_to_temperature_rate(a[0], a[1], a[2])
    });
    // em_absorption_to_heat_source(q)
    let v = tuples(&[vec![500.0], vec![-100.0]], 200, 701, &|r| {
        vec![r.signed_log_uniform(-3.0, 9.0)]
    });
    rows(m, "bridge::em_absorption_to_heat_source", Class::B, &v, &|a| {
        br::em_absorption_to_heat_source(a[0])
    });
}

fn kimiya(m: &mut Module) {
    // reaction_heat_rate(ΔH, ṅ)
    let v = tuples(
        &[vec![-100_000.0, 0.1], vec![50_000.0, 0.01]],
        300,
        702,
        &|r| vec![r.signed_log_uniform(1.0, 7.0), r.signed_log_uniform(-6.0, 2.0)],
    );
    rows(m, "bridge::reaction_heat_rate", Class::B, &v, &|a| {
        br::reaction_heat_rate(a[0], a[1])
    });
    // equilibrium_temperature(ΔG°, K), with K across 1 (|ln K| < 1e-15) and on either side of it
    let mut v = tuples(
        &[vec![-10_000.0, 100.0], vec![10_000.0, 0.01]],
        400,
        703,
        &|r| vec![r.signed_log_uniform(1.0, 6.0), r.log_uniform(-10.0, 10.0)],
    );
    for k in around(1.0, 12) {
        v.push(vec![-10_000.0, k]);
        v.push(vec![10_000.0, k]);
    }
    v.extend([
        vec![0.0, 100.0],
        vec![-0.0, 100.0],
        vec![-10_000.0, f64::MAX],
        vec![-10_000.0, f64::from_bits(1)],
    ]);
    rows(m, "bridge::equilibrium_temperature", Class::C, &v, &|a| {
        br::equilibrium_temperature(a[0], a[1])
    });
}

fn dravya(m: &mut Module) {
    // temperature_to_thermal_strain(T, α, T_ref)
    let v = tuples(&[vec![393.15, 12e-6, 293.15]], 300, 704, &|r| {
        vec![
            r.uniform(1.0, 2000.0),
            r.log_uniform(-7.0, -4.0),
            r.uniform(250.0, 350.0),
        ]
    });
    rows(m, "bridge::temperature_to_thermal_strain", Class::B, &v, &|a| {
        br::temperature_to_thermal_strain(a[0], a[1], a[2])
    });
    // thermal_gradient_to_stress(dT/dx, E, α, L)
    let v = tuples(&[vec![10.0, 200e9, 12e-6, 1.0]], 300, 705, &|r| {
        vec![
            r.signed_log_uniform(-1.0, 4.0),
            r.log_uniform(9.0, 12.0),
            r.log_uniform(-7.0, -4.0),
            r.log_uniform(-3.0, 1.0),
        ]
    });
    rows(m, "bridge::thermal_gradient_to_stress", Class::B, &v, &|a| {
        br::thermal_gradient_to_stress(a[0], a[1], a[2], a[3])
    });
}

fn badal(m: &mut Module) {
    // altitude_to_temperature(T_surface, h, Γ), with T reaching 0
    let mut v = tuples(
        &[vec![288.15, 1000.0, 0.0065], vec![288.15, 0.0, 0.0065]],
        300,
        706,
        &|r| {
            vec![
                r.uniform(200.0, 320.0),
                r.uniform(-500.0, 20_000.0),
                r.uniform(0.002, 0.01),
            ]
        },
    );
    for h in around(288.15 / 0.0065, 2) {
        v.push(vec![288.15, h, 0.0065]);
    }
    v.extend([
        vec![0.0, 0.0, 0.0065],
        vec![-0.0, 0.0, 0.0065],
        vec![288.15, 1e6, 0.0065],
    ]);
    rows(m, "bridge::altitude_to_temperature", Class::B, &v, &|a| {
        br::altitude_to_temperature(a[0], a[1], a[2])
    });
    // humidity_to_wet_bulb(T, RH): the fitted range, an extrapolation, and the clamp near saturation
    let bases = [
        vec![300.0, 0.2],
        vec![300.0, 1.0],
        vec![310.0, 0.99],
        vec![273.15, 1.0],
        vec![293.15, 0.5],
    ];
    let mut v = tuples(&bases, 500, 707, &|r| {
        vec![r.uniform(253.15, 323.15), r.uniform(0.05, 0.99)]
    });
    let mut r = Rng::new(708);
    for _ in 0..200 {
        v.push(vec![r.uniform(150.0, 400.0), r.uniform(-0.2, 1.3)]);
    }
    for k in 0..=12 {
        let t = 253.15 + 6.0 * f64::from(k);
        for rh in [0.9, 0.95, 0.97, 0.98, 0.99, 0.995, 1.0, 1.01] {
            v.push(vec![t, rh]);
        }
    }
    rows(m, "bridge::humidity_to_wet_bulb", Class::C, &v, &|a| {
        br::humidity_to_wet_bulb(a[0], a[1])
    });
}

fn pravash(m: &mut Module) {
    // velocity_to_convection_coefficient(v, ρ, μ, k, c_p, D)
    let bases = [
        vec![1.0, 1000.0, 0.001, 0.6, 4186.0, 0.1],
        vec![1.0, 998.0, 1e-3, 0.6, 4180.0, 0.05],
        vec![0.0, 1000.0, 0.001, 0.6, 4186.0, 0.1],
        vec![-2.0, 1.2, 1.8e-5, 0.026, 1005.0, 0.2],
    ];
    let v = tuples(&bases, 400, 709, &|r| {
        vec![
            r.signed_log_uniform(-3.0, 2.0),
            r.log_uniform(-1.0, 4.0),
            r.log_uniform(-6.0, 0.0),
            r.log_uniform(-2.0, 2.5),
            r.log_uniform(2.5, 4.2),
            r.log_uniform(-3.0, 0.0),
        ]
    });
    rows(m, "bridge::velocity_to_convection_coefficient", Class::C, &v, &|a| {
        br::velocity_to_convection_coefficient(a[0], a[1], a[2], a[3], a[4], a[5])
    });
    // tke_to_eddy_thermal_diffusivity(k, L)
    let v = tuples(&[vec![1.0, 0.1]], 300, 710, &|r| {
        vec![r.log_uniform(-6.0, 3.0), r.log_uniform(-4.0, 1.0)]
    });
    rows(m, "bridge::tke_to_eddy_thermal_diffusivity", Class::B, &v, &|a| {
        br::tke_to_eddy_thermal_diffusivity(a[0], a[1])
    });
}

fn sharira(m: &mut Module) {
    // muscle_work_to_heat(P, η), with η at 0, at 1 ± 2 ulps, above 1 and NaN
    let bases = [
        vec![100.0, 0.25],
        vec![100.0, 0.0],
        vec![100.0, 1.0],
        vec![100.0, 1.5],
        vec![-100.0, 0.25],
    ];
    let mut v = tuples(&bases, 300, 711, &|r| {
        vec![r.signed_log_uniform(0.0, 3.0), r.uniform(0.0, 1.2)]
    });
    for e in around(1.0, 2) {
        v.push(vec![100.0, e]);
    }
    v.push(vec![100.0, f64::from_bits(1)]);
    rows(m, "bridge::muscle_work_to_heat", Class::B, &v, &|a| {
        br::muscle_work_to_heat(a[0], a[1])
    });
    // basal_metabolic_rate(m)
    let v = tuples(&[vec![70.0]], 300, 712, &|r| vec![r.log_uniform(-3.0, 5.0)]);
    rows(m, "bridge::basal_metabolic_rate", Class::C, &v, &|a| {
        br::basal_metabolic_rate(a[0])
    });
    // body_surface_area(m, h)
    let v = tuples(&[vec![70.0, 1.75]], 300, 713, &|r| {
        vec![r.uniform(1.0, 200.0), r.uniform(0.3, 2.3)]
    });
    rows(m, "bridge::body_surface_area", Class::C, &v, &|a| {
        br::body_surface_area(a[0], a[1])
    });
    // skin_radiation_loss(A, T_skin, T_env), T⁴ overflowing among the edges
    let mut v = tuples(&[vec![1.85, 310.0, 293.0]], 300, 714, &|r| {
        vec![
            r.uniform(0.1, 3.0),
            r.uniform(290.0, 315.0),
            r.uniform(200.0, 330.0),
        ]
    });
    v.extend([vec![1.85, 310.0, 310.0], vec![1.0, 1e100, 1.0]]);
    rows(m, "bridge::skin_radiation_loss", Class::B, &v, &|a| {
        br::skin_radiation_loss(a[0], a[1], a[2])
    });
    // muscle_activation_heat(a, F), with a at and across 0 and 1
    let mut v = tuples(&[vec![0.5, 1000.0]], 300, 715, &|r| {
        vec![r.uniform(-0.2, 1.2), r.signed_log_uniform(0.0, 4.0)]
    });
    for a in around(1.0, 2) {
        v.push(vec![a, 1000.0]);
    }
    v.push(vec![f64::from_bits(1), 1000.0]);
    rows(m, "bridge::muscle_activation_heat", Class::B, &v, &|a| {
        br::muscle_activation_heat(a[0], a[1])
    });
    // evaporative_cooling(ṁ)
    let v = tuples(&[vec![0.5 / 3600.0], vec![-0.001]], 200, 716, &|r| {
        vec![r.signed_log_uniform(-6.0, -2.0)]
    });
    rows(m, "bridge::evaporative_cooling", Class::B, &v, &|a| {
        br::evaporative_cooling(a[0])
    });
}

pub fn emit() -> Module {
    let mut m = Module::new("bridge");
    bijli(&mut m);
    kimiya(&mut m);
    dravya(&mut m);
    badal(&mut m);
    pravash(&mut m);
    sharira(&mut m);
    m
}
