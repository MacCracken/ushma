//! `transfer.tsv` (M1c onward): every public function of Rust 1.4.0's `transfer`, over seeded physical grids and the
//! edges every argument can take.
//!
//! For each function: a seeded sweep over its arguments' physical ranges, then, for each argument in turn, a physical
//! base case with that argument replaced by each edge value (0, −0, −1, NaN, ±∞, the 1e-30 divisor threshold and its
//! neighbours, 1e-300, 1e300). Outputs are the f64 result or `E:<Variant>`. Classes: B for + − × ÷ √ and `powi`; C
//! where `exp` or `powf` enters (`lumped_capacitance`, `nusselt_dittus_boelter`, `nusselt_natural_vertical`,
//! `thermal_boundary_layer`), compared against the corrected oracle (`transfer.cr.tsv`, ADR 0004).

use ushma::transfer as t;

use crate::grids::Rng;
use crate::tsv::{Class, Module, Tok, err, f};

/// A result as its tokens.
fn res(r: ushma::error::Result<f64>) -> Vec<Tok> {
    match r {
        Ok(x) => vec![f(x)],
        Err(e) => vec![err(&e)],
    }
}

/// The edge values every argument is tried at.
fn edges() -> Vec<f64> {
    let tiny = 1e-30_f64;
    vec![
        0.0,
        -0.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        tiny,
        f64::from_bits(tiny.to_bits() - 1),
        f64::from_bits(tiny.to_bits() + 1),
        -tiny,
        1e-300,
        1e300,
    ]
}

/// The argument tuples for a function of `bases.len()` arguments: `n` seeded draws from `draw`, then every base
/// case with each argument at each edge.
fn tuples(
    bases: &[Vec<f64>],
    n: usize,
    seed: u64,
    draw: &dyn Fn(&mut Rng) -> Vec<f64>,
) -> Vec<Vec<f64>> {
    let mut v: Vec<Vec<f64>> = bases.to_vec();
    let mut r = Rng::new(seed);
    for _ in 0..n {
        v.push(draw(&mut r));
    }
    for base in bases {
        for k in 0..base.len() {
            for e in edges() {
                let mut c = base.clone();
                c[k] = e;
                v.push(c);
            }
        }
    }
    v
}

fn row(m: &mut Module, name: &str, class: Class, args: &[f64], out: Vec<Tok>) {
    let ins: Vec<Tok> = args.iter().map(|&x| f(x)).collect();
    m.row(name, class, &ins, &out);
}

fn part1(m: &mut Module) {
    // conduction(k, A, T_hot, T_cold, L)
    for a in tuples(
        &[
            vec![401.0, 1.0, 400.0, 300.0, 0.1],
            vec![0.0263, 2.5, 293.15, 313.15, 0.005],
        ],
        400,
        1,
        &|r| {
            vec![
                r.log_uniform(-3.0, 3.5),
                r.log_uniform(-4.0, 2.0),
                r.uniform(1.0, 2000.0),
                r.uniform(1.0, 2000.0),
                r.signed_log_uniform(-6.0, 1.0),
            ]
        },
    ) {
        row(
            m,
            "transfer::conduction",
            Class::B,
            &a,
            res(t::conduction(a[0], a[1], a[2], a[3], a[4])),
        );
    }
    // convection(h, A, T_s, T_f)
    for a in tuples(&[vec![25.0, 2.0, 350.0, 300.0]], 400, 2, &|r| {
        vec![
            r.log_uniform(0.0, 5.0),
            r.log_uniform(-4.0, 2.0),
            r.uniform(1.0, 2000.0),
            r.uniform(1.0, 2000.0),
        ]
    }) {
        row(
            m,
            "transfer::convection",
            Class::B,
            &a,
            vec![f(t::convection(a[0], a[1], a[2], a[3]))],
        );
    }
    // radiation(ε, A, T_s, T_surr)
    let mut rad = tuples(&[vec![0.9, 1.0, 500.0, 300.0]], 400, 3, &|r| {
        vec![
            r.unit(),
            r.log_uniform(-4.0, 2.0),
            r.uniform(0.0, 3000.0),
            r.uniform(0.0, 3000.0),
        ]
    });
    rad.extend([
        vec![1.0, 1.0, 300.0, 300.0],
        vec![1.0 + f64::EPSILON, 1.0, 300.0, 0.0],
        vec![-0.0, 1.0, 1e77, 0.0],
    ]);
    for a in rad {
        row(
            m,
            "transfer::radiation",
            Class::B,
            &a,
            res(t::radiation(a[0], a[1], a[2], a[3])),
        );
    }
    // thermal_resistance_conduction(k, A, L)
    for a in tuples(&[vec![401.0, 1.0, 0.1]], 300, 4, &|r| {
        vec![
            r.log_uniform(-3.0, 3.5),
            r.log_uniform(-4.0, 2.0),
            r.log_uniform(-6.0, 1.0),
        ]
    }) {
        row(
            m,
            "transfer::thermal_resistance_conduction",
            Class::B,
            &a,
            res(t::thermal_resistance_conduction(a[0], a[1], a[2])),
        );
    }
    // thermal_resistance_convection(h, A)
    for a in tuples(&[vec![25.0, 2.0]], 300, 5, &|r| {
        vec![r.log_uniform(-1.0, 5.0), r.log_uniform(-4.0, 2.0)]
    }) {
        row(
            m,
            "transfer::thermal_resistance_convection",
            Class::B,
            &a,
            res(t::thermal_resistance_convection(a[0], a[1])),
        );
    }
    // thermal_resistance_series / parallel over lists (`v:n` then the elements)
    let mut lists: Vec<Vec<f64>> = vec![
        vec![],
        vec![1.0],
        vec![0.5, 0.25, 0.125],
        vec![1e300, 1e300],
        vec![-0.0],
        vec![0.0, 0.0],
    ];
    for e in edges() {
        lists.push(vec![2.0, e, 3.0]);
    }
    let mut r = Rng::new(6);
    for _ in 0..200 {
        let n = r.int(1, 8) as usize;
        lists.push((0..n).map(|_| r.log_uniform(-4.0, 4.0)).collect());
    }
    for l in &lists {
        let mut ins = vec![Tok::V(l.len())];
        ins.extend(l.iter().map(|&x| f(x)));
        m.row(
            "transfer::thermal_resistance_series",
            Class::B,
            &ins,
            &[f(t::thermal_resistance_series(l))],
        );
        m.row(
            "transfer::thermal_resistance_parallel",
            Class::B,
            &ins,
            &res(t::thermal_resistance_parallel(l)),
        );
    }
    // heat_stored(m, c, ΔT)
    for a in tuples(&[vec![2.0, 4186.0, 25.0]], 300, 7, &|r| {
        vec![
            r.log_uniform(-3.0, 4.0),
            r.log_uniform(2.0, 4.0),
            r.signed_log_uniform(-2.0, 3.0),
        ]
    }) {
        row(
            m,
            "transfer::heat_stored",
            Class::B,
            &a,
            vec![f(t::heat_stored(a[0], a[1], a[2]))],
        );
    }
    // thermal_diffusivity(k, ρ, c_p)
    for a in tuples(&[vec![401.0, 8960.0, 385.0]], 300, 8, &|r| {
        vec![
            r.log_uniform(-2.0, 3.5),
            r.log_uniform(0.0, 4.5),
            r.log_uniform(2.0, 4.0),
        ]
    }) {
        row(
            m,
            "transfer::thermal_diffusivity",
            Class::B,
            &a,
            res(t::thermal_diffusivity(a[0], a[1], a[2])),
        );
    }
    // biot_number(h, L, k)
    for a in tuples(&[vec![25.0, 0.01, 50.0]], 300, 9, &|r| {
        vec![
            r.log_uniform(0.0, 4.0),
            r.log_uniform(-4.0, 0.0),
            r.log_uniform(-2.0, 3.5),
        ]
    }) {
        row(
            m,
            "transfer::biot_number",
            Class::B,
            &a,
            res(t::biot_number(a[0], a[1], a[2])),
        );
    }
    // lumped_capacitance(T0, T_env, t, τ)
    for a in tuples(&[vec![500.0, 300.0, 60.0, 120.0]], 400, 10, &|r| {
        vec![
            r.uniform(1.0, 2000.0),
            r.uniform(1.0, 2000.0),
            r.log_uniform(-2.0, 5.0),
            r.log_uniform(-1.0, 4.0),
        ]
    }) {
        row(
            m,
            "transfer::lumped_capacitance",
            Class::C,
            &a,
            vec![f(t::lumped_capacitance(a[0], a[1], a[2], a[3]))],
        );
    }
    // reynolds_number(V, L, ν)
    for a in tuples(&[vec![2.0, 0.05, 1.5e-5]], 300, 11, &|r| {
        vec![
            r.log_uniform(-3.0, 2.0),
            r.log_uniform(-4.0, 1.0),
            r.log_uniform(-7.0, -3.0),
        ]
    }) {
        row(
            m,
            "transfer::reynolds_number",
            Class::B,
            &a,
            res(t::reynolds_number(a[0], a[1], a[2])),
        );
    }
    // prandtl_number(ν, α)
    for a in tuples(&[vec![1.5e-5, 2.2e-5]], 300, 12, &|r| {
        vec![r.log_uniform(-7.0, -3.0), r.log_uniform(-8.0, -3.0)]
    }) {
        row(
            m,
            "transfer::prandtl_number",
            Class::B,
            &a,
            res(t::prandtl_number(a[0], a[1])),
        );
    }
    // nusselt_number(h, L, k)
    for a in tuples(&[vec![25.0, 0.1, 0.026]], 300, 13, &|r| {
        vec![
            r.log_uniform(0.0, 4.0),
            r.log_uniform(-4.0, 1.0),
            r.log_uniform(-2.0, 3.0),
        ]
    }) {
        row(
            m,
            "transfer::nusselt_number",
            Class::B,
            &a,
            res(t::nusselt_number(a[0], a[1], a[2])),
        );
    }
    // nusselt_dittus_boelter(Re, Pr)
    for a in tuples(&[vec![50_000.0, 0.7]], 600, 14, &|r| {
        vec![r.log_uniform(3.0, 7.0), r.log_uniform(-1.0, 3.0)]
    }) {
        row(
            m,
            "transfer::nusselt_dittus_boelter",
            Class::C,
            &a,
            res(t::nusselt_dittus_boelter(a[0], a[1])),
        );
    }
    // nusselt_natural_vertical(Ra, Pr)
    for a in tuples(&[vec![1e9, 0.71], vec![0.0, 7.0]], 600, 15, &|r| {
        vec![r.log_uniform(-2.0, 14.0), r.log_uniform(-2.0, 3.0)]
    }) {
        row(
            m,
            "transfer::nusselt_natural_vertical",
            Class::C,
            &a,
            res(t::nusselt_natural_vertical(a[0], a[1])),
        );
    }
    // boundary_layer_thickness(x, Re_x)
    for a in tuples(&[vec![0.5, 1e5]], 300, 16, &|r| {
        vec![r.log_uniform(-3.0, 1.0), r.log_uniform(1.0, 6.0)]
    }) {
        row(
            m,
            "transfer::boundary_layer_thickness",
            Class::B,
            &a,
            res(t::boundary_layer_thickness(a[0], a[1])),
        );
    }
    // thermal_boundary_layer(x, Re_x, Pr)
    for a in tuples(&[vec![0.5, 1e5, 0.7]], 600, 17, &|r| {
        vec![
            r.log_uniform(-3.0, 1.0),
            r.log_uniform(1.0, 6.0),
            r.log_uniform(-2.0, 3.0),
        ]
    }) {
        row(
            m,
            "transfer::thermal_boundary_layer",
            Class::C,
            &a,
            res(t::thermal_boundary_layer(a[0], a[1], a[2])),
        );
    }
}

pub fn emit() -> Module {
    let mut m = Module::new("transfer");
    part1(&mut m);
    m
}
