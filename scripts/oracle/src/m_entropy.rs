//! `entropy.tsv` (M2c): every public function of Rust 1.4.0's `entropy`.
//!
//! For each function: a seeded sweep over its arguments' physical ranges and each argument of each base case at every
//! edge value (`sweep::edges`), then its switches: T₁ = T₂ and V₁ = V₂ ± 2 ulps for the entropy changes; the Carnot
//! limits with T_cold at T_hot, within 2 ulps of it and above it. Mixing runs lists of fractions: zero, −0, NaN,
//! negative and out-of-range fractions, sums across the 1e-10 tolerance, subnormal fractions and seeded mixtures.
//! Classes: B for + − × ÷, A for `is_spontaneous`'s bool, C where `ln` enters (`ideal_gas_entropy_change`,
//! `isothermal_entropy_change`, `entropy_of_mixing`).

use ushma::entropy as en;

use crate::grids::Rng;
use crate::sweep::{around, res, row, tuples};
use crate::tsv::{Class, Module, Tok, b, f};

fn changes(m: &mut Module) {
    // ideal_gas_entropy_change(n, C_v, T₁, T₂, V₁, V₂), with T₁ = T₂ and V₁ = V₂
    let mut ig = tuples(
        &[vec![1.0, 20.8, 300.0, 400.0, 0.01, 0.02]],
        500,
        60,
        &|r| {
            vec![
                r.log_uniform(-2.0, 2.0),
                r.uniform(10.0, 60.0),
                r.uniform(1.0, 3000.0),
                r.uniform(1.0, 3000.0),
                r.log_uniform(-5.0, 1.0),
                r.log_uniform(-5.0, 1.0),
            ]
        },
    );
    for t2 in around(300.0, 2) {
        ig.push(vec![1.0, 20.8, 300.0, t2, 0.01, 0.01]);
    }
    for v2 in around(0.01, 2) {
        ig.push(vec![1.0, 20.8, 300.0, 300.0, 0.01, v2]);
    }
    for a in ig {
        row(
            m,
            "entropy::ideal_gas_entropy_change",
            Class::C,
            &a,
            res(en::ideal_gas_entropy_change(
                a[0], a[1], a[2], a[3], a[4], a[5],
            )),
        );
    }
    // isothermal_entropy_change(n, V₁, V₂)
    let mut iso = tuples(&[vec![1.0, 0.01, 0.02]], 400, 61, &|r| {
        vec![
            r.log_uniform(-2.0, 2.0),
            r.log_uniform(-5.0, 1.0),
            r.log_uniform(-5.0, 1.0),
        ]
    });
    for v2 in around(0.01, 2) {
        iso.push(vec![1.0, 0.01, v2]);
    }
    for a in iso {
        row(
            m,
            "entropy::isothermal_entropy_change",
            Class::C,
            &a,
            res(en::isothermal_entropy_change(a[0], a[1], a[2])),
        );
    }
    // heat_transfer_entropy(Q, T)
    for a in tuples(
        &[vec![1000.0, 300.0], vec![-1000.0, 300.0]],
        300,
        62,
        &|r| vec![r.signed_log_uniform(-2.0, 7.0), r.uniform(1.0, 3000.0)],
    ) {
        row(
            m,
            "entropy::heat_transfer_entropy",
            Class::B,
            &a,
            res(en::heat_transfer_entropy(a[0], a[1])),
        );
    }
}

fn carnot(m: &mut Module) {
    // carnot_efficiency / carnot_cop_refrigeration(T_hot, T_cold), with T_cold at, near and above T_hot
    let draw = |r: &mut Rng| vec![r.uniform(1.0, 3000.0), r.uniform(1.0, 3000.0)];
    let mut v = tuples(&[vec![500.0, 300.0], vec![300.0, 250.0]], 400, 63, &draw);
    for th in [300.0, 1.0, 1e-300, 1e300] {
        for tc in around(th, 2) {
            v.push(vec![th, tc]);
        }
    }
    v.extend([
        vec![200.0, 300.0],
        vec![f64::MAX, 1.0],
        vec![1.0, f64::MIN_POSITIVE],
    ]);
    for a in &v {
        row(
            m,
            "entropy::carnot_efficiency",
            Class::B,
            a,
            res(en::carnot_efficiency(a[0], a[1])),
        );
        row(
            m,
            "entropy::carnot_cop_refrigeration",
            Class::B,
            a,
            res(en::carnot_cop_refrigeration(a[0], a[1])),
        );
    }
}

fn potentials(m: &mut Module) {
    // helmholtz(U, T, S) and gibbs(H, T, S)
    let draw = |r: &mut Rng| {
        vec![
            r.signed_log_uniform(-2.0, 7.0),
            r.uniform(0.0, 3000.0),
            r.signed_log_uniform(-3.0, 4.0),
        ]
    };
    for a in tuples(&[vec![1000.0, 300.0, 5.0]], 300, 64, &draw) {
        row(
            m,
            "entropy::helmholtz",
            Class::B,
            &a,
            vec![f(en::helmholtz(a[0], a[1], a[2]))],
        );
    }
    for a in tuples(&[vec![2000.0, 300.0, 5.0]], 300, 65, &draw) {
        row(
            m,
            "entropy::gibbs",
            Class::B,
            &a,
            vec![f(en::gibbs(a[0], a[1], a[2]))],
        );
    }
    // is_spontaneous(ΔS)
    let mut sp = tuples(&[vec![0.1]], 100, 66, &|r| {
        vec![r.signed_log_uniform(-310.0, 300.0)]
    });
    sp.extend([vec![f64::from_bits(1)], vec![-f64::from_bits(1)]]);
    for a in sp {
        row(
            m,
            "entropy::is_spontaneous",
            Class::A,
            &a,
            vec![b(en::is_spontaneous(a[0]))],
        );
    }
}

fn mixing(m: &mut Module) {
    // entropy_of_mixing(n, xs): `f:n` then `v:len` and the fractions
    let nan = f64::NAN;
    let mut cases: Vec<(f64, Vec<f64>)> = vec![
        (1.0, vec![0.5, 0.5]),
        (1.0, vec![0.2, 0.3, 0.5]),
        (1.0, vec![1.0]),
        (1.0, vec![0.0, 1.0]),
        (1.0, vec![-0.0, 1.0]),
        (1.0, vec![0.3, 0.3]),
        (1.0, vec![]),
        (1.0, vec![nan, 1.0]),
        (1.0, vec![1.5, -0.5]),
        (1.0, vec![f64::from_bits(1), 1.0]),
        (1.0, vec![1e-300, 1.0]),
        (1.0, vec![f64::MIN_POSITIVE, 1.0]),
        (2.5, vec![0.78, 0.21, 0.01]),
    ];
    for n in crate::sweep::edges() {
        cases.push((n, vec![0.5, 0.5]));
    }
    for x in around(1.0 - 1e-10, 3).into_iter().chain(around(1.0, 2)) {
        cases.push((1.0, vec![x]));
    }
    for x in around(0.5 + 1e-10, 3) {
        cases.push((1.0, vec![0.5, x]));
    }
    let mut r = Rng::new(67);
    for _ in 0..300 {
        let len = r.int(1, 8) as usize;
        let mut xs: Vec<f64> = (0..len).map(|_| r.unit()).collect();
        if r.int(0, 3) != 0 {
            let s: f64 = xs.iter().sum();
            for x in &mut xs {
                *x /= s;
            }
        }
        cases.push((r.log_uniform(-2.0, 3.0), xs));
    }
    for (n, xs) in &cases {
        let mut ins = vec![f(*n), Tok::V(xs.len())];
        ins.extend(xs.iter().map(|&x| f(x)));
        m.row(
            "entropy::entropy_of_mixing",
            Class::C,
            &ins,
            &res(en::entropy_of_mixing(*n, xs)),
        );
    }
}

pub fn emit() -> Module {
    let mut m = Module::new("entropy");
    changes(&mut m);
    carnot(&mut m);
    potentials(&mut m);
    mixing(&mut m);
    m
}
