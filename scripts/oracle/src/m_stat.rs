//! `stat.tsv` (M5b): every public function of Rust 1.4.0's `stat`.
//!
//! Each function seeded over its physical ranges, then each argument of its base cases at the edges, and every switch
//! either side: Einstein's x = Θ_E/T at 500 ± 1 ulp and where (eˣ − 1)² falls below 1e-30 (bisected) ± 2 ulps; Debye's
//! Θ_D/T at 0.01 and 50 ± 1 ulp, the reference values of the 1.4.0 tests, and a sweep across the integrated range;
//! the vibrational partition function where hν ≪ kT and where the exponent underflows to 0 (`DivisionByZero`). The
//! rotational symmetry number, Rust's `u32`, is an `i:` token (0, 1, 2, 3, 12, u32::MAX). Classes: B for + − × ÷ √,
//! C where exp, expm1, ln or powf enters, D for `debye_cv` (hisab's Simpson sums in another order; probe 1).

use ushma::stat as st;

use crate::sweep::{around, res, row, tuples};
use crate::tsv::{Class, Module, f, i};

const N2_MASS: f64 = 4.651_754e-26;

/// The smallest positive double in (lo, hi] where `up` holds, for `up` false at lo and true at hi, switching once.
fn first_true(lo: f64, hi: f64, up: &dyn Fn(f64) -> bool) -> f64 {
    let (mut a, mut b) = (lo.to_bits(), hi.to_bits());
    while b - a > 1 {
        let m = a + (b - a) / 2;
        if up(f64::from_bits(m)) {
            b = m;
        } else {
            a = m;
        }
    }
    f64::from_bits(b)
}

fn boltzmann(m: &mut Module) {
    for a in tuples(
        &[vec![1e-21, 300.0], vec![0.0, 300.0], vec![-0.0, 300.0]],
        300,
        110,
        &|r| vec![r.uniform(-1e-19, 1e-19), r.log_uniform(-1.0, 4.0)],
    ) {
        row(
            m,
            "stat::boltzmann_probability",
            Class::C,
            &a,
            res(st::boltzmann_probability(a[0], a[1])),
        );
    }
    for a in tuples(
        &[vec![0.0, 1e-20, 300.0], vec![1.0, 1.0, 300.0]],
        300,
        111,
        &|r| {
            vec![
                r.uniform(-1e-19, 1e-19),
                r.uniform(-1e-19, 1e-19),
                r.log_uniform(-1.0, 4.0),
            ]
        },
    ) {
        row(
            m,
            "stat::boltzmann_ratio",
            Class::C,
            &a,
            res(st::boltzmann_ratio(a[0], a[1], a[2])),
        );
    }
    let mut w = tuples(&[vec![100.0]], 300, 112, &|r| {
        vec![r.log_uniform(0.0, 300.0)]
    });
    w.extend(around(1.0, 2).into_iter().map(|x| vec![x]));
    for a in w {
        row(
            m,
            "stat::entropy_from_microstates",
            Class::C,
            &a,
            res(st::entropy_from_microstates(a[0])),
        );
    }
    for a in tuples(&[vec![3.0, 300.0]], 200, 113, &|r| {
        vec![r.uniform(0.0, 12.0), r.log_uniform(-1.0, 4.0)]
    }) {
        row(
            m,
            "stat::equipartition_energy",
            Class::B,
            &a,
            vec![f(st::equipartition_energy(a[0], a[1]))],
        );
    }
    for a in tuples(&[vec![3.0]], 200, 114, &|r| vec![r.uniform(0.0, 12.0)]) {
        row(
            m,
            "stat::equipartition_cv",
            Class::B,
            &a,
            vec![f(st::equipartition_cv(a[0]))],
        );
    }
}

fn partitions(m: &mut Module) {
    for a in tuples(&[vec![N2_MASS, 0.02241, 300.0]], 300, 115, &|r| {
        vec![
            r.log_uniform(-27.0, -24.0),
            r.log_uniform(-6.0, 1.0),
            r.log_uniform(0.0, 4.0),
        ]
    }) {
        row(
            m,
            "stat::partition_translational",
            Class::C,
            &a,
            res(st::partition_translational(a[0], a[1], a[2])),
        );
    }
    // rotational: I and T seeded and at the edges, σ over its own values
    let sigmas: [u32; 6] = [0, 1, 2, 3, 12, u32::MAX];
    let mut rot: Vec<(f64, f64, u32)> = Vec::new();
    for a in tuples(&[vec![1.4e-46, 300.0]], 200, 116, &|r| {
        vec![r.log_uniform(-48.0, -44.0), r.log_uniform(0.0, 4.0)]
    }) {
        rot.push((a[0], a[1], 2));
    }
    for &s in &sigmas {
        rot.push((1.4e-46, 300.0, s));
        rot.push((2.0e-47, 50.0, s));
    }
    for (inertia, t, sigma) in rot {
        m.row(
            "stat::partition_rotational",
            Class::B,
            &[f(inertia), f(t), i(sigma)],
            &res(st::partition_rotational(inertia, t, sigma)),
        );
    }
    // vibrational: seeded, the edges, hν ≪ kT, and an exponent that underflows to 0
    let mut vib = tuples(&[vec![1e13, 300.0], vec![1.0, 1e6]], 300, 117, &|r| {
        vec![r.log_uniform(9.0, 15.0), r.log_uniform(0.0, 4.0)]
    });
    vib.extend([
        vec![1e-300, 1e300],
        vec![1e-200, 1e200],
        vec![1e15, 1.0],
        vec![5e13, 3.0],
    ]);
    for a in vib {
        row(
            m,
            "stat::partition_vibrational",
            Class::C,
            &a,
            res(st::partition_vibrational(a[0], a[1])),
        );
    }
}

fn speeds(m: &mut Module) {
    for a in tuples(
        &[vec![500.0, N2_MASS, 300.0], vec![0.0, N2_MASS, 300.0]],
        300,
        118,
        &|r| {
            vec![
                r.uniform(0.0, 3000.0),
                r.log_uniform(-27.0, -24.0),
                r.log_uniform(0.0, 4.0),
            ]
        },
    ) {
        row(
            m,
            "stat::maxwell_boltzmann_speed_pdf",
            Class::C,
            &a,
            res(st::maxwell_boltzmann_speed_pdf(a[0], a[1], a[2])),
        );
    }
    let mt = tuples(&[vec![N2_MASS, 300.0]], 200, 119, &|r| {
        vec![r.log_uniform(-27.0, -24.0), r.log_uniform(0.0, 4.0)]
    });
    for a in &mt {
        row(
            m,
            "stat::most_probable_speed",
            Class::B,
            a,
            res(st::most_probable_speed(a[0], a[1])),
        );
        row(
            m,
            "stat::mean_speed",
            Class::B,
            a,
            res(st::mean_speed(a[0], a[1])),
        );
        row(
            m,
            "stat::rms_speed",
            Class::B,
            a,
            res(st::rms_speed(a[0], a[1])),
        );
    }
}

fn solids(m: &mut Module) {
    // Einstein: seeded, the edges, x at 500 ± 1 ulp, and where (eˣ − 1)² < 1e-30 (Rust's own test) ± 2 ulps
    let mut ein = tuples(&[vec![300.0, 300.0]], 300, 120, &|r| {
        vec![r.uniform(10.0, 2000.0), r.log_uniform(-1.0, 4.0)]
    });
    for x in around(500.0, 1) {
        ein.push(vec![x, 1.0]);
    }
    let limit = |th: f64| {
        let ex = th.exp();
        ((ex - 1.0) * (ex - 1.0)).abs() >= 1e-30
    };
    for x in around(first_true(1e-300, 1.0, &limit), 2) {
        ein.push(vec![x, 1.0]);
    }
    ein.extend([vec![200.0, 5000.0], vec![1000.0, 10.0], vec![300.0, 50.0]]);
    for a in ein {
        row(
            m,
            "stat::einstein_cv",
            Class::C,
            &a,
            res(st::einstein_cv(a[0], a[1])),
        );
    }
    // Debye: seeded, the edges, Θ/T across the integrated range and at both switches ± 1 ulp, 1.4.0's reference cases
    let mut deb = tuples(&[vec![300.0, 300.0]], 300, 121, &|r| {
        vec![r.uniform(50.0, 2500.0), r.log_uniform(-1.0, 4.0)]
    });
    for x in around(0.01, 1).into_iter().chain(around(50.0, 1)) {
        deb.push(vec![x, 1.0]);
    }
    for &x in &[
        0.005, 0.02, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 49.999_999, 50.000_001,
        100.0, 709.0, 1000.0,
    ] {
        deb.push(vec![x * 100.0, 100.0]);
    }
    deb.extend([
        vec![2230.0, 3.0],
        vec![300.0, 0.42],
        vec![428.0, 1.0],
        vec![300.0, 3000.0],
        vec![300.0, 5.0],
    ]);
    for a in deb {
        row(
            m,
            "stat::debye_cv",
            Class::D,
            &a,
            res(st::debye_cv(a[0], a[1])),
        );
    }
}

pub fn emit() -> Module {
    let mut m = Module::new("stat");
    boltzmann(&mut m);
    partitions(&mut m);
    speeds(&mut m);
    solids(&mut m);
    m
}
