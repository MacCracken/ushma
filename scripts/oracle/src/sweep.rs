//! The argument sweeps every module's function rows share (M1c onward): a function's seeded draws over its physical
//! ranges, then each argument of each base case at every edge value, and the neighbours of a branch switch.

use crate::grids::Rng;
use crate::tsv::{Class, Module, Tok, err, f};

/// A result as its tokens.
pub fn res(r: ushma::error::Result<f64>) -> Vec<Tok> {
    match r {
        Ok(x) => vec![f(x)],
        Err(e) => vec![err(&e)],
    }
}

/// The edge values every argument is tried at.
pub fn edges() -> Vec<f64> {
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
pub fn tuples(
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

/// One row of f64 inputs.
pub fn row(m: &mut Module, name: &str, class: Class, args: &[f64], out: Vec<Tok>) {
    let ins: Vec<Tok> = args.iter().map(|&x| f(x)).collect();
    m.row(name, class, &ins, &out);
}

/// x's neighbours `k` ulps either side (positive finite x), and x itself.
pub fn around(x: f64, k: u64) -> Vec<f64> {
    let b = x.to_bits();
    ((b - k)..=(b + k)).map(f64::from_bits).collect()
}
