//! `fmtf.tsv` — the f64 text ushma writes (M1a): Rust's `Display` (`{}`) and serde_json 1.0.154's form
//! (`serde_json::to_string`), for `src/fmtf.cyr`, the shortest-digit printer vendored from impetus (ADR 0005).
//!
//! The doubles: the specials (±0, ±∞, NaN, the extremes, the smallest subnormal), every power of two (2^-1074 …
//! 2^1023, which include the exact ties the two Rust printers break differently: core up, zmij to even), the powers
//! of ten 1e-323 … 1e308 and their neighbours, the form thresholds (serde_json switches to an exponent below 1e-5 and
//! from 1e16), 0 … 1000 by tenths and hundredths, and a seeded sweep of 4,000 random bit patterns and 2,000 short
//! decimals. Every row is class A: the bytes.

use crate::grids::Rng;
use crate::tsv::{Class, Module, f, s};

fn doubles() -> Vec<f64> {
    let mut xs = vec![
        0.0,
        -0.0,
        f64::MAX,
        f64::MIN,
        f64::MIN_POSITIVE,
        f64::EPSILON,
        f64::NAN,
        f64::INFINITY,
    ];
    xs.extend([
        f64::NEG_INFINITY,
        f64::from_bits(1),
        f64::from_bits(0x000F_FFFF_FFFF_FFFF),
    ]);
    for k in -1074..=1023 {
        xs.push(2f64.powi(k));
    }
    for e in -323..=308 {
        let p: f64 = format!("1e{e}").parse().expect("power of ten");
        xs.extend([
            p,
            f64::from_bits(p.to_bits() + 1),
            f64::from_bits(p.to_bits() - 1),
        ]);
    }
    for t in [
        "1e15",
        "1e16",
        "1e17",
        "9007199254740992",
        "9007199254740993",
        "1e-5",
        "1e-6",
        "1e-7",
        "0.0001",
    ] {
        let x: f64 = t.parse().expect("threshold");
        xs.extend([
            x,
            f64::from_bits(x.to_bits() + 1),
            f64::from_bits(x.to_bits() - 1),
            -x,
        ]);
    }
    for i in 0..=1000 {
        xs.push(f64::from(i));
        xs.push(f64::from(i) / 10.0);
        xs.push(f64::from(i) / 100.0);
    }
    let mut r = Rng::new(0x6d74_6674_7573_6875);
    for _ in 0..4000 {
        xs.push(f64::from_bits(r.bits()));
    }
    for _ in 0..2000 {
        let digits = (r.bits() % 17) as u32 + 1;
        let m = r.bits() % 10u64.pow(digits);
        let e = (r.bits() % 60) as i32 - 30;
        let x: f64 = format!("{m}e{e}").parse().expect("decimal");
        xs.push(if r.bits() & 1 == 1 { -x } else { x });
    }
    xs
}

pub fn emit() -> Module {
    let mut m = Module::new("fmtf");
    for x in doubles() {
        m.row("fmtf::display", Class::A, &[f(x)], &[s(format!("{x}"))]);
        let j = serde_json::to_string(&x).expect("serde_json f64");
        m.row("fmtf::json", Class::A, &[f(x)], &[s(j)]);
    }
    m
}
