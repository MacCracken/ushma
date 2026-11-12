//! `math.tsv` — what Rust **std** returns for every f64 method ushma 1.4.0's library code
//! calls: the reference for the M0 step 8 baselines (port plan §4 step 5, §7.1).
//!
//! The calls, counted in `rust-old/src` outside the test modules, whole crate (2026-11-11): abs 30,
//! powf 29, sqrt 23, ln 23, is_finite 21, exp 14, max 11, atan 11, clamp 9, powi 7, min 6,
//! is_nan 5, tanh 2, exp_m1 2, ceil 1. (Some `max`/`min`/`clamp` calls are on `usize` and are
//! not recorded here.)
//!
//! Where each result comes from (rustc 1.99.1, x86_64-unknown-linux-gnu; the oracle binary's
//! imports, `nm -D`, are listed in scripts/oracle/README.md):
//! - **glibc's libm** (dynamic): exp, expm1, log, pow, atan, tanh. These rows are glibc's bits
//!   (`glibc.rs`): the header records the glibc, and `--check` refuses another.
//! - **inline code or compiler-builtins**: sqrt (`sqrtsd`), abs, ceil, min, max, clamp,
//!   is_finite, is_nan, and `__powidf2` for powi. min/max lower to `minsd`/`maxsd`, so on
//!   (+0, −0) they return the SECOND operand; Rust leaves that case unspecified, and the rows
//!   record what 1.99.1 does.
//!
//! Every input passes through `black_box`, so each row is a run-time library call, never a
//! compile-time fold. ushma's `powf` calls with a literal exponent (`re.powf(0.8)`,
//! `ra.powf(1.0 / 6.0)`, …) are recorded as `f64::powf` rows with that exponent's exact Rust
//! value (a literal, or a constant-folded `a / b`, correctly rounded), over the ranges the call
//! sites feed it.
//!
//! Classes:
//! - **A**: exact on both sides — sqrt, abs, ceil, min, max, clamp, is_finite, is_nan.
//! - **B**: powi (repeated multiplication in `__powidf2`'s order).
//! - **C**: exp, exp_m1, ln, powf, atan, tanh. ADR 0004: the port's primitives are correctly
//!   rounded, so a C row's expected value is the corrected oracle (`<module>.cr.tsv`, M0 part 2c)
//!   wherever glibc misrounds.

use std::hint::black_box;

use crate::grids::{EXP_SUBNORMAL_X, Rng, decades, pow2, specials};
use crate::tsv::{Class, Module, b, f, i};

fn unary(m: &mut Module, name: &str, class: Class, fun: fn(f64) -> f64, xs: &[f64]) {
    for &x in xs {
        let y = fun(black_box(x));
        m.row(name, class, &[f(x)], &[f(y)]);
    }
}

fn binary(m: &mut Module, name: &str, class: Class, fun: fn(f64, f64) -> f64, xy: &[(f64, f64)]) {
    for &(x, y) in xy {
        let r = fun(black_box(x), black_box(y));
        m.row(name, class, &[f(x), f(y)], &[f(r)]);
    }
}

/// Every ordered pair of `xs`.
fn grid(xs: &[f64]) -> Vec<(f64, f64)> {
    xs.iter()
        .flat_map(|&x| xs.iter().map(move |&y| (x, y)))
        .collect()
}

/// Specials, ±2^k at every 16th binade, and a wide signed log-uniform sweep.
fn wide(seed: u64, n: usize) -> Vec<f64> {
    let mut v = specials();
    for k in (-1074..=1023).step_by(16) {
        let x = pow2(k);
        v.push(x);
        v.push(-x);
    }
    let mut rng = Rng::new(seed);
    v.extend((0..n).map(|_| rng.signed_log_uniform(-300.0, 300.0)));
    v
}

fn sweep(seed: u64, n: usize, lo: f64, hi: f64) -> Vec<f64> {
    let mut rng = Rng::new(seed);
    (0..n).map(|_| rng.uniform(lo, hi)).collect()
}

fn log_sweep(seed: u64, n: usize, lo_exp10: f64, hi_exp10: f64) -> Vec<f64> {
    let mut rng = Rng::new(seed);
    (0..n)
        .map(|_| rng.log_uniform(lo_exp10, hi_exp10))
        .collect()
}

fn pairs(seed: u64, n: usize, lo_exp10: f64, hi_exp10: f64) -> Vec<(f64, f64)> {
    let mut rng = Rng::new(seed);
    (0..n)
        .map(|_| {
            let x = rng.signed_log_uniform(lo_exp10, hi_exp10);
            (x, rng.signed_log_uniform(lo_exp10, hi_exp10))
        })
        .collect()
}

/// `xs` paired with the one exponent `y`.
fn with_exp(xs: Vec<f64>, y: f64) -> Vec<(f64, f64)> {
    xs.into_iter().map(|x| (x, y)).collect()
}

fn clamp_row(m: &mut Module, x: f64, lo: f64, hi: f64) {
    let r = black_box(x).clamp(black_box(lo), black_box(hi));
    m.row("f64::clamp", Class::A, &[f(x), f(lo), f(hi)], &[f(r)]);
}

#[allow(clippy::too_many_lines)] // one grid per function, in the module comment's order
pub fn emit() -> Module {
    let mut m = Module::new("math");

    // ── exact (class A) ──
    // sqrt: view factors √(1 + x²), RMS / mean speeds √(kT/m), fin parameters √(hP/kA), the IF97
    // Region 4 roots, Kay's-rule and Peng–Robinson terms; squares exactly.
    let mut xs = wide(0x5157, 300);
    xs.extend(log_sweep(0x5158, 300, -40.0, 40.0));
    xs.extend((1..=100).map(|k| f64::from(k * k)));
    xs.extend(sweep(0x5159, 100, 0.0, 1e4));
    unary(&mut m, "f64::sqrt", Class::A, f64::sqrt, &xs);
    unary(&mut m, "f64::abs", Class::A, f64::abs, &wide(0xab5, 100));

    // ceil: the adaptive step count ratio·(1 − roundoff) for ratios 1 … 1e3 (and the integers
    // themselves), ties and near-integers both ways, then a wide sweep.
    let mut ce: Vec<f64> = (-20..=20).map(|k| f64::from(k) + 0.5).collect();
    ce.extend((0..=1000).step_by(7).map(f64::from));
    ce.extend((1..=200).map(|k| f64::from(k) + f64::EPSILON * f64::from(k)));
    ce.extend((1..=200).map(|k| f64::from(k) * (1.0 - f64::EPSILON)));
    ce.extend(sweep(0xce11, 300, 0.0, 1e3));
    ce.extend(wide(0xce12, 100));
    unary(&mut m, "f64::ceil", Class::A, f64::ceil, &ce);

    // min / max: every special pair (NaN on either side, ±0 both ways, ±∞), equal pairs, a value
    // against 0.0 as ushma writes `x.max(0.0)`, a value against 1.0 (`.max(1.0)`), and random pairs.
    let mut mm = grid(&specials());
    let mut rng = Rng::new(0x3a3a);
    for _ in 0..50 {
        let x = rng.signed_log_uniform(-30.0, 30.0);
        mm.push((x, x));
        mm.push((x, 0.0));
        mm.push((0.0, x));
        mm.push((x, 1.0));
        mm.push((x, f64::NAN));
        mm.push((f64::NAN, x));
    }
    mm.extend(pairs(0x3a3b, 200, -300.0, 300.0));
    binary(&mut m, "f64::max", Class::A, f64::max, &mm);
    binary(&mut m, "f64::min", Class::A, f64::min, &mm);

    // clamp: ushma's bounds — (0, 1) (quality, efficiencies), (0, 100) (relative humidity), the
    // IF97 Region 4 range (273.15, 647.096) — then other valid bounds (lo ≤ hi, neither NaN;
    // Rust panics otherwise).
    let mut cl: Vec<f64> = specials();
    cl.extend(sweep(0xc1a0, 200, -0.5, 1.5));
    cl.extend([1.0 + f64::EPSILON, -f64::EPSILON, 1.0 - f64::EPSILON / 2.0]);
    for &x in &cl {
        clamp_row(&mut m, x, 0.0, 1.0);
    }
    for &x in specials().iter().chain(&sweep(0xc1a3, 60, -20.0, 120.0)) {
        clamp_row(&mut m, x, 0.0, 100.0);
    }
    for &x in specials().iter().chain(&sweep(0xc1a4, 60, 200.0, 700.0)) {
        clamp_row(&mut m, x, 273.15, 647.096);
    }
    let mut rng = Rng::new(0xc1a2);
    let mut bounds = vec![(0.0, 0.0), (-0.0, 0.0), (-0.0, -0.0), (0.0, f64::INFINITY)];
    for _ in 0..40 {
        let a = rng.signed_log_uniform(-10.0, 10.0);
        let c = rng.signed_log_uniform(-10.0, 10.0);
        bounds.push(if a <= c { (a, c) } else { (c, a) });
    }
    for &(lo, hi) in &bounds {
        for x in [
            lo,
            hi,
            0.0,
            -0.0,
            f64::NAN,
            rng.signed_log_uniform(-10.0, 10.0),
        ] {
            clamp_row(&mut m, x, lo, hi);
        }
    }

    let mut fin = specials();
    fin.extend([
        f64::from_bits(0xfff8_0000_0000_0000), // -NaN
        f64::from_bits(0x7ff0_0000_0000_0001), // signalling-pattern NaN
        -f64::MAX,
        f64::from_bits(0x7fef_ffff_ffff_ffff),
    ]);
    fin.extend(log_sweep(0xf1f1, 20, -300.0, 300.0));
    for &x in &fin {
        m.row(
            "f64::is_finite",
            Class::A,
            &[f(x)],
            &[b(black_box(x).is_finite())],
        );
        m.row(
            "f64::is_nan",
            Class::A,
            &[f(x)],
            &[b(black_box(x).is_nan())],
        );
    }

    // ── arithmetic (class B) ──
    // powi: `t.powi(4)` (Stefan–Boltzmann, radiation; T from 1e-3 K to 1e5 K and beyond) and
    // `.powi(2)`; n is black-boxed too, so every row is a `__powidf2` call.
    let mut pw = specials();
    pw.extend(sweep(0x9101, 200, 0.0, 6000.0));
    pw.extend(log_sweep(0x9102, 100, -3.0, 77.0));
    pw.extend(sweep(0x9103, 60, -10.0, 10.0));
    pw.extend(log_sweep(0x9104, 60, -150.0, 150.0));
    for n in [4, 2, -2, -1, 0, 1, 3] {
        let set: &[f64] = if n == 4 || n == 2 { &pw } else { &pw[..80] };
        for &x in set {
            let r = black_box(x).powi(black_box(n));
            m.row("f64::powi", Class::B, &[f(x), i(n)], &[f(r)]);
        }
    }

    // ── transcendental (class C) ──
    // exp: Boltzmann factors −E/kT and ε-NTU's exp(−NTU·(1 + C_r)) in [−700, 0] above all, the
    // equilibrium constant exp(−ΔG/RT) across the whole finite range and past it, Clausius–
    // Clapeyron and the vapor curves (ln P vs 1/T), the subnormal-result rows, integers, near 0.
    let mut ex = specials();
    ex.extend(sweep(0xe001, 400, -700.0, 0.0));
    ex.extend(sweep(0xe002, 200, -745.2, 709.8));
    ex.extend(sweep(0xe003, 200, -50.0, 50.0));
    ex.extend(sweep(0xe004, 100, -1e-6, 1e-6));
    ex.extend((-745..=709).step_by(7).map(f64::from));
    ex.extend(EXP_SUBNORMAL_X);
    ex.extend([
        709.782_712_893_384,
        709.782_712_893_384_1,
        -708.396_418_532_264_1,
        -745.133_219_101_941_1,
        -745.133_219_101_941_2,
        837.0, // propane ln K: past the f64 range (ln_equilibrium_constant)
        1e-300,
        -1e-300,
    ]);
    unary(&mut m, "f64::exp", Class::C, f64::exp, &ex);

    // exp_m1: both calls take −x — Debye's integrand x⁴e⁻ˣ/expm1(−x)² for x in (0, 50] and
    // partition_vibrational's 1 − e^(−hν/kT) for hν/kT from 1e-10 to past 745.
    let mut em = specials();
    em.extend(sweep(0xe101, 100, -1e-6, 1e-6));
    em.extend(sweep(0xe102, 200, -1.0, 1.0));
    em.extend(sweep(0xe103, 300, -50.0, 0.0));
    em.extend(sweep(0xe104, 100, -745.2, 709.8));
    em.extend(log_sweep(0xe105, 300, -10.0, 3.0).into_iter().map(|x| -x));
    unary(&mut m, "f64::exp_m1", Class::C, f64::exp_m1, &em);

    // ln: pressure, volume and temperature ratios (entropy, isothermal work, adiabatic exponents)
    // close to 1 and far from it, LMTD's ln(ΔT₁/ΔT₂) as ΔT₁ → ΔT₂, mole fractions x ln x in
    // (0, 1], vapor pressures, equilibrium constants, a wide sweep, decades.
    let mut lg = wide(0x1060, 300);
    lg.extend(sweep(0x1061, 200, 0.5, 2.0));
    lg.extend(sweep(0x1062, 100, 1.0 - 1e-9, 1.0 + 1e-9));
    lg.extend(sweep(0x1063, 200, 0.0, 1.0));
    lg.extend(log_sweep(0x1064, 300, -300.0, 300.0));
    lg.extend(decades(-40, 30, &["1", "2", "5"]));
    unary(&mut m, "f64::ln", Class::C, f64::ln, &lg);

    // powf, both operands run-time: random positive pairs, the special grid, negative bases with
    // integer exponents.
    let mut pf: Vec<(f64, f64)> = Vec::new();
    let mut rng = Rng::new(0x9002);
    for _ in 0..300 {
        let x = rng.log_uniform(-30.0, 30.0);
        pf.push((x, rng.uniform(-10.0, 10.0)));
    }
    pf.extend(grid(&specials()));
    for _ in 0..50 {
        let x = -rng.uniform(0.1, 10.0);
        #[allow(clippy::cast_precision_loss)]
        let y = rng.int(-5, 5) as f64;
        pf.push((x, y));
    }
    // The γ − 1 family: (V₁/V₂)^(γ−1) and (P₂/P₁)^((γ−1)/γ) for compression and pressure ratios
    // 1 … 30, γ 1.01 … 1.67 (Otto, Diesel, Brayton, adiabatic relations); T_a·(T_b/T_a)^f along a
    // diagram path, f in [0, 1].
    let mut rng = Rng::new(0x9003);
    for _ in 0..300 {
        let r = rng.uniform(1.0, 30.0);
        let gamma = rng.uniform(1.01, 1.67);
        pf.push((r, gamma - 1.0));
        pf.push((r, (gamma - 1.0) / gamma));
        pf.push((1.0 / r, gamma - 1.0));
    }
    for _ in 0..200 {
        let ratio = rng.uniform(0.2, 10.0);
        pf.push((ratio, rng.unit()));
    }
    binary(&mut m, "f64::powf", Class::C, f64::powf, &pf);
    // The literal exponents, each over the range its call site feeds it. The exponent is Rust's
    // value of the literal as written (`9.0 / 16.0` is a correctly rounded constant fold).
    let literal: [(f64, Vec<f64>); 13] = [
        (1.5, log_sweep(0x9010, 200, -6.0, 32.0)), // Maxwell–Boltzmann, translational partition, Stull
        (2.5, sweep(0x9011, 150, 5.0, 1000.0)),    // Redlich–Kwong Tc^2.5
        (0.8, log_sweep(0x9012, 150, 2.0, 7.5)),   // Dittus–Boelter Re^0.8
        (0.4, log_sweep(0x9013, 150, -1.0, 3.5)),  // Dittus–Boelter Pr^0.4
        (1.6, sweep(0x9014, 150, 0.2, 5.0)),       // Pitzer B0: Tr^1.6
        (4.2, sweep(0x9015, 150, 0.2, 5.0)),       // Pitzer B1: Tr^4.2
        (0.75, log_sweep(0x9016, 150, -3.0, 4.0)), // Kleiber: mass^0.75
        (0.425, sweep(0x9017, 150, 1.0, 300.0)),   // DuBois: mass^0.425
        (0.725, sweep(0x9018, 150, 40.0, 250.0)),  // DuBois: height^0.725
        (9.0 / 16.0, log_sweep(0x9019, 150, -3.5, 0.5)), // Churchill–Chu: (0.492/Pr)^(9/16)
        (-8.0 / 27.0, sweep(0x901a, 150, 1.0, 3.0)), // Churchill–Chu: [1 + …]^(−8/27)
        (1.0 / 6.0, log_sweep(0x901b, 150, -2.0, 13.0)), // Churchill–Chu: Ra^(1/6)
        (1.0 / 3.0, log_sweep(0x901c, 150, -1.0, 4.0)), // Blasius: Pr^(1/3)
    ];
    for (y, xs) in literal {
        let mut xs = xs;
        xs.extend([0.0, 1.0, 2.0, 10.0, f64::INFINITY]);
        binary(&mut m, "f64::powf", Class::C, f64::powf, &with_exp(xs, y));
    }

    // atan: the view-factor terms atan(x/√(1+y²)), x·atan(x), w·atan(1/w) for aspect ratios 1e-6
    // … 1e6, Stull's wet-bulb terms over RH 5 … 99 %, a general sweep, near 0, wide.
    let mut at = wide(0xa7a0, 200);
    at.extend(log_sweep(0xa7a1, 400, -6.0, 6.0));
    at.extend(sweep(0xa7a2, 300, -1e3, 1e3));
    at.extend(sweep(0xa7a3, 200, -2.0, 2.0));
    at.extend(sweep(0xa7a4, 100, -1e-6, 1e-6));
    at.extend(sweep(0xa7a5, 150, 0.0, 120.0)); // t + RH, RH − 1.676…, 0.023101·RH, 0.151977·√(RH + 8.31)
    unary(&mut m, "f64::atan", Class::C, f64::atan, &at);

    // tanh: fin heat transfer tanh(mL) and fin efficiency tanh(mL)/mL, mL in [0, 50] above all,
    // tiny mL (the efficiency's limit at 1), large, negative, specials.
    let mut th = specials();
    th.extend(sweep(0x7a00, 500, 0.0, 50.0));
    th.extend(log_sweep(0x7a01, 200, -12.0, 0.0));
    th.extend(sweep(0x7a02, 100, -50.0, 0.0));
    th.extend(log_sweep(0x7a03, 50, 1.5, 300.0));
    unary(&mut m, "f64::tanh", Class::C, f64::tanh, &th);

    m
}
