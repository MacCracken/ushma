//! `numerical.tsv` (M6a onward): Rust 1.4.0's `numerical`.
//!
//! A boundary condition is written as its variant (`e:BoundaryCondition.Fixed`, `.Insulated`, `.Convective`) and the
//! values it carries; a 1D grid as `v:n` and its nodes, then Δx, α and its two conditions. A step's row is a grid state
//! and the step in, the grid after it out (or the refusal), so every step is checked from the oracle's own state.
//!
//! M6a: `ThermalGrid1D::new` over node counts, lengths, diffusivities, initial temperatures and every condition at its
//! edges; `fourier_number` and `max_stable_dt`; `step_explicit`, `step_crank_nicolson` and `step_adaptive` on seven
//! setups (fixed, insulated and convective ends, Biot numbers from 0 to 10, 3 to 50 nodes), each along a trajectory of
//! consecutive states, at fractions and multiples of the stable step, at the explicit stability limit (bisected) ± 2 ulps,
//! at the adaptive exact multiples ± 1 ulp and its sub-step cap (bisected, on the 10-node setup) ± 2 ulps, and on malformed grids (fewer than
//! 3 nodes, Δx or α not positive, a condition out of range, a NaN node). JSON: `BoundaryCondition`, `Side` and
//! `ThermalGrid1D`. Class B throughout (+ − × ÷ and ceil), A for JSON.

use ushma::numerical::{BoundaryCondition as Bc, Side, ThermalGrid1D};

use crate::json_cases::variants;
use crate::sweep::{around, edges};
use crate::tsv::{Class, Module, Tok, err, f, i, s};

fn bc_toks(b: &Bc) -> Vec<Tok> {
    match *b {
        Bc::Fixed(t) => vec![Tok::E("BoundaryCondition", "Fixed".into()), f(t)],
        Bc::Insulated => vec![Tok::E("BoundaryCondition", "Insulated".into())],
        Bc::Convective { h, t_fluid, k } => vec![
            Tok::E("BoundaryCondition", "Convective".into()),
            f(h),
            f(t_fluid),
            f(k),
        ],
        _ => unreachable!("1.4.0 has three conditions"),
    }
}

fn g1_toks(g: &ThermalGrid1D) -> Vec<Tok> {
    let mut v = vec![Tok::V(g.nodes.len())];
    v.extend(g.nodes.iter().map(|&x| f(x)));
    v.extend([f(g.dx), f(g.alpha)]);
    v.extend(bc_toks(&g.bc_left));
    v.extend(bc_toks(&g.bc_right));
    v
}

fn conv(h: f64, t_fluid: f64, k: f64) -> Bc {
    Bc::Convective { h, t_fluid, k }
}

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

/// The seven step setups: (nodes, length, α, T₀, left, right), and a profile to start from (None: uniform).
fn setups() -> Vec<ThermalGrid1D> {
    let mut v = Vec::new();
    let mk = |n: usize, l: f64, a: f64, t: f64, bl: Bc, br: Bc| {
        ThermalGrid1D::new(n, l, a, t, bl, br).expect("setup")
    };
    v.push(mk(10, 1.0, 1e-4, 300.0, Bc::Fixed(400.0), Bc::Fixed(300.0)));
    v.push(mk(11, 0.5, 1e-5, 300.0, Bc::Fixed(400.0), Bc::Insulated));
    v.push(mk(
        21,
        0.2,
        4e-5,
        300.0,
        Bc::Fixed(500.0),
        conv(25.0, 300.0, 50.0),
    ));
    v.push(mk(
        8,
        0.01,
        1e-6,
        350.0,
        conv(1000.0, 500.0, 1.0),
        conv(10.0, 250.0, 40.0),
    ));
    let mut e = mk(12, 1.0, 2e-5, 300.0, Bc::Insulated, Bc::Insulated);
    for (k, x) in e.nodes.iter_mut().enumerate() {
        *x = 300.0 + 50.0 * ((k as f64) * 0.7).sin();
    }
    v.push(e);
    v.push(mk(
        3,
        0.1,
        1e-4,
        280.0,
        conv(0.0, 300.0, 10.0),
        Bc::Fixed(320.0),
    ));
    v.push(mk(
        50,
        2.0,
        1e-4,
        300.0,
        Bc::Fixed(373.15),
        Bc::Fixed(273.15),
    ));
    v
}

/// Grids the constructor would refuse, built from Rust's public fields (as a deserialized or edited grid can be).
fn malformed() -> Vec<ThermalGrid1D> {
    let base = setups().remove(0);
    let edit = |f: &dyn Fn(&mut ThermalGrid1D)| {
        let mut g = base.clone();
        f(&mut g);
        g
    };
    vec![
        edit(&|g| g.nodes.truncate(2)),
        edit(&|g| g.nodes.clear()),
        edit(&|g| g.dx = 0.0),
        edit(&|g| g.dx = f64::NAN),
        edit(&|g| g.dx = -1.0),
        edit(&|g| g.alpha = 0.0),
        edit(&|g| g.alpha = f64::INFINITY),
        edit(&|g| g.bc_left = Bc::Fixed(f64::NAN)),
        edit(&|g| g.bc_right = conv(25.0, 300.0, 0.0)),
        edit(&|g| g.bc_right = conv(-1.0, 300.0, 50.0)),
        edit(&|g| g.bc_left = conv(25.0, f64::INFINITY, 50.0)),
        edit(&|g| g.nodes[4] = f64::NAN),
        edit(&|g| g.nodes[0] = f64::INFINITY),
    ]
}

fn grid_rows(m: &mut Module) {
    // new(n, L, α, T₀, left, right)
    let bcs = [
        Bc::Fixed(400.0),
        Bc::Insulated,
        conv(25.0, 300.0, 50.0),
        conv(0.0, 300.0, 50.0),
        Bc::Fixed(f64::NAN),
        Bc::Fixed(f64::INFINITY),
        Bc::Fixed(-0.0),
        conv(-1.0, 300.0, 50.0),
        conv(f64::NAN, 300.0, 50.0),
        conv(25.0, 300.0, 0.0),
        conv(25.0, 300.0, -0.0),
        conv(25.0, f64::NAN, 50.0),
        conv(f64::INFINITY, 300.0, 50.0),
    ];
    let mut news: Vec<(usize, f64, f64, f64, Bc, Bc)> = Vec::new();
    for n in [0, 1, 2, 3, 4, 10] {
        news.push((n, 1.0, 1e-4, 300.0, Bc::Fixed(400.0), Bc::Fixed(300.0)));
    }
    for e in edges() {
        news.push((10, e, 1e-4, 300.0, Bc::Insulated, Bc::Insulated));
        news.push((10, 1.0, e, 300.0, Bc::Insulated, Bc::Insulated));
        news.push((10, 1.0, 1e-4, e, Bc::Insulated, Bc::Insulated));
    }
    for b in &bcs {
        news.push((5, 1.0, 1e-4, 300.0, *b, Bc::Insulated));
        news.push((5, 1.0, 1e-4, 300.0, Bc::Insulated, *b));
    }
    news.push((
        7,
        0.3,
        1e-4,
        300.0,
        Bc::Fixed(f64::NAN),
        conv(-1.0, 1.0, 1.0),
    ));
    for (n, l, a, t, bl, br) in news {
        let mut ins = vec![i(i64::try_from(n).expect("small")), f(l), f(a), f(t)];
        ins.extend(bc_toks(&bl));
        ins.extend(bc_toks(&br));
        let out = match ThermalGrid1D::new(n, l, a, t, bl, br) {
            Ok(g) => g1_toks(&g),
            Err(z) => vec![err(&z)],
        };
        m.row("numerical::ThermalGrid1D::new", Class::B, &ins, &out);
    }
    // fourier_number and max_stable_dt on every setup and malformed grid
    let mut grids = setups();
    grids.extend(malformed());
    for g in &grids {
        for dt in [1.0, 0.5, 1e-3, 0.0, -1.0, f64::NAN, f64::INFINITY] {
            let mut ins = g1_toks(g);
            ins.push(f(dt));
            m.row(
                "numerical::ThermalGrid1D::fourier_number",
                Class::B,
                &ins,
                &[f(g.fourier_number(dt))],
            );
        }
        m.row(
            "numerical::ThermalGrid1D::max_stable_dt",
            Class::B,
            &g1_toks(g),
            &[f(g.max_stable_dt())],
        );
    }
}

/// A step row: the grid and dt in, the grid after the step (and, for the adaptive step, the time advanced) out.
fn step_row(m: &mut Module, which: u8, g: &ThermalGrid1D, dt: f64) -> ThermalGrid1D {
    let mut after = g.clone();
    let (name, r) = match which {
        0 => (
            "numerical::ThermalGrid1D::step_explicit",
            after.step_explicit(dt).map(|()| None),
        ),
        1 => (
            "numerical::ThermalGrid1D::step_crank_nicolson",
            after.step_crank_nicolson(dt).map(|()| None),
        ),
        _ => (
            "numerical::ThermalGrid1D::step_adaptive",
            after.step_adaptive(dt).map(Some),
        ),
    };
    let out = match r {
        Ok(None) => g1_toks(&after),
        Ok(Some(t)) => {
            let mut v = vec![f(t)];
            v.extend(g1_toks(&after));
            v
        }
        Err(z) => vec![err(&z)],
    };
    let mut ins = g1_toks(g);
    ins.push(f(dt));
    m.row(name, Class::B, &ins, &out);
    after
}

fn step_rows(m: &mut Module) {
    for g0 in setups() {
        let dts = g0.max_stable_dt();
        // trajectories: 20 explicit steps at half the stable step, 20 Crank–Nicolson at 10×, 10 adaptive at 3.7×
        let mut g = g0.clone();
        for _ in 0..20 {
            g = step_row(m, 0, &g, 0.5 * dts);
        }
        let mut g = g0.clone();
        for _ in 0..20 {
            g = step_row(m, 1, &g, 10.0 * dts);
        }
        let mut g = g0.clone();
        for _ in 0..10 {
            g = step_row(m, 2, &g, 3.7 * dts);
        }
        // fractions and multiples of the stable step, and dt at the edges
        for k in [0.1, 1.0, 2.0, 100.0, 1e6] {
            // 10⁶ adaptive sub-steps only on the small grids (validate.py recomputes each in Python)
            if k > 1e5 && g0.nodes.len() > 10 {
                continue;
            }
            for which in 0..3 {
                step_row(m, which, &g0, k * dts);
            }
        }
        for e in edges() {
            for which in 0..3 {
                step_row(m, which, &g0, e);
            }
        }
        // the explicit stability limit, bisected, ± 2 ulps
        let unstable = |dt: f64| g0.clone().step_explicit(dt).is_err();
        for dt in around(first_true(dts, 2.0 * dts, &unstable), 2) {
            step_row(m, 0, &g0, dt);
        }
        // the adaptive step at exact multiples of the stable step ± 1 ulp, and its sub-step cap, bisected, ± 2 ulps
        for k in [1.0, 2.0, 3.0] {
            for dt in around(k * dts, 1) {
                step_row(m, 2, &g0, dt);
            }
        }
        let capped = |dt: f64| g0.clone().step_adaptive(dt).is_err();
        if g0.nodes.len() == 10 && !capped(1e5 * dts) && capped(1e7 * dts) {
            for dt in around(first_true(1e5 * dts, 1e7 * dts, &capped), 2) {
                step_row(m, 2, &g0, dt);
            }
        }
        // a huge Crank–Nicolson step: its pivots stop being finite
        step_row(m, 1, &g0, 1e300);
    }
    for g in malformed() {
        let dts = setups()[0].max_stable_dt();
        for which in 0..3 {
            step_row(m, which, &g, 0.5 * dts);
        }
    }
}

fn decode<T: for<'de> serde::Deserialize<'de>>(
    doc: &str,
    toks: impl Fn(&T) -> Vec<Tok>,
) -> Vec<Tok> {
    let mut de = serde_json::Deserializer::from_str(doc);
    match T::deserialize(&mut de).and_then(|x| de.end().map(|()| x)) {
        Ok(x) => toks(&x),
        Err(_) => vec![Tok::Err("Deserialize".to_string())],
    }
}

fn add_variants(docs: &mut Vec<String>, base: &str) {
    for (_kind, bytes) in variants(base) {
        if let Ok(text) = String::from_utf8(bytes)
            && !text.contains('\0')
        {
            docs.push(text);
        }
    }
}

fn json_rows(m: &mut Module) {
    // BoundaryCondition
    let bcs = [
        Bc::Fixed(400.0),
        Bc::Fixed(f64::NAN),
        Bc::Fixed(-0.0),
        Bc::Insulated,
        conv(25.0, 300.0, 50.0),
        conv(f64::INFINITY, f64::NAN, 0.0),
    ];
    let mut docs = Vec::new();
    for b in &bcs {
        let doc = serde_json::to_string(b).expect("serialize");
        m.row(
            "numerical::BoundaryCondition::to_json",
            Class::A,
            &bc_toks(b),
            &[s(doc.clone())],
        );
        docs.push(doc);
    }
    let base = docs[4].clone();
    add_variants(&mut docs, &base);
    docs.extend(
        [
            "\"Fixed\"",
            r#"{"Fixed":null}"#,
            r#"{"Fixed":"1"}"#,
            r#"{"Insulated":null}"#,
            r#"{"Insulated":1}"#,
            "\"insulated\"",
            r#"{"Convective":[1,2,3]}"#,
            r#"{"Convective":[1,2]}"#,
            r#"{"Convective":{"h":1,"t_fluid":2}}"#,
            r#"{"Convective":{"h":1,"t_fluid":2,"k":3,"x":4}}"#,
            r#"{"Fixed":1,"Insulated":null}"#,
            r#"{"Radiative":1}"#,
            "[]",
            "{}",
            "null",
        ]
        .map(String::from),
    );
    for doc in &docs {
        m.row(
            "numerical::BoundaryCondition::from_json",
            Class::A,
            &[s(doc.clone())],
            &decode::<Bc>(doc, bc_toks),
        );
    }
    // Side
    for sd in [Side::Left, Side::Right, Side::Top, Side::Bottom] {
        m.row(
            "numerical::Side::to_json",
            Class::A,
            &[Tok::E("Side", format!("{sd:?}"))],
            &[s(serde_json::to_string(&sd).expect("serialize"))],
        );
    }
    for doc in [
        "\"Left\"",
        "\"Bottom\"",
        "\"left\"",
        "\"Front\"",
        r#"{"Top":null}"#,
        "1",
        "",
    ] {
        m.row(
            "numerical::Side::from_json",
            Class::A,
            &[s(doc)],
            &decode::<Side>(doc, |x: &Side| vec![Tok::E("Side", format!("{x:?}"))]),
        );
    }
    // ThermalGrid1D
    let mut grids = setups();
    grids.truncate(4);
    let mut odd = grids[0].clone();
    odd.nodes[2] = f64::NAN;
    odd.nodes.truncate(3);
    grids.push(odd);
    let mut empty = grids[1].clone();
    empty.nodes.clear();
    grids.push(empty);
    let mut gdocs = Vec::new();
    for g in &grids {
        let doc = serde_json::to_string(g).expect("serialize");
        m.row(
            "numerical::ThermalGrid1D::to_json",
            Class::A,
            &g1_toks(g),
            &[s(doc.clone())],
        );
        gdocs.push(doc);
    }
    let base = gdocs[0].clone();
    add_variants(&mut gdocs, &base);
    gdocs.extend(
        [
            r#"[[1,2,3],0.5,0.0001,"Insulated",{"Fixed":1}]"#,
            r#"{"nodes":[1,2,3],"dx":0.5,"alpha":0.0001,"bc_left":"Insulated"}"#,
            r#"{"nodes":3,"dx":0.5,"alpha":0.0001,"bc_left":"Insulated","bc_right":"Insulated"}"#,
            r#"{"nodes":[1,"2",3],"dx":0.5,"alpha":0.0001,"bc_left":"Insulated","bc_right":"Insulated"}"#,
            r#"{"nodes":[1,2,3],"dx":0.5,"alpha":0.0001,"bc_left":"Open","bc_right":"Insulated"}"#,
            r#"{"nodes":[],"dx":-1,"alpha":0,"bc_left":{"Fixed":1},"bc_right":{"Convective":[1,2,3]}}"#,
            "[]",
            "null",
        ]
        .map(String::from),
    );
    for doc in &gdocs {
        m.row(
            "numerical::ThermalGrid1D::from_json",
            Class::A,
            &[s(doc.clone())],
            &decode::<ThermalGrid1D>(doc, g1_toks),
        );
    }
}

pub fn emit() -> Module {
    let mut m = Module::new("numerical");
    grid_rows(&mut m);
    step_rows(&mut m);
    json_rows(&mut m);
    m
}
