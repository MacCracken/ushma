//! `chem.tsv`'s function and JSON rows (M5a), added to the table module's preset rows.
//!
//! A reaction side is written `v:n` and, per term, its coefficient then its species by value (`s:name`, ΔH°f, ΔG°f,
//! C°p), so edge species (non-finite, negative, zero C°p) sit beside the presets. Rows: Hess's law and the Gibbs sum over
//! named reactions, their reverses, empty sides, seeded random reactions, and each coefficient and species field of a
//! base reaction at the edges; ΔH − TΔS, ln K and Van 't Hoff seeded and at the edges, ΔG° ±0 among them; K and
//! Van 't Hoff where exp leaves the f64 range (found by bisection on the bit patterns) ± 2 ulps; the flame temperature
//! over named combustions and initial temperatures, every coefficient and species at the edges (reactants checked before
//! products), an empty product side, and where T_flame crosses 0 (bisected) ± 2 ulps; `Species`'s JSON. Classes: B for
//! the sums and ln K, C where exp enters, A for JSON.

use std::borrow::Cow;

use ushma::chem::{self as ch, Species};

use crate::grids::Rng;
use crate::json_cases::variants;
use crate::sweep::{around, edges, res, row, tuples};
use crate::tsv::{Class, Module, Tok, f, s};

fn sp_toks(x: &Species) -> Vec<Tok> {
    vec![s(x.name.to_string()), f(x.delta_hf), f(x.delta_gf), f(x.cp)]
}

fn side(terms: &[(f64, &Species)]) -> Vec<Tok> {
    let mut v = vec![Tok::V(terms.len())];
    for (n, x) in terms {
        v.push(f(*n));
        v.extend(sp_toks(x));
    }
    v
}

fn species(name: &str, hf: f64, gf: f64, cp: f64) -> Species {
    Species {
        name: Cow::Owned(name.to_string()),
        delta_hf: hf,
        delta_gf: gf,
        cp,
    }
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

type Side<'a> = Vec<(f64, &'a Species)>;

/// Named reactions as (reactants, products).
fn reactions() -> Vec<(Side<'static>, Side<'static>)> {
    vec![
        (
            vec![(1.0, &ch::CH4), (2.0, &ch::O2)],
            vec![(1.0, &ch::CO2), (2.0, &ch::H2O_GAS)],
        ),
        (
            vec![(1.0, &ch::H2), (0.5, &ch::O2)],
            vec![(1.0, &ch::H2O_LIQUID)],
        ),
        (
            vec![(1.0, &ch::C3H8), (5.0, &ch::O2)],
            vec![(3.0, &ch::CO2), (4.0, &ch::H2O_GAS)],
        ),
        (vec![(1.0, &ch::CO), (0.5, &ch::O2)], vec![(1.0, &ch::CO2)]),
        (vec![(1.0, &ch::N2), (3.0, &ch::H2)], vec![(2.0, &ch::NH3)]),
        (vec![(1.0, &ch::N2), (1.0, &ch::O2)], vec![(2.0, &ch::NO)]),
        (
            vec![(1.0, &ch::C2H6), (3.5, &ch::O2)],
            vec![(2.0, &ch::CO2), (3.0, &ch::H2O_LIQUID)],
        ),
        (vec![], vec![(1.0, &ch::CO2)]),
        (vec![(1.0, &ch::H2O_GAS)], vec![]),
        (vec![], vec![]),
        (vec![(1.0, &ch::H2), (1.0, &ch::O2)], vec![(1.0, &ch::N2)]),
    ]
}

fn hess_rows(m: &mut Module) {
    let edge_species: Vec<Species> = vec![
        species("nan", f64::NAN, f64::NAN, f64::NAN),
        species("inf", f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY),
        species("negzero", -0.0, -0.0, -0.0),
        species("huge", f64::MAX, -f64::MAX, f64::MAX),
        species("tiny", f64::from_bits(1), -f64::from_bits(1), 1e-300),
    ];
    let mut cases: Vec<(Side, Side)> = Vec::new();
    for (r, p) in reactions() {
        cases.push((r.clone(), p.clone()));
        cases.push((p, r));
    }
    // each coefficient of methane's combustion, and each of its species, at the edges
    let (base_r, base_p) = reactions().remove(0);
    for k in 0..base_r.len() + base_p.len() {
        for e in edges() {
            let (mut r, mut p) = (base_r.clone(), base_p.clone());
            if k < r.len() {
                r[k].0 = e;
            } else {
                p[k - r.len()].0 = e;
            }
            cases.push((r, p));
        }
        for x in &edge_species {
            let (mut r, mut p) = (base_r.clone(), base_p.clone());
            if k < r.len() {
                r[k].1 = x;
            } else {
                p[k - r.len()].1 = x;
            }
            cases.push((r, p));
        }
    }
    // seeded random reactions over the presets
    let mut rng = Rng::new(100);
    let draw_side = |rng: &mut Rng| -> Side<'static> {
        let n = (rng.unit() * 5.0) as usize;
        (0..n)
            .map(|_| {
                let i = ((rng.unit() * 12.0) as usize).min(11);
                (rng.uniform(0.0, 10.0), ch::ALL_SPECIES[i])
            })
            .collect()
    };
    for _ in 0..200 {
        let r = draw_side(&mut rng);
        let p = draw_side(&mut rng);
        cases.push((r, p));
    }
    for (r, p) in &cases {
        let mut ins = side(p);
        ins.extend(side(r));
        m.row(
            "chem::reaction_enthalpy",
            Class::B,
            &ins,
            &[f(ch::reaction_enthalpy(p, r))],
        );
        m.row(
            "chem::reaction_gibbs",
            Class::B,
            &ins,
            &[f(ch::reaction_gibbs(p, r))],
        );
    }
}

fn scalar_rows(m: &mut Module) {
    // ΔG = ΔH − TΔS
    for a in tuples(&[vec![-100_000.0, -200.0, 298.15]], 200, 101, &|r| {
        vec![
            r.uniform(-3e6, 3e6),
            r.uniform(-1e3, 1e3),
            r.uniform(1.0, 3000.0),
        ]
    }) {
        row(
            m,
            "chem::gibbs_from_enthalpy_entropy",
            Class::B,
            &a,
            vec![f(ch::gibbs_from_enthalpy_entropy(a[0], a[1], a[2]))],
        );
    }
    // ln K and K: seeded, the edges, ΔG° ±0, and where K leaves the f64 range
    let mut lk = tuples(
        &[
            vec![-50_000.0, 298.15],
            vec![0.0, 298.15],
            vec![-0.0, 298.15],
        ],
        200,
        102,
        &|r| vec![r.uniform(-3e6, 3e6), r.uniform(1.0, 3000.0)],
    );
    for &t in &[298.15, 1000.0, 10.0] {
        // overflow (ΔG° < 0) and underflow (ΔG° > 0): the switch from Ok to Err, ± 2 ulps
        let over = first_true(1.0, 1e9, &|x| ch::equilibrium_constant(-x, t).is_err());
        let under = first_true(1.0, 1e9, &|x| ch::equilibrium_constant(x, t).is_err());
        for x in around(over, 2) {
            lk.push(vec![-x, t]);
        }
        for x in around(under, 2) {
            lk.push(vec![x, t]);
        }
    }
    lk.push(vec![-2_073_900.0, 298.15]);
    lk.push(vec![2_073_900.0, 298.15]);
    for a in &lk {
        row(
            m,
            "chem::ln_equilibrium_constant",
            Class::B,
            a,
            res(ch::ln_equilibrium_constant(a[0], a[1])),
        );
        row(
            m,
            "chem::equilibrium_constant",
            Class::C,
            a,
            res(ch::equilibrium_constant(a[0], a[1])),
        );
    }
    // Van 't Hoff: seeded, the edges, and where K₂ overflows
    let mut vh = tuples(&[vec![100.0, -50_000.0, 298.15, 400.0]], 300, 103, &|r| {
        vec![
            r.log_uniform(-10.0, 10.0),
            r.uniform(-3e5, 3e5),
            r.uniform(100.0, 3000.0),
            r.uniform(100.0, 3000.0),
        ]
    });
    let over = first_true(1.0, 1e12, &|dh| {
        ch::vant_hoff_k(1.0, dh, 298.15, 400.0).is_err()
    });
    for x in around(over, 2) {
        vh.push(vec![1.0, x, 298.15, 400.0]);
    }
    vh.push(vec![100.0, 0.0, 298.15, 400.0]);
    vh.push(vec![100.0, -0.0, 298.15, 400.0]);
    vh.push(vec![100.0, -50_000.0, 298.15, 298.15]);
    for a in vh {
        row(
            m,
            "chem::vant_hoff_k",
            Class::C,
            &a,
            res(ch::vant_hoff_k(a[0], a[1], a[2], a[3])),
        );
    }
}

fn flame_rows(m: &mut Module) {
    let air = |o2: f64| vec![(1.0, &ch::CH4), (o2, &ch::O2), (3.76 * o2, &ch::N2)];
    let methane_r = air(2.0);
    let methane_p = vec![(1.0, &ch::CO2), (2.0, &ch::H2O_GAS), (7.52, &ch::N2)];
    let mut cases: Vec<(Side, Side, f64)> = Vec::new();
    let named: Vec<(Side, Side)> = vec![
        (methane_r.clone(), methane_p.clone()),
        (
            air(3.0),
            vec![
                (1.0, &ch::CO2),
                (2.0, &ch::H2O_GAS),
                (1.0, &ch::O2),
                (11.28, &ch::N2),
            ],
        ),
        (
            vec![(1.0, &ch::H2), (0.5, &ch::O2), (1.88, &ch::N2)],
            vec![(1.0, &ch::H2O_GAS), (1.88, &ch::N2)],
        ),
        (
            vec![(1.0, &ch::C3H8), (5.0, &ch::O2), (18.8, &ch::N2)],
            vec![(3.0, &ch::CO2), (4.0, &ch::H2O_GAS), (18.8, &ch::N2)],
        ),
        (vec![(1.0, &ch::CO), (0.5, &ch::O2)], vec![(1.0, &ch::CO2)]),
        (vec![(1.0, &ch::CH4)], vec![(1.0, &ch::CO2)]),
        (methane_p.clone(), methane_r.clone()),
        (methane_r.clone(), vec![]),
        (vec![], methane_p.clone()),
    ];
    for (r, p) in &named {
        for &t in &[298.15, 200.0, 500.0, 1000.0, 2000.0] {
            cases.push((r.clone(), p.clone(), t));
        }
        for e in edges() {
            cases.push((r.clone(), p.clone(), e));
        }
    }
    // every coefficient of methane-air at the edges (the reactants are checked first), and edge species
    let zero_cp = species("zero_cp", -100.0, -50.0, 0.0);
    let neg_cp = species("neg_cp", -100.0, -50.0, -40.0);
    let nan_cp = species("nan_cp", -100.0, -50.0, f64::NAN);
    let inf_h = species("inf_h", f64::NEG_INFINITY, 0.0, 30.0);
    for k in 0..methane_r.len() + methane_p.len() {
        for e in edges() {
            let (mut r, mut p) = (methane_r.clone(), methane_p.clone());
            if k < r.len() {
                r[k].0 = e;
            } else {
                p[k - r.len()].0 = e;
            }
            cases.push((r, p, 298.15));
        }
    }
    for x in [&zero_cp, &neg_cp, &nan_cp, &inf_h] {
        cases.push((methane_r.clone(), vec![(1.0, x)], 298.15));
        cases.push((vec![(1.0, x)], methane_p.clone(), 298.15));
    }
    // a bad coefficient on both sides: the reactants' is the one reported (the code is the same)
    cases.push((vec![(-1.0, &ch::CH4)], vec![(f64::NAN, &ch::CO2)], 298.15));
    // where T_flame crosses 0 for the reverse of methane's combustion, ± 2 ulps
    let (rr, rp) = (methane_p.clone(), methane_r.clone());
    let ok = |t: f64| ch::adiabatic_flame_temperature(&rr, &rp, t).is_ok();
    if !ok(298.15) && ok(1e5) {
        for x in around(first_true(298.15, 1e5, &ok), 2) {
            cases.push((rr.clone(), rp.clone(), x));
        }
    }
    // seeded random combustions over the presets
    let mut rng = Rng::new(104);
    for _ in 0..200 {
        let n_r = 1 + (rng.unit() * 3.0) as usize;
        let n_p = 1 + (rng.unit() * 3.0) as usize;
        let pick = |rng: &mut Rng| ch::ALL_SPECIES[((rng.unit() * 12.0) as usize).min(11)];
        let r: Side = (0..n_r)
            .map(|_| (rng.uniform(0.0, 5.0), pick(&mut rng)))
            .collect();
        let p: Side = (0..n_p)
            .map(|_| (rng.uniform(0.0, 5.0), pick(&mut rng)))
            .collect();
        let t = rng.uniform(200.0, 1500.0);
        cases.push((r, p, t));
    }
    for (r, p, t) in &cases {
        let mut ins = side(r);
        ins.extend(side(p));
        ins.push(f(*t));
        m.row(
            "chem::adiabatic_flame_temperature",
            Class::B,
            &ins,
            &res(ch::adiabatic_flame_temperature(r, p, *t)),
        );
    }
}

fn decoded(doc: &str) -> Vec<Tok> {
    let mut de = serde_json::Deserializer::from_str(doc);
    match Species::deserialize(&mut de).and_then(|x| de.end().map(|()| x)) {
        Ok(x) => sp_toks(&x),
        Err(_) => vec![Tok::Err("Deserialize".to_string())],
    }
}

use serde::Deserialize;

fn json_rows(m: &mut Module) {
    let mut values: Vec<Species> = ch::ALL_SPECIES.iter().map(|x| (*x).clone()).collect();
    for name in [
        "",
        "a\"b\\c/d",
        "tab\tnew\nline\r",
        "ctl\u{1}\u{1f}\u{7f}",
        "\u{e9}\u{1f600} ok",
        "\u{0}nul",
    ] {
        values.push(species(name, 1.0, 2.0, 3.0));
    }
    values.push(species("specials", f64::NAN, f64::INFINITY, -0.0));
    values.push(species(
        "extremes",
        f64::MAX,
        f64::MIN_POSITIVE,
        f64::from_bits(1),
    ));
    values.push(species("forms", 1e16, 1e-5, 0.1 + 0.2));
    let mut docs = Vec::new();
    for x in &values {
        let doc = serde_json::to_string(x).expect("serialize");
        m.row(
            "chem::Species::to_json",
            Class::A,
            &sp_toks(x),
            &[s(doc.clone())],
        );
        docs.push(doc);
    }
    let base = docs[7].clone();
    for (_kind, bytes) in variants(&base) {
        if let Ok(text) = String::from_utf8(bytes)
            && !text.contains('\0')
        {
            docs.push(text);
        }
    }
    for doc in [
        r#"{"name":"CH\u0000x","delta_hf":1,"delta_gf":2.5,"cp":-0.0}"#,
        r#"["CH4",-74810,-50720,35.31]"#,
        r#"["CH4",-74810,-50720]"#,
        r#"{"cp":3,"delta_gf":2,"delta_hf":1,"name":"rev"}"#,
        r#"{"name":"n","delta_hf":1,"delta_gf":2}"#,
        r#"{"name":null,"delta_hf":1,"delta_gf":2,"cp":3}"#,
        r#"{"name":"n","delta_hf":[1],"delta_gf":2,"cp":3}"#,
        r#"{"name":"n","delta_hf":1e400,"delta_gf":2,"cp":3}"#,
        "[]",
        "null",
        "\"CH4\"",
        "",
    ] {
        docs.push(doc.to_string());
    }
    for doc in &docs {
        m.row(
            "chem::Species::from_json",
            Class::A,
            &[s(doc.clone())],
            &decoded(doc),
        );
    }
}

/// Adds the module's function and JSON rows to `m` (the `chem` module, after its preset rows).
pub fn rows(m: &mut Module) {
    hess_rows(m);
    scalar_rows(m);
    flame_rows(m);
    json_rows(m);
}
