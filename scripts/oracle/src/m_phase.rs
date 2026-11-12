//! `phase.tsv`'s function and JSON rows (M3a), added to the module the table presets seed (`t_tables.rs`).
//!
//! A row's substance is its eleven fields in declaration order (`s:name`, then triple_point_t … molar_mass), placed
//! where Rust's signature puts it (first for the `SubstanceData` methods and `heat_of_*`, last for
//! `heat_for_phase_change`). The substances are the 8 presets and edge substances (each property at 0, −1, NaN, ±∞;
//! a boiling point outside the triple–critical span; a sublimer whose boiling point is above its triple point; zero
//! molar mass). Every join is crossed: `phase_at` runs each substance at temperatures across the triple point, the
//! melting and boiling points and the critical point (each ± 1 ulp), and at each temperature at the curve pressure
//! Rust computes there ± 1 ulp, the critical pressure, 1 atm and fixed pressures; the curves run their own ranges and
//! past them. Classes: A for the bools, the phases and JSON; B for + − × ÷; C where `ln` or `exp` enters.

use std::borrow::Cow;

use serde::Deserialize as _;
use ushma::phase::{self as ph, ALL_SUBSTANCES, Phase, SubstanceData};
use ushma::state::ATM;

use crate::grids::Rng;
use crate::json_cases::variants;
use crate::sweep::{around, edges, res, tuples};
use crate::tsv::{Class, Module, Tok, b, e, err, f, s};

fn fields(x: &SubstanceData) -> Vec<Tok> {
    vec![
        s(x.name.to_string()),
        f(x.triple_point_t),
        f(x.triple_point_p),
        f(x.critical_t),
        f(x.critical_p),
        f(x.critical_density),
        f(x.latent_heat_fusion),
        f(x.latent_heat_vaporization),
        f(x.melting_point),
        f(x.boiling_point),
        f(x.molar_mass),
    ]
}

fn with(base: &SubstanceData, name: &str, k: usize, v: f64) -> SubstanceData {
    let mut x = base.clone();
    x.name = Cow::Owned(name.to_string());
    match k {
        0 => x.triple_point_t = v,
        1 => x.triple_point_p = v,
        2 => x.critical_t = v,
        3 => x.critical_p = v,
        4 => x.critical_density = v,
        5 => x.latent_heat_fusion = v,
        6 => x.latent_heat_vaporization = v,
        7 => x.melting_point = v,
        8 => x.boiling_point = v,
        _ => x.molar_mass = v,
    }
    x
}

/// The presets, then each property of Water and of CO₂ at the edge values, then the shapes the branches turn on.
fn substances() -> Vec<SubstanceData> {
    let mut v: Vec<SubstanceData> = ALL_SUBSTANCES.iter().map(|x| (*x).clone()).collect();
    let water = ph::WATER_PHASE.clone();
    let co2 = ph::CARBON_DIOXIDE.clone();
    for k in 0..10 {
        for x in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            v.push(with(&water, "edge-water", k, x));
            v.push(with(&co2, "edge-co2", k, x));
        }
    }
    v.push(with(&water, "boil-above-critical", 8, 700.0));
    v.push(with(&water, "boil-below-triple", 8, 200.0));
    v.push(with(&co2, "sublimer-boil-above-triple", 8, 250.0));
    v.push(with(&co2, "sublimer-boil-at-triple", 8, co2.triple_point_t));
    v.push(with(&water, "triple-at-1atm", 1, ATM));
    v.push(with(&water, "triple-above-1atm", 1, ATM * 2.0));
    v
}

fn sub_rows(
    m: &mut Module,
    name: &str,
    class: Class,
    x: &SubstanceData,
    extra: &[f64],
    out: Vec<Tok>,
) {
    let mut ins = fields(x);
    ins.extend(extra.iter().map(|&a| f(a)));
    m.row(name, class, &ins, &out);
}

/// The temperatures across a substance's joins: the triple, melting, boiling and critical points ± 1 ulp, points
/// between them, and below and above.
fn joins(x: &SubstanceData) -> Vec<f64> {
    let mut t = Vec::new();
    for j in [
        x.triple_point_t,
        x.melting_point,
        x.boiling_point,
        x.critical_t,
    ] {
        if j.is_finite() && j > 0.0 {
            t.extend(around(j, 1));
        }
    }
    let lo = x.triple_point_t;
    let hi = x.critical_t;
    if lo.is_finite() && hi.is_finite() && lo > 0.0 && hi > lo {
        for k in 1..8 {
            t.push(lo + (hi - lo) * f64::from(k) / 8.0);
        }
        t.extend([lo * 0.5, lo * 0.9, hi * 1.1, hi * 2.0]);
    }
    t
}

fn methods(m: &mut Module) {
    let subs = substances();
    for x in &subs {
        sub_rows(
            m,
            "phase::SubstanceData::sublimes_at_atmospheric",
            Class::A,
            x,
            &[],
            vec![b(x.sublimes_at_atmospheric())],
        );
        sub_rows(
            m,
            "phase::SubstanceData::latent_heat_sublimation",
            Class::C,
            x,
            &[],
            vec![f(x.latent_heat_sublimation())],
        );
        let mut ts = joins(x);
        ts.extend(edges());
        for &t in &ts {
            sub_rows(
                m,
                "phase::SubstanceData::vapor_pressure",
                Class::C,
                x,
                &[t],
                res(x.vapor_pressure(t)),
            );
            sub_rows(
                m,
                "phase::SubstanceData::sublimation_pressure",
                Class::C,
                x,
                &[t],
                res(x.sublimation_pressure(t)),
            );
        }
    }
    // phase_at(T, P): the presets on the full grid, the edge substances on a light one
    for (i, x) in subs.iter().enumerate() {
        let mut ts = joins(x);
        if i >= ALL_SUBSTANCES.len() {
            ts.truncate(6);
            ts.push(300.0);
        }
        let mut pts: Vec<(f64, f64)> = Vec::new();
        for &t in &ts {
            let mut ps = vec![ATM, 1.0, 1e3, 1e5, 1e7, 1e9];
            if x.critical_p.is_finite() && x.critical_p > 0.0 {
                ps.extend(around(x.critical_p, 1));
            }
            for curve in [x.vapor_pressure(t), x.sublimation_pressure(t)] {
                if let Ok(p) = curve
                    && p.is_finite()
                    && p > 0.0
                {
                    ps.extend(around(p, 1));
                }
            }
            for p in ps {
                pts.push((t, p));
            }
        }
        if i < ALL_SUBSTANCES.len() {
            for e in edges() {
                pts.push((e, ATM));
                pts.push((300.0, e));
            }
        }
        for (t, p) in pts {
            let out = match x.phase_at(t, p) {
                Ok(v) => vec![e("Phase", v)],
                Err(z) => vec![err(&z)],
            };
            sub_rows(
                m,
                "phase::SubstanceData::phase_at",
                Class::C,
                x,
                &[t, p],
                out,
            );
        }
    }
}

fn free(m: &mut Module) {
    let subs = substances();
    let masses = [
        1.0,
        0.0,
        -0.0,
        2.5,
        1e-300,
        1e300,
        -1.0,
        f64::NAN,
        f64::INFINITY,
    ];
    for x in &subs {
        for &mass in &masses {
            sub_rows(
                m,
                "phase::heat_of_fusion",
                Class::B,
                x,
                &[mass],
                vec![f(ph::heat_of_fusion(x, mass))],
            );
            sub_rows(
                m,
                "phase::heat_of_vaporization",
                Class::B,
                x,
                &[mass],
                vec![f(ph::heat_of_vaporization(x, mass))],
            );
        }
    }
    // heat_for_phase_change(m, c_p, T_start, T_end, substance): ranges across the melting and boiling points, both
    // ways, starting and ending on them; then each argument at the edges
    for (i, x) in subs.iter().enumerate() {
        let mut ts = vec![10.0, 100.0, 250.0, 300.0, 400.0, 1000.0];
        for j in [x.melting_point, x.boiling_point] {
            if j.is_finite() && j > 0.0 {
                ts.extend(around(j, 1));
            }
        }
        let n = if i < ALL_SUBSTANCES.len() {
            ts.len()
        } else {
            4
        };
        for &a in &ts[..n] {
            for &c in &ts {
                let mut ins = vec![f(1.0), f(2000.0), f(a), f(c)];
                ins.extend(fields(x));
                m.row(
                    "phase::heat_for_phase_change",
                    Class::C,
                    &ins,
                    &res(ph::heat_for_phase_change(1.0, 2000.0, a, c, x)),
                );
            }
        }
    }
    let water = &ph::WATER_PHASE;
    for k in 0..4 {
        for e in edges() {
            let mut a = [1.0, 4186.0, 250.0, 400.0];
            a[k] = e;
            let mut ins: Vec<Tok> = a.iter().map(|&v| f(v)).collect();
            ins.extend(fields(water));
            m.row(
                "phase::heat_for_phase_change",
                Class::C,
                &ins,
                &res(ph::heat_for_phase_change(a[0], a[1], a[2], a[3], water)),
            );
        }
    }
    // clausius_clapeyron_slope(L, T, Δv) and clausius_clapeyron_pressure(P₁, T₁, T₂, L)
    for a in tuples(&[vec![2.26e6, 373.15, 1.672]], 300, 70, &|r: &mut Rng| {
        vec![
            r.log_uniform(3.0, 7.0),
            r.uniform(1.0, 2000.0),
            r.signed_log_uniform(-4.0, 3.0),
        ]
    }) {
        crate::sweep::row(
            m,
            "phase::clausius_clapeyron_slope",
            Class::B,
            &a,
            res(ph::clausius_clapeyron_slope(a[0], a[1], a[2])),
        );
    }
    let mut cc = tuples(
        &[vec![101_325.0, 373.15, 393.15, 40_650.0]],
        400,
        71,
        &|r: &mut Rng| {
            vec![
                r.log_uniform(1.0, 8.0),
                r.uniform(1.0, 2000.0),
                r.uniform(1.0, 2000.0),
                r.log_uniform(2.0, 6.0),
            ]
        },
    );
    cc.extend([
        vec![101_325.0, 373.15, 1.0, 40_650.0],
        vec![101_325.0, 1.0, 1e300, 1e9],
        vec![1e300, 300.0, 1e6, 1e6],
        vec![101_325.0, 373.15, 373.15, 40_650.0],
    ]);
    for a in cc {
        crate::sweep::row(
            m,
            "phase::clausius_clapeyron_pressure",
            Class::C,
            &a,
            res(ph::clausius_clapeyron_pressure(a[0], a[1], a[2], a[3])),
        );
    }
}

fn substance_values() -> Vec<SubstanceData> {
    let mut v: Vec<SubstanceData> = ALL_SUBSTANCES.iter().map(|x| (*x).clone()).collect();
    let water = ph::WATER_PHASE.clone();
    for name in [
        "",
        "a\"b\\c/d",
        "tab\tnew\nline\r",
        "ctl\u{1}\u{1f}\u{7f}",
        "\u{e9}\u{1f600} ok",
        "\u{0}nul",
    ] {
        let mut x = water.clone();
        x.name = Cow::Owned(name.to_string());
        v.push(x);
    }
    let mut sp = with(&water, "specials", 0, f64::NAN);
    sp.triple_point_p = f64::INFINITY;
    sp.critical_t = f64::NEG_INFINITY;
    sp.critical_p = -0.0;
    sp.molar_mass = f64::MAX;
    sp.boiling_point = f64::from_bits(1);
    v.push(sp);
    v
}

fn json_rows(m: &mut Module) {
    // Phase
    for p in [
        Phase::Solid,
        Phase::Liquid,
        Phase::Gas,
        Phase::SupercriticalFluid,
    ] {
        m.row(
            "phase::Phase::to_json",
            Class::A,
            &[e("Phase", p)],
            &[s(serde_json::to_string(&p).expect("serialize"))],
        );
    }
    for doc in [
        "\"Solid\"",
        "\"Liquid\"",
        "\"Gas\"",
        "\"SupercriticalFluid\"",
        "\"liquid\"",
        "\"Plasma\"",
        "\"Liquid \"",
        "\"\\u004ciquid\"",
        "{\"Liquid\":null}",
        "{\"Liquid\":1}",
        "{\"Liquid\":[]}",
        "{\"Liquid\":null,\"Gas\":null}",
        "{}",
        "[\"Liquid\"]",
        "0",
        "1",
        "null",
        "true",
        "\"\"",
        "",
        " \"Gas\" ",
        "\"Gas\"x",
    ] {
        let out = match serde_json::from_str::<Phase>(doc) {
            Ok(p) => vec![e("Phase", p)],
            Err(_) => vec![Tok::Err("Deserialize".to_string())],
        };
        m.row("phase::Phase::from_json", Class::A, &[s(doc)], &out);
    }
    // SubstanceData
    let mut docs = Vec::new();
    for x in substance_values() {
        let doc = serde_json::to_string(&x).expect("serialize");
        m.row(
            "phase::SubstanceData::to_json",
            Class::A,
            &fields(&x),
            &[s(doc.clone())],
        );
        docs.push(doc);
    }
    let water = docs[0].clone();
    for (_kind, bytes) in variants(&water) {
        if let Ok(text) = String::from_utf8(bytes)
            && !text.contains('\0')
        {
            docs.push(text);
        }
    }
    docs.extend(
        [
            r#"["Water",273.16,611.657,647.096,22064000,322,334000,2260000,273.15,373.15,0.018015]"#,
            r#"["Water",273.16,611.657,647.096,22064000,322,334000,2260000,273.15,373.15]"#,
            r#"{"name":"n","triple_point_t":1,"triple_point_p":2,"critical_t":3,"critical_p":4,"critical_density":5,"latent_heat_fusion":6,"latent_heat_vaporization":7,"melting_point":8,"boiling_point":9}"#,
            r#"{"name":"n","triple_point_t":1e400,"triple_point_p":2,"critical_t":3,"critical_p":4,"critical_density":5,"latent_heat_fusion":6,"latent_heat_vaporization":7,"melting_point":8,"boiling_point":9,"molar_mass":10}"#,
            "[]",
            "null",
            "",
        ]
        .map(String::from),
    );
    for doc in &docs {
        let mut de = serde_json::Deserializer::from_str(doc);
        let out = match SubstanceData::deserialize(&mut de).and_then(|x| de.end().map(|()| x)) {
            Ok(x) => fields(&x),
            Err(_) => vec![Tok::Err("Deserialize".to_string())],
        };
        m.row(
            "phase::SubstanceData::from_json",
            Class::A,
            &[s(doc.clone())],
            &out,
        );
    }
}

/// Adds the module's function and JSON rows to `m` (the `phase` module).
pub fn rows(m: &mut Module) {
    methods(m);
    free(m);
    json_rows(m);
}
