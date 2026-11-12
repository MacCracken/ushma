//! `state.tsv`'s function and JSON rows (M2a onward), added to the module the table presets seed (`t_tables.rs`).
//!
//! For each function: a seeded sweep over its arguments' physical ranges, each argument of each base case at every
//! edge value (`sweep::edges`), then the function's own switches (V − nb across 0 for van der Waals, γ just above 1,
//! V₁ = V₂, nRT at the 1e-30 divisor threshold). The reduced coordinates run every gas over a grid of T and P. Classes:
//! B for + − × ÷, C where `ln` or `powf` enters (`isothermal_work`, `adiabatic_temperature`).
//!
//! `GasData`'s JSON (class A), as `material.tsv`'s: serde_json's document for the presets and for names and values
//! that exercise its escapes and number forms, and serde_json's verdict on each of those documents, on the variants of
//! Nitrogen's (`json_cases::variants`) and on hand-written edges. Each row's input is the four fields in declaration
//! order (`s:name`, then critical_t, critical_p, acentric_factor), or for `from_json` the document.

use std::borrow::Cow;

use serde::Deserialize as _;
use ushma::state::{self as st, ALL_GASES, GasData};

use crate::json_cases::variants;
use crate::sweep::{around, res, row, tuples};
use crate::tsv::{Class, Module, Tok, f, s};

fn part1(m: &mut Module) {
    // ideal_gas_pressure(n, T, V)
    for a in tuples(
        &[vec![1.0, 273.15, 0.02241], vec![2.0, 300.0, 0.01]],
        400,
        40,
        &|r| {
            vec![
                r.log_uniform(-3.0, 3.0),
                r.uniform(0.0, 3000.0),
                r.log_uniform(-5.0, 2.0),
            ]
        },
    ) {
        row(
            m,
            "state::ideal_gas_pressure",
            Class::B,
            &a,
            res(st::ideal_gas_pressure(a[0], a[1], a[2])),
        );
    }
    // ideal_gas_volume(n, T, P)
    for a in tuples(&[vec![1.0, 273.15, 101_325.0]], 400, 41, &|r| {
        vec![
            r.log_uniform(-3.0, 3.0),
            r.uniform(0.0, 3000.0),
            r.log_uniform(2.0, 8.0),
        ]
    }) {
        row(
            m,
            "state::ideal_gas_volume",
            Class::B,
            &a,
            res(st::ideal_gas_volume(a[0], a[1], a[2])),
        );
    }
    // ideal_gas_temperature(P, V, n)
    for a in tuples(&[vec![101_325.0, 0.02241, 1.0]], 400, 42, &|r| {
        vec![
            r.log_uniform(2.0, 8.0),
            r.log_uniform(-5.0, 2.0),
            r.log_uniform(-3.0, 3.0),
        ]
    }) {
        row(
            m,
            "state::ideal_gas_temperature",
            Class::B,
            &a,
            res(st::ideal_gas_temperature(a[0], a[1], a[2])),
        );
    }
    // van_der_waals_pressure(n, T, V, a, b), with V − nb across 0
    let mut vdw = tuples(
        &[
            vec![1.0, 273.15, 0.02241, 0.3658, 4.286e-5],
            vec![1.0, 300.0, 0.0224, 0.1364, 3.219e-5],
        ],
        500,
        43,
        &|r| {
            vec![
                r.log_uniform(-2.0, 2.0),
                r.uniform(1.0, 2000.0),
                r.log_uniform(-4.0, 1.0),
                r.signed_log_uniform(-3.0, 1.0),
                r.log_uniform(-6.0, -3.0),
            ]
        },
    );
    for v in around(4.286e-5, 3) {
        vdw.push(vec![1.0, 300.0, v, 0.3658, 4.286e-5]);
    }
    vdw.push(vec![0.0, 300.0, 1e-6, 0.3658, 4.286e-5]);
    vdw.push(vec![1.0, 300.0, 1e-6, -0.3658, 0.0]);
    for a in vdw {
        row(
            m,
            "state::van_der_waals_pressure",
            Class::B,
            &a,
            res(st::van_der_waals_pressure(a[0], a[1], a[2], a[3], a[4])),
        );
    }
    // isothermal_work(n, T, V₁, V₂), with V₁ = V₂
    let mut iso = tuples(&[vec![1.0, 300.0, 0.01, 0.02]], 500, 44, &|r| {
        vec![
            r.log_uniform(-2.0, 2.0),
            r.uniform(0.0, 2000.0),
            r.log_uniform(-5.0, 1.0),
            r.log_uniform(-5.0, 1.0),
        ]
    });
    for v in around(0.01, 2) {
        iso.push(vec![1.0, 300.0, 0.01, v]);
    }
    for a in iso {
        row(
            m,
            "state::isothermal_work",
            Class::C,
            &a,
            res(st::isothermal_work(a[0], a[1], a[2], a[3])),
        );
    }
    // isobaric_work(P, V₁, V₂)
    for a in tuples(&[vec![101_325.0, 0.01, 0.02]], 300, 45, &|r| {
        vec![
            r.log_uniform(2.0, 8.0),
            r.log_uniform(-5.0, 1.0),
            r.log_uniform(-5.0, 1.0),
        ]
    }) {
        row(
            m,
            "state::isobaric_work",
            Class::B,
            &a,
            vec![f(st::isobaric_work(a[0], a[1], a[2]))],
        );
    }
    // adiabatic_temperature(T₁, V₁, V₂, γ), with γ just above 1
    let mut adi = tuples(&[vec![300.0, 0.02, 0.01, 1.4]], 500, 46, &|r| {
        vec![
            r.uniform(1.0, 2000.0),
            r.log_uniform(-5.0, 1.0),
            r.log_uniform(-5.0, 1.0),
            r.uniform(1.0, 1.7),
        ]
    });
    for g in around(1.0, 2)
        .into_iter()
        .chain([1.0 + 1e-12, 5.0 / 3.0, 1.4, 1.3])
    {
        adi.push(vec![300.0, 0.02, 0.01, g]);
    }
    for a in adi {
        row(
            m,
            "state::adiabatic_temperature",
            Class::C,
            &a,
            res(st::adiabatic_temperature(a[0], a[1], a[2], a[3])),
        );
    }
    // compressibility_factor(P, V, n, T)
    for a in tuples(&[vec![101_325.0, 0.02241, 1.0, 273.15]], 400, 47, &|r| {
        vec![
            r.log_uniform(2.0, 8.0),
            r.log_uniform(-5.0, 2.0),
            r.signed_log_uniform(-3.0, 3.0),
            r.uniform(1.0, 3000.0),
        ]
    }) {
        row(
            m,
            "state::compressibility_factor",
            Class::B,
            &a,
            res(st::compressibility_factor(a[0], a[1], a[2], a[3])),
        );
    }
    // reduced_temperature(T, T_c) and reduced_pressure(P, P_c): every gas over a grid, and the edges
    let ts = [
        0.0, 1.0, 100.0, 273.15, 300.0, 500.0, 1000.0, 3000.0, -0.0, -1.0,
    ];
    let ps = [0.0, 1.0, 1e3, 101_325.0, 1e6, 5e6, 2e7, 1e8, -0.0, -1.0];
    let mut specials = vec![f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 1e-300, 1e300];
    specials.push(0.0);
    for g in ALL_GASES {
        for &t in ts.iter().chain(specials.iter()) {
            row(
                m,
                "state::reduced_temperature",
                Class::B,
                &[t, g.critical_t],
                vec![f(st::reduced_temperature(t, g.critical_t))],
            );
        }
        for &p in ps.iter().chain(specials.iter()) {
            row(
                m,
                "state::reduced_pressure",
                Class::B,
                &[p, g.critical_p],
                vec![f(st::reduced_pressure(p, g.critical_p))],
            );
        }
    }
    for &c in &specials {
        row(
            m,
            "state::reduced_temperature",
            Class::B,
            &[300.0, c],
            vec![f(st::reduced_temperature(300.0, c))],
        );
        row(
            m,
            "state::reduced_pressure",
            Class::B,
            &[1e6, c],
            vec![f(st::reduced_pressure(1e6, c))],
        );
    }
}

fn fields(g: &GasData) -> Vec<Tok> {
    vec![
        s(g.name.to_string()),
        f(g.critical_t),
        f(g.critical_p),
        f(g.acentric_factor),
    ]
}

fn gas(name: &str, tc: f64, pc: f64, om: f64) -> GasData {
    GasData {
        name: Cow::Owned(name.to_string()),
        critical_t: tc,
        critical_p: pc,
        acentric_factor: om,
    }
}

fn json_values() -> Vec<GasData> {
    let mut v: Vec<GasData> = ALL_GASES.iter().map(|g| (*g).clone()).collect();
    for name in [
        "",
        "a\"b\\c/d",
        "tab\tnew\nline\r",
        "ctl\u{1}\u{1f}\u{7f}",
        "\u{e9}\u{1f600} ok",
        "\u{0}nul",
    ] {
        v.push(gas(name, 1.0, 2.0, 3.0));
    }
    v.push(gas("specials", f64::NAN, f64::INFINITY, -0.0));
    v.push(gas("negative", f64::NEG_INFINITY, -1.0, 0.0));
    v.push(gas(
        "extremes",
        f64::MAX,
        f64::MIN_POSITIVE,
        f64::from_bits(1),
    ));
    v.push(gas("forms", 1e16, 1e-5, 0.1 + 0.2));
    v
}

fn decoded(doc: &str) -> Vec<Tok> {
    let mut de = serde_json::Deserializer::from_str(doc);
    match GasData::deserialize(&mut de).and_then(|g| de.end().map(|()| g)) {
        Ok(x) => fields(&x),
        Err(_) => vec![Tok::Err("Deserialize".to_string())],
    }
}

fn json_rows(m: &mut Module) {
    let mut docs = Vec::new();
    for x in json_values() {
        let doc = serde_json::to_string(&x).expect("serialize");
        m.row(
            "state::GasData::to_json",
            Class::A,
            &fields(&x),
            &[s(doc.clone())],
        );
        docs.push(doc);
    }
    let nitrogen = docs[0].clone();
    for (_kind, bytes) in variants(&nitrogen) {
        if let Ok(text) = String::from_utf8(bytes)
            && !text.contains('\0')
        {
            docs.push(text);
        }
    }
    for doc in [
        r#"{"name":"N\u0000x","critical_t":1,"critical_p":2.5,"acentric_factor":-0.0}"#,
        r#"{"name":"é😀","critical_t":1,"critical_p":2,"acentric_factor":3}"#,
        r#"["Nitrogen",126.19,3390000,0.037]"#,
        r#"["Nitrogen",126.19,3390000]"#,
        r#"{"acentric_factor":3,"critical_p":2,"critical_t":1,"name":"rev"}"#,
        r#"{"name":"n","critical_t":1,"critical_p":2}"#,
        r#"{"name":null,"critical_t":1,"critical_p":2,"acentric_factor":3}"#,
        r#"{"name":"n","critical_t":[1],"critical_p":2,"acentric_factor":3}"#,
        r#"{"name":"n","critical_t":1,"critical_p":2,"acentric_factor":3,"extra":{"a":[1,{"b":"c"}]}}"#,
        r#"{"name":"n","critical_t":1e400,"critical_p":2,"acentric_factor":3}"#,
        "[]",
        "null",
        "\"Nitrogen\"",
        "",
    ] {
        docs.push(doc.to_string());
    }
    for doc in &docs {
        m.row(
            "state::GasData::from_json",
            Class::A,
            &[s(doc.clone())],
            &decoded(doc),
        );
    }
}

/// Adds the module's function and JSON rows to `m` (the `state` module).
pub fn rows(m: &mut Module) {
    part1(m);
    json_rows(m);
}
