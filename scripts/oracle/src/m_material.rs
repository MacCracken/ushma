//! `material.tsv`'s function and JSON rows (M1b), added to the module the table presets seed (`t_tables.rs`).
//!
//! Every row's input is a `ThermalMaterial` written as its six fields in declaration order (`s:name`, then
//! conductivity, specific_heat, density, melting_point, boiling_point), or for `from_json` the document:
//! - `material::diffusivity`, `material::volumetric_heat_capacity` (class B): the 10 presets, a grid over every
//!   property's physical values and edges (0, −0, negative, NaN, ±∞, 1e-300, 1e300), and 300 seeded materials;
//! - `material::to_json` (class A): serde_json's document for the presets and for names and values that exercise its
//!   escapes and number forms (control bytes, quotes, backslashes, non-ASCII, NaN and ±∞ as `null`, −0.0, the
//!   extremes);
//! - `material::from_json` (class A): serde_json's verdict on every `to_json` document and on the 60-odd variants of
//!   Copper's (`json_cases::variants`: malformed, duplicate, unknown and missing fields, the sequence form, number and
//!   string forms): the six fields read, or `E:Deserialize`. Documents that aren't UTF-8 are left to the unit tests,
//!   since a golden `s:` token is text.

use serde::Deserialize as _;
use ushma::material::{ALL_MATERIALS, ThermalMaterial};

use crate::grids::Rng;
use crate::json_cases::variants;
use crate::tsv::{Class, Module, Tok, f, s};

fn fields(m: &ThermalMaterial) -> Vec<Tok> {
    vec![
        s(m.name.to_string()),
        f(m.conductivity),
        f(m.specific_heat),
        f(m.density),
        f(m.melting_point),
        f(m.boiling_point),
    ]
}

fn mat(name: &str, k: f64, cp: f64, rho: f64, mp: f64, bp: f64) -> ThermalMaterial {
    ThermalMaterial {
        name: name.to_string().into(),
        conductivity: k,
        specific_heat: cp,
        density: rho,
        melting_point: mp,
        boiling_point: bp,
    }
}

fn presets() -> Vec<ThermalMaterial> {
    ALL_MATERIALS.iter().map(|m| (*m).clone()).collect()
}

fn method_rows(m: &mut Module) {
    let edges = [0.0, -0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 1e-300, 1e300];
    let mut ms = presets();
    for k in [401.0, 0.0263, 2200.0].iter().chain(edges.iter()) {
        for cp in [385.0, 4186.0].iter().chain(edges.iter()) {
            for rho in [8960.0, 1.1614].iter().chain(edges.iter()) {
                ms.push(mat("x", *k, *cp, *rho, 0.0, 0.0));
            }
        }
    }
    let mut r = Rng::new(0x6d61_7465_7269_616c);
    for _ in 0..300 {
        let k = r.log_uniform(-3.0, 4.0);
        let cp = r.log_uniform(1.0, 4.5);
        let rho = r.log_uniform(-1.0, 4.5);
        ms.push(mat("seeded", k, cp, rho, r.uniform(0.0, 4000.0), r.uniform(0.0, 6000.0)));
    }
    for x in &ms {
        let ins = fields(x);
        m.row("material::diffusivity", Class::B, &ins, &[f(x.diffusivity())]);
        m.row("material::volumetric_heat_capacity", Class::B, &ins, &[f(x.volumetric_heat_capacity())]);
    }
}

fn json_values() -> Vec<ThermalMaterial> {
    let mut v = presets();
    for name in ["", "a\"b\\c/d", "tab\tnew\nline\r", "ctl\u{1}\u{1f}\u{7f}", "\u{e9}\u{1f600} ok", "\u{0}nul"] {
        v.push(mat(name, 1.0, 2.0, 3.0, 4.0, 5.0));
    }
    v.push(mat("specials", f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.0, 0.0));
    v.push(mat("extremes", f64::MAX, f64::MIN_POSITIVE, f64::from_bits(1), 1e16, 1e-5));
    v.push(mat("forms", 1e15, 123_456_789_012_345_680.0, 0.000_012_5, 9.999_999_999_999_999e-6, 0.1 + 0.2));
    v
}

fn decoded(doc: &str) -> Vec<Tok> {
    let mut de = serde_json::Deserializer::from_str(doc);
    match ThermalMaterial::deserialize(&mut de).and_then(|m| de.end().map(|()| m)) {
        Ok(x) => fields(&x),
        Err(_) => vec![Tok::Err("Deserialize".to_string())],
    }
}

fn json_rows(m: &mut Module) {
    let mut docs = Vec::new();
    for x in json_values() {
        let doc = serde_json::to_string(&x).expect("serialize");
        m.row("material::to_json", Class::A, &fields(&x), &[s(doc.clone())]);
        docs.push(doc);
    }
    let copper = docs[0].clone();
    for (_kind, bytes) in variants(&copper) {
        if let Ok(text) = String::from_utf8(bytes)
            && !text.contains('\0')
        {
            docs.push(text);
        }
    }
    for doc in [
        r#"{"name":"Cu\u0000x","conductivity":1,"specific_heat":2.5,"density":3,"melting_point":-0.0,"boiling_point":1e3}"#,
        r#"{"name":"é😀","conductivity":1,"specific_heat":2,"density":3,"melting_point":4,"boiling_point":5}"#,
        r#"["Copper",401,385,8960,1358,2835]"#,
        r#"{"boiling_point":5,"melting_point":4,"density":3,"specific_heat":2,"conductivity":1,"name":"rev"}"#,
        r#"{"name":"n","conductivity":1,"specific_heat":2,"density":3,"melting_point":4}"#,
        r#"{"name":null,"conductivity":1,"specific_heat":2,"density":3,"melting_point":4,"boiling_point":5}"#,
        r#"{"name":"n","conductivity":[1],"specific_heat":2,"density":3,"melting_point":4,"boiling_point":5}"#,
        r#"{"name":"n","conductivity":1,"specific_heat":2,"density":3,"melting_point":4,"boiling_point":5,"extra":{"a":[1,{"b":"c"}]}}"#,
        "[]",
        "null",
        "\"Copper\"",
        "",
    ] {
        docs.push(doc.to_string());
    }
    for doc in &docs {
        m.row("material::from_json", Class::A, &[s(doc.clone())], &decoded(doc));
    }
}

/// Adds the module's function and JSON rows to `m` (the `material` module).
pub fn rows(m: &mut Module) {
    method_rows(m);
    json_rows(m);
}
