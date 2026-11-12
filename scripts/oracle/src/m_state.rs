//! `state.tsv`'s function and JSON rows (M2a onward), added to the module the table presets seed (`t_tables.rs`).
//!
//! For each function: a seeded sweep over its arguments' physical ranges, each argument of each base case at every
//! edge value (`sweep::edges`), then the function's own switches (V − nb across 0 for van der Waals, γ just above 1,
//! V₁ = V₂, nRT at the 1e-30 divisor threshold). The reduced coordinates run every gas over a grid of T and P, and the
//! cubic equations every gas × T across T_c × V_m across the co-volume b; Pitzer's Z runs T_r × P_r out to where
//! Z ≤ 0; Kay's rules run mixtures with zero, NaN and negative fractions, unequal lengths and sums across the 1e-10
//! tolerance. Classes: B for + − × ÷ and √, C where `ln` or `powf` enters (`isothermal_work`,
//! `adiabatic_temperature`, `redlich_kwong_*`, `pitzer_second_virial`, `compressibility_pitzer`).
//!
//! `GasData`'s JSON (class A), as `material.tsv`'s: serde_json's document for the presets and for names and values
//! that exercise its escapes and number forms, and serde_json's verdict on each of those documents, on the variants of
//! Nitrogen's (`json_cases::variants`) and on hand-written edges. Each row's input is the four fields in declaration
//! order (`s:name`, then critical_t, critical_p, acentric_factor), or for `from_json` the document.

use std::borrow::Cow;

use serde::Deserialize as _;
use ushma::state::{self as st, ALL_GASES, GasData};

use crate::grids::Rng;
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

/// T across the critical point, for a gas of critical temperature `tc`.
fn temps(tc: f64) -> Vec<f64> {
    [0.3, 0.5, 0.8, 0.95, 0.99, 1.0, 1.01, 1.05, 1.5, 3.0, 10.0]
        .iter()
        .map(|k| k * tc)
        .collect()
}

/// Molar volumes from just below the co-volume `b` (V_m − b across 0) to the ideal-gas range.
fn volumes(b: f64) -> Vec<f64> {
    let mut v = around(b, 2);
    v.extend(
        [1.01, 1.5, 2.0, 3.0, 4.0, 10.0, 100.0]
            .iter()
            .map(|k| k * b),
    );
    v.extend([1e-3, 0.0224, 1.0]);
    v
}

/// Kay's-rule inputs: each list as `v:n` then its elements, fractions first.
fn lists(ys: &[f64], xs: &[f64]) -> Vec<Tok> {
    let mut v = vec![Tok::V(ys.len())];
    v.extend(ys.iter().map(|&y| f(y)));
    v.push(Tok::V(xs.len()));
    v.extend(xs.iter().map(|&x| f(x)));
    v
}

fn part2(m: &mut Module) {
    let edge_gases = [vec![304.13, 7_375_000.0], vec![126.19, 3_390_000.0]];
    // redlich_kwong_params(T_c, P_c) -> (a, b)
    let mut rkp: Vec<Vec<f64>> = ALL_GASES
        .iter()
        .map(|g| vec![g.critical_t, g.critical_p])
        .collect();
    rkp.extend(tuples(&edge_gases, 300, 50, &|r| {
        vec![r.log_uniform(0.0, 3.5), r.log_uniform(4.0, 8.0)]
    }));
    for a in rkp {
        let (ra, rb) = st::redlich_kwong_params(a[0], a[1]);
        row(
            m,
            "state::redlich_kwong_params",
            Class::C,
            &a,
            vec![f(ra), f(rb)],
        );
    }
    // peng_robinson_params(T_c, P_c, ω) -> (a, b, κ)
    let mut prp: Vec<Vec<f64>> = ALL_GASES
        .iter()
        .map(|g| vec![g.critical_t, g.critical_p, g.acentric_factor])
        .collect();
    prp.extend(tuples(
        &[
            vec![304.13, 7_375_000.0, 0.224],
            vec![126.19, 3_390_000.0, 0.037],
        ],
        300,
        51,
        &|r| {
            vec![
                r.log_uniform(0.0, 3.5),
                r.log_uniform(4.0, 8.0),
                r.uniform(-0.5, 1.5),
            ]
        },
    ));
    for a in prp {
        let (pa, pb, pk) = st::peng_robinson_params(a[0], a[1], a[2]);
        row(
            m,
            "state::peng_robinson_params",
            Class::B,
            &a,
            vec![f(pa), f(pb), f(pk)],
        );
    }
    // redlich_kwong_pressure(T, V_m, T_c, P_c) and peng_robinson_pressure(T, V_m, T_c, P_c, ω): every gas × T × V_m
    let mut rk: Vec<Vec<f64>> = Vec::new();
    let mut pr: Vec<Vec<f64>> = Vec::new();
    for g in ALL_GASES {
        let (_, brk) = st::redlich_kwong_params(g.critical_t, g.critical_p);
        let (_, bpr, _) = st::peng_robinson_params(g.critical_t, g.critical_p, g.acentric_factor);
        for t in temps(g.critical_t) {
            for vm in volumes(brk) {
                rk.push(vec![t, vm, g.critical_t, g.critical_p]);
            }
            for vm in volumes(bpr) {
                pr.push(vec![t, vm, g.critical_t, g.critical_p, g.acentric_factor]);
            }
        }
    }
    rk.extend(tuples(
        &[vec![273.15, 0.02241, 304.13, 7_375_000.0]],
        300,
        52,
        &|r| {
            vec![
                r.uniform(50.0, 2000.0),
                r.log_uniform(-5.0, 0.0),
                r.uniform(100.0, 700.0),
                r.log_uniform(6.0, 7.5),
            ]
        },
    ));
    pr.extend(tuples(
        &[vec![273.15, 0.02241, 304.13, 7_375_000.0, 0.224]],
        300,
        53,
        &|r| {
            vec![
                r.uniform(50.0, 2000.0),
                r.log_uniform(-5.0, 0.0),
                r.uniform(100.0, 700.0),
                r.log_uniform(6.0, 7.5),
                r.uniform(-0.3, 1.0),
            ]
        },
    ));
    for a in rk {
        row(
            m,
            "state::redlich_kwong_pressure",
            Class::C,
            &a,
            res(st::redlich_kwong_pressure(a[0], a[1], a[2], a[3])),
        );
    }
    for a in pr {
        row(
            m,
            "state::peng_robinson_pressure",
            Class::B,
            &a,
            res(st::peng_robinson_pressure(a[0], a[1], a[2], a[3], a[4])),
        );
    }
    // pitzer_second_virial(T, T_c, P_c, ω): every gas over T, seeded, and the edges
    let mut pz: Vec<Vec<f64>> = Vec::new();
    for g in ALL_GASES {
        for t in temps(g.critical_t) {
            pz.push(vec![t, g.critical_t, g.critical_p, g.acentric_factor]);
        }
    }
    pz.extend(tuples(
        &[vec![300.0, 304.13, 7_375_000.0, 0.224]],
        300,
        54,
        &|r| {
            vec![
                r.uniform(50.0, 2000.0),
                r.uniform(100.0, 700.0),
                r.log_uniform(6.0, 7.5),
                r.uniform(-0.3, 1.0),
            ]
        },
    ));
    for a in pz {
        row(
            m,
            "state::pitzer_second_virial",
            Class::C,
            &a,
            vec![f(st::pitzer_second_virial(a[0], a[1], a[2], a[3]))],
        );
    }
    // virial_pressure_2nd(T, V_m, B) and virial_pressure_3rd(T, V_m, B, C)
    for a in tuples(&[vec![300.0, 0.02241, -1e-4]], 300, 55, &|r| {
        vec![
            r.uniform(1.0, 2000.0),
            r.log_uniform(-5.0, 1.0),
            r.signed_log_uniform(-7.0, -2.0),
        ]
    }) {
        row(
            m,
            "state::virial_pressure_2nd",
            Class::B,
            &a,
            res(st::virial_pressure_2nd(a[0], a[1], a[2])),
        );
    }
    for a in tuples(&[vec![300.0, 0.02241, -1e-4, 1e-8]], 300, 56, &|r| {
        vec![
            r.uniform(1.0, 2000.0),
            r.log_uniform(-5.0, 1.0),
            r.signed_log_uniform(-7.0, -2.0),
            r.signed_log_uniform(-12.0, -6.0),
        ]
    }) {
        row(
            m,
            "state::virial_pressure_3rd",
            Class::B,
            &a,
            res(st::virial_pressure_3rd(a[0], a[1], a[2], a[3])),
        );
    }
    // compressibility_pitzer(P_r, T_r, ω): T_r × P_r out to where Z ≤ 0, for four ω
    let mut cp: Vec<Vec<f64>> = Vec::new();
    for om in [0.0, 0.037, 0.224, 0.344, -0.2] {
        for tr in [0.3, 0.5, 0.7, 0.9, 1.0, 1.1, 1.5, 2.0, 3.0, 5.0] {
            for prr in [0.0, -0.0, 0.01, 0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0] {
                cp.push(vec![prr, tr, om]);
            }
        }
    }
    cp.extend(tuples(
        &[vec![0.5, 1.5, 0.1], vec![5.0, 1.0, 0.0]],
        400,
        57,
        &|r| {
            vec![
                r.uniform(0.0, 10.0),
                r.uniform(0.1, 5.0),
                r.uniform(-0.3, 1.0),
            ]
        },
    ));
    for a in cp {
        row(
            m,
            "state::compressibility_pitzer",
            Class::C,
            &a,
            res(st::compressibility_pitzer(a[0], a[1], a[2])),
        );
    }
    // kays_rule_tc / pc / omega(ys, xs): pure, binary, ternary, sums across the 1e-10 tolerance, and the refusals
    let nan = f64::NAN;
    let mut mix: Vec<(Vec<f64>, Vec<f64>)> = vec![
        (vec![1.0], vec![304.13]),
        (vec![0.5, 0.5], vec![126.19, 154.58]),
        (vec![0.78, 0.21, 0.01], vec![126.19, 154.58, 150.86]),
        (vec![0.0, 1.0], vec![126.19, 154.58]),
        (vec![-0.0, 1.0], vec![126.19, 154.58]),
        (vec![0.3, 0.3], vec![126.19, 154.58]),
        (vec![0.5, 0.5], vec![126.19]),
        (vec![], vec![]),
        (vec![nan], vec![300.0]),
        (vec![1.5, -0.5], vec![300.0, 200.0]),
        (vec![0.5, 0.5], vec![4e6, nan]),
        (vec![0.5, 0.5], vec![f64::INFINITY, 1.0]),
        (vec![0.5, 0.5], vec![1e308, 1e308]),
        (vec![1.0], vec![-0.0]),
    ];
    for y in around(1.0 - 1e-10, 3).into_iter().chain(around(1.0, 2)) {
        mix.push((vec![y], vec![300.0]));
    }
    for y in around(0.5 + 1e-10, 3) {
        mix.push((vec![0.5, y], vec![300.0, 200.0]));
    }
    let mut r = Rng::new(58);
    for _ in 0..200 {
        let n = r.int(1, 6) as usize;
        let mut ys: Vec<f64> = (0..n).map(|_| r.unit()).collect();
        let s: f64 = ys.iter().sum();
        if r.int(0, 3) != 0 {
            for y in &mut ys {
                *y /= s;
            }
        }
        let xs: Vec<f64> = (0..n).map(|_| r.log_uniform(-1.0, 7.5)).collect();
        mix.push((ys, xs));
    }
    for (ys, xs) in &mix {
        let ins = lists(ys, xs);
        m.row(
            "state::kays_rule_tc",
            Class::B,
            &ins,
            &res(st::kays_rule_tc(ys, xs)),
        );
        m.row(
            "state::kays_rule_pc",
            Class::B,
            &ins,
            &res(st::kays_rule_pc(ys, xs)),
        );
        m.row(
            "state::kays_rule_omega",
            Class::B,
            &ins,
            &res(st::kays_rule_omega(ys, xs)),
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
    part2(m);
    json_rows(m);
}
