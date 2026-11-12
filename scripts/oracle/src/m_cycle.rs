//! `cycle.tsv` (M4a onward): every public function of Rust 1.4.0's `cycle`.
//!
//! A `CycleResult` is written as its kind (`e:CycleKind.X`), `v:n` and the state points' five words each, `v:m` and
//! the processes (`e:ProcessKind.X`), then its six numbers in declaration order.
//!
//! M4a: the process heats and works and `adiabatic_temp_from_pressure` over seeded sweeps and the edges; the Otto,
//! Diesel and Brayton cycles over r, γ and T grids (r and r_p just above 1, γ just above 1 and up to 5/3), Diesel's
//! cutoff at r ± 1 ulp, Brayton's turbine inlet at the compressor outlet T₂ ± 1 ulp (computed as Rust computes it), and
//! each argument of a base case at the edges. JSON: `StatePoint`, `ProcessKind`, `CycleKind` and `CycleResult`, from
//! the cycles' own results and edges (empty lists, unequal lists, a bad element). Classes: B for + − × ÷, C where `ln`
//! or `powf` enters, A for JSON.

use ushma::cycle::{self as cy, CycleKind, CycleResult, ProcessKind, StatePoint};

use crate::grids::Rng;
use crate::json_cases::variants;
use crate::sweep::{around, res, row, tuples};
use crate::tsv::{Class, Module, Tok, e, err, f, s};

fn sp_toks(p: &StatePoint) -> Vec<Tok> {
    vec![
        f(p.temperature),
        f(p.pressure),
        f(p.volume),
        f(p.entropy),
        f(p.enthalpy),
    ]
}

/// A CycleResult as its tokens.
pub fn cycle_toks(r: &CycleResult) -> Vec<Tok> {
    let mut v = vec![e("CycleKind", r.kind), Tok::V(r.state_points.len())];
    for p in &r.state_points {
        v.extend(sp_toks(p));
    }
    v.push(Tok::V(r.processes.len()));
    for k in &r.processes {
        v.push(e("ProcessKind", k));
    }
    v.extend([
        f(r.heat_in),
        f(r.heat_out),
        f(r.net_work),
        f(r.efficiency),
        f(r.back_work_ratio),
        f(r.gamma),
    ]);
    v
}

fn cyc_res(r: ushma::error::Result<CycleResult>) -> Vec<Tok> {
    match r {
        Ok(c) => cycle_toks(&c),
        Err(z) => vec![err(&z)],
    }
}

fn helpers(m: &mut Module) {
    let draw = |r: &mut Rng| {
        vec![
            r.log_uniform(-2.0, 2.0),
            r.uniform(10.0, 60.0),
            r.uniform(1.0, 3000.0),
            r.uniform(1.0, 3000.0),
        ]
    };
    for (name, seed) in [
        ("cycle::isochoric_heat", 90),
        ("cycle::isobaric_heat", 91),
        ("cycle::adiabatic_work", 92),
    ] {
        for a in tuples(&[vec![1.0, 20.8, 300.0, 600.0]], 200, seed, &draw) {
            let out = match name {
                "cycle::isochoric_heat" => cy::isochoric_heat(a[0], a[1], a[2], a[3]),
                "cycle::isobaric_heat" => cy::isobaric_heat(a[0], a[1], a[2], a[3]),
                _ => cy::adiabatic_work(a[0], a[1], a[2], a[3]),
            };
            row(m, name, Class::B, &a, vec![f(out)]);
        }
    }
    let mut at = tuples(&[vec![300.0, 101_325.0, 1_013_250.0, 1.4]], 400, 93, &|r| {
        vec![
            r.uniform(1.0, 2000.0),
            r.log_uniform(3.0, 7.0),
            r.log_uniform(3.0, 8.0),
            r.uniform(1.0, 1.7),
        ]
    });
    for g in around(1.0, 2) {
        at.push(vec![300.0, 101_325.0, 1_013_250.0, g]);
    }
    for a in at {
        row(
            m,
            "cycle::adiabatic_temp_from_pressure",
            Class::C,
            &a,
            res(cy::adiabatic_temp_from_pressure(a[0], a[1], a[2], a[3])),
        );
    }
}

fn gas_cycles(m: &mut Module) {
    let rs = [1.5, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0, 30.0];
    let gs = [1.1, 1.3, 1.4, 5.0 / 3.0];
    let ts = [250.0, 300.0, 800.0];
    // otto_cycle(T₁, P₁, r, q_in, γ, n)
    let mut otto: Vec<Vec<f64>> = Vec::new();
    for &r in rs.iter().chain(around(1.0, 1).iter().skip(2)) {
        for &g in gs.iter().chain(around(1.0, 1).iter().skip(2)) {
            for &t in &ts {
                otto.push(vec![t, 101_325.0, r, 50_000.0, g, 1.0]);
            }
        }
    }
    otto.extend(tuples(
        &[vec![300.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0]],
        300,
        94,
        &|r| {
            vec![
                r.uniform(200.0, 1000.0),
                r.log_uniform(4.0, 7.0),
                r.uniform(1.0, 25.0),
                r.log_uniform(3.0, 6.0),
                r.uniform(1.05, 1.7),
                r.log_uniform(-2.0, 2.0),
            ]
        },
    ));
    for a in otto {
        row(
            m,
            "cycle::otto_cycle",
            Class::C,
            &a,
            cyc_res(cy::otto_cycle(a[0], a[1], a[2], a[3], a[4], a[5])),
        );
    }
    // diesel_cycle(T₁, P₁, r, r_c, γ, n), with r_c at r ± 1 ulp
    let mut diesel: Vec<Vec<f64>> = Vec::new();
    for &r in &[4.0, 12.0, 16.0, 20.0, 30.0] {
        for &rc in &[1.5, 2.0, 3.0] {
            for &g in &gs {
                diesel.push(vec![300.0, 101_325.0, r, rc, g, 1.0]);
            }
        }
        for rc in around(r, 1) {
            diesel.push(vec![300.0, 101_325.0, r, rc, 1.4, 1.0]);
        }
    }
    diesel.extend(
        around(1.0, 1)
            .into_iter()
            .map(|rc| vec![300.0, 101_325.0, 20.0, rc, 1.4, 1.0]),
    );
    diesel.extend(tuples(
        &[vec![300.0, 101_325.0, 20.0, 2.0, 1.4, 1.0]],
        300,
        95,
        &|r| {
            let cr = r.uniform(5.0, 25.0);
            vec![
                r.uniform(200.0, 1000.0),
                r.log_uniform(4.0, 7.0),
                cr,
                r.uniform(1.01, cr),
                r.uniform(1.05, 1.7),
                r.log_uniform(-2.0, 2.0),
            ]
        },
    ));
    for a in diesel {
        row(
            m,
            "cycle::diesel_cycle",
            Class::C,
            &a,
            cyc_res(cy::diesel_cycle(a[0], a[1], a[2], a[3], a[4], a[5])),
        );
    }
    // brayton_cycle(T₁, P₁, r_p, T₃, γ, n), with T₃ at T₂ ± 1 ulp
    let mut bray: Vec<Vec<f64>> = Vec::new();
    for &rp in &[1.5, 4.0, 8.0, 10.0, 15.0, 25.0] {
        for &g in &gs {
            for &t3 in &[1000.0, 1400.0, 1800.0] {
                bray.push(vec![300.0, 101_325.0, rp, t3, g, 1.0]);
            }
            let t2 = 300.0 * rp.powf((g - 1.0) / g);
            for t3 in around(t2, 1) {
                bray.push(vec![300.0, 101_325.0, rp, t3, g, 1.0]);
            }
        }
    }
    bray.extend(
        around(1.0, 1)
            .into_iter()
            .map(|rp| vec![300.0, 101_325.0, rp, 1400.0, 1.4, 1.0]),
    );
    bray.extend(tuples(
        &[vec![300.0, 101_325.0, 10.0, 1400.0, 1.4, 1.0]],
        300,
        96,
        &|r| {
            vec![
                r.uniform(200.0, 600.0),
                r.log_uniform(4.0, 7.0),
                r.uniform(1.0, 40.0),
                r.uniform(500.0, 2200.0),
                r.uniform(1.05, 1.7),
                r.log_uniform(-2.0, 2.0),
            ]
        },
    ));
    for a in bray {
        row(
            m,
            "cycle::brayton_cycle",
            Class::C,
            &a,
            cyc_res(cy::brayton_cycle(a[0], a[1], a[2], a[3], a[4], a[5])),
        );
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

fn json_rows(m: &mut Module) {
    let otto = cy::otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0).expect("otto");
    let bray = cy::brayton_cycle(300.0, 101_325.0, 10.0, 1400.0, 1.4, 1.0).expect("brayton");
    // StatePoint
    let mut sps: Vec<StatePoint> = otto.state_points.clone();
    sps.push(StatePoint {
        temperature: f64::NAN,
        pressure: f64::INFINITY,
        volume: -0.0,
        entropy: f64::MAX,
        enthalpy: f64::from_bits(1),
    });
    let mut docs = Vec::new();
    for p in &sps {
        let doc = serde_json::to_string(p).expect("serialize");
        m.row(
            "cycle::StatePoint::to_json",
            Class::A,
            &sp_toks(p),
            &[s(doc.clone())],
        );
        docs.push(doc);
    }
    let base = docs[0].clone();
    for (_kind, bytes) in variants(&base) {
        if let Ok(text) = String::from_utf8(bytes)
            && !text.contains('\0')
        {
            docs.push(text);
        }
    }
    docs.extend(["[1,2,3,4,5]", "[1,2,3,4]", "{}", "null"].map(String::from));
    for doc in &docs {
        m.row(
            "cycle::StatePoint::from_json",
            Class::A,
            &[s(doc.clone())],
            &decode::<StatePoint>(doc, sp_toks),
        );
    }
    // ProcessKind and CycleKind
    let procs = [
        ProcessKind::Isentropic,
        ProcessKind::Isochoric,
        ProcessKind::Isobaric,
        ProcessKind::Isothermal,
        ProcessKind::Isenthalpic,
    ];
    for p in procs {
        m.row(
            "cycle::ProcessKind::to_json",
            Class::A,
            &[e("ProcessKind", p)],
            &[s(serde_json::to_string(&p).expect("serialize"))],
        );
    }
    let kinds = [
        CycleKind::Otto,
        CycleKind::Diesel,
        CycleKind::Brayton,
        CycleKind::Rankine,
        CycleKind::Refrigeration,
    ];
    for k in kinds {
        m.row(
            "cycle::CycleKind::to_json",
            Class::A,
            &[e("CycleKind", k)],
            &[s(serde_json::to_string(&k).expect("serialize"))],
        );
    }
    for doc in [
        "\"Isentropic\"",
        "\"Isenthalpic\"",
        "\"isobaric\"",
        "\"Adiabatic\"",
        "{\"Isochoric\":null}",
        "{\"Isochoric\":0}",
        "[\"Isobaric\"]",
        "null",
        "",
    ] {
        m.row(
            "cycle::ProcessKind::from_json",
            Class::A,
            &[s(doc)],
            &decode::<ProcessKind>(doc, |k| vec![e("ProcessKind", k)]),
        );
    }
    for doc in [
        "\"Otto\"",
        "\"Refrigeration\"",
        "\"otto\"",
        "\"Stirling\"",
        "{\"Diesel\":null}",
        "{\"Diesel\":null,\"Otto\":null}",
        "1",
        "",
    ] {
        m.row(
            "cycle::CycleKind::from_json",
            Class::A,
            &[s(doc)],
            &decode::<CycleKind>(doc, |k| vec![e("CycleKind", k)]),
        );
    }
    // CycleResult
    let mut results = vec![otto.clone(), bray.clone()];
    let mut empty = otto.clone();
    empty.state_points.clear();
    empty.processes.clear();
    results.push(empty);
    let mut uneven = bray.clone();
    uneven.processes.truncate(2);
    uneven.heat_in = f64::NAN;
    uneven.gamma = f64::INFINITY;
    results.push(uneven);
    let mut rdocs = Vec::new();
    for r in &results {
        let doc = serde_json::to_string(r).expect("serialize");
        m.row(
            "cycle::CycleResult::to_json",
            Class::A,
            &cycle_toks(r),
            &[s(doc.clone())],
        );
        rdocs.push(doc);
    }
    let base = rdocs[0].clone();
    for (_kind, bytes) in variants(&base) {
        if let Ok(text) = String::from_utf8(bytes)
            && !text.contains('\0')
        {
            rdocs.push(text);
        }
    }
    rdocs.extend(
        [
            r#"{"kind":"Otto","state_points":[],"processes":[],"heat_in":1,"heat_out":2,"net_work":3,"efficiency":4,"back_work_ratio":5,"gamma":6}"#,
            r#"["Brayton",[[1,2,3,4,5]],["Isobaric"],1,2,3,4,5,6]"#,
            r#"{"kind":"Otto","state_points":[{"temperature":1,"pressure":2,"volume":3,"entropy":4}],"processes":[],"heat_in":1,"heat_out":2,"net_work":3,"efficiency":4,"back_work_ratio":5,"gamma":6}"#,
            r#"{"kind":"Otto","state_points":[],"processes":["Adiabatic"],"heat_in":1,"heat_out":2,"net_work":3,"efficiency":4,"back_work_ratio":5,"gamma":6}"#,
            r#"{"kind":"Otto","state_points":{},"processes":[],"heat_in":1,"heat_out":2,"net_work":3,"efficiency":4,"back_work_ratio":5,"gamma":6}"#,
            r#"{"kind":"Jet","state_points":[],"processes":[],"heat_in":1,"heat_out":2,"net_work":3,"efficiency":4,"back_work_ratio":5,"gamma":6}"#,
            "[]",
            "null",
        ]
        .map(String::from),
    );
    for doc in &rdocs {
        m.row(
            "cycle::CycleResult::from_json",
            Class::A,
            &[s(doc.clone())],
            &decode::<CycleResult>(doc, cycle_toks),
        );
    }
}

pub fn emit() -> Module {
    let mut m = Module::new("cycle");
    helpers(&mut m);
    gas_cycles(&mut m);
    json_rows(&mut m);
    m
}
