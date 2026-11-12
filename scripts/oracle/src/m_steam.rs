//! `steam.tsv` (M3b onward): every public function of Rust 1.4.0's `steam`.
//!
//! The private tables are reached through the public lookups only. Region 4 runs both ways over its range: a dense
//! sweep of 273.15 … 647.096 K, every table row's temperature ± 1 ulp and the limits ± 2 ulps, and the inverse at each
//! of those pressures (and its neighbours) back to T. The saturated lookups run every table row by temperature and by
//! pressure, each ± 1 ulp, the midpoints between rows, seeded points and the edges. Quality runs saturated states (rows,
//! midpoints and the critical row, where v_g = v_f and h_fg = 0) at the liquid and vapor values, midway, past them by
//! the 1e-9 round-off tolerance ± 2 ulps, at x = 1 + 2e-16 and at the edges; wet steam at x across [0, 1] and just
//! outside. Superheated steam (M3c): the lookup at every tier pressure ± 1 ulp and between tiers, at every table
//! temperature ± 1 ulp, the midpoints and around each pressure's saturation temperature (the tier's anchor);
//! `superheated_by_entropy` at each of those pressures over the entropies of the breakpoints ± 1 ulp, the midpoints and
//! past the range. JSON: `SaturatedEntry`, `WetSteamProperties` and `SuperheatedEntry`, as `material.tsv`'s. Classes:
//! B (Region 4 is + − × ÷ √), C for superheated steam (ln and exp between tiers), A for JSON.

use serde::Deserialize as _;
use ushma::steam::{self as st, SaturatedEntry, SuperheatedEntry, WetSteamProperties};

use crate::grids::Rng;
use crate::json_cases::variants;
use crate::sweep::{around, edges, res};
use crate::tsv::{Class, Module, Tok, err, f, s};

/// Every saturated table row's (temperature, pressure), parsed from the frozen source as `steam_data.cyr`'s columns
/// are (`t_steam::saturated_rows`; the table is private).
fn table_points() -> Vec<(f64, f64)> {
    let src = include_str!("../../../rust-old/src/steam.rs");
    crate::t_steam::saturated_rows(src)
        .expect("SATURATED_TABLE parses")
        .iter()
        .map(|r| (r[0], r[1]))
        .collect()
}

fn sat_fields(e: &SaturatedEntry) -> Vec<Tok> {
    vec![
        f(e.temperature),
        f(e.pressure),
        f(e.v_f),
        f(e.v_g),
        f(e.h_f),
        f(e.h_fg),
        f(e.h_g),
        f(e.s_f),
        f(e.s_g),
    ]
}

fn wet_fields(w: &WetSteamProperties) -> Vec<Tok> {
    vec![
        f(w.quality),
        f(w.specific_volume),
        f(w.specific_enthalpy),
        f(w.specific_entropy),
    ]
}

fn sat_res(r: ushma::error::Result<SaturatedEntry>) -> Vec<Tok> {
    match r {
        Ok(e) => sat_fields(&e),
        Err(z) => vec![err(&z)],
    }
}

fn region4(m: &mut Module) {
    let pts = table_points();
    let mut ts: Vec<f64> = Vec::new();
    for &(t, _) in &pts {
        ts.extend(around(t, 1));
    }
    ts.extend(around(273.15, 2));
    ts.extend(around(647.096, 2));
    let mut r = Rng::new(80);
    for _ in 0..1500 {
        ts.push(r.uniform(273.15, 647.096));
    }
    ts.extend(edges());
    ts.extend([273.0, 700.0, 300.0]);
    let mut ps: Vec<f64> = Vec::new();
    for &t in &ts {
        m.row(
            "steam::saturation_pressure",
            Class::B,
            &[f(t)],
            &res(st::saturation_pressure(t)),
        );
        if let Ok(p) = st::saturation_pressure(t) {
            ps.push(p);
        }
    }
    let mut back: Vec<f64> = Vec::new();
    for &p in &ps {
        back.extend(around(p, 1));
    }
    if let Ok(pmin) = st::saturation_pressure(273.15) {
        back.extend(around(pmin, 2));
    }
    back.extend(around(22_064_000.0, 2));
    for &(_, p) in &pts {
        back.push(p);
    }
    back.extend(edges());
    back.extend([611.0, 3e7, 101_325.0]);
    for &p in &back {
        m.row(
            "steam::saturation_temperature",
            Class::B,
            &[f(p)],
            &res(st::saturation_temperature(p)),
        );
    }
}

fn lookups(m: &mut Module) -> Vec<SaturatedEntry> {
    let pts = table_points();
    let mut ts: Vec<f64> = Vec::new();
    let mut ps: Vec<f64> = Vec::new();
    for (k, &(t, p)) in pts.iter().enumerate() {
        ts.extend(around(t, 1));
        ps.extend(around(p, 1));
        if let Some(&(t2, p2)) = pts.get(k + 1) {
            ts.push(0.5 * (t + t2));
            ps.push((p * p2).sqrt());
        }
    }
    let mut r = Rng::new(81);
    for _ in 0..400 {
        ts.push(r.uniform(273.16, 647.096));
        ps.push(r.log_uniform(2.8, 7.34));
    }
    ts.extend(edges());
    ps.extend(edges());
    ts.extend([273.15, 647.1, 373.15]);
    ps.extend([611.0, 22_100_000.0, 101_325.0]);
    let mut states = Vec::new();
    for &t in &ts {
        let r = st::saturated_by_temperature(t);
        if let Ok(e) = r {
            states.push(e);
        }
        m.row(
            "steam::saturated_by_temperature",
            Class::B,
            &[f(t)],
            &sat_res(r),
        );
    }
    for &p in &ps {
        m.row(
            "steam::saturated_by_pressure",
            Class::B,
            &[f(p)],
            &sat_res(st::saturated_by_pressure(p)),
        );
    }
    states
}

fn quality(m: &mut Module, all: &[SaturatedEntry]) {
    // The states: one in 8 of the lookups' (rows and midpoints), the critical row, and an edge state.
    let mut states: Vec<SaturatedEntry> = all.iter().step_by(8).copied().collect();
    if let Ok(c) = st::saturated_by_temperature(647.096) {
        states.push(c);
    }
    if let Ok(mut z) = st::saturated_by_temperature(373.15) {
        z.v_g = f64::NAN;
        z.h_fg = f64::INFINITY;
        z.s_g = z.s_f;
        states.push(z);
    }
    let tol = 1e-9_f64;
    for e in &states {
        let probe = |lo: f64, hi: f64| -> Vec<f64> {
            let d = hi - lo;
            let mut v = vec![lo, hi, 0.5 * (lo + hi), lo + 0.25 * d];
            for x in [1.0 + tol, 1.0 + 2e-16, -tol] {
                let target = lo + x * d;
                if target.is_finite() && target != 0.0 {
                    v.extend(
                        around(target.abs(), 2)
                            .into_iter()
                            .map(|a| a.copysign(target)),
                    );
                }
            }
            v.extend([f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 0.0]);
            v
        };
        for v in probe(e.v_f, e.v_g) {
            let mut ins = vec![f(v)];
            ins.extend(sat_fields(e));
            m.row(
                "steam::quality_from_volume",
                Class::B,
                &ins,
                &res(st::quality_from_volume(v, e)),
            );
        }
        for h in probe(e.h_f, e.h_f + e.h_fg) {
            let mut ins = vec![f(h)];
            ins.extend(sat_fields(e));
            m.row(
                "steam::quality_from_enthalpy",
                Class::B,
                &ins,
                &res(st::quality_from_enthalpy(h, e)),
            );
        }
        for s in probe(e.s_f, e.s_g) {
            let mut ins = vec![f(s)];
            ins.extend(sat_fields(e));
            m.row(
                "steam::quality_from_entropy",
                Class::B,
                &ins,
                &res(st::quality_from_entropy(s, e)),
            );
        }
        let mut xs = vec![0.0, -0.0, 1.0, 0.5, 0.25, 0.9];
        xs.extend(around(1.0, 1));
        xs.extend([
            -f64::from_bits(1),
            f64::from_bits(1),
            1.0 + 2e-16,
            f64::NAN,
            f64::INFINITY,
            -1.0,
        ]);
        for x in xs {
            let mut ins = vec![f(x)];
            ins.extend(sat_fields(e));
            let out = match st::wet_steam_properties(x, e) {
                Ok(w) => wet_fields(&w),
                Err(z) => vec![err(&z)],
            };
            m.row("steam::wet_steam_properties", Class::B, &ins, &out);
        }
    }
}

fn json_rows(m: &mut Module, states: &[SaturatedEntry]) {
    let mut sats: Vec<SaturatedEntry> = states.iter().step_by(40).copied().collect();
    if let Some(&base) = states.first() {
        let mut sp = base;
        sp.temperature = f64::NAN;
        sp.pressure = f64::INFINITY;
        sp.v_f = -0.0;
        sp.v_g = f64::MAX;
        sp.h_f = f64::from_bits(1);
        sats.push(sp);
    }
    let mut docs = Vec::new();
    for e in &sats {
        let doc = serde_json::to_string(e).expect("serialize");
        m.row(
            "steam::SaturatedEntry::to_json",
            Class::A,
            &sat_fields(e),
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
    docs.extend(
        [
            "[1,2,3,4,5,6,7,8,9]",
            "[1,2,3,4,5,6,7,8]",
            r#"{"temperature":1,"pressure":2,"v_f":3,"v_g":4,"h_f":5,"h_fg":6,"h_g":7,"s_f":8}"#,
            r#"{"temperature":1e400,"pressure":2,"v_f":3,"v_g":4,"h_f":5,"h_fg":6,"h_g":7,"s_f":8,"s_g":9}"#,
            "[]",
            "null",
            "",
        ]
        .map(String::from),
    );
    for doc in &docs {
        let mut de = serde_json::Deserializer::from_str(doc);
        let out = match SaturatedEntry::deserialize(&mut de).and_then(|x| de.end().map(|()| x)) {
            Ok(x) => sat_fields(&x),
            Err(_) => vec![Tok::Err("Deserialize".to_string())],
        };
        m.row(
            "steam::SaturatedEntry::from_json",
            Class::A,
            &[s(doc.clone())],
            &out,
        );
    }
    // WetSteamProperties
    let mut wets = Vec::new();
    for e in sats.iter().take(4) {
        for x in [0.0, 0.5, 1.0] {
            if let Ok(w) = st::wet_steam_properties(x, e) {
                wets.push(w);
            }
        }
    }
    wets.push(WetSteamProperties {
        quality: f64::NAN,
        specific_volume: f64::INFINITY,
        specific_enthalpy: -0.0,
        specific_entropy: 1e300,
    });
    let mut wdocs = Vec::new();
    for w in &wets {
        let doc = serde_json::to_string(w).expect("serialize");
        m.row(
            "steam::WetSteamProperties::to_json",
            Class::A,
            &wet_fields(w),
            &[s(doc.clone())],
        );
        wdocs.push(doc);
    }
    let wbase = wdocs[1].clone();
    for (_kind, bytes) in variants(&wbase) {
        if let Ok(text) = String::from_utf8(bytes)
            && !text.contains('\0')
        {
            wdocs.push(text);
        }
    }
    wdocs.extend(["[0.5,1,2,3]", "[0.5,1,2]", "{}", "null"].map(String::from));
    for doc in &wdocs {
        let mut de = serde_json::Deserializer::from_str(doc);
        let out = match WetSteamProperties::deserialize(&mut de).and_then(|x| de.end().map(|()| x))
        {
            Ok(x) => wet_fields(&x),
            Err(_) => vec![Tok::Err("Deserialize".to_string())],
        };
        m.row(
            "steam::WetSteamProperties::from_json",
            Class::A,
            &[s(doc.clone())],
            &out,
        );
    }
}

fn sh_fields(e: &SuperheatedEntry) -> Vec<Tok> {
    vec![
        f(e.temperature),
        f(e.pressure),
        f(e.specific_volume),
        f(e.specific_enthalpy),
        f(e.specific_entropy),
    ]
}

fn sh_res(r: ushma::error::Result<SuperheatedEntry>) -> Vec<Tok> {
    match r {
        Ok(e) => sh_fields(&e),
        Err(z) => vec![err(&z)],
    }
}

/// The superheated table's rows, parsed as `steam_data.cyr`'s columns are.
fn sh_table() -> Vec<Vec<f64>> {
    let src = include_str!("../../../rust-old/src/steam.rs");
    crate::t_steam::superheated_rows(src).expect("SUPERHEATED_DATA parses")
}

/// The pressures to try: every tier ± 1 ulp, the geometric midpoint of each pair of tiers, and the range's ends.
fn sh_pressures(rows: &[Vec<f64>]) -> Vec<f64> {
    let mut tiers: Vec<f64> = rows.iter().map(|r| r[1]).collect();
    tiers.dedup();
    let mut ps = Vec::new();
    for (k, &p) in tiers.iter().enumerate() {
        ps.extend(around(p, 1));
        if let Some(&p2) = tiers.get(k + 1) {
            ps.push((p * p2).sqrt());
            ps.push(p + 0.25 * (p2 - p));
        }
    }
    ps.extend([5_000.0, 2e7]);
    ps
}

fn superheated(m: &mut Module) {
    let rows = sh_table();
    let ps = sh_pressures(&rows);
    // superheated_lookup(T, P): at each pressure, every table temperature ± 1 ulp, the midpoints between rows, and
    // around each tier's saturation temperature (the anchor) ± 1 ulp
    let mut temps: Vec<f64> = Vec::new();
    for (k, r) in rows.iter().enumerate() {
        temps.extend(around(r[0], 1));
        if let Some(r2) = rows.get(k + 1)
            && r2[1] == r[1]
        {
            temps.push(0.5 * (r[0] + r2[0]));
        }
    }
    let mut sats: Vec<f64> = Vec::new();
    for &p in &ps {
        if let Ok(e) = st::saturated_by_pressure(p) {
            sats.extend(around(e.temperature, 1));
            sats.push(e.temperature + 1.0);
        }
    }
    temps.sort_by(f64::total_cmp);
    temps.dedup();
    for &p in &ps {
        for &t in temps.iter().chain(sats.iter()) {
            m.row(
                "steam::superheated_lookup",
                Class::C,
                &[f(t), f(p)],
                &sh_res(st::superheated_lookup(t, p)),
            );
        }
    }
    for e in edges() {
        for (t, p) in [(e, 500_000.0), (573.15, e)] {
            m.row(
                "steam::superheated_lookup",
                Class::C,
                &[f(t), f(p)],
                &sh_res(st::superheated_lookup(t, p)),
            );
        }
    }
    let mut r = Rng::new(82);
    for _ in 0..400 {
        let p = r.log_uniform(4.0, 7.0);
        let t = r.uniform(300.0, 1100.0);
        m.row(
            "steam::superheated_lookup",
            Class::C,
            &[f(t), f(p)],
            &sh_res(st::superheated_lookup(t, p)),
        );
    }
    // superheated_by_entropy(P, s): at each pressure, the entropy of every breakpoint the lookup covers there ± 1 ulp,
    // the midpoints, and below and above the range; then the edges and seeded points
    for &p in &ps {
        let mut ss: Vec<f64> = Vec::new();
        for &t in temps.iter().chain(sats.iter()) {
            if let Ok(e) = st::superheated_lookup(t, p) {
                ss.extend(around(e.specific_entropy, 1));
            }
        }
        ss.sort_by(f64::total_cmp);
        ss.dedup();
        let mids: Vec<f64> = ss.windows(2).map(|w| 0.5 * (w[0] + w[1])).collect();
        ss.extend(mids);
        ss.extend([1000.0, 20_000.0]);
        for s0 in ss.iter().step_by(2) {
            m.row(
                "steam::superheated_by_entropy",
                Class::C,
                &[f(p), f(*s0)],
                &sh_res(st::superheated_by_entropy(p, *s0)),
            );
        }
    }
    for e in edges() {
        for (p, s0) in [(e, 8000.0), (100_000.0, e)] {
            m.row(
                "steam::superheated_by_entropy",
                Class::C,
                &[f(p), f(s0)],
                &sh_res(st::superheated_by_entropy(p, s0)),
            );
        }
    }
    for _ in 0..300 {
        let p = r.log_uniform(4.0, 7.0);
        let s0 = r.uniform(6000.0, 9500.0);
        m.row(
            "steam::superheated_by_entropy",
            Class::C,
            &[f(p), f(s0)],
            &sh_res(st::superheated_by_entropy(p, s0)),
        );
    }
    // JSON
    let mut vals: Vec<SuperheatedEntry> = Vec::new();
    for &(t, p) in &[
        (573.15, 500_000.0),
        (573.15, 750_000.0),
        (373.15, 10_000.0),
        (873.15, 10_000_000.0),
    ] {
        if let Ok(e) = st::superheated_lookup(t, p) {
            vals.push(e);
        }
    }
    vals.push(SuperheatedEntry {
        temperature: f64::NAN,
        pressure: f64::INFINITY,
        specific_volume: -0.0,
        specific_enthalpy: f64::MAX,
        specific_entropy: f64::from_bits(1),
    });
    let mut docs = Vec::new();
    for e in &vals {
        let doc = serde_json::to_string(e).expect("serialize");
        m.row(
            "steam::SuperheatedEntry::to_json",
            Class::A,
            &sh_fields(e),
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
        let mut de = serde_json::Deserializer::from_str(doc);
        let out = match SuperheatedEntry::deserialize(&mut de).and_then(|x| de.end().map(|()| x)) {
            Ok(x) => sh_fields(&x),
            Err(_) => vec![Tok::Err("Deserialize".to_string())],
        };
        m.row(
            "steam::SuperheatedEntry::from_json",
            Class::A,
            &[s(doc.clone())],
            &out,
        );
    }
}

pub fn emit() -> Module {
    let mut m = Module::new("steam");
    region4(&mut m);
    let states = lookups(&mut m);
    quality(&mut m, &states);
    json_rows(&mut m, &states);
    superheated(&mut m);
    m
}
