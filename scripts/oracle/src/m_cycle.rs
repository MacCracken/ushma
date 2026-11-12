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
//!
//! M4b: `rankine_cycle` (its `Option<f64>` superheat as `n:` or `f:`) over condenser × boiler pressures × superheats,
//! with every branch switch tried either side: the boiler at the condenser pressure ± 1 ulp, the table's ends ± 1 ulp,
//! each tier temperature and each tier's saturation temperature ± 1 ulp, and the superheat where the turbine exhaust
//! leaves the wet region (s₃ = s_g at the condenser, found by bisection) ± 2 ulps. `refrigeration_cycle` over
//! evaporator × condenser pressures, with the condenser at the evaporator ± 1 ulp and the evaporator where the
//! compressor exit enters the superheated table (bisected) ± 2 ulps. `heat_pump_cop` seeded and at the edges.
//! `RefrigerationResult`'s JSON. The steam cycles are class C (the superheated lookups' ln and exp between tiers).

use ushma::cycle::{self as cy, CycleKind, CycleResult, ProcessKind, RefrigerationResult, StatePoint};
use ushma::steam;

use crate::grids::Rng;
use crate::json_cases::variants;
use crate::sweep::{around, edges, res, row, tuples};
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

/// A RefrigerationResult as its tokens: its cycle's, then the two COPs.
fn refr_toks(r: &RefrigerationResult) -> Vec<Tok> {
    let mut v = cycle_toks(&r.cycle);
    v.push(f(r.cop_refrigeration));
    v.push(f(r.cop_heat_pump));
    v
}

fn refr_res(r: ushma::error::Result<RefrigerationResult>) -> Vec<Tok> {
    match r {
        Ok(c) => refr_toks(&c),
        Err(z) => vec![err(&z)],
    }
}

/// The superheated table's tier pressures (1.4.0's private `SUPERHEAT_PRESSURES`, as `t_steam` parses them).
fn tiers() -> Vec<f64> {
    let src = include_str!("../../../rust-old/src/steam.rs");
    let mut t: Vec<f64> = crate::t_steam::superheated_rows(src)
        .expect("SUPERHEATED_DATA parses")
        .iter()
        .map(|r| r[1])
        .collect();
    t.dedup();
    t
}

/// The saturated table's first and last pressures (its range).
fn sat_range() -> (f64, f64) {
    let src = include_str!("../../../rust-old/src/steam.rs");
    let rows = crate::t_steam::saturated_rows(src).expect("SATURATED_TABLE parses");
    (rows[0][1], rows[rows.len() - 1][1])
}

/// The smallest positive double in (lo, hi] where `up` holds, for `up` false at lo, true at hi and switching once
/// between: a bisection on the bit patterns, which order the positive doubles.
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

fn rankine_row(m: &mut Module, pc: f64, pb: f64, t: Option<f64>) {
    let sh = match t {
        Some(x) => f(x),
        None => Tok::None,
    };
    m.row(
        "cycle::rankine_cycle",
        Class::C,
        &[f(pc), f(pb), sh],
        &cyc_res(cy::rankine_cycle(pc, pb, t)),
    );
}

fn steam_cycles(m: &mut Module) {
    let tiers = tiers();
    let (p_min, p_max) = sat_range();
    let sat = |p: f64| steam::saturated_by_pressure(p).expect("in the table");
    // rankine_cycle(P_cond, P_boil, T_superheat): the grid
    let mut pbs: Vec<f64> = Vec::new();
    for (k, &p) in tiers.iter().enumerate() {
        pbs.push(p);
        if let Some(&p2) = tiers.get(k + 1) {
            pbs.push((p * p2).sqrt());
        }
    }
    pbs.extend([2e7, p_max]);
    let pcs = [p_min, 2_000.0, 7_384.0, 10_000.0, 47_390.0, 100_000.0, 500_000.0];
    for &pc in &pcs {
        for &pb in &pbs {
            rankine_row(m, pc, pb, None);
            for &t in &[373.15, 473.15, 573.15, 673.15, 873.15, 1073.15] {
                rankine_row(m, pc, pb, Some(t));
            }
        }
    }
    // each tier temperature ± 1 ulp, on three tiers; each tier's saturation temperature ± 1 ulp
    let src = include_str!("../../../rust-old/src/steam.rs");
    let mut temps: Vec<f64> = crate::t_steam::superheated_rows(src)
        .expect("SUPERHEATED_DATA parses")
        .iter()
        .map(|r| r[0])
        .collect();
    temps.sort_by(f64::total_cmp);
    temps.dedup();
    for &pb in &[500_000.0, 2_000_000.0, 10_000_000.0] {
        for &t in &temps {
            for x in around(t, 1) {
                rankine_row(m, 10_000.0, pb, Some(x));
            }
        }
    }
    for &pb in &tiers {
        let pc = if pb > 10_000.0 { 10_000.0 } else { 5_000.0 };
        for x in around(sat(pb).temperature, 1) {
            rankine_row(m, pc, pb, Some(x));
        }
    }
    // the boiler at the condenser pressure ± 1 ulp; the table's ends ± 1 ulp
    for &pc in &[7_384.0, 100_000.0] {
        for pb in around(pc, 1) {
            rankine_row(m, pc, pb, None);
            rankine_row(m, pc, pb, Some(873.15));
        }
    }
    for pc in around(p_min, 1) {
        rankine_row(m, pc, 1e6, None);
    }
    for pb in around(p_max, 1) {
        rankine_row(m, 10_000.0, pb, None);
    }
    // the turbine exhaust leaving the wet region: s₃(T) = s_g(P_cond), ± 2 ulps
    for &(pc, pb) in &[
        (10_000.0, 2e6),
        (100_000.0, 500_000.0),
        (47_390.0, 1e6),
        (20_000.0, 5e6),
        (7_384.0, 500_000.0),
        (100_000.0, 1e7),
    ] {
        let sg = sat(pc).s_g;
        let wet = |t: f64| match steam::superheated_lookup(t, pb) {
            Ok(e) => e.specific_entropy > sg,
            Err(_) => false,
        };
        let lo = sat(pb).temperature;
        // the highest table temperature the lookup takes at P_boil (a tier's range, or the narrower of two)
        let hi = temps
            .iter()
            .rev()
            .copied()
            .find(|&t| steam::superheated_lookup(t, pb).is_ok())
            .unwrap_or(lo);
        if !wet(lo) && wet(hi) {
            for x in around(first_true(lo, hi, &wet), 2) {
                rankine_row(m, pc, pb, Some(x));
            }
        }
    }
    // seeded draws, then each argument of two base cases at the edges
    let mut r = crate::grids::Rng::new(97);
    for _ in 0..300 {
        let pc = r.log_uniform(2.8, 6.3);
        let pb = pc * 10f64.powf(r.uniform(0.0, 2.5));
        let t = if r.unit() < 0.33 {
            None
        } else {
            Some(r.uniform(300.0, 1100.0))
        };
        rankine_row(m, pc, pb, t);
    }
    for (pc, pb, t) in [(10_000.0, 2e6, Some(573.15)), (10_000.0, 1e6, None)] {
        for e in edges() {
            rankine_row(m, e, pb, t);
            rankine_row(m, pc, e, t);
            if t.is_some() {
                rankine_row(m, pc, pb, Some(e));
            }
        }
    }
    // refrigeration_cycle(P_evap, P_cond)
    let mut rf: Vec<Vec<f64>> = Vec::new();
    let pes = [p_min, 1_000.0, 2_000.0, 3_169.0, 5_000.0, 7_384.0, 10_000.0, 20_000.0, 47_390.0];
    let mut pcs2: Vec<f64> = pbs[..pbs.len() - 2].to_vec();
    pcs2.extend([3_169.0, 7_384.0, 20_000.0, 47_390.0, 101_325.0, 2e7]);
    for &pe in &pes {
        for &pc in &pcs2 {
            rf.push(vec![pe, pc]);
        }
    }
    for &pe in &[3_169.0, 7_384.0] {
        for pc in around(pe, 1) {
            rf.push(vec![pe, pc]);
        }
    }
    for pe in around(p_min, 1) {
        rf.push(vec![pe, 100_000.0]);
    }
    // the evaporator where the compressor exit enters the superheated table, ± 2 ulps
    for &pc in &[100_000.0, 200_000.0, 500_000.0, 1e6] {
        let ok = |pe: f64| cy::refrigeration_cycle(pe, pc).is_ok();
        let (lo, hi) = (p_min, 0.5 * pc);
        if !ok(lo) && ok(hi) {
            for x in around(first_true(lo, hi, &ok), 2) {
                rf.push(vec![x, pc]);
            }
        }
    }
    rf.extend(tuples(
        &[vec![7_384.0, 47_390.0], vec![3_169.0, 100_000.0]],
        300,
        98,
        &|r| {
            let pe = r.log_uniform(2.8, 5.0);
            vec![pe, pe * 10f64.powf(r.uniform(0.0, 3.0))]
        },
    ));
    for a in rf {
        m.row(
            "cycle::refrigeration_cycle",
            Class::C,
            &[f(a[0]), f(a[1])],
            &refr_res(cy::refrigeration_cycle(a[0], a[1])),
        );
    }
    // heat_pump_cop(COP_ref)
    for a in tuples(&[vec![3.0]], 200, 99, &|r| vec![r.uniform(0.0, 20.0)]) {
        row(m, "cycle::heat_pump_cop", Class::B, &a, vec![f(cy::heat_pump_cop(a[0]))]);
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
    // RefrigerationResult
    let r1 = cy::refrigeration_cycle(7_384.0, 47_390.0).expect("refrigeration");
    let r2 = cy::refrigeration_cycle(3_169.0, 100_000.0).expect("refrigeration");
    let mut odd = r1.clone();
    odd.cycle.state_points.clear();
    odd.cycle.processes.truncate(1);
    odd.cop_refrigeration = f64::NAN;
    odd.cop_heat_pump = f64::NEG_INFINITY;
    let mut fdocs = Vec::new();
    for r in [&r1, &r2, &odd] {
        let doc = serde_json::to_string(r).expect("serialize");
        m.row(
            "cycle::RefrigerationResult::to_json",
            Class::A,
            &refr_toks(r),
            &[s(doc.clone())],
        );
        fdocs.push(doc);
    }
    let base = fdocs[0].clone();
    for (_kind, bytes) in variants(&base) {
        if let Ok(text) = String::from_utf8(bytes)
            && !text.contains('\0')
        {
            fdocs.push(text);
        }
    }
    let seq = r#"["Refrigeration",[[1,2,3,4,5]],["Isenthalpic"],1,2,-3,4,1,0]"#;
    fdocs.extend([
        format!(r#"{{"cycle":{seq},"cop_refrigeration":1.5,"cop_heat_pump":2.5}}"#),
        format!(r#"[{seq},1.5,2.5]"#),
        format!(r#"{{"cop_heat_pump":2.5,"cop_refrigeration":1.5,"cycle":{seq}}}"#),
        format!(r#"{{"cycle":{seq},"cop_refrigeration":1.5}}"#),
        format!(r#"[{seq},1.5]"#),
        r#"{"cycle":null,"cop_refrigeration":1.5,"cop_heat_pump":2.5}"#.to_string(),
        format!(r#"{{"cycle":{},"cop_refrigeration":1.5,"cop_heat_pump":2.5}}"#, rdocs[0]),
        "[]".to_string(),
        "null".to_string(),
    ]);
    for doc in &fdocs {
        m.row(
            "cycle::RefrigerationResult::from_json",
            Class::A,
            &[s(doc.clone())],
            &decode::<RefrigerationResult>(doc, refr_toks),
        );
    }
}

pub fn emit() -> Module {
    let mut m = Module::new("cycle");
    helpers(&mut m);
    gas_cycles(&mut m);
    steam_cycles(&mut m);
    json_rows(&mut m);
    m
}
