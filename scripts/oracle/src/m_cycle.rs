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
//!
//! M4c: `compare_cycles` (its entries as `v:m`, then each entry's label, result and limits), each limit at the edges and
//! T_H at T_L ± 1 ulp, on every cycle's result and on results relabelled; the T-s and P-v diagrams (a result and
//! `i:points_per_process` in, `v:count` and each point's x, y out) on every cycle's results and on results edited to
//! reach each path's guard: short and empty process lists, an empty and a one-point cycle, isothermal and isenthalpic
//! processes, γ 0 and NaN, zero, negative and non-finite temperatures, volumes and pressures, equal volumes and volumes
//! 1 ulp apart (the isentrope's |ln(V_a/V_b)| > 1e-12 guard). `DiagramPoint`'s and `CycleComparisonEntry`'s JSON.
//! Classes: B for the comparison, C for the diagrams (pow and ln).

use ushma::cycle::{
    self as cy, CycleComparisonEntry, CycleKind, CycleResult, DiagramPoint, ProcessKind, RefrigerationResult, StatePoint,
};
use ushma::steam;

use crate::grids::Rng;
use crate::json_cases::variants;
use crate::sweep::{around, edges, res, row, tuples};
use crate::tsv::{Class, Module, Tok, e, err, f, i, s};

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
    // DiagramPoint
    let pts = [
        DiagramPoint { x: 6_000.5, y: 101_325.0 },
        DiagramPoint { x: f64::NAN, y: f64::INFINITY },
        DiagramPoint { x: -0.0, y: f64::from_bits(1) },
        DiagramPoint { x: 1e300, y: -1e-300 },
    ];
    let mut pdocs = Vec::new();
    for p in &pts {
        let doc = serde_json::to_string(p).expect("serialize");
        m.row("cycle::DiagramPoint::to_json", Class::A, &dp_toks(p), &[s(doc.clone())]);
        pdocs.push(doc);
    }
    let base = pdocs[0].clone();
    for (_kind, bytes) in variants(&base) {
        if let Ok(text) = String::from_utf8(bytes)
            && !text.contains('\0')
        {
            pdocs.push(text);
        }
    }
    pdocs.extend(
        ["[1,2]", "[1]", "[1,2,3]", r#"{"x":1}"#, r#"{"y":2,"x":1}"#, r#"{"x":"1","y":2}"#, "{}", "null"].map(String::from),
    );
    for doc in &pdocs {
        m.row(
            "cycle::DiagramPoint::from_json",
            Class::A,
            &[s(doc.clone())],
            &decode::<DiagramPoint>(doc, dp_toks),
        );
    }
    // CycleComparisonEntry
    let mut cmp = cy::compare_cycles(&[
        (CycleKind::Otto, &otto, 1_800.0, 300.0),
        (CycleKind::Refrigeration, &r1.cycle, 353.0, 313.0),
        (CycleKind::Brayton, &bray, 250.0, 300.0),
    ]);
    let mut odd = cmp[0].clone();
    odd.efficiency = f64::NAN;
    odd.second_law_efficiency = f64::NEG_INFINITY;
    odd.net_work = -0.0;
    cmp.push(odd);
    let mut cdocs = Vec::new();
    for c in &cmp {
        let doc = serde_json::to_string(c).expect("serialize");
        m.row("cycle::CycleComparisonEntry::to_json", Class::A, &cc_toks(c), &[s(doc.clone())]);
        cdocs.push(doc);
    }
    let base = cdocs[0].clone();
    for (_kind, bytes) in variants(&base) {
        if let Ok(text) = String::from_utf8(bytes)
            && !text.contains('\0')
        {
            cdocs.push(text);
        }
    }
    cdocs.extend(
        [
            r#"["Otto",1,2,3,4,5]"#,
            r#"["Otto",1,2,3,4]"#,
            r#"{"kind":"Jet","efficiency":1,"net_work":2,"back_work_ratio":3,"carnot_efficiency":4,"second_law_efficiency":5}"#,
            r#"{"second_law_efficiency":5,"carnot_efficiency":4,"back_work_ratio":3,"net_work":2,"efficiency":1,"kind":"Diesel"}"#,
            r#"{"kind":"Otto","efficiency":1,"net_work":2,"back_work_ratio":3,"carnot_efficiency":4}"#,
            "[]",
            "null",
        ]
        .map(String::from),
    );
    for doc in &cdocs {
        m.row(
            "cycle::CycleComparisonEntry::from_json",
            Class::A,
            &[s(doc.clone())],
            &decode::<CycleComparisonEntry>(doc, cc_toks),
        );
    }
}

fn cc_toks(c: &CycleComparisonEntry) -> Vec<Tok> {
    vec![
        e("CycleKind", c.kind),
        f(c.efficiency),
        f(c.net_work),
        f(c.back_work_ratio),
        f(c.carnot_efficiency),
        f(c.second_law_efficiency),
    ]
}

fn dp_toks(p: &DiagramPoint) -> Vec<Tok> {
    vec![f(p.x), f(p.y)]
}

/// The results the comparison and diagram rows run on: every cycle's, then results edited to reach each guard.
fn sample_results() -> Vec<CycleResult> {
    let otto = cy::otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0).expect("otto");
    let diesel = cy::diesel_cycle(300.0, 101_325.0, 20.0, 2.0, 1.4, 1.0).expect("diesel");
    let bray = cy::brayton_cycle(300.0, 101_325.0, 10.0, 1400.0, 1.4, 1.0).expect("brayton");
    let mono = cy::otto_cycle(300.0, 101_325.0, 10.0, 40_000.0, 5.0 / 3.0, 2.0).expect("otto, monatomic");
    let rk = cy::rankine_cycle(10_000.0, 2_000_000.0, Some(623.15)).expect("rankine");
    let rk_sat = cy::rankine_cycle(10_000.0, 1_000_000.0, None).expect("rankine, saturated");
    let rk_sh = cy::rankine_cycle(100_000.0, 500_000.0, Some(573.15)).expect("rankine, superheated exhaust");
    let rf = cy::refrigeration_cycle(7_384.0, 47_390.0).expect("refrigeration").cycle;
    let mut v = vec![otto.clone(), diesel, bray.clone(), mono, rk, rk_sat, rk_sh, rf];
    let edit = |base: &CycleResult, f: &dyn Fn(&mut CycleResult)| {
        let mut c = base.clone();
        f(&mut c);
        c
    };
    v.push(edit(&otto, &|c| c.processes.truncate(1)));
    v.push(edit(&otto, &|c| c.processes.clear()));
    v.push(edit(&otto, &|c| c.state_points.clear()));
    v.push(edit(&otto, &|c| c.state_points.truncate(1)));
    v.push(edit(&bray, &|c| {
        c.processes = vec![
            ProcessKind::Isothermal,
            ProcessKind::Isenthalpic,
            ProcessKind::Isothermal,
            ProcessKind::Isenthalpic,
        ];
    }));
    v.push(edit(&otto, &|c| c.gamma = 0.0));
    v.push(edit(&otto, &|c| c.gamma = f64::NAN));
    v.push(edit(&otto, &|c| {
        c.state_points[1].temperature = 0.0;
        c.state_points[2].temperature = f64::NAN;
        c.state_points[3].temperature = -5.0;
    }));
    v.push(edit(&bray, &|c| {
        c.state_points[0].volume = 0.0;
        c.state_points[2].pressure = -1.0;
        c.processes[2] = ProcessKind::Isothermal;
    }));
    v.push(edit(&otto, &|c| {
        let v0 = c.state_points[0].volume;
        c.state_points[1].volume = v0;
        c.state_points[3].volume = f64::from_bits(v0.to_bits() + 1);
    }));
    v.push(edit(&bray, &|c| {
        c.state_points[1].volume = f64::INFINITY;
        c.state_points[3].pressure = f64::NAN;
    }));
    v
}

fn comparison_rows(m: &mut Module) {
    let results = sample_results();
    let mut refrig = results[0].clone();
    refrig.kind = CycleKind::Refrigeration;
    refrig.efficiency = 5.0;
    let row_of = |m: &mut Module, entries: &[(CycleKind, &CycleResult, f64, f64)]| {
        let mut ins = vec![Tok::V(entries.len())];
        for (k, r, th, tl) in entries {
            ins.push(e("CycleKind", *k));
            ins.extend(cycle_toks(r));
            ins.extend([f(*th), f(*tl)]);
        }
        let out = cy::compare_cycles(entries);
        let mut outs = vec![Tok::V(out.len())];
        for c in &out {
            outs.extend(cc_toks(c));
        }
        m.row("cycle::compare_cycles", Class::B, &ins, &outs);
    };
    // each result between its own hottest and coldest states, and between 300 and 250 K
    for r in results.iter().take(8).chain([&refrig]) {
        let (th, tl) = match (r.state_points.get(2), r.state_points.first()) {
            (Some(h), Some(c)) => (h.temperature, c.temperature),
            _ => (1000.0, 300.0),
        };
        row_of(m, &[(r.kind, r, th, tl)]);
        row_of(m, &[(r.kind, r, 300.0, 250.0)]);
    }
    // the limits at the edges, and T_H at T_L ± 1 ulp, for a heat engine and for a refrigeration result
    for r in [&results[0], &refrig] {
        for x in crate::sweep::edges() {
            row_of(m, &[(r.kind, r, x, 300.0)]);
            row_of(m, &[(r.kind, r, 1000.0, x)]);
        }
        for th in around(300.0, 1) {
            row_of(m, &[(r.kind, r, th, 300.0)]);
        }
    }
    // a label that differs from the result's own kind; no entries; several at once
    row_of(m, &[(CycleKind::Otto, &results[2], 1400.0, 300.0)]);
    row_of(m, &[]);
    let all: Vec<(CycleKind, &CycleResult, f64, f64)> =
        results.iter().take(8).map(|r| (r.kind, r, 1200.0, 290.0)).collect();
    row_of(m, &all);
}

fn diagram_rows(m: &mut Module) {
    let results = sample_results();
    for (k, r) in results.iter().enumerate() {
        let mut ptss: Vec<usize> = vec![0, 1, 2, 3, 5, 10];
        if k < 8 {
            ptss.push(50);
        }
        for pts in ptss {
            let mut ins = cycle_toks(r);
            ins.push(i(i64::try_from(pts).expect("small")));
            for (name, out) in [
                ("cycle::cycle_ts_diagram", cy::cycle_ts_diagram(r, pts)),
                ("cycle::cycle_pv_diagram", cy::cycle_pv_diagram(r, pts)),
            ] {
                let mut outs = vec![Tok::V(out.len())];
                for p in &out {
                    outs.extend(dp_toks(p));
                }
                m.row(name, Class::C, &ins, &outs);
            }
        }
    }
}

pub fn emit() -> Module {
    let mut m = Module::new("cycle");
    helpers(&mut m);
    gas_cycles(&mut m);
    steam_cycles(&mut m);
    comparison_rows(&mut m);
    diagram_rows(&mut m);
    json_rows(&mut m);
    m
}
