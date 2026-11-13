//! `soorat.tsv` (M7b): Rust 1.4.0's `integration::soorat`, the visualization data built from the M4 and M6 states.
//!
//! A heatmap is written as `v:n` and its values, nx, ny, its origin, spacing and bounds; a profile as `v:n` and its
//! temperatures, Δx, origin, direction and bounds; a cycle diagram as `v:n` T-s points (x, y each), `v:n` P-v points,
//! `v:m` state points (T, P, V, S each), its kind `s:` and efficiency; a network graph as `v:n` node temperatures, `v:m`
//! edges (a, b each) and `v:m` conductances; a flux field as `v:n` vectors (qx, qy each), nx, ny, Δx, Δy and the
//! largest magnitude. Each constructor's row is its source (the M6 grid or network tokens, the M4 cycle result) and
//! arguments in, the value out.
//!
//! The 2D grids are M6b's seven setups before and after their solves, and edits a deserialized grid can hold: missing
//! and short last rows, an extra row, ny past or short of the rows, one-cell-wide axes, the empty grid, NaN, ±∞, all-NaN
//! and signed-zero fields, and Δx, Δy at their edges (the flux). A ragged row other than the last cannot be held by the
//! port's flat grid (ADR 0002 row 2), so none is generated. The 1D grids are M6a's setups, stepped and not, and its
//! malformed grids; the cycles M4c's sample results at 0 to 50 points a process; the networks Rust's test networks and
//! every R edge with temperature lists short, exact and long. JSON: the six types, malformed documents among them.
//! Class B, C for the cycle diagram (its pows and ln), A for JSON.

use ushma::cycle::CycleResult;
use ushma::integration::soorat::{
    CycleDiagramData, CycleStatePoint, HeatFluxField, TemperatureProfile,
    ThermalGridVisualization, ThermalNetworkVisualization,
};
use ushma::numerical::{BoundaryCondition as Bc, ThermalGrid1D, ThermalGrid2D, ThermalNetwork};

use crate::json_cases::variants;
use crate::m_cycle::{cycle_toks, sample_results};
use crate::m_numerical::{g1_toks, g2_toks, malformed, net, net_toks, setups, setups2, usize_tok};
use crate::sweep::edges;
use crate::tsv::{Class, Module, Tok, f, i, s};

fn gv_toks(v: &ThermalGridVisualization) -> Vec<Tok> {
    let mut t = vec![Tok::V(v.values.len())];
    t.extend(v.values.iter().map(|&x| f(x)));
    t.extend([
        usize_tok(v.dimensions[0]),
        usize_tok(v.dimensions[1]),
        f(v.origin[0]),
        f(v.origin[1]),
        f(v.spacing[0]),
        f(v.spacing[1]),
        f(v.min_temp),
        f(v.max_temp),
    ]);
    t
}

fn tp_toks(p: &TemperatureProfile) -> Vec<Tok> {
    let mut t = vec![Tok::V(p.temperatures.len())];
    t.extend(p.temperatures.iter().map(|&x| f(x)));
    t.push(f(p.dx));
    t.extend(p.origin.iter().map(|&x| f(x)));
    t.extend(p.direction.iter().map(|&x| f(x)));
    t.extend([f(p.min_temp), f(p.max_temp)]);
    t
}

fn csp_toks(p: &CycleStatePoint) -> Vec<Tok> {
    vec![
        f(p.temperature),
        f(p.pressure),
        f(p.volume),
        f(p.entropy),
    ]
}

fn pairs(t: &mut Vec<Tok>, ps: &[[f64; 2]]) {
    t.push(Tok::V(ps.len()));
    for p in ps {
        t.extend([f(p[0]), f(p[1])]);
    }
}

fn cd_toks(d: &CycleDiagramData) -> Vec<Tok> {
    let mut t = Vec::new();
    pairs(&mut t, &d.ts_points);
    pairs(&mut t, &d.pv_points);
    t.push(Tok::V(d.state_points.len()));
    for p in &d.state_points {
        t.extend(csp_toks(p));
    }
    t.extend([s(d.kind.clone()), f(d.efficiency)]);
    t
}

fn nv_toks(v: &ThermalNetworkVisualization) -> Vec<Tok> {
    let mut t = vec![Tok::V(v.node_temperatures.len())];
    t.extend(v.node_temperatures.iter().map(|&x| f(x)));
    t.push(Tok::V(v.edges.len()));
    for e in &v.edges {
        t.extend([usize_tok(e[0]), usize_tok(e[1])]);
    }
    t.push(Tok::V(v.conductances.len()));
    t.extend(v.conductances.iter().map(|&x| f(x)));
    t
}

fn hf_toks(h: &HeatFluxField) -> Vec<Tok> {
    let mut t = Vec::new();
    pairs(&mut t, &h.fluxes);
    t.extend([
        usize_tok(h.dimensions[0]),
        usize_tok(h.dimensions[1]),
        f(h.spacing[0]),
        f(h.spacing[1]),
        f(h.max_magnitude),
    ]);
    t
}

fn f64s(xs: &[f64]) -> Vec<Tok> {
    let mut t = vec![Tok::V(xs.len())];
    t.extend(xs.iter().map(|&x| f(x)));
    t
}

/// The 2D grids: M6b's setups, each also solved, then the edits a deserialized grid can hold.
fn grids2() -> Vec<ThermalGrid2D> {
    let mut v = Vec::new();
    for g in setups2() {
        let tol = if g.nx == 11 {
            1e-12
        } else if g.nx == 20 {
            1e-4
        } else {
            1e-6
        };
        let mut solved = g.clone();
        solved.solve_steady_state(tol, 1_000_000).expect("converges");
        v.push(g);
        v.push(solved);
    }
    let base = setups2().remove(1);
    let edit = |f: &dyn Fn(&mut ThermalGrid2D)| {
        let mut g = base.clone();
        f(&mut g);
        g
    };
    v.extend([
        edit(&|g| {
            g.nodes.pop();
        }),
        edit(&|g| g.nodes.truncate(4)),
        edit(&|g| {
            g.nodes.last_mut().expect("rows").truncate(3);
        }),
        edit(&|g| g.nodes.push(vec![123.0; 7])),
        edit(&|g| g.nodes.clear()),
        edit(&|g| g.ny = 9),
        edit(&|g| g.ny = 3),
        edit(&|g| g.ny = 0),
        edit(&|g| {
            g.nx = 1;
            for r in &mut g.nodes {
                r.truncate(1);
            }
        }),
        edit(&|g| {
            g.ny = 1;
            g.nodes.truncate(1);
        }),
        edit(&|g| {
            g.nx = 1;
            g.ny = 1;
            g.nodes = vec![vec![321.0]];
        }),
        edit(&|g| {
            g.nx = 0;
            g.ny = 0;
            g.nodes.clear();
        }),
        edit(&|g| g.nodes[3][3] = f64::NAN),
        edit(&|g| g.nodes[2][4] = f64::INFINITY),
        edit(&|g| g.nodes[0][0] = f64::NEG_INFINITY),
        edit(&|g| {
            for r in &mut g.nodes {
                r.fill(f64::NAN);
            }
        }),
        edit(&|g| {
            for (k, r) in g.nodes.iter_mut().enumerate() {
                r.fill(if k % 2 == 0 { 0.0 } else { -0.0 });
            }
        }),
        edit(&|g| {
            for r in &mut g.nodes {
                r.fill(-0.0);
            }
            g.nodes[3][3] = 0.0;
        }),
        edit(&|g| g.dx = 0.0),
        edit(&|g| g.dx = f64::NAN),
        edit(&|g| g.dx = -0.5),
        edit(&|g| g.dy = f64::INFINITY),
        edit(&|g| g.dy = -0.0),
    ]);
    // Rust's two tests: the 4×4 hot corner, the 5×5 gradient in x
    let mut hot = ThermalGrid2D::new(4, 4, 1.0, 1.0, 300.0).expect("4x4");
    hot.nodes[0][0] = 400.0;
    v.push(hot);
    let mut ramp = ThermalGrid2D::new(5, 5, 1.0, 1.0, 300.0).expect("5x5");
    for row in &mut ramp.nodes {
        for (ix, t) in row.iter_mut().enumerate() {
            *t = 300.0 + ix as f64 * 25.0;
        }
    }
    v.push(ramp);
    v
}

fn heatmap_rows(m: &mut Module) {
    for g in grids2() {
        m.row(
            "soorat::ThermalGridVisualization::from_grid_2d",
            Class::B,
            &g2_toks(&g),
            &gv_toks(&ThermalGridVisualization::from_grid_2d(&g)),
        );
    }
}

fn flux_rows(m: &mut Module) {
    let row = |m: &mut Module, g: &ThermalGrid2D, k: f64| {
        let mut ins = g2_toks(g);
        ins.push(f(k));
        m.row(
            "soorat::HeatFluxField::from_grid_2d",
            Class::B,
            &ins,
            &hf_toks(&HeatFluxField::from_grid_2d(g, k)),
        );
    };
    let gs = grids2();
    for g in &gs {
        row(m, g, 1.0);
        row(m, g, 50.0);
    }
    // the conductivity at its edges, on the gradient grid and a solved setup
    for g in [&gs[gs.len() - 1], &gs[3]] {
        for k in edges() {
            row(m, g, k);
        }
    }
}

fn profile_rows(m: &mut Module) {
    let mut grids: Vec<ThermalGrid1D> = Vec::new();
    for g in setups() {
        let mut stepped = g.clone();
        let dt = 0.5 * stepped.max_stable_dt();
        for _ in 0..5 {
            stepped.step_explicit(dt).expect("steps");
        }
        grids.push(g);
        grids.push(stepped);
    }
    grids.extend(malformed());
    let base = setups().remove(0);
    let edit = |f: &dyn Fn(&mut ThermalGrid1D)| {
        let mut g = base.clone();
        f(&mut g);
        g
    };
    grids.extend([
        edit(&|g| g.nodes.fill(f64::NAN)),
        edit(&|g| {
            for (k, x) in g.nodes.iter_mut().enumerate() {
                *x = if k % 3 == 0 { -0.0 } else { 0.0 };
            }
        }),
        edit(&|g| {
            g.nodes.fill(-0.0);
            g.nodes[7] = 0.0;
        }),
        edit(&|g| g.nodes = vec![f64::INFINITY, 5.0, f64::NEG_INFINITY]),
    ]);
    // Rust's tests: a 5-node grid between 400 K and 300 K, and a 3-node insulated one
    grids.push(
        ThermalGrid1D::new(5, 1.0, 1e-5, 350.0, Bc::Fixed(400.0), Bc::Fixed(300.0)).expect("5"),
    );
    grids.push(
        ThermalGrid1D::new(3, 1.0, 1e-5, 300.0, Bc::Insulated, Bc::Insulated).expect("3"),
    );
    let frames: [([f64; 3], [f64; 3]); 3] = [
        ([0.0; 3], [1.0, 0.0, 0.0]),
        ([1.5, -2.0, 0.25], [0.0, 0.6, 0.8]),
        ([f64::NAN, f64::INFINITY, -0.0], [-0.0, f64::NAN, 1e300]),
    ];
    for (k, g) in grids.iter().enumerate() {
        for (o, d) in frames.iter().take(if k < 2 { 3 } else { 1 }) {
            let mut ins = g1_toks(g);
            ins.extend(o.iter().map(|&x| f(x)));
            ins.extend(d.iter().map(|&x| f(x)));
            m.row(
                "soorat::TemperatureProfile::from_grid_1d",
                Class::B,
                &ins,
                &tp_toks(&TemperatureProfile::from_grid_1d(g, *o, *d)),
            );
        }
    }
}

fn diagram_rows(m: &mut Module) {
    let mut results: Vec<CycleResult> = sample_results();
    // Rust's test: Otto at r = 8 with 2 kJ in
    results.push(ushma::cycle::otto_cycle(300.0, 101_325.0, 8.0, 2000.0, 1.4, 1.0).expect("otto"));
    for (k, r) in results.iter().enumerate() {
        let mut ptss: Vec<usize> = vec![0, 1, 2, 3, 5, 20];
        if k < 8 {
            ptss.push(50);
        }
        for pts in ptss {
            let mut ins = cycle_toks(r);
            ins.push(i(i64::try_from(pts).expect("small")));
            m.row(
                "soorat::CycleDiagramData::from_cycle_result",
                Class::C,
                &ins,
                &cd_toks(&CycleDiagramData::from_cycle_result(r, pts)),
            );
        }
    }
}

fn network_rows(m: &mut Module) {
    let mut nets: Vec<ThermalNetwork> = Vec::new();
    // Rust's tests
    let mut a = ThermalNetwork::new(3);
    a.add_resistance(0, 1, 10.0).expect("r");
    a.add_resistance(1, 2, 20.0).expect("r");
    nets.push(a);
    let mut b = ThermalNetwork::new(3);
    b.add_resistance(0, 1, 2.0).expect("r");
    nets.push(b);
    // a solved chain with fixed ends, an empty network, one with resistances but no nodes, indices past num_nodes,
    // and every R edge
    nets.push(net(
        5,
        &[(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 4, 1.0)],
        &[(0, 400.0), (4, 300.0)],
    ));
    nets.push(ThermalNetwork::new(0));
    nets.push(net(0, &[(0, 1, 1.0)], &[]));
    nets.push(net(2, &[(5, 7, 1.0), (1, 1, 3.0)], &[(9, 1.0)]));
    let mut rs: Vec<(usize, usize, f64)> = Vec::new();
    for (k, r) in edges()
        .into_iter()
        .chain([f64::MAX, f64::MIN_POSITIVE, f64::from_bits(1), 2.2e-308])
        .enumerate()
    {
        rs.push((k % 4, (k + 1) % 4, r));
    }
    nets.push(net(4, &rs, &[]));
    let temps: Vec<Vec<f64>> = vec![
        vec![300.0, 350.0, 400.0],
        vec![300.0],
        vec![],
        vec![400.0, 375.0, 350.0, 325.0, 300.0, 275.0, 250.0],
        vec![f64::NAN, f64::INFINITY, -0.0, 0.0],
    ];
    for x in &nets {
        let mut lists = temps.clone();
        if let Ok(t) = x.solve() {
            lists.push(t);
        }
        for t in &lists {
            let mut ins = net_toks(x);
            ins.extend(f64s(t));
            m.row(
                "soorat::ThermalNetworkVisualization::from_network",
                Class::B,
                &ins,
                &nv_toks(&ThermalNetworkVisualization::from_network(x, t)),
            );
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

/// A type's JSON rows: each value's document, then the decoder over those documents, the variants of `base` (an index
/// into them) and the hand-written ones.
fn json_type<T: serde::Serialize + for<'de> serde::Deserialize<'de>>(
    m: &mut Module,
    name: &str,
    values: &[T],
    toks: &dyn Fn(&T) -> Vec<Tok>,
    base: usize,
    extra: &[String],
) {
    let mut docs = Vec::new();
    for x in values {
        let doc = serde_json::to_string(x).expect("serialize");
        m.row(&format!("soorat::{name}::to_json"), Class::A, &toks(x), &[s(doc.clone())]);
        docs.push(doc);
    }
    let b = docs[base].clone();
    add_variants(&mut docs, &b);
    docs.extend(extra.iter().cloned());
    docs.extend(["[]", "{}", "null", ""].map(String::from));
    for doc in &docs {
        m.row(
            &format!("soorat::{name}::from_json"),
            Class::A,
            &[s(doc.clone())],
            &decode::<T>(doc, toks),
        );
    }
}

fn json_rows(m: &mut Module) {
    let strs = |xs: &[&str]| xs.iter().map(|x| (*x).to_string()).collect::<Vec<String>>();
    // ThermalGridVisualization: Rust's test value, from grids (NaN bounds among them)
    let gs = grids2();
    let mut gvs = vec![ThermalGridVisualization {
        values: vec![300.0; 4],
        dimensions: [2, 2],
        origin: [0.0, 0.0],
        spacing: [1.0, 1.0],
        min_temp: 300.0,
        max_temp: 300.0,
    }];
    for k in [4, 14, 18, 25, 29] {
        gvs.push(ThermalGridVisualization::from_grid_2d(&gs[k]));
    }
    let gv_body = |vals: &str, dims: &str| {
        format!(
            r#"{{"values":{vals},"dimensions":{dims},"origin":[0.0,0.0],"spacing":[1.0,1.0],"min_temp":1.0,"max_temp":2.0}}"#
        )
    };
    let mut extra = vec![
        gv_body("[1,2,3,4]", "[2,2]"),
        gv_body("[1,2,3]", "[2,2]"),
        gv_body("[]", "[0,0]"),
        gv_body("[1,2,3,4]", "[2]"),
        gv_body("[1,2,3,4]", "[2,2,1]"),
        gv_body("[1,2,3,4]", "[2.0,2]"),
        gv_body("[1,2,3,4]", "[-2,2]"),
        gv_body("[1,null,3,4]", "[2,2]"),
        gv_body("[1,2,3,4]", "{}"),
        r#"{"values":[1],"dimensions":[1,1],"origin":[0.0],"spacing":[1.0,1.0],"min_temp":1,"max_temp":1}"#.into(),
        r#"{"values":[1],"dimensions":[1,1],"origin":[0,0],"spacing":[1,1],"min_temp":null,"max_temp":1}"#.into(),
        r#"{"values":[1],"dimensions":[1,1],"origin":[0,0],"spacing":[1,1],"min_temp":1}"#.into(),
        r#"[[1,2],[2,1],[0,0],[1,1],1,2]"#.into(),
        r#"[[1,2],[2,1],[0,0],[1,1],1]"#.into(),
    ];
    json_type(m, "ThermalGridVisualization", &gvs, &gv_toks, 0, &extra);
    // TemperatureProfile
    let g1 = setups();
    let mut tps = vec![TemperatureProfile {
        temperatures: vec![300.0, 310.0, 320.0],
        dx: 0.1,
        origin: [0.0; 3],
        direction: [1.0, 0.0, 0.0],
        min_temp: 300.0,
        max_temp: 320.0,
    }];
    for g in [&g1[2], &g1[5]] {
        tps.push(TemperatureProfile::from_grid_1d(g, [1.0, 2.0, 3.0], [0.0, 0.0, -1.0]));
    }
    let mut nan_g = g1[0].clone();
    nan_g.nodes.fill(f64::NAN);
    tps.push(TemperatureProfile::from_grid_1d(&nan_g, [0.0; 3], [1.0, 0.0, 0.0]));
    let tp = |t: &str, o: &str, d: &str| {
        format!(r#"{{"temperatures":{t},"dx":0.1,"origin":{o},"direction":{d},"min_temp":1,"max_temp":2}}"#)
    };
    extra = vec![
        tp("[1,2]", "[0,0,0]", "[1,0,0]"),
        tp("[]", "[0,0,0]", "[1,0,0]"),
        tp("[1,2]", "[0,0]", "[1,0,0]"),
        tp("[1,2]", "[0,0,0,0]", "[1,0,0]"),
        tp("[1,2]", "[0,0,0]", "[1,0,null]"),
        tp("[1,2]", "[0,0,\"0\"]", "[1,0,0]"),
        tp("{}", "[0,0,0]", "[1,0,0]"),
        r#"[[1,2],0.1,[0,0,0],[1,0,0],1,2]"#.into(),
        r#"{"temperatures":[1],"dx":0.1,"origin":[0,0,0],"direction":[1,0,0],"min_temp":1}"#.into(),
    ];
    json_type(m, "TemperatureProfile", &tps, &tp_toks, 0, &extra);
    // CycleStatePoint
    let csps = [
        CycleStatePoint {
            temperature: 300.0,
            pressure: 100_000.0,
            volume: 0.001,
            entropy: 0.0,
        },
        CycleStatePoint {
            temperature: f64::NAN,
            pressure: -0.0,
            volume: 1e300,
            entropy: f64::from_bits(1),
        },
    ];
    extra = strs(&[
        r#"{"temperature":1,"pressure":2,"volume":3}"#,
        r#"{"temperature":1,"pressure":2,"volume":3,"entropy":"4"}"#,
        r#"{"temperature":1,"pressure":2,"volume":3,"entropy":4,"enthalpy":5}"#,
        r#"[1,2,3,4]"#,
        r#"[1,2,3]"#,
        r#"[1,2,3,4,5]"#,
    ]);
    json_type(m, "CycleStatePoint", &csps, &csp_toks, 0, &extra);
    // CycleDiagramData: Rust's test value, two cycles' data, and kinds needing escapes
    let res = sample_results();
    let mut cds = vec![CycleDiagramData {
        ts_points: vec![[0.0, 300.0], [100.0, 600.0]],
        pv_points: vec![[0.001, 100_000.0], [0.01, 50_000.0]],
        state_points: vec![csps[0]],
        kind: "Otto".to_string(),
        efficiency: 0.56,
    }];
    cds.push(CycleDiagramData::from_cycle_result(&res[0], 3));
    cds.push(CycleDiagramData::from_cycle_result(&res[4], 2));
    let mut odd = cds[0].clone();
    odd.kind = "tab\there \"quoted\" \\ back\u{1}slash é → 𝄞".to_string();
    odd.ts_points.clear();
    cds.push(odd);
    let cd = |ts: &str, sp: &str, kind: &str| {
        format!(r#"{{"ts_points":{ts},"pv_points":[],"state_points":{sp},"kind":{kind},"efficiency":0.5}}"#)
    };
    let sp1 = r#"[{"temperature":1,"pressure":2,"volume":3,"entropy":4}]"#;
    extra = vec![
        cd("[[1,2]]", sp1, "\"Otto\""),
        cd("[[1,2]]", "[[1,2,3,4]]", "\"Carnot\""),
        cd("[[1,2]]", sp1, "\"\""),
        cd("[[1,2]]", sp1, r#""é𝄞\n""#),
        cd("[[1,2]]", sp1, r#""\ud834""#),
        cd("[[1,2]]", sp1, "1"),
        cd("[[1,2]]", sp1, "null"),
        cd("[[1]]", sp1, "\"Otto\""),
        cd("[[1,2,3]]", sp1, "\"Otto\""),
        cd("[1,2]", sp1, "\"Otto\""),
        cd("[[1,2]]", r#"[{"temperature":1,"pressure":2,"volume":3}]"#, "\"Otto\""),
        r#"[[[1,2]],[],[],"Otto",0.5]"#.into(),
        r#"[[[1,2]],[],[],"Otto"]"#.into(),
        r#"{"ts_points":[],"pv_points":[],"state_points":[],"kind":"Otto"}"#.into(),
    ];
    json_type(m, "CycleDiagramData", &cds, &cd_toks, 0, &extra);
    // ThermalNetworkVisualization: Rust's test value, from networks (NaN temperatures, ∞ conductance)
    let mut nvs = vec![ThermalNetworkVisualization {
        node_temperatures: vec![300.0, 350.0, 400.0],
        edges: vec![[0, 1], [1, 2]],
        conductances: vec![10.0, 5.0],
    }];
    let tn = net(3, &[(0, 1, 2.0), (2, 0, 2.2e-308)], &[]);
    nvs.push(ThermalNetworkVisualization::from_network(&tn, &[300.0]));
    nvs.push(ThermalNetworkVisualization::from_network(&ThermalNetwork::new(0), &[]));
    let nv = |t: &str, e: &str, c: &str| {
        format!(r#"{{"node_temperatures":{t},"edges":{e},"conductances":{c}}}"#)
    };
    extra = vec![
        nv("[1,2]", "[[0,1]]", "[0.5]"),
        nv("[1,2]", "[[0,1],[1,0]]", "[0.5]"),
        nv("[1,2]", "[[0,1,2]]", "[0.5]"),
        nv("[1,2]", "[[0]]", "[0.5]"),
        nv("[1,2]", "[[0,-1]]", "[0.5]"),
        nv("[1,2]", "[[0,1.0]]", "[0.5]"),
        nv("[1,2]", "[[0,9223372036854775807]]", "[0.5]"),
        nv("[1,2]", "[[0,9223372036854775808]]", "[0.5]"),
        nv("[1,null]", "[[0,1]]", "[0.5]"),
        nv("[1,2]", "[{\"0\":1}]", "[0.5]"),
        r#"[[1,2],[[0,1]],[0.5]]"#.into(),
        r#"{"node_temperatures":[1,2],"edges":[[0,1]]}"#.into(),
    ];
    json_type(m, "ThermalNetworkVisualization", &nvs, &nv_toks, 0, &extra);
    // HeatFluxField: Rust's test value, from grids
    let mut hfs = vec![HeatFluxField {
        fluxes: vec![[100.0, 0.0], [-50.0, 50.0]],
        dimensions: [2, 1],
        spacing: [0.1, 0.1],
        max_magnitude: 100.0,
    }];
    hfs.push(HeatFluxField::from_grid_2d(&gs[gs.len() - 1], 1.0));
    hfs.push(HeatFluxField::from_grid_2d(&gs[27], 2.0));
    let hf = |fl: &str, d: &str| {
        format!(r#"{{"fluxes":{fl},"dimensions":{d},"spacing":[0.1,0.1],"max_magnitude":1}}"#)
    };
    extra = vec![
        hf("[[1,2],[3,4]]", "[2,1]"),
        hf("[[1,2],[3]]", "[2,1]"),
        hf("[[1,2],[3,4,5]]", "[2,1]"),
        hf("[[1,2],[3,null]]", "[2,1]"),
        hf("[[1,2]]", "[1,1,1]"),
        hf("[[1,2]]", "[1]"),
        hf("[]", "[0,0]"),
        r#"[[[1,2]],[1,1],[0.1,0.1],1]"#.into(),
        r#"{"fluxes":[[1,2]],"dimensions":[1,1],"spacing":[0.1,0.1]}"#.into(),
    ];
    json_type(m, "HeatFluxField", &hfs, &hf_toks, 0, &extra);
}

pub fn emit() -> Module {
    let mut m = Module::new("soorat");
    heatmap_rows(&mut m);
    profile_rows(&mut m);
    diagram_rows(&mut m);
    network_rows(&mut m);
    flux_rows(&mut m);
    json_rows(&mut m);
    m
}
