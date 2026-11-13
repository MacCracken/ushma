//! `numerical.tsv` (M6a onward): Rust 1.4.0's `numerical`.
//!
//! A boundary condition is written as its variant (`e:BoundaryCondition.Fixed`, `.Insulated`, `.Convective`) and the
//! values it carries; a 1D grid as `v:n` and its nodes, then Δx, α and its two conditions. A step's row is a grid state
//! and the step in, the grid after it out (or the refusal), so every step is checked from the oracle's own state.
//!
//! M6a: `ThermalGrid1D::new` over node counts, lengths, diffusivities, initial temperatures and every condition at its
//! edges; `fourier_number` and `max_stable_dt`; `step_explicit`, `step_crank_nicolson` and `step_adaptive` on seven
//! setups (fixed, insulated and convective ends, Biot numbers from 0 to 10, 3 to 50 nodes), each along a trajectory of
//! consecutive states, at fractions and multiples of the stable step, at the explicit stability limit (bisected) ± 2 ulps,
//! at the adaptive exact multiples ± 1 ulp and its sub-step cap (bisected, on the 10-node setup) ± 2 ulps, and on malformed grids (fewer than
//! 3 nodes, Δx or α not positive, a condition out of range, a NaN node). JSON: `BoundaryCondition`, `Side` and
//! `ThermalGrid1D`. Class B throughout (+ − × ÷ and ceil), A for JSON.
//!
//! M6b: a 2D grid is written as its rows (`v:ny`, each `v:len` and its nodes), then nx, ny, Δx, Δy and its four
//! conditions (left, right, bottom, top). `ThermalGrid2D::new` over sizes, extents and initial temperatures at the
//! edges; `set_boundary` on every side with each kind of condition, on a usable and a malformed grid;
//! `solve_steady_state` on seven setups (square and not, every condition, 3×3 to 20×20), each at a generous iteration
//! limit and at exactly the converging count and one fewer (the switch to "did not converge"), the tolerance at the
//! edges, and malformed grids (a missing row, a short last row, nx < 3, Δx not positive, a bad condition, a NaN node).
//! A ragged row other than the last cannot be held by the port's flat grid, so none is generated. JSON:
//! `ThermalGrid2D`, its ragged and mis-shaped documents among them (ADR 0002 row 2).

use ushma::numerical::{
    BoundaryCondition as Bc, Side, ThermalGrid1D, ThermalGrid2D, ThermalNetwork,
};

use crate::grids::Rng;
use crate::json_cases::variants;
use crate::sweep::{around, edges};
use crate::tsv::{Class, Module, Tok, err, f, i, s};

pub(crate) fn bc_toks(b: &Bc) -> Vec<Tok> {
    match *b {
        Bc::Fixed(t) => vec![Tok::E("BoundaryCondition", "Fixed".into()), f(t)],
        Bc::Insulated => vec![Tok::E("BoundaryCondition", "Insulated".into())],
        Bc::Convective { h, t_fluid, k } => vec![
            Tok::E("BoundaryCondition", "Convective".into()),
            f(h),
            f(t_fluid),
            f(k),
        ],
        _ => unreachable!("1.4.0 has three conditions"),
    }
}

pub(crate) fn g1_toks(g: &ThermalGrid1D) -> Vec<Tok> {
    let mut v = vec![Tok::V(g.nodes.len())];
    v.extend(g.nodes.iter().map(|&x| f(x)));
    v.extend([f(g.dx), f(g.alpha)]);
    v.extend(bc_toks(&g.bc_left));
    v.extend(bc_toks(&g.bc_right));
    v
}

pub(crate) fn conv(h: f64, t_fluid: f64, k: f64) -> Bc {
    Bc::Convective { h, t_fluid, k }
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

/// The seven step setups: (nodes, length, α, T₀, left, right), and a profile to start from (None: uniform).
pub(crate) fn setups() -> Vec<ThermalGrid1D> {
    let mut v = Vec::new();
    let mk = |n: usize, l: f64, a: f64, t: f64, bl: Bc, br: Bc| {
        ThermalGrid1D::new(n, l, a, t, bl, br).expect("setup")
    };
    v.push(mk(10, 1.0, 1e-4, 300.0, Bc::Fixed(400.0), Bc::Fixed(300.0)));
    v.push(mk(11, 0.5, 1e-5, 300.0, Bc::Fixed(400.0), Bc::Insulated));
    v.push(mk(
        21,
        0.2,
        4e-5,
        300.0,
        Bc::Fixed(500.0),
        conv(25.0, 300.0, 50.0),
    ));
    v.push(mk(
        8,
        0.01,
        1e-6,
        350.0,
        conv(1000.0, 500.0, 1.0),
        conv(10.0, 250.0, 40.0),
    ));
    let mut e = mk(12, 1.0, 2e-5, 300.0, Bc::Insulated, Bc::Insulated);
    for (k, x) in e.nodes.iter_mut().enumerate() {
        *x = 300.0 + 50.0 * ((k as f64) * 0.7).sin();
    }
    v.push(e);
    v.push(mk(
        3,
        0.1,
        1e-4,
        280.0,
        conv(0.0, 300.0, 10.0),
        Bc::Fixed(320.0),
    ));
    v.push(mk(
        50,
        2.0,
        1e-4,
        300.0,
        Bc::Fixed(373.15),
        Bc::Fixed(273.15),
    ));
    v
}

/// Grids the constructor would refuse, built from Rust's public fields (as a deserialized or edited grid can be).
pub(crate) fn malformed() -> Vec<ThermalGrid1D> {
    let base = setups().remove(0);
    let edit = |f: &dyn Fn(&mut ThermalGrid1D)| {
        let mut g = base.clone();
        f(&mut g);
        g
    };
    vec![
        edit(&|g| g.nodes.truncate(2)),
        edit(&|g| g.nodes.clear()),
        edit(&|g| g.dx = 0.0),
        edit(&|g| g.dx = f64::NAN),
        edit(&|g| g.dx = -1.0),
        edit(&|g| g.alpha = 0.0),
        edit(&|g| g.alpha = f64::INFINITY),
        edit(&|g| g.bc_left = Bc::Fixed(f64::NAN)),
        edit(&|g| g.bc_right = conv(25.0, 300.0, 0.0)),
        edit(&|g| g.bc_right = conv(-1.0, 300.0, 50.0)),
        edit(&|g| g.bc_left = conv(25.0, f64::INFINITY, 50.0)),
        edit(&|g| g.nodes[4] = f64::NAN),
        edit(&|g| g.nodes[0] = f64::INFINITY),
    ]
}

fn grid_rows(m: &mut Module) {
    // new(n, L, α, T₀, left, right)
    let bcs = [
        Bc::Fixed(400.0),
        Bc::Insulated,
        conv(25.0, 300.0, 50.0),
        conv(0.0, 300.0, 50.0),
        Bc::Fixed(f64::NAN),
        Bc::Fixed(f64::INFINITY),
        Bc::Fixed(-0.0),
        conv(-1.0, 300.0, 50.0),
        conv(f64::NAN, 300.0, 50.0),
        conv(25.0, 300.0, 0.0),
        conv(25.0, 300.0, -0.0),
        conv(25.0, f64::NAN, 50.0),
        conv(f64::INFINITY, 300.0, 50.0),
    ];
    let mut news: Vec<(usize, f64, f64, f64, Bc, Bc)> = Vec::new();
    for n in [0, 1, 2, 3, 4, 10] {
        news.push((n, 1.0, 1e-4, 300.0, Bc::Fixed(400.0), Bc::Fixed(300.0)));
    }
    for e in edges() {
        news.push((10, e, 1e-4, 300.0, Bc::Insulated, Bc::Insulated));
        news.push((10, 1.0, e, 300.0, Bc::Insulated, Bc::Insulated));
        news.push((10, 1.0, 1e-4, e, Bc::Insulated, Bc::Insulated));
    }
    for b in &bcs {
        news.push((5, 1.0, 1e-4, 300.0, *b, Bc::Insulated));
        news.push((5, 1.0, 1e-4, 300.0, Bc::Insulated, *b));
    }
    news.push((
        7,
        0.3,
        1e-4,
        300.0,
        Bc::Fixed(f64::NAN),
        conv(-1.0, 1.0, 1.0),
    ));
    for (n, l, a, t, bl, br) in news {
        let mut ins = vec![i(i64::try_from(n).expect("small")), f(l), f(a), f(t)];
        ins.extend(bc_toks(&bl));
        ins.extend(bc_toks(&br));
        let out = match ThermalGrid1D::new(n, l, a, t, bl, br) {
            Ok(g) => g1_toks(&g),
            Err(z) => vec![err(&z)],
        };
        m.row("numerical::ThermalGrid1D::new", Class::B, &ins, &out);
    }
    // fourier_number and max_stable_dt on every setup and malformed grid
    let mut grids = setups();
    grids.extend(malformed());
    for g in &grids {
        for dt in [1.0, 0.5, 1e-3, 0.0, -1.0, f64::NAN, f64::INFINITY] {
            let mut ins = g1_toks(g);
            ins.push(f(dt));
            m.row(
                "numerical::ThermalGrid1D::fourier_number",
                Class::B,
                &ins,
                &[f(g.fourier_number(dt))],
            );
        }
        m.row(
            "numerical::ThermalGrid1D::max_stable_dt",
            Class::B,
            &g1_toks(g),
            &[f(g.max_stable_dt())],
        );
    }
}

/// A step row: the grid and dt in, the grid after the step (and, for the adaptive step, the time advanced) out.
fn step_row(m: &mut Module, which: u8, g: &ThermalGrid1D, dt: f64) -> ThermalGrid1D {
    let mut after = g.clone();
    let (name, r) = match which {
        0 => (
            "numerical::ThermalGrid1D::step_explicit",
            after.step_explicit(dt).map(|()| None),
        ),
        1 => (
            "numerical::ThermalGrid1D::step_crank_nicolson",
            after.step_crank_nicolson(dt).map(|()| None),
        ),
        _ => (
            "numerical::ThermalGrid1D::step_adaptive",
            after.step_adaptive(dt).map(Some),
        ),
    };
    let out = match r {
        Ok(None) => g1_toks(&after),
        Ok(Some(t)) => {
            let mut v = vec![f(t)];
            v.extend(g1_toks(&after));
            v
        }
        Err(z) => vec![err(&z)],
    };
    let mut ins = g1_toks(g);
    ins.push(f(dt));
    m.row(name, Class::B, &ins, &out);
    after
}

fn step_rows(m: &mut Module) {
    for g0 in setups() {
        let dts = g0.max_stable_dt();
        // trajectories: 20 explicit steps at half the stable step, 20 Crank–Nicolson at 10×, 10 adaptive at 3.7×
        let mut g = g0.clone();
        for _ in 0..20 {
            g = step_row(m, 0, &g, 0.5 * dts);
        }
        let mut g = g0.clone();
        for _ in 0..20 {
            g = step_row(m, 1, &g, 10.0 * dts);
        }
        let mut g = g0.clone();
        for _ in 0..10 {
            g = step_row(m, 2, &g, 3.7 * dts);
        }
        // fractions and multiples of the stable step, and dt at the edges
        for k in [0.1, 1.0, 2.0, 100.0, 1e6] {
            // 10⁶ adaptive sub-steps only on the small grids (validate.py recomputes each in Python)
            if k > 1e5 && g0.nodes.len() > 10 {
                continue;
            }
            for which in 0..3 {
                step_row(m, which, &g0, k * dts);
            }
        }
        for e in edges() {
            for which in 0..3 {
                step_row(m, which, &g0, e);
            }
        }
        // the explicit stability limit, bisected, ± 2 ulps
        let unstable = |dt: f64| g0.clone().step_explicit(dt).is_err();
        for dt in around(first_true(dts, 2.0 * dts, &unstable), 2) {
            step_row(m, 0, &g0, dt);
        }
        // the adaptive step at exact multiples of the stable step ± 1 ulp, and its sub-step cap, bisected, ± 2 ulps
        for k in [1.0, 2.0, 3.0] {
            for dt in around(k * dts, 1) {
                step_row(m, 2, &g0, dt);
            }
        }
        let capped = |dt: f64| g0.clone().step_adaptive(dt).is_err();
        if g0.nodes.len() == 10 && !capped(1e5 * dts) && capped(1e7 * dts) {
            for dt in around(first_true(1e5 * dts, 1e7 * dts, &capped), 2) {
                step_row(m, 2, &g0, dt);
            }
        }
        // a huge Crank–Nicolson step: its pivots stop being finite
        step_row(m, 1, &g0, 1e300);
    }
    for g in malformed() {
        let dts = setups()[0].max_stable_dt();
        for which in 0..3 {
            step_row(m, which, &g, 0.5 * dts);
        }
    }
}

/// A usize as an `i:` token, or, past i64 (only a document's 2^63 variant reaches that; bayan refuses it, ADR 5), its
/// exact decimal as an `s:` token.
pub(crate) fn usize_tok(n: usize) -> Tok {
    match i64::try_from(n) {
        Ok(x) => i(x),
        Err(_) => s(n.to_string()),
    }
}

pub(crate) fn g2_toks(g: &ThermalGrid2D) -> Vec<Tok> {
    let mut v = vec![Tok::V(g.nodes.len())];
    for row in &g.nodes {
        v.push(Tok::V(row.len()));
        v.extend(row.iter().map(|&x| f(x)));
    }
    v.extend([usize_tok(g.nx), usize_tok(g.ny), f(g.dx), f(g.dy)]);
    for b in [&g.bc_left, &g.bc_right, &g.bc_bottom, &g.bc_top] {
        v.extend(bc_toks(b));
    }
    v
}

fn side_tok(sd: Side) -> Tok {
    Tok::E("Side", format!("{sd:?}"))
}

/// The seven solve setups.
pub(crate) fn setups2() -> Vec<ThermalGrid2D> {
    let mk = |nx: usize, ny: usize, lx: f64, ly: f64, t: f64, bcs: &[(Side, Bc)]| {
        let mut g = ThermalGrid2D::new(nx, ny, lx, ly, t).expect("setup");
        for &(sd, b) in bcs {
            g.set_boundary(sd, b);
        }
        g
    };
    let fixed4 = [
        (Side::Left, Bc::Fixed(400.0)),
        (Side::Right, Bc::Fixed(400.0)),
        (Side::Top, Bc::Fixed(400.0)),
        (Side::Bottom, Bc::Fixed(400.0)),
    ];
    vec![
        mk(5, 5, 1.0, 1.0, 300.0, &fixed4),
        mk(
            7,
            7,
            1.0,
            1.0,
            350.0,
            &[
                (Side::Left, Bc::Fixed(400.0)),
                (Side::Right, Bc::Fixed(300.0)),
            ],
        ),
        mk(
            11,
            5,
            0.5,
            0.2,
            350.0,
            &[
                (Side::Left, Bc::Fixed(400.0)),
                (Side::Right, conv(10.0, 300.0, 2.0)),
            ],
        ),
        mk(
            6,
            8,
            0.3,
            0.7,
            320.0,
            &[
                (Side::Right, Bc::Fixed(280.0)),
                (Side::Bottom, conv(15.0, 350.0, 1.5)),
                (Side::Top, Bc::Fixed(310.0)),
            ],
        ),
        mk(3, 3, 1.0, 1.0, 300.0, &[(Side::Left, Bc::Fixed(500.0))]),
        mk(
            20,
            20,
            1.0,
            1.0,
            300.0,
            &[
                (Side::Left, Bc::Fixed(400.0)),
                (Side::Right, Bc::Fixed(300.0)),
            ],
        ),
        mk(
            6,
            6,
            0.1,
            0.1,
            300.0,
            &[
                (Side::Left, conv(100.0, 400.0, 2.0)),
                (Side::Right, conv(5.0, 250.0, 0.5)),
                (Side::Bottom, conv(0.0, 300.0, 1.0)),
                (Side::Top, Bc::Fixed(360.0)),
            ],
        ),
    ]
}

fn solve_row(m: &mut Module, g: &ThermalGrid2D, tol: f64, max: usize) {
    let mut after = g.clone();
    let out = match after.solve_steady_state(tol, max) {
        Ok(n) => {
            let mut v = vec![i(i64::try_from(n).expect("small"))];
            v.extend(g2_toks(&after));
            v
        }
        Err(z) => vec![err(&z)],
    };
    let mut ins = g2_toks(g);
    ins.extend([f(tol), i(i64::try_from(max).expect("small"))]);
    m.row(
        "numerical::ThermalGrid2D::solve_steady_state",
        Class::B,
        &ins,
        &out,
    );
}

fn grid2d_rows(m: &mut Module) {
    // new(nx, ny, lx, ly, T₀)
    let mut news: Vec<(usize, usize, f64, f64, f64)> = Vec::new();
    for nx in [0, 2, 3, 4, 10] {
        for ny in [2, 3, 7] {
            news.push((nx, ny, 1.0, 1.0, 300.0));
        }
    }
    for e in edges() {
        news.push((5, 5, e, 1.0, 300.0));
        news.push((5, 5, 1.0, e, 300.0));
        news.push((5, 5, 1.0, 1.0, e));
    }
    for (nx, ny, lx, ly, t) in news {
        let ins = [
            i(i64::try_from(nx).expect("small")),
            i(i64::try_from(ny).expect("small")),
            f(lx),
            f(ly),
            f(t),
        ];
        let out = match ThermalGrid2D::new(nx, ny, lx, ly, t) {
            Ok(g) => g2_toks(&g),
            Err(z) => vec![err(&z)],
        };
        m.row("numerical::ThermalGrid2D::new", Class::B, &ins, &out);
    }
    // set_boundary(side, bc) on a usable grid and on one missing a row
    let base = ThermalGrid2D::new(5, 4, 1.0, 0.6, 300.0).expect("base");
    let mut popped = base.clone();
    popped.nodes.pop();
    let bcs = [
        Bc::Fixed(400.0),
        Bc::Insulated,
        conv(25.0, 350.0, 2.0),
        Bc::Fixed(f64::NAN),
    ];
    for g in [&base, &popped] {
        for sd in [Side::Left, Side::Right, Side::Top, Side::Bottom] {
            for b in &bcs {
                let mut after = g.clone();
                after.set_boundary(sd, *b);
                let mut ins = g2_toks(g);
                ins.push(side_tok(sd));
                ins.extend(bc_toks(b));
                m.row(
                    "numerical::ThermalGrid2D::set_boundary",
                    Class::B,
                    &ins,
                    &g2_toks(&after),
                );
            }
        }
    }
    // solve_steady_state: each setup at a generous limit, at its converging count and one fewer, tolerances at the edges
    for g in setups2() {
        let tol = if g.nx == 11 {
            1e-12
        } else if g.nx == 20 {
            1e-4
        } else {
            1e-6
        };
        let n = g
            .clone()
            .solve_steady_state(tol, 1_000_000)
            .expect("converges");
        solve_row(m, &g, tol, 1_000_000);
        solve_row(m, &g, tol, n);
        solve_row(m, &g, tol, n - 1);
        solve_row(m, &g, tol, 0);
        solve_row(m, &g, tol, 1);
        if g.nx <= 7 {
            for e in edges() {
                solve_row(m, &g, e, 50);
            }
        }
    }
    // malformed grids
    let ok = setups2().remove(1);
    let edit = |f: &dyn Fn(&mut ThermalGrid2D)| {
        let mut g = ok.clone();
        f(&mut g);
        g
    };
    for g in [
        edit(&|g| {
            g.nodes.pop();
        }),
        edit(&|g| {
            g.nodes.last_mut().expect("rows").pop();
        }),
        edit(&|g| g.nodes.push(vec![300.0; 7])),
        edit(&|g| {
            g.nx = 2;
            for r in &mut g.nodes {
                r.truncate(2);
            }
        }),
        edit(&|g| g.dx = 0.0),
        edit(&|g| g.dy = f64::NAN),
        edit(&|g| g.bc_left = Bc::Fixed(f64::NAN)),
        edit(&|g| g.bc_top = conv(25.0, 300.0, 0.0)),
        edit(&|g| g.bc_bottom = conv(-1.0, 300.0, 1.0)),
        edit(&|g| g.nodes[3][3] = f64::NAN),
        edit(&|g| g.nodes[2][4] = f64::INFINITY),
    ] {
        solve_row(m, &g, 1e-6, 1000);
    }
}

pub(crate) fn net_toks(x: &ThermalNetwork) -> Vec<Tok> {
    let mut v = vec![usize_tok(x.num_nodes), Tok::V(x.resistances.len())];
    for &(a, b, r) in &x.resistances {
        v.extend([usize_tok(a), usize_tok(b), f(r)]);
    }
    v.push(Tok::V(x.fixed.len()));
    for &(node, t) in &x.fixed {
        v.extend([usize_tok(node), f(t)]);
    }
    v
}

/// A network with these lists as they stand, unchecked (Rust's fields are public).
pub(crate) fn net(n: usize, res: &[(usize, usize, f64)], fixed: &[(usize, f64)]) -> ThermalNetwork {
    let mut x = ThermalNetwork::new(n);
    x.resistances.extend_from_slice(res);
    x.fixed.extend_from_slice(fixed);
    x
}

fn below(rng: &mut Rng, k: usize) -> usize {
    usize::try_from(rng.int(0, i64::try_from(k).expect("small") - 1)).expect("small")
}

/// A connected network as probe 9 draws them: a chain through every node plus extra edges, R log-uniform over
/// 1e-3 … 1e14, one to three fixed nodes at 200 … 1000 K. `loops` also draws self-loops and repeated edges.
fn seeded(rng: &mut Rng, max_n: usize, loops: bool) -> ThermalNetwork {
    let n = 2 + below(rng, max_n - 1);
    let mut x = ThermalNetwork::new(n);
    for k in 1..n {
        let r = rng.log_uniform(-3.0, 14.0);
        x.resistances.push((k - 1, k, r));
    }
    for _ in 0..below(rng, 2 * n) {
        let a = below(rng, n);
        let b = below(rng, n);
        if a != b || loops {
            let r = rng.log_uniform(-3.0, 14.0);
            x.resistances.push((a, b, r));
            if loops && below(rng, 4) == 0 {
                x.resistances.push((b, a, r));
            }
        }
    }
    let nfix = 1 + below(rng, 3.min(n - 1));
    for _ in 0..nfix {
        let node = below(rng, n);
        if loops || x.fixed.iter().all(|&(q, _)| q != node) {
            let t = rng.uniform(200.0, 1000.0);
            x.fixed.push((node, t));
        }
    }
    x
}

fn net_solve_row(m: &mut Module, x: &ThermalNetwork) {
    let out = match x.solve() {
        Ok(t) => {
            let mut v = vec![Tok::V(t.len())];
            v.extend(t.iter().map(|&y| f(y)));
            v
        }
        Err(z) => vec![err(&z)],
    };
    m.row(
        "numerical::ThermalNetwork::solve",
        Class::B,
        &net_toks(x),
        &out,
    );
}

fn network_rows(m: &mut Module) {
    let big = usize::try_from(i64::MAX).expect("64-bit");
    for n in [0, 1, 2, 3, 40, 1000, big] {
        m.row(
            "numerical::ThermalNetwork::new",
            Class::A,
            &[usize_tok(n)],
            &net_toks(&ThermalNetwork::new(n)),
        );
    }
    // add_resistance(a, b, R) and set_fixed_temperature(node, T) at every index edge and every value edge
    let bases = [
        ThermalNetwork::new(0),
        ThermalNetwork::new(3),
        net(4, &[(0, 1, 2.0)], &[(0, 300.0)]),
    ];
    let mut vals = edges();
    vals.extend([
        1.0,
        300.0,
        1e14,
        5e-324,
        f64::MIN_POSITIVE,
        f64::MAX,
        f64::MIN,
    ]);
    for x in &bases {
        let n = x.num_nodes;
        let mut idx = vec![0, 1, n.saturating_sub(1), n, n + 1, big];
        idx.sort_unstable();
        idx.dedup();
        let mut adds: Vec<(usize, usize, f64)> = Vec::new();
        for &a in &idx {
            for &b in &idx {
                adds.push((a, b, 1.0));
            }
        }
        for &r in &vals {
            adds.push((0, n.saturating_sub(1), r));
        }
        for (a, b, r) in adds {
            let mut after = x.clone();
            let out = match after.add_resistance(a, b, r) {
                Ok(()) => net_toks(&after),
                Err(z) => vec![err(&z)],
            };
            let mut ins = net_toks(x);
            ins.extend([usize_tok(a), usize_tok(b), f(r)]);
            m.row(
                "numerical::ThermalNetwork::add_resistance",
                Class::A,
                &ins,
                &out,
            );
        }
        let mut sets: Vec<(usize, f64)> = idx.iter().map(|&k| (k, 300.0)).collect();
        sets.extend(vals.iter().map(|&t| (n.saturating_sub(1), t)));
        for (node, t) in sets {
            let mut after = x.clone();
            let out = match after.set_fixed_temperature(node, t) {
                Ok(()) => net_toks(&after),
                Err(z) => vec![err(&z)],
            };
            let mut ins = net_toks(x);
            ins.extend([usize_tok(node), f(t)]);
            m.row(
                "numerical::ThermalNetwork::set_fixed_temperature",
                Class::A,
                &ins,
                &out,
            );
        }
    }
    // solve: the Rust tests' networks, garbage lists, structure edges, overflow, then seeded networks
    let series =
        |r1: f64, r2: f64, t0: f64, t2: f64| net(3, &[(0, 1, r1), (1, 2, r2)], &[(0, t0), (2, t2)]);
    let mut nets = vec![
        series(1.0, 1.0, 400.0, 300.0),
        net(
            3,
            &[(0, 1, 2.0), (1, 2, 2.0), (0, 2, 2.0)],
            &[(0, 400.0), (2, 300.0)],
        ),
        net(3, &[(0, 1, 1e-3), (1, 2, 1e-3), (0, 2, 3e-4)], &[]),
        net(4, &[(0, 1, 1.0), (2, 3, 1.0)], &[(0, 300.0)]),
        series(1e14, 1e14, 400.0, 300.0),
        ThermalNetwork::new(0),
        net(2, &[(0, 5, 1.0)], &[]),
        net(2, &[(0, 1, 1.0)], &[(0, 300.0), (5, 300.0)]),
        net(2, &[(0, 1, 1.0)], &[(0, f64::NAN)]),
        net(2, &[(0, 1, 1.0)], &[(0, f64::INFINITY)]),
        net(2, &[(1, 1, 0.0)], &[(0, 300.0)]),
        net(2, &[(0, 1, -1.0)], &[(0, 300.0)]),
        net(2, &[(0, 1, f64::NAN)], &[(0, 300.0)]),
        net(2, &[(0, 1, f64::INFINITY)], &[(0, 300.0)]),
        net(2, &[(0, 1, 1.0)], &[(big, 300.0)]),
        // every node fixed, no free unknowns; a later entry for a node wins
        net(1, &[], &[(0, 300.0)]),
        net(
            3,
            &[(0, 1, 1.0)],
            &[(0, 300.0), (1, 310.0), (2, 320.0), (1, 330.0)],
        ),
        net(3, &[], &[(0, 300.0), (1, 310.0), (2, 320.0)]),
        // a lone free node, unanchored; a free node with no edge beside an anchored one
        ThermalNetwork::new(1),
        net(3, &[(0, 1, 1.0)], &[(0, 300.0)]),
        // self-loops (free and fixed), repeated edges, an edge between two fixed nodes
        net(
            3,
            &[(0, 1, 1.0), (1, 1, 3.0), (1, 2, 1.0)],
            &[(0, 400.0), (2, 300.0)],
        ),
        net(
            3,
            &[(0, 0, 1.0), (0, 1, 1.0), (1, 2, 1.0), (2, 2, 7.0)],
            &[(0, 400.0), (2, 300.0)],
        ),
        net(3, &[(0, 1, 1e20), (1, 2, 1e20), (1, 1, 1.0)], &[(0, 400.0)]),
        net(3, &[(0, 1, 0.1), (1, 1, 1e-17), (1, 2, 0.3)], &[(0, 400.0)]),
        net(
            4,
            &[
                (0, 1, 2.0),
                (0, 1, 2.0),
                (1, 0, 4.0),
                (1, 2, 1.0),
                (2, 3, 1.0),
            ],
            &[(0, 400.0), (3, 300.0)],
        ),
        net(
            4,
            &[(0, 3, 1.0), (0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0)],
            &[(0, 400.0), (3, 300.0)],
        ),
        net(3, &[(1, 0, 1.0), (2, 1, 1.0)], &[(2, 300.0), (0, 400.0)]),
        net(
            5,
            &[(4, 3, 1.0), (3, 2, 2.0), (2, 1, 3.0), (1, 0, 4.0)],
            &[(4, 500.0)],
        ),
        // overflow: a conductance 1/R past f64::MAX, a right-hand side g·T past it, a subnormal conductance
        series(1e-310, 1.0, 400.0, 300.0),
        series(5e-324, 1.0, 400.0, 300.0),
        series(f64::MIN_POSITIVE, 1.0, 400.0, 300.0),
        series(0.5, 1.0, f64::MAX, 300.0),
        series(0.5, 0.5, f64::MAX, f64::MAX),
        series(1e308, 1.0, 400.0, 300.0),
        series(f64::MAX, f64::MAX, 400.0, 300.0),
        net(
            3,
            &[(0, 1, 6e-309), (0, 1, 6e-309), (1, 2, 1.0)],
            &[(0, 400.0), (2, 300.0)],
        ),
        net(2, &[(0, 1, 1e-308)], &[(0, 400.0)]),
        net(2, &[(0, 1, 1e-308)], &[(0, 1e300)]),
        net(3, &[(0, 1, 1.0), (1, 2, 1.0)], &[(0, -0.0), (2, -0.0)]),
        net(3, &[(0, 1, 1.0), (1, 2, 1.0)], &[(0, -0.0), (2, 0.0)]),
    ];
    let mut vals2 = edges();
    vals2.extend([1e14, 1e15, 1e16, 1e-300, 5e-324, f64::MAX]);
    for &r in &vals2 {
        nets.push(series(r, 1.0, 400.0, 300.0));
        nets.push(series(1.0, 1.0, r, 300.0));
    }
    let mut rng = Rng::new(0x006e_6574_776f_726b);
    for k in 0..400 {
        nets.push(seeded(
            &mut rng,
            if k % 8 == 0 { 120 } else { 40 },
            k % 4 == 3,
        ));
    }
    for x in &nets {
        net_solve_row(m, x);
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

fn json_rows(m: &mut Module) {
    // BoundaryCondition
    let bcs = [
        Bc::Fixed(400.0),
        Bc::Fixed(f64::NAN),
        Bc::Fixed(-0.0),
        Bc::Insulated,
        conv(25.0, 300.0, 50.0),
        conv(f64::INFINITY, f64::NAN, 0.0),
    ];
    let mut docs = Vec::new();
    for b in &bcs {
        let doc = serde_json::to_string(b).expect("serialize");
        m.row(
            "numerical::BoundaryCondition::to_json",
            Class::A,
            &bc_toks(b),
            &[s(doc.clone())],
        );
        docs.push(doc);
    }
    let base = docs[4].clone();
    add_variants(&mut docs, &base);
    docs.extend(
        [
            "\"Fixed\"",
            r#"{"Fixed":null}"#,
            r#"{"Fixed":"1"}"#,
            r#"{"Insulated":null}"#,
            r#"{"Insulated":1}"#,
            "\"insulated\"",
            r#"{"Convective":[1,2,3]}"#,
            r#"{"Convective":[1,2]}"#,
            r#"{"Convective":{"h":1,"t_fluid":2}}"#,
            r#"{"Convective":{"h":1,"t_fluid":2,"k":3,"x":4}}"#,
            r#"{"Fixed":1,"Insulated":null}"#,
            r#"{"Radiative":1}"#,
            "[]",
            "{}",
            "null",
        ]
        .map(String::from),
    );
    for doc in &docs {
        m.row(
            "numerical::BoundaryCondition::from_json",
            Class::A,
            &[s(doc.clone())],
            &decode::<Bc>(doc, bc_toks),
        );
    }
    // Side
    for sd in [Side::Left, Side::Right, Side::Top, Side::Bottom] {
        m.row(
            "numerical::Side::to_json",
            Class::A,
            &[Tok::E("Side", format!("{sd:?}"))],
            &[s(serde_json::to_string(&sd).expect("serialize"))],
        );
    }
    for doc in [
        "\"Left\"",
        "\"Bottom\"",
        "\"left\"",
        "\"Front\"",
        r#"{"Top":null}"#,
        "1",
        "",
    ] {
        m.row(
            "numerical::Side::from_json",
            Class::A,
            &[s(doc)],
            &decode::<Side>(doc, |x: &Side| vec![Tok::E("Side", format!("{x:?}"))]),
        );
    }
    // ThermalGrid1D
    let mut grids = setups();
    grids.truncate(4);
    let mut odd = grids[0].clone();
    odd.nodes[2] = f64::NAN;
    odd.nodes.truncate(3);
    grids.push(odd);
    let mut empty = grids[1].clone();
    empty.nodes.clear();
    grids.push(empty);
    let mut gdocs = Vec::new();
    for g in &grids {
        let doc = serde_json::to_string(g).expect("serialize");
        m.row(
            "numerical::ThermalGrid1D::to_json",
            Class::A,
            &g1_toks(g),
            &[s(doc.clone())],
        );
        gdocs.push(doc);
    }
    let base = gdocs[0].clone();
    add_variants(&mut gdocs, &base);
    gdocs.extend(
        [
            r#"[[1,2,3],0.5,0.0001,"Insulated",{"Fixed":1}]"#,
            r#"{"nodes":[1,2,3],"dx":0.5,"alpha":0.0001,"bc_left":"Insulated"}"#,
            r#"{"nodes":3,"dx":0.5,"alpha":0.0001,"bc_left":"Insulated","bc_right":"Insulated"}"#,
            r#"{"nodes":[1,"2",3],"dx":0.5,"alpha":0.0001,"bc_left":"Insulated","bc_right":"Insulated"}"#,
            r#"{"nodes":[1,2,3],"dx":0.5,"alpha":0.0001,"bc_left":"Open","bc_right":"Insulated"}"#,
            r#"{"nodes":[],"dx":-1,"alpha":0,"bc_left":{"Fixed":1},"bc_right":{"Convective":[1,2,3]}}"#,
            "[]",
            "null",
        ]
        .map(String::from),
    );
    for doc in &gdocs {
        m.row(
            "numerical::ThermalGrid1D::from_json",
            Class::A,
            &[s(doc.clone())],
            &decode::<ThermalGrid1D>(doc, g1_toks),
        );
    }
    // ThermalGrid2D
    let mut g2s = setups2();
    g2s.truncate(5);
    let mut solved = g2s[1].clone();
    solved.solve_steady_state(1e-6, 100_000).expect("solves");
    g2s.push(solved);
    let mut odd = g2s[0].clone();
    odd.nodes[2][2] = f64::NAN;
    g2s.push(odd);
    let mut hdocs = Vec::new();
    for g in &g2s {
        let doc = serde_json::to_string(g).expect("serialize");
        m.row(
            "numerical::ThermalGrid2D::to_json",
            Class::A,
            &g2_toks(g),
            &[s(doc.clone())],
        );
        hdocs.push(doc);
    }
    let base = hdocs[4].clone();
    add_variants(&mut hdocs, &base);
    let bcs = r#""bc_left":"Insulated","bc_right":"Insulated","bc_bottom":"Insulated","bc_top":{"Fixed":1}"#;
    for body in [
        r#""nodes":[[1,2,3],[4,5,6],[7,8,9]],"nx":3,"ny":3"#,
        r#""nodes":[[1,2,3],[4,5],[7,8,9]],"nx":3,"ny":3"#,
        r#""nodes":[[1,2,3],[4,5,6]],"nx":3,"ny":3"#,
        r#""nodes":[[1,2,3,4],[5,6,7,8],[9,10,11,12]],"nx":3,"ny":4"#,
        r#""nodes":[[1,2,3],[4,5,6],[7,8,9]],"nx":3.0,"ny":3"#,
        r#""nodes":[[1,2,3],[4,5,6],[7,8,9]],"nx":-3,"ny":3"#,
        r#""nodes":[[1,2,3],[4,5,6],[7,8,9]],"ny":3"#,
        r#""nodes":[[1,"2",3],[4,5,6],[7,8,9]],"nx":3,"ny":3"#,
        r#""nodes":[1,2,3],"nx":3,"ny":1"#,
        r#""nodes":[],"nx":0,"ny":0"#,
    ] {
        hdocs.push(format!(r#"{{{body},"dx":0.5,"dy":0.5,{bcs}}}"#));
    }
    hdocs.push(r#"[[[1,2,3],[4,5,6],[7,8,9]],3,3,0.5,0.5,"Insulated","Insulated","Insulated","Insulated"]"#.to_string());
    hdocs.extend(["[]", "null"].map(String::from));
    for doc in &hdocs {
        m.row(
            "numerical::ThermalGrid2D::from_json",
            Class::A,
            &[s(doc.clone())],
            &decode::<ThermalGrid2D>(doc, g2_toks),
        );
    }
    // ThermalNetwork
    let nets = [
        ThermalNetwork::new(0),
        net(3, &[(0, 1, 1.0), (1, 2, 1.0)], &[(0, 400.0), (2, 300.0)]),
        net(4, &[(0, 1, 2.5e-7), (3, 3, 1e14)], &[(1, -0.0)]),
        net(
            2,
            &[(0, 5, f64::NAN), (1, 0, f64::INFINITY)],
            &[(7, f64::NEG_INFINITY)],
        ),
    ];
    let mut ndocs = Vec::new();
    for x in &nets {
        let doc = serde_json::to_string(x).expect("serialize");
        m.row(
            "numerical::ThermalNetwork::to_json",
            Class::A,
            &net_toks(x),
            &[s(doc.clone())],
        );
        ndocs.push(doc);
    }
    let base = ndocs[1].clone();
    add_variants(&mut ndocs, &base);
    for body in [
        r#""num_nodes":3,"resistances":[[0,1,1]],"fixed":[[0,300]]"#,
        r#""num_nodes":3,"resistances":[[0,1]],"fixed":[]"#,
        r#""num_nodes":3,"resistances":[[0,1,1.0,2]],"fixed":[]"#,
        r#""num_nodes":3,"resistances":[[0,-1,1.0]],"fixed":[]"#,
        r#""num_nodes":3,"resistances":[[0,1.0,1.0]],"fixed":[]"#,
        r#""num_nodes":3,"resistances":[[0,1,"1"]],"fixed":[]"#,
        r#""num_nodes":3,"resistances":[[0,1,null]],"fixed":[]"#,
        r#""num_nodes":3,"resistances":[],"fixed":[[0]]"#,
        r#""num_nodes":3,"resistances":[],"fixed":[[0,300,1]]"#,
        r#""num_nodes":3,"resistances":[],"fixed":[{"0":300}]"#,
        r#""num_nodes":3,"resistances":{},"fixed":[]"#,
        r#""num_nodes":-1,"resistances":[],"fixed":[]"#,
        r#""num_nodes":3.0,"resistances":[],"fixed":[]"#,
        r#""num_nodes":3,"resistances":[]"#,
        r#""num_nodes":3,"resistances":[],"fixed":[],"extra":1"#,
        r#""num_nodes":9223372036854775807,"resistances":[[9223372036854775806,0,1e-300]],"fixed":[]"#,
    ] {
        ndocs.push(format!("{{{body}}}"));
    }
    ndocs.extend(
        [
            r#"[3,[[0,1,1.0]],[[0,300.0]]]"#,
            r#"[3,[[0,1,1.0]]]"#,
            "[]",
            "null",
        ]
        .map(String::from),
    );
    for doc in &ndocs {
        m.row(
            "numerical::ThermalNetwork::from_json",
            Class::A,
            &[s(doc.clone())],
            &decode::<ThermalNetwork>(doc, net_toks),
        );
    }
}

pub fn emit() -> Module {
    let mut m = Module::new("numerical");
    grid_rows(&mut m);
    step_rows(&mut m);
    grid2d_rows(&mut m);
    network_rows(&mut m);
    json_rows(&mut m);
    m
}
