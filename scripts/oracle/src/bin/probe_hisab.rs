//! `probe_hisab` — what Rust hisab 1.4.0 itself returns, for M0 probes 1 and 9 (port plan §4 step 9; compared by
//! `scripts/probes/p1_p9_hisab.cyr`). It writes `<out>/simpson.txt` and `<out>/gauss.txt` (default `build/probe-hisab`,
//! gitignored), every value as its 16-hex-digit bit pattern:
//! - `simpson a b n result`: `hisab::calc::integral_simpson` of x⁴/(1 + x)² (+ − × ÷ only, so the Cyrius integrand
//!   gives the same bits and any difference is the integrator's) over seeded intervals and step counts, debye_cv's
//!   200 among them;
//! - `gauss nf a… x…`: ushma's `ThermalNetwork::solve` system for a seeded network (the conductance matrix and
//!   right-hand side, each row scaled by its diagonal, as 1.4.0 builds it), row-major nf × (nf + 1), then
//!   `hisab::num::gaussian_elimination`'s solution. The same network is also solved through ushma's public API and the
//!   two solutions must be bit-identical, which proves this file's copy of the construction;
//! - `gref nf a…`: a system hisab 1.4.0 refuses (its pivot test), with no solution, for the singularity decisions.
//!
//! ```sh
//! cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml --bin probe_hisab
//! ```

use std::fmt::Write as _;
use std::process::ExitCode;

/// xorshift64*, as the oracle's grids (fixed seeds: the same data on every run).
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }

    fn unit(&mut self) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let top = (self.next() >> 11) as f64;
        top * (1.0 / 9_007_199_254_740_992.0)
    }

    fn below(&mut self, n: usize) -> usize {
        usize::try_from(self.next() % (n as u64)).unwrap_or(0)
    }
}

fn hex(x: f64) -> String {
    format!("{:016x}", x.to_bits())
}

fn integrand(x: f64) -> f64 {
    let y = 1.0 + x;
    (x * x * x * x) / (y * y)
}

/// A network's nodes, resistances (a, b, R) and fixed temperatures.
struct Net {
    n: usize,
    res: Vec<(usize, usize, f64)>,
    fixed: Vec<(usize, f64)>,
}

/// ushma 1.4.0's `ThermalNetwork::solve` system, row-scaled (numerical.rs), and the free nodes' indices.
fn system(net: &Net) -> (Vec<Vec<f64>>, Vec<usize>) {
    let mut is_fixed = vec![false; net.n];
    let mut fixed_t = vec![0.0; net.n];
    for &(node, t) in &net.fixed {
        is_fixed[node] = true;
        fixed_t[node] = t;
    }
    let free: Vec<usize> = (0..net.n).filter(|i| !is_fixed[*i]).collect();
    let nf = free.len();
    let mut map = vec![usize::MAX; net.n];
    for (fi, &ni) in free.iter().enumerate() {
        map[ni] = fi;
    }
    let mut m = vec![vec![0.0; nf + 1]; nf];
    for &(a, b, r) in &net.res {
        let g = 1.0 / r;
        if !is_fixed[a] && !is_fixed[b] {
            let (fa, fb) = (map[a], map[b]);
            m[fa][fa] += g;
            m[fa][fb] -= g;
            m[fb][fb] += g;
            m[fb][fa] -= g;
        } else if !is_fixed[a] && is_fixed[b] {
            let fa = map[a];
            m[fa][fa] += g;
            m[fa][nf] += g * fixed_t[b];
        } else if is_fixed[a] && !is_fixed[b] {
            let fb = map[b];
            m[fb][fb] += g;
            m[fb][nf] += g * fixed_t[a];
        }
    }
    for (fi, row) in m.iter_mut().enumerate() {
        let diag = row[fi];
        if diag.is_finite() && diag > 0.0 {
            row.iter_mut().for_each(|v| *v /= diag);
        }
    }
    (m, free)
}

/// A connected network: a chain through every node plus extra edges, resistances log-uniform over 1e-3 … 1e14
/// (the 1.4.0 regression's R ≈ 1e14 included), one to three fixed nodes at 200 … 1000 K.
fn network(rng: &mut Rng) -> Net {
    let n = 2 + rng.below(39);
    let mut res = Vec::new();
    let r_of = |rng: &mut Rng| {
        let e = -3.0 + 17.0 * rng.unit();
        let s = format!("{:.6}e{}", 1.0 + 9.0 * rng.unit(), e.floor());
        s.parse::<f64>().unwrap_or(1.0)
    };
    for k in 1..n {
        let r = r_of(rng);
        res.push((k - 1, k, r));
    }
    for _ in 0..rng.below(2 * n) {
        let a = rng.below(n);
        let b = rng.below(n);
        if a != b {
            let r = r_of(rng);
            res.push((a, b, r));
        }
    }
    let mut fixed = Vec::new();
    let nfix = 1 + rng.below(3.min(n - 1));
    for _ in 0..nfix {
        let node = rng.below(n);
        if fixed.iter().all(|&(f, _)| f != node) {
            fixed.push((node, 200.0 + 800.0 * rng.unit()));
        }
    }
    Net { n, res, fixed }
}

fn run() -> Result<(), String> {
    let out = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "build/probe-hisab".to_string());
    std::fs::create_dir_all(&out).map_err(|e| format!("{out}: {e}"))?;

    let mut rng = Rng(0x51_3950_0a1d);
    let mut s = String::new();
    for k in 0..3000 {
        let a = if k % 3 == 0 { 0.0 } else { 10.0 * rng.unit() };
        let b = a + 50.0 * rng.unit() + 1e-3;
        let n = [200, 200, 200, 2, 4, 10, 7, 1000][rng.below(8)];
        let r = hisab::calc::integral_simpson(integrand, a, b, n)
            .map_err(|e| format!("simpson: {e}"))?;
        let _ = writeln!(s, "simpson {} {} {n} {}", hex(a), hex(b), hex(r));
    }
    std::fs::write(format!("{out}/simpson.txt"), &s).map_err(|e| e.to_string())?;

    let mut g = String::new();
    let mut rows = 0;
    let mut mirrored = 0;
    let mut refused = 0;
    while rows < 2000 {
        let net = network(&mut rng);
        let (mut m, free) = system(&net);
        if free.is_empty() {
            continue;
        }
        let aug = m.clone();
        let Ok(x) = hisab::num::gaussian_elimination(&mut m) else {
            // A system 1.4.0 refuses (its pivot test): recorded for the singularity decisions, without a solution.
            let _ = write!(g, "gref {}", free.len());
            for row in &aug {
                for &v in row {
                    let _ = write!(g, " {}", hex(v));
                }
            }
            g.push('\n');
            refused += 1;
            continue;
        };
        let mut tn = ushma::numerical::ThermalNetwork::new(net.n);
        for &(a, b, r) in &net.res {
            tn.add_resistance(a, b, r)
                .map_err(|e| format!("add_resistance: {e}"))?;
        }
        for &(node, t) in &net.fixed {
            tn.set_fixed_temperature(node, t)
                .map_err(|e| format!("set_fixed: {e}"))?;
        }
        let api = tn.solve().map_err(|e| format!("solve: {e}"))?;
        for (fi, &ni) in free.iter().enumerate() {
            if api[ni].to_bits() != x[fi].to_bits() {
                return Err(format!(
                    "network {rows}: the mirrored system's solution differs from solve()"
                ));
            }
        }
        mirrored += 1;
        let _ = write!(g, "gauss {}", free.len());
        for row in &aug {
            for &v in row {
                let _ = write!(g, " {}", hex(v));
            }
        }
        for &v in &x {
            let _ = write!(g, " {}", hex(v));
        }
        g.push('\n');
        rows += 1;
    }
    std::fs::write(format!("{out}/gauss.txt"), &g).map_err(|e| e.to_string())?;
    println!(
        "probe_hisab: {out}/simpson.txt 3000 rows, {out}/gauss.txt {rows} rows ({mirrored} mirrored bit-for-bit against \
         ThermalNetwork::solve) and {refused} systems hisab 1.4.0 refuses"
    );
    Ok(())
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("probe_hisab: {e}");
            ExitCode::from(2)
        }
    }
}
