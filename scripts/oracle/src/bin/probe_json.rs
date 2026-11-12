//! `probe_json` — serde_json 1.0.154 on ushma 1.4.0's 28 serializable types, for M0 probe 11 (port plan §4 step 9;
//! read back by `scripts/probes/p11_json.cyr`). It writes, under `<out>` (default `build/probe-json`, gitignored):
//! - `cases.txt`, one line per document: `type kind syntax typed rt doc stream back`, tab-separated:
//!   - `kind`: `base` (serde_json's own text for a value), or a variant made from a base document (malformed,
//!     duplicate, unknown and missing fields, other number, string and shape forms);
//!   - `syntax`: 1 when serde_json accepts the text as JSON at all (`from_str::<IgnoredAny>`), else 0;
//!   - `typed`: serde_json's verdict reading the text as the type: `ok`, or `err:<syntax|data|eof>:<message>`;
//!   - `rt`: for `ok`, 1 when the value read writes back as the base document's bytes, else 0 (`-` otherwise);
//!   - `doc`: the document bytes, hex;
//!   - `stream`: the document as a stream of space-separated tokens when `syntax` is 1 (`{` `}` `[` `]`, `k:` and
//!     `s:` keys and strings percent-encoded, `f:` a JSON float's binary64 bits as Rust parses it, `i:` an
//!     integer's text, `b:0`/`b:1`, `n`), else `-`. The Cyrius probe flattens bayan's tree the same way;
//!   - `back`: for `ok`, what serde_json wrote back for the value it read, hex (the decoded value), else `-`.
//! - `f64.txt`, one line per double: `bits serde_json display back` (serde_json's text, Rust's `{}`, and the bits
//!   serde_json 1.0.154 reads its own text back as, `-` for `null`: its default float parser is best-effort, not
//!   correctly rounded, so `back` is not always `bits`), over every f64
//!   in the base documents, every power of two, the powers of ten and their neighbours, the printers' form
//!   thresholds, small integers and a seeded sweep of random bit patterns and short decimals.
//!
//! ```sh
//! cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml --bin probe_json
//! ```

use std::fmt::Write as _;
use std::process::ExitCode;

use serde::Serialize;
use serde::de::{DeserializeOwned, IgnoredAny};

use ushma::ai::{DaimonConfig, HooshConfig};
use ushma::chem::ALL_SPECIES;
use ushma::cycle::{
    self, CycleComparisonEntry, CycleKind, CycleResult, DiagramPoint, ProcessKind,
    RefrigerationResult, StatePoint,
};
use ushma::integration::soorat::{
    CycleDiagramData, CycleStatePoint, HeatFluxField, TemperatureProfile, ThermalGridVisualization,
    ThermalNetworkVisualization,
};
use ushma::material::{ALL_MATERIALS, ThermalMaterial};
use ushma::numerical::{BoundaryCondition, Side, ThermalGrid1D, ThermalGrid2D, ThermalNetwork};
use ushma::phase::{ALL_SUBSTANCES, Phase, SubstanceData};
use ushma::state::{ALL_GASES, GasData};
use ushma::steam::{self, SaturatedEntry, SuperheatedEntry, WetSteamProperties};

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
}

#[path = "../json_cases.rs"]
mod json_cases;
use json_cases::{hex, pct, stream, variants};

#[derive(Default)]
struct Out {
    cases: String,
    count: usize,
    floats: Vec<f64>,
}

/// serde_json's verdict on `doc` as a T: (`ok` or the error, the rt flag, the value read written back as hex).
fn verdict<T: DeserializeOwned + Serialize>(doc: &[u8], base: &str) -> (String, String, String) {
    match serde_json::from_slice::<T>(doc) {
        Ok(v) => {
            let back = serde_json::to_string(&v).unwrap_or_default();
            let rt = if back == base { "1" } else { "0" };
            ("ok".into(), rt.into(), hex(back.as_bytes()))
        }
        Err(e) => (err(&e), "-".into(), "-".into()),
    }
}

fn err(e: &serde_json::Error) -> String {
    let c = match e.classify() {
        serde_json::error::Category::Io => "io",
        serde_json::error::Category::Syntax => "syntax",
        serde_json::error::Category::Data => "data",
        serde_json::error::Category::Eof => "eof",
    };
    format!("err:{c}:{}", pct(e.to_string().as_bytes()))
}

fn line<T: DeserializeOwned + Serialize>(
    o: &mut Out,
    ty: &str,
    kind: &str,
    doc: &[u8],
    base: &str,
) {
    let syntax = serde_json::from_slice::<IgnoredAny>(doc).is_ok();
    let st = if syntax {
        match stream(doc) {
            Some(s) => s,
            None => panic!("{ty} {kind}: serde_json accepts a document this tokenizer refuses"),
        }
    } else {
        "-".into()
    };
    let (typed, rtv, back) = verdict::<T>(doc, base);
    let _ = writeln!(
        o.cases,
        "{ty}\t{kind}\t{}\t{typed}\t{rtv}\t{}\t{st}\t{back}",
        u8::from(syntax),
        hex(doc)
    );
    o.count += 1;
}

fn emit<T: DeserializeOwned + Serialize>(o: &mut Out, ty: &str, values: &[T]) {
    for (n, v) in values.iter().enumerate() {
        let base = serde_json::to_string(v).expect("serialize");
        line::<T>(o, ty, "base", base.as_bytes(), &base);
        if let Some(s) = stream(base.as_bytes()) {
            for t in s.split(' ') {
                if let Some(b) = t.strip_prefix("f:") {
                    o.floats
                        .push(f64::from_bits(u64::from_str_radix(b, 16).expect("hex")));
                }
            }
        }
        if n == 0 {
            for (kind, doc) in variants(&base) {
                line::<T>(o, ty, &kind, &doc, &base);
            }
        }
    }
}

fn floats(seen: &[f64]) -> String {
    let mut xs: Vec<f64> = seen.to_vec();
    xs.extend([
        0.0,
        -0.0,
        f64::MAX,
        f64::MIN,
        f64::MIN_POSITIVE,
        f64::EPSILON,
        f64::NAN,
        f64::INFINITY,
    ]);
    xs.extend([
        f64::NEG_INFINITY,
        f64::from_bits(1),
        f64::from_bits(0x000F_FFFF_FFFF_FFFF),
    ]);
    for k in -1074..=1023 {
        xs.push(2f64.powi(k));
    }
    for e in -324..=308 {
        let p: f64 = format!("1e{e}").parse().expect("pow10");
        xs.extend([
            p,
            f64::from_bits(p.to_bits() + 1),
            f64::from_bits(p.to_bits().saturating_sub(1)),
        ]);
        for m in [1.5, 2.5, 5.0, 9.999_999_999_999_998] {
            xs.push(m * p);
        }
    }
    for t in [
        "1e15",
        "1e16",
        "1e17",
        "9007199254740992",
        "9007199254740993",
        "1e-5",
        "1e-6",
        "1e-7",
        "0.0001",
    ] {
        let x: f64 = t.parse().expect("threshold");
        xs.extend([
            x,
            f64::from_bits(x.to_bits() + 1),
            f64::from_bits(x.to_bits() - 1),
            -x,
        ]);
    }
    for i in 0..=1000 {
        xs.push(f64::from(i));
        xs.push(f64::from(i) / 10.0);
        xs.push(f64::from(i) / 100.0);
    }
    let mut r = Rng(0x9e37_79b9_7f4a_7c15);
    for _ in 0..100_000 {
        xs.push(f64::from_bits(r.next()));
    }
    for _ in 0..50_000 {
        let digits = (r.next() % 17) as u32 + 1;
        let m = r.next() % 10u64.pow(digits);
        let e = (r.next() % 60) as i32 - 30;
        let x: f64 = format!("{m}e{e}").parse().expect("decimal");
        xs.push(if r.next() & 1 == 1 { -x } else { x });
    }
    let mut s = String::new();
    for x in xs {
        let j = serde_json::to_string(&x).expect("f64");
        // serde_json's own reading of its text (its default, best-effort float parser), `-` for `null`.
        let back = match serde_json::from_str::<f64>(&j) {
            Ok(b) => format!("{:016x}", b.to_bits()),
            Err(_) => "-".to_string(),
        };
        let _ = writeln!(s, "{:016x} {j} {x} {back}", x.to_bits());
    }
    s
}

fn main() -> ExitCode {
    let dir = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "build/probe-json".to_string());
    let mut o = Out::default();

    let otto = cycle::otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0).expect("otto");
    let diesel = cycle::diesel_cycle(300.0, 101_325.0, 20.0, 2.0, 1.4, 1.0).expect("diesel");
    let brayton = cycle::brayton_cycle(300.0, 101_325.0, 10.0, 1400.0, 1.4, 1.0).expect("brayton");
    let rankine = cycle::rankine_cycle(10_000.0, 2_000_000.0, Some(573.15)).expect("rankine");
    let rankine_sat = cycle::rankine_cycle(10_000.0, 2_000_000.0, None).expect("rankine");
    let refr = cycle::refrigeration_cycle(7_384.0, 47_390.0).expect("refrigeration");
    let cycles = vec![
        otto.clone(),
        diesel.clone(),
        brayton.clone(),
        rankine.clone(),
        rankine_sat,
        refr.cycle.clone(),
    ];

    let mut g1 = ThermalGrid1D::new(
        6,
        0.1,
        1.0e-5,
        300.0,
        BoundaryCondition::Fixed(400.0),
        BoundaryCondition::Convective {
            h: 10.0,
            t_fluid: 293.15,
            k: 50.0,
        },
    )
    .expect("grid1d");
    for _ in 0..5 {
        g1.step_explicit(1.0).expect("step");
    }
    let mut g2 = ThermalGrid2D::new(4, 3, 0.1, 0.1, 300.0).expect("grid2d");
    g2.set_boundary(Side::Left, BoundaryCondition::Fixed(400.0));
    g2.set_boundary(Side::Top, BoundaryCondition::Insulated);
    g2.set_boundary(
        Side::Right,
        BoundaryCondition::Convective {
            h: 25.0,
            t_fluid: 290.0,
            k: 45.0,
        },
    );
    g2.solve_steady_state(1.0e-9, 10_000).expect("steady");
    let mut net = ThermalNetwork::new(4);
    net.add_resistance(0, 1, 0.5).expect("r");
    net.add_resistance(1, 2, 0.25).expect("r");
    net.add_resistance(2, 3, 1.0e-3).expect("r");
    net.set_fixed_temperature(0, 400.0).expect("fix");
    net.set_fixed_temperature(3, 300.0).expect("fix");
    let temps = net.solve().expect("solve");

    let sat: Vec<SaturatedEntry> = [273.16, 373.15, 450.0, 600.0]
        .iter()
        .map(|&t| steam::saturated_by_temperature(t).expect("sat"))
        .collect();
    let wet: Vec<WetSteamProperties> = [0.0, 0.5, 1.0]
        .iter()
        .map(|&q| steam::wet_steam_properties(q, &sat[1]).expect("wet"))
        .collect();
    let sup: Vec<SuperheatedEntry> = [(573.15, 1.0e6), (773.15, 8.0e6)]
        .iter()
        .map(|&(t, p)| steam::superheated_lookup(t, p).expect("sup"))
        .collect();

    emit(
        &mut o,
        "DaimonConfig",
        &[
            DaimonConfig::default(),
            DaimonConfig {
                endpoint: "https://h\u{e9}st:9090/\"x\"".into(),
                api_key: Some("k\\e\ny\u{1}".into()),
            },
        ],
    );
    emit(&mut o, "HooshConfig", &[HooshConfig::default()]);
    emit(
        &mut o,
        "Species",
        &ALL_SPECIES.iter().map(|s| (*s).clone()).collect::<Vec<_>>(),
    );
    emit::<StatePoint>(
        &mut o,
        "StatePoint",
        &cycles
            .iter()
            .flat_map(|c| c.state_points.clone())
            .collect::<Vec<_>>(),
    );
    emit(
        &mut o,
        "ProcessKind",
        &[
            ProcessKind::Isentropic,
            ProcessKind::Isochoric,
            ProcessKind::Isobaric,
            ProcessKind::Isothermal,
            ProcessKind::Isenthalpic,
        ],
    );
    emit(
        &mut o,
        "CycleKind",
        &[
            CycleKind::Otto,
            CycleKind::Diesel,
            CycleKind::Brayton,
            CycleKind::Rankine,
            CycleKind::Refrigeration,
        ],
    );
    emit::<CycleResult>(&mut o, "CycleResult", &cycles);
    let mut dp = cycle::cycle_ts_diagram(&otto, 3);
    dp.push(DiagramPoint {
        x: f64::NAN,
        y: f64::INFINITY,
    });
    emit(&mut o, "DiagramPoint", &dp);
    emit::<RefrigerationResult>(&mut o, "RefrigerationResult", &[refr]);
    emit::<CycleComparisonEntry>(
        &mut o,
        "CycleComparisonEntry",
        &cycle::compare_cycles(&[
            (CycleKind::Otto, &otto, 2000.0, 300.0),
            (CycleKind::Diesel, &diesel, 2000.0, 300.0),
            (CycleKind::Brayton, &brayton, 1400.0, 300.0),
            (CycleKind::Rankine, &rankine, 573.15, 318.95),
        ]),
    );
    emit(
        &mut o,
        "ThermalGridVisualization",
        &[ThermalGridVisualization::from_grid_2d(&g2)],
    );
    emit(
        &mut o,
        "TemperatureProfile",
        &[TemperatureProfile::from_grid_1d(
            &g1,
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        )],
    );
    let cdd = CycleDiagramData::from_cycle_result(&otto, 4);
    emit(&mut o, "CycleDiagramData", std::slice::from_ref(&cdd));
    emit::<CycleStatePoint>(&mut o, "CycleStatePoint", &cdd.state_points);
    emit(
        &mut o,
        "ThermalNetworkVisualization",
        &[ThermalNetworkVisualization::from_network(&net, &temps)],
    );
    emit(
        &mut o,
        "HeatFluxField",
        &[HeatFluxField::from_grid_2d(&g2, 50.0)],
    );
    emit::<ThermalMaterial>(
        &mut o,
        "ThermalMaterial",
        &ALL_MATERIALS
            .iter()
            .map(|m| (*m).clone())
            .collect::<Vec<_>>(),
    );
    emit(
        &mut o,
        "BoundaryCondition",
        &[
            BoundaryCondition::Fixed(300.0),
            BoundaryCondition::Insulated,
            BoundaryCondition::Convective {
                h: 10.0,
                t_fluid: 293.15,
                k: 50.0,
            },
        ],
    );
    // The variants of the other two payload forms (base documents 2 and 3 get only their base line above).
    for (kind, doc) in [
        ("e_obj_null", &b"{\"Insulated\":null}"[..]),
        ("e_fixed_bare", b"\"Fixed\""),
        ("e_fixed_int", b"{\"Fixed\":300}"),
        ("e_fixed_seq", b"{\"Fixed\":[300.0]}"),
        ("e_conv_seq", b"{\"Convective\":[10.0,293.15,50.0]}"),
        (
            "e_conv_missing",
            b"{\"Convective\":{\"h\":10.0,\"t_fluid\":293.15}}",
        ),
        (
            "e_conv_unknown",
            b"{\"Convective\":{\"h\":10.0,\"t_fluid\":293.15,\"k\":50.0,\"zz\":1}}",
        ),
        (
            "e_conv_dup",
            b"{\"Convective\":{\"h\":10.0,\"h\":11.0,\"t_fluid\":293.15,\"k\":50.0}}",
        ),
        ("e_two_tags", b"{\"Fixed\":300.0,\"Insulated\":null}"),
        ("e_unknown", b"{\"Radiative\":1.0}"),
    ] {
        line::<BoundaryCondition>(&mut o, "BoundaryCondition", kind, doc, "{\"Fixed\":300.0}");
    }
    emit(
        &mut o,
        "Side",
        &[Side::Left, Side::Right, Side::Top, Side::Bottom],
    );
    emit(&mut o, "ThermalGrid1D", &[g1]);
    emit(&mut o, "ThermalGrid2D", &[g2]);
    emit(&mut o, "ThermalNetwork", &[net]);
    emit(
        &mut o,
        "Phase",
        &[
            Phase::Solid,
            Phase::Liquid,
            Phase::Gas,
            Phase::SupercriticalFluid,
        ],
    );
    emit::<SubstanceData>(
        &mut o,
        "SubstanceData",
        &ALL_SUBSTANCES
            .iter()
            .map(|s| (*s).clone())
            .collect::<Vec<_>>(),
    );
    emit::<GasData>(
        &mut o,
        "GasData",
        &ALL_GASES.iter().map(|g| (*g).clone()).collect::<Vec<_>>(),
    );
    emit(&mut o, "SaturatedEntry", &sat);
    emit(&mut o, "WetSteamProperties", &wet);
    emit(&mut o, "SuperheatedEntry", &sup);

    if let Err(e) = std::fs::create_dir_all(&dir) {
        eprintln!("probe_json: {dir}: {e}");
        return ExitCode::FAILURE;
    }
    let f = floats(&o.floats);
    let types: std::collections::BTreeSet<&str> = o
        .cases
        .lines()
        .filter_map(|l| l.split('\t').next())
        .collect();
    for (name, body) in [("cases.txt", &o.cases), ("f64.txt", &f)] {
        if let Err(e) = std::fs::write(format!("{dir}/{name}"), body) {
            eprintln!("probe_json: {dir}/{name}: {e}");
            return ExitCode::FAILURE;
        }
    }
    println!(
        "probe_json: {} documents over {} types, {} doubles -> {dir}",
        o.count,
        types.len(),
        f.lines().count()
    );
    ExitCode::SUCCESS
}
