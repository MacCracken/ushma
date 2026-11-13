//! `oracle` — writes everything the Cyrius port derives from the frozen Rust ushma 1.4.0 in
//! `rust-old/` (read-only; it is never edited):
//! - the golden corpus, `tests/golden/<module>.tsv` (port plan §4 step 5);
//! - the generated Cyrius sources, `src/constants.cyr` and `src/<module>_data.cyr`
//!   (port plan §4 step 6).
//!
//! ```sh
//! cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml            # write
//! cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check # verify
//! ```
//!
//! `--check` regenerates every output in memory and compares it byte for byte with the
//! committed file. It exits 1 if any output differs or is missing, if a `.tsv` in
//! `tests/golden/` has no emitter, or if the corpus was recorded on another glibc platform
//! (`glibc.rs`): the transcendental rows are glibc's bits, so they can only be checked where
//! they were made. A `<module>.cr.tsv` beside a module's file is the corrected oracle (ADR 0004),
//! written and checked by `validate.py`; here it is only stale when its module has no emitter.
//!
//! Adapted from kimiya's `scripts/oracle/` (same format, same rules; tanmatra's originally).
//! Modules are added at their milestone, with the golden grid from the port plan §6: `constants`
//! and `math` (M0 part 2a); the table data and the rows that pin it (part 2b: `chem`, `material`,
//! `phase`, `state`, plus the private `steam_data` and `constants_private`); each module's own rows
//! from M1 on.

mod cyr;
mod glibc;
mod grids;
mod json_cases;
mod m_bridge;
mod m_chem;
mod m_constants;
mod m_cycle;
mod m_entropy;
mod m_fmtf;
mod m_material;
mod m_math;
mod m_numerical;
mod m_phase;
mod m_stat;
mod m_state;
mod m_steam;
mod m_transfer;
mod rustsrc;
mod sweep;
mod t_constants;
mod t_private;
mod t_steam;
mod t_tables;
mod tsv;

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

/// Every golden emitter, sorted by file name: the modules emitted whole (`bridge`, `constants`, `cycle`, `entropy`, `fmtf`, `math`,
/// `numerical`, `stat`, `steam`, `transfer`) and the four the public tables seed (`chem`, `material`, `phase`, `state`; M0 part 2b), with their
/// function and JSON rows added (`chem`, `material`, `phase`, `state`).
fn golden_modules() -> Result<Vec<tsv::Module>, String> {
    let mut v = vec![
        m_bridge::emit(),
        m_constants::emit(),
        m_cycle::emit(),
        m_entropy::emit(),
        m_fmtf::emit(),
        m_math::emit(),
        m_numerical::emit(),
        m_stat::emit(),
        m_steam::emit(),
        m_transfer::emit(),
    ];
    let mut tables = t_tables::golden()?;
    for m in &mut tables {
        if m.name == "material" {
            m_material::rows(m);
        }
        if m.name == "state" {
            m_state::rows(m);
        }
        if m.name == "phase" {
            m_phase::rows(m);
        }
        if m.name == "chem" {
            m_chem::rows(m);
        }
    }
    v.extend(tables);
    v.sort_by_key(|m| m.name);
    Ok(v)
}

/// One generated file: its path from the repo root, its contents, and a summary.
struct Output {
    rel: String,
    text: String,
    summary: String,
}

fn read(path: &Path) -> Result<String, String> {
    fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))
}

/// Every generated Cyrius source, in file-name order. Each is built from the oracle's API; the
/// Rust source is read for comments (and, for private data, literals with a completeness proof).
fn cyrius_outputs(root: &Path, prov: &str) -> Result<Vec<Output>, String> {
    let steam_src = read(&root.join("rust-old/src/steam.rs"))?;
    let mut files = vec![
        (
            "src/constants.cyr".to_string(),
            t_constants::emit(root, prov)?,
        ),
        (
            "src/constants_private.cyr".to_string(),
            t_private::emit(root, prov)?,
        ),
        (
            "src/steam_data.cyr".to_string(),
            t_steam::emit(&steam_src, prov)?,
        ),
    ];
    for t in t_tables::all() {
        files.push((format!("src/{}.cyr", t.file), t.emit_cyr(root, prov)?));
    }
    files.sort();
    Ok(files
        .into_iter()
        .map(|(rel, text)| {
            let lines = text.lines().count();
            Output {
                rel,
                text,
                summary: format!("{lines} lines"),
            }
        })
        .collect())
}

/// `version = "x"` from the first line that starts with it (the `[package]` one).
fn toml_version(text: &str) -> Option<&str> {
    let line = text.lines().find(|l| l.starts_with("version = \""))?;
    line.strip_prefix("version = \"")?.strip_suffix('"')
}

/// The locked version of `crate_name` in a Cargo.lock.
fn lock_version<'a>(lock: &'a str, crate_name: &str) -> Option<&'a str> {
    let needle = format!("name = \"{crate_name}\"\n");
    let at = lock.find(&needle)? + needle.len();
    let line = lock[at..].lines().next()?;
    line.strip_prefix("version = \"")?.strip_suffix('"')
}

/// `ushma 1.4.0 (rust-old/), hisab 1.4.0`: the oracle and its numerical dependency.
fn oracle_provenance(root: &Path) -> Result<String, String> {
    let cargo = read(&root.join("rust-old/Cargo.toml"))?;
    let lock = read(&root.join("scripts/oracle/Cargo.lock"))?;
    let us = toml_version(&cargo).ok_or("rust-old/Cargo.toml: no version")?;
    let hs = lock_version(&lock, "hisab").ok_or("scripts/oracle/Cargo.lock: no hisab")?;
    Ok(format!("ushma {us} (rust-old/), hisab {hs}"))
}

/// Write `text` to `root/rel`, or in check mode compare it. Returns false on a mismatch.
fn emit(root: &Path, check: bool, out: &Output) -> Result<bool, String> {
    let path: PathBuf = root.join(&out.rel);
    if !check {
        fs::write(&path, &out.text).map_err(|e| format!("{}: {e}", path.display()))?;
        println!("wrote    {} ({})", out.rel, out.summary);
        return Ok(true);
    }
    match fs::read_to_string(&path) {
        Ok(old) if old == out.text => {
            println!("ok       {} ({})", out.rel, out.summary);
            Ok(true)
        }
        Ok(_) => {
            println!("DIFFERS  {}", out.rel);
            Ok(false)
        }
        Err(_) => {
            println!("MISSING  {}", out.rel);
            Ok(false)
        }
    }
}

/// The libm platform in a TSV's `# oracle:` header line: the text from `libm: ` to the next
/// `. `, e.g. `glibc 2.44 x86_64 fma`.
fn recorded_platform(text: &str) -> Option<String> {
    let line = text
        .lines()
        .take(5)
        .find_map(|l| l.strip_prefix("# oracle: "))?;
    let (_, rest) = line.split_once("libm: ")?;
    Some(rest.split_once(". ").map_or(rest, |(p, _)| p).to_string())
}

/// In check mode: every committed `.tsv` must have been recorded on the running platform.
fn check_platform(dir: &Path, platform: &str) -> bool {
    let mut ok = true;
    let Ok(entries) = fs::read_dir(dir) else {
        return true; // every file then reports MISSING
    };
    let mut names: Vec<String> = entries
        .flatten()
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| n.ends_with(".tsv"))
        .collect();
    names.sort();
    for name in names {
        let Ok(text) = fs::read_to_string(dir.join(&name)) else {
            continue;
        };
        match recorded_platform(&text) {
            Some(p) if p == platform => {}
            Some(p) => {
                println!(
                    "PLATFORM tests/golden/{name}: recorded on {p}, running {platform}. Its \
                     transcendental rows are glibc's bits: check it on the recorded platform, \
                     or regenerate deliberately"
                );
                ok = false;
            }
            None => {
                println!("PLATFORM tests/golden/{name}: no recorded libm platform in its header");
                ok = false;
            }
        }
    }
    ok
}

fn run() -> Result<bool, String> {
    let check = std::env::args().skip(1).any(|a| a == "--check");
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let root = root
        .canonicalize()
        .map_err(|e| format!("{}: {e}", root.display()))?;
    let oracle = oracle_provenance(&root)?;
    let platform = glibc::platform();
    // Every TSV row may hold glibc's bits; the generated .cyr files hold only the oracle's
    // constants, so they name the oracle alone and don't change with the platform.
    let tsv_prov = format!("{oracle}; libm: {platform}");
    let dir = root.join("tests/golden");
    if !check {
        fs::create_dir_all(&dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    }

    let mut ok = true;
    if check {
        ok &= check_platform(&dir, &platform);
    }

    let mut outputs = Vec::new();
    let mut emitted = BTreeSet::new();
    for module in golden_modules()? {
        let name = format!("{}.tsv", module.name);
        emitted.insert(name.clone());
        outputs.push(Output {
            rel: format!("tests/golden/{name}"),
            text: module.render(&tsv_prov),
            summary: format!("{} rows", module.len()),
        });
    }
    // The harness reads this list, so a module added here can't go untested: without a
    // dispatcher it is reported as pending, row by row.
    let mut index = String::from(
        "# ushma golden corpus: one module per line; tests/parity.tcyr replays each\n",
    );
    for name in &emitted {
        index.push_str(name.trim_end_matches(".tsv"));
        index.push('\n');
    }
    outputs.push(Output {
        rel: "tests/golden/modules.txt".to_string(),
        summary: format!("{} modules", emitted.len()),
        text: index,
    });
    outputs.extend(cyrius_outputs(&root, &oracle)?);

    for out in &outputs {
        ok &= emit(&root, check, out)?;
    }

    // A .tsv with no emitter is stale: the harness would test rows nobody regenerates. A
    // `<module>.cr.tsv` (ADR 0004, written by validate.py) belongs to its module's emitter.
    if let Ok(entries) = fs::read_dir(&dir) {
        for entry in entries.flatten() {
            let name = entry.file_name().to_string_lossy().into_owned();
            let module_file = name
                .strip_suffix(".cr.tsv")
                .map_or_else(|| name.clone(), |m| format!("{m}.tsv"));
            if name.ends_with(".tsv") && !emitted.contains(&module_file) {
                println!("STALE    tests/golden/{name} (no emitter)");
                ok = false;
            }
        }
    }
    Ok(ok)
}

fn main() -> ExitCode {
    match run() {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::from(1),
        Err(e) => {
            eprintln!("oracle: {e}");
            ExitCode::from(2)
        }
    }
}
