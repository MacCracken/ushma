//! `src/constants.cyr` — ushma's public scalar constants, as `US_` consts.
//!
//! ushma 1.4.0 has 9 public scalar consts over 6 modules (8 `f64`, 1 `usize`). `ALL` lists every
//! one (`complete` proves it against `rust-old/src`), and `tests/golden/constants.tsv` pins every
//! bit (`m_constants.rs`). They all go into `src/constants.cyr`, which comes second in `[lib]`
//! because the modules reach across for them (`GAS_CONSTANT` from state is used by entropy,
//! phase, cycle, chem, stat and bridge).
//!
//! The struct-valued consts (`material::COPPER`, `phase::WATER_PHASE`, `state::GAS_NITROGEN`,
//! `chem::CO`, …) and the private constants (steam's IF97 coefficients, bridge's Kleiber factor,
//! …) are table data: they are generated with their module's data at its milestone (M0 part 2b).
//!
//! Values come from the oracle's API; the comments are the Rust doc comments, which carry the
//! citations (CODATA 2018, SI 2019).

use crate::cyr::{doc_comment, f64_lit, header};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

/// A constant's value: an `f64`, or an integer (`usize` in Rust).
#[derive(Clone, Copy)]
pub enum Value {
    F(f64),
    I(i64),
}

/// One public scalar const of the oracle.
pub struct Constant {
    pub module: &'static str,
    pub name: &'static str,
    pub value: Value,
}

impl Constant {
    /// `state::GAS_CONSTANT`, the golden row's input.
    pub fn path(&self) -> String {
        format!("{}::{}", self.module, self.name)
    }
}

/// `MAX_ADAPTIVE_SUBSTEPS` as an i64 (a small usize; checked, not cast).
fn max_substeps() -> i64 {
    i64::try_from(ushma::numerical::MAX_ADAPTIVE_SUBSTEPS).unwrap_or(-1)
}

/// Every public scalar const of ushma 1.4.0: modules in file-name order, then source order.
fn all_constants() -> Vec<Constant> {
    vec![
        Constant {
            module: "chem",
            name: "STANDARD_T",
            value: Value::F(ushma::chem::STANDARD_T),
        },
        Constant {
            module: "numerical",
            name: "MAX_ADAPTIVE_SUBSTEPS",
            value: Value::I(max_substeps()),
        },
        Constant {
            module: "stat",
            name: "PLANCK",
            value: Value::F(ushma::stat::PLANCK),
        },
        Constant {
            module: "stat",
            name: "AVOGADRO",
            value: Value::F(ushma::stat::AVOGADRO),
        },
        Constant {
            module: "state",
            name: "GAS_CONSTANT",
            value: Value::F(ushma::state::GAS_CONSTANT),
        },
        Constant {
            module: "state",
            name: "ATM",
            value: Value::F(ushma::state::ATM),
        },
        Constant {
            module: "state",
            name: "STANDARD_TEMP",
            value: Value::F(ushma::state::STANDARD_TEMP),
        },
        Constant {
            module: "transfer",
            name: "STEFAN_BOLTZMANN",
            value: Value::F(ushma::transfer::STEFAN_BOLTZMANN),
        },
        Constant {
            module: "transfer",
            name: "BOLTZMANN_K",
            value: Value::F(ushma::transfer::BOLTZMANN_K),
        },
    ]
}

/// Every constant, for the golden emitter too.
pub fn all() -> Vec<Constant> {
    all_constants()
}

/// The scalar types a public const can have in the oracle; anything else is table data.
const SCALARS: [&str; 9] = [
    "f64", "f32", "usize", "u64", "u32", "i64", "i32", "u8", "bool",
];

/// Proves `all()` is every top-level `pub const NAME: <scalar>` item in `rust-old/src` (no more,
/// no fewer, `integration/` included), and that the Cyrius names are unique.
pub fn complete(root: &Path) -> Result<(), String> {
    let mut found = BTreeSet::new();
    let mut dirs = vec![root.join("rust-old/src")];
    while let Some(dir) = dirs.pop() {
        let entries = std::fs::read_dir(&dir).map_err(|e| format!("{}: {e}", dir.display()))?;
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                dirs.push(path);
                continue;
            }
            if path.extension().is_none_or(|x| x != "rs") {
                continue;
            }
            let module = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or_default()
                .to_string();
            let text =
                std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
            for line in text.lines() {
                if let Some(rest) = line.strip_prefix("pub const ")
                    && let Some((name, ty)) = rest.split_once(':')
                {
                    let ty = ty.trim_start();
                    if SCALARS.iter().any(|s| {
                        ty.starts_with(&format!("{s} ")) || ty.starts_with(&format!("{s}="))
                    }) {
                        found.insert(format!("{module}::{name}"));
                    }
                }
            }
        }
    }
    let consts = all_constants();
    let listed: BTreeSet<String> = consts.iter().map(Constant::path).collect();
    if listed.len() != consts.len() {
        return Err("t_constants: a constant is listed twice".to_string());
    }
    let missing: Vec<_> = found.difference(&listed).collect();
    let extra: Vec<_> = listed.difference(&found).collect();
    if !missing.is_empty() || !extra.is_empty() {
        return Err(format!(
            "t_constants: not listed {missing:?}; listed but not in rust-old/src {extra:?}"
        ));
    }
    let mut names = BTreeMap::new();
    for c in &consts {
        if let Some(prev) = names.insert(c.name, c.module) {
            return Err(format!(
                "US_{} would name both {prev}::{} and {}",
                c.name,
                c.name,
                c.path()
            ));
        }
    }
    Ok(())
}

/// `src/constants.cyr`: each constant as `const US_<NAME> = <literal>;` under its Rust doc comment
/// and a `# Rust: ushma::<module>::<NAME>` line (which `validate.py` and `probe_tables.py` read to
/// tie it to its golden row).
pub fn emit(root: &Path, provenance: &str) -> Result<String, String> {
    complete(root)?;
    let mut out = header(
        "public constants",
        provenance,
        &[
            "Each f64 is the shortest decimal that round-trips to the oracle's value, written out\n\
             positionally: the Cyrius lexer has no exponent notation, and decimal literals are\n\
             correctly rounded (6.6.10+). tests/golden/constants.tsv pins every bit.",
            "Struct-valued and private constants are table data, generated with their module's\n\
             data file (port plan §4 step 6).",
        ],
    );
    let mut sources = BTreeMap::new();
    for c in all_constants() {
        if !sources.contains_key(c.module) {
            let p = root.join(format!("rust-old/src/{}.rs", c.module));
            let text = std::fs::read_to_string(&p).map_err(|e| format!("{}: {e}", p.display()))?;
            sources.insert(c.module, text);
        }
        let src = &sources[c.module];
        out.push_str(&doc_comment(src, &format!("pub const {}:", c.name))?);
        out.push_str(&format!("# Rust: ushma::{}\n", c.path()));
        let lit = match c.value {
            Value::F(x) => f64_lit(x)?,
            Value::I(n) => n.to_string(),
        };
        out.push_str(&format!("const US_{} = {lit};\n\n", c.name));
    }
    Ok(out.trim_end().to_string() + "\n")
}
