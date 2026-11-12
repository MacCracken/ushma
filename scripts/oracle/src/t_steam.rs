//! `src/steam_data.cyr` — the private IAPWS-IF97 steam tables and constants (port plan §2.1).
//!
//! `steam.rs`'s tables are private: the API exposes interpolated values, not the rows. So every
//! cell comes from its literal in `rust-old/src/steam.rs`, parsed with Rust's `f64::from_str`
//! (`rustsrc.rs`), never re-typed. The rows are proven by the steam module's golden rows (M3),
//! and `probe_tables.py` re-parses every literal independently in Python and checks each
//! compiled cell. What the Rust code relies on is asserted here: every entry carries every field
//! in declaration order; `SUPERHEATED_DATA.len() == SUPERHEAT_PRESSURES.len() * TEMPS_PER_TIER`
//! (the source's own `const _` assertion); and row k of `SUPERHEATED_DATA` is at pressure
//! `SUPERHEAT_PRESSURES[k / TEMPS_PER_TIER]`.

use crate::cyr::{comment_line, f64_column, f64_lit, header};
use crate::rustsrc::{ConstItem, array_items, const_items, eval_f64, eval_int, struct_entries};
use std::fmt::Write as _;

/// The steam items this file generates, and the one it excuses (the compile-time assertion it
/// re-checks).
pub const GENERATED: &[&str] = &[
    "SATURATED_TABLE",
    "R4_N",
    "R4_T_MIN",
    "T_CRIT",
    "P_CRIT",
    "QUALITY_ROUNDOFF",
    "SUPERHEAT_PRESSURES",
    "TEMPS_PER_TIER",
    "SUPERHEATED_DATA",
];
pub const EXCUSED: &[(&str, &str)] = &[(
    "_",
    "the compile-time assertion SUPERHEATED_DATA.len() == SUPERHEAT_PRESSURES.len() * \
     TEMPS_PER_TIER, re-checked by t_steam",
)];

const SAT_FIELDS: [&str; 9] = [
    "temperature",
    "pressure",
    "v_f",
    "v_g",
    "h_f",
    "h_fg",
    "h_g",
    "s_f",
    "s_g",
];
const SH_FIELDS: [&str; 5] = [
    "temperature",
    "pressure",
    "specific_volume",
    "specific_enthalpy",
    "specific_entropy",
];

/// The `//` / `///` comment block directly above the line that starts with `item`, as `#` lines.
fn comments_before(src: &str, item: &str) -> String {
    let lines: Vec<&str> = src.lines().collect();
    let Some(at) = lines.iter().position(|l| l.trim_start().starts_with(item)) else {
        return String::new();
    };
    let mut block = Vec::new();
    let mut k = at;
    while k > 0 {
        k -= 1;
        let l = lines[k].trim_start();
        let text = l
            .strip_prefix("///")
            .or_else(|| l.strip_prefix("//"))
            .map(|r| r.strip_prefix(' ').unwrap_or(r).trim_end());
        match text {
            Some(t) => block.push(t.to_string()),
            None if l.starts_with("#[") => {}
            None => break,
        }
    }
    block.reverse();
    let mut out = String::new();
    for b in block {
        let _ = writeln!(out, "{}", comment_line(&b));
    }
    out
}

fn item<'a>(items: &'a [ConstItem], name: &str) -> Result<&'a ConstItem, String> {
    items
        .iter()
        .find(|c| c.name == name)
        .ok_or_else(|| format!("steam.rs: no const {name}"))
}

/// A parsed table: each entry's values, and their literal text (for row labels).
type Rows = Vec<(Vec<f64>, Vec<String>)>;

/// Parse a table of struct literals and check each entry's fields against `fields`.
fn table(c: &ConstItem, ty: &str, fields: &[&str]) -> Result<Rows, String> {
    let mut out = Vec::new();
    for (k, entry) in struct_entries(&c.expr, ty)?.into_iter().enumerate() {
        let names: Vec<&str> = entry.iter().map(|(f, _)| f.as_str()).collect();
        if names != fields {
            return Err(format!(
                "{}[{k}]: fields {names:?}, expected {fields:?}",
                c.name
            ));
        }
        let mut vals = Vec::new();
        let mut lits = Vec::new();
        for (_, v) in &entry {
            vals.push(eval_f64(v).map_err(|e| format!("{}[{k}]: {e}", c.name))?);
            lits.push(v.clone());
        }
        out.push((vals, lits));
    }
    Ok(out)
}

/// The literal text of `v` without digit separators, for a row label.
fn label_lit(v: &str) -> String {
    v.chars().filter(|&ch| ch != '_').collect()
}

/// The saturated table's rows (each entry's nine values, in `SaturatedEntry` order) from `rust-old/src/steam.rs`'s
/// source, for the golden rows that run every row through the lookups (`m_steam.rs`).
pub fn saturated_rows(src: &str) -> Result<Vec<Vec<f64>>, String> {
    let items = const_items("steam", src)?;
    let fields = [
        "temperature",
        "pressure",
        "v_f",
        "v_g",
        "h_f",
        "h_fg",
        "h_g",
        "s_f",
        "s_g",
    ];
    Ok(
        table(item(&items, "SATURATED_TABLE")?, "SaturatedEntry", &fields)?
            .into_iter()
            .map(|(vals, _)| vals)
            .collect(),
    )
}

pub fn emit(src: &str, provenance: &str) -> Result<String, String> {
    let items = const_items("steam", src)?;
    let mut out = header(
        "steam data: the IAPWS-IF97 saturated and superheated tables and the Region 4 constants",
        provenance,
        &[
            "Private to the steam module (`_US_STEAM_`): ushma 1.4.0's tables are private, so every\n\
             cell is its literal in rust-old/src/steam.rs, parsed with Rust's f64::from_str and written\n\
             as a round-trip decimal. scripts/oracle/probe_tables.py re-parses each literal\n\
             independently and checks every compiled cell; the steam golden rows (M3) prove them\n\
             through the API.",
        ],
    );

    // Region 4.
    let r4 = item(&items, "R4_N")?;
    let r4_vals: Vec<f64> = array_items(&r4.expr)?
        .iter()
        .map(|v| eval_f64(v))
        .collect::<Result<_, _>>()?;
    if r4_vals.len() != 10 {
        return Err(format!("R4_N: {} coefficients, expected 10", r4_vals.len()));
    }
    out.push_str(&comments_before(src, "const R4_N:"));
    out.push_str(
        "# Rust: ushma::steam::R4_N (private), as ten scalar consts named by IF97's n1 … n10: a const is an immediate,\n\
         # about 13 ns a Region 4 evaluation cheaper than ten indexed loads (M3b's review).\n",
    );
    for (k, v) in r4_vals.iter().enumerate() {
        let _ = writeln!(out, "const _US_STEAM_R4_N{} = {};", k + 1, f64_lit(*v)?);
    }
    out.push('\n');
    for name in ["R4_T_MIN", "T_CRIT", "P_CRIT", "QUALITY_ROUNDOFF"] {
        let c = item(&items, name)?;
        if c.ty != "f64" {
            return Err(format!("steam::{name}: type {}", c.ty));
        }
        out.push_str(&comments_before(src, &format!("const {name}:")));
        let _ = writeln!(out, "# Rust: ushma::steam::{name} (private)");
        let _ = writeln!(
            out,
            "const _US_STEAM_{name} = {};\n",
            f64_lit(eval_f64(&c.expr)?)?
        );
    }

    // The saturated table.
    let sat_item = item(&items, "SATURATED_TABLE")?;
    let sat = table(sat_item, "SaturatedEntry", &SAT_FIELDS)?;
    out.push_str(&comments_before(src, "const SATURATED_TABLE:"));
    out.push_str("# Rust: ushma::steam::SATURATED_TABLE (private; its length)\n");
    let _ = writeln!(out, "const _US_STEAM_SATURATED_COUNT = {};\n", sat.len());
    let sat_labels: Vec<String> = sat
        .iter()
        .map(|(_, lits)| format!("T = {} K", label_lit(&lits[0])))
        .collect();
    for (j, field) in SAT_FIELDS.iter().enumerate() {
        out.push_str(&crate::cyr::doc_comment_after(
            src,
            "pub struct SaturatedEntry {",
            &format!("pub {field}:"),
        )?);
        let _ = writeln!(
            out,
            "# Indexed by SATURATED_TABLE row; mirrors `SaturatedEntry::{field}`."
        );
        let col: Vec<f64> = sat.iter().map(|(v, _)| v[j]).collect();
        out.push_str(&f64_column(
            &format!("_US_STEAM_SAT_{}", field.to_uppercase()),
            &col,
            &sat_labels,
        )?);
        out.push('\n');
    }

    // The superheated table.
    let tiers_item = item(&items, "SUPERHEAT_PRESSURES")?;
    let tiers: Vec<f64> = array_items(&tiers_item.expr)?
        .iter()
        .map(|v| eval_f64(v))
        .collect::<Result<_, _>>()?;
    let per = eval_int(&item(&items, "TEMPS_PER_TIER")?.expr)?;
    let sh_item = item(&items, "SUPERHEATED_DATA")?;
    let sh = table(sh_item, "SuperheatedEntry", &SH_FIELDS)?;
    let per_u = usize::try_from(per).map_err(|e| format!("TEMPS_PER_TIER: {e}"))?;
    if sh.len() != tiers.len() * per_u {
        return Err(format!(
            "SUPERHEATED_DATA: {} rows, but {} tiers × {per} per tier",
            sh.len(),
            tiers.len()
        ));
    }
    for (k, (v, _)) in sh.iter().enumerate() {
        let tier = tiers[k / per_u];
        if v[1].to_bits() != tier.to_bits() {
            return Err(format!(
                "SUPERHEATED_DATA[{k}]: pressure {} is not its tier's {tier}",
                v[1]
            ));
        }
    }
    out.push_str(&comments_before(src, "const SUPERHEAT_PRESSURES:"));
    out.push_str("# Rust: ushma::steam::SUPERHEAT_PRESSURES (private)\n");
    let _ = writeln!(
        out,
        "const _US_STEAM_SUPERHEAT_TIER_COUNT = {};",
        tiers.len()
    );
    let tier_labels: Vec<String> = (0..tiers.len()).map(|k| format!("tier {k}")).collect();
    out.push_str(&f64_column(
        "_US_STEAM_SUPERHEAT_PRESSURES",
        &tiers,
        &tier_labels,
    )?);
    out.push('\n');
    out.push_str(&comments_before(src, "const TEMPS_PER_TIER:"));
    out.push_str("# Rust: ushma::steam::TEMPS_PER_TIER (private)\n");
    let _ = writeln!(out, "const _US_STEAM_TEMPS_PER_TIER = {per};\n");
    out.push_str(&comments_before(src, "const SUPERHEATED_DATA:"));
    out.push_str("# Rust: ushma::steam::SUPERHEATED_DATA (private; its length)\n");
    let _ = writeln!(out, "const _US_STEAM_SUPERHEATED_COUNT = {};\n", sh.len());
    let sh_labels: Vec<String> = sh
        .iter()
        .map(|(_, lits)| format!("{} Pa, {} K", label_lit(&lits[1]), label_lit(&lits[0])))
        .collect();
    for (j, field) in SH_FIELDS.iter().enumerate() {
        out.push_str(&crate::cyr::doc_comment_after(
            src,
            "pub struct SuperheatedEntry {",
            &format!("pub {field}:"),
        )?);
        let _ = writeln!(
            out,
            "# Indexed by SUPERHEATED_DATA row; mirrors `SuperheatedEntry::{field}`."
        );
        let col: Vec<f64> = sh.iter().map(|(v, _)| v[j]).collect();
        out.push_str(&f64_column(
            &format!("_US_STEAM_SH_{}", field.to_uppercase()),
            &col,
            &sh_labels,
        )?);
        out.push('\n');
    }
    Ok(out.trim_end().to_string() + "\n")
}
