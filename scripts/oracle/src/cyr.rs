//! Cyrius source emission for the generated modules (port plan §4 step 6). Adapted from
//! kimiya's `scripts/oracle/src/cyr.rs` (itself from tanmatra's).
//!
//! Rules the output must keep. Each was probed on the 6.7.18 toolchain by the sibling ports,
//! and is re-checked here by `probe_tables.py` (M0 part 2b):
//! - **f64 values are decimal literals, written out positionally.** The lexer has no
//!   exponent notation. Decimals have been correctly rounded since 6.6.10, so the shortest
//!   round-trip digits give back the oracle's exact bits.
//! - **Hex is never used in a `const`**, where it is an integer and `CONST * x` then does
//!   integer arithmetic.
//! - **±∞ are `(1.0 / 0.0)` and `-(1.0 / 0.0)`.** NaN can't be a compile-time f64 (the
//!   compiler refuses it: its bits depend on the target), so a NaN is an error here.
//! - **Every public global gets the `US_` prefix** (ADR 0001, P6). A bundle is concatenated
//!   into its consumers, so names share one namespace; enum variants are global too.

use std::fmt::Write as _;

/// The Cyrius literal for `x`, checked to parse back to the same bits.
pub fn f64_lit(x: f64) -> Result<String, String> {
    if x.is_nan() {
        return Err("NaN cannot be a compile-time f64 in Cyrius".to_string());
    }
    if x.is_infinite() {
        return Ok(if x > 0.0 {
            "(1.0 / 0.0)"
        } else {
            "-(1.0 / 0.0)"
        }
        .to_string());
    }
    let lit = positional(x);
    if lit.contains(['e', 'E', 'x', 'X']) {
        return Err(format!("{lit}: not a positional decimal"));
    }
    let back: f64 = lit.parse().map_err(|e| format!("{lit}: {e}"))?;
    if back.to_bits() != x.to_bits() {
        return Err(format!("{lit} does not round-trip to {:016x}", x.to_bits()));
    }
    Ok(lit)
}

/// The shortest round-trip digits of `x` (Rust's `{:e}`), written without an exponent.
fn positional(x: f64) -> String {
    let sci = format!("{:e}", x.abs());
    let (mant, exp) = sci.split_once('e').unwrap_or((&sci, "0"));
    let exp: i64 = exp.parse().unwrap_or(0);
    let digits: String = mant.chars().filter(char::is_ascii_digit).collect();
    let point = exp + 1; // digits before the decimal point
    let n = i64::try_from(digits.len()).unwrap_or(i64::MAX);
    let body = if point <= 0 {
        let zeros = "0".repeat(usize::try_from(-point).unwrap_or(0));
        format!("0.{zeros}{digits}")
    } else if point >= n {
        let zeros = "0".repeat(usize::try_from(point - n).unwrap_or(0));
        format!("{digits}{zeros}.0")
    } else {
        let at = usize::try_from(point).unwrap_or(0);
        format!("{}.{}", &digits[..at], &digits[at..])
    };
    if x.is_sign_negative() {
        format!("-{body}")
    } else {
        body
    }
}

/// The `///` doc block directly above the first line of `src` that starts with `item`
/// (attributes in between are skipped), as `# ` comment lines. Empty if there is none.
/// Only comments are read from the Rust source; values always come from the oracle's API.
pub fn doc_comment(src: &str, item: &str) -> Result<String, String> {
    doc_comment_after(src, "", item)
}

/// As [`doc_comment`], for the first line starting with `item` (after leading whitespace) that
/// follows the first line starting with `anchor`: a field inside `pub struct GasData {`.
pub fn doc_comment_after(src: &str, anchor: &str, item: &str) -> Result<String, String> {
    let lines: Vec<&str> = src.lines().collect();
    let from = lines
        .iter()
        .position(|l| l.trim_start().starts_with(anchor))
        .ok_or_else(|| format!("no `{anchor}` in the Rust source"))?;
    let at = from
        + lines[from..]
            .iter()
            .position(|l| l.trim_start().starts_with(item))
            .ok_or_else(|| format!("no `{item}` after `{anchor}` in the Rust source"))?;
    let mut docs = Vec::new();
    let mut k = at;
    while k > 0 {
        k -= 1;
        let l = lines[k].trim_start();
        if l.starts_with("#[") {
            continue;
        }
        match l.strip_prefix("///") {
            Some(rest) => docs.push(rest.strip_prefix(' ').unwrap_or(rest).to_string()),
            None => break,
        }
    }
    docs.reverse();
    let mut out = String::new();
    for d in docs {
        let _ = writeln!(out, "{}", comment_line(d.trim_end()));
    }
    Ok(out)
}

/// `var <name>: f64[n] = { … };`, one value per line with its row label as a comment.
pub fn f64_column(name: &str, values: &[f64], labels: &[String]) -> Result<String, String> {
    if values.len() != labels.len() {
        return Err(format!(
            "{name}: {} values, {} labels",
            values.len(),
            labels.len()
        ));
    }
    let mut out = format!("var {name}: f64[{}] = {{\n", values.len());
    for (v, l) in values.iter().zip(labels) {
        let _ = writeln!(out, "    {},  # {l}", f64_lit(*v)?);
    }
    out.push_str("};\n");
    Ok(out)
}

/// A NUL-terminated string blob and its offsets: `var <prefix>_NAME_OFF: u16[n]` (row i's offset)
/// and `var <prefix>_NAME_BYTES: u8[m]` (every string's UTF-8 bytes and a NUL, in row order).
pub fn name_blob(prefix: &str, names: &[String], what: &str) -> Result<String, String> {
    let mut off = Vec::new();
    let mut bytes: Vec<(String, Vec<u8>)> = Vec::new();
    let mut at = 0_usize;
    for n in names {
        if n.as_bytes().contains(&0) {
            return Err(format!("{prefix}: `{n}` holds a NUL"));
        }
        off.push(u16::try_from(at).map_err(|_| format!("{prefix}: blob past 64 KiB"))?);
        let mut b = n.as_bytes().to_vec();
        b.push(0);
        at += b.len();
        bytes.push((n.clone(), b));
    }
    let mut out = format!("# Offset of row i's name in {prefix}_NAME_BYTES ({what}).\n");
    let _ = writeln!(out, "var {prefix}_NAME_OFF: u16[{}] = {{", off.len());
    for (o, n) in off.iter().zip(names) {
        let _ = writeln!(out, "    {o},  # {n}");
    }
    out.push_str("};\n\n");
    out.push_str("# Every name, NUL-terminated, in row order (UTF-8).\n");
    let _ = writeln!(out, "var {prefix}_NAME_BYTES: u8[{at}] = {{");
    for (n, b) in &bytes {
        let list: Vec<String> = b.iter().map(u8::to_string).collect();
        let _ = writeln!(out, "    {},  # {n}", list.join(", "));
    }
    out.push_str("};\n");
    Ok(out)
}

/// `line` as a `#` comment, with no trailing space when it is empty (lint-clean).
pub fn comment_line(line: &str) -> String {
    if line.is_empty() {
        "#".to_string()
    } else {
        format!("# {line}")
    }
}

/// The fixed header of a generated `.cyr` file.
pub fn header(title: &str, provenance: &str, notes: &[&str]) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "# ushma — {title}");
    out.push_str("#\n");
    let _ = writeln!(
        out,
        "# GENERATED by scripts/oracle from {provenance}. Do not edit."
    );
    out.push_str(
        "# Regenerate: cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml\n",
    );
    for n in notes {
        out.push_str("#\n");
        for line in n.lines() {
            let _ = writeln!(out, "{}", comment_line(line));
        }
    }
    out.push('\n');
    out
}
