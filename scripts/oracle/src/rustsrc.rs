//! Reading `const` items out of the frozen Rust source (port plan §2.1, §4 step 6).
//!
//! A private constant isn't reachable through the oracle's API, so its value comes from its
//! literal in `rust-old/src`, parsed with Rust's own `f64::from_str` (correctly rounded), and is
//! proven by the golden rows that depend on it. `const_items` lists every `const` item a module
//! defines outside its test module (top level and fn-local), so the generators can prove they
//! emit or excuse each one (`t_private::complete`). Only comments and literals are read here;
//! nothing is evaluated beyond what Rust's const evaluation does for these items (IEEE `*` and
//! `/` left to right, integer `<<`).

/// One `const` item of the Rust source.
pub struct ConstItem {
    pub module: String,
    pub name: String,
    pub public: bool,
    pub ty: String,
    /// The initializer's text, `//` comments removed, without the trailing `;`.
    pub expr: String,
    /// The enclosing `fn`, for a fn-local const.
    pub in_fn: Option<String>,
    /// The `//` comment at the end of the item's line, if any.
    pub trailing: Option<String>,
}

impl ConstItem {
    pub fn path(&self) -> String {
        format!("{}::{}", self.module, self.name)
    }
}

fn strip_line_comment(line: &str) -> (&str, Option<&str>) {
    match line.find("//") {
        Some(at) => (&line[..at], Some(line[at + 2..].trim())),
        None => (line, None),
    }
}

/// The name of a `fn` declared on `line`, if it declares one.
fn fn_name(line: &str) -> Option<String> {
    let t = line.trim_start();
    let t = t.strip_prefix("pub ").unwrap_or(t);
    let t = t.strip_prefix("pub(crate) ").unwrap_or(t);
    let t = t.strip_prefix("const ").unwrap_or(t);
    let rest = t.strip_prefix("fn ")?;
    let end = rest.find(['(', '<'])?;
    Some(rest[..end].to_string())
}

/// Every `const NAME: TYPE = EXPR;` of `src` before its `#[cfg(test)]` module, in source order.
pub fn const_items(module: &str, src: &str) -> Result<Vec<ConstItem>, String> {
    let mut items = Vec::new();
    let mut lines = src.lines();
    // The enclosing fns: (name, brace depth at the line that declares it).
    let mut fns: Vec<(String, i64)> = Vec::new();
    let mut depth: i64 = 0;
    while let Some(raw) = lines.next() {
        if raw.trim_start().starts_with("#[cfg(test)]") {
            break;
        }
        let (code, trailing) = strip_line_comment(raw);
        let t = code.trim_start();
        let (public, rest) = if let Some(r) = t.strip_prefix("pub const ") {
            (true, Some(r))
        } else if let Some(r) = t.strip_prefix("pub(crate) const ") {
            (false, Some(r))
        } else {
            (false, t.strip_prefix("const "))
        };
        let item = rest
            .and_then(|r| r.split_once(':'))
            .filter(|(n, _)| !n.contains([' ', '(']));
        if let Some((name, after)) = item {
            let (ty, init) = after
                .split_once('=')
                .ok_or_else(|| format!("{module}::{name}: no `=`"))?;
            let mut expr = init.trim().to_string();
            let mut trail = trailing.map(str::to_string);
            // A multi-line initializer runs to the `;` that closes it. Its braces always balance,
            // so they are not counted.
            while !expr.trim_end().ends_with(';') {
                let next = lines
                    .next()
                    .ok_or_else(|| format!("{module}::{name}: unterminated"))?;
                let (c, _) = strip_line_comment(next);
                expr.push('\n');
                expr.push_str(c.trim_end());
                trail = None;
            }
            let expr = expr.trim_end().trim_end_matches(';').trim().to_string();
            let in_fn = fns
                .last()
                .filter(|(_, d)| depth > *d)
                .map(|(n, _)| n.clone());
            items.push(ConstItem {
                module: module.to_string(),
                name: name.trim().to_string(),
                public,
                ty: ty.trim().to_string(),
                expr,
                in_fn,
                trailing: trail,
            });
            continue;
        }
        if let Some(name) = fn_name(code) {
            fns.push((name, depth));
        }
        for ch in code.chars() {
            match ch {
                '{' => depth += 1,
                '}' => depth -= 1,
                _ => {}
            }
        }
        while fns
            .last()
            .is_some_and(|(_, d)| depth <= *d && code.contains('}'))
        {
            fns.pop();
        }
    }
    Ok(items)
}

/// One numeric literal (`1.380_649e-23`, `22_064_000.0`, `-0.238`), by Rust's `f64::from_str`.
pub fn literal_f64(tok: &str) -> Result<f64, String> {
    let t: String = tok.trim().chars().filter(|&c| c != '_').collect();
    let t = t.strip_suffix("f64").unwrap_or(&t);
    if t.is_empty() || !t.chars().all(|c| c.is_ascii_digit() || ".eE+-".contains(c)) {
        return Err(format!("`{tok}` is not a numeric literal"));
    }
    t.parse::<f64>().map_err(|e| format!("`{tok}`: {e}"))
}

/// An f64 initializer: a literal, or literals joined by `*` and `/`, evaluated left to right in
/// IEEE double as Rust's const evaluation does (`70.0 * 4184.0 / 86_400.0`).
pub fn eval_f64(expr: &str) -> Result<f64, String> {
    let mut acc: Option<f64> = None;
    let mut op = '*';
    let mut tok = String::new();
    let flush = |acc: &mut Option<f64>, op: char, tok: &str| -> Result<(), String> {
        let v = literal_f64(tok)?;
        *acc = Some(match (*acc, op) {
            (None, _) => v,
            (Some(a), '*') => a * v,
            (Some(a), '/') => a / v,
            _ => return Err(format!("operator {op}")),
        });
        Ok(())
    };
    for ch in expr.chars() {
        if (ch == '*' || ch == '/') && !tok.trim().is_empty() {
            flush(&mut acc, op, &tok)?;
            tok.clear();
            op = ch;
        } else if !ch.is_whitespace() {
            tok.push(ch);
        }
    }
    flush(&mut acc, op, &tok)?;
    acc.ok_or_else(|| format!("`{expr}`: empty"))
}

/// An integer initializer: a literal, or `a << b`.
pub fn eval_int(expr: &str) -> Result<i64, String> {
    let lit = |s: &str| -> Result<i64, String> {
        let t: String = s.trim().chars().filter(|&c| c != '_').collect();
        t.parse::<i64>().map_err(|e| format!("`{s}`: {e}"))
    };
    match expr.split_once("<<") {
        Some((a, b)) => {
            let shift = u32::try_from(lit(b)?).map_err(|e| format!("`{expr}`: {e}"))?;
            lit(a)?
                .checked_shl(shift)
                .ok_or_else(|| format!("`{expr}` overflows"))
        }
        None => lit(expr),
    }
}

/// The comma-separated items of an array initializer `[a, b, …]` or `&[a, b, …]`.
pub fn array_items(expr: &str) -> Result<Vec<String>, String> {
    let t = expr.trim().trim_start_matches('&').trim();
    let inner = t
        .strip_prefix('[')
        .and_then(|r| r.strip_suffix(']'))
        .ok_or_else(|| format!("not an array initializer: `{}`", &t[..t.len().min(40)]))?;
    Ok(inner
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(str::to_string)
        .collect())
}

/// The entries of a slice of struct literals `&[Ty { f: v, … }, …]`: each entry's
/// `(field, expr)` pairs in source order.
pub fn struct_entries(expr: &str, ty: &str) -> Result<Vec<Vec<(String, String)>>, String> {
    let mut out = Vec::new();
    let mut rest = expr;
    let open = format!("{ty} {{");
    while let Some(at) = rest.find(&open) {
        let body_start = at + open.len();
        let close = rest[body_start..]
            .find('}')
            .ok_or_else(|| format!("{ty}: unterminated entry"))?;
        let body = &rest[body_start..body_start + close];
        let mut fields = Vec::new();
        for part in body.split(',') {
            let part = part.trim();
            if part.is_empty() {
                continue;
            }
            let (f, v) = part
                .split_once(':')
                .ok_or_else(|| format!("{ty}: `{part}` is not `field: value`"))?;
            fields.push((f.trim().to_string(), v.trim().to_string()));
        }
        out.push(fields);
        rest = &rest[body_start + close + 1..];
    }
    Ok(out)
}
