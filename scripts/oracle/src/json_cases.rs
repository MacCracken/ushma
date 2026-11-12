//! The serde_json document machinery shared by the M0 probe (`src/bin/probe_json.rs`, probe 11) and the golden
//! corpus's JSON rows (`m_material.rs` and each later module's): a strict reader of what serde_json accepts (its
//! token stream, and a compact document's top-level members), and the variants of a base document, malformed,
//! duplicate, unknown and missing fields, number, string and shape forms, that every type's decoder rows replay.
//! Each user takes the parts it needs, hence the module-level `dead_code` allowance.
#![allow(dead_code)]

use std::fmt::Write as _;

pub fn pct(bytes: &[u8]) -> String {
    let mut s = String::new();
    for &b in bytes {
        if b.is_ascii_alphanumeric() || b == b'.' || b == b'_' || b == b'-' {
            s.push(b as char);
        } else {
            let _ = write!(s, "%{b:02X}");
        }
    }
    if s.is_empty() {
        s.push('%');
    }
    s
}

pub fn hex(bytes: &[u8]) -> String {
    let mut s = String::with_capacity(bytes.len() * 2);
    for &b in bytes {
        let _ = write!(s, "{b:02x}");
    }
    s
}

/// A strict reader of the JSON serde_json accepts, emitting the token stream and recording the top level's members.
pub struct Tok<'a> {
    pub b: &'a [u8],
    pub i: usize,
    pub out: Vec<String>,
}

impl Tok<'_> {
    fn ws(&mut self) {
        while self.i < self.b.len() && matches!(self.b[self.i], b' ' | b'\t' | b'\n' | b'\r') {
            self.i += 1;
        }
    }

    fn string(&mut self) -> Option<Vec<u8>> {
        self.i += 1;
        let mut v = Vec::new();
        loop {
            let c = *self.b.get(self.i)?;
            self.i += 1;
            match c {
                b'"' => return Some(v),
                b'\\' => {
                    let e = *self.b.get(self.i)?;
                    self.i += 1;
                    match e {
                        b'"' | b'\\' | b'/' => v.push(e),
                        b'b' => v.push(8),
                        b'f' => v.push(12),
                        b'n' => v.push(10),
                        b'r' => v.push(13),
                        b't' => v.push(9),
                        b'u' => {
                            // A lone surrogate (which serde_json skips over as JSON but refuses as a String) is
                            // U+FFFD here, as RFC 8259 readers commonly decode it.
                            let mut cp = self.hex4()?;
                            if (0xD800..0xDC00).contains(&cp) {
                                let pair = self.b.get(self.i) == Some(&b'\\') && self.b.get(self.i + 1) == Some(&b'u');
                                let lo = if pair {
                                    let save = self.i;
                                    self.i += 2;
                                    let lo = self.hex4()?;
                                    if !(0xDC00..0xE000).contains(&lo) {
                                        self.i = save;
                                    }
                                    lo
                                } else {
                                    0
                                };
                                cp = if (0xDC00..0xE000).contains(&lo) {
                                    0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00)
                                } else {
                                    0xFFFD
                                };
                            } else if (0xDC00..0xE000).contains(&cp) {
                                cp = 0xFFFD;
                            }
                            let ch = char::from_u32(cp)?;
                            let mut buf = [0u8; 4];
                            v.extend_from_slice(ch.encode_utf8(&mut buf).as_bytes());
                        }
                        _ => return None,
                    }
                }
                _ => v.push(c),
            }
        }
    }

    fn hex4(&mut self) -> Option<u32> {
        let s = std::str::from_utf8(self.b.get(self.i..self.i + 4)?).ok()?;
        self.i += 4;
        u32::from_str_radix(s, 16).ok()
    }

    fn value(&mut self) -> Option<()> {
        self.ws();
        match *self.b.get(self.i)? {
            b'{' => {
                self.i += 1;
                self.out.push("{".into());
                self.ws();
                if self.b.get(self.i) == Some(&b'}') {
                    self.i += 1;
                    self.out.push("}".into());
                    return Some(());
                }
                loop {
                    self.ws();
                    if self.b.get(self.i) != Some(&b'"') {
                        return None;
                    }
                    let k = self.string()?;
                    self.out.push(format!("k:{}", pct(&k)));
                    self.ws();
                    if self.b.get(self.i) != Some(&b':') {
                        return None;
                    }
                    self.i += 1;
                    self.value()?;
                    self.ws();
                    match *self.b.get(self.i)? {
                        b',' => self.i += 1,
                        b'}' => {
                            self.i += 1;
                            self.out.push("}".into());
                            return Some(());
                        }
                        _ => return None,
                    }
                }
            }
            b'[' => {
                self.i += 1;
                self.out.push("[".into());
                self.ws();
                if self.b.get(self.i) == Some(&b']') {
                    self.i += 1;
                    self.out.push("]".into());
                    return Some(());
                }
                loop {
                    self.value()?;
                    self.ws();
                    match *self.b.get(self.i)? {
                        b',' => self.i += 1,
                        b']' => {
                            self.i += 1;
                            self.out.push("]".into());
                            return Some(());
                        }
                        _ => return None,
                    }
                }
            }
            b'"' => {
                let s = self.string()?;
                self.out.push(format!("s:{}", pct(&s)));
                Some(())
            }
            b't' if self.b[self.i..].starts_with(b"true") => {
                self.i += 4;
                self.out.push("b:1".into());
                Some(())
            }
            b'f' if self.b[self.i..].starts_with(b"false") => {
                self.i += 5;
                self.out.push("b:0".into());
                Some(())
            }
            b'n' if self.b[self.i..].starts_with(b"null") => {
                self.i += 4;
                self.out.push("n".into());
                Some(())
            }
            _ => {
                let s = self.i;
                while self.i < self.b.len() && matches!(self.b[self.i], b'-' | b'+' | b'.' | b'e' | b'E' | b'0'..=b'9') {
                    self.i += 1;
                }
                let t = std::str::from_utf8(&self.b[s..self.i]).ok()?;
                if t.contains(['.', 'e', 'E']) {
                    let x: f64 = t.parse().ok()?;
                    self.out.push(format!("f:{:016x}", x.to_bits()));
                } else {
                    if t.is_empty() {
                        return None;
                    }
                    self.out.push(format!("i:{t}"));
                }
                Some(())
            }
        }
    }
}

pub fn stream(doc: &[u8]) -> Option<String> {
    let mut t = Tok { b: doc, i: 0, out: Vec::new() };
    t.value()?;
    t.ws();
    if t.i != doc.len() {
        return None;
    }
    Some(t.out.join(" "))
}

/// The top-level object's members of a compact serde_json document: (key, value text) in order.
pub fn members(doc: &str) -> Option<Vec<(String, String)>> {
    let b = doc.as_bytes();
    if b.first() != Some(&b'{') {
        return None;
    }
    let mut out = Vec::new();
    let mut i = 1;
    while i < b.len() - 1 {
        let ks = i;
        let mut t = Tok { b, i, out: Vec::new() };
        t.string()?;
        let key = doc[ks..t.i].to_string();
        t.i += 1;
        let vs = t.i;
        t.value()?;
        out.push((key, doc[vs..t.i].to_string()));
        i = t.i + 1;
    }
    Some(out)
}

pub fn obj(m: &[(String, String)]) -> String {
    let body: Vec<String> = m.iter().map(|(k, v)| format!("{k}:{v}")).collect();
    format!("{{{}}}", body.join(","))
}

pub fn is_float(v: &str) -> bool {
    v.contains(['.', 'e', 'E']) && !v.starts_with(['"', '[', '{'])
}

pub fn is_int(v: &str) -> bool {
    !v.is_empty() && v.bytes().all(|c| c.is_ascii_digit())
}

/// Variants of one base document. Each is (kind, bytes).
pub fn variants(base: &str) -> Vec<(String, Vec<u8>)> {
    let mut v: Vec<(String, Vec<u8>)> = Vec::new();
    let mut add = |k: &str, s: Vec<u8>| v.push((k.to_string(), s));
    add("ws", format!("\t\n {base} \r\n").into_bytes());
    add("trunc", base.as_bytes()[..base.len() - 1].to_vec());
    add("trailing", format!("{base}x").into_bytes());
    add("two_docs", format!("{base} {base}").into_bytes());
    add("bom", [&[0xEF, 0xBB, 0xBF][..], base.as_bytes()].concat());
    if let Some(m) = members(base) {
        let ins = |pre: &str| format!("{{{pre}{}", &base[1..]).into_bytes();
        add("trailing_comma", format!("{},}}", &base[..base.len() - 1]).into_bytes());
        add("dup_first", ins(&format!("{}:{},", m[0].0, m[0].1)));
        add("unknown", ins("\"zz_unknown\":{\"a\":[1,2,{\"b\":null}],\"c\":\"x\"},"));
        add("unknown_dup", ins("\"zz\":1,\"zz\":2,"));
        add("unknown_nested_dup", ins("\"zz\":{\"a\":1,\"a\":2},"));
        add("deep_100", ins(&format!("\"zz\":{}{},", "[".repeat(100), "]".repeat(100))));
        add("deep_200", ins(&format!("\"zz\":{}{},", "[".repeat(200), "]".repeat(200))));
        add("missing_first", obj(&m[1..]).into_bytes());
        add("empty_obj", b"{}".to_vec());
        let mut r = m.clone();
        r.reverse();
        add("reorder", obj(&r).into_bytes());
        let vals: Vec<&str> = m.iter().map(|(_, x)| x.as_str()).collect();
        add("seq_form", format!("[{}]", vals.join(",")).into_bytes());
        add("seq_short", format!("[{}]", vals[..vals.len() - 1].join(",")).into_bytes());
        add("seq_long", format!("[{},0]", vals.join(",")).into_bytes());
        let with = |i: usize, val: &str| {
            let mut c = m.clone();
            c[i].1 = val.to_string();
            obj(&c).into_bytes()
        };
        if let Some(i) = m.iter().position(|(_, x)| is_float(x)) {
            for (k, val) in [
                ("f_null", "null"),
                ("f_int", "1"),
                ("f_neg_int", "-7"),
                ("f_huge", "1e400"),
                ("f_tiny", "1e-400"),
                ("f_neg_zero_int", "-0"),
                ("f_neg_zero", "-0.0"),
                ("f_exp_upper", "1E2"),
                ("f_exp_plus", "1e+2"),
                ("f_lead_zero", "01.5"),
                ("f_plus", "+1.5"),
                ("f_dot_first", ".5"),
                ("f_dot_last", "1."),
                ("f_nan", "NaN"),
                ("f_inf", "Infinity"),
                ("f_string", "\"1.5\""),
                ("f_bool", "true"),
                ("f_long", "3.14159265358979323846264338327950288419716939937510582097494459"),
                ("f_halfway", "9007199254740993"),
                ("f_i64_over", "9223372036854775808"),
            ] {
                add(k, with(i, val));
            }
        }
        if let Some(i) = m.iter().position(|(_, x)| is_int(x)) {
            for (k, val) in [
                ("u_max", "18446744073709551615"),
                ("u_over", "18446744073709551616"),
                ("u_i64_over", "9223372036854775808"),
                ("u_neg", "-1"),
                ("u_float", "3.0"),
                ("u_exp", "3e0"),
                ("u_null", "null"),
            ] {
                add(k, with(i, val));
            }
        }
        if let Some(i) = m.iter().position(|(_, x)| x.starts_with('"')) {
            for (k, val) in [
                ("s_escapes", "\"\\u00e9\\ud83d\\ude00\\n\\t\\\"\\\\\\/\\b\\f\\r\\u0000\""),
                ("s_lone_high", "\"\\ud800\""),
                ("s_lone_low", "\"\\udc00\""),
                ("s_bad_escape", "\"\\x41\""),
                ("s_null", "null"),
                ("s_number", "42"),
            ] {
                add(k, with(i, val));
            }
            let mut raw = m.clone();
            raw[i].1 = "\"\u{e9}\"".to_string();
            add("s_raw_utf8", obj(&raw).into_bytes());
            let mut b = with(i, "\"A\"");
            if let Some(p) = b.iter().position(|&c| c == b'A') {
                b[p] = 0xFF;
                add("s_bad_utf8", b.clone());
                b[p] = 0x01;
                add("s_raw_ctrl", b);
            }
        }
        if let Some(i) = m.iter().position(|(_, x)| x.starts_with('[')) {
            let a = &m[i].1;
            let inner = &a[1..a.len() - 1];
            if let Some(cut) = inner.rfind(',') {
                add("a_short", with(i, &format!("[{}]", &inner[..cut])));
            }
            add("a_long", with(i, &format!("[{inner},0.5]")));
            add("a_null", with(i, "null"));
        }
    } else if base.starts_with('"') {
        let name = &base[1..base.len() - 1];
        add("e_obj_null", format!("{{{base}:null}}").into_bytes());
        add("e_unknown", b"\"Nope\"".to_vec());
        add("e_lower", format!("\"{}\"", name.to_lowercase()).into_bytes());
        add("e_int", b"0".to_vec());
        add("e_obj_two", format!("{{{base}:null,\"Nope\":null}}").into_bytes());
        add("e_seq", format!("[{base}]").into_bytes());
    }
    v
}
