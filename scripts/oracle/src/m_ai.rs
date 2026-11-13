//! `ai.tsv` (M7c): Rust 1.4.0's `ai`: the configs' defaults and JSON, and `DaimonClient::register_agent`'s requests and
//! outcomes against a loopback stub.
//!
//! A config is its endpoint `s:` and its key (`s:`, or `n:` for None). A registration row is the endpoint (`STUB` for
//! the stub's `127.0.0.1:<port>`, `REFUSED` for a port nothing listens on), the key, and the stub's scripted answer:
//! its status, the body it gives the oracle and the body it gives the port. The two bodies differ only where ADR 0002
//! row 6 does: daimon 2.5.6 sends an integer `id` where 1.4.0 reads a string `agent_id`. The stub answers at the
//! client's own path (1.4.0's `/v1/agents/register`, the port's `/v1/agents`) and 404 anywhere else. The row's outputs
//! are the request the stub received (`n:` when none came) and the outcome: the agent id `s:`, or the error.
//!
//! A request is its method, its path, the five headers 1.4.0 sends by value (host with the stub's address as `STUB`,
//! content-type, content-length, authorization, accept; `n:` for one not sent), the names of any others, and the body.
//! Class A throughout.

use std::io::{ErrorKind, Read, Write};
use std::net::TcpListener;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::thread;
use std::time::{Duration, Instant};

use ushma::ai::{DaimonClient, DaimonConfig, HooshConfig};

use crate::json_cases::variants;
use crate::tsv::{Class, Module, Tok, err, i, s};

/// 1.4.0's registration path; the port's is `/v1/agents` (ADR 0002 row 6).
const PATH: &str = "/v1/agents/register";

fn key_tok(k: Option<&str>) -> Tok {
    k.map_or(Tok::None, s)
}

fn dc_toks(c: &DaimonConfig) -> Vec<Tok> {
    vec![s(c.endpoint.clone()), key_tok(c.api_key.as_deref())]
}

fn hc_toks(c: &HooshConfig) -> Vec<Tok> {
    vec![s(c.endpoint.clone())]
}

/// The position just past the blank line ending a request's head, if it has one.
fn head_end(b: &[u8]) -> Option<usize> {
    b.windows(4).position(|w| w == b"\r\n\r\n").map(|p| p + 4)
}

/// The value of header `name` (case-insensitive) in a head's lines.
fn header<'a>(lines: &'a [(String, String)], name: &str) -> Option<&'a str> {
    lines
        .iter()
        .find(|(k, _)| k == name)
        .map(|(_, v)| v.as_str())
}

/// Reads one request from the stub's connection: its head and the body its content-length gives.
fn read_request(c: &mut std::net::TcpStream) -> Vec<u8> {
    let mut b = Vec::new();
    let mut chunk = [0u8; 4096];
    while head_end(&b).is_none() {
        match c.read(&mut chunk) {
            Ok(0) | Err(_) => return b,
            Ok(n) => b.extend_from_slice(&chunk[..n]),
        }
    }
    let end = head_end(&b).expect("head");
    let lines = parse_head(&b[..end]).1;
    let len: usize = header(&lines, "content-length")
        .and_then(|v| v.trim().parse().ok())
        .unwrap_or(0);
    while b.len() < end + len {
        match c.read(&mut chunk) {
            Ok(0) | Err(_) => break,
            Ok(n) => b.extend_from_slice(&chunk[..n]),
        }
    }
    b
}

/// A request head as its request line and its (lowercase name, value) header lines.
fn parse_head(head: &[u8]) -> (String, Vec<(String, String)>) {
    let text = String::from_utf8_lossy(head);
    let mut lines = text.split("\r\n");
    let first = lines.next().unwrap_or_default().to_string();
    let hs = lines
        .filter(|l| !l.is_empty())
        .filter_map(|l| l.split_once(':'))
        .map(|(k, v)| (k.trim().to_ascii_lowercase(), v.trim().to_string()))
        .collect();
    (first, hs)
}

/// A captured request as its tokens, the stub's address written `STUB`.
fn request_toks(raw: &[u8], addr: &str) -> Vec<Tok> {
    let end = head_end(raw).unwrap_or(raw.len());
    let (first, hs) = parse_head(&raw[..end]);
    let mut parts = first.split(' ');
    let method = parts.next().unwrap_or_default().to_string();
    let path = parts.next().unwrap_or_default().to_string();
    let opt = |name: &str| header(&hs, name).map_or(Tok::None, |v| s(v.to_string()));
    let host = header(&hs, "host").map_or(Tok::None, |v| s(v.replace(addr, "STUB")));
    let len = header(&hs, "content-length").map_or(Tok::None, |v| {
        i(v.trim().parse::<i64>().expect("content-length"))
    });
    let known = [
        "host",
        "content-type",
        "content-length",
        "authorization",
        "accept",
    ];
    let mut others: Vec<String> = hs
        .iter()
        .map(|(k, _)| k.clone())
        .filter(|k| !known.contains(&k.as_str()))
        .collect();
    others.sort();
    let mut t = vec![
        s(method),
        s(path),
        host,
        opt("content-type"),
        len,
        opt("authorization"),
        opt("accept"),
        Tok::V(others.len()),
    ];
    t.extend(others.into_iter().map(s));
    t.push(s(String::from_utf8_lossy(&raw[end..]).into_owned()));
    t
}

/// One registration of 1.4.0's client with endpoint `endpoint` (`STUB` and `REFUSED` replaced) and key `key`, against
/// a stub answering `status` and `body` at `path_prefix` + 1.4.0's path: the request it received and the outcome.
fn register(endpoint: &str, key: Option<&str>, status: u16, body: &str) -> Vec<Tok> {
    let listener = TcpListener::bind("127.0.0.1:0").expect("bind");
    let addr = listener.local_addr().expect("addr").to_string();
    let refused = endpoint.contains("REFUSED");
    let url = endpoint.replace("STUB", &addr).replace("REFUSED", &addr);
    let prefix = url
        .trim_end_matches('/')
        .split_once(&addr)
        .map_or(String::new(), |(_, p)| p.to_string());
    let want = format!("{prefix}{PATH}");
    let stop = Arc::new(AtomicBool::new(false));
    let stub = if refused {
        drop(listener);
        None
    } else {
        listener.set_nonblocking(true).expect("nonblocking");
        let stop2 = Arc::clone(&stop);
        let body = body.to_string();
        Some(thread::spawn(move || -> Option<Vec<u8>> {
            let deadline = Instant::now() + Duration::from_secs(30);
            loop {
                match listener.accept() {
                    Ok((mut c, _)) => {
                        c.set_nonblocking(false).expect("blocking");
                        c.set_read_timeout(Some(Duration::from_secs(5)))
                            .expect("timeout");
                        let raw = read_request(&mut c);
                        let end = head_end(&raw).unwrap_or(raw.len());
                        let (first, _) = parse_head(&raw[..end]);
                        let path = first.split(' ').nth(1).unwrap_or_default();
                        let (st, b) = if path == want {
                            (status, body.as_str())
                        } else {
                            (404, "")
                        };
                        let resp = format!(
                            "HTTP/1.1 {st} X\r\nContent-Type: application/json\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{b}",
                            b.len()
                        );
                        let _ = c.write_all(resp.as_bytes());
                        return Some(raw);
                    }
                    Err(e) if e.kind() == ErrorKind::WouldBlock => {
                        if stop2.load(Ordering::SeqCst) || Instant::now() > deadline {
                            return None;
                        }
                        thread::sleep(Duration::from_millis(1));
                    }
                    Err(_) => return None,
                }
            }
        }))
    };
    let config = DaimonConfig {
        endpoint: url,
        api_key: key.map(str::to_string),
    };
    let client = DaimonClient::new(config).expect("client");
    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()
        .expect("runtime");
    let outcome = match rt.block_on(client.register_agent()) {
        Ok(id) => s(id),
        Err(e) => err(&e),
    };
    stop.store(true, Ordering::SeqCst);
    let raw = stub.and_then(|h| h.join().expect("stub"));
    let mut out = match raw {
        Some(r) => request_toks(&r, &addr),
        None => vec![Tok::None],
    };
    out.push(outcome);
    out
}

/// The scripted answers: (status, 1.4.0's body, the port's body). Where the two differ, the oracle's has 1.4.0's
/// `agent_id` string and the port's daimon 2.5.6's integer `id` (ADR 0002 row 6).
fn answers() -> Vec<(u16, &'static str, &'static str)> {
    vec![
        (
            201,
            r#"{"agent_id":"7","name":"ushma"}"#,
            r#"{"id":7,"name":"ushma","type":"User","status":0}"#,
        ),
        (200, r#"{"agent_id":"12"}"#, r#"{"id":12}"#),
        (
            201,
            r#"{"agent_id":"9223372036854775807"}"#,
            r#"{"id":9223372036854775807}"#,
        ),
        (201, r#"{"agent_id":"-3"}"#, r#"{"id":-3}"#),
        (201, r#"{"agent_id":"0","agent_id2":"x"}"#, r#"{"id":0,"idx":5}"#),
        (401, r#"{"error":"unauthorized"}"#, r#"{"error":"unauthorized"}"#),
        (404, "", ""),
        (405, r#"{"error":"method not allowed"}"#, r#"{"error":"method not allowed"}"#),
        (409, r#"{"error":"agent limit reached"}"#, r#"{"error":"agent limit reached"}"#),
        (500, r#"{"error":"internal"}"#, r#"{"error":"internal"}"#),
        (302, "", ""),
        (201, r#"{"agent_id":"7""#, r#"{"id":7"#),
        (201, "ok", "ok"),
        (201, "", ""),
        (201, r#"{"name":"ushma"}"#, r#"{"name":"ushma"}"#),
        (201, r#"{"agent_id":7}"#, r#"{"id":"7"}"#),
        (201, r#"{"agent_id":null}"#, r#"{"id":null}"#),
        (201, r#"["7"]"#, "[7]"),
        (201, r#"{"agent_id":7.5}"#, r#"{"id":7.5}"#),
        (201, r#"{"agent_id":"7"} x"#, r#"{"id":7} x"#),
    ]
}

fn register_rows(m: &mut Module) {
    let row = |m: &mut Module, ep: &str, key: Option<&str>, a: (u16, &str, &str)| {
        let ins = [
            s(ep),
            key_tok(key),
            i(i64::from(a.0)),
            s(a.1),
            s(a.2),
        ];
        m.row(
            "ai::DaimonClient::register_agent",
            Class::A,
            &ins,
            &register(ep, key, a.0, a.1),
        );
    };
    // every answer, with Rust's probe key
    for a in answers() {
        row(m, "http://STUB", Some("probe-key"), a);
    }
    // the endpoint's trailing slashes and a path prefix, the key's forms, and a refused connection
    let ok = answers()[0];
    for ep in ["http://STUB/", "http://STUB///", "http://STUB/api", "http://STUB/api//"] {
        row(m, ep, Some("probe-key"), ok);
    }
    for key in [None, Some(""), Some("tab\tkey"), Some("k\u{e9}y"), Some("bad\r\nkey"), Some("del\u{7f}")] {
        row(m, "http://STUB", key, ok);
    }
    row(m, "http://REFUSED", Some("probe-key"), ok);
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

fn config_rows(m: &mut Module) {
    m.row(
        "ai::DaimonConfig::default",
        Class::A,
        &[],
        &dc_toks(&DaimonConfig::default()),
    );
    m.row(
        "ai::HooshConfig::default",
        Class::A,
        &[],
        &hc_toks(&HooshConfig::default()),
    );
    let dcs = [
        DaimonConfig::default(),
        DaimonConfig {
            endpoint: "http://example.com".into(),
            api_key: Some("key123".into()),
        },
        DaimonConfig {
            endpoint: "https://custom.host:9090".into(),
            api_key: Some("test-key".into()),
        },
        DaimonConfig {
            endpoint: String::new(),
            api_key: Some(String::new()),
        },
        DaimonConfig {
            endpoint: "a\"b\\c\t\u{1}\u{e9}\u{1d11e}".into(),
            api_key: None,
        },
    ];
    let mut docs = Vec::new();
    for c in &dcs {
        let doc = serde_json::to_string(c).expect("serialize");
        m.row("ai::DaimonConfig::to_json", Class::A, &dc_toks(c), &[s(doc.clone())]);
        docs.push(doc);
    }
    let base = docs[1].clone();
    add_variants(&mut docs, &base);
    docs.extend(
        [
            r#"{"endpoint":"x"}"#,
            r#"{"endpoint":"x","api_key":null}"#,
            r#"{"endpoint":"x","api_key":5}"#,
            r#"{"endpoint":"x","api_key":["k"]}"#,
            r#"{"api_key":"k"}"#,
            r#"{"endpoint":null,"api_key":"k"}"#,
            r#"{"endpoint":1}"#,
            r#"{"endpoint":"x","api_key":"k","extra":1}"#,
            r#"{"endpoint":"é𝄞","api_key":"\n"}"#,
            r#"{"endpoint":"\ud834"}"#,
            r#"["x","k"]"#,
            r#"["x",null]"#,
            r#"["x"]"#,
            r#"["x","k",1]"#,
            "[]",
            "{}",
            "null",
            "",
        ]
        .map(String::from),
    );
    for doc in &docs {
        m.row(
            "ai::DaimonConfig::from_json",
            Class::A,
            &[s(doc.clone())],
            &decode::<DaimonConfig>(doc, dc_toks),
        );
    }
    let hcs = [
        HooshConfig::default(),
        HooshConfig {
            endpoint: "http://h:1/x".into(),
        },
        HooshConfig {
            endpoint: "\u{e9}\"\u{7f}".into(),
        },
    ];
    let mut hdocs = Vec::new();
    for c in &hcs {
        let doc = serde_json::to_string(c).expect("serialize");
        m.row("ai::HooshConfig::to_json", Class::A, &hc_toks(c), &[s(doc.clone())]);
        hdocs.push(doc);
    }
    let base = hdocs[1].clone();
    add_variants(&mut hdocs, &base);
    hdocs.extend(
        [
            r#"{"endpoint":"x","api_key":"k"}"#,
            r#"{"endpoint":1}"#,
            r#"{"endpoint":null}"#,
            r#"["x"]"#,
            r#"["x","y"]"#,
            "[]",
            "{}",
            "null",
            "",
        ]
        .map(String::from),
    );
    for doc in &hdocs {
        m.row(
            "ai::HooshConfig::from_json",
            Class::A,
            &[s(doc.clone())],
            &decode::<HooshConfig>(doc, hc_toks),
        );
    }
}

pub fn emit() -> Module {
    let mut m = Module::new("ai");
    config_rows(&mut m);
    register_rows(&mut m);
    m
}
