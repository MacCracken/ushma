# 0005 — serde-compatible JSON codecs on bayan

**Status**: Accepted; amended 2026-11-11 (probe 11): the decoders' rules
**Date**: 2026-11-11

## Context

28 ushma types derive serde `Serialize`/`Deserialize`, and the crate has 30 serde round-trip test
sites. soorat consumes the visualization types as JSON, and `ai` sends JSON bodies. Everything
ports (ADR 0001), so the JSON does too.

Cyrius 6.7.18 offers `#derive(Serialize)` on bayan (guide L4548–4660), but its semantics are not
serde's:
- a missing key decodes to 0 instead of an error;
- u64, `Vec<Str>`, tuple and payload-enum fields are refused;
- it prints f64 with bayan's printer, which round-trips but isn't always shortest and formats
  differently from serde_json (zmij).

Its decoder also has three open issues (2026-11-07): it reads past a NUL, memsets a refused
allocation, and hangs on a stray brace. impetus reached the same decision (its ADR 0010).

## Decision

The port reads and writes **the same JSON documents as serde_json 1.0.154**, byte for byte:
- **Readers** work on bayan's value tree (`bayan_json_v_parse_buf`; `bayan_json_parse` is the flat
  parser). serde's semantics apply: required fields, absent `Option` fields, unknown fields
  ignored, duplicate fields refused (bayan 1.11.1 already refuses a duplicate key), an integer
  accepted for an f64, `null` refused for an f64. A malformed document is
  `Err(US_ERR_DESERIALIZE)`.
- **Writers** use `str_builder`.
- **f64 text** comes from a shortest-digit printer, `src/fmtf.cyr`, vendored from impetus, under
  zmij's rules for JSON. The same printer, under Rust's `Display` rules, prints `UshmaError`
  messages.
- **Parity.** Every document the oracle writes decodes to the oracle's values bit for bit, and
  the port's document for every golden value is serde_json's text byte for byte. Both are
  class A.
- **Hardening.** Every decoder is fuzzed (`cyrius fuzz --poison --bounds`).

## Amendment, 2026-11-11 (M0 probe 11): what the decoders add to bayan

Probe 11 read serde_json's 1,148 documents for the 28 types through `bayan_json_v_parse_buf` (`port-parity.md` §3).
Wherever both read a document, the tree is serde's, keys, strings and float bits alike, so the readers stay on
bayan's tree. Where the two differ, the decoders do this:
- **They match serde, by their own checks:**
  - a FLOAT that parsed to ±inf (`1e400`) is refused, as serde's `number out of range`;
  - a decoded string that isn't valid UTF-8 is refused (bayan passes raw bytes through);
  - a raw control character (U+0000–U+001F) inside a string is refused. bayan accepts it, and its decoded bytes
    can't tell it from an escape, so each decode pre-scans the text: a local stopgap until bayan refuses them
    (port plan §8).
- **The port refuses what serde reads.** Each document goes to `Err(US_ERR_DESERIALIZE)`, a mapped parity row
  (impetus's S3 and S4):
  - a duplicate key in an unknown field, or nested inside one (serde ignores it);
  - nesting past 128 levels inside an unknown field;
  - an integer ≥ 2^63, whether in an f64 field or a usize up to u64::MAX.
- **One open divergence: the integer `-0`.** serde reads it into an f64 as −0.0, and bayan as +0, a sign the tree
  can't recover. serde_json never writes it (−0.0 is `-0.0`), so only hand-written documents carry it. It is a mapped
  row until bayan keeps the sign (port plan §8).
- **The f64 oracle.** serde_json's default float parse is best-effort: it reads 35,970 of 160,320 of its own shortest
  texts back 1–2 ulp off. The parse precision is the consumer's serde_json choice, not ushma's, which only derives
  `Deserialize`. So the oracle crate builds serde_json with `float_roundtrip` (as impetus's does), and every decode
  agrees with bayan's correctly rounded `bayan_f64_parse`. Writing is the same in both modes.
- **The printer.** `fmtf.cyr` matched serde_json's text and Rust's Display on all 160,378 doubles; bayan's Grisu2
  printer differs on 6,276. `fmtf.cyr` needs one bit-cast helper from impetus's `mathx.cyr`, vendored with it under
  ushma's name.

## Consequences

- **Positive** — Rust consumers' JSON (soorat, stored grids) reads unchanged, and untrusted input
  never reaches the derive scanner's open defects.
- **Negative** — 28 hand-written codecs to own, and a second private copy of `fmtf.cyr` until
  bayan has a serde-compatible printer (plan §8).
- **Neutral** — M0 probe 11 records exactly where `#derive(Serialize)` differs, in case a later
  bayan closes the gap: it accepts 49 kinds of document serde refuses, refuses serde's sequence form, writes bayan's
  f64 text, faults re-encoding a struct whose `Str` was missing, and names enum variants by their `US_` constants.

## Alternatives considered

- **`#derive(Serialize)`, compared by parsed value only.** Rejected: its decoder accepts
  documents serde refuses, its text differs, and its open defects are reachable from untrusted
  input.
