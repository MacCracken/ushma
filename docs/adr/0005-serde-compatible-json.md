# 0005 — serde-compatible JSON codecs on bayan

**Status**: Accepted
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
- **Readers** work on `bayan_json_parse`'s tree. serde's semantics apply: required fields, absent
  `Option` fields, unknown fields ignored, duplicate fields refused (bayan 1.11.1 already refuses
  a duplicate key), an integer accepted for an f64, `null` refused for an f64. A malformed
  document is `Err(US_ERR_DESERIALIZE)`.
- **Writers** use `str_builder`.
- **f64 text** comes from a shortest-digit printer, `src/fmtf.cyr`, vendored from impetus, under
  zmij's rules for JSON. The same printer, under Rust's `Display` rules, prints `UshmaError`
  messages.
- **Parity.** Every document the oracle writes decodes to the oracle's values bit for bit, and
  the port's document for every golden value is serde_json's text byte for byte. Both are
  class A.
- **Hardening.** Every decoder is fuzzed (`cyrius fuzz --poison --bounds`).

## Consequences

- **Positive** — Rust consumers' JSON (soorat, stored grids) reads unchanged, and untrusted input
  never reaches the derive scanner's open defects.
- **Negative** — 28 hand-written codecs to own, and a second private copy of `fmtf.cyr` until
  bayan has a serde-compatible printer (plan §8).
- **Neutral** — M0 probe 11 records exactly where `#derive(Serialize)` differs, in case a later
  bayan closes the gap.

## Alternatives considered

- **`#derive(Serialize)`, compared by parsed value only.** Rejected: its decoder accepts
  documents serde refuses, its text differs, and its open defects are reachable from untrusted
  input.
