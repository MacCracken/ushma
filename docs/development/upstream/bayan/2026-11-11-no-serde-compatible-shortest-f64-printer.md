# `bayan_f64_to_json` isn't serde_json's (or Rust's) f64 text: not always the closest shortest digits, other exponent thresholds — open

**Status:** 🟡 **OPEN**. Verified 2026-11-11 against bayan 1.11.1 (`965298d`), as bundled in the cyrius 6.7.18 stdlib.
The repro exits 1. This is a stdlib-surface request with a verified implementation to seed it.
**Placement:** unpinned.
**Discovered:** impetus's port (its ADR 0010, port plan §8 "S2", not filed) and ushma's (M0 probe 11,
`ushma/docs/development/port-parity.md` §3). Filed by ushma on 2026-11-11.
**Reporter:** ushma 2.0.0 port (`ushma/docs/development/ushma-port-plan.md` §8; ushma ADR 0005)
**Severity:** Low–Medium. bayan's text round-trips, so no value is lost. But a port that must write the same JSON
bytes as the Rust program it replaces can't use bayan's printer. Each one carries a private copy of a shortest-digit
printer instead: impetus's `src/fmtf.cyr`, which ushma vendors at M1.
**Affects:** bayan 1.11.1, `bayan_f64_to_json` (Grisu2, the miloyip/rapidjson port), and with it `#derive(Serialize)`'s
f64 fields.

## Summary

serde_json 1.0.154 (zmij) and Rust's `Display` both print the **shortest** digits that read back to the same double,
and the **closest** of those. They differ only on an exact tie: zmij rounds to even, Rust's core rounds up. serde_json
writes positionally while the decimal exponent is in −5..=15, with `.0` on an integer, and otherwise as `d.ddde±x`.

`bayan_f64_to_json` is Grisu2. It always round-trips, but it isn't always shortest or closest, and its format
thresholds differ. ushma's probe compared it on 160,378 doubles: every f64 in ushma's 131 serde documents, every power
of two, the powers of ten and their neighbours, the form thresholds, 100,000 random bit patterns and 50,000 short
decimals. It differs from serde_json's text on **6,276**, with 1,402 of them carrying more digits. impetus's
`fmtf.cyr` (Burger & Dybvig on a small bignum, no tables) matched serde_json's text and Rust's Display on all 160,378.

## Reproduction

`repros/2026-11-11-no-serde-compatible-shortest-f64-printer.cyr` prints six doubles against serde_json's text:

```
serde_json 2.9545531576914354e-126   bayan 2.9545531576914355e-126   (2^-417: not the closest)
serde_json 2.9802322387695312e-8     bayan 2.9802322387695313e-8     (2^-25: an exact tie, rounded up)
serde_json 1e-6                      bayan 0.000001
serde_json 1.5e-6                    bayan 0.0000015
serde_json 1.8014398509481984e+16    bayan 18014398509481984.0
serde_json 5.902958103587057e+20     bayan 590295810358705700000.0
```

It exits 1.

## Root cause

Grisu2 without a fallback stops at "round-trips", not "shortest and closest". The positional range is a format choice.

## Proposed fix

Add the two Rust forms alongside the existing printer, which stays for current callers:
- `bayan_f64_to_json_serde(bits, buf)`: serde_json's text, with ties to even;
- `bayan_f64_to_display(bits, buf)`: Rust's `{}`, with ties up, always positional, no `.0`, and `NaN` / `inf` /
  `-inf`.

impetus's `src/fmtf.cyr` (ADR 0010, 400 lines, pinned by `tests/golden/fmtf.tsv` over every binade, the powers of ten
and a seeded sweep) implements both and can seed them, under bayan's names. `#derive(Serialize)` could then switch its
f64 fields to the serde form.

## Consumer-side workaround

Each port carries `fmtf.cyr`: impetus's, and ushma's vendored copy from M1, two copies already.
