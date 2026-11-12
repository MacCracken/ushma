# The JSON integer `-0` parses as `JTAG_INT` 0, so a reader that wants an f64 gets +0.0, not −0.0 — open

**Status:** 🟡 **OPEN**. Verified 2026-11-11 against bayan 1.11.1 (`965298d`), as bundled in the cyrius 6.7.18 stdlib.
The repro exits 1.
**Placement:** unpinned.
**Discovered:** 2026-11-11, by ushma's Rust → Cyrius port (M0 probe 11, `ushma/docs/development/port-parity.md` §3),
while replaying serde_json 1.0.154's documents through bayan's tree parser.
**Reporter:** ushma 2.0.0 port (`ushma/docs/development/ushma-port-plan.md` §8; ushma ADR 0005 amendment)
**Severity:** Low. Only a hand-written document carries `-0`: serde_json writes −0.0 as `-0.0`. But no caller can
recover the sign: the tree holds +0, and the text is gone.
**Affects:** bayan 1.11.1, `_jp_parse_number_a` (`src/json.cyr:1579`): the integer path, through `_jp_atoi` (:1643).

## Summary

`{"x":-0}` parses with `x` as `JTAG_INT` and payload 0. `bayan_json_v_float(x)` then gives `f64_from(0)`, +0.0.
serde_json reads `-0` into an f64 as −0.0, as JavaScript's `JSON.parse` does. (Python's `json` returns the integer 0.)
`-0.0` is fine: it has a `.`, so it is a FLOAT and keeps the sign. ushma's probe replayed 20 documents with an f64
field written `-0`, and each decoded to +0.0 where serde_json gives −0.0. Every other number form agreed.

## Reproduction

`repros/2026-11-11-json-integer-minus-zero-loses-its-sign.cyr` parses `{"x":-0,"y":-0.0}`:

```
-0   tag 2
-0   f64 bits negative (1 = sign kept): 0
-0.0 f64 bits negative (1 = sign kept): 1
```

It exits 1.

## Root cause

An integer token goes through `_jp_atoi` into an i64, which has no negative zero. The sign of a zero integer is
dropped there.

## Proposed fix

Parse an integer token that is a zero with a leading `-` as `JTAG_FLOAT` −0.0, so the token stays a number of the same
value. Integer readers are unaffected: `bayan_json_v_int` coerces a FLOAT with `f64_to`, which gives 0. Only the tag
changes (INT to FLOAT) for that one spelling. If the tag must stay INT, an alternative is a flag in the cell's
reserved word (+16), with `bayan_json_v_float` returning −0.0 when it is set.

## Consumer-side workaround

None from the tree. ushma's ADR 0005 maps the 20 rows as an open divergence until this lands.
