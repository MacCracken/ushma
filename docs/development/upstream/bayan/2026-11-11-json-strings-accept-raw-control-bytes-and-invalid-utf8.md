# JSON strings accept raw control bytes (U+0000–U+001F unescaped) and bytes that aren't UTF-8 — open

**Status:** 🟡 **OPEN**. Verified 2026-11-11 against bayan 1.11.1 (`965298d`), as bundled in the cyrius 6.7.18 stdlib.
The repro exits 1.
**Placement:** unpinned.
**Discovered:** 2026-11-11, by ushma's Rust → Cyrius port (M0 probe 11, `ushma/docs/development/port-parity.md` §3),
while replaying serde_json 1.0.154's verdicts through bayan's tree parser.
**Reporter:** ushma 2.0.0 port (`ushma/docs/development/ushma-port-plan.md` §8; ushma ADR 0005 amendment)
**Severity:** Medium. bayan accepts documents that RFC 8259 calls invalid and that serde_json refuses. A raw control
byte can't be refused by the caller either: after decoding it is indistinguishable from its escape.
**Affects:** bayan 1.11.1, `_jp_parse_string_a` (`src/json.cyr:1396`). A byte that is neither `"` nor `\` is copied
as is (:1496), so the tree, streaming and flat parsers all take it.

## Summary

RFC 8259 §7: "All Unicode characters may be placed within the quotation marks, except for the characters that MUST be
escaped: quotation mark, reverse solidus, and the control characters (U+0000 through U+001F)." §8.1: JSON text
exchanged between systems "MUST be encoded using UTF-8".

bayan accepts both kinds of violation:
- **Raw control bytes.** A raw 0x01, or even a raw tab or newline, inside a string. serde_json refuses: "control
  character (\u0000-\u001F) found while parsing a string".
- **Invalid UTF-8.** A byte such as 0xFF inside a string. serde_json refuses it when reading a `String`: "invalid
  unicode code point".

ushma's probe replayed 10 documents of each kind, and bayan parsed all 20.

For a caller that must match serde, invalid UTF-8 can at least be caught after the parse, because a valid escape
never produces it. A raw control byte can't be: the decoded `"a\u0001b"` and a raw `a<0x01>b` are the same bytes.

## Reproduction

`repros/2026-11-11-json-strings-accept-raw-control-bytes-and-invalid-utf8.cyr` parses `{"s":"a?b"}` with a raw 0x01,
then a raw tab, then 0xFF in place of the `?`. All three parse, and it exits 1.

## Root cause

`_jp_parse_string_a` checks only for `"` and `\`. It applies no test for bytes below 0x20 and no UTF-8 validation.

## Proposed fix

- **Control bytes.** In the string loop, refuse any byte below 0x20 with "control character in string" and its
  offset. This is the RFC rule, and serde_json's.
- **UTF-8.** Validate the decoded bytes, or the raw bytes between escapes, as UTF-8, rejecting overlong forms and
  surrogates as RFC 3629 does.

Both are refusals of invalid JSON, so they need no opt-in. If compatibility with existing lenient readers matters,
the UTF-8 check can sit behind the same kind of option the duplicate-key request
(`2026-11-09-json-duplicate-key-has-no-opt-in-for-ignored-members.md`) proposes, with strict as the default.

## Consumer-side workaround

ushma's decoders validate decoded strings as UTF-8. Raw control bytes need a pre-scan of the text for any byte below
0x20 inside a string (ushma ADR 0005 amendment), which repeats work the parser has already done.
