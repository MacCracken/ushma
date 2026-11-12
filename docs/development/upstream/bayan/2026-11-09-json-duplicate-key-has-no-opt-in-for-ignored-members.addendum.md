<!-- Draft addendum: appended to bayan/docs/development/issues/2026-11-09-json-duplicate-key-has-no-opt-in-for-ignored-members.md -->

## Also measured: ushma (2026-11-11, bayan 1.11.1)

ushma's port replayed serde_json 1.0.154's documents for its 24 struct types through `bayan_json_v_parse_buf`
(`ushma/docs/development/port-parity.md` §3, probe 11). serde_json read each of these, and bayan refused it:

| Document | serde_json | bayan | Documents |
|---|---|---|---|
| `"zz":1,"zz":2` as an unknown top-level member | read (ignored) | `duplicate key` | 23 |
| `"zz":{"a":1,"a":2}`: a duplicate inside an unknown member | read (ignored) | `duplicate key` | 23 |
| an unknown member nested 200 deep | read (ignored, no depth limit) | `nesting too deep` (128) | 23 |

A duplicate of a known field was refused by both, and so was one inside a payload variant
(`{"Convective":{"h":10.0,"h":11.0,…}}`): 24 + 1 documents.

The third row is the same shape of difference: bayan's cap applies inside members the caller ignores, where
serde_json's `IgnoredAny` skip has none. ushma would keep the depth cap regardless (it is a resource guard), and
accepts both rows as mapped divergences in its ADR 0005 amendment. The opt-in this issue proposes would cover the
first two.
