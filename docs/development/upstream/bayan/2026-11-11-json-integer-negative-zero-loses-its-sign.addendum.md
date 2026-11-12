<!-- Addendum: appended to bayan/docs/development/issues/2026-11-11-json-integer-negative-zero-loses-its-sign.md -->

## Also measured: ushma (2026-11-11, bayan 1.11.1)

ushma's port hit the same case (`ushma/docs/development/port-parity.md` §3, probe 11). It replayed serde_json
1.0.154's verdicts on documents for its 20 struct types that have an f64 field, each with that field written as `-0`.
serde_json read all 20 as −0.0, and bayan's tree gave +0.0 every time. It was the only number spelling where the two
disagreed on the value: `-0.0`, `1E2`, `1e+2`, `1e-400`, long and integer forms all matched. ushma maps the 20 rows as
an open divergence until the sign is kept (its ADR 0005 amendment). Either of the two fixes above would close them.
