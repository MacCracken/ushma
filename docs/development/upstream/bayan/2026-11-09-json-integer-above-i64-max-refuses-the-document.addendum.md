<!-- Draft addendum: appended to bayan/docs/development/issues/2026-11-09-json-integer-above-i64-max-refuses-the-document.md -->

## Also measured: ushma (2026-11-11, bayan 1.11.1)

ushma's port hit both halves of the proposed fix (`ushma/docs/development/port-parity.md` §3, probe 11):
- **(1) u64.** serde_json reads `usize` fields at 2^63 and at u64::MAX (`18446744073709551615`): ThermalNetwork's
  `num_nodes` and a grid's dimensions. bayan refuses both documents (4).
- **(2) beyond i64, for an f64 field.** serde_json reads the integer `9223372036854775808` in an f64 field as
  9.223372036854776e+18, its correctly rounded value. bayan refuses each document (20).

The second half matters to every Rust port with f64 fields: serde accepts any JSON integer for an f64. (2)'s
`JTAG_FLOAT` with the correctly rounded value is exactly what an f64 reader needs; `_jp_atof` already gives it. ushma
maps both as divergences until this lands (its ADR 0005 amendment).
