<!-- Draft addendum: appended to cyrius/docs/development/issues/2026-11-11-derive-json-no-vector-array-or-tagged-variant-shapes.md -->

## Also measured: ushma (2026-11-11, cyrius 6.7.18)

ushma's port ran `#derive(Serialize)` on three of its 28 serde types (StatePoint, ThermalMaterial, Phase) against
serde_json 1.0.154's own documents (`ushma/docs/development/port-parity.md` §3, probe 11). It found two more shapes
a port needs, besides the vector-array and tagged-variant ones above.

- **Variant names differ from the constants.** Enum members are global, so a library prefixes them
  (`US_PHASE_SOLID`), and the derived codec writes and reads that name: `"US_PHASE_SOLID"`. serde writes `"Solid"`,
  and the derived reader can't read it. A per-variant rename, or a prefix to strip (`#derive(Serialize,
  strip_prefix = "US_PHASE_")`), would cover every prefixed library.
- **The sequence form.** serde accepts a struct written as an array, its fields by position
  (`[300.0,101325.0,…]`); the derived reader returns 0.
- **Smaller:**
  - a variant name over 31 bytes is a compile error (`#derive field name too long (max 31 bytes)`);
  - the f64 text is bayan's printer, not serde's (`9223372036854776000.0` against `9.223372036854776e+18`);
  - the decoders accept documents serde refuses: truncated, trailing data, a BOM, trailing commas, duplicate and
    missing fields, `null`, NaN and Infinity, `1e400`, leading zeros, `+1.5`, `.5`, `1.`, and a string or bool for an
    f64. That is 49 of 51 disagreements; `2026-11-11-derive-writer-faults-on-a-zero-str-field.md` files the one
    that faults.

ushma hand-writes its codecs on bayan's tree for now (its ADR 0005), as impetus and pravash do.
