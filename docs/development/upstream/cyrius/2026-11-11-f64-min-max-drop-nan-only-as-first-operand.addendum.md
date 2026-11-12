<!-- Draft addendum: appended to cyrius/docs/development/issues/2026-11-11-f64-min-max-drop-nan-only-as-first-operand.md -->

## Also needed: ushma (2026-11-11, cyrius 6.7.18)

ushma's port measured the same asymmetry on its golden corpus. `f64_max` and `f64_min` differ from Rust on 59 of 600
rows each, and every one has a NaN second operand (`ushma/docs/development/port-parity.md` §2a). `f64_clamp` matches
Rust on all 617 of its rows with valid bounds.

ushma 1.4.0's library code (tests excluded) calls this surface at:

| Rust | Call sites |
|---|---|
| `f64::max` / `f64::min` | 11 / 6 |
| `is_finite()` / `is_nan()` | 21 / 5 |
| `f64::NAN` / `f64::INFINITY` / `f64::NEG_INFINITY` / `f64::EPSILON` | 11 / 1 / 1 / 1 |
| `clamp` | 9 (`f64_clamp` is fine) |

Its stand-in is kimiya's `_km_min`, `_km_max`, `_km_is_finite`, `_km_is_nan` and `_km_powi`, vendored into ushma's
`src/mathx.cyr` as `_us_*`. That is a fourth private copy, after impetus's, kimiya's and pravash's.
