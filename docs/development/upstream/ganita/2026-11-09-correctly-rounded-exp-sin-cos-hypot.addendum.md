<!-- Draft addendum: appended to ganita/docs/development/issues/2026-11-09-correctly-rounded-exp-sin-cos-hypot.md -->

## A fourth consumer: ushma (2026-11-11)

ushma's Rust → Cyrius port needs the same tier, by its ADR 0004 ("correct rounding is the accuracy target"). That
makes it the fourth consumer after impetus, kimiya and pravash.

- **Two more functions: tanh and expm1.**
  - **tanh.** `ganita_f64_tanh` is not correctly rounded.
  - **expm1.** ganita has no public expm1 (or log1p) at all; `_gn_expm1` is private. ushma's `partition_vibrational`
    and `debye_cv` need it, and so do kimiya's decay chains.
- **Its measurements** cover ushma's golden corpus (`ushma/tests/golden/math.tsv`, the 9,836 class-C rows), against
  the correctly rounded value (glibc's result, or MPFR's where glibc 2.44 misrounds, on 26 rows: exp_m1 23, ln 1,
  pow 2). They are x86_64 ELF numbers; on aarch64 the builtins are the polyfills.

  | Function (rows) | Counterpart | Rows off (all 1 ulp) |
  |---|---|---|
  | exp (1,140) | `f64_exp` builtin / `_f64_exp_polyfill` | 45 / 99 |
  | ln (1,587) | `f64_ln` builtin / `_f64_ln_polyfill` | 0 / 29 |
  | pow (3,615) | `ganita_f64_pow` | 375, plus `pow(-1, DBL_MAX)` = +∞ (this repo's open defect) |
  | atan (1,624) | `f64_atan` builtin / `_f64_atan_polyfill` / `ganita_f64_atan2(x, 1)` | 0 / 11 / 0 |
  | **tanh** (860) | `ganita_f64_tanh` | **52** |
  | **exp_m1** (1,010) | none in ganita or the stdlib | — |
- **The stand-in** (ushma `src/mathx.cyr`, M0 parts 2c-2b and 2c-3):
  - kimiya's correctly rounded exp, expm1, ln, pow and atan2, vendored mechanically at kimiya `5b14c85`;
  - **ushma's own `_us_cr_tanh`**: tanh a = E / (E + 2) with E = expm1(2a), a double-double fast path with a proven
    rounding test, and a triple-double accurate path. `TANH_ONE`, the threshold past which tanh rounds to ±1, is
    generated from MPFR.

  ushma's `scripts/oracle/cr_verify.py` (kimiya's, plus a tanh family) checks them against MPFR 4.2.2. Each has 0
  mismatches, with identical bits and flags on x86_64 ELF, x86_64 PE and aarch64:

  | Primitive | Arguments | Accurate path forced |
  |---|---|---|
  | exp | 1,096,454 | 100,045 |
  | expm1 | 1,099,052 | 100,045 |
  | ln | 1,001,201 | 100,039 |
  | pow | 1,306,099 | 100,537 |
  | atan2 | 1,014,417 | 100,411 |
  | tanh | 1,082,051 | 100,045 |

  So tanh and expm1 for the tier can be seeded from ushma's and kimiya's verified code, as exp, pow and atan2 are
  from impetus's.
