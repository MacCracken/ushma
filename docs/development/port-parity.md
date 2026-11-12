# ushma port parity — measured baselines

Measurements the port plan's parity policy (§5) rests on, recorded as they are taken. Each
section names the probe or step that produced it and how to reproduce it.

## 1. Math counterparts and their semantics (M0 step 4, 2026-11-11)

Probe: `scripts/probes/p0_math.cyr`
(`cyrius build scripts/probes/p0_math.cyr build/p0_math && ./build/p0_math`; the same with
`--aarch64` under `qemu-aarch64`). Rust's side: `tests/golden/math.tsv`.

### Every Rust method has a Cyrius 6.7.18 counterpart

| Rust | Cyrius present and called | Note |
|---|---|---|
| `sqrt`, `abs`, `ceil` | `f64_sqrt`, `f64_abs`, `f64_ceil` (builtins) | exact |
| `exp`, `ln` | `f64_exp`, `f64_ln` (builtins); `_f64_exp_polyfill`, `_f64_ln_polyfill` | same bits as each other at the probed point; the port uses correctly rounded `_us_cr_*` (ADR 0004) |
| `atan` | `f64_atan` (builtin), `_f64_atan_polyfill`, `ganita_f64_atan2(x, 1)` | the three agree at the probed point |
| `powf` | `ganita_f64_pow` | **not correctly rounded**: `pow(2, 0.5)` = `0x3ff6a09e667f3bcc`, while √2 correctly rounded (`f64_sqrt(2)`) is `…bcd`. One ulp off on a trivial input, more evidence for ADR 0004 |
| `tanh` | `ganita_f64_tanh` | `tanh(0.5)` = `0x3fdd9353d7568af3`; correct rounding is measured against MPFR in step 8 |
| `exp_m1` | none in math or ganita | the port's `_us_cr_expm1` (ADR 0004) |
| `powi` | none | repeated multiplication in `__powidf2`'s order |
| `max`, `min`, `clamp` | `f64_max`, `f64_min`, `f64_clamp` (math.cyr) | semantics below |
| `is_nan`, `is_finite` | none; the idioms `x != x` and `f64_eq(x - x, 0.0)` | both give the right answer on NaN, ±∞ and finite values |

### min / max / clamp: where the stdlib and Rust differ

| Case | Cyrius 6.7.18 (x86_64) | Rust 1.99 (golden rows) | Port |
|---|---|---|---|
| `max(NaN, 1)` / `min(NaN, 1)` | 1 | 1 | — |
| `max(1, NaN)` / `min(1, NaN)` | **NaN** | **1** | `_us_max`/`_us_min` wrappers ignore a NaN in either slot |
| `max(+0, −0)` / `min(+0, −0)` | −0 (the second operand) | −0 (the second operand, `maxsd`/`minsd`) | same |
| `max(−0, +0)` / `min(−0, +0)` | +0 (the second operand) | +0 | same |
| `clamp(NaN, 0, 1)` | NaN | NaN | same |
| `clamp(−0, +0, 1)` | −0 | −0 | same |

On aarch64 the probe's output is identical except for the sign of the NaNs:
`0x7ff8000000000000` there, `0xfff8000000000000` on x86_64. Class A treats any NaN as
matching any NaN, so this doesn't affect parity.

`validate.py` independently confirms Rust's semantics on all 600 `f64::max` and 600 `f64::min`
golden rows: a NaN operand is ignored, and (±0, ∓0) gives the second operand.

## 1b. Generated tables on three targets (M0 step 6, 2026-11-11)

`python3 -I -B scripts/oracle/probe_tables.py . --all-targets` compiles all 1,052 generated
values with the Cyrius compiler: the constants, the four public tables, steam's private
tables and the other private constants. Each value is checked against the golden rows or an
independent Python parse of the Rust literals. The output is byte-identical on x86_64 ELF,
aarch64 (qemu) and x86_64 PE (wine). Decimal literals convert to the same bits on every
target.

## 2. Baselines on three targets (M0 step 8)

### 2a. The stdlib against glibc (2026-11-11, `tests/parity.tcyr`)

The parity harness dispatches each class-C math row to its stdlib counterpart and records the
ulp distance from the oracle (glibc 2.44) without gating it. These numbers measure how far the
stdlib is from **glibc**, not from the correct result. Which side is right on each differing
row comes from the corrected oracle (ADR 0004, MPFR, M0 part 2c), and the port replaces every
one of these with a correctly rounded `_us_cr_*` primitive.

Rows at 1 ulp from glibc (none is further, except the one `pow` row below):

| Function (rows) | Counterpart | x86_64 ELF | x86_64 PE | aarch64 |
|---|---|---|---|---|
| exp (1,140) | `f64_exp` builtin | 45 | 45 | 99 |
| | `_f64_exp_polyfill` | 99 | 99 | 99 |
| ln (1,587) | `f64_ln` builtin | 1 | 1 | 28 |
| | `_f64_ln_polyfill` | 28 | 28 | 28 |
| powf (3,615) | `ganita_f64_pow` | 375 + 1 non-finite | same | same |
| atan (1,624) | `f64_atan` builtin | 0 | 0 | 11 |
| | `_f64_atan_polyfill` | 11 | 11 | 11 |
| | `ganita_f64_atan2(x, 1)` | 0 | 0 | 11 |
| tanh (860) | `ganita_f64_tanh` | 52 | 52 | 52 |
| exp_m1 (1,010) | none | — | — | — |

- **Builtins vs polyfills.** On aarch64 the builtins are the polyfills, so the builtins'
  bits differ by target. The polyfills are target-independent but further from glibc. No
  stdlib option is both target-independent and as close to glibc as the x87 builtins: ADR
  0004's case.
- **The non-finite `pow` row** is `pow(−1, DBL_MAX)`. glibc gives 1 and `ganita_f64_pow`
  gives +∞: ganita's known defect (tanmatra ADR 0007, kimiya ADR 0004).
- **Class A.** `f64_sqrt`, `f64_abs`, `f64_ceil` and `f64_clamp` match Rust on every golden
  row (1,074 / 374 / 1,258 / 617), so `f64_clamp` has Rust's semantics for valid bounds.
  `f64_max` / `f64_min` differ on 59 of 600 rows each, every one a NaN second operand (the
  shadow lines). M1's `_us_max` / `_us_min` fix that.

The harness's summary is identical on x86_64 ELF, `--bounds`, PE (wine) and aarch64 (qemu):
`6 modules, 0 failing | rows: 3374 ok, 0 failed, 8826 measured, 3538 pending`.

### 2b. The correctly rounded primitives against MPFR

*Pending (M0 part 2c-2).*

## 3. Probes 1–12 (M0 step 9)

*Pending (M0 part 3).*
