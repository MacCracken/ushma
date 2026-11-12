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

*Pending (M0 part 3).*

## 3. Probes 1–12 (M0 step 9)

*Pending (M0 part 3).*
