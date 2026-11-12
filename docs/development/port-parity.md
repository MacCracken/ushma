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

### 2b. The corrected oracle (ADR 0004, 2026-11-11)

`python3 -I -B scripts/oracle/validate.py . --write-corrections` recomputes every class-C math row
correctly rounded with MPFR 4.2.2 (`scripts/oracle/mpfr_ref.py`: 53 bits, binary64's exponent
range, subnormalized). It writes `tests/golden/math.cr.tsv` for the rows where glibc 2.44's result
differs: **26 of 9,836**.

| Function | Rows | glibc misrounds | So the stdlib's differences from glibc (§2a) are |
|---|---|---|---|
| exp | 1,140 | 0 | the builtin's and polyfill's own errors (45 / 99 rows on x86) |
| exp_m1 | 1,010 | **23** | — (no stdlib function) |
| ln | 1,587 | 1 | the x87 builtin is correctly rounded on every row, the one glibc misrounds included; the polyfill is off on 29 |
| powf | 3,615 | 2 | ganita's own errors (373 of its 375), plus `pow(−1, DBL_MAX)`; on the 2 corrected rows it matches the correct value once and glibc's error once |
| atan | 1,624 | 0 | the polyfill's own errors (11); the x87 builtin and `ganita_f64_atan2(x, 1)` are correct on every row on x86 |
| tanh | 860 | 0 | **ganita's own errors (52 rows)**: `ganita_f64_tanh` is not correctly rounded, so the port needs `_us_cr_tanh` (M0 part 2c-3) |

The harness loads the file and compares each corrected row against the corrected outputs. While
the class-C rows are measured, it reports per function how many corrected rows the counterpart
gets right (`port = corrected`) and how many repeat glibc's error (`= glibc's`).

### 2c. The correctly rounded primitives on the corpus (M0 part 2c-2b, 2026-11-11)

`src/mathx.cyr` is vendored from kimiya (commit 5b14c85): 55 functions, the closure of ushma's 11 roots and their
fast and accurate entry points. `src/mathx_data.cyr` holds its 46 tables and constants, regenerated by
`scripts/oracle/mathx_tables.py` from MPFR with kimiya's assertions; `--check-kimiya` confirms every value is
bit-identical to kimiya's. The harness now compares them at **0 ulp** against the corrected oracle:

| Function (rows) | Port | Result |
|---|---|---|
| exp (1,140) | `_us_cr_exp` | 1,140 ok |
| exp_m1 (1,010) | `_us_cr_expm1` | 987 ok + 23 accepted (ADR 4): correct on every row glibc misrounds |
| ln (1,587) | `_us_cr_ln` | 1,586 ok + 1 accepted (ADR 4) |
| powf (3,615) | `_us_cr_pow` | 3,613 ok + 2 accepted (ADR 4); `pow(−1, DBL_MAX)` = 1, as C99 and glibc |
| atan (1,624) | `_us_cr_atan2(x, 1)` | 1,624 ok |
| max, min (600 each), is_finite, is_nan (34 each), powi (1,260) | `_us_max`, `_us_min`, `_us_is_finite`, `_us_is_nan`, `_us_powi` | all ok (class A / B, 0 ulp) |

No row matches glibc's misrounding instead of the correct value. The gated lines, and the summary
`14852 ok, 0 failed, 860 measured (tanh), 0 pending, 26 accepted`, are byte-identical on x86_64 ELF, `--bounds`,
PE and aarch64. The corpus is about 9,000 transcendental arguments. The dense MPFR sweep (kimiya's
`cr_verify.py`: ≥ 10⁶ arguments per primitive, hard cases, the accurate path forced, on every target) comes with
the correctly rounded tanh in M0 part 2c-3.

### 2d. tanh, and the dense MPFR sweep (M0 part 2c-3, 2026-11-11)

`_us_cr_tanh` (ADR 0004 amendment) is exact on all 860 golden rows, where `ganita_f64_tanh` is off on 52. With it,
every golden row is gated: `15712 ok, 0 failed, 0 measured, 0 pending, 26 accepted (ADR 4: 26)`.

`python3 -I -B scripts/oracle/cr_verify.py` (vendored from kimiya, plus ushma's tanh family) checks every
primitive against MPFR on about 10^6 seeded arguments per target, with the accurate path forced on 100,000 of
them. It also checks that the fast path agrees with the full result wherever its test passes, and that every target
gives the same bits and flags. The results are in §2e.

### 2e. The MPFR sweep on three targets (2026-11-11)

`python3 -I -B scripts/oracle/cr_verify.py` (MPFR 4.2.2; x86_64 ELF, x86_64 PE under wine, aarch64 under qemu):

| Primitive | Arguments | Mismatches | Accurate path forced | Fast path falls back on |
|---|---|---|---|---|
| exp | 1,096,454 | 0 | 100,045, 0 mismatches | 2.86·10⁻⁴ (mostly tiny \|x\|) |
| expm1 | 1,099,052 | 0 | 100,045, 0 | 5.78·10⁻³ (mostly the overflow end, which the fast path leaves to the accurate one) |
| ln | 1,001,201 | 0 | 100,039, 0 | 6.19·10⁻⁵ (mostly near 1) |
| pow | 1,306,099 | 0 | 100,537, 0 | 1.91·10⁻⁴ |
| atan2 | 1,014,417 | 0 | 100,411, 0 | 1.87·10⁻⁵ |
| tanh | 1,082,051 | 0 | 100,045, 0 | 1.29·10⁻⁵ (8 of the 14 at the 1 end) |

- **Fast path vs full result:** wherever the fast path's test passes, its result equals the full one, on every
  argument.
- **Across targets:** bits and flags are identical on all three.
- **x^2.5 (Redlich–Kwong's Tc^2.5):** 300,007 arguments, 0 mismatches.
- **Special values:** pow 537, atan2 361 and tanh 31, all equal to MPFR's. glibc and MPFR differ on two of the
  published hard pow cases (`pow(0.5951224589825083, 0.8256884680195185)`,
  `pow(2.6430276526426257, 0.46875)`), where glibc misrounds and the port matches MPFR.

So every class-C function ushma calls is correctly rounded and target-independent: exp, exp_m1, ln, powf, atan
(through atan2) and tanh. The class-C bound is therefore 0 ulp against the corrected oracle everywhere (ADR 0004).

## 3. Probes (M0 step 9)

Each probe is a small program in `scripts/probes/`, run on 6.7.18 (x86_64) on 2026-11-11. Batch 3a holds probes
2, 3, 4, 5, 6, 8 and 12; probes 1, 7, 9, 10 and 11 follow in batch 3b.

### Probe 12: f64 hygiene and the field-store miscompile (`p12_f64_hygiene.cyr`)

| Trap | Diagnosed |
|---|---|
| `0 - x`, `2 * x` | yes: *integer arithmetic with an f64 right operand* |
| an integer argument to an f64 parameter (`f(3)`) | yes: *an integer stored into an f64/f32 slot keeps its integer bits* |
| an integer stored into an f64 local | yes, the same |
| an untyped builtin result compared as an f64 | yes: *comparison mixes f64 and integer operands* |
| an f64 compared with an integer | yes, the same |
| `if (x)` on an f64 (true for −0.0) | **no**: CLAUDE.md now forbids it |

So the warning-free `cyrius test` gate catches six of the seven. The open
`field-store-of-struct-call-copies-a-field-argument` miscompile still reproduces at the pin (7 of 12 forms wrong, no
diagnostic). The workaround (bind the call to a local, then store) is a CLAUDE.md rule.

### Probe 6: names (`p6_names/`)

- **`mod us;` is unusable for a bundle.** Its scope runs to the end of the compilation unit, past the included file,
  so it mangled the consumer's own `main` into `us_main` and the build failed. No form ends it, and the compiler's
  own tests only ever use it at the end of a file. Inside the module a bare call to another of its fns doesn't
  resolve without a `use`, and globals aren't mangled (`us_GLOBALV` is undefined).
- **A collision with kimiya's bundle** (`include "../kimiya/dist/kimiya.cyr"` plus a ushma-style
  `ideal_gas_pressure`): *duplicate fn 'ideal_gas_pressure' (last definition wins)*, and here a hard error because
  kimiya's returns a `Result` and this one an f64. With matching signatures it would have compiled, and the last
  definition would silently have won.
- **The stdlib shares the namespace too.** cr_verify's probe had to rename its `run`, which collided with
  `lib/process.cyr`'s `run` (pulled in by sandhi's dependencies).

So P6 is a **written `us_` prefix** on every public function, with `US_` on constants, enum members and error codes
(ADR 0001 amendment).

### Probe 8: span cost (`p8_span_cost.bcyr`)

`sakshi_span_enter` / `_exit` have **no level check**: they read the clock, push the span stack and emit every call
(by default a formatted line to stderr). Around `n·R·T/P`:

| Variant | ns per call |
|---|---|
| plain | 8.1 (Rust 1.4.0's `state/ideal_gas_volume`: 3.7) |
| level-gated span (`if (level >= SK_DEBUG)` around enter and exit), DEBUG filtered out | 9.7 |
| sakshi span into its ring buffer (its cheapest output) | 192 |

So ushma's spans go through a level-gated `_us_span_enter` / `_us_span_exit`, the Cyrius form of Rust's filtered
`#[tracing::instrument(level = "debug")]`: about 1.6 ns when DEBUG is off.

### Probe 3: wide results (`p3_wide_results.bcyr`)

A CycleResult-shaped record (a kind and six f64s), 10^6 calls:

| Shape | ns per call | heap per call | error channel |
|---|---|---|---|
| `fn f(…, out: *R): Result`, Ok(0) / Err(code) | 14 | 0 | yes |
| `fn f(…): Result`, Ok(a fresh record) | 25–29 | **56 B, never freed** | yes |
| `fn f(…): R` by value | 13–14 | 0 | no |

So results come back through a caller's out-pointer with a value-form `Result`. All three gave correct values. (An
untyped word becomes a typed pointer through a one-line cast helper, impetus's idiom; a bare assignment warns.)

### Probe 4: point lists (`p4_point_lists.bcyr`)

A 200-point diagram, built 10^4 times:

| Shape | µs per diagram | heap per diagram |
|---|---|---|
| vec of 16-byte point records | 7.0 | 7,192 B |
| two f64 vecs, fresh | 6.6 | 7,984 B |
| two f64 vecs, reused (`vec_truncate(v, 0)`) | 5.3 | 0 |
| caller's flat f64 buffer, count returned | 2.0 | 0 |

So `cycle_ts_diagram` / `cycle_pv_diagram` take a caller buffer (with a too-small error) or caller vecs they
refill; a fresh allocation per call leaks about 8 KB.

### Probe 2: memory (`p2_memory.cyr`)

| Per call | Heap |
|---|---|
| hisab `calc_integral_simpson`, 200 steps (debye_cv's integral) | 0 B |
| hisab `calc_lerp` | 0 B |
| a 100-node explicit grid step, fresh scratch | 800 B, never freed |
| the same, a growable arena reset each step | 0 B |
| the same, scratch allocated once | 0 B |

So the open `no-scoped-allocator-for-alloc` issue (hisab allocating per call) doesn't reach ushma's two hisab calls,
which are scalar. Stepped state keeps its scratch in the object, or in an arena it resets each step.

### Probe 5: steam lookups (`p5_table_lookup.bcyr`)

On the generated column tables, the search and interpolation alone: saturated by T (binary search over 42 rows, 9
lerps) 77 ns; superheated by (T, P) (tier and row search, 3 lerps) 43 ns. Rust's whole functions take 18.2 ns and
31.4 ns. Cyrius doesn't inline, so the lerp calls dominate; this is M3's optimisation target (inline the lerps,
read each row pair once). Speed is measured, not gated.
