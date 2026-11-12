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
2, 3, 4, 5, 6, 8 and 12. Batch 3b holds probes 1, 7, 9, 10 and 11. Their Rust halves are binaries of the oracle crate
(`scripts/oracle/src/bin/`: `probe_hisab`, `probe_ai`, `probe_json`), which drive the frozen 1.4.0 and Rust hisab 1.4.0
read-only.

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

### Probe 1: Simpson against Rust hisab (`p1_p9_hisab.cyr`, `probe_hisab`)

hisab 3.4.1's `calc_integral_simpson`, called through a function pointer, against `hisab::calc::integral_simpson` 1.4.0
on 3,000 seeded intervals. The integrand is x⁴/(1 + x)², with + − × ÷ only, so any difference is the integrator's:

| Rows | Bit-identical |
|---|---|
| all, 3,000 | 1,942 |
| n = 200 (debye_cv's), 1,152 | 763 |
| odd n, 371 | 248 |

The worst distance is 5 ulp. The two integrators sum in a different order, so `debye_cv` is **class D** (relative
≤ 1e-12, port plan §5), not B. The callback form works, and probe 2 measured it at 0 B per call.

### Probe 9: network solves (`p1_p9_hisab.cyr`, `probe_hisab`)

The probe generated 2,178 seeded connected networks: a chain through every node plus extra edges, R log-uniform over
1e-3 … 1e14, and one to three fixed nodes. Each system was built and row-scaled as 1.4.0's `ThermalNetwork::solve`
builds it, and the copy was asserted bit-identical to `solve()` itself:

| Systems | hisab 1.4.0 | `ganita_mat_gaussian_elim` | `ganita_mat_lu` + `_lu_solve` |
|---|---|---|---|
| 2,000 that 1.4.0 solves | solved | 1,741 bit-identical, **259 refused** | **2,000 bit-identical** (37,424 unknowns, 0 ulp) |
| 178 that 1.4.0 refuses | `SingularPivot` | 178 refused | **178 refused** |

So LU reproduces 1.4.0's answers and its singularity decisions on every system:
- **Why LU agrees.** LU refuses when |pivot| ≤ 1e-12 · max|aᵢⱼ|, and 1.4.0 when |pivot| < 1e-12 (`EPSILON_F64`). After
  row scaling, every diagonal is exactly 1 and every other entry is at most 1 in magnitude. So on ushma's matrices the
  two are the same test, apart from a pivot of exactly 1e-12.
- **Why elimination refuses more.** Its scale also spans the right-hand side, temperatures up to 1,000 K, so it refuses
  259 systems that 1.4.0 solves.
- **The leak.** `ganita_mat_lu_solve` allocates its forward-substitution vector on every call (`alloc(n * 8)`: 320 B at
  n = 40, never freed), even when the factors, pivots and solution are reused.
- **The decision.** P4 moves to LU, class B (ADR 0003 amendment). The scratch allocation is filed upstream as an `_a`
  variant (§8).

### Probe 7: the daimon registration on sandhi (`p7_sandhi_post.cyr`, `p7_daimon_stub.py`, `probe_ai`)

Both clients ran against a loopback stub that answers as daimon 2.5.6 does and records each raw request. The oracle's
`DaimonClient::register_agent` (`probe_ai`) and sandhi 1.13.7 compare like this:

| | Rust 1.4.0 (reqwest) | sandhi |
|---|---|---|
| request line | `POST /v1/agents/register` (the P8 defect: the stub's 404) | `POST /v1/agents` (P8) |
| body | 96 B, `{"capabilities":[…],"name":"ushma"}` (`json!` sorts the keys) | the same bytes, once the port writes serde's key order |
| host, content-type, content-length, authorization (`Bearer …`), accept (`*/*`) | sent | the same values; Accept is set explicitly |
| other headers | none | `User-Agent: sandhi/1.13.7`, `Accept-Encoding: identity`, `TE: trailers`, `Connection: close` |

Header names differ only in case (reqwest lowercases them). A 201 comes back with its body. A refused connection is
status 0 with `err_kind` 2.

**Timeouts.** The test server accepts the connection, then sends nothing for 3 s:

| Options | Result |
|---|---|
| `total_ms` 1000 | 200 after 3 s: the deadline doesn't bound a read that is blocked waiting for the first byte |
| `total_ms` + `read_ms` 1000 | `err_kind` 4 (timeout) after 1 s |
| `total_ms` + `idle_ms` 1000 | 200 after 3 s |

So the client sets `read_ms` as well as `total_ms`, both 30 s (reqwest's `timeout(30 s)`), and the gap is filed
upstream (§8). The results are the same on x86_64, PE (wine) and aarch64 (qemu).

### Probe 10: the manifest and both bundles (`p10_manifest.sh`)

The probe ran in a scratch copy (`build/p10/`):
- `[lib]` holds the generated `src/` modules plus a stand-in that calls hisab (Simpson, lerp), ganita (LU), sakshi, log
  and bayan.
- `[lib.ai]` holds the same, plus the sandhi client.

**sakshi twice.** Every resolve warns `refusing to overwrite stdlib leaf 'sakshi'`: ushma's `cyrius deps`, `distlib
--all` and `distlib --check`, and every consumer's own `cyrius deps`. The cause is hisab 3.4.1's `[deps.sakshi]` package,
which arrives beside the stdlib leaf that log and sandhi need.
- The three copies are byte-identical (sha256 `8a819554…`), and the stdlib copy wins: the vendored leaf, the 6.7.18
  snapshot and the 2.5.10 package.
- hisab is a production dependency here (debye_cv, the steam lerps), so prakash's way out (`scope = "test"`) doesn't
  apply. Only hisab can retire it (§8).

**The bundles.** `distlib --all` writes both, each with its compile-verified requires block, and `--check` reads both
as current:

| Bundle | Size | Sidecar leaves |
|---|---|---|
| `dist/ushma.cyr` | 235,305 B | `log bayan` |
| `dist/ushma-ai.cyr` | 236,341 B | `log bayan sandhi` |

The sidecars list only ushma's own leaves. math, ganita, sakshi and the rest come with hisab, which a consumer gets
through ushma's `[deps.hisab]` without declaring it (vendored as `lib/hisab.cyr`, commit-pinned).

**Consumers.** A consumer project declares `[deps.ushma]` with the bundle as its module, and either the sidecar's
leaves or `stdlib = []`. Both forms build and run on x86_64, PE and aarch64:
- the core consumer: all 7 dependency calls answer;
- the ai consumer: its POST gets a 201 (6 of 6 at the stub).

One PE warning comes with `[lib.ai]`, from the vendored stdlib: `lib/dynlib.cyr:674: syscall 5 … is not routed on
CYRIUS_TARGET_WIN=1`. The gates already filter `^warning:lib/`, and it is noted upstream.

| Consumer (DCE) | x86_64 | PE | aarch64 |
|---|---|---|---|
| `[lib]` | 786,920 B | 3,457,024 B | 3,588,240 B |
| `[lib.ai]` | 1,617,096 B | 5,151,744 B | 5,290,744 B |

### Probe 11: JSON (`p11_json.cyr`, `p11_json.sh`, `probe_json`)

serde_json 1.0.154 wrote 1,148 documents over the 28 types:
- 131 base documents: every preset, every cycle, grids, networks, steam entries, both configs and every enum variant;
- 1,017 variants of one value per type, in 71 kinds: malformed, duplicate, unknown and missing fields, numbers,
  strings, shapes and enum forms.

Each document carries serde_json's verdict and the value it read. bayan 1.11.1's `bayan_json_v_parse_buf` read every
one, and its tree was flattened into the oracle's token stream: keys in order, decoded strings, every float's bits.

**Where both read a document, bayan's tree is serde_json's**, with one exception: bayan reads the integer `-0` as +0,
where serde_json reads −0.0 into an f64 (20 documents). Their verdicts on what is JSON differ in six cases:

| Document | serde_json | bayan |
|---|---|---|
| a duplicate key in an unknown field, or nested inside one | ignored | refused: `duplicate key` |
| nesting past 128 levels inside an unknown field | read | refused: `nesting too deep` |
| an integer ≥ 2^63 (an f64 field, or a usize up to u64::MAX) | read | refused: `integer is out of range for i64` |
| a raw control character in a string | refused | parsed |
| invalid UTF-8 in a string | refused (typed) | parsed |
| `1e400` | refused: `number out of range` | parsed as +inf |

Everything else agrees:
- **Both refuse:** a truncated document, trailing data, two documents, a BOM, trailing commas, leading zeros, `+1.5`,
  `.5`, `1.`, NaN, Infinity, bad escapes and lone surrogates.
- **Both read:** whitespace, unknown fields, 100-deep nesting, exponent forms, underflow to 0 and `-0.0`.
- **The semantic cases parse in bayan,** and the decoders apply serde's rules to them: missing fields, `null`, wrong
  types, enum and array shapes, and a struct written as an array (serde's sequence form, which serde accepts).

**Floats** (160,378 doubles: every f64 in the base documents, every power of two, the powers of ten and their
neighbours, the printers' form thresholds, small integers, 100,000 random bit patterns and 50,000 short decimals):

| | Differs on |
|---|---|
| impetus's `fmtf.cyr` (sha256 `1d925be3…`), serde_json form vs serde_json's text | 0 |
| `fmtf.cyr`, Display form vs Rust's `{}` | 0 |
| bayan's own printer (Grisu2) vs serde_json's text | 6,276 (1,402 of them longer than shortest) |
| `bayan_f64_parse(serde_json's text)` vs the bits written | 0 |
| serde_json's own default parse of its text vs the bits written | **35,970 of 160,320**, by 1–2 ulp |

serde_json's default float parse is best-effort, so the oracle crate enables `float_roundtrip`, as impetus's does. All
160,320 texts and all 130 encodable base documents then read back exactly. Writing is the same in both modes.

**`#derive(Serialize)`** (P7's alternative), on StatePoint, ThermalMaterial and Phase:
- **It accepts what serde refuses** (49 of its 51 disagreements): truncated documents, trailing data, a BOM, trailing
  commas, duplicate fields, missing fields (read as 0), `{}`, `null`, NaN and Infinity (read as NaN), `1e400`, leading
  zeros, `+1.5`, `.5`, `1.`, and a string or a bool for an f64 (0).
- **It refuses** serde's sequence form (`[300.0,…]` returns 0).
- **Its f64 text is bayan's printer:** `9223372036854776000.0`, where serde writes `9.223372036854776e+18`.
- **It faults on missing strings.** A document without a `Str` field leaves the field 0, and the derived writer then
  faults on it.
- **Enum names are the constants'.** An enum writes `"US_PHASE_SOLID"` where serde writes `"Solid"`, and it can't read
  `"Solid"`. A variant name over 31 bytes is a compile error (`#derive field name too long`).

So P7 stands: hand-written codecs on bayan's tree, with the ADR 0005 amendment's rules. The summaries are identical on
x86_64, PE and aarch64.

## 4. M1 (2026-11-11 to 2026-11-12)

The harness's summary at the end of M1, identical on x86_64 ELF, `--bounds`, PE (wine) and aarch64 (qemu):
`8 modules, 0 failing | rows: 59803 ok, 0 failed, 0 measured, 0 pending, 39 accepted (ADR 4: 34)`.

- **fmtf** (M1a): 26,088 rows, every one ok. The vendored printer writes Rust's `{}` and serde_json's text for 13,044
  doubles. `validate.py` checks each row independently (read-back, shortest digits, closest, form).
- **material** (M1b): 2,927 rows.
  - The presets: 11.
  - `diffusivity` and `volumetric_heat_capacity`: 1,410 each, 0 ulp.
  - `to_json`: 19 documents, byte for byte.
  - `from_json`: 77 documents. 72 are ok, and 5 are accepted under ADR 5, which are exactly its mapped divergences:
    - a duplicate key inside an unknown field, and another nested in one: refused, where serde_json reads them;
    - an unknown field nested 200 deep: refused (bayan's cap is 128);
    - an integer ≥ 2^63 for an f64: refused;
    - an integer `-0`: read as +0, where serde_json gives −0.0.

    The harness accepts such a row only when bayan's own error is one of the three mapped refusals, or the document
    holds an integer `-0`.
- **transfer, part 1** (M1c): 7,217 rows over 18 functions.
  - **Class B** (conduction to boundary layers): every row bit-identical, the 233 refusals included (each the same
    Rust variant).
  - **Class C** (lumped capacitance, Dittus–Boelter, Churchill–Chu, the thermal boundary layer): 0 ulp against the
    corrected oracle. glibc 2.44's `pow` misrounds inside 4 rows, 2 in Dittus–Boelter and 2 in Churchill–Chu
    (`transfer.cr.tsv`, each call listed). The port matches MPFR there, and the oracle's own results are 1–2 ulp off.
  - **`validate.py`'s models** recompute every transfer row in Rust's order. With glibc's primitives all 7,217 match;
    with MPFR's, they write the corrections.
- **transfer, part 2** (M1d): 7,883 rows over 14 functions (fins, LMTD, ε-NTU, view factors), so 15,100 in all.
  - **The grid** adds each function's branch switches to the edges: LMTD at ΔT₁ = ΔT₂, 1 and 2 ulps apart, and at
    |ΔT₁/ΔT₂ − 1| = 1e-10 and its neighbours; ε-NTU over NTU 0…100 × C_r 0…1, with C_r at 1e-30 and 1 − 1e-10 and
    their neighbours; the parallel-plate factor across max(X, Y) = 0.01; fin efficiency at m·L = 1e-30; the
    perpendicular factor at aspect ratios where its clamp acts.
  - **Class B** (fin parameter, fin effectiveness, the two heat rates, NTU, coaxial disks): every row bit-identical,
    the 99 refusals included.
  - **Class C** (the fin heat and efficiency, both LMTDs, both effectivenesses, the parallel and perpendicular
    factors): 0 ulp against the corrected oracle, the 209 refusals included. glibc misrounds inside 4 rows, 3 `exp`
    (two in parallel-flow effectiveness, one in counter-flow) and 1 `atan` (the parallel-plate factor), so
    `transfer.cr.tsv` has 8 rows. The port matches MPFR there.
  - **`validate.py`** models all 32 functions: 15,100 rows reproduced with glibc's primitives.
- **The review** (M1d) changed two things, and every row above is unchanged by them:
  - `_us_finite_gt` / `_us_finite_ge` are one branch-free expression, `(x > bound) & (x < +∞)`, which `#inline` takes
    (three checks: 21 → 6 ns);
  - Rust's `powi(2)` and `powi(4)` are written as the products `__powidf2` forms, x·x and (x·x)·(x·x).

  The numbers are in `docs/benchmarks-rust-v-cyrius.md`.
- **Surface parity** (`scripts/surface_parity.py`): every public item of error (19), logging (2), material (16) and
  transfer (34) has its Cyrius name, and the 7 public constants of later modules have theirs. A missing name fails
  it, as checked by renaming one.
- **Tests** (`scripts/rust_test_ledger.py`, `docs/development/rust-test-ledger.md`): all 96 unit tests of error,
  validate, logging, material and transfer are ported, and the 6 integration tests that call only transfer and
  material (`tests/integration.tcyr`). The other 377 are counted under their milestones.

## 5. M2 (2026-11-12)

The harness's summary at the end of M2 is `9 modules, 0 failing | rows: 74377 ok, 0 failed, 0 measured, 0 pending, 46
accepted (ADR 4: 36)`, identical on x86_64 ELF, `--bounds`, PE and aarch64.

- **state, part 1** (M2a): 4,250 rows.
  - The presets and `ALL_GASES`: 9.
  - **Class B** (the ideal gas law three ways, van der Waals, isobaric work, Z, the reduced coordinates): every row
    bit-identical, refusals included. Van der Waals runs V across nb ± 3 ulps, so the V − nb > 0 refusal switches
    on the row where Rust's does.
  - **Class C** (isothermal work's ln, the adiabatic temperature's pow): 0 ulp against MPFR, and glibc misrounds
    none of these rows, so there is no `state.cr.tsv`.
  - **`GasData`'s JSON:** 18 `to_json` documents byte for byte; 78 `from_json` documents, 73 ok and 5 accepted under
    ADR 5, the same mapped divergences as material's (two duplicate keys in an unknown field, deep nesting, an
    integer `-0`, an integer ≥ 2^63).
  - **`validate.py`** reproduces all 4,145 function rows with glibc's primitives.
- **state, part 2** (M2b): 6,886 more rows, 11,136 in all.
  - **The cubic equations:** each gas at 11 temperatures from 0.3 T_c to 10 T_c (0.99, 1.0, 1.01 T_c among them) ×
    molar volumes from the co-volume b ± 2 ulps to 1 m³/mol, so the V_m > b refusal switches on the row where Rust's
    does. Redlich–Kwong is class C (T_c^2.5), Peng–Robinson class B (√ only): all 3,350 rows of the two pass, and so
    do their parameter rows (Redlich–Kwong's pair, Peng–Robinson's triple).
  - **Pitzer–Abbott Z** over T_r 0.3…5 × P_r 0…10 × five ω, out to where the correlation gives Z ≤ 0 and refuses:
    0 ulp against MPFR. glibc's `pow` misrounds T_r^4.2 inside one row, the first entry of `state.cr.tsv`; the port
    matches MPFR there.
  - **Kay's rules:** 233 mixtures × 3 rules, class B, every row bit-identical, the refusals included (unequal and
    empty lists, NaN, negative and out-of-range fractions, sums across the 1e-10 tolerance, non-finite values).
  - **Surface parity:** state 36 of 36.
- **entropy** (M2c): 3,454 rows over 9 functions (class A 115, B 1,970, C 1,369; 661 refusals).
  - **Class B and A** (heat transfer's entropy, both Carnot limits, Helmholtz, Gibbs, the Clausius check):
    bit-identical. The Carnot rows put T_cold at T_hot and 1–2 ulps either side at 300 K, 1 K, 1e-300 and 1e300, so
    the T_cold < T_hot refusal switches where Rust's does.
  - **Class C** (the ideal-gas and isothermal changes, mixing): 0 ulp against MPFR. glibc's `log` misrounds inside
    one ideal-gas row, `entropy.cr.tsv`'s only entry; the port matches MPFR there. Mixing runs 344 lists, subnormal
    fractions and sums across the 1e-10 tolerance among them.
  - **Surface parity:** entropy 9 of 9. **Tests:** all 23 of entropy's and the 2 entropy integration tests; every
    test due at M2 is ported (ledger: 180 of 479).

## 6. M3 (2026-11-12)

The harness's summary at the end of M3 is `10 modules, 0 failing | rows: 132784 ok, 0 failed, 0 measured, 0 pending,
76 accepted (ADR 4: 46)`, identical on x86_64 ELF, `--bounds`, PE and aarch64.

- **phase** (M3a): 27,333 rows (class A 232, B 2,389, C 24,712; 6,054 refusals).
  - **The substances:** the 8 presets and 106 edge substances: each of Water's and CO₂'s ten properties at 0, −1,
    NaN, ±∞, a boiling point above the critical point or below the triple point, a sublimer whose boiling point is
    above or at its triple point, and a triple point at and above 1 atm.
  - **The curves** (vapor and sublimation pressure, the sublimation latent heat): each substance at temperatures
    across every join, the triple, melting, boiling and critical points ± 1 ulp, and past both ends. 0 ulp against
    MPFR; glibc misrounds none of these rows.
  - **`phase_at`:** 10,797 points. At every join temperature each preset is asked at the pressure Rust's own curve
    gives there and 1 ulp either side, so each boundary decision (solid or gas, liquid or gas, supercritical) is
    tested on the row where it flips. Every row matches.
  - **The heats:** `heat_for_phase_change` runs every substance over ranges that start, end and cross on the
    melting and boiling points ± 1 ulp, both ways.
  - **JSON:** `Phase` as serde writes it (a string) and reads it (a string, or `{"Variant":null}`): 26 rows.
    `SubstanceData`: 15 `to_json` and 68 `from_json` documents, 5 accepted under ADR 5's mapped divergences.
  - **A predicate found in porting:** `vapor_pressure` checks T against the substance's own triple point, which a
    caller can set to −∞. The fast `_us_finite_ge` answers 1 for x = bound = −∞ where Rust's `finite_ge` answers
    false, so a data-dependent bound uses `_us_finite_ge_any`, exact for every bound.
  - **Surface parity:** phase 29 of 29, its enum members included.
- **steam, Region 4 and the saturated table** (M3b): 14,025 rows (class B 13,880, A 145; 2,478 refusals).
  - **Region 4:** the pressure equation over 273.15 … 647.096 K (1,500 seeded points, every table row's temperature
    ± 1 ulp, the limits ± 2 ulps) and the backward equation at each resulting pressure ± 1 ulp, at the range's ends ± 2
    ulps and at every row's pressure. Bit-identical: it is algebraic with √.
  - **The lookups:** by temperature and by pressure at every table row ± 1 ulp and midway between rows, seeded points,
    and the edges. Bit-identical, hisab's `calc_lerp` included.
  - **Quality:** each of 74 saturated states (the critical row among them, where v_g = v_f and h_fg = 0 refuse) at
    its liquid and vapor values, midway, past them by the 1e-9 round-off ± 2 ulps and at x = 1 + 2e-16. Wet steam at
    x across [0, 1] and just outside. Bit-identical.
  - **JSON:** `SaturatedEntry` and `WetSteamProperties`, 145 rows, 10 accepted under ADR 5's mapped divergences.
- **steam, superheated** (M3c): 17,088 rows more, 31,113 in all.
  - **The lookup** (10,999 rows) at every tier pressure ± 1 ulp, between tiers (the geometric midpoint and a
    quarter of the way), and the range's ends; at every table temperature ± 1 ulp, the midpoints between rows, and
    around each pressure's saturation temperature, where the tier's saturated-vapor anchor takes over. 0 ulp against
    MPFR; glibc misrounds a ln or an exp inside 6 rows.
  - **`superheated_by_entropy`** (6,037 rows) at each of those pressures over the entropy of every breakpoint ± 1 ulp,
    the midpoints, and below s_g (wet) and above the range. 0 ulp against MPFR; 4 rows inherit a glibc misround from
    an intermediate lookup. The breakpoints are sorted by insertion, the same order as Rust's `sort_unstable_by
    (total_cmp)` for these finite positive temperatures.
  - **`steam.cr.tsv`:** 10 rows (log 5, exp 5). **Surface parity:** steam 19 of 19. **Tests:** all 46 of steam's and
    M3's 5 integration tests (ledger: 275 of 479, every one due).
