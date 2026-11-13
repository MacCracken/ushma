# scripts/oracle — ushma's golden-value generator

Drives the frozen Rust ushma 1.4.0 in `rust-old/` (by path, read-only; it is never edited)
and writes everything the Cyrius port derives from it:
- the golden corpus, `tests/golden/<module>.tsv`, which the parity harness replays;
- the generated Cyrius sources, `src/constants.cyr` and `src/<module>_data.cyr`.

Adapted from kimiya's `scripts/oracle/` (format and rules unchanged; tanmatra's originally).

```sh
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml             # write
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check  # verify
python3 -I -B scripts/oracle/validate.py .                                          # independent check
```

`--check` regenerates every output in memory and compares it byte for byte. It exits 1 on:
- a difference or a missing file;
- a `.tsv` with no emitter;
- a corpus recorded on another glibc platform.

`Cargo.lock` is committed. It was seeded from `rust-old/Cargo.lock`, so every shared
dependency, hisab 1.4.0 above all, is the version the oracle was tested with.

## What's generated (M0 parts 2a and 2b, M1)

| Output | Rows / lines | What it is |
|---|---|---|
| `tests/golden/constants.tsv` | 9 | every public scalar const of ushma 1.4.0 (8 `f64`, 1 `usize`), class A |
| `tests/golden/fmtf.tsv` | 26,088 | Rust's `{}` and `serde_json::to_string` of 13,044 doubles (`m_fmtf.rs`: the specials, every power of two, the powers of ten ± 1 ulp, the form thresholds, 0 … 1000 in tenths and hundredths, 4,000 random bit patterns, 2,000 short decimals), class A; `validate.py` checks each text by its properties |
| `tests/golden/steam.tsv` | 31,113 | M3 (`m_steam.rs`): Region 4 both ways, both saturated lookups at every row ± 1 ulp and between rows, quality and wet steam around their tolerances (class B); the superheated lookup at every tier ± 1 ulp and between tiers, and `superheated_by_entropy` over every breakpoint's entropy (class C; `steam.cr.tsv`: 10 rows); the three steam types' JSON (class A). The tables are private: their rows come from `t_steam::saturated_rows` / `superheated_rows` |
| `tests/golden/cycle.tsv` | 5,607 | M4 (`m_cycle.rs`): the process helpers (class B; `adiabatic_temp_from_pressure` class C); Otto, Diesel and Brayton over r, γ and T grids, with Diesel's r_c at r ± 1 ulp and Brayton's T₃ at T₂ ± 1 ulp, each row the whole result; Rankine (its `Option<f64>` superheat as `n:` or `f:`) and refrigeration over the steam tables, each branch switch tried either side (the exhaust's and the compressor exit's found by bisection on the bit patterns); `heat_pump_cop` and `compare_cycles` (class B); the T-s and P-v diagrams on every cycle's results and on results edited to reach each path's guard (a result and `i:points_per_process` in, `v:count` and the points out). Class C for the cycles and diagrams (`cycle.cr.tsv`: 3 rows). JSON for the seven cycle types (class A) |
| `tests/golden/transfer.tsv` | 15,100 | M1c–M1d (`m_transfer.rs`): each of the 32 public transfer functions over a seeded physical sweep, every argument at each edge, and each branch switch (LMTD at ΔT₁ = ΔT₂ and across the 1e-10 band, NTU 0…100 × C_r 0…1 with both switches, view factors across 0.01); class B, or C where `exp`, `ln`, `powf`, `tanh` or `atan` enters (`transfer.cr.tsv`: 8 rows) |
| `tests/golden/math.tsv` | 15,687 | every std f64 method ushma's library calls, over ushma's ranges (port plan §7.1) |
| `tests/golden/entropy.tsv` | 3,454 | M2c (`m_entropy.rs`): the 9 entropy functions over seeded sweeps, every edge and each switch (T₁ = T₂, V₁ = V₂, T_cold at and around T_hot, 344 mixtures); classes A, B and C (`entropy.cr.tsv`: 1 row) |
| `tests/golden/phase.tsv` | 27,333 | the 8 presets and `ALL_SUBSTANCES` (part 2b), then M3a (`m_phase.rs`): the 10 phase functions over the presets and 106 edge substances across every join, `phase_at` at Rust's own curve pressures ± 1 ulp (classes A, B, C); 4 + 22 `Phase` and 15 + 68 `SubstanceData` JSON rows (class A) |
| `tests/golden/state.tsv` | 11,136 | the 8 presets and `ALL_GASES` (part 2b), then M2 (`m_state.rs`): all 21 state functions over seeded sweeps, every edge and each switch (the cubic equations: every gas × T across T_c × V_m across b; Pitzer–Abbott: T_r × P_r to Z ≤ 0; Kay's rules: 233 mixtures), classes B and C (`state.cr.tsv`: 1 row); 18 `GasData::to_json` and 78 `GasData::from_json` documents (class A) |
| `tests/golden/material.tsv` | 2,927 | M1b (`m_material.rs`): the presets, then 1,410 rows per method (`diffusivity`, `volumetric_heat_capacity`, class B), 19 `to_json` and 77 `from_json` documents (class A; the variants from `json_cases.rs`, shared with probe 11) |
| `tests/golden/{state,phase,chem}.tsv` | 9 / 9 / 13 | every public table preset (name + fields) and its `ALL_*` order, class A; each module's function rows join it at its milestone |
| `tests/golden/chem.tsv` | 2,340 | M5a (`m_chem.rs`), after the 13 preset rows: Hess's law and the Gibbs sum over named, reversed, empty and seeded reactions (a side is `v:n` then per term its coefficient and the species by value), each coefficient and species at the edges (class B); ln K and ΔH − TΔS (class B); K and Van 't Hoff, with where exp leaves the f64 range bisected ± 2 ulps (class C; `chem.cr.tsv`: 3 rows); the flame temperature over named and seeded combustions, the edges, and where T_flame crosses 0 (class B); `Species`'s JSON (class A) |
| `tests/golden/stat.tsv` | 4,170 | M5b (`m_stat.rs`): every stat function seeded and at the edges, and every switch either side: Einstein's x = 500 ± 1 ulp and where (eˣ − 1)² < 1e-30 (bisected) ± 2 ulps; Debye's Θ/T at 0.01 and 50 ± 1 ulp, a sweep and 1.4.0's reference cases; the vibrational partition function where hν ≪ kT and where its exponent underflows; σ (Rust's `u32`) as `i:` tokens. Class B (+ − × ÷ √), C (exp, expm1, ln, pow; `stat.cr.tsv`: 11 rows), D for `debye_cv` (hisab's Simpson sums in another order; probe 1) |
| `tests/golden/numerical.tsv` | 2,211 | M6 (`m_numerical.rs`): `ThermalGrid1D::new` over its arguments and every condition at its edges; `fourier_number`, `max_stable_dt`; the three 1D steps on seven setups, each row one step from the oracle's own state (a grid is `v:n` and its nodes, Δx, α and its two conditions, each its variant token and values), along trajectories and at the explicit stability limit, the adaptive exact multiples and sub-step cap either side (bisected), dt at the edges and malformed grids; `ThermalGrid2D::new`, `set_boundary` and `solve_steady_state` (a 2D grid is its rows, `v:ny` then each `v:len`, then nx, ny, Δx, Δy and four conditions), the solve at its converging count and one fewer; `ThermalNetwork` (a network is num_nodes, `v:` resistances a, b, R and `v:` fixed nodes node, T): `new`, `add_resistance` and `set_fixed_temperature` at every index and value edge, `solve` on Rust's test networks, garbage lists, every structure edge (self-loops, repeated edges, every node fixed, −0 K), the overflow cases (ADR 0002 row 16) and 400 seeded networks as probe 9 draws them; the JSON of `BoundaryCondition`, `Side`, both grids and the network. Class B, A for JSON and the network's edits |
| `tests/golden/bridge.tsv` | 6,439 | M7a (`m_bridge.rs`): the 16 bridge conversions, each seeded over its physical range, each argument of its base cases (the 1.4.0 tests' values) at every edge, and each fallback's switch either side (m·c_p at 0, K = 1 ± 12 ulps across the |ln K| < 1e-15 test, T reaching 0 with altitude, η at 1 ± 2 ulps, the activation at 0 and 1, the wet bulb's clamp near saturation). Class B (+ − × ÷ √, powi), C where ln, powf or atan enters (`bridge.cr.tsv`: 3 rows) |
| `tests/golden/modules.txt` | 14 | the module list the harness reads |
| `src/constants.cyr` | 47 lines | the 9 public constants as `US_` consts, under their Rust doc comments |
| `src/{material,state,phase,chem}_data.cyr` | 101–202 lines | the public tables: `US_<TABLE>_COUNT`, one `US_<TABLE>_<PRESET>` index constant per preset, one `f64[n]` column per field, the names as a NUL-separated blob with `u16` offsets (`t_tables.rs`) |
| `src/steam_data.cyr` | 877 lines | steam's private IAPWS-IF97 data: `R4_N`, the Region 4 scalars, the 42 × 9 saturated and 72 × 5 superheated tables, the 9 tier pressures (`t_steam.rs`) |
| `src/constants_private.cyr` | 41 lines | the 8 other private constants, fn-local ones included (bridge, cycle, numerical, stat; `t_private.rs`) |

**Private data** isn't reachable through the API. It is read from its literal in
`rust-old/src` (`rustsrc.rs`, Rust's `f64::from_str`; `KLEIBER_W`'s `70.0 * 4184.0 / 86_400.0`
folded as Rust's const evaluation folds it), never re-typed. `t_private::complete` proves that
every `const` item in every module is generated by exactly one emitter or excused with a reason.
The one excuse is steam's `const _` assertion, whose property `t_steam` re-checks. `t_steam` also
asserts what the steam code relies on: each entry carries every field in declaration order, and
each superheated row's pressure is its tier's.

Added at their milestones (port plan §4 step 6, §6):
- each module's rows (M1 on);
- the JSON rows (each module's milestone);
- the request rows (M7);
- the corrected-oracle files `<module>.cr.tsv` (ADR 0004, M0 part 2c).

## Format 1

One row per call, five tab-separated columns:

```text
fn <TAB> class <TAB> inputs <TAB> outputs <TAB> note
```

- **fn** — the Rust function, e.g. `f64::powf` or `transfer::conduction`; `const` for a constant.
- **class** — the parity class (port plan §5):
  - `A` exact (NaN matches any NaN);
  - `B` arithmetic (≤ 2 ulp, expected 0);
  - `C` transcendental (0 ulp against the corrected oracle, ADR 0004);
  - `D` iterative (relative ≤ 1e-12).
- **inputs**, **outputs** — space-separated typed tokens (`-` when empty):

  | Token | Meaning |
  |---|---|
  | `f:<16 hex>` | an f64 by its IEEE-754 bit pattern, lowercase |
  | `i:<decimal>` | an integer; booleans are 0 / 1 |
  | `e:<Enum>.<Variant>` | a Rust enum variant by name |
  | `s:<percent-encoded UTF-8>` | a string; every byte outside printable ASCII, and `%`, is `%XX` |
  | `n:` | `Option::None` |
  | `E:<Variant>` | `Err(UshmaError::<Variant>)`, the variant name up to its payload |
  | `v:<len>` | a sequence; its elements' tokens follow in order |

- **note** — a readable rendering of the call, for review. The harness ignores it, except for
  the ADR 0002 notes that end `(<tokens>)]` and name the port's value.

Each file starts with five `#` lines: the format and module, the oracle provenance (ushma and
hisab versions; `libm: glibc <version> <arch> <fma|no-fma>`), the regenerate command, the
column legend, and the row count.

## Where the math rows come from

The release-built oracle binary's libm imports (`nm -D target/release/oracle`, rustc 1.99.1,
x86_64): `atan@GLIBC_2.2.5`, `exp@GLIBC_2.29`, `expm1@GLIBC_2.2.5`, `log@GLIBC_2.29`,
`pow@GLIBC_2.29`, `tanh@GLIBC_2.2.5`. `powi` is the statically linked compiler-builtins
`__powidf2`, and `sqrt`, `abs`, `ceil`, `min`, `max`, `clamp`, `is_finite` and `is_nan` are
inline. So the class C rows are glibc 2.44's bits, and `--check` refuses another platform.

## validate.py

It shares no code with the generator. It re-derives from first principles:
- the format of every file, and that `modules.txt` lists exactly the module files;
- every constant, against its Rust literal parsed by Python's correctly rounded `float()`;
- every `src/constants.cyr` literal, against its golden bits;
- every class A/B math row, recomputed in binary64 (`powi` in `__powidf2`'s order, `min`/`max`
  with Rust 1.99's recorded semantics: a NaN operand is ignored, and (±0, ∓0) gives the second
  operand).

Since part 2b it also checks every public table row (name and fields) against the preset's
struct literal in the Rust source, and the `ALL_*` row against the slice's entries. Class C rows
are checked against MPFR by the corrected-oracle tooling (M0 part 2c). Run it with
`python3 -I -B`. A tampered row fails it, as checked when it was written.

From M1 it also models each ported module's functions (`MODELS`: transfer's 32, state's 21, entropy's 9, phase's 10, steam's 10, cycle's 13, chem's 7, stat's 14 with Rust hisab 1.4.0's Simpson for `debye_cv`, numerical's 9 so far; the comparison, diagrams and term-list functions read their raw tokens, `RAW_ARGS`): each is recomputed in Rust's
operation order with Python binary64 arithmetic, its transcendentals (`exp`, `ln`, `pow`, `tanh`, `atan`) taken from
glibc through ctypes, and every row must come out token for token. With MPFR's primitives instead, the same models
give the module's corrected oracle (`<module>.cr.tsv`), each misrounded call recorded.

## mpfr_ref.py and the corrected oracle (ADR 0004)

`mpfr_ref.py` calls the system MPFR (4.x, through ctypes): 53 bits, binary64's exponent range,
subnormalized. That makes every result the IEEE round-to-nearest-even value of the exact
function, exact ties and subnormals included. `python3 -I -B scripts/oracle/mpfr_ref.py` runs its
self-checks.

`validate.py . --write-corrections` writes `tests/golden/math.cr.tsv`: one row per class-C math
row whose correctly rounded output differs from the oracle's (glibc's):

```text
fn <TAB> inputs <TAB> corrected outputs <TAB> note
```

- **fn** and **inputs** are the golden row's, byte for byte (the key the harness matches).
- **note** is `oracle <outputs>; <fn>(<args>) glibc f:<g> correctly rounded f:<c>`.
- The header is the golden file's provenance plus the MPFR version.

The plain `validate.py .` run recomputes the file and checks it. On the recorded libm platform
it also checks that each note's glibc value is this machine's libm result, and the exact text.
`math.cr.tsv` has 26 rows (exp_m1 23, ln 1, powf 2), `transfer.cr.tsv` 8 (pow 4, exp 3, atan 1),
`state.cr.tsv` 1 (pow), `entropy.cr.tsv` 1 (log), `steam.cr.tsv` 10 (log 5, exp 5), `cycle.cr.tsv` 3 (log 1, pow 2), `chem.cr.tsv` 3 (exp) and `stat.cr.tsv` 11 (expm1 10, exp 1).

## mathx_tables.py

Generates `src/mathx_data.cyr`, the tables and constants of `src/mathx.cyr`'s correctly rounded exp, expm1, ln,
pow and atan2. It is vendored from kimiya's generator with the same algorithms (MPFR at 400 bits, exact
`Fraction`s) and the same assertions on what each error analysis relies on, emitting only the 46 items ushma's
closure reads. `--check` verifies the file is current. `--check-kimiya` also confirms every value is bit-identical
to kimiya's `src/mathx_data.cyr` (the gates run it when kimiya's checkout is beside this one).

## cr_verify.py

Checks `src/mathx.cyr`'s correctly rounded exp, expm1, ln, pow, atan2 and tanh against MPFR, bit for bit, on
x86_64 ELF, PE (wine) and aarch64 (qemu). Each primitive gets about 10^6 seeded arguments:
- binade-uniform over its domain, plus the dense regions where the code changes path;
- the known hard-to-round cases;
- the accurate path forced on 100,000 of them;
- a check that the fast path agrees with the full result wherever its test passes, and that every target gives the
  same bits and flags.

The families for the five vendored primitives are kimiya's, unchanged, and tanh's is ushma's own. It also checks
x^2.5 over every binade, and the special values of pow, atan2 and tanh against MPFR and glibc. Options:
`--only tanh,exp`, `--targets elf`, `--work DIR` (default `build/cr-verify`).

## probe_tables.py

Compiles every generated constant and cell with the Cyrius compiler (`build/oracle-probe/`)
and checks the bits:
- **Public data:** against the golden rows (the oracle's API).
- **Private data:** against the script's own Python parse of the Rust literals, independent
  of the Rust generator.

It fails on:
- a mismatch;
- a `var` or `const` it didn't probe;
- a compiler diagnostic other than the vendored ones recorded in `state.md`.

`--all-targets` also builds the probe for aarch64 (qemu) and PE (wine) and requires the same
output. 1,052 values at part 2b. Tampering with a cell, a name byte, an index constant or a
private constant, or adding an unprobed array, fails it, as checked when it was written.

## The probe binaries (M0 part 3b)

Three more binaries of this crate give the M0 probes their Rust halves. Each drives the frozen 1.4.0, and Rust hisab
1.4.0 itself, read-only, and writes under `build/` (gitignored). Nothing they write is golden data:

| Binary | Probe | Writes |
|---|---|---|
| `probe_hisab` | 1, 9 | `build/probe-hisab/simpson.txt` (3,000 `integral_simpson` rows) and `gauss.txt` (2,000 network systems with 1.4.0's solutions, each mirrored bit for bit against `ThermalNetwork::solve`, plus `gref` rows for the systems 1.4.0 refuses) |
| `probe_ai` | 7 | nothing: it runs `DaimonClient::register_agent` once against the loopback stub (`scripts/probes/p7_daimon_stub.py`), which records the request |
| `probe_json` | 11 | `build/probe-json/cases.txt` (1,148 serde_json documents over the 28 types, with serde's verdicts, decoded values and token streams) and `f64.txt` (160,378 doubles: serde_json's text, Rust's Display, serde's read-back) |

```sh
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml --bin probe_hisab
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml --bin probe_json
```

For these, the crate also depends on `hisab`, `serde`, `serde_json` and `tokio`, each pinned to the version the
lockfile already resolved, so no package was added. serde_json is built with **`float_roundtrip`**: its default
float parse is best-effort, and reads 22 % of its own shortest texts back 1–2 ulp off (probe 11). The parse precision
belongs to the consumer, not to ushma, which only derives `Deserialize`, and the port's parse is correctly rounded
(bayan). Writing is identical in both modes, and nothing in the golden corpus parses JSON. tokio's `rt`, `net` and
`time` features add no crate; `probe_ai` uses a current-thread runtime.
