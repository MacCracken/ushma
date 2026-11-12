# ushma Rust → Cyrius port plan

> **Plan of record** for porting ushma to Cyrius. Rust **1.4.0** is the final Rust feature
> release and the **reference implementation** (the oracle); the Cyrius port ships as **2.0.0**.
> The method is tanmatra's ([its ADR 0006](../../../tanmatra/docs/adr/0006-cyrius-port-with-rust-oracle.md),
> [its port plan](../../../tanmatra/docs/development/tanmatra-port-plan.md)), as adapted by
> kimiya ([its port plan](../../../kimiya/docs/development/kimiya-port-plan.md)) and impetus
> ([its port plan](../../../impetus/docs/development/impetus-port-plan.md)). Where this plan says
> "as tanmatra", "as kimiya" or "as impetus", follow that document. §2 lists what is different
> here and why.
>
> ▶ **Start on Cyrius 6.7.18 with hisab 3.4.1**, after 1.4.0 is tagged (M0 step 0).
>
> 🟡 **Status:** planned 2026-11-11. The proposals in §0.2 are confirmed when ADRs 0001–0003 are
> written (M0 step 2).

| Thing | Value | Source |
|---|---|---|
| Rust source (oracle) | ushma **1.4.0** at `9528e1f` (committed, **not yet tagged**): ≈5,100 code lines in `src/` (no tests, comments or blank lines), 155 public functions and methods, 18 files, 10 domain features plus `soorat-compat`, `ai`, `logging` | `Cargo.toml`, `VERSION`, `git tag` |
| Rust health, 2026-11-11 on rustc 1.99.1 | ✅ 449 unit + 30 integration tests (no doc tests); `cargo fmt --check`, clippy `--all-features --all-targets -D warnings`, rustdoc `-D warnings`, `cargo audit` (224 crates, clean), `cargo deny` (advisories, bans, licenses, sources ok) | local run |
| Correctness baseline | The 1.4.0 audit: IAPWS-IF97 steam tables, exact Region 4, real process paths in the diagrams, NaN/±∞ rejected by every guard, 80 regression tests with IF97, Incropera and NIST reference values | `CHANGELOG.md` 1.4.0 |
| Cyrius toolchain | **6.7.18** (installed, latest) | `~/.cyrius/current` |
| ganita | **1.2.15**, bundled in the 6.7.18 stdlib | `~/.cyrius/lib/ganita.cyr` |
| hisab | **3.4.1**, Cyrius, tag `3.4.1`, bundle `dist/hisab.cyr` | `/home/macro/Repos/hisab` |
| Oracle's math | Rust **std** → host **glibc 2.44**. Calls in the in-scope modules: `powf` (23), `ln` (22), `exp` (14), `atan` (7), `tanh` (2), `exp_m1` (2). Rust **hisab 1.4.0** (pinned in the lockfile): `calc::lerp` (3), `calc::integral_simpson` (1), `num::gaussian_elimination` (1) | `grep` of non-test `src/`, `Cargo.lock` |
| Dependents | Local: **soorat 1.0.0** (Rust), optional feature `thermo` → `ushma = "1"` with `soorat-compat`. kiran, joshua, badal and dravya have no local repo. crates.io reverse dependencies are checked at M0 step 0. **No Cyrius consumer yet.** | `~/Repos/*/Cargo.toml`, `cyrius.cyml` |
| Target | ushma **2.0.0**, Cyrius | this document |

## 0. Decisions

### 0.1 Carried from the sibling ports

The maintainer settled these during the tanmatra, kimiya and impetus ports. They apply here
unchanged.

- **The oracle is frozen.** Rust 1.4.0 moves to `rust-old/` and is never edited. A generator
  crate in `scripts/oracle/` drives it by path. A wrong oracle value is fixed in Rust as 1.4.x
  and the corpus is regenerated.
- **No Rust 2.0.** API shape changes land in the Cyrius API only, each accepted by an ADR row
  (§3.3). Rust stays on 1.x for soorat and any crates.io dependent.
- **hisab-first** (impetus ADR 0007). The port calls hisab for everything hisab provides. Where
  hisab or the stdlib falls short, the port files the gap upstream (§8) with the probe that
  shows it. A stopgap is local and named as one.
- **Accuracy over glibc parity** (kimiya ADRs 0005 and 0006, impetus ADR 0011). Every
  transcendental primitive is correctly rounded and gives the same bits on every target.
  "Parity" means functional parity: the same results, `Err` codes, check order and branches.
  It does not mean copying glibc's last-ulp errors. The class-C reference is the *corrected
  oracle*: Rust's operation sequence evaluated with correctly rounded primitives (§5).
- **Performance is measured, not gated.** Every in-scope Rust benchmark gets a Cyrius twin, and
  both columns go in `docs/benchmarks-rust-v-cyrius.md`. A slow row is a finding to file, not a
  reason to change the port.
- **New physics stays out.** The roadmap backlog (the 2D ghost-node walls, SOR, the Dittus–Boelter
  cooling exponent, compressed liquid, the sloped melting curve, …) comes after 2.0, demand-gated.
- **The maintainer does all git.** Each verified milestone part ends with a summary and a request
  to commit before the next part starts.

### 0.2 Proposed for ushma (confirm in M0 step 2)

| # | Proposal | Why | Alternative |
|---|---|---|---|
| P1 | **Scope.** Port the 10 domain modules plus `error` and `validate`. Defer `bridge` and `integration/soorat`. Drop `ai` and `logging`. | No Cyrius code calls the bridges. soorat, the only known dependent, is Rust and keeps using Rust 1.x. `ai` is a reqwest client that registers with daimon, not thermodynamics. `logging` wraps tracing-subscriber (§3.2). | Port `bridge` too: 16 thin functions with no dependencies, about one sub-milestone. |
| P2 | **Range-edge sweep inside each milestone**, not a library-wide audit afterwards. Each milestone runs kimiya's overflow/exact-constant checks on its own functions before its gate. | kimiya ran its audit after M4 and had to reopen every module (batches B1–B12). Doing it per module costs the same work once. | kimiya's order: port first, audit later. |
| P3 | **Correctly rounded primitives.** Vendor the subset of kimiya's `src/mathx.cyr` that ushma calls into `src/mathx.cyr` (`_us_cr_exp`, `_ln`, `_pow`, `_expm1`, and atan as atan2(x, 1)). Write a new correctly rounded `tanh`. File a request to host the shared set upstream. | ushma would be the **third** private copy, after impetus's `crmath.cyr` (1,994 lines) and kimiya's `mathx.cyr` (5,384 lines). CLAUDE.md's third-instance rule says extract it. Extracting it first would block ushma on kimiya's open audit batches B10 and B11, which still change `mathx.cyr`. | Extract to ganita or hisab first, then port. |
| P4 | **Dense Gaussian elimination** (`ThermalNetwork::solve`): transcribe Rust hisab 1.4.0's `gaussian_elimination` into `src/numerics.cyr` as a private `_us_*` stopgap (class B), and file a hisab request. | hisab 3.4.1 has `num_tridiag_solve` but no dense direct solver. kimiya's M5 needs the same routine (its ADR 0003), so the request has two consumers. | Wait for hisab before M6. |
| P5 | **One bundle.** `[lib]` holds every in-scope module, steam included. | Rust's 10 features are additive compile-time gates. Cyrius `[features]` only select optional dependencies, and profiles are for exclusive variants (impetus's 2D/3D). Rankine and refrigeration need steam anyway. | `[lib]` without steam plus a `[lib.steam]` profile (as tanmatra's optics). |
| P6 | **Public names get a `us_` prefix** (`us_ideal_gas_pressure`), with `US_` constants and `US_ERR_*` codes. | Cyrius has one global namespace. A duplicate non-private `fn` is only a warning, and **the last definition wins**. 10 ushma functions already collide with kimiya's Cyrius bundle: `adiabatic_flame_temperature`, `equilibrium_constant`, `ideal_gas_pressure`, `ideal_gas_volume`, `is_spontaneous`, `peng_robinson_pressure`, `reaction_enthalpy`, `redlich_kwong_pressure`, `van_der_waals_pressure`, `vant_hoff_k`. kiran uses both crates in Rust. | Keep Rust's names (as kimiya and tanmatra) and rename only the 10. |

---

## 1. Porting models and references

- **kimiya** (`/home/macro/Repos/kimiya`, in progress, M0–M4 done) is the **closest model**. Both
  are pure-function science libraries of similar size (kimiya 1.2.0: 6,824 lines, 202 functions,
  662 tests) with a glibc oracle. Copy and adapt:
  - `scripts/oracle/`: the generator crate, `validate.py` (and `--sensitivity`), `mpfr_ref.py`,
    `cr_verify.py`, `probe_tables.py` and `mathx_tables.py`;
  - `tests/parity.tcyr`, `tests/parity_core.cyr` and `tests/ulp.cyr`, with the ADR-mapping table
    and the "accepted under ADR 6" count;
  - `src/mathx.cyr` (P3);
  - ADR 0002's divergence format, one row per function family, with the maintainer's decisions
    quoted verbatim;
  - `handoff/gates/gates.sh`, the full three-target gate run. `handoff/` is untracked and will be
    deleted, so copy it into `scripts/gates.sh` at M0;
  - `handoff/HANDOFF.md`'s per-batch loop. An agent ports one or two families while independent
    Python models built from `rust-old/src` are written. Then a probe compares the port with the
    models bit for bit on all three targets, and every difference is investigated.
- **impetus 2.0.0** (`/home/macro/Repos/impetus`, done 2026-11-10) is the model for:
  - hisab-first (ADR 0007) and its gap list (plan §8);
  - **single-step parity rows** for stateful objects (plan §2.1, §5): ushma's thermal grids and
    network (§2.2);
  - the M0 **probes** (plan §4 step 9, `scripts/probes/`): memory, struct values, per-test
    module lists;
  - private constants extracted from the Rust source with a completeness proof
    (`scripts/oracle/src/t_private.rs`);
  - `scripts/bench_record.py` and `scripts/doc_examples.py`.
- **tanmatra 2.0.2** (`/home/macro/Repos/tanmatra`) is the original template. Use it for
  `scripts/rust_test_ledger.py`, `scripts/surface_parity.py`, the port-parity rule for class-C
  bounds, and the CLAUDE.md translation (its plan §2 step 2).
- **hisab 3.4.1** (`/home/macro/Repos/hisab`) is the dependency. Read its `CLAUDE.md`,
  `docs/guides/migration-3.4.md` (the Result convention) and the doc comment of every function
  ushma calls (§7.2).
- **Language and toolchain** (`/home/macro/Repos/cyrius`): `docs/guides/cyrius-guide.md`,
  `docs/stdlib-reference.md`, `docs/development/ecosystem-migration.md` (6.7.x entries),
  `cyrius help manifest`.

---

## 2. What is different from the siblings

1. **The steam tables are private and the largest data in the crate.** `SATURATED_TABLE` has 42
   rows × 9 fields; `SUPERHEATED_DATA` has 72 rows (9 pressure tiers × 8 temperatures); `R4_N`
   holds the 10 IF97 Region 4 coefficients. Their API exposes interpolated values, not the rows.
   They are extracted from `rust-old/src/steam.rs` with Rust's own `f64::from_str`, plus a proof
   that every `const` in the module is emitted or excused (impetus's `t_private.rs`). The golden
   rows that read them then prove them. They are never re-typed.
2. **Stateful objects, but with public state.** `ThermalGrid1D`, `ThermalGrid2D` and
   `ThermalNetwork` step or solve in place. Their fields are public, so the oracle can dump the
   full state without impetus's hidden `oracle` feature (confirm at M0 step 0). They get impetus's
   second parity layer: load the state, step once, compare the whole state. A 1-ulp difference
   doesn't compound here the way a stacked pile does (diffusion damps it), but each step is still
   compared, not a long trajectory.
3. **ushma's own iterations.** These are ushma's loops, transcribed exactly, expected class B (D
   only where an iteration count can legitimately differ):
   - Rankine's and refrigeration's s(T, P) = s₁ inversions (bisection over the table);
   - `superheated_by_entropy`'s piecewise-linear inversion;
   - Gauss–Seidel steady state;
   - adaptive sub-stepping;
   - the adiabatic flame temperature.
4. **Three Rust hisab 1.4.0 calls change library** (§7.2). `calc_lerp` keeps the fast form's
   bits (class B). `calc_integral_simpson` changed its framing in 3.3.6, so probe it (class B or
   D). `gaussian_elimination` has no hisab 3.4.1 counterpart (P4).
5. **Results are wide structs and point lists.** Each cycle function returns a struct of state
   points, heats and efficiencies; the diagrams return `Vec<DiagramPoint>`. The open
   small-struct `Result` payload bug
   (`cyrius/docs/development/issues/2026-11-08-small-struct-result-payload-read-as-address.md`)
   means these go through out-pointers or handles. Probe 3 checks this at the 6.7.18 pin.

---

## 3. Scope

### 3.1 In scope for 2.0.0

| Milestone | Rust module | Code lines | Public fns | Unit tests | Cyrius module(s) |
|---|---|---|---|---|---|
| M0/M1 | *(correctly rounded primitives, P3)* | — | — | — | `src/mathx.cyr` (private `_us_cr_*`) |
| M1 | `error.rs` | 37 | — | 18 | `src/error.cyr` |
| M1 | `validate.rs` | 19 | (private) | 2 | `src/validate.cyr` |
| M1 | *(8 public constants in 4 modules)* | — | (consts) | — | `src/constants.cyr` (generated) |
| M1 | `material.rs` | 106 | 2 | 8 | `src/material_data.cyr` (generated), `src/material.cyr` |
| M1 | `transfer.rs` | 443 | 32 | 67 | `src/transfer.cyr` |
| M2 | `state.rs` | 473 | 21 | 47 | `src/state_data.cyr` (generated), `src/state.cyr` |
| M2 | `entropy.rs` | 132 | 9 | 23 | `src/entropy.cyr` |
| M3 | `phase.rs` | 353 | 10 | 44 | `src/phase_data.cyr` (generated), `src/phase.cyr` |
| M3 | `steam.rs` | 1,320 | 10 | 46 | `src/steam_data.cyr` (generated, §2.1), `src/steam.cyr` |
| M4 | `cycle.rs` | 730 | 13 | 59 | `src/cycle.cyr` |
| M5 | `chem.rs` | 203 | 7 | 23 | `src/chem_data.cyr` (generated), `src/chem.cyr` |
| M5 | `stat.rs` | 214 | 14 | 33 | `src/stat.cyr` |
| M6 | *(hisab gap, P4)* | — | — | — | `src/numerics.cyr` (private `_us_*`) |
| M6 | `numerical.rs` | 599 | 13 | 27 | `src/numerical.cyr` |

These counts were measured on 2026-11-11 against 1.4.0. In scope: 131 of the 155 public
functions, 397 of the 449 unit tests, and all 30 integration tests (none touches a deferred
module).

### 3.2 Deferred or dropped (P1)

| Rust surface | Size | Status | Unblocked when |
|---|---|---|---|
| `bridge.rs`: bijli, kimiya, dravya, badal, pravash, sharira conversions | 168 lines, 16 fns, 31 tests | **Deferred** | A Cyrius consumer calls one. kimiya 2.0 deferred its own ushma bridges until ushma ports, so the two sets should port together |
| `integration/soorat.rs` (`soorat-compat`) | 213 lines, 5 fns, 13 tests | **Deferred** | soorat ports. Until then soorat's `thermo` feature uses Rust 1.x |
| `ai.rs` (`ai`: daimon registration over reqwest) | 86 lines, 1 fn, 7 tests | **Dropped** | A consumer asks (Cyrius would need an HTTP client) |
| `logging.rs` (`logging`: tracing-subscriber init) | 14 lines, 2 fns, 1 test | **Dropped** | — (the Cyrius analogue is sakshi) |

`rust-old/` keeps every deferred function and its tests as the reference.

### 3.3 API changes in the Cyrius version (accepted divergences, ADR 0002)

These are the starting list. Rows found during the port are added there, one per function
family, and anything that moves results noticeably goes to the maintainer first.

| # | Rust 1.4.0 | Cyrius 2.0.0 | Golden comparison |
|---|---|---|---|
| 1 | `UshmaError` variants with `String` payloads (`DivisionByZero`, `InvalidParameter`, `SubstanceNotFound`) | `US_ERR_*` code enum (append-only) + numeric detail, as kimiya | the code is class A; message text is not compared |
| 2 | `ThermalGrid2D` stores `Vec<Vec<f64>>`; Crank–Nicolson allocates per step (roadmap backlog: breaking in Rust) | flat row-major buffer + dims; caller-reusable scratch | the state dump is flattened the same way, class B |
| 3 | `#[tracing::instrument]` on 28 functions | no spans (as tanmatra, kimiya and impetus) | n/a |
| 4 | Lookups by `&str` name (`SubstanceNotFound`) | table index (−1 when absent) plus a name search | same index set, class A |
| — | Exponents written as rounded decimals (`powf(1.0/3.0)`, `powf(1.0/6.0)`) and any similar exact-constant miss | fixed under the accuracy policy (§0.1, P2) as each module ports | mapped per row, as kimiya's ADR 0002 #12 |

### 3.4 Not carried over

- **serde.** There are 28 `Serialize`/`Deserialize` derives, but only `soorat-compat` and `ai`
  depend on them, and both are out (§3.2).
- **Rust-only mechanics.** `#[inline]`, `#[cold]` error helpers, closures (→ function pointer +
  context pointer, §7.3), `&'static` lifetimes, `Option<&'static T>`.
- **The 5 Rust examples** (`basic_thermo`, `heat_engine`, `building_thermal`, `cooling_system`,
  `phase_change`) stay in `rust-old/examples/`. Cyrius versions are an M7 nice-to-have for
  `docs/examples/`.

---

## 4. M0: prep, then initialize with `cyrius port`

`cyrius port` moves the Rust project into `rust-old/` and writes the Cyrius skeleton, so don't
hand-build the structure. The scaffold still overwrites project files on 6.7.18 (issue
`cyrius/docs/development/issues/2026-11-08-cyrius-port-clobbers-project-docs-and-omits-bench.md`,
**open**), so steps 2–3 repair what it writes.

0. **Freeze a green oracle (Rust, before the move).** *Maintainer.*
   - Tag `1.4.0` at `9528e1f`; the release workflow requires the tag to equal `VERSION`.
     Publishing to crates.io is optional, but it makes 1.4.0 the version Rust dependents see.
   - Re-run the P(-1) checks on the tagged commit. They are clean on rustc 1.99.1 (header
     table). After this point nothing in `src/` changes.
   - Run `./scripts/bench-history.sh` on the tag. The last row in `bench-history.csv` is from
     `2f35058` (1.3.0), and `docs/benchmarks-rust-v-cyrius.md` needs a 1.4.0 Rust column.
   - Record glibc (**2.44** today); the oracle's transcendental bits depend on it.
   - Check crates.io reverse dependencies (as kimiya did) and record them in the header table.
     The CLAUDE.md consumer list (kiran, joshua) may be stale.
   - Confirm that every state the generator needs is reachable through the public API (§2.2).
     If something isn't, add it as **unreleased port tooling**, as impetus did: a hidden,
     read-only `oracle` feature, `VERSION` stays 1.4.0, no tag.
1. **Scaffold.**
   ```sh
   cd /home/macro/Repos/ushma
   cyrius port --dry-run .
   cyrius port .
   ```
   Expect these to move into `rust-old/`: `Cargo.toml`, `Cargo.lock`, `rust-toolchain.toml`,
   `src/`, `tests/`, `benches/`, `examples/`, `Makefile`, `codecov.yml`, `deny.toml` and `target/`
   (`--dry-run` doesn't list them by name). Expect the scaffold to write `src/main.cyr`,
   `cyrius.cyml`, `tests/ushma.{tcyr,bcyr,fcyr}`, `docs/{adr,architecture,guides,examples}/` and
   `docs/development/state.md`.
   - **Track `rust-old/Cargo.lock`.** It's untracked today (the library-crate rule), but the
     oracle has to rebuild against exactly hisab 1.4.0 and its other pins. kimiya committed its
     lockfile.
2. **Repair what the scaffold overwrites.** *(Claude; maintainer reviews.)*
   - **`CLAUDE.md`.** Merge `git show HEAD:CLAUDE.md` with the scaffold's Cyrius rules. Keep the
     identity, process, task sizing and DO NOT list. Translate the principles as impetus did:
     - `#[non_exhaustive]` → append-only enum values;
     - `#[must_use]` → `#must_use` / `#pure`;
     - no `unwrap`/`panic!` → no `+?`, `-?`, `*?`, `assert_fatal`, `result_unwrap` or `panic` in
       `src/`;
     - "Vec arena over HashMap" → sorted arrays and typed stride arenas, never `hashmap` (its
       seed is random per process);
     - "`write!` over `format!`" and "Cow over clone" → caller-provided buffers and borrowed
       pointers;
     - "tracing on all operations" → sakshi spans behind a level, added only when a consumer
       asks;
     - "Feature-gate optional deps" → `[features]` for optional dependencies only.

     Replace "Do not commit `Cargo.lock`" with the lockfile rule for `rust-old/`. Refresh the
     consumer list from step 0.
   - **Roadmap** restored from git, with a "2.0 — Cyrius port" section. **`state.md`** rewritten
     for 1.4.0 → 2.0.0-dev.
   - **Workflows.** Keep the scaffold's Cyrius CI and add `test --bounds`, `deps --locked` and
     `audit`. Move the Rust workflows to `rust-old/.github/workflows/`. Append `/rust-old/target/`
     and `/build/` to `.gitignore`. The scaffold says it appends `/rust-old/` but writes
     `/rust-old/target/`; `rust-old/` must stay tracked.
   - **Move Rust-only material** to `rust-old/`: `scripts/bench-history.sh`,
     `scripts/version-bump.sh`, `benchmarks.md`, `bench-history.csv`,
     `docs/architecture/overview.md`. Restore `rust-old/README.md` (the crate's README, as kimiya)
     and add `rust-old/LINES_OF_SOURCE.txt`.
   - **ADRs.** ushma has none yet. Use the scaffold's four-digit form:
     - `0001-cyrius-port-with-rust-oracle.md`: this approach, hisab-first, P1, P5, P6;
     - `0002-port-api-divergences.md`: §3.3;
     - `0003-transcribe-dense-gaussian-elimination.md`: P4;
     - `0004-correct-rounding-is-the-accuracy-target.md`: P2 and P3, adopting kimiya's ADR 0006
       from the start rather than by amendment.
   - Check that `cyrius.cyml` reads `version = "${file:VERSION}"`, not `0.1.0`, and that `bench`
     is in `[test] stdlib`. Document the scaffold's stub functions so `cyrius audit` passes.
     `VERSION` stays 1.4.0 until M7.
   - Copy kimiya's `handoff/gates/gates.sh` to `scripts/gates.sh` and point it at ushma.
3. **`cyrius.cyml`.** Confirm the keys with `cyrius help manifest`. There is no `[lib]` until M1:
   `modules = []` is invalid.
   ```toml
   [package]
   name = "ushma"
   version = "${file:VERSION}"
   license = "GPL-3.0-only"
   language = "cyrius"
   cyrius = "6.7.18"

   [deps]
   stdlib = ["syscalls", "string", "alloc", "str", "fmt", "io", "vec", "math", "ganita",
             "tagged", "result", "overflow", "fnptr"]   # sakshi arrives via hisab: don't list it

   [deps.hisab]
   git = "https://github.com/MacCracken/hisab.git"
   tag = "3.4.1"
   modules = ["dist/hisab.cyr"]

   [test]
   stdlib = ["assert", "bench"]

   [lib]   # from M1, in dependency order
   modules = ["src/error.cyr", "src/constants.cyr", "src/validate.cyr", "src/mathx.cyr", ...]
   ```
   - `[lib]` order: error → constants → validate → mathx → material_data → material → transfer →
     state_data → state → entropy → phase_data → phase → steam_data → steam → cycle → chem_data →
     chem → stat → numerics → numerical.
   - Module dependencies, from Rust's features: entropy and phase need state; steam needs phase;
     cycle needs state, entropy and steam; numerical needs transfer; chem needs state and
     entropy; stat needs state and transfer.
   - hisab's leaves (`dist/hisab.deps`: math, alloc, vec, ganita, result, str, fmt) must be in
     `[deps] stdlib`.
4. **Toolchain and math check** (§7.1). Confirm a Cyrius counterpart for every Rust f64 method
   ushma calls, and re-probe what impetus found on 6.7.18:
   - `f64_min`/`f64_max` return the second operand when the first is NaN (Rust ignores a NaN
     either way). The in-scope modules have 5 `max`, 5 `min` and 7 `clamp` calls, and 1.4.0's
     "f64::max drops NaN" fix depends on these semantics.
   - The stdlib has no public `is_nan`, `is_finite`, `copysign`, `EPSILON`, `INFINITY` or `NAN`
     (18 `is_finite` and 3 `is_nan` calls), so they become local `_us_*` helpers.
   - `powi` (5 calls) becomes repeated multiplication in Rust's order.
5. **Golden-value generator** `scripts/oracle/`: copy kimiya's crate and adapt it.
   - `ushma-oracle` depends on `rust-old/` by path, with the 10 domain features. Its
     `Cargo.lock` is seeded from `rust-old/Cargo.lock`.
   - Use TSV format 1 (`fn<TAB>class<TAB>inputs<TAB>outputs<TAB>note`, f64 values as bit
     patterns), deterministic output, a `--check` mode and glibc in the header (`--check` refuses
     another glibc).
   - Start with `constants` and the **math sweep**: every std method ushma calls, over special
     values and ushma's ranges, including the actual `powf` exponents: the computed ones (γ − 1
     and its expressions, a polytropic k, an interpolation fraction) and the literals 1.5, 2.5,
     0.8, 0.4, 1.6, 4.2, 9/16, −8/27, 1/3 and 1/6.
6. **Table and constant generator** (same crate). Tables come through the oracle's public API or,
   for private data, from the Rust literals with the completeness proof (§2.1):
   - `constants.cyr`: `GAS_CONSTANT`, `ATM`, `STANDARD_TEMP`, `STANDARD_T`, `STEFAN_BOLTZMANN`,
     `BOLTZMANN_K`, `PLANCK`, `AVOGADRO` and `MAX_ADAPTIVE_SUBSTEPS`, plus the private ones
     (`T_CRIT`, `P_CRIT`, `R4_T_MIN`, `QUALITY_ROUNDOFF`, …);
   - `material_data.cyr`: 10 materials;
   - `state_data.cyr`: the gas critical data and acentric factors;
   - `phase_data.cyr`: 8 substances;
   - `chem_data.cyr`: 12 species;
   - `steam_data.cyr`: 42 saturated rows, 72 superheated rows, 9 tier pressures, `R4_N`.

   Write baked typed column arrays with decimal literals in positional form (the lexer has no
   exponent notation), and never a hex pattern in a `const`. `probe_tables.py` compiles every
   cell and compares it bit for bit.
7. **Parity harness.** Copy kimiya's `tests/parity.tcyr`, `tests/parity_core.cyr` and
   `tests/ulp.cyr`:
   - per-module dispatchers;
   - outcomes ok / FAILED / measured / pending / accepted;
   - exit code = the number of failing modules;
   - the ADR 0002 mapping table.

   For M6, add impetus's **state loader** and **step comparator**.
8. **Baseline measurements** → `docs/development/port-parity.md`. Run the math sweep through every
   Cyrius counterpart on x86_64 ELF, x86_64 PE (wine, one process at a time) and aarch64 (qemu).
   The corrected-oracle references come from `mpfr_ref.py`. Confirm the vendored `_us_cr_*`
   primitives give identical bits on all three targets.
9. **Probes.** Each one is a small program with numbers, run before M1 writes code. They go in
   `scripts/probes/` and `port-parity.md`:
   1. **Callback into hisab.** `debye_cv`'s integrand passed to `calc_integral_simpson` as a
      function pointer + context. Compare its bits against Rust hisab 1.4.0 at 200 steps over
      Θ/T ∈ (0, 50]: equal means class B, different means class D (§7.2).
   2. **Memory.** The stdlib `alloc` is a bump allocator with no free. Step a 1D grid 10⁵ times
      (explicit, Crank–Nicolson, adaptive) and run `solve_steady_state` on a 2D grid. The heap
      high-water mark must go flat after warm-up, which needs reusable scratch (§3.3 #2).
   3. **Wide results.** Return a cycle result (state points, heats, efficiencies): out-pointer vs
      handle vs `Result` payload, at the 6.7.18 pin (the open small-struct issue).
   4. **Point lists.** `Vec<DiagramPoint>` as a vec of struct values vs parallel f64 vecs, timed
      for `cycle_ts_diagram` at its point counts.
   5. **Table lookup.** Steam saturated lookups by T and by P, and superheated by (T, P): binary
      search over baked columns, compared with Rust's timings.
   6. **Name collisions** (P6). Include kimiya's `dist/kimiya.cyr` beside a ushma stub and record
      the `duplicate fn` warnings and which definition wins.
10. **File upstream** (§8). Do this before M1, so the requests carry M0's numbers.

**Gate:** `cyrius build` and `cyrius bench` pass; `cyrius test` runs the harness (0 modules) with
no warnings; `cyrius lint`, `fmt --check` and `audit` are clean. Committed: the generator,
`tests/ulp.cyr`, `src/mathx.cyr` with its verification, ADRs 0001–0004, the probe results and
the upstream issues.

---

## 5. Parity policy

tanmatra's classes apply, with kimiya's ADR 0006 reading of class C:

| Class | Tolerance | ushma examples |
|---|---|---|
| **A, exact** | bit-identical; equal integers and codes | constants, every table cell, materials/substances/gases/species, `phase_at` classifications, enum properties, `US_ERR_*` codes, `sublimes_at_atmospheric` |
| **B, arithmetic** | ≤ 2 ulp, expected 0 | + − × ÷ √ only: Fourier, Newton cooling, resistances in series/parallel, Q = mcΔT, diffusivity, Biot, Stefan–Boltzmann (T⁴ by `powi`), ideal gas, Van der Waals, Kay's rule, Carnot, Helmholtz/Gibbs, Hess sums, IF97 Region 4 (algebraic + √), `calc_lerp` interpolation, explicit FD steps, Gauss–Seidel sweeps, transcribed elimination |
| **C, transcendental** | **0 ulp against the corrected oracle** once a function uses the `_us_cr_*` primitives. Rows where glibc misrounds are counted "accepted under ADR 0004" | anything through exp/ln/pow/atan/tanh/expm1: LMTD, ε-NTU, fins, Nusselt correlations, view factors, adiabatic relations, isothermal work, entropy changes, Clausius–Clapeyron, vapor/sublimation curves, ln-P superheated interpolation, K and ln K, van 't Hoff, Boltzmann, partition functions, Maxwell–Boltzmann, Einstein |
| **D, iterative** | relative ≤ 1e-12 | only where the iteration path can differ: the hisab 3.4.1 Simpson (if probe 1 says so), and any solver path an ADR names |

ushma specifics:
- **Never widen a class to make a row pass.** A failure is one of three things:
  - a port defect;
  - a Cyrius/hisab defect (file it upstream);
  - a Rust oracle defect (fix it in 1.4.x and regenerate).

  The only exception is an ADR-accepted divergence, which the harness still reports by name.
- **The reference-value tests travel with the code.** Parity proves "same as Rust"; the 1.4.0
  regression tests prove "right". They cover:
  - IAPWS-IF97 saturation to 1e-8;
  - Incropera Table A.4 air;
  - the Rankine 13.2 % case;
  - the unit-square perpendicular view factor of 0.2000;
  - diamond's Debye Cv at 3 K;
  - propane's ln K ≈ 837.
- **Test branch points on both sides.** Every switch gets golden rows just below and just above
  it:
  - steam tier boundaries and exact tier rows;
  - T_sat and the first tabulated row of each tier;
  - the quality round-off clamp (1e-9);
  - `phase_at`'s piecewise curve joins (triple, boiling and critical points) and the
    sublimation branch;
  - the view-factor asymptote at X, Y < 0.01;
  - Debye's T³ switch at Θ/T > 50, and Θ/T > 709;
  - the explicit stability limit Fo·(1 + Bi) = 0.5;
  - the adaptive sub-step cap;
  - Pitzer's Z ≤ 0 error;
  - K leaving the f64 range;
  - the refrigeration Carnot limits (T_L ≤ 0).

---

## 6. Milestones

Every milestone ends with the siblings' gates:
- `cyrius build`;
- `cyrius test` and `test --bounds` **with no warnings** (the f64 type diagnostics appear only
  there);
- the full corpus on ELF, PE and aarch64 via `scripts/gates.sh`;
- `cyrius bench` recorded against the Rust numbers;
- `vet`, `lint`, `fmt --check` and `audit` on the module's files;
- `distlib --all` and `distlib --check`;
- `deps --locked`;
- **surface parity**: every in-scope Rust `pub fn` has a Cyrius counterpart or a written reason;
- the milestone's **range-edge sweep** (P2);
- CHANGELOG and `state.md`.

| Milestone | Modules | Key work | Golden grid (minimum) |
|---|---|---|---|
| **M1 Foundation & transfer** | error, constants, validate, mathx, material, transfer | `US_ERR_*` codes + detail; `finite_gt`/`finite_ge` with Rust's NaN rejection; vendored correctly rounded primitives + new `tanh`; 10-material table; 32 transfer functions: conduction/convection/radiation, resistance networks, fins (tanh), LMTD (equal-ΔT limit), ε-NTU (C_r = 1 branch), Re/Pr/Nu correlations (`powf`), Blasius, three view factors (atan, ln, the asymptote switch) | every material × field; each function over seeded physical grids plus NaN/±∞/0/negative edges; LMTD with ΔT₁ = ΔT₂ and 1 ulp apart; NTU 0…100 × C_r 0…1; view factors across the 0.01 switch and at the unit-square reciprocity case |
| **M2 Equations of state & entropy** | state, entropy | Ideal gas, Van der Waals, Redlich–Kwong, Peng–Robinson pressures; virial; Pitzer–Abbott Z (Err at Z ≤ 0); Kay's rule; isothermal/isobaric/adiabatic work and temperature; ideal-gas, isothermal, heat-transfer and mixing entropy; Carnot efficiency and COP; Helmholtz/Gibbs; Clausius inequality | every gas × T × V_m grid incl. near-critical; Pitzer Tr × Pr incl. the Z ≤ 0 edge; mixtures with zero, NaN and negative fractions; Carnot at T_c ≥ T_h |
| **M3 Phase & steam** | phase, steam | Substance table; `phase_at` with the piecewise ln P vs 1/T curves and sublimation; latent heats incl. the CO₂ sublimation case; Clausius–Clapeyron (overflow → Err). Steam: generated tables; exact IF97 Region 4 both ways; saturated by T and by P (P·v_g/T interpolation); quality from v/h/s with the round-off clamp; superheated ln-P tier interpolation with the saturated-vapor lower row; `superheated_by_entropy` | `phase_at` on a T × P grid per substance with rows on both sides of every join; every steam row through the lookups; Region 4 over 273.16…647.096 K and its inverse; between-row and between-tier interpolation; every tier through `superheated_by_entropy`; x = 1 + 2e-16 |
| **M4 Cycles** | cycle | Otto, Diesel, Brayton, Rankine (superheated exhaust via the s-inversion), refrigeration (compressor exit via the s-inversion, negative net work), heat-pump COP, `compare_cycles` (Carnot COP for refrigeration), T-s and P-v diagrams along real process paths; results through out-pointers (probe 3) | r, γ, T grids for the gas cycles; Rankine/refrigeration over the table range incl. the 13.2 % case and COP ≤ Carnot; every cycle kind's diagram point lists; `processes` shorter than `state_points` |
| **M5 Chemical & statistical** | chem, stat | 12-species table; Hess's law; reaction Gibbs; K and `ln_equilibrium_constant` (range-exit Err); van 't Hoff; flame temperature (ushma's loop). Boltzmann probability/ratio, S = k ln W, equipartition, three partition functions (vibrational via expm1), Maxwell–Boltzmann pdf and speeds, Einstein Cv, Debye Cv via hisab Simpson (function-pointer callback) with the T³ switch | every species × function; propane ln K ≈ 837; partitions at hν ≪ kT; Debye over Θ/T ∈ (0, 50] ∪ (50, 10⁴] incl. diamond at 3 K |
| **M6 Numerical** | numerics, numerical | `ThermalGrid1D`: explicit (Fo·(1 + Bi) limit), Crank–Nicolson with caller-reusable scratch, adaptive sub-stepping, `max_stable_dt`. `ThermalGrid2D`: flat storage (§3.3 #2), Gauss–Seidel with convective BCs, convergence on non-finite fields. `ThermalNetwork`: connectivity check, row scaling, transcribed elimination (P4). Malformed public state → Err | single-step rows: every BC kind × explicit/CN/adaptive incl. Fo·(1 + Bi) exactly 0.5; 2D steady states with their iteration counts; network solves incl. a floating node and R ≈ 1e14; the 10⁵-step soak from probe 2 |
| **M7 Release 2.0.0** | all | Rust-test ledger (tanmatra's script): each of the 397 + 30 tests ported, folded into golden rows, or marked n/a / deferred; full surface-parity diff; `docs/benchmarks-rust-v-cyrius.md` (63 criterion benchmarks); `cyrius distlib --all`; `cyrius audit`; migration guide (`docs/guides/migration-2.0.md`, covering §3.3 and P6's names); VERSION 2.0.0 | full corpus green on x86_64 ELF, x86_64 PE and aarch64 |

Order: M1 → M2 → M3 → M4 → M5 → M6 → M7. M6 is last because it carries the hisab gap, the
memory design and the stateful parity layer. M5 could move earlier if probe 1 is clean. Under the
CLAUDE.md task-sizing rule, M3 (steam) and M4 (cycles) are large: split each into verified
sub-parts (M3a phase, M3b steam tables and Region 4, M3c superheated; M4a gas cycles, M4b steam
cycles, M4c diagrams).

---

## 7. Rust → Cyrius mapping (ushma specifics)

### 7.1 Math functions

Every std f64 method in the non-test code of the in-scope modules (§3.1), 1.4.0. The deferred
`bridge` adds its own exponents (Kleiber's 0.75; DuBois's 0.425 and 0.725), which port with it.

| Rust (calls) | Cyrius | Parity |
|---|---|---|
| `sqrt` (20), `abs` (24), `ceil` (1) | `f64_sqrt`, `f64_abs`, `f64_ceil` (builtins) | exact |
| `exp` (14), `ln` (22), `exp_m1` (2) | `_us_cr_exp`, `_us_cr_ln`, `_us_cr_expm1` (P3) | class C, same bits on every target |
| `powf` (23) | `_us_cr_pow`; exact-constant exponents reviewed under the accuracy policy (1/3 → a correctly rounded cbrt; 1/6, 9/16, −8/27) | class C |
| `atan` (7) | `_us_cr_atan2(x, 1.0)` (atan2 with x = 1 > 0 is atan exactly) | class C |
| `tanh` (2) | **new** `_us_cr_tanh` (fins) | class C |
| `powi` (5) | repeated multiplication in Rust's order | class B |
| `max` (5), `min` (5), `clamp` (7) | local helpers with Rust's NaN semantics (step 4) | exact |
| `is_finite` (18), `is_nan` (3) | `x - x == 0`, `x != x` helpers | exact |

### 7.2 hisab: what the oracle used vs what hisab 3.4.1 provides

| Rust hisab 1.4.0 (used by) | hisab 3.4.1 Cyrius | Plan |
|---|---|---|
| `calc::lerp` (steam, 3 sites) | `calc_lerp`: same `a + (b − a)·t` fast form; a half-scale repair only where the fast form overflows | use hisab; class B |
| `calc::integral_simpson` (stat `debye_cv`, n = 200) | `calc_integral_simpson(f, a, b, steps, out)`, function pointer + out-pointer; reframed in 3.3.6 for subnormal widths | use hisab; class B or D per probe 1 |
| `num::gaussian_elimination` (`ThermalNetwork::solve`) | **none**: `num_tridiag_solve` and iterative solvers only | **gap** → transcribe into `numerics.cyr` (P4) and file it |
| *(the roadmap's "bisection")* | `num_bisection` | not called in 1.4.0; ushma's s-inversions are its own loops (§2.3) |

### 7.3 Types and idioms

| Rust | Cyrius |
|---|---|
| `Result<T, UshmaError>` | value-form `Result` (stdlib `result`) + `US_ERR_*` code enum and numeric detail |
| `ThermalMaterial`, `SubstanceData`, `GasData`, `Species` `pub const`s and `ALL_*` slices | baked column tables + index; named index constants for the presets |
| `Phase`, `ProcessKind`, `CycleKind`, `Side` | append-only integer enums |
| `BoundaryCondition` (with `Convective { h, t_inf }` data) | tag + payload fields (stdlib `tagged`) |
| cycle result structs, `Vec<DiagramPoint>` | caller out-pointers / handles (probe 3); a vec or parallel f64 vecs (probe 4) |
| `Vec<Vec<f64>>` 2D field | flat row-major buffer + `nx`, `ny` (§3.3 #2) |
| local closures passed to hisab | function pointer + context pointer (stdlib `fnptr`) |
| `usize` counts, `num_nodes − 1` | i64 with explicit checks; stdlib `overflow` operators where Rust used checked arithmetic |
| `#[must_use]` / `#[non_exhaustive]` | `#must_use` / append-only enum values |
| `#[tracing::instrument]` | none (§3.3 #3) |

The siblings have already probed these language traps on 6.7.18:
- there are no exponent literals;
- `(0 - x)` on an f64 is wrong;
- a hex literal in a `const` is an integer;
- `f64_round` ties to even;
- `+?` aborts the process;
- `cyrius test` must be warning-free;
- the compiler refuses a source over 24 MB, includes counted (split large probes);
- `cyrius audit` rejects comment lines over 120 characters;
- `cyrius bench` stops any file after 300 s.

---

## 8. Upstream (file at M0, then as found)

| To | Request | Evidence |
|---|---|---|
| ganita (cyrius stdlib) or hisab | Host the shared correctly rounded primitives (exp, expm1, exp10, ln, log10, pow, sin, cos, acos, atan2, tanh) so ports stop vendoring them | three private copies: impetus `crmath.cyr`, kimiya `mathx.cyr`, ushma `mathx.cyr` |
| hisab | Dense direct solver (Gaussian elimination with partial pivoting, or LU) | ushma `ThermalNetwork::solve`; kimiya `simulate_stiff` (its ADR 0003) |
| cyrius | `cyrius port` still clobbers project docs (existing issue) | add ushma as another reproduction if 6.7.18 still does it |
| cyrius | A duplicate non-private `fn` across bundles is a warning and the last one wins | probe 6; ushma × kimiya, 10 names |
| cyrius stdlib | `is_nan`, `is_finite`, `copysign`, `EPSILON`, `INFINITY`, `NAN`; Rust-semantics `min`/`max` | check impetus's filing first and add to it rather than duplicate it |
| kimiya | Its `KM_GAS_CONSTANT` is 8.314462618, truncated at 1e-11 relative; ushma 1.4.0 uses the exact N_A·k_B = 8.31446261815324 | for its bridge phase; values crossing the ushma ↔ kimiya bridges disagree in the 10th digit |

---

## 9. After 2.0.0

- Bridges (with kimiya's ushma bridges) once a Cyrius consumer exists; `soorat-compat` once soorat
  ports.
- Switch `mathx.cyr` and `numerics.cyr` to the upstream versions when §8 lands, then re-measure.
- The roadmap backlog, demand-gated.
- Rust 1.4.x stays maintained for soorat and crates.io dependents: oracle fixes only.

---

## 10. Tracking

| Milestone | Status | Golden rows | Reference tests | Benches | Notes |
|---|---|---|---|---|---|
| M0 Prep & initialize | ☐ planned 2026-11-11 | — | — | — | Blocked on step 0: tag 1.4.0. §0.2 to confirm |
| M1 Foundation & transfer | ☐ | | | | |
| M2 Equations of state & entropy | ☐ | | | | |
| M3 Phase & steam | ☐ | | | | Large: M3a/b/c |
| M4 Cycles | ☐ | | | | Large: M4a/b/c |
| M5 Chemical & statistical | ☐ | | | | Probe 1 decides Simpson's class |
| M6 Numerical | ☐ | | | | ADR 0003 (P4); probe 2 (memory) |
| M7 Release 2.0.0 | ☐ | | | | |
