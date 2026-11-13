# ushma Rust → Cyrius port plan

> **Plan of record** for porting ushma to Cyrius. Rust **1.4.0** is the final Rust feature
> release and the **reference implementation** (the oracle); the Cyrius port ships as **2.0.0**.
> **Everything ports**: every module, every feature, every test, the examples and the benchmarks.
>
> The parity method (frozen oracle, golden rows, classes A–D) is tanmatra's
> ([its ADR 0006](../../../tanmatra/docs/adr/0006-cyrius-port-with-rust-oracle.md)), as kimiya
> and impetus refined it. The **scope and the Rust → Cyrius mapping come from the Cyrius 6.7.18
> language guide and stdlib** (§7), not from the sibling ports.
>
> ▶ **Start on Cyrius 6.7.18 with hisab 3.4.1**, after 1.4.0 is tagged (M0 step 0).
>
> 🟡 **Status:** planned 2026-11-11, revised the same day (everything in scope; stdlib-first
> mapping). The maintainer **accepted P2–P8** (§0.2) on 2026-11-11; ADRs 0001–0005 record them at
> M0 step 2. M0 starts once 1.4.0 is tagged as it stands (step 0). **The Rust code is never edited**,
> and every fix goes into the Cyrius port.

| Thing | Value | Source |
|---|---|---|
| Rust source (oracle) | ushma **1.4.0** at `9528e1f` (committed, **not yet tagged**): ≈5,400 code lines in `src/` (no tests, comments or blank lines), 155 public functions and methods, 18 files; 10 domain features plus `soorat-compat`, `ai`, `logging` (`full` = all) | `Cargo.toml`, `VERSION`, `git tag` |
| Rust health, 2026-11-11 on rustc 1.99.1 | ✅ 449 unit + 30 integration tests (no doc tests); `cargo fmt --check`, clippy `--all-features --all-targets -D warnings`, rustdoc `-D warnings`, `cargo audit` (224 crates, clean), `cargo deny` (advisories, bans, licenses, sources ok) | local run |
| Correctness baseline | The 1.4.0 audit: IAPWS-IF97 steam tables, exact Region 4, real process paths in the diagrams, NaN/±∞ rejected by every guard, 80 regression tests with IF97, Incropera and NIST reference values | `CHANGELOG.md` 1.4.0 |
| Cyrius toolchain | **6.7.18** (installed, latest) | `~/.cyrius/current` |
| Cyrius libraries the port uses | **math** (f64 ops, `F64_*`, `f64_lerp`/`_min`/`_max`/`_clamp`, the same-bits exp/ln/atan polyfills, `f64_parse`); **ganita 1.2.15** (`ganita_f64_pow`/`_tanh`/`_atan2`/`_cbrt`, `ganita_mat_gaussian_elim`, `ganita_mat_lu`/`_lu_solve`); **sakshi 2.5.10** + **log** (spans, levels, outputs); **sandhi 1.13.7** (HTTP/1.1 + HTTP/2 client, headers, TLS policy) with its companions tls, async, random, fdlopen, dynlib; **bayan 1.11.1** (JSON parse/build); **io** (`getenv`); result, tagged, boxed, vec, fnptr, alloc arenas, freelist, overflow | `~/.cyrius/lib`, `cyrius/docs/stdlib-modules.md`, `stdlib-reference.md` |
| hisab | **3.4.1**, Cyrius, tag `3.4.1`, bundle `dist/hisab.cyr` (`calc_lerp`, `calc_integral_simpson`) | `/home/macro/Repos/hisab` |
| Oracle's math | Rust **std** → host **glibc 2.44**: `powf` (29 calls), `ln` (23), `exp` (14), `atan` (11), `tanh` (2), `exp_m1` (2). Rust **hisab 1.4.0** (pinned in the lockfile): `calc::lerp` (3), `calc::integral_simpson` (1), `num::gaussian_elimination` (1) | `grep` of non-test `src/`, `Cargo.lock` |
| Services `ai` talks to | **daimon 2.5.6** and **hoosh 2.12.3**, both Cyrius. daimon registers agents at `POST /v1/agents` and answers `{"id":N,…}` (§2.6) | `/home/macro/Repos/daimon/src/router.cyr`, `api_agent.cyr` |
| Dependents | Local: **soorat 1.0.0** (Rust), optional feature `thermo` → `ushma = "1"` with `soorat-compat`. kiran, joshua, badal and dravya have no local repo. crates.io reverse dependencies are checked at M0 step 0. **No Cyrius consumer yet.** | `~/Repos/*/Cargo.toml`, `cyrius.cyml` |
| Target | ushma **2.0.0**, Cyrius | this document |

## 0. Decisions

### 0.1 Decided

- **Everything ports** (maintainer, 2026-11-11). All 18 Rust files and every feature:
  - the 10 domain modules;
  - `bridge` and `integration/soorat`;
  - `logging`, on `lib/log.cyr` and sakshi;
  - `ai`, on sandhi;
  - serde, as JSON codecs on bayan for all 28 types that derive it;
  - all 449 unit and 30 integration tests, the 63 benchmarks and the 5 examples.

  Nothing is dropped or deferred. If something has no Cyrius counterpart, the port shows the
  evidence and asks. It never leaves the thing out.
- **The Rust code is frozen and never touched** (maintainer, 2026-11-11). Rust 1.4.0 moves to
  `rust-old/` as it stands, with no fixes, no 1.4.x releases and no port tooling added. A
  generator crate in `scripts/oracle/` drives it by path, read-only. When the oracle is wrong,
  the port does it right, and an ADR 0002 row maps the oracle's rows to the corrected value.
- **No further Rust releases.** API shape changes and fixes land in the Cyrius API only, each
  accepted by an ADR row (§3.3).
- **Stdlib- and hisab-first.** The port calls the Cyrius stdlib (math, ganita, sakshi, log,
  sandhi, bayan, …) and hisab for everything they provide. Where they fall short, the port files
  the gap upstream (§8) with the probe that shows it. A stopgap is local and named as one.
- **Accuracy over glibc parity** (the maintainer's decision in kimiya, its ADRs 0005 and 0006).
  Every transcendental primitive is correctly rounded and gives the same bits on every target.
  "Parity" means functional parity: the same results, `Err` codes, check order and branches.
  It does not mean copying glibc's last-ulp errors. The class-C reference is the *corrected
  oracle* (§5).
- **Performance is measured, not gated.** Every Rust benchmark gets a Cyrius twin, and both
  columns go in `docs/benchmarks-rust-v-cyrius.md`.
- **New physics stays out.** The roadmap backlog comes after 2.0, demand-gated. It is new work,
  not existing features.
- **The maintainer does all git.** Each verified milestone part ends with a summary and a request
  to commit before the next part starts.

### 0.2 Accepted 2026-11-11 (maintainer: "Accept all"; ADRs 0001–0005 at M0 step 2)

| # | Decision | Why | Rejected alternative |
|---|---|---|---|
| P2 | **Range-edge sweep inside each milestone**, not a library-wide audit afterwards. | kimiya ran its audit after M4 and had to reopen every module (batches B1–B12). | Port first, audit later. |
| P3 | **Correctly rounded primitives** in `src/mathx.cyr` (`_us_cr_exp`, `_ln`, `_expm1`, `_pow`, atan as atan2(x, 1), and tanh), vendored from kimiya's `src/mathx.cyr`, plus an upstream request to host one shared set. | The stdlib's exp/ln/atan polyfills and ganita's pow/tanh are ≤1 ulp with the same bits everywhere, but they are **not correctly rounded**, and the accuracy policy (§0.1) requires correct rounding. The stdlib has no expm1. ushma would be the third private copy, after impetus `crmath.cyr` and kimiya `mathx.cyr`. | Use the stdlib polyfills and ganita, and accept ≤1 ulp under an ADR. |
| P4 | **`ThermalNetwork::solve` on `ganita_mat_lu` + `ganita_mat_lu_solve`**, after 1.4.0's diagonal row scaling (amended by probe 9; first decided as `ganita_mat_gaussian_elim`). | The stdlib already has the solver; nothing to transcribe. Probe 9 (2,178 networks): LU is bit-identical on all 2,000 systems 1.4.0 solves and refuses exactly the 178 it refuses, because after row scaling LU's relative threshold is 1.4.0's absolute `ε`. Elimination's threshold spans [A\|b] and refuses 259 systems 1.4.0 solves. ADR 0003 and its amendment. `_lu_solve` allocates 8·n B per call (§8). | `ganita_mat_gaussian_elim` (rejected by probe 9). |
| P5 | **Two bundles, as prakash:** `[lib]` (every module except the HTTP client, logging included) and `[lib.ai]` (`[lib]` + `src/ai.cyr`). | sandhi pulls the whole TLS stack (tls, async, random, fdlopen, dynlib). A consumer that wants only the physics shouldn't link it, which is what Rust's `ai` feature did. The 10 domain features are additive, so they need no profiles. | A single bundle with ai included. |
| P6 | **Public names live under `us`.** Functions become `us_ideal_gas_pressure`, `us_ai_register_agent`, …. Enum variants, constants and error codes get `US_` names (`US_PHASE_SOLID`, `US_ERR_INVALID_TEMPERATURE`). Probe 6 picks the mechanism: the language's own `mod us;` mangling (a consumer then writes `use us.ideal_gas_pressure;` and calls the bare name), or a written prefix. | Cyrius has **one global namespace**. A duplicate non-private `fn` is a warning, and **the last definition wins** (guide L5267–5279). A bare enum variant name also resolves to the last definition, retroactively (L391–396, L10297). Collisions already exist: 10 functions with kimiya's bundle (`ideal_gas_pressure`, `reaction_enthalpy`, `equilibrium_constant`, …), and prakash's `ai_*` client names. kiran uses both kimiya and ushma. `mod`/`use` (6.7.17, L625) is documented in one sentence and no library uses it yet, hence the probe. | Keep Rust's names and rename only the collisions. |
| P7 | **JSON is byte-identical to serde_json**, with serde's field semantics. Readers work on bayan's value tree (`bayan_json_v_parse_buf`), writers on `str_builder`. The f64 text comes from a shortest-digit printer (`src/fmtf.cyr`, vendored from impetus), which also prints `UshmaError`'s Display text the way Rust does. | `#derive(Serialize)` exists (guide L4548–4660, on bayan), but it isn't serde's semantics: a missing key leaves 0 instead of an error; u64, `Vec<Str>`, tuple and payload-enum fields are refused; its f64 text is bayan's printer, which round-trips but isn't always shortest and formats differently from serde_json. Its scanner also has three open issues (2026-11-07: reads past a NUL, memsets a refused allocation, hangs on a stray brace). impetus reached the same decision (its ADR 0010). | `#derive(Serialize)`, compared by parsed value only, with its decode differences accepted by ADR. |
| P8 | **The port talks to daimon's real API: `POST /v1/agents`, reading `id`.** Rust 1.4.0 posts to `/v1/agents/register` and reads `agent_id`, which matches no daimon that exists. daimon is Cyrius-only (2.5.6), and there is no Rust daimon. prakash made the same correction in its 2.5.1 (`prakash/src/ai.cyr`). | The Rust code is frozen (§0.1), so the fix lives only in the Cyrius port, as ADR 0002 row #6. | Mirror Rust's path, which would never register. |

---

## 1. References

- **The Cyrius language and stdlib** (`/home/macro/Repos/cyrius`), reviewed in full for this
  plan on 2026-11-11:
  - `docs/guides/cyrius-guide.md` (all 10,725 lines; its newest notes are for 6.7.17);
  - `docs/stdlib-modules.md` and `docs/stdlib-reference.md`;
  - `docs/development/ecosystem-migration.md` (6.6–6.7.18);
  - `docs/ecosystem.md` (fold pins);
  - `cyrius help` and `cyrius help manifest`.

  §7 maps every Rust construct and dependency ushma uses onto them, with guide line numbers.
- **sandhi 1.13.7** (`/home/macro/Repos/sandhi`): `README.md` (the client and its companion
  modules; "use bayan, not json"), `docs/examples/`.
- **sakshi 2.5.10** (`lib/sakshi.cyr`) and `lib/log.cyr`: spans (`sakshi_span_enter`/`_exit`),
  levels, outputs, the `SAKSHI_LEVEL` compile-time ceiling.
- **bayan 1.11.1** (`lib/bayan.cyr`): `bayan_json_parse`, `bayan_json_v_*`, `bayan_json_build`.
- **hisab 3.4.1** (`/home/macro/Repos/hisab`): its `CLAUDE.md`, `docs/guides/migration-3.4.md`,
  and the doc comments of `calc_lerp` and `calc_integral_simpson`.
- **Process models** (method only, never scope):
  - **kimiya**, in progress: the generator crate `scripts/oracle/` (`validate.py`, `mpfr_ref.py`,
    `cr_verify.py`, `probe_tables.py`, `mathx_tables.py`), the parity harness
    (`tests/parity.tcyr`, `parity_core.cyr`, `ulp.cyr`), `src/mathx.cyr`, ADR 0002's divergence
    format, and `handoff/` (gates script and the per-batch verification loop);
  - **impetus 2.0.0**: single-step parity rows, M0 probes (`scripts/probes/`), private-constant
    extraction (`t_private.rs`), serde-compatible JSON (ADR 0010, `src/json.cyr`,
    `src/fmtf.cyr`);
  - **prakash 2.15.2**: the Cyrius AI client on sandhi (`src/ai.cyr`, `tests/ai.tcyr`,
    `tests/ai_bundle.tcyr`), the `[lib.ai]` profile and the sakshi resolution note in
    `cyrius.cyml`, and `src/serialize.cyr`;
  - **tanmatra 2.0.2**: `scripts/rust_test_ledger.py`, `scripts/surface_parity.py`, the class-C
    port-parity rule.

---

## 2. What this port has to handle

1. **Private steam tables.** `SATURATED_TABLE` has 42 rows × 9 fields; `SUPERHEATED_DATA` has 72
   rows (9 pressure tiers × 8 temperatures); `R4_N` holds the 10 IF97 Region 4 coefficients. The
   API exposes interpolated values, not the rows. They are extracted from `rust-old/src/steam.rs`
   with Rust's `f64::from_str`, plus a proof that every `const` in the module is emitted or
   excused (impetus's `t_private.rs`). They are never re-typed.
2. **Stateful objects with public state.** `ThermalGrid1D`, `ThermalGrid2D` and `ThermalNetwork`
   step or solve in place. Their fields are public, so the generator dumps the full state through
   the public API without touching the Rust code (confirm at M0 step 0). They get single-step
   parity rows: load the state, step once, compare the whole state.
3. **ushma's own iterations**, transcribed exactly and expected class B:
   - Rankine's and refrigeration's s(T, P) = s₁ inversions;
   - `superheated_by_entropy`;
   - Gauss–Seidel;
   - adaptive sub-stepping;
   - the flame temperature.
4. **Wide results and point lists.** Each cycle function returns a struct of state points, heats
   and efficiencies; the diagrams return `Vec<DiagramPoint>`. The open small-struct `Result`
   payload issue
   (`cyrius/docs/development/issues/2026-11-08-small-struct-result-payload-read-as-address.md`)
   means these go through out-pointers or handles (probe 3).
5. **All three optional stacks at once.** No sibling port has carried tracing spans, an HTTP
   client and serde JSON together:
   - **Logging:** `logging::init()` / `init_with_level()` become `us_log_init()` /
     `us_log_init_with_level()` on `lib/log.cyr`, reading `USHMA_LOG` through `getenv`. All 28
     `#[tracing::instrument]` spans are `level = "debug"`. They become `sakshi_span_enter`/`_exit`
     with the same names, at DEBUG. Probe 8 measures the cost when spans are filtered out, against
     functions that take 1–2 ns.
   - **HTTP:** `DaimonClient::register_agent` becomes `us_ai_register_agent` on
     `sandhi_http_post_opts`: serde's body bytes (`capabilities` before `name`), a bearer header and
     `Accept: */*` through `sandhi_headers_*`, `total_ms` and `read_ms` at 30 s each (probe 7:
     `total_ms` alone doesn't bound a silent server), and the status check. `DaimonConfig`'s redacted `Debug` becomes a formatter that
     prints `<redacted>`. Endpoint and key strings are refused when they hold control bytes (the
     request-splitting guard prakash added).
   - **JSON:** 28 types (§3.1), 30 serde round-trip test sites, documents byte-identical to
     serde_json (P7).
6. **An oracle defect in `ai`** (P8). The Rust client posts to the wrong path and reads a field
   daimon never sends. The Cyrius client uses daimon 2.5.6's API, and the frozen Rust stays as it
   is.
7. **sakshi arrives twice.** hisab 3.4.1's `cyrius.cyml` declares `[deps.sakshi]` (tag 2.5.10),
   although its bundle makes no `sakshi_` call. ushma needs sakshi as a stdlib leaf for
   `lib/log.cyr` and sandhi. Whenever both are present, the resolver warns: "refusing to overwrite
   stdlib leaf 'sakshi' … declare it in [deps] stdlib OR depend on the package, not both"
   (`cbt/deps/stdlib.cyr:251–259`). The stdlib copy wins, and the two are byte-identical at 2.5.10,
   so this is benign but can't be silenced from ushma. The fix is upstream: hisab drops
   `[deps.sakshi]`, as patra did at 1.13.0 (§8). Probe 10 confirmed it with real resolves: every one warns (ushma's
   `deps`, `distlib --all`, `distlib --check` and each consumer's `deps`), and all three copies are byte-identical
   (sha256 `8a819554…`).
8. **Open toolchain issues the port must work around** (`cyrius/docs/development/issues/`, open
   at 6.7.18):

   | Issue | Severity | Hits ushma | Workaround in the port |
   |---|---|---|---|
   | `field-store-of-struct-call-copies-a-field-argument` | 🔴 High, **silent wrong code** | `s.f = g(t.e)` stores `t.e`, not `g`'s result: the shape of stepping a grid in place | bind the call to a local, then store; a CLAUDE.md rule and a probe-12 case |
   | `small-struct-result-payload-read-as-address` | 🟡 Medium | a struct ≤ 8 B in `Ok()` crashes when bound back | error/result records go through out-pointers (probe 3) |
   | `no-scoped-allocator-for-alloc` | 🟡 High | hisab allocates from the global heap on every call and never frees; `alloc_reset` clobbers the caller | probe 2 measures heap growth per `debye_cv` (Simpson) and per steam lookup (`calc_lerp`); a hisab call that leaks stays off per-call paths, or the formula is evaluated locally with a stated reason |
   | `arena-grow-keeps-a-chunk-per-larger-allocation` | 🔴 Medium | growable scratch arenas for the grids | size the arena's chunk ≥ the largest per-step need |
   | `exp-polyfill-double-rounds-subnormal-results` | 🟡 Low | aarch64 `f64_exp`/`f64_exp2` misround (up to 39 ulp measured downstream); x87 misrounds about 5 % | the correctly rounded primitives (P3) |
   | `derive-serialize-*` (three) | 🔴/🟡 | the derived JSON decoder | hand-written codecs on bayan's tree (P7) |
   | `f64v-builtins-no-scalar-tail` | 🔴 High | any SIMD sweep over a grid row | no `f64v_*` without an explicit scalar tail |
   | `an-attribute-before-a-var-attaches-to-the-next-fn` | 🟡 Medium | `#must_use`/`#inline` near the table globals | never put an attribute before a `var` or `const` (CLAUDE.md rule) |
   | `manifest-build-profiles-and-target-unread`, `constant-if-arms-not-folded` | 🟡 Low | feature gating | `[lib.<profile>]` plus `#ifdef`, never `if (CONST)` |
   | `cyrius-port-clobbers-project-docs-and-omits-bench` | 🟡 Medium | M0 | step 2 |

   No open issue concerns sandhi's client or sakshi.

---

## 3. Scope

### 3.1 Everything, by milestone

| Milestone | Rust module | Code lines | Public fns | Unit tests | Cyrius module(s) |
|---|---|---|---|---|---|
| M1 | *(correctly rounded primitives, P3)* | — | — | — | `src/mathx.cyr` (private `_us_cr_*`) |
| M1 | *(JSON infrastructure, P7)* | — | — | — | `src/fmtf.cyr` (vendored), `src/json.cyr` (shared reader/writer helpers) |
| M1 | `error.rs` | 37 | — | 18 | `src/error.cyr` (codes, names, messages) |
| M1 | `validate.rs` | 19 | (private) | 2 | `src/validate.cyr` |
| M1 | `logging.rs` | 14 | 2 | 1 | `src/logging.cyr` (log + sakshi, `USHMA_LOG`) |
| M1 | *(9 public constants in 5 modules)* | — | (consts) | — | `src/constants.cyr` (generated) |
| M1 | `material.rs` | 106 | 2 | 8 | `src/material_data.cyr` (generated), `src/material.cyr` |
| M1 | `transfer.rs` | 443 | 32 | 67 | `src/transfer.cyr` |
| M2 | `state.rs` | 473 | 21 | 47 | `src/state_data.cyr` (generated), `src/state.cyr` |
| M2 | `entropy.rs` | 132 | 9 | 23 | `src/entropy.cyr` |
| M3 | `phase.rs` | 353 | 10 | 44 | `src/phase_data.cyr` (generated), `src/phase.cyr` |
| M3 | `steam.rs` | 1,320 | 10 | 46 | `src/steam_data.cyr` (generated, §2.1), `src/steam.cyr` |
| M4 | `cycle.rs` | 730 | 13 | 59 | `src/cycle.cyr` |
| M5 | `chem.rs` | 203 | 7 | 23 | `src/chem_data.cyr` (generated), `src/chem.cyr` |
| M5 | `stat.rs` | 214 | 14 | 33 | `src/stat.cyr` |
| M6 | `numerical.rs` | 599 | 13 | 27 | `src/numerical.cyr` (on ganita, P4) |
| M7 | `bridge.rs` | 168 | 16 | 31 | `src/bridge.cyr` |
| M7 | `integration/soorat.rs` | 215 | 5 (+6 types) | 13 | `src/soorat.cyr` |
| M7 | `ai.rs` | 86 | 1 (+3 types) | 7 | `src/ai.cyr` (`[lib.ai]` only) |

**Totals:** 155 of 155 public functions, 449 of 449 unit tests, 30 of 30 integration tests, 63
benchmarks and 5 examples (M8). `lib.rs` (29 lines of module wiring and re-exports) becomes the
`[lib]` / `[lib.ai]` module lists.

Each module's types get their serde-compatible JSON codecs in that module's milestone. The 28
types (`UshmaError` derives no serde) are:
- `ThermalMaterial`;
- `GasData`;
- `Phase`, `SubstanceData`;
- `SaturatedEntry`, `SuperheatedEntry`, `WetSteamProperties`;
- `StatePoint`, `ProcessKind`, `CycleKind`, `CycleResult`, `DiagramPoint`, `RefrigerationResult`,
  `CycleComparisonEntry`;
- `Species`;
- `BoundaryCondition`, `Side`, `ThermalGrid1D`, `ThermalGrid2D`, `ThermalNetwork`;
- the 6 soorat types;
- `DaimonConfig`, `HooshConfig`.

### 3.2 Rust features → Cyrius

| Rust feature | Cyrius |
|---|---|
| `transfer`, `state`, `entropy`, `phase`, `steam`, `cycle`, `numerical`, `chem`, `stat`, `material` | modules in `[lib]`, in dependency order (§4 step 3). They were additive compile-time gates. Cyrius `[features]` only select optional dependencies, so nothing needs a profile |
| `soorat-compat` | `src/soorat.cyr` in `[lib]`; it has no extra dependencies |
| `logging` | `src/logging.cyr` in `[lib]` (log + sakshi). Spans compile to nothing above the `SAKSHI_LEVEL` ceiling |
| `ai` | `src/ai.cyr` in `[lib.ai]` (sandhi + its TLS companions + bayan), as prakash (P5) |
| `full` | `[lib.ai]` |

### 3.3 API changes in the Cyrius version (accepted divergences, ADR 0002)

This is the starting list. Rows found during the port are added there, one per function family,
and anything that moves results noticeably goes to the maintainer first.

| # | Rust 1.4.0 | Cyrius 2.0.0 | Golden comparison |
|---|---|---|---|
| 1 | `UshmaError` variants (thiserror) with `String` and f64 payloads, and formatted Display text | `US_ERR_*` codes (append-only) in the value-form `Result`'s payload; `us_error_name` / `us_error_message` give each code a fixed message; the failing values go to a level-gated sakshi event at the failure site. Rust's formatted message strings are Rust machinery and are not reproduced | the error code, class A; message text not compared |
| 2 | `ThermalGrid2D` stores `Vec<Vec<f64>>`; Crank–Nicolson allocates per step (roadmap backlog: breaking in Rust) | flat row-major buffer + dims; caller-reusable scratch; the JSON keeps Rust's nested-array shape | state dump flattened the same way, class B; JSON class A |
| 3 | lookups by `&str` name | table index (−1 when absent) plus a name search | same index set, class A |
| 4 | `register_agent` is `async` (tokio + reqwest) | synchronous `sandhi_http_post_opts` with the same 30 s timeout | request bytes (method, headers, body) class A against the oracle's captured request; the path is mapped by #6 |
| 5 | `EnvFilter` directives in `USHMA_LOG` (`ushma=debug,…`) | the level names trace / debug / info / warn / error, and `ushma=<level>`; other directive syntax is refused with a warning | parsing rows class A; log line format is sakshi's, not compared |
| 6 | `POST /v1/agents/register`, reads `agent_id` (matches no daimon that exists) | `POST /v1/agents`, reads `id`, as daimon 2.5.6 serves (P8) | the oracle's path is mapped to `/v1/agents`; response handling is checked against a stub answering as daimon 2.5.6 (`{"id":N}`, 401, 500, malformed, missing `id`) |
| 7 | reqwest's `timeout(30 s)` bounds the whole request | `total_ms` **and** `read_ms`, 30 s each: `total_ms` alone doesn't bound a server silent before its first byte (probe 7) | stub rows: a silent server gives the timeout error |
| 8 | reqwest sends exactly host, content-type, content-length, authorization, `accept: */*` | sandhi adds `User-Agent`, `Accept-Encoding: identity`, `TE: trailers`, `Connection: close`; the port sets `Accept: */*` | body bytes class A; the oracle's five headers by value (names case-insensitive); sandhi's four listed as extra (probe 7) |
| — | exponents written as rounded decimals (`powf(1.0/3.0)`, `powf(1.0/6.0)`) and similar exact-constant misses | fixed under the accuracy policy as each module ports (P2) | mapped per row, as kimiya's ADR 0002 #12 |

### 3.4 Not carried

Only Rust mechanics that have no meaning in Cyrius: `#[inline]`, `#[cold]`, `&'static` lifetimes,
and `Option<&'static T>`. Each one's effect is kept, as §7.3 shows. Nothing observable is dropped.

---

## 4. M0: prep, then initialize with `cyrius port`

`cyrius port` moves the Rust project into `rust-old/` and writes the Cyrius skeleton, so don't
hand-build the structure. On 6.7.18 the scaffold still overwrites project files (issue
`cyrius/docs/development/issues/2026-11-08-cyrius-port-clobbers-project-docs-and-omits-bench.md`,
**open**), so step 2 repairs what it writes.

0. **Freeze the oracle as it stands.** *Maintainer: git only, no code changes.*
   - Tag `1.4.0` at `9528e1f` (the release workflow requires the tag to equal `VERSION`).
     Publishing is optional.
   - The P(-1) checks are already clean on rustc 1.99.1 (header table). From here on, nothing in
     the Rust code changes.
   - Run `./scripts/bench-history.sh` on the tag. The last row in `bench-history.csv` is from
     `2f35058` (1.3.0).
   - Record glibc (**2.44** today); the oracle's transcendental bits depend on it.
   - Check crates.io reverse dependencies and record them in the header table. The CLAUDE.md
     consumer list (kiran, joshua) may be stale.
   - Confirm that every state the generator needs is reachable through the public API (§2.2). If
     something isn't, the generator extracts it from the Rust source (as the steam tables, §2.1).
     The Rust code is not edited.
1. **Scaffold.**
   ```sh
   cd /home/macro/Repos/ushma
   cyrius port --dry-run .
   cyrius port .
   ```
   Expect these to move into `rust-old/`: `Cargo.toml`, `Cargo.lock`, `rust-toolchain.toml`,
   `src/`, `tests/`, `benches/`, `examples/`, `Makefile`, `codecov.yml`, `deny.toml` and `target/`.
   Expect the scaffold to write `src/main.cyr`, `cyrius.cyml`, `tests/ushma.{tcyr,bcyr,fcyr}`,
   `docs/{adr,architecture,guides,examples}/` and `docs/development/state.md`.
   - **Track `rust-old/Cargo.lock`.** It's untracked today (the library-crate rule), but the
     oracle has to rebuild against exactly hisab 1.4.0, reqwest 0.13, serde_json 1.0.154 and the
     rest of its pins.
2. **Repair what the scaffold overwrites.** *(Claude; maintainer reviews.)*
   - **`CLAUDE.md`.** Merge `git show HEAD:CLAUDE.md` with the scaffold's Cyrius rules. Keep the
     identity, process, task sizing and DO NOT list. Translate the principles to their Cyrius
     mechanisms:
     - `#[non_exhaustive]` → append-only enum values;
     - `#[must_use]` → `#must_use` / `#pure`;
     - no `unwrap`/`panic!` → no `+?`, `-?`, `*?`, `assert_fatal`, `result_unwrap`, `unwrap` or
       `panic` in `src/`;
     - "Vec arena over HashMap" → sorted arrays and typed stride arenas; never `hashmap` for
       results (its seed is random per process);
     - "`write!` over `format!`" and "Cow over clone" → `str_builder` into caller-provided buffers,
       and borrowed pointers;
     - "tracing on all operations" → a sakshi span on every operation Rust instruments, and
       `log_*` events;
     - "Feature-gate optional deps" → the `[lib.ai]` profile for the HTTP/TLS stack;
     - new rules from the open issues (§2.8): bind a struct-returning call to a local before
       storing it into a field; never put an attribute before a `var` or `const`; grids and other
       stepped state are passed as `p: *T` and called with `&x`, because a by-value struct over
       8 B is a copy (6.6.16).

     Replace "Do not commit `Cargo.lock`" with the lockfile rule for `rust-old/`. Refresh the
     consumer list from step 0.
   - **Roadmap** restored from git, with its 2.0 section. **`state.md`** rewritten for 1.4.0 →
     2.0.0-dev.
   - **Workflows.** Keep the scaffold's Cyrius CI and add `test --bounds`, `deps --locked`,
     `audit` and a `[lib.ai]` profile job. Move the Rust workflows to
     `rust-old/.github/workflows/`. Append `/rust-old/target/` and `/build/` to `.gitignore`.
     `rust-old/` itself must stay tracked.
   - **Move Rust-only material** to `rust-old/`: `scripts/bench-history.sh`,
     `scripts/version-bump.sh`, `benchmarks.md`, `bench-history.csv`,
     `docs/architecture/overview.md`. Restore `rust-old/README.md` and add
     `rust-old/LINES_OF_SOURCE.txt`.
   - **ADRs** (scaffold's four-digit form):
     - `0001-cyrius-port-with-rust-oracle.md`: everything ports, stdlib- and hisab-first, P5, P6;
     - `0002-port-api-divergences.md`: §3.3;
     - `0003-ganita-gaussian-elimination.md`: P4 and the singularity-threshold divergence;
     - `0004-correct-rounding-is-the-accuracy-target.md`: P2 and P3;
     - `0005-serde-compatible-json.md`: P7.
   - Check that `cyrius.cyml` reads `version = "${file:VERSION}"`, not `0.1.0`, and that `bench`
     is in `[test] stdlib`. Document the scaffold's stub functions so `cyrius audit` passes.
     `VERSION` stays 1.4.0 until M8.
   - Copy kimiya's `handoff/gates/gates.sh` to `scripts/gates.sh`, adapted for ushma and both
     profiles. `handoff/` is untracked there and will be deleted.
3. **`cyrius.cyml`.** Confirm the keys with `cyrius help manifest`, and model it on prakash's
   (`[lib]` + `[lib.ai]`, sandhi's companions, the sakshi note). The migration notes govern:
   - The port template is shaped for a binary (`entry = "src/main.cyr"`, no `[lib]`), so add
     `[lib]` yourself. There is no `[lib]` until M1: `modules = []` is invalid.
   - sandhi, sakshi and bayan are **folds**: they go in `[deps] stdlib`, never `[deps.x]`. List
     every leaf the code calls (6.6.18). A fold's own Requires block brings its leaves (EM899),
     but sandhi's README asks consumers to list the companions, so list them.
   - Only hisab is a `[deps.X]`, with an https URL; 6.7.18 refuses `git://`/`http://`.
   - The manifest has no feature-gated dependencies, so the HTTP stack is always vendored. The
     `[lib.ai]` profile decides what a consumer links. sandhi's TLS stack is refused on the cx
     bytecode target, so `[lib.ai]` builds for ELF, PE and aarch64 only, the same set as the
     gates.
   - Use `[test]` (files, stdlib, timeout), `[deps.X] scope = "test"` where a dependency is
     test-only, `[build] dce = true`, and `cyrius deps --locked` in CI.
   ```toml
   [package]
   name = "ushma"
   version = "${file:VERSION}"
   license = "GPL-3.0-only"
   language = "cyrius"
   cyrius = "6.7.18"

   [deps]
   stdlib = ["syscalls", "string", "alloc", "str", "fmt", "io", "vec", "math", "ganita",
             "tagged", "result", "boxed", "overflow", "fnptr", "sakshi", "log", "bayan",
             "sandhi", "tls", "async", "random", "fdlopen", "dynlib"]   # sakshi warns twice: probe 10

   [deps.hisab]
   git = "https://github.com/MacCracken/hisab.git"
   tag = "3.4.1"
   modules = ["dist/hisab.cyr"]

   [test]
   stdlib = ["assert", "bench"]

   [lib]      # from M1, in dependency order; dist/ushma.cyr, sidecar `log bayan` (probe 10)
   modules = ["src/error.cyr", "src/constants.cyr", "src/validate.cyr", "src/mathx.cyr", ...]

   [lib.ai]   # [lib] + src/ai.cyr; dist/ushma-ai.cyr, sidecar `log bayan sandhi` (probe 10)
   modules = [...]
   ```
   - `[lib]` order: error → constants → validate → mathx → fmtf → json → logging → material_data →
     material → transfer → state_data → state → entropy → phase_data → phase → steam_data → steam →
     cycle → chem_data → chem → stat → numerical → bridge → soorat. `[lib.ai]` appends `ai`.
   - Module dependencies, from Rust's features: entropy and phase need state; steam needs phase;
     cycle needs state, entropy and steam; numerical needs transfer; chem needs state and
     entropy; stat needs state and transfer.
4. **Toolchain and math check** (§7.1). Confirm the stdlib counterpart of every Rust f64 method
   ushma calls:
   - **Rust NaN semantics.** `f64_min`/`f64_max`/`f64_clamp` (math.cyr) disagree with Rust when an
     operand is NaN (Rust ignores a NaN operand). Probe the exact behaviour; ushma has 11 `max`,
     6 `min` and 9 `clamp` calls. Write `_us_*` wrappers where it matters (1.4.0's
     "`f64::max` drops NaN" fix depends on it).
   - **`is_nan` / `is_finite`** (5 and 21 calls). math.cyr has no such fns and no
     `INFINITY`/`NAN`/`EPSILON` constants among its `F64_*` set, so define local helpers (`x != x`,
     `x - x == 0`) and constants.
   - **`powi`** (7 calls) becomes repeated multiplication in Rust's order.
5. **Golden-value generator** `scripts/oracle/` (copy kimiya's crate and adapt it):
   - `ushma-oracle` depends on `rust-old/` by path with `features = ["full"]`. Its `Cargo.lock` is
     seeded from `rust-old/Cargo.lock`.
   - Use TSV format 1 (f64 as bit patterns), deterministic output, a `--check` mode and glibc in
     the header.
   - Start with `constants` and the **math sweep**: every std method ushma calls, over special
     values and ushma's ranges, including the actual `powf` exponents. Those are the computed ones
     (γ − 1 and its expressions, a polytropic k, an interpolation fraction) and the literals 1.5,
     2.5, 0.8, 0.4, 1.6, 4.2, 0.75, 0.425, 0.725, 9/16, −8/27, 1/3 and 1/6.
   - **JSON rows:** for every value the corpus produces of a serde type, serde_json's document,
     byte for byte (P7).
   - **Request rows:** the generator runs a loopback listener that answers as daimon 2.5.6, and
     records the exact request `register_agent` sends (method, path, headers that matter, body).
     It also records the client's outcome for each scripted response: 201 with `id`, 401, 500, a
     malformed body, a missing `id`.
   - **Logging rows:** `USHMA_LOG` values → the resulting level.
6. **Table and constant generator** (same crate). Tables come through the oracle's public API or,
   for private data, from the Rust literals with the completeness proof (§2.1):
   - `constants.cyr`: `GAS_CONSTANT`, `ATM`, `STANDARD_TEMP`, `STANDARD_T`, `STEFAN_BOLTZMANN`,
     `BOLTZMANN_K`, `PLANCK`, `AVOGADRO`, `MAX_ADAPTIVE_SUBSTEPS`, the private ones (`T_CRIT`,
     `P_CRIT`, `R4_T_MIN`, `QUALITY_ROUNDOFF`, …) and the bridge constants (`KLEIBER_W`, …);
   - `material_data.cyr`: 10 materials;
   - `state_data.cyr`: the gas critical data and acentric factors;
   - `phase_data.cyr`: 8 substances;
   - `chem_data.cyr`: 12 species;
   - `steam_data.cyr`: 42 saturated rows, 72 superheated rows, 9 tier pressures, `R4_N`.

   Write baked typed column arrays and decimal literals in positional form. `probe_tables.py`
   compiles every cell and compares it bit for bit.
7. **Parity harness.** Copy kimiya's `tests/parity.tcyr`, `tests/parity_core.cyr` and
   `tests/ulp.cyr`:
   - per-module dispatchers;
   - outcomes ok / FAILED / measured / pending / accepted;
   - exit code = the number of failing modules;
   - the ADR 0002 mapping table.

   Add a JSON comparator (document bytes, and decode → values), impetus's state loader and step
   comparator for M6, and a request comparator for M7.
8. **Baseline measurements** → `docs/development/port-parity.md`. Run the math sweep on x86_64
   ELF, x86_64 PE (wine, one process at a time) and aarch64 (qemu). Measure these against MPFR:
   - the `_us_cr_*` primitives;
   - the stdlib alternatives (`_f64_exp_polyfill`, `_f64_ln_polyfill`, `_f64_atan_polyfill`,
     `ganita_f64_pow`, `ganita_f64_tanh`).

   That records exactly what P3 buys over the stdlib.
9. **Probes.** Each one is a small program with numbers, run before M1 writes code. They go in
   `scripts/probes/` and `port-parity.md`:
   1. **Callback into hisab.** `debye_cv`'s integrand passed to `calc_integral_simpson` as a
      function pointer + context. Compare its bits against Rust hisab 1.4.0 at 200 steps: equal
      means class B, different means class D.
   2. **Memory.** The global allocator is a bump with no free. Arenas (`arena_new_growable`,
      `arena_reset`) and `freelist` are the tools. Step a 1D grid 10⁵ times (explicit,
      Crank–Nicolson, adaptive) and solve a 2D steady state. The heap high-water mark must go
      flat after warm-up. Also measure the heap growth of one `calc_integral_simpson` call and
      one `calc_lerp` call (the `no-scoped-allocator-for-alloc` issue). Every steam lookup and
      every `debye_cv` makes them, so they have to be flat.
   3. **Wide results.** Return a cycle result: out-pointer vs handle vs `Result` payload, at the
      6.7.18 pin (the open small-struct issue).
   4. **Point lists.** `Vec<DiagramPoint>` as a vec of struct pointers vs parallel f64 vecs, timed
      at `cycle_ts_diagram`'s point counts.
   5. **Table lookup.** Steam lookups by T, by P and by (T, P) over baked columns, against Rust's
      timings.
   6. **Names** (P6). Include `kimiya/dist/kimiya.cyr` and a prakash ai bundle beside a ushma stub,
      and record the `duplicate fn` warnings and which definition wins. Then try `mod us;` in a
      `[lib]` module:
      - does it mangle fns only, or also globals, enums and structs?
      - does `distlib` keep the mangling?
      - does a consumer's `use us.f;` bind across a bundle include?
      - does `cyrius doc` document the result?

      If yes on all four, P6 uses `mod`. Otherwise, a written prefix.
   7. **sandhi POST.** `sandhi_http_post_opts` with a bearer header and a JSON body to a loopback
      stub. Check the request bytes against the oracle's captured request, the timeout and the
      status handling. Build `[lib.ai]` for all three targets.
   8. **Span cost.** `sakshi_span_enter`/`_exit` around `ideal_gas_volume` (1.6 ns in Rust) at a
      level that filters it out, vs no span, vs a `SAKSHI_LEVEL` compile-time ceiling. This
      decides whether spans are runtime-filtered or compiled out by default.
   9. **Gaussian elimination.** `ganita_mat_gaussian_elim` vs Rust hisab 1.4.0 on ushma's network
      matrices after 1.4.0's row scaling (bits and singularity decisions), and
      `ganita_mat_lu`/`_lu_solve` as the alternative (P4).
   10. **Manifest.** sakshi via hisab's `[deps.sakshi]` and via the stdlib leaf (log, sandhi):
       `cyrius deps` warnings, byte-identity at the pin, and `distlib` for both profiles.
   11. **JSON.** `bayan_json_v_parse_buf` on serde_json's documents for all 28 types (malformed,
       duplicate and unknown fields), and the vendored `fmtf.cyr` printer against serde_json's f64
       text and Rust's Display. Also run `#derive(Serialize)` on two of the types, to record
       exactly where its output and decoding differ from serde (P7's alternative).
   12. **f64 hygiene and silent miscompiles.** Build a small module with every trap from §7.3 and
       confirm that `cyrius test`/`lint` flag each one, so the warning-free gate catches them in
       real code. Add the open `s.f = g(t.e)` miscompile (§2.8), which nothing flags: record it and
       write the grid code to the workaround from the start.
10. **File upstream** (§8) before M1, so the requests carry M0's numbers.

**Gate:** `cyrius build` (both profiles) and `cyrius bench` pass; `cyrius test` runs the harness
(0 modules) with no warnings; `cyrius lint`, `fmt --check` and `audit` are clean. Committed: the
generator, `tests/ulp.cyr`, `src/mathx.cyr` with its verification, ADRs 0001–0005, the probe
results and the upstream issues.

---

## 5. Parity policy

tanmatra's classes apply, with kimiya's ADR 0006 reading of class C:

| Class | Tolerance | ushma examples |
|---|---|---|
| **A, exact** | bit-identical; equal integers, codes and bytes | constants, every table cell, `phase_at` classifications, enum properties, `US_ERR_*` codes, **JSON documents**, **request bytes**, `USHMA_LOG` parsing, `DaimonConfig`'s redacted text |
| **B, arithmetic** | ≤ 2 ulp, expected 0 | + − × ÷ √ only: Fourier, Newton cooling, resistances, Q = mcΔT, diffusivity, Biot, Stefan–Boltzmann (T⁴ by `powi`), ideal gas, Van der Waals, Kay's rule, Carnot, Helmholtz/Gibbs, Hess sums, IF97 Region 4 (algebraic + √), lerp interpolation, explicit FD steps, Gauss–Seidel sweeps, the network solve (ganita LU: 0 ulp, probe 9) |
| **C, transcendental** | **0 ulp against the corrected oracle** once a function uses the `_us_cr_*` primitives; rows where glibc misrounds are counted "accepted under ADR 0004" | LMTD, ε-NTU, fins, Nusselt correlations, view factors, adiabatic relations, isothermal work, entropy changes, Clausius–Clapeyron, vapor/sublimation curves, ln-P interpolation, K and ln K, van 't Hoff, Boltzmann, partition functions, Maxwell–Boltzmann, Einstein, the bridges' power laws |
| **D, iterative** | relative ≤ 1e-12 | only where an iteration path can differ: the hisab 3.4.1 Simpson (`debye_cv`: 1,942 of 3,000 rows bit-identical, ≤ 5 ulp, probe 1), and any solver path an ADR names |

- **A zero carries the oracle's sign** in every class (M4b). The ulp line puts +0 and −0 at one point, so the harness
  also compares a zero's sign bit. Without that check, ports of Rust's `sum` that started at +0.0 had passed.
- **Never widen a class to make a row pass.** A failure is one of three things:
  - a port defect;
  - a Cyrius stdlib or hisab defect (file it upstream);
  - a Rust oracle defect (fixed in the port, with an ADR 0002 row; the Rust code is never
    edited).

  The only exception is an ADR-accepted divergence, which the harness still reports by name.
- **The reference-value tests travel with the code.** They cover IAPWS-IF97 saturation to 1e-8,
  Incropera A.4 air, the Rankine 13.2 % case, the unit-square view factor of 0.2000, diamond's
  Debye Cv at 3 K and propane's ln K ≈ 837.
- **Test branch points on both sides.** Every switch gets golden rows just below and just above
  it:
  - steam tier boundaries and exact tier rows;
  - T_sat and the first tabulated row of each tier;
  - the quality clamp (1e-9);
  - `phase_at`'s curve joins and the sublimation branch;
  - the view-factor asymptote at X, Y < 0.01;
  - Debye's switch at Θ/T > 50, and Θ/T > 709;
  - Fo·(1 + Bi) = 0.5;
  - the adaptive sub-step cap;
  - Pitzer Z ≤ 0;
  - K leaving the f64 range;
  - the refrigeration Carnot limits;
  - the bridges' fallbacks for NaN inputs;
  - `muscle_work_to_heat` at η = 1.

---

## 6. Milestones

Every milestone ends with these gates:
- `cyrius build` for both profiles;
- `cyrius test` and `test --bounds` **with no warnings**;
- the full corpus on ELF, PE and aarch64 via `scripts/gates.sh`;
- `cyrius bench` against the Rust numbers;
- `vet`, `lint`, `fmt --check`, `doc --check` and `audit`;
- `coverage --exec` at ≥ 80 % (CLAUDE.md's target), with `coverage --min` gating CI;
- `fuzz --poison --bounds` on every reader of untrusted input: the 28 JSON decoders,
  `USHMA_LOG` parsing, and daimon response handling;
- `distlib --all` and `distlib --check`;
- `deps --locked`;
- **surface parity**: every Rust `pub fn` has a Cyrius counterpart, and the port has no
  exemptions;
- the milestone's range-edge sweep (P2), its sakshi spans and its JSON codecs;
- CHANGELOG and `state.md`.

| Milestone | Modules | Key work | Golden grid (minimum) |
|---|---|---|---|
| **M1 Foundation, logging & transfer** | error, constants, validate, mathx, fmtf, json, logging, material, transfer | `US_ERR_*` codes, names and fixed messages, failure detail to sakshi; `finite_gt`/`finite_ge`; correctly rounded primitives; the shortest-digit printer and JSON helpers; `us_log_init`/`_with_level` (`USHMA_LOG`); the span macro pattern every later module uses; 10-material table; 32 transfer functions (conduction, convection, radiation, resistance networks, fins, LMTD, ε-NTU, Re/Pr/Nu, Blasius, three view factors) | every material × field + JSON; each function over seeded physical grids plus NaN/±∞/0/negative edges; LMTD at ΔT₁ = ΔT₂ and 1 ulp apart; NTU 0…100 × C_r 0…1; view factors across the 0.01 switch; every `USHMA_LOG` form |
| **M2 Equations of state & entropy** | state, entropy | Ideal gas, Van der Waals, Redlich–Kwong, Peng–Robinson; virial; Pitzer–Abbott Z; Kay's rule; work and adiabatic relations; entropy changes; Carnot; Helmholtz/Gibbs; Clausius inequality; `GasData` JSON | every gas × T × V_m incl. near-critical; Pitzer Tr × Pr incl. Z ≤ 0; mixtures with zero/NaN/negative fractions; Carnot at T_c ≥ T_h |
| **M3 Phase & steam** (M3a phase, M3b tables + Region 4, M3c superheated) | phase, steam | substance table; `phase_at`'s piecewise curves and sublimation; latent heats (CO₂ sublimation); Clausius–Clapeyron; generated steam tables; exact IF97 Region 4 both ways; saturated by T/P; quality with the clamp; ln-P tier interpolation; `superheated_by_entropy`; JSON for `Phase`, `SubstanceData` and the three steam types | `phase_at` T × P per substance across every join; every steam row through the lookups; Region 4 over 273.16…647.096 K and back; between-row/-tier points; every tier through `superheated_by_entropy`; x = 1 + 2e-16 |
| **M4 Cycles** (M4a gas, M4b steam, M4c diagrams + comparison) | cycle | Otto, Diesel, Brayton, Rankine (s-inversion), refrigeration (negative net work), heat-pump COP, `compare_cycles`, T-s/P-v diagrams; results through out-pointers (probe 3); JSON for the 7 cycle types | r, γ, T grids; Rankine/refrigeration over the table range incl. 13.2 % and COP ≤ Carnot; every diagram's point list; `processes` shorter than `state_points` |
| **M5 Chemical & statistical** | chem, stat | species table and JSON; Hess; reaction Gibbs; K and `ln_equilibrium_constant`; van 't Hoff; flame temperature; Boltzmann; S = k ln W; equipartition; partition functions; Maxwell–Boltzmann; Einstein; Debye via hisab Simpson (function-pointer callback) | every species × function; propane ln K ≈ 837; partitions at hν ≪ kT; Debye over Θ/T ∈ (0, 50] ∪ (50, 10⁴] incl. diamond at 3 K |
| **M6 Numerical** | numerical | `ThermalGrid1D` explicit / Crank–Nicolson (reusable scratch) / adaptive / `max_stable_dt`; `ThermalGrid2D` flat storage, Gauss–Seidel with convective BCs; `ThermalNetwork` connectivity, row scaling, ganita LU (P4); malformed state → Err; JSON for the 5 types (nested 2D shape kept) | single-step rows for every BC kind × method incl. Fo·(1 + Bi) = 0.5; 2D steady states with iteration counts; networks incl. a floating node and R ≈ 1e14; the 10⁵-step soak (probe 2) |
| **M7 Integrations** | bridge, soorat, ai | 16 bridge conversions (bijli, kimiya, dravya, badal, pravash, sharira) with their NaN fallbacks; the 6 soorat visualization types (heatmap, profile, cycle diagram, network graph, flux field) with JSON; the daimon/hoosh configs (JSON, redacted text), `us_ai_daimon_client_new`, `us_ai_register_agent` on sandhi, in `[lib.ai]` | every bridge over its documented ranges and NaN inputs; soorat from every M3–M6 grid/network/cycle golden state, incl. empty and all-NaN fields; ai against the oracle's captured requests and every scripted daimon response |
| **M8 Release 2.0.0** | all | Rust-test ledger: every one of the 449 + 30 tests is ported or folded into golden rows, with no deferrals; full surface-parity diff (155/155); `docs/benchmarks-rust-v-cyrius.md` (63 benchmarks); **the 5 examples as Cyrius programs** in `docs/examples/` (`basic_thermo`, `heat_engine`, `building_thermal`, `cooling_system`, `phase_change`), built and run in CI; `cyrius distlib --all`; `cyrius audit`; `docs/guides/migration-2.0.md` (§3.3, P6's names, `[lib]` vs `[lib.ai]`); VERSION 2.0.0 | full corpus green on x86_64 ELF, x86_64 PE and aarch64, both profiles |

Order: M1 → … → M8. M6 is late because it carries the memory design and the stateful parity
layer, and M7 because its modules sit on everything else. M3 and M4 are large: each splits into
verified sub-parts, per CLAUDE.md task sizing.

---

## 7. Rust → Cyrius mapping

### 7.1 Math functions

Every std f64 method in non-test `src/`, 1.4.0, whole crate:

| Rust (calls) | Stdlib counterpart (6.7.18) | Port uses | Parity |
|---|---|---|---|
| `sqrt` (23), `abs` (30), `ceil` (1) | builtins `f64_sqrt`, `f64_abs`, `f64_ceil` | the builtins | exact |
| `exp` (14), `ln` (23) | builtins `f64_exp`/`f64_ln` (x87 on x86, polyfill on aarch64: bits differ by target, FAQ #9); `_f64_exp_polyfill`/`_f64_ln_polyfill` (≤1 ulp, same bits everywhere). The open `exp-polyfill-double-rounds-subnormal-results` issue: up to 39 ulp on aarch64 subnormal results, and x87 misrounds about 5 % | `_us_cr_exp`, `_us_cr_ln` (P3) | class C |
| `exp_m1` (2) | **none** in math or ganita | `_us_cr_expm1` (P3) | class C |
| `powf` (29) | `ganita_f64_pow` (note its `pow(−1, huge)` defect, tanmatra ADR 0007) | `_us_cr_pow`; exact-constant exponents reviewed: the correlation exponents (0.8, 0.4, 1/6, 9/16, −8/27, 1/3) stay Rust's doubles (ADR 0002 row 9) | class C |
| `atan` (11) | builtin `f64_atan`; `_f64_atan_polyfill` (≤1 ulp, same bits); `ganita_f64_atan2` | `_us_cr_atan2(x, 1.0)` (atan2 with x = 1 > 0 is atan exactly) | class C |
| `tanh` (2) | `ganita_f64_tanh`: **not correctly rounded** (52 of 860 golden rows off, port-parity.md §2b) | `_us_cr_tanh`, ushma's own (ADR 0004 amendment) | class C |
| `powi` (7) | — | repeated multiplication in Rust's order | class B |
| `max` (11), `min` (6), `clamp` (9) | `f64_max`, `f64_min`, `f64_clamp` (math.cyr) | wrappers with Rust's NaN semantics (step 4) | exact |
| `is_finite` (21), `is_nan` (5) | none | `x - x == 0`, `x != x` helpers | exact |

### 7.2 Libraries

| Rust dependency (use) | Cyrius | Parity |
|---|---|---|
| hisab 1.4.0 `calc::lerp` (steam ×3) | hisab `calc_lerp` (same `a + (b − a)·t` fast form; stdlib `f64_lerp` rounds identically) | class B |
| hisab 1.4.0 `calc::integral_simpson` (`debye_cv`, n = 200) | hisab `calc_integral_simpson(f, a, b, steps, out)` | class D (probe 1) |
| hisab 1.4.0 `num::gaussian_elimination` (`ThermalNetwork::solve`) | ganita `ganita_mat_lu` + `ganita_mat_lu_solve` (P4, ADR 0003 amendment) | class B, 0 ulp and the same refusals (probe 9) |
| serde + serde_json (28 types; `ai` bodies) | bayan (`bayan_json_parse`, `bayan_json_v_*`, `bayan_json_build`) + `str_builder` writers + vendored `fmtf.cyr` (P7) | class A |
| thiserror (`UshmaError` Display) | `US_ERR_*` codes, `us_error_name` / `us_error_message` (fixed text per code); the values at the failure site go to sakshi | the code, class A |
| tracing (`#[instrument]` ×28, all DEBUG) | sakshi `sakshi_span_enter(name, len)` / `sakshi_span_exit()` | span names and levels: reference tests |
| tracing-subscriber (`EnvFilter`, `USHMA_LOG`) | `log_init(level)`, `log_set_level`, io `getenv("USHMA_LOG")` | class A on parsing (§3.3 #5) |
| reqwest 0.13 (POST, JSON, bearer, 30 s, `error_for_status`) | sandhi `sandhi_http_post_opts(url, headers, body, len, opts)`, `sandhi_headers_new`/`_set`, `sandhi_http_options_*` (timeouts), `sandhi_http_status`/`_body`/`_err_kind` | request bytes class A |
| tokio (async runtime for `register_agent`) | none needed: sandhi's call is synchronous. Cyrius async (`CYRIUS_ASYNC=1`) runs coroutines only on x86_64 Linux, macOS and Windows; aarch64 refuses them (guide L8529–8675), so an async client wouldn't cover every target | §3.3 #4 |

### 7.3 Types and idioms

| Rust | Cyrius |
|---|---|
Each row cites the language guide (`cyrius/docs/guides/cyrius-guide.md`, "L") where it decides
the mapping.

| Rust | Cyrius |
|---|---|
| `Result<T, UshmaError>` + `?` | value-form `Result` (`Ok`/`Err` return a (tag, payload) register pair with zero allocation; bind both halves, `var t, v = f();`) + `?`, in fns declared `: Result` (L6995–7131, L7244–7294). The payload is one word: the code (§3.3 #1) |
| structs of f64 fields | `struct S { cp: f64; … }`, with typed fields doing float arithmetic (L1676–1753). Every field must be initialized, and a struct literal is legal only as a declaration's initializer (L1514–1573) |
| `pub const` struct instances (`WATER`, `GAS_NITROGEN`, …) and `ALL_*: &[&T]` slices | `const` holds only an integer, f64 or string (L7737–7772). An instance becomes a top-level `var` with constant fields, which is baked into the image (L8711–8761). Array initializers refuse struct elements (L429–432), so the tables are generated typed **column arrays** plus named index constants, and `ALL_*` is an index range |
| `Phase`, `ProcessKind`, `CycleKind`, `Side` | integer enums with `US_`-prefixed variants (P6); variant names are global (L391–396); `match` with `_ =>` where Rust's `#[non_exhaustive]` mattered (L7451–7465) |
| `BoundaryCondition` (`Convective { h, t_inf }`) | a struct with a kind tag and `h`, `t_inf` fields. A `: stack` enum carries one payload word, and a boxed payload enum allocates on every construction, never reclaimed (L7009–7051) |
| cycle result structs | caller out-pointers (`p: *T`). A struct over 16 B returns through a hidden pointer, and a by-value struct returned through a fn pointer held in a variable faults (L1944–1966, L2270–2276). Probe 3 |
| `Vec<f64>`, `Vec<DiagramPoint>` | `vec` handles (`vec_get` bound as `var x: f64 = …`, L128), a vec of struct pointers, or parallel f64 vecs (probe 4) |
| `Vec<Vec<f64>>` 2D field | flat row-major `*f64` buffer + `nx`, `ny` (`a[i][j]` is refused, L324–536; §3.3 #2) |
| the integrand closure passed to hisab | a top-level fn + context pointer. Closure parameters are untyped words that refuse `: f64` (rebind them), and a capturing closure heap-allocates its environment on every construction (L7574–7703) |
| `String`, `&str`, `format!`, `write!` | `Str` handles and cstrings, `str_builder`, `fmt` (L3919–4079). There is no format and no public f64 text API, hence `fmtf.cyr` |
| `Option<String>` (`api_key`) | a cstring or 0 |
| custom redacting `Debug` | `us_ai_daimon_config_fmt` printing `<redacted>`. There is no Debug derive |
| `Default` impls | `us_*_default()` constructors (`default` itself is a reserved word, L696–790) |
| `impl` methods | `impl T { fn f(self, …) }` → `T_f`, or plain prefixed fns (L672–694) |
| `usize` counts, checked arithmetic | i64 (there is no usize, L199–212), with explicit range checks. `+? -? *?` **exit the process with code 57** (L1074–1135), so they are banned in `src/` |
| 7-parameter fns (`fin_rectangular_heat`) | fine: 7+ parameters go on the stack (5+ on Windows) (L547) |
| `#[must_use]` / `#[non_exhaustive]` | `#must_use` (warns on a statement-level discard) and `#pure` (checked, L5422–5483) / append-only enum values |
| `#[inline]` / `#[cold]` | `#inline` where a bench shows it helps; it is ignored on fns with defaults, defer or async (L887). There is no `#cold` |
| doc comments, `missing_docs`, rustdoc `-D warnings` | `#` doc lines; `cyrius doc --check` fails an undocumented public fn (L6677–6682) |

**Language traps.** The f64 ones are **warnings, not errors**, so the build succeeds. That's why
the gate is a warning-free `cyrius test`:
- **An operator is typed by its left operand.** `0 - x` and `2 * x` are integer operations on an
  f64's bits; write `0.0 - x` or `x * 2.0` (L145–167).
- **An integer constant stored in an f64 slot keeps its bits**, and that includes arguments:
  `f(3)` passes a subnormal, so write `3.0` (L120–128).
- **`==` and `<` on builtin results compare bit patterns.** Use `f64_lt`/`f64_eq`, or declare the
  result `: f64` (L104–118). `if (x)` and `!x` test the raw bits (L7851).
- **No exponent literals.** `1e-9` is written `0.000000001` (L53–60). A hex literal in a `const`
  is an integer.
- **`f64_round` ties to even.** Rust's `round` ties away from zero.
- **Large local arrays.** A local array over about 120 KB silently gets static, shared storage
  (L2797–2811).
- **Limits:** at most 256 distinct includes (L5076); a 2 MiB string pool (L5811); a 24 MB
  source, includes counted.
- **Gate quirks (from the siblings):** `cyrius audit` rejects comment lines over 120 characters,
  and `cyrius bench` stops any file after 300 s.

---

## 8. Upstream (file at M0, then as found)

**Filed 2026-11-11** (uncommitted in each target, for the maintainer): [`upstream/README.md`](upstream/README.md)
lists the 11 new issues (each with a repro that reproduced from its own repository), the 7 addenda to issues other
ports filed, and the note for kimiya.

| To | Request | Evidence |
|---|---|---|
| ganita (cyrius stdlib) | correctly rounded exp, expm1, ln, log10, pow, atan2, tanh (and the rest kimiya has), so ports stop vendoring them | three private copies (impetus `crmath.cyr`, kimiya `mathx.cyr`, ushma `mathx.cyr`); step 8's MPFR measurements |
| ganita | an `expm1` (and `log1p`) at all | ushma `partition_vibrational`, `debye_cv`; kimiya's decay chains |
| bayan | a serde_json-compatible shortest f64 printer | impetus ADR 0010 and ushma P7: two copies of `fmtf.cyr` |
| cyrius stdlib (math) | `f64_is_nan`, `f64_is_finite`, `F64_INF`, `F64_NAN`, `F64_EPSILON`; Rust-semantics min/max | 26 ushma call sites; check impetus's filing first and add to it |
| hisab | drop `[deps.sakshi]` from `cyrius.cyml`: the bundle makes no `sakshi_` call, and the package makes every consumer that uses log, sandhi or sakshi resolve with a warning | §2.7, probe 10: every resolve warns, ushma's and each consumer's; the copies are byte-identical; patra did the same at 1.13.0 |
| ganita | an allocator-taking `ganita_mat_lu_solve_a` (or a caller scratch): `_lu_solve` allocates its forward-substitution vector on every call and never frees it | probe 9: 320 B per call at n = 40 with every other buffer reused. Until it lands, ushma runs the substitution itself over ganita's factors, measured bit-identical on 600 networks (M6c; ADR 0003, second amendment) |
| ganita | `ganita_mat_gaussian_elim`'s threshold scales with the right-hand side, so it refuses well-posed systems whose b is large | probe 9: 259 of 2,000 row-scaled network systems that LU and Rust hisab 1.4.0 solve |
| sandhi | `total_ms` doesn't bound a read that blocks before the first response byte (nor does `idle_ms`); only `read_ms` does | probe 7: a server silent for 3 s answers 200 under `total_ms` 1000 |
| bayan | keep the sign of the integer `-0` (serde_json reads it as −0.0 into an f64; the tree has +0) | probe 11: 20 documents; an addendum to pravash's issue, filed the same evening |
| bayan | refuse a raw control character (U+0000–U+001F) inside a string, as RFC 8259 and serde_json do | probe 11: 10 documents parsed that serde_json refuses |
| cyrius stdlib (fdlopen, dynlib) | six raw syscalls with no PE route (`fdlopen.cyr:472, 530, 626, 954`; `dynlib.cyr:674, 950`) warn on every PE build that links sandhi | probe 10: every `[lib.ai]` consumer's PE build keeps `dynlib.cyr:674` after DCE |
| sakshi | a level-taking span entry (`sakshi_span_enter_at`): `sakshi_span_enter` / `_exit` never test the level | probe 8: 192 ns per filtered-out span into the ring buffer, 9.7 ns level-gated, around an 8.1 ns function |
| cyrius | `mod NAME;` in an included file stays in force to the end of the unit, so it mangles the includer's `main` | probe 6 |
| cyrius | `#derive(Serialize)`'s `T_to_json` faults on a 0 `Str`, which a document missing that field decodes to | probe 11 |
| cyrius | `cyrius port` still clobbers project docs (existing issue) | it did at 6.7.18: items 1–4 and pravash's `Cargo.lock` point (addendum) |
| cyrius | the open `field-store-of-struct-call-copies-a-field-argument` miscompile | nothing to add: probe 12 re-ran the issue's own repro (7 of 12 forms wrong), and ushma's grid code is written to the workaround |
| cyrius | a duplicate non-private `fn` across bundles is a warning, and the last one wins | probe 6; ushma × kimiya (10 names), × prakash (`ai_*`) |
| cyrius | `#inline` declines a small body for one `if` or one local ("body has control flow or a var declaration"); filed 2026-11-12 as an addendum to pravash's `#inline` issue | M1d's review: three range checks 21 → 6 ns once rewritten as one expression; the level-gated spans, 7 ns a call, have no expression form |
| kimiya, pravash, impetus, tanmatra (**later**: not filed, as those ports are in progress) | their `tests/parity_core.cyr` compares floats by ulp distance alone, and `tests/ulp.cyr` puts +0 and −0 at one point, so a zero of the wrong sign passes. Once each port is done, check it with ushma's one-line zero-sign check in `_par_compare` | M4b: ushma's copy of the same harness hid 5 wrong-sign rows (`sum` starting at +0.0). None of the four uses ushma's `0.0 - 0.0` start, so what they would find is unknown |
| cyrius stdlib (string) | `memcpy` and `memset` copy one byte per iteration. This is already open (ranga 2026-11-02, with impetus's small-copy and pravash's numbers); ushma's numbers are not yet added | M4a: 160 B (four state points) took 274 ns, against 45 ns for a word loop. ushma copies words itself |
| kimiya | ganita's LU reproduces Rust hisab 1.4.0's elimination bit for bit on row-scaled systems, refusals included, so its planned transcription (ADR 0003) may not be needed; `mat_gaussian_elim` refuses more | probe 9 |
| kimiya | its `KM_GAS_CONSTANT` is 8.314462618, truncated at 1e-11 relative; ushma 1.4.0 uses the exact N_A·k_B = 8.31446261815324 | values crossing the ushma ↔ kimiya bridges disagree in the 10th digit |

---

## 9. After 2.0.0

- Switch `mathx.cyr` and `fmtf.cyr` to upstream versions when §8 lands, then re-measure.
- The roadmap backlog, demand-gated.
- Rust 1.4.0 stays frozen in `rust-old/` as the parity reference. It gets no fixes and no
  releases.

---

## 10. Tracking

| Milestone | Status | Golden rows | Reference tests | Benches | Notes |
|---|---|---|---|---|---|
| M0 Prep & initialize | 🟡 part 1 done 2026-11-11 (steps 0–3: tag verified, Rust bench baseline, `cyrius port`, repairs, manifest, ADRs 0001–0005; gates green on ELF, `--bounds`, PE, aarch64) | 15,738 (constants 9, math 15,687, table presets 42) | — | Rust baseline: 63 | Parts 2a (math check, generator, validate.py), 2b (table data, probe_tables.py), 2c-1 (parity harness: 3,374 ok, 0 failed) 2c-2a (MPFR corrected oracle: 26 rows), 2c-2b (mathx vendored) and 2c-3 (CR tanh, cr_verify: every golden row gated, 15,712 ok + 26 accepted) done. Part 3a (probes 2, 3, 4, 5, 6, 8, 12) done; P6 settled as a written prefix. Part 3b (probes 1, 7, 9, 10, 11) done: Simpson class D, P4 moved to LU (class B), the ai client sets `read_ms`, both bundles link on three targets, P7 stands with the ADR 0005 amendment's decoder rules. Upstream filed (step 10, §8: 11 issues with verified repros, 7 addenda, a kimiya note; `upstream/README.md`), uncommitted in the targets. Next: M1 |
| M1 Foundation, logging & transfer | ✅ done 2026-11-12 (M1a–M1c 2026-11-11, M1d 2026-11-12) | fmtf 26,088; material 2,927; transfer 15,100 (32 fns; `transfer.cr.tsv` 8) | foundation 120, material 54, transfer 260, integration 19; Rust tests 102 of 102 due (96 unit, 6 integration) | 21 (material 4, span, transfer 16) | Parts: **M1a** error (codes, names, messages; detail to sakshi), validate, fmtf (vendored), logging (`USHMA_LOG`, level-gated spans), the `[lib]` / `[lib.ai]` profiles; **M1b** json helpers, material + its JSON, the first fuzz target; **M1c** transfer 1 (conduction to boundary layers, 18 fns); **M1d** transfer 2 (fins, exchangers, view factors, 14 fns), the review (inlined finite checks, powi as products: diffusivity 33 → 16 ns, radiation 50 → 25), coverage 100 % of public fns executed, `scripts/surface_parity.py` and `scripts/rust_test_ledger.py` in the gates and CI, the spans under Rust's span names |
| M2 Equations of state & entropy | ✅ done 2026-11-12 | state 11,136 (presets 9, 21 fns, GasData JSON 96; `state.cr.tsv` 1); entropy 3,454 (9 fns; `entropy.cr.tsv` 1) | state 181, entropy 85, integration 54 (14 of M1–M2); Rust tests 180 of 180 due (state 47, entropy 23, integration 8) | 14 (state 8, entropy 6) | Parts: **M2a** state 1 (ideal gas, van der Waals, the works, adiabatic, Z, reduced coordinates; 10 fns), `UsGas` and its JSON, the GasData fuzz target; **M2b** state 2 (Redlich–Kwong, Peng–Robinson, virial, Pitzer Z, Kay's rules; 11 fns), the 6 state integration tests; **M2c** entropy (9 fns), the 2 entropy integration tests, the review (nothing to change: inlined checks, no allocation), surface parity 36/36 and 9/9, coverage 79/79 |
| M3 Phase & steam | ✅ done 2026-11-12 | phase 27,333 (presets 9, 10 fns, JSON 109); steam 31,113 (10 fns, JSON 197; `steam.cr.tsv` 10) | phase 171, steam 169, integration 74; Rust tests 275 of 275 due (phase 44, steam 46, integration 5 of M3's) | 14 (phase 6, steam 8) | Parts: **M3a** phase (`UsPhase`, `UsSubstance`, the curves, `phase_at`, latent heats, Clausius–Clapeyron, both JSON types, the two fuzz targets); **M3b** steam: Region 4 both ways, the saturated lookups, quality and wet steam, both JSON types, the review (R4 coefficients as generated consts: 174 → 141 ns); **M3c** superheated (the lookup, the by-entropy inverse, `SuperheatedEntry`), the Rankine integration test, surface 29/29 and 19/19, coverage 117/117 |
| M4 Cycles | ✅ done 2026-11-12 | cycle 5,607 (13 fns, JSON 312; `cycle.cr.tsv` 3) | cycle 192, integration 85; Rust tests 339 of 339 due (cycle 59 of 59, integration 24) | 5 (cycle/otto, diesel, brayton, rankine, refrigeration) and 3 Cyrius-only (two diagrams, a comparison) | Parts: **M4a** the four types and their JSON, the process helpers, Otto, Diesel and Brayton, `US_ERR_OUT_OF_MEMORY`, ADR 0002 row 10 (the result filled in place: 746 → 510 ns), 4 fuzz targets; **M4b** Rankine (ADR 0002 row 11), refrigeration with `UsRefrigerationResult`, `heat_pump_cop`, the superheated inversion's entropy-only probe, the zero-sign harness rule and the six `sum` starts it found; **M4c** the comparison and both diagrams into caller buffers (ADR 0002 row 12), their two types, the review (`push_process` refuses an unknown code), cycle PORTED (surface 44/44), coverage 166/166 |
| M5 Chemical & statistical | ✅ done 2026-11-12 | chem 2,340 (presets 13, 7 fns, JSON 100; `chem.cr.tsv` 3); stat 4,170 (14 fns; class D 352; `stat.cr.tsv` 11) | chem 114, stat 107, integration 100; Rust tests 401 of 401 due (chem 23, stat 33, integration 30 of 30) | 7 (chem 3, stat 4) | Parts: **M5a** chem, `UsSpecies`, term lists (ADR 0002 row 13), surface 24/24; **M5b** stat, Debye on hisab's Simpson (class D, 0–2 ulp), σ (ADR 0002 row 14), surface 16/16, coverage 192/192 |
| M6 Numerical | ✅ done 2026-11-12 | numerical 2,211 (13 fns, JSON 256) | numerical 347; Rust tests 428 of 428 due (numerical 27 of 27) | 5 (numerical/explicit_1d, crank_nicolson_1d, adaptive_1d_4substeps, gauss_seidel_2d_20x20, thermal_network_5) | Parts: **M6a** boundary conditions (ADR 0002 row 15), `Side`, the 1D grid and its three steps, scratch reuse, the review (rolling loops: 822 → 591 ns); **M6b** the 2D grid (flat nodes, ADR 0002 row 2, its mis-shaped documents accepted by name), Gauss–Seidel, the size guards; **M6c** the network (ADR 0002 row 16: its lists and solve scratch in `UsNetwork`, a system that overflows f64 refused, 9 rows accepted by name), LU with the `_lu_solve` stopgap measured bit-identical (ADR 0003, second amendment), the review (949 → 829 ns), the 10⁵-step soak, numerical PORTED (surface 36/36), coverage 242/242 |
| M7 Integrations | 🟡 M7a done 2026-11-12 | bridge 6,439 (16 fns; `bridge.cr.tsv` 3) | bridge 60; Rust tests 459 of 459 due (bridge 31 of 31) | 16 Cyrius-only (bridge) | Parts: **M7a** bridge, the 16 conversions with their fallbacks in 1.4.0's order, the constant review (ADR 0002, under the table: no miss), surface 16/16; **M7b** soorat's 6 types and their JSON; **M7c** the ai client in `[lib.ai]`; `total_ms` + `read_ms`, serde's body bytes (probe 7) |
| M8 Release 2.0.0 | ☐ | | | | 155/155 fns, 479/479 tests, 5 examples |
