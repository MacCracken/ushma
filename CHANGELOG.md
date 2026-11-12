# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

The Rust → Cyrius port has begun; it ships as **2.0.0** with every module, feature, test,
benchmark and example ported. Plan of record: `docs/development/ushma-port-plan.md`, with
ADRs 0001–0005.

### Changed
- **project** — the Rust 1.4.0 crate (tag `1.4.0`) is frozen in `rust-old/` as the parity
  oracle. It is never edited and never re-released, and every fix lands in the Cyrius port
  (ADR 0001). Moved there with it: `Cargo.toml`, `Cargo.lock` (now tracked, so the oracle
  rebuilds against its exact pins), `src/`, `tests/`, `benches/`, `examples/`, `Makefile`,
  `deny.toml`, `codecov.yml`, `rust-toolchain.toml`, the Rust CI workflows,
  `scripts/bench-history.sh`, `scripts/version-bump.sh`, `bench-history.csv`, `benchmarks.md`,
  `docs/architecture/overview.md` and the 1.4.0 README
- **project** — scaffolded with `cyrius port` on Cyrius 6.7.18: `cyrius.cyml` (stdlib leaves
  incl. ganita, sakshi, log, bayan and sandhi; hisab 3.4.1), Cyrius CI and release workflows,
  `src/main.cyr` smoke entry, test/bench/fuzz stubs
- **docs** — CLAUDE.md merged for Cyrius (Rust principles translated, new rules for the
  toolchain's open issues); README, CONTRIBUTING, getting-started guide and `state.md`
  rewritten for the port; ADRs 0001–0005 (port with a frozen oracle, API divergences, ganita
  Gaussian elimination, correct rounding, serde-compatible JSON)

- **CLAUDE.md** — a new rule. `p + n` on a typed pointer (a `p: *T` parameter, or `&s` of a struct-typed local) steps `n · sizeof(T)` bytes, so a byte offset goes through an untyped copy (`var q = p;`). M4a's first build wrote state points at `out + 8·k` on a `*UsStatePoint`, which gave a SIGSEGV and wrong points. The fuzz harness's `fz_substance` read the same way, out of bounds, and is fixed
- **transfer** — `us_lmtd_parallel` and `us_lmtd_counter` return `Err(v)` / `Ok(v)` from their helper's result instead of the pair `(t, v)`. The pair form worked where it was the only return, but cyrius warns that beside an `Err(..)` path it is read as a single value, so the port uses one form throughout (found writing `superheated_by_entropy`)
- **transfer** — `us_conduction` and `us_radiation` now open their spans as `conduction` and `radiation`, Rust's span names (CLAUDE.md); M1c had named them `us_conduction` and `us_radiation`
- **scripts/gates.sh, CI** — the gates now run `cyrius coverage --exec --min 80`, `scripts/surface_parity.py` and `scripts/rust_test_ledger.py --check` (both were SKIP "from M1"). CI gains `distlib --check`, `fuzz --poison --bounds`, the coverage gate and both parity scripts
- **tests/parity_core.cyr** — the report's name column is 43 wide and always ends in a space, so a 42-character name (`transfer::view_factor_perpendicular_plates`) no longer runs into its class
- **ADR 0002 row 9** — a correlation's fitted exponents (Dittus–Boelter's 0.8, 0.4; Churchill–Chu's 1/6, 9/16, −8/27; Pohlhausen's 1/3) stay the doubles Rust wrote, computed correctly rounded; an exponent that is a physical constant is restored exactly. The plan had planned `cbrt` for Pohlhausen's 1/3
- **ADR 0002 row 1** — errors are Cyrius-native: `US_ERR_*` codes in the value-form `Result`, `us_error_name` / `us_error_message` with a fixed text per code, and the failing values to a level-gated sakshi event at the failure site. thiserror's formatted messages are not reproduced, and parity is on the code. `UshmaError` derives no serde, so it has no JSON codec (the plan had listed one)
- **ADR 0003** (amended, probe 9) — `ThermalNetwork::solve` uses `ganita_mat_lu` + `ganita_mat_lu_solve`, not `ganita_mat_gaussian_elim`: LU is bit-identical to Rust hisab 1.4.0 on all 2,000 network systems 1.4.0 solves and refuses exactly the 178 it refuses, while elimination refuses 259 that 1.4.0 solves. `debye_cv`'s Simpson is class D (probe 1: ≤ 5 ulp from 1.4.0's)
- **ADR 0005** (amended, probe 11) — the JSON decoders' rules beyond bayan: refuse non-finite numbers, invalid UTF-8 and raw control characters as serde_json does; refuse duplicate keys and deep nesting in unknown fields and integers ≥ 2^63 (mapped rows); the integer `-0` is an open divergence until bayan keeps its sign
- **ADR 0002** (rows 7, 8, probe 7) — the ai client sets `read_ms` as well as `total_ms` (30 s each), and sends serde's body bytes and `Accept: */*`; sandhi's four transport headers are extra
- **scripts/gates.sh** — adds `cyrius fuzz --poison --bounds` and `vendor_fmtf.py --check`
- **scripts/gates.sh** — the run now names every step that failed and exits 1: a nonzero exit, a test log without a clean parity or assert summary, an ushma warning line, a fmt / lint / doc failure, or a staged file. It used to print each step's exit code and end 0, which hid a broken `--check` once the oracle crate gained its probe binaries (now `default-run = "oracle"`)
- **scripts/oracle** — serde_json is built with `float_roundtrip`: its default float parse reads 35,970 of 160,320 of its own shortest texts back 1–2 ulp off (probe 11). Writing is unchanged, and the golden corpus is unaffected

### Added
- **M6a, boundary conditions and the 1D grid** (port plan §6): `src/numerical.cyr`'s first part.
  - **Types:** `UsBoundary` holds the condition kinds (`US_BOUNDARY_FIXED`, `_INSULATED`, `_CONVECTIVE`), and `UsBoundaryCondition` holds a condition's value, filled by `us_boundary_fixed` / `_insulated` / `_convective` (ADR 0002 row 15). `UsSide` holds the four sides. `UsGrid1D` is the 1D transient grid.
  - **Functions:** `us_grid1d_new` / `_new_a`, `us_grid1d_fourier_number`, `us_grid1d_max_stable_dt`, `us_grid1d_step_explicit`, `us_grid1d_step_crank_nicolson` and `us_grid1d_step_adaptive`, with accessors for the nodes (`us_grid1d_len`, `_node`, `_set_node`) and the two conditions (`us_grid1d_left`, `_right`). All three types have JSON.
  - **Behaviour:** Rust's operation order, with unary minus as `_us_neg`. Every step checks the grid as Rust's `check` does, so a grid edited by hand or read from JSON is refused, not trusted. The explicit and Crank–Nicolson steps open spans under Rust's names.
  - **Memory:** the grid keeps its Crank–Nicolson scratch, sized on the first step, so no step allocates after it (3,000 steps, 0 heap bytes; ADR 0002 row 2). A refused allocation is `Err(US_ERR_OUT_OF_MEMORY)`.
  - **Golden rows:** `tests/golden/numerical.tsv` (new, 1,190 rows: B 1,073, A 117; 336 refusals). A step's row runs from the oracle's own state to the next, so every step is checked bit for bit:
    - `new` over node counts, lengths, diffusivities, initial temperatures and every condition at its edges;
    - the three steps on seven setups (fixed, insulated and convective ends, Biot numbers from 0 to 10, 3 to 50 nodes), along trajectories;
    - fractions and multiples of the stable step, and dt at the edges;
    - the explicit stability limit (bisected) ± 2 ulps, the adaptive exact multiples ± 1 ulp, and its 10⁶ sub-step cap (bisected) ± 2 ulps;
    - Crank–Nicolson's pivots at dt = 1e300 (`DivisionByZero`);
    - malformed grids;
    - the three types' JSON.
    - Every row passes, 5 JSON rows under ADR 5.
  - **`validate.py`** models all six from their raw tokens and reproduces their 1,073 rows.
  - **Tests:** `tests/numerical.tcyr` (83 assertions) ports 15 of Rust's 27 numerical tests, every one due at M6a. Reference tests cover the spans, the struct layout, the scratch's reuse, refused allocations and the documents. The fuzz harness reads `BoundaryCondition`, `Side` and `ThermalGrid1D` documents (964 read and round-tripped).
  - **Benches:** Rust's three 1D benchmarks, recorded against Rust's.
  - **The review:** each loop loads every node once now, rolling the old values forward, and Crank–Nicolson's boundary rows sit outside the interior loop. Every row is unchanged. `explicit_1d` went 822 → 591 ns, `crank_nicolson_1d` 2,665 → 2,104 and `adaptive_1d_4substeps` 1,739 → 1,295 (the minimum of the runs; the machine was loaded).
  - **`_us_vec_size`** moved from `src/cycle.cyr` to `src/validate.cyr`, since the grids use it too.
  - **`scripts/surface_parity.py`:** the three types are in TYPES, and `BoundaryCondition` maps to its kind enum and its value struct. numerical stays scheduled until M6c; probed as ported, it lacks only M6b's and M6c's 13 items (23 of 36)
- **M5b, stat and the M5 gates** (port plan §6): `src/stat.cyr`, all 14 functions, so M5 is complete. They are the Boltzmann factor and population ratio, S = k·ln W, both equipartition forms, the three partition functions, the Maxwell–Boltzmann density and three speeds, and the Einstein and Debye heat capacities.
  - **Behaviour:** Rust's operation order, with unary minus as `_us_neg`. exp, expm1, ln and pow are correctly rounded: the vibrational partition function's −expm1(−x) keeps its precision where hν ≪ kT, as 1.4.0's does. `debye_cv` integrates on hisab's `calc_integral_simpson` (200 intervals, a plain integrand function; Rust's closure captures nothing) and allocates nothing (1,000 integrals, 0 heap bytes). Rust's `u32` symmetry number is an integer whose negatives are refused (ADR 0002 row 14).
  - **Golden rows:** `tests/golden/stat.tsv` (new, 4,170 rows: B 1,350, C 2,468, D 352; 219 refusals). Every function is seeded and run at the edges. Every switch runs either side:
    - Einstein's x = 500 ± 1 ulp, and where (eˣ − 1)² drops below 1e-30 (bisected) ± 2 ulps;
    - Debye's Θ/T at 0.01 and 50 ± 1 ulp, a sweep across the integrated range, and 1.4.0's reference cases;
    - the vibrational denominator where hν ≪ kT, and where it underflows to 0 (`DivisionByZero`);
    - σ at 0, 1, 2, 3, 12 and u32::MAX.
  - **Parity:** every row passes. `debye_cv` is class D (relative ≤ 1e-12, probe 1). Its 352 rows are 0–2 ulp from 1.4.0 (280 bit-identical), since hisab 3.4.1's Simpson sums in another order and its expm1 is correctly rounded. glibc misrounds an expm1 inside 10 vibrational rows and an exp inside 1 Boltzmann row (`tests/golden/stat.cr.tsv`, new).
  - **`validate.py`** models all 14, with Rust hisab 1.4.0's Neumaier-compensated Simpson transcribed for `debye_cv`, and reproduces all 4,170 rows bit for bit. `Prim` gains expm1, and `model_args` reads `i:` tokens.
  - **Tests:** `tests/stat.tcyr` (107 assertions) ports all 33 of Rust's stat tests. Reference tests cover the switches, a negative σ, the denominator's underflow and the Debye heap. The 3 stat integration tests are ported, so all 30 of Rust's integration tests are (100 assertions); the ledger has 401 of 479.
  - **Benches:** Rust's four stat benchmarks, recorded against Rust's.
  - **The M5 gates:** stat is PORTED in both parity scripts (16 of 16). Coverage is 192 of 192
- **M5a, chem** (port plan §6): `src/chem.cyr`, all 7 functions. `UsSpecies` holds a `Species`, with `us_species_preset` / `us_species_name` over the generated table and its JSON. The functions are `us_reaction_enthalpy`, `us_reaction_gibbs`, `us_gibbs_from_enthalpy_entropy`, `us_ln_equilibrium_constant`, `us_equilibrium_constant`, `us_vant_hoff_k` and `us_adiabatic_flame_temperature`.
  - **Shapes** (ADR 0002 row 13): a reaction side is a term list, n word pairs of (coefficient, `*UsSpecies`) with its length beside it, the flat form of Rust's `&[(f64, &Species)]`.
  - **Behaviour:** Rust's operation order. Every sum starts at −0.0, and Rust's unary minus (`-ΔG°`, `-(ΔH°/R)`, `-ΔH°rxn`) is `_us_neg`, so ln K at ΔG° = +0 is −0.0 as in 1.4.0. exp is correctly rounded. `equilibrium_constant` and `adiabatic_flame_temperature` open spans under Rust's names. The flame temperature checks the reactants' coefficients before the products', as Rust's chained iterator does.
  - **Golden rows:** `tests/golden/chem.tsv` grows from the 13 preset rows to 2,340 (A 113, B 1,563, C 664; 380 refusals):
    - Hess's law and the Gibbs sum over named reactions, their reverses, empty sides, 200 seeded reactions, and each coefficient and species of methane's combustion at the edges. Edge species are written by value: NaN, ±∞, −0, f64::MAX, subnormal.
    - ΔH − TΔS, ln K and Van 't Hoff, seeded and at the edges, ΔG° ±0 among them.
    - K and Van 't Hoff where exp leaves the f64 range, bisected, ± 2 ulps.
    - The flame temperature over named combustions and initial temperatures, every coefficient at the edges, edge species, an empty product side, and where T_flame crosses 0 (bisected) ± 2 ulps.
    - `Species`'s JSON.
    - Every row passes. glibc misrounds an exp inside 3 `equilibrium_constant` rows (`tests/golden/chem.cr.tsv`, new), and 5 JSON rows are accepted under ADR 5.
  - **`validate.py`** models all 7 (the term-list functions from their raw tokens) and reproduces their 2,227 rows.
  - **Tests:** `tests/chem.tcyr` (114 assertions) ports all 23 of Rust's chem tests, with propane's ln K ≈ 837 among them. Reference tests cover the spans, ln K's signed zero and the preset refusals. chem's 3 integration tests are ported (93 assertions in `tests/integration.tcyr`). The fuzz harness reads `Species` documents: `fz_gas` became `fz_named` for both name-and-three-numbers types.
  - **Benches:** Rust's three chem benchmarks, recorded against Rust's.
  - **Surface parity:** chem is PORTED, 24 of 24, with `ALL_SPECIES` mapped to `US_SPECIES_COUNT` and `us_species_preset`. The ledger has 365 of 479
- **M4c, comparison, diagrams and the M4 gates** (port plan §6): `us_compare_cycles` / `us_compare_cycle`, `us_cycle_ts_diagram`, `us_cycle_pv_diagram` and `us_cycle_diagram_len`, with `UsCycleComparison` and `UsDiagramPoint` and both types' JSON. cycle is 13 of 13, 44 of its 44 public items, and M4 is complete.
  - **Shapes** (ADR 0002 row 12): a diagram writes `UsDiagramPoint` records into the caller's buffer and returns Ok(the count), or Err(`US_ERR_INVALID_PARAMETER`) when they would not fit. The check divides, so no point count overflows it. `us_compare_cycles(kinds, results, t_highs, t_lows, n, out)` takes Rust's tuples as four parallel lists and writes n records; `us_compare_cycle` does one.
  - **Behaviour:** Rust's operation order, with ln and pow correctly rounded. The T-s isochores and isobars of an ideal gas follow T_a·(T_b/T_a)^f, and the P-v isentropes the polytropic through both state points. A refrigeration result is compared with the Carnot COP. Non-physical limits give 0 and 0.
  - **Golden rows:** `tests/golden/cycle.tsv` grows to 5,607. M4c adds 434 (A 115, B 75, C 244):
    - the comparison on every cycle's results and on relabelled ones, each limit at the edges and T_H at T_L ± 1 ulp;
    - both diagrams at 0, 1, 2, 3, 5, 10 and 50 points a process, on every cycle's results and on results edited to reach each path's guard: short and empty process lists, isothermal and isenthalpic processes, γ 0 and NaN, non-positive and non-finite states, and volumes equal or 1 ulp apart;
    - both types' JSON.
    - Every row passes. glibc misrounds a pow inside one T-s and one P-v row (`cycle.cr.tsv` has 3 rows now), and 10 JSON rows are accepted under ADR 5.
  - **`validate.py`** models the comparison and both diagrams from their raw tokens (`RAW_ARGS`) and reproduces all 319 rows.
  - **Tests:** `tests/cycle.tcyr` (192 assertions) ports all 59 of Rust's cycle tests and checks a too-small buffer's refusal. The last 3 cycle integration tests are ported (85 assertions in `tests/integration.tcyr`): every Rust test due through M4, 339 of 479. The fuzz harness reads `DiagramPoint` and `CycleComparisonEntry` documents.
  - **Benches:** Rust benches neither the diagrams nor the comparison. The record has them as Cyrius-only rows: a 200-point T-s diagram takes 16.5 µs and a P-v diagram 17.6 µs (100 correctly rounded pows each, about 146 ns apiece), and one comparison 20 ns.
  - **The review** (work loop step 5): `us_cycle_push_process` now refuses a code that is not a `US_PROCESS_*` (Err(`US_ERR_INVALID_PARAMETER`)). Rust's enum holds no such value, and it would write a document that never reads back.
  - **The M4 gates:** cycle is PORTED in `scripts/surface_parity.py` (44 of 44) and `scripts/rust_test_ledger.py`. Coverage is 166 of 166 public functions
- **M4b, the steam cycles** (port plan §6): `us_rankine_cycle`, `us_refrigeration_cycle` and `us_heat_pump_cop`, with `UsRefrigerationResult` and its JSON. cycle has 10 of its 13 functions.
  - **Rankine:** Rust's `Option<f64>` superheat is a flag and a value (`superheated`, `t_superheat`; ADR 0002 row 11). The boiler exit is the saturated vapor or the superheated state. The turbine exhaust is wet, or superheated by entropy inversion where s₃ > s_g at the condenser.
  - **Refrigeration:** `UsRefrigerationResult` nests its `UsCycleResult` as the first field and adds the two COPs. `us_refrigeration_result_init` / `_init_a` make one, and `us_refrigeration_result_cycle` gives the cycle's address for the cycle accessors. Its `net_work` is Rust's `-w_comp`, a sign flip (`_us_neg`, new).
  - **Behaviour:** Rust's operation order. A steam lookup's refusal is returned as Rust's `?` returns it. Both cycles open spans under Rust's names, around the inversion's own.
  - **Golden rows:** `tests/golden/cycle.tsv` grows to 5,173. M4b adds 2,334 (A 53, B 213, C 2,068; 1,140 refusals):
    - Rankine: 1,482 rows over condenser × boiler pressures × superheats. Each tier temperature and each tier's saturation temperature runs ± 1 ulp, the boiler at the condenser pressure ± 1 ulp, the table's ends ± 1 ulp, and the superheat where the exhaust leaves the wet region (bisected) ± 2 ulps.
    - Refrigeration: 586 rows, with the condenser at the evaporator ± 1 ulp and the evaporator where the compressor exit enters the superheated table (bisected) ± 2 ulps.
    - `heat_pump_cop`: 213 rows. `RefrigerationResult`'s JSON: 53 rows.
    - Every row passes. 5 JSON rows are accepted under ADR 5, and glibc misrounds none.
  - **`validate.py`** models the three functions (`model_args` reads `n:` as `None`) and reproduces their 2,281 rows.
  - **Tests:** `tests/cycle.tcyr` grows to 162 assertions. It ports 15 more Rust tests, 45 of 59 in all (every one due), with the 13.2 % turbine-exhaust and the 677 K / 2.877 COP reference cases among them. The 2 steam-cycle integration tests are ported (82 assertions in `tests/integration.tcyr`). Reference tests cover a refrigeration result's document, the steam spans, reuse (1,000 analyses of each, 0 heap bytes) and refused allocations. The fuzz harness reads `RefrigerationResult` documents (618 read and round-tripped).
  - **Benches:** Rust's `cycle/rankine` and `cycle/refrigeration`, recorded against Rust's.
  - **`scripts/surface_parity.py`:** `RefrigerationResult` is in TYPES. Probed as ported, cycle has 35 of its 44 items, and the 9 missing are M4c's.
- **M4a, the cycle types and the ideal-gas cycles** (port plan §6): `src/cycle.cyr`'s first 15 functions and its first 4 types.
  - **Types:** `UsProcess` and `UsCycle` hold Rust's `ProcessKind` and `CycleKind` (`US_PROCESS_ISENTROPIC` …, `US_CYCLE_OTTO` …, append-only), with `us_process_name` / `us_cycle_name`. `UsStatePoint` holds a `StatePoint`, and `UsCycleResult` a `CycleResult`. All four have JSON.
  - **Functions:** the process helpers (`us_isochoric_heat`, `us_isobaric_heat`, `us_adiabatic_work`, `us_adiabatic_temp_from_pressure`) and the ideal-gas cycles (`us_otto_cycle`, `us_diesel_cycle`, `us_brayton_cycle`).
  - **Behaviour:** Rust's operation order, with ln and pow correctly rounded. The three cycles open spans under Rust's names.
  - **Results:** the caller owns a `UsCycleResult`, and each analysis fills it (ADR 0002 row 10). Its state points and processes are two stdlib vecs on the result's own allocator, sized and stored in place, so a reused result allocates nothing after its first analysis (1,000 analyses, 0 heap bytes). `us_cycle_state` / `us_cycle_process` read a result; `us_cycle_push_state` / `us_cycle_push_process` build one by hand.
  - **Golden rows:** `tests/golden/cycle.tsv` (new, 2,839 rows: A 144, B 747, C 1,948; 395 refusals).
    - The helpers on seeded grids and at their edges.
    - The three cycles over r, γ and T grids, with Diesel's cutoff ratio at r ± 1 ulp and Brayton's T₃ at T₂ ± 1 ulp.
    - 144 JSON rows.
    - Every row passes. glibc's `log` misrounds inside one Otto row (`tests/golden/cycle.cr.tsv`), and 10 JSON rows are accepted under ADR 5's mapped divergences.
  - **`validate.py`** models all 7 functions and reproduces their 2,695 rows.
  - **Tests:** `tests/cycle.tcyr` (104 assertions) ports 30 of Rust's 59 cycle tests, every one due at M4a, and checks each refusal's code. M4b and M4c bring the rest. Reference tests cover:
    - a result's document;
    - the three spans;
    - a reused result;
    - every allocation refused in turn.
  - **Fuzzing:** the harness reads `StatePoint`, `ProcessKind`, `CycleKind` and `CycleResult` documents (691 results read and round-tripped in 50,000 inputs). Its targets are grouped by shape (`fz_enum`, `fz_f64`), and its buffer is 2,048 bytes for the 1.4 KB result documents.
  - **Benches:** Rust's three ideal-gas cycle benchmarks, recorded against Rust's.
  - **The review:** storing the state points in place cut Otto from 746 to 510 ns. 24 `vec_push` calls had cost 266 ns, and the stdlib's byte-loop `memcpy` was no faster (274 ns for 160 B, against 45 ns for a word loop; cyrius's open `memcpy` issue).
  - **`US_ERR_OUT_OF_MEMORY`** (new, `src/error.cyr`): code 18, `"OutOfMemory"`, for a refused allocation where Rust would abort the process. The cycle analyses return it, and so do the result's init, its push accessors and its decoder. No library code calls `vec_push`, which aborts on a refused grow.
  - **`src/json.cyr`:** `_us_jr_unit_enum` reads a unit enum through its name fn, and phase's decoder uses it too.
  - **`scripts/surface_parity.py`:** the four cycle types are in TYPES. cycle stays scheduled until M4c. Probed as ported, it lacks only M4b's and M4c's 15 items (29 of 44 present).
- **M3c, superheated steam and the M3 gates** (port plan §6): `us_superheated_lookup` (filling a `UsSuperheated`: within a tier linearly in T from the saturated-vapor anchor; between tiers h and s linearly in ln P and v log-linearly) and `us_superheated_by_entropy` (the exact inverse for the table's model: its breakpoints sorted, bisected, the segment solved), with `SuperheatedEntry`'s JSON. steam is 10 of 10, and M3 is complete.
  - **Behaviour:** Rust's operation order, ln and exp correctly rounded; every `?` of Rust's propagates the inner refusal; `superheated_by_entropy` opens its span under Rust's name. Where Rust's `SuperheatOutOfRange` carries the requested (T, P), the refusal's DEBUG event does (`T=…, P=…`).
  - **Golden rows:** `tests/golden/steam.tsv` grows to 31,113 (B 13,880, C 17,036, A 197; 8,806 refusals). The lookup runs at every tier pressure ± 1 ulp and between tiers, at every table temperature ± 1 ulp, the midpoints and around each pressure's saturation temperature (10,999 rows); `superheated_by_entropy` at those pressures over every breakpoint's entropy ± 1 ulp, the midpoints and past the range (6,037). Every row passes; glibc misrounds a ln or exp inside 10 of them (`tests/golden/steam.cr.tsv`: 6 lookups, 4 inversions), and the port matches MPFR there.
  - **`validate.py`** parses the tier pressures and the 72-row table itself and reproduces all 17,036 superheated rows.
  - **Tests:** `tests/steam.tcyr` (169 assertions): all 46 of Rust's steam tests, the span test among them. `simple_rankine_cycle_energy_balance` is ported, so all 5 of M3's integration tests are. The fuzz harness reads `SuperheatedEntry` documents too.
  - **Benches:** the three superheated benchmarks, against Rust's.
  - **The M3 gates:** surface parity checks steam (19 of 19), the ledger has 275 of 479 ported (every one due), coverage is 117 of 117
- **M3b, steam: Region 4, the saturated table, quality** (port plan §6): `src/steam.cyr`'s first 8 functions: `us_saturation_pressure` and `us_saturation_temperature` (IAPWS-IF97 Region 4 both ways), `us_saturated_by_temperature` / `us_saturated_by_pressure` (filling a `UsSaturated`), the three `us_quality_from_*` and `us_wet_steam_properties` (filling a `UsWetSteam`), with both types' JSON.
  - **Behaviour:** Rust's operation order (class B: Region 4 is algebraic with √); between table rows, hisab's `calc_lerp` in T, with v_g through P·v_g/T as Rust interpolates it; quality's 1e-9 round-off clamp; results through the caller's struct (Ok(0)).
  - **Golden rows:** `tests/golden/steam.tsv` (new, 14,025 rows: B 13,880, A 145; 2,478 refusals). Region 4 over 273.15 … 647.096 K (1,500 seeded, every table row ± 1 ulp, the limits ± 2 ulps) and back at every resulting pressure ± 1 ulp; both lookups at every row ± 1 ulp and between rows; quality at the liquid and vapor values, past them by the tolerance ± 2 ulps and at x = 1 + 2e-16, the critical row (v_g = v_f, h_fg = 0) included; 16 + 61 `SaturatedEntry` and 13 + 55 `WetSteamProperties` JSON rows. Every row passes (10 JSON rows under ADR 5's mapped divergences).
  - **`validate.py`** parses the table and R4's coefficients from `rust-old/src/steam.rs` itself (`load_steam`) and reproduces all 13,880 function rows.
  - **Tests:** `tests/steam.tcyr` (111 assertions): 32 of Rust's 46 steam tests (the superheated 14 come with M3c). Three steam integration tests are ported early (`clausius_clapeyron_vs_steam_table`, `quality_consistency_across_methods`, `phase_lookup_agrees_with_steam_table`). The fuzz harness reads both steam documents.
  - **Benches:** five, against Rust's.
  - **`src/json.cyr`:** the codec of a struct of f64 fields only (`_us_jw_f64_struct` / `_us_jr_f64_struct`, by a field-name fn), which the steam types use and `SuperheatedEntry` and cycle's types will.
- **M3a, phase** (port plan §6): `src/phase.cyr`. `UsPhase` holds Rust's `Phase` (`US_PHASE_SOLID` … `US_PHASE_SUPERCRITICAL_FLUID`, append-only) with `us_phase_name` and its JSON; `UsSubstance` holds `SubstanceData`, with `us_substance_preset` / `us_substance_name` over the generated table and its JSON. Rust's five `SubstanceData` methods are `us_substance_sublimes_at_atmospheric`, `_vapor_pressure`, `_sublimation_pressure`, `_latent_heat_sublimation` and `_phase_at` (Ok(US_PHASE_*)); the free functions are `us_heat_of_fusion`, `us_heat_of_vaporization`, `us_heat_for_phase_change` and both Clausius–Clapeyron forms.
  - **Behaviour:** Rust's operation order; ln and exp correctly rounded. `phase_at` and `heat_for_phase_change` open spans under Rust's names, and Rust's span field (`substance = <name>`) is a DEBUG event inside the span (`_us_span_enter_kv`, new in `src/logging.cyr`), since sakshi's spans take no fields.
  - **Golden rows:** `tests/golden/phase.tsv` grows to 27,333 (A 232, B 2,389, C 24,712; 6,054 refusals). The 8 presets and 106 edge substances (each property at 0, −1, NaN, ±∞; a boiling point outside the triple–critical span; sublimers either side of their triple point) run the curves across every join (the triple, melting, boiling and critical points ± 1 ulp), and `phase_at` runs 10,797 points, each preset at Rust's own curve pressure ± 1 ulp at every join temperature. 4 `Phase::to_json`, 22 `Phase::from_json` (serde's string and `{"Variant":null}` forms), 15 `SubstanceData::to_json` and 68 `from_json` documents. Every row passes; glibc misrounds none of them, and 5 `SubstanceData` documents are accepted under ADR 5's mapped divergences.
  - **`validate.py`** models all 10 functions (`model_args` passes a substance's `s:` name through) and reproduces all 27,215 of their rows.
  - **Tests:** `tests/phase.tcyr` (171 assertions): all 44 of Rust's phase tests, with the refusal codes and the span reference test (the substance field included). The integration test on material and phase (`water_boiling_point_consistency`) is ported. The fuzz harness reads `SubstanceData` and `Phase` documents too, and its seed buffer is 512 bytes and clamped (the 300-byte Water seed overran the old 192, caught by the poisoning allocator).
  - **Benches:** Rust's six phase benchmarks, against Rust's.
  - **`_us_finite_ge_any`** (new, `src/validate.cyr`): `_us_finite_ge` for a bound read from the caller's data, exact at a −∞ bound (`vapor_pressure` checks T against the substance's own triple point). `_us_finite_ge` keeps the faster two-comparison form for its 115 literal bounds, about 1 ns a check cheaper; `_us_finite_gt` is exact for every bound as it is.
  - **`scripts/surface_parity.py`:** phase is checked (29 of 29), and a public enum's variants must be members of its Cyrius enum (`Phase::Gas` → `US_PHASE_GAS` in `UsPhase`; a renamed member fails it)
- **docs/development/upstream/** — an addendum to cyrius's open `#inline` issue (pravash's `2026-11-11-inline-refuses-small-numeric-helpers.md`), filed 2026-11-12 into the cyrius tree (uncommitted) with a repro that exits 1: `#inline` also declines a 2-parameter, 12-token body for one `if` or one local. It costs ushma's range checks 21 ns per three where one expression costs 6, and its level-gated spans 7 ns a call, with no expression form. Ask: inline `if`, early returns and locals within the size limit
- **M2c, entropy and the M2 gates** (port plan §6): `src/entropy.cyr`'s 9 functions: the ideal-gas and isothermal entropy changes, heat transfer's entropy, the Carnot efficiency and refrigeration COP, the Helmholtz and Gibbs energies, the Clausius check and the entropy of mixing.
  - **Behaviour:** Rust's operation order, refusals as `US_ERR_*`, spans for `carnot_efficiency` and `entropy_of_mixing`. `us_is_spontaneous` returns 1 or 0, and `us_entropy_of_mixing(moles, xs, n)` reads the caller's list.
  - **Golden rows:** `tests/golden/entropy.tsv` (new, 3,454 rows): a seeded sweep and every argument at each edge, T₁ = T₂ and V₁ = V₂ ± 2 ulps, Carnot with T_cold at, within 2 ulps of and above T_hot (at 300 K, 1 K, 1e-300 and 1e300), and 344 mixtures (zero, −0, NaN, negative and out-of-range fractions, subnormal fractions, sums across the 1e-10 tolerance). Every row passes: class B and A bit-identical, class C 0 ulp against MPFR. glibc's `log` misrounds inside one ideal-gas entropy row (`tests/golden/entropy.cr.tsv`).
  - **`validate.py`** models all 9 and reproduces all 3,454 rows.
  - **Tests:** `tests/entropy.tcyr` (85 assertions): all 23 of Rust's entropy tests, with the refusal codes and the span reference test. `tests/integration.tcyr` gains the 2 integration tests on entropy (54 assertions in all); every integration test due at M2 is ported.
  - **Benches:** Rust's six entropy benchmarks, against Rust's.
  - **The M2 gates:** surface parity checks entropy (9 of 9), the ledger has 180 of 479 ported (every one due), and coverage is 79 of 79 public functions
- **M2b, state part 2** (port plan §6): the other 11 state functions, so all 21 are ported: Redlich–Kwong's parameters and pressure, Peng–Robinson's parameters and pressure, Pitzer's second virial coefficient, the two- and three-term virial pressures, the Pitzer–Abbott compressibility, and Kay's rules for T_c, P_c and ω.
  - **Behaviour:** Rust's operation order; refusals as `US_ERR_*`, Pitzer–Abbott refusing Z ≤ 0 as 1.4.0 does; spans for `redlich_kwong_pressure` and `peng_robinson_pressure`. The parameter functions return their pair or triple as multiple values (`var a, b = us_redlich_kwong_params(tc, pc);`). Kay's rules take each list with its own length (`us_kays_rule_tc(ys, ny, tcs, nv)`), so a mismatch stays a refusal. T_c^2.5 and Pitzer's T_r^1.6 and T_r^4.2 are the correctly rounded pow on the exponents Rust writes (ADR 0002 row 9).
  - **Golden rows:** `tests/golden/state.tsv` grows to 11,136. The cubic equations run every gas × 11 temperatures across T_c × V_m from b ± 2 ulps to 1 m³/mol; Pitzer–Abbott runs T_r × P_r out to Z ≤ 0 for five ω; Kay's rules run 233 mixtures (zero, −0, NaN, negative and out-of-range fractions, unequal and empty lists, sums across the 1e-10 tolerance, non-finite values). Every row passes: class B bit-identical, class C 0 ulp against MPFR. glibc's `pow` misrounds T_r^4.2 in one Pitzer–Abbott row, the first entry of `tests/golden/state.cr.tsv`.
  - **`validate.py`** models all 21 state functions (`model_args` now reads a row of two lists) and reproduces all 11,031 function rows.
  - **Tests:** `tests/state.tcyr` has 181 assertions: all 47 of Rust's state tests, and the span test covers all four instrumented functions. `tests/integration.tcyr` ports the 6 integration tests that call only state (47 assertions). The ledger has 155 of 479 ported, every one due.
  - **Benches:** Redlich–Kwong and Peng–Robinson, against Rust's.
  - **Surface parity:** state is checked now, 36 of 36
- **M2a, state part 1 and the reference gases** (port plan §6): `src/state.cyr`'s first 10 functions: the ideal gas law (pressure, volume, temperature), van der Waals, isothermal and isobaric work, the adiabatic temperature, the compressibility factor and the reduced coordinates. `UsGas` holds `GasData`, with `us_gas_preset` / `us_gas_name` over the generated table and `us_gas_to_json` / `us_gas_from_json_a` / `us_gas_from_json` (serde's document, ADR 0005).
  - **Behaviour:** Rust's operation order, refusals as `US_ERR_*`, spans for `ideal_gas_pressure` and `van_der_waals_pressure` under Rust's names; ln and pow correctly rounded.
  - **Golden rows:** `tests/golden/state.tsv` grows to 4,250: a seeded sweep and every argument at each edge for each function, van der Waals with V − nb across 0, γ at 1 and its neighbours, V₁ = V₂ ± 2 ulps, every gas against a grid of T and P; 18 `to_json` and 78 `from_json` documents (Nitrogen's variants among them). Class B is bit-identical, class C 0 ulp against MPFR (glibc misrounds none of these rows), and 5 JSON rows are accepted under ADR 5, the same mapped divergences as material's.
  - **`validate.py`** models the 10 functions and reproduces all 4,145 of their rows; the presets stay `check_tables`' and the JSON rows the harness's (`unmodelled`).
  - **Tests:** `tests/state.tcyr` (101 assertions): 28 of Rust's 47 state tests, each with its refusal's code, and the span reference test; the other 19 come with M2b. The fuzz harness reads `GasData` documents too (2,247 of 50,000 inputs accepted and round-tripped).
  - **Benches:** six, against Rust's (`docs/benchmarks-rust-v-cyrius.md`, M2).
  - **`scripts/oracle/src/sweep.rs`** (new): the argument sweeps transfer's rows used, now shared with state's (`edges`, `tuples`, `res`, `row`, `around`); `transfer.tsv` is unchanged.
  - **`scripts/surface_parity.py`:** `GasData`'s names (`UsGas`, `us_gas_*`, the presets' `US_GAS_*`), and a mapping may now precede its module's milestone. State stays scheduled until M2b; with it counted, the 11 part-2 functions are all that is missing.
  - **`_US_TINY`** (the 1e-30 divisor threshold) moves from `src/transfer.cyr` to `src/validate.cyr`, now that state uses it
- **M1d, transfer part 2 and the M1 gates** (port plan §6): the 14 remaining transfer functions, so all 32 are ported. They are the fin parameter, rectangular-fin heat, efficiency and effectiveness; parallel- and counter-flow LMTD, the LMTD heat rate, NTU, both ε-NTU effectivenesses and the ε-NTU heat rate; and the parallel-plate, perpendicular-plate and coaxial-disk view factors.
  - **Behaviour:** Rust's operation order, refusals as `US_ERR_*`, and spans for the five Rust instruments (`fin_rectangular_heat`, `lmtd_parallel`, `lmtd_counter`, `view_factor_parallel_plates`, `view_factor_perpendicular_plates`). tanh, ln and atan are the correctly rounded `_us_cr_tanh`, `_us_cr_ln` and `_us_cr_atan2(x, 1)`; π is the stdlib's `F64_PI`, Rust's `consts::PI` bit for bit. A view factor's refusal reports the first dimension that isn't positive (ADR 0002 row 1).
  - **Golden rows:** 7,883 more, 15,100 in `tests/golden/transfer.tsv`. The grid adds each branch switch: LMTD at ΔT₁ = ΔT₂, 1–2 ulps apart and across the 1e-10 band; NTU 0…100 × C_r 0…1 with C_r at both switches; view factors across 0.01; fin efficiency at m·L = 1e-30. Class B is bit-identical. Class C is 0 ulp against the corrected oracle, and `transfer.cr.tsv` gains the 4 rows where glibc misrounds (3 `exp`, 1 `atan`).
  - **`validate.py`** models all 32 transfer functions (`Prim` gains glibc's and MPFR's `log`, `tanh`, `atan`) and reproduces all 15,100 rows.
  - **Tests:** `tests/transfer.tcyr` has 260 assertions: all 67 of Rust's transfer tests, plus the refusal codes and a reference test that each of the seven spans opens under Rust's name and closes on every path. `tests/integration.tcyr` (new) ports the 6 Rust integration tests that call only transfer and material. `tests/foundation.tcyr` ports Rust's 18 error tests one for one: each refusal's DEBUG event names its subject and value, as Rust's Display text did. `_us_fail_s`, a text-valued refusal, is new for the reason and name fields.
  - **Benches:** fin heat, counter-flow LMTD and effectiveness, the disk and perpendicular view factors, recorded against Rust's.
  - **`scripts/surface_parity.py`** (new): every public item of each ported module has its ADR 0001 Cyrius name (error 19, logging 2, material 16, transfer 34), and so does every public constant. Items with no single name are mapped to the names that carry them, and those are checked too: there are no exemptions. Every other Rust file is listed with its milestone.
  - **`scripts/rust_test_ledger.py`** (new): each Cyrius test names the Rust tests it ports (`# rust <file>: <test>`), and the ledger (`docs/development/rust-test-ledger.md`, `--check` in the gates) fails on a due test with no port. 102 of the 479 are ported, every one that is due; the other 377 are counted under their milestones.
  - **Coverage:** `cyrius coverage --exec`: 44 of 44 public functions executed (100 %)
- **M1c, transfer part 1** (port plan §6): `src/transfer.cyr`'s first 18 functions. They cover conduction, convection and radiation; conduction, convection, series and parallel resistances; heat stored; diffusivity; Biot; lumped capacitance; Reynolds, Prandtl and Nusselt; Dittus–Boelter; Churchill–Chu; and the Blasius and thermal boundary layers.
  - **Behaviour:** each computes in Rust's operation order, with refusals as `US_ERR_*` and conduction and radiation in level-gated spans. The correlation exponents are Rust's doubles, through the correctly rounded `pow` (ADR 0002 row 9).
  - **Golden rows:** `tests/golden/transfer.tsv` holds 7,217 rows, a seeded physical sweep plus every argument at each edge (0, −0, −1, NaN, ±∞, the 1e-30 divisor threshold and its neighbours, 1e-300, 1e300). Every class-B row is bit-identical. The class-C rows are 0 ulp against the corrected oracle, and `tests/golden/transfer.cr.tsv` lists the 4 rows where glibc's `pow` misrounds (accepted under ADR 4).
  - **`validate.py`:** it now models each transfer function in Rust's order, with glibc's or MPFR's primitives. It reproduces all 7,217 rows and writes each module's corrected oracle (`check_corrections` is per module now).
  - **Tests and benches:** `tests/transfer.tcyr` has 100 assertions (Rust's 32 tests for these functions, each also checking the refusal's code); 11 transfer benchmarks are recorded against Rust's
- **M1b, JSON and materials** (port plan §6):
  - `src/json.cyr`: the shared serde-compatible helpers (ADR 0005 and its amendment). Writers cover keys, escaped strings and f64. Readers work on bayan's tree on the caller's allocator: object or exact-length sequence form, fields by name or position, f64 (refusing `null` and ±∞), UTF-8-checked strings, and a pre-scan that refuses a raw control byte inside a string;
  - `src/material.cyr`: `UsMaterial`, with the name as bytes and length so a decoded `\u0000` survives. It adds `us_material_preset` / `us_material_name` over the generated table, `us_material_diffusivity`, `us_material_volumetric_heat_capacity`, `us_material_to_json`, and `us_material_from_json_a` / `us_material_from_json`, which read the tree and name on the caller's allocator, so a decode on a reset arena never grows the heap;
  - golden rows from the oracle (`m_material.rs`, with `json_cases.rs` now shared with probe 11): 1,410 per method, 19 `to_json` documents, and 77 `from_json` documents (serde's own, Copper's 60-odd variants, and hand-written edges). All pass, 5 of them accepted under ADR 5: exactly the mapped divergences (duplicate keys in an unknown field, deep nesting, an integer `-0` and one ≥ 2^63). The harness accepts them only when the port did what the ADR says;
  - `tests/material.tcyr` (54 assertions: Rust's 8 material tests, presets, UTF-8 and control-byte refusals, the sequence form, a 1,000-decode heap check at 0 bytes); the fuzz harness reads material documents too and checks each accepted one round-trips (3,308 of 50,000 inputs);
  - `tests/ushma.bcyr`'s first real benchmarks and `docs/benchmarks-rust-v-cyrius.md`: the material methods at 7 ns (harness floor 4 ns) against Rust's 0.25 ns inlined, and the Cyrius-only JSON and span costs
- **M1a, the foundation** (port plan §6; the `[lib]` modules so far):
  - `src/error.cyr`: `US_ERR_*` codes, one per Rust 1.4.0 `UshmaError` variant in Rust's order, plus `US_ERR_DESERIALIZE`; append-only. `us_error_name` gives the Rust variant name (the golden corpus's `E:` token) and `us_error_message` a fixed message (ADR 0002 row 1);
  - `src/validate.cyr`: `_us_finite_gt` / `_us_finite_ge`, which refuse NaN and ±∞, and the `_us_fail*` refusals. A refusal returns its code and, when DEBUG is on, reports the function, the error and the failing value in one sakshi event, rendered into a stack buffer (no allocation);
  - `src/fmtf.cyr`: impetus's shortest-digit printer (Rust's Display and serde_json's form), vendored mechanically by `scripts/vendor_fmtf.py` (`--check` in the gates) and pinned by `tests/golden/fmtf.tsv`: 26,088 rows over 13,044 doubles, every one ok. `validate.py` checks the same rows independently (read-back, shortest, closest, form);
  - `src/logging.cyr`: `us_log_init` / `us_log_init_with_level` (USHMA_LOG: a level name or `ushma=<level>`, any case, `off`; anything else warns and keeps the default), `us_log_parse`, and the level-gated `_us_span_enter` / `_us_span_exit` every instrumented function uses;
  - `[lib]` and `[lib.ai]` in `cyrius.cyml`, with `dist/ushma.cyr` and `dist/ushma-ai.cyr` from `cyrius distlib --all` (sidecar `str log io`). `src/main.cyr` includes the same modules;
  - `tests/foundation.tcyr` (68 assertions: codes, predicates, refusals and their events, every USHMA_LOG form, init, spans) and `tests/ushma.fcyr`, the fuzz harness, with USHMA_LOG parsing as its first target (50,000 seeded inputs under `--poison --bounds`)
- **scripts/oracle** — the golden-value generator: drives `rust-old/` read-only and writes `tests/golden/constants.tsv` (9 rows), `tests/golden/math.tsv` (15,687 rows: every std f64 method ushma calls, over its call sites' ranges) and `src/constants.cyr` (the 9 public constants as `US_` consts). `--check` refuses a corpus from another glibc. `validate.py` independently re-derives the constants and every class A/B math row (5,869 checks)
- **table data** (M0 part 2b), generated by scripts/oracle:
  - `src/material_data.cyr`, `src/state_data.cyr`, `src/phase_data.cyr`, `src/chem_data.cyr`: the 38 public presets as typed column arrays with `US_<TABLE>_<PRESET>` index constants and name blobs, pinned by 42 golden rows;
  - `src/steam_data.cyr`: steam's private IAPWS-IF97 tables (42 saturated × 9, 72 superheated × 5, 9 tiers, the 10 Region 4 coefficients and scalars);
  - `src/constants_private.cyr`: the 8 other private constants.

  Private data comes from its Rust literals, with a proof that every `const` item in `rust-old/src` is generated or excused. `scripts/oracle/probe_tables.py` compiles all 1,052 values with the Cyrius compiler and checks them against the golden rows or an independent Python parse, with identical output on ELF, aarch64 and PE
- **docs/development/upstream/** (M0 step 10): the port plan's §8 in each target's issue format, filed into the targets on 2026-11-11 (uncommitted there) and kept here as the record. It holds 11 new issues, each with a repro under `repros/` that exits 1 while the defect is present (all run and reproduced on 6.7.18, from ushma and from the target):
  - hisab: the unused `[deps.sakshi]`;
  - ganita: `lu_solve`'s per-call allocation, and elimination's threshold scaling with b;
  - sandhi: `total_ms` not bounding a silent server;
  - sakshi: spans with no level;
  - bayan: raw control bytes and invalid UTF-8 in strings, and a serde-compatible printer;
  - cyrius: `mod` scope leaking into the includer, duplicate fns across bundles, the derive writer faulting on a 0 `Str`, and fdlopen / dynlib's PE warnings.

  It also holds 7 addenda to issues impetus, kimiya and pravash filed (the `-0` sign among them: pravash filed it minutes earlier, so ushma's duplicate was withdrawn), and a note for kimiya
- **scripts/probes/** (M0 part 3b): probes 1 (Simpson against Rust hisab), 7 (the daimon registration on sandhi, against a loopback stub), 9 (network solves against Rust hisab), 10 (`[lib]` / `[lib.ai]` bundles and consumers, `p10_manifest.sh`) and 11 (JSON against serde_json, `p11_json.sh`), with results in `docs/development/port-parity.md` §3; their Rust halves are oracle binaries (`scripts/oracle/src/bin/probe_hisab.rs`, `probe_ai.rs`, `probe_json.rs`). Each bundle links and runs a consumer on x86_64, PE and aarch64, and the probe-7 and probe-11 summaries are identical on all three
- **scripts/probes/** (M0 part 3a): probes 2 (memory), 3 (wide results), 4 (point lists), 5 (steam lookups), 6 (names and `mod`), 8 (span cost), 12 (f64 hygiene and the open field-store miscompile), with their results in `docs/development/port-parity.md` §3. They settle P6 as a written `us_` prefix (`mod`'s scope leaks into consumers; ADR 0001 amendment). They also set the design rules now in CLAUDE.md: level-gated spans (sakshi's spans have no level check, 192 ns each), out-pointer results, caller buffers for point lists, arenas for per-step scratch, and no truth test on an f64
- **`_us_cr_tanh`** (M0 part 2c-3, ADR 0004 amendment): ushma's own correctly rounded tanh, needed because `ganita_f64_tanh` is not correctly rounded (off on 52 of 860 golden rows). It computes tanh a = E/(E + 2) with E = expm1(2a), a double-double fast path with a proven rounding test and a triple-double accurate path, with `TANH_ONE` generated from MPFR. Every golden row is now gated (15,712 ok + 26 accepted)
- **scripts/oracle/cr_verify.py** (M0 part 2c-3): the dense MPFR sweep, vendored from kimiya (its families for exp, expm1, ln, pow and atan2, with the published hard cases) plus a tanh family. It runs about 6.6·10⁶ arguments per target, the accurate paths forced on 600,000, with 0 mismatches and identical bits and flags on ELF, PE and aarch64
- **src/mathx.cyr, src/mathx_data.cyr** (M0 part 2c-2b, ADR 0004, port plan P3): the vendored subset of kimiya's correctly rounded primitives (kimiya commit 5b14c85, mechanically extracted, identifiers renamed `_us_`). That's 55 functions: the Rust-semantics `min`, `max`, `powi`, `is_finite`, `is_nan`, and the correctly rounded `exp`, `expm1`, `ln`, `pow`, `atan2` with their fast and accurate entry points, all on + − × ÷, √ and integer operations, so they give the same bits on every target. `scripts/oracle/mathx_tables.py` regenerates the 46 tables from MPFR with kimiya's assertions, bit-identical to kimiya's (`--check-kimiya`). The parity harness gates them at 0 ulp against the corrected oracle: every row exact, including the 26 glibc misrounds, identically on ELF, `--bounds`, PE and aarch64
- **corrected oracle** (M0 part 2c-2a, ADR 0004): `scripts/oracle/mpfr_ref.py` gives correctly rounded binary64 results through the system MPFR. `validate.py --write-corrections` writes `tests/golden/math.cr.tsv`, the 26 class-C math rows where glibc 2.44 misrounds (exp_m1 23, ln 1, powf 2), and `validate.py` checks every listed call against MPFR and glibc. It shows `ganita_f64_tanh` is not correctly rounded (52 rows off where glibc is right), so the port gets its own correctly rounded tanh
- **tests/parity.tcyr** (M0 part 2c-1): the parity harness that replays the golden corpus against the port (`tests/parity_core.cyr` and `tests/ulp.cyr` from kimiya, renumbered for ADR 0004). Every constant and table preset matches the oracle; so do sqrt, abs, ceil and clamp on every row. The stdlib's transcendentals are measured against glibc (port-parity.md §2a), and the summary is identical on ELF, `--bounds`, PE and aarch64
- **scripts/probes/p0_math.cyr** — the math-counterpart probe: `f64_max`/`f64_min` return NaN for a NaN second operand where Rust returns the other value, so the port wraps them; `ganita_f64_pow(2, 0.5)` is 1 ulp off √2 (`docs/development/port-parity.md` §1)
- **scripts/gates.sh** — the full gate run (oracle, ELF / `--bounds` / PE / aarch64 tests, bench, fmt / lint / doc, vet, lock, audit)

### Fixed
- **transfer, state, entropy** — every port of Rust's f64 `sum` now starts at −0.0, as Rust's does (`_US_NEG_ZERO`, new in `src/validate.cyr`). Six sums had started at `0.0 - 0.0`, which is +0.0. So an empty series resistance, and Kay's rules over a −0.0 product, gave +0.0 where 1.4.0 gives −0.0. The other four (the parallel resistance, the mixture check and `entropy_of_mixing`'s two sums) differ only where no output shows it. Found at M4b
- **tests/parity_core.cyr** — a zero must now carry the oracle's sign in every class. +0 and −0 are one point on the ulp line, so the distance alone had passed the 5 rows above as 0 ulp, and the port's "class B bit-identical" was not quite true. The 8 decoder rows that read the integer `-0` stay accepted under ADR 5 (bayan's open `-0` issue), as before
- **tests/parity.tcyr** — the ADR 5 flag is reset for every row. It had been reset only for the material, state and phase routes, so a failing transfer, entropy, steam or cycle row after a JSON row would have been accepted as ADR 5. The tightened check showed it: it accepted the two series rows until this fix
- **tests/foundation.tcyr** — the predicates' "−0" cases tested +0 (`0.0 - 0.0`), and now use −0.0. `_us_neg` (new: Rust's unary minus as a sign flip, which `0.0 - x` is not at +0.0) is tested beside it

### Performance
- **steam** (M4b's review) — `superheated_by_entropy` probes the entropy alone during its bisection: the ln-P weight is computed once per inversion, and the volume's ln, ln and exp are skipped until the final lookup. Every bit is unchanged. Between tiers the inversion went from 4.6 to 2.2 µs, and `cycle/refrigeration` from 5.1 to 2.6 µs (Rust 890 ns). On a tier, `steam/superheated_by_entropy` went from 1,498 ns to 1,005 (Rust 342.05). All 6,037 inversion rows and every steam-cycle row are unchanged
- **steam** (M3b's review) — Region 4's coefficients are generated as scalar consts (immediates) instead of an indexed array, and the saturated interpolation reads its columns directly: `saturated_by_temp` 174 → 141 ns (Rust 18.19), `saturated_by_pressure` 189 → 159 (43.89), `saturation_pressure` 68 → 59 (6.45). Every golden row is unchanged; hisab's `calc_lerp` stays (an inlined lerp measured 3 ns faster)
- **validate** — `_us_finite_gt` / `_us_finite_ge` are one branch-free expression, `(x > bound) & (x < +∞)`, so `#inline` takes them; they were an `if` around a call. Three checks went from 21 ns to 6. Transfer, before → after, against Rust 1.4.0: diffusivity 33 → 16 ns (Rust 4.49), Biot 29 → 15 (3.99), parallel resistance 47 → 31 (9.28), conduction 31 → 23 (13.1), Reynolds 16 → 10 (2.57), fin heat 180 → 157 (89.76), coaxial disks 65 → 44 (8.19). Every golden row is unchanged (`docs/benchmarks-rust-v-cyrius.md`)
- **transfer** — Rust's `powi(4)` and `powi(2)` are written as the products `__powidf2` forms, (x·x)·(x·x) and x·x, instead of the loop. Radiation went from 40 to 25 ns (Rust 9.95), with the same bits on every row
- Rust 1.4.0 baseline: 63 benchmarks on the tag (`rust-old/bench-history.csv`, commit
  `2688099`), the Rust column for `docs/benchmarks-rust-v-cyrius.md`

## [1.4.0]

Toolchain/dependency update and a full correctness audit. Several fixes change
numerical results — see **Changed** for anything a consumer may observe.

### Changed
- **toolchain** — MSRV 1.89 → 1.99 (latest stable; edition 2024 is still the newest edition). CI MSRV job and release builds use 1.99
- **deps** — hisab 1.3 → 1.4, reqwest 0.12 → 0.13 (`ai` feature; default TLS is now rustls/aws-lc-rs instead of native-tls), criterion 0.5 → 0.8 (dev; `criterion::black_box` → `std::hint::black_box`), serde 1.0.229, thiserror 2.0.21, tracing 0.1.44, tracing-subscriber 0.3.23, serde_json 1.0.154, tokio 1.54.1
- **deps** — `ai` feature no longer pulls in `tokio` with `features = ["full"]`; nothing in the crate used it (reqwest still brings its own tokio dependency)
- **ci** — actions/checkout v4 → v7, codecov-action v4 → v7, action-gh-release v2 → v3 (Node 24 runtimes)
- **deny** — allow `ISC` and `CDLA-Permissive-2.0` (rustls / webpki-root-certs licenses, `ai` feature only)
- **validation (all modules)** — every positivity/range guard now also rejects NaN and ±∞. `x <= 0.0` let NaN through (all NaN comparisons are false), so e.g. `carnot_efficiency(NaN, 300.0)`, `otto_cycle(NaN, …)`, `phase_at(NaN, …)` and `saturated_by_temperature(NaN)` returned `Ok(NaN…)`; they now return errors. Implemented with crate-private `validate::finite_gt` / `finite_ge`
- **steam** — saturated and superheated tables regenerated from IAPWS-IF97 on the same grid. Fixes s_g 4–26 J/(kg·K) low at 523–633 K, h_g ≠ h_f + h_fg by ±1 kJ/kg in 16 rows (made `quality_from_enthalpy(h_g)` fail with x = 1.0005), the critical-point v_c (0.003155 → 1/322 m³/kg), and superheated rows off by up to 7 kJ/kg / 16 J/(kg·K)
- **steam** — saturation pressure/temperature now come from the exact IF97 Region 4 equations (linear-in-T interpolation was up to 3.3 % off between rows); v_g is interpolated through the smooth ratio P·v_g/T instead of linearly (leave-one-out test on the table: mean error 4.3 % → 0.2 %, worst 15 % → 3.5 %)
- **steam** — `superheated_lookup` interpolates between pressure tiers in ln P (log-log for v). Linear-in-P interpolation gave v 80 % high at 30 kPa and 12 % high at 150 kPa
- **steam** — `superheated_lookup` uses the saturated-vapor state as the lowest point of each tier, so states between T_sat and the first tabulated row are covered
- **cycle** — refrigeration `net_work` is now negative (−W_compressor), matching the documented "positive = work produced" and the first-law identity `heat_in − heat_out = net_work` that every other cycle satisfies. **Behavior change** for consumers reading `net_work` of a refrigeration result
- **cycle** — `compare_cycles` uses the Carnot COP T_L/(T_H − T_L) for refrigeration results (it divided a COP by the Carnot *efficiency*: COP 5 at 300/250 K gave η_II = 30 instead of 1.0), and reports 0 for non-physical limits (T_L ≤ 0 used to give η_Carnot = 1.33)
- **cycle** — `cycle_ts_diagram` traces ideal-gas isochoric/isobaric paths as T = T_a·(T_b/T_a)^f (was a straight line: Otto 2→3 midpoint 1892 K vs 1460 K); `cycle_pv_diagram` draws isentropes through both end points with k = ln(P_b/P_a)/ln(V_a/V_b), fixing flat, discontinuous steam-cycle isentropes (γ = 0), and draws isotherms as hyperbolas
- **phase** — `phase_at` vapor-pressure curve now runs through the triple point, normal boiling point and critical point (piecewise ln P vs 1/T). The single-reference curve overshot P_c for water (26 MPa at T_c) and was 85 % off for CO₂, so e.g. CO₂ at (300 K, 3 MPa) was "Liquid" (now Gas) and water at (640 K, 23 MPa) was "Gas" (now Liquid)
- **phase** — below the triple point `phase_at` uses a sublimation curve instead of always returning Solid (CO₂ at 200 K / 1 atm is Gas; water at 250 K / 10 Pa is Gas)
- **phase** — `heat_for_phase_change` treats substances that sublime at 1 atm (CO₂) as a single sublimation event instead of adding fusion + vaporization at unrelated temperatures
- **state** — `GAS_CONSTANT` = 8.314 462 618 153 24 (exact N_A·k_B; was truncated at 1e-11 relative)
- **state** — `compressibility_pitzer` documented as the Pitzer–Abbott virial correlation (not Lee-Kesler) and returns an error for non-physical Z ≤ 0 (e.g. Pr = 5, Tr = 1 gave Z = −0.695)
- **stat** — `debye_cv` integrand rewritten as x⁴e⁻ˣ/expm1(−x)² and switched to the exact T³ law for Θ_D/T > 50; it returned NaN for Θ_D/T > 709 (e.g. diamond at 3 K) and lost accuracy from Θ_D/T ≈ 300
- **chem** — `equilibrium_constant` and `vant_hoff_k` return an error instead of `Ok(∞)` / `Ok(0)` when K leaves the f64 range (propane combustion: ln K ≈ 837)
- **transfer** — `view_factor_parallel_plates` uses the asymptote (XY/π)(1 − (X²+Y²)/3) for X, Y < 0.01, where the exact formula cancels to 0
- **numerical** — `ThermalGrid1D::step_explicit` stability limit is Fo·(1 + Bi) ≤ 0.5 with a convective boundary (accepting Fo = 0.4 at Bi = 10 drove T to −400 K in one step)
- **numerical** — `ThermalGrid2D` applies `Convective` boundaries (T_b = (T_in + Bi·T∞)/(1 + Bi)); they were silently held at the initial temperature
- **numerical** — `ThermalNetwork::solve` errors when a free node has no path to a fixed node (a floating triangle "solved" to 0 K) and scales rows by their diagonal, so very small conductances (R ≈ 1e14 K/W) no longer trip hisab's absolute pivot threshold
- **bridge** — `basal_metabolic_rate` uses Kleiber's 70 kcal/day = 3.39 W (was 3.5 W, 3 % high: 70 kg → 82.0 W, was 84.7 W)
- **bridge** — `humidity_to_wet_bulb` never exceeds the dry-bulb temperature; docs state the Stull fit range. `velocity_to_convection_coefficient` documented as Dittus-Boelter internal flow (not flat plate) and returns 0 for invalid k/ρ/c_p (k = 0 gave NaN). `equilibrium_temperature` uses the exact R
- **material** — `AIR` at 300 K per Incropera Table A.4: ρ 1.225 → 1.1614 kg/m³ (1.225 is the 288 K ISA value), c_p 1005 → 1007, k 0.026 → 0.0263
- **integration/soorat** — empty or all-NaN fields report `min_temp`/`max_temp` = NaN instead of ±f64::MAX; missing network temperatures are NaN instead of 0 K

### Added
- **steam** — `saturation_pressure(T)` and `saturation_temperature(P)` (IAPWS-IF97 Region 4, verified against the IF97 tables to 1e-8)
- **steam** — `superheated_by_entropy(P, s)`: superheated state with a given entropy (turbine/compressor exit). Exact inversion of the piecewise-linear table model (~7 lookups, no iteration tolerance)
- **phase** — `SubstanceData::vapor_pressure`, `sublimation_pressure`, `latent_heat_sublimation`, `sublimes_at_atmospheric`
- **chem** — `ln_equilibrium_constant` for reactions whose K overflows f64
- **numerical** — `ThermalGrid1D::max_stable_dt()` and `MAX_ADAPTIVE_SUBSTEPS`
- **docs** — `[package.metadata.docs.rs] all-features = true`, so docs.rs documents `steam`, `ai`, `logging` and `soorat-compat`; `soorat-compat` added to the feature tables
- **tests** — 80 new regression tests (399 → 479), with reference values from IAPWS-IF97, Incropera and NIST
- **benches** — `transfer/view_factor_perpendicular`, `steam/saturation_pressure`, `steam/superheated_between_tiers`, `steam/superheated_by_entropy`, `stat/debye_cv_low_t`, `numerical/adaptive_1d_4substeps`

### Fixed
- **steam** — `superheated_lookup` rejected exact table rows at tier pressures, e.g. (423.15 K, 200 kPa) and (523.15 K, 2 MPa), because it always interpolated toward the next tier up, whose range starts higher. It now uses the single tier. Errors also report the requested pressure, not the tier's
- **steam** — `validate_quality` tolerates round-off (x = 1 + 2e-16 was rejected) and clamps to [0, 1]
- **cycle** — `rankine_cycle` treated a superheated turbine exhaust as saturated vapor. h4 came out too low and η too high: 14.7 % instead of 13.2 % at 100 kPa / 500 kPa / 573 K. It now solves s(T, P_cond) = s4
- **cycle** — `refrigeration_cycle` silently replaced a failed compressor-exit bisection with the bracket midpoint. With a bracket that was too narrow, this gave COP 6.49, above the Carnot 4.0. It now solves s(T, P_cond) = s1 over the whole table range and propagates errors
- **cycle** — `cycle_pv_diagram` panicked when `processes` was shorter than `state_points` (both fields are public and deserializable)
- **transfer** — `view_factor_perpendicular_plates` used a garbled formula. Unit squares gave 0.298 instead of 0.2000, and other inputs were clamped to 0 or 1. It now uses the full Incropera/Howell C-14 expression, and reciprocity holds exactly
- **transfer** — `view_factor_coaxial_disks` is rewritten without cancellation. It was 25 % low at r = 1e-4·d, and 0 at 3e-5·d
- **transfer** — negative k, ρ, c_p or R are now rejected consistently:
  - `thermal_resistance_conduction` accepted k < 0 and A < 0 when their product was positive
  - `thermal_resistance_parallel` accepted R < 0; [2, −1] gave −2 K/W
  - also covered: `thermal_diffusivity`, `biot_number`, `nusselt_number`, `fin_parameter`, `fin_efficiency_rectangular` and `prandtl_number`
- **numerical** — `step_adaptive` rejected about 5 % of valid steps through round-off (Fo = 0.5000000000000001). It now caps the sub-step count, and it validates dt: NaN dt returned Ok without stepping
- **numerical** — `step_explicit` and `step_crank_nicolson` accepted negative or NaN dt, which ran diffusion backwards
- **numerical** — `ThermalGrid2D::solve_steady_state` reported convergence on an all-NaN field, because `f64::max` drops NaN
- **numerical** — `ThermalNetwork::new(0).add_resistance(0, 0, …)` underflowed `num_nodes - 1`
- **numerical** — grids or networks with inconsistent public fields (deserialized, or edited by hand) now return errors instead of panicking on out-of-bounds indexing
- **integration/soorat** — `HeatFluxField::from_grid_2d` and `ThermalGridVisualization::from_grid_2d` no longer panic or report dimensions that don't match `values` when the grid's rows disagree with nx/ny
- **phase** — `clausius_clapeyron_pressure` returns an error instead of `Ok(∞)` on overflow
- **state** — `ideal_gas_*`, `van_der_waals_pressure` and `isothermal_work` reject negative moles: `ideal_gas_temperature(ATM, 0.0224, −1)` gave −273 K
- **state** — `adiabatic_temperature` validates T₁
- **state** — `redlich_kwong_pressure` and `peng_robinson_pressure` validate Tc, Pc, ω and V_m; a negative Tc gave `Ok(NaN)`
- **state** — Kay's-rule mixture validation rejects NaN and negative mole fractions
- **entropy** — `ideal_gas_entropy_change`, `isothermal_entropy_change` and `entropy_of_mixing` validate moles; `ideal_gas_entropy_change` also validates Cv
- **chem** — `adiabatic_flame_temperature` rejects negative or NaN stoichiometric coefficients
- **stat** — `partition_vibrational` uses `exp_m1`, so it stays accurate for hν ≪ kT
- **bridge** — `muscle_work_to_heat` with η = 1 gave 300 % heat; it now gives 0
- **bridge** — NaN inputs fall back to the documented values instead of propagating

### Security
- **ai** — `DaimonConfig`'s `Debug` output redacts `api_key`; a derived `Debug` would have leaked it through `{:?}` or tracing fields
- **ai** — `register_agent` checks the HTTP status with `error_for_status`, so a 401 or 500 response is reported as such rather than as "missing agent_id"
- **ai** — a trailing `/` on the endpoint no longer produces `//v1/…`

### Performance
Measured as an interleaved A/B run of 1.3.0 against 1.4.0, both under the criterion 0.8 harness. Each figure is the minimum of the per-round medians over 2 rounds. Note that the criterion 0.5 → 0.8 switch by itself shifts some `#[tracing::instrument]`-ed functions by several ns, so `bench-history.csv` runs from before and after this release are not directly comparable.

- **`InvalidParameter` errors** are built out of line through a `#[cold]` helper (`validate::invalid_param`), and small validated functions are `#[inline]`. Before this, adding an error branch that formats a message turned leaf functions into non-inlinable calls: `ideal_gas_volume` went from 1.6 ns to 9.4 ns in an intermediate build.
- **`ThermalGrid1D`:** `step_explicit` updates in place with no clone, and Crank-Nicolson uses 2 scratch vectors instead of 6.
- **steam**
  - `wet_steam_properties`: 9.38 → 1.74 ns (−81 %)
  - `quality_from_enthalpy`: 8.23 → 2.08 ns (−75 %)
- **stat**
  - `rms_speed`: 14.78 → 3.88 ns (−74 %)
  - `einstein_cv`: 12.35 → 8.81 ns (−29 %)
- **transfer**
  - `reynolds`: 8.96 → 2.65 ns (−70 %)
  - `biot_number`: 8.59 → 4.09 ns (−52 %)
  - `diffusivity`: 7.79 → 4.65 ns (−40 %)
  - `eff_counter`: 11.90 → 8.33 ns (−30 %)
- **phase**
  - `clausius_clapeyron_slope`: 8.32 → 3.27 ns (−61 %)
  - `clausius_clapeyron_pressure`: 19.74 → 10.72 ns (−46 %)
- **entropy** — `carnot_cop`: 8.01 → 3.40 ns (−58 %)
- **state** — `compressibility_factor`: 8.25 → 4.44 ns (−46 %)
- **chem** — `reaction_enthalpy`: 4.56 → 2.48 ns (−46 %)
- **numerical**
  - `explicit_1d`: 67.70 → 43.93 ns (−35 %)
  - `crank_nicolson_1d`: 1059.7 → 838.2 ns (−21 %)
  - `gauss_seidel_2d_20x20`: 1.44 → 1.34 ms (−7 %)
- **Costlier, for correctness:**
  - `steam/saturated_by_pressure`: 12.2 → 35.7 ns (exact IF97 backward equation)
  - `steam/saturated_by_temp`: 14.7 → 17.9 ns
  - `cycle/rankine`: 67 → 114 ns
  - `cycle/refrigeration`: 422 → 708 ns (exact compressor-exit solve)
  - `phase/phase_lookup`: 14 → 22 ns (three-anchor vapor-pressure curve)
  - `chem/adiabatic_flame`: 13 → 27 ns (coefficient validation)
  - `stat/debye_cv`: 1.52 → 2.04 µs (overflow-proof integrand)
  - `state/ideal_gas_volume`: 1.6 → 3.7 ns (NaN-safe guards)
- **New benchmarks**
  - `transfer/view_factor_perpendicular`: 50 ns
  - `steam/saturation_pressure`: 6.0 ns
  - `steam/superheated_between_tiers`: 60 ns
  - `steam/superheated_by_entropy`: 0.33 µs
  - `stat/debye_cv_low_t`: 8.9 ns
  - `numerical/adaptive_1d_4substeps`: 117 ns

## [1.3.0]

### Added
- **bridge** — sharira: muscle work to heat generation, basal metabolic rate (Kleiber's law), body surface area (Du Bois formula), skin radiation loss (Stefan-Boltzmann), muscle activation heat, evaporative cooling from sweat rate

## [1.2.0]

### Added
- **bridge** — cross-crate primitive-value bridges for bijli (Joule heating rate, EM absorption to heat source), kimiya (reaction heat rate, equilibrium temperature), dravya (thermal strain, thermal gradient to stress), badal (altitude to temperature, humidity to wet-bulb), pravash (velocity to convection coefficient, TKE to eddy diffusivity)
- **integration/soorat** — feature-gated `soorat-compat` module with visualization data structures: `ThermalGridVisualization` (2D heatmap from `ThermalGrid2D`), `TemperatureProfile` (1D line from `ThermalGrid1D`), `CycleDiagramData` (T-s and P-v plots from `CycleResult`), `ThermalNetworkVisualization` (node-link graph), `HeatFluxField` (gradient vectors from `ThermalGrid2D`)

### Updated
- hisab 1.0.0 -> 1.3.0, zerocopy 0.8.47 -> 0.8.48, iri-string 0.7.11 -> 0.7.12

## [1.1.0] - 2026-03-26

### Changed
- **state** — `van_der_waals_pressure()` now rejects T=0 K (`<= 0.0` instead of `< 0.0`) for consistency with real gas EOS functions
- **state** — `ideal_gas_temperature()` now validates pressure and volume parameters
- **entropy** — `entropy_of_mixing()` now accepts zero mole fractions (absent components contribute nothing: lim x→0+ of x·ln(x) = 0)

### Fixed
- **numerical** — removed `unwrap()` in `ThermalGrid1D::new()` boundary condition application; replaced with safe index access
- **ai** — `DaimonClient::register_agent()` now returns proper error instead of silently defaulting when `agent_id` field is missing from API response
- **phase** — `heat_for_phase_change()` now validates `specific_heat > 0` before computation

### Removed
- **deps** — removed unused `impetus` dependency (was declared but never imported)

### Added
- **tracing** — added `#[tracing::instrument]` to 23 core functions across transfer, state, entropy, phase, cycle, numerical, and chem modules for structured debug logging

## [1.0.0] - 2026-03-25

### Added

#### transfer — Heat Transfer
- Fourier's law: conduction through solid materials
- Newton's law of cooling: convective heat transfer
- Stefan-Boltzmann law: radiative heat transfer
- Thermal resistance: conduction, convection, series, parallel networks
- Heat stored (Q = mcΔT), thermal diffusivity, Biot number
- Lumped capacitance transient cooling model
- Reynolds, Prandtl, Nusselt number calculations
- Dittus-Boelter (turbulent pipe) and Churchill-Chu (natural convection) correlations
- Laminar and thermal boundary layer thickness (Blasius)
- Fin heat transfer: rectangular fins, efficiency, effectiveness
- Heat exchangers: LMTD method (parallel/counter flow)
- Heat exchangers: ε-NTU method with special cases (Cr=0, Cr=1)
- Radiation view factors: parallel plates, perpendicular plates, coaxial disks

#### state — Equations of State
- Ideal gas law: pressure, volume, temperature (PV = nRT)
- Van der Waals equation for real gases
- Redlich-Kwong equation of state
- Peng-Robinson equation of state (with acentric factor)
- Virial equation of state (2nd and 3rd order)
- Pitzer correlation for second virial coefficient
- Generalized compressibility (Pitzer Z⁰ + ωZ¹)
- Isothermal, isobaric, and adiabatic work
- Compressibility factor
- GasData struct with 8 built-in gases (N₂, O₂, CO₂, CH₄, H₂O, NH₃, C₂H₆, C₃H₈)
- Kay's rule mixture rules for Tc, Pc, ω

#### entropy — Thermodynamic Potentials
- Ideal gas entropy change, isothermal entropy change
- Heat transfer entropy (ΔS = Q/T)
- Carnot efficiency and COP for refrigeration
- Helmholtz and Gibbs free energy
- Clausius inequality check (spontaneity)
- Entropy of mixing for ideal gases

#### material — Thermal Properties
- ThermalMaterial struct with Cow<'static, str> name
- 10 built-in materials: copper, aluminum, iron, steel, glass, water, air, wood (oak), concrete, diamond
- Thermal diffusivity and volumetric heat capacity methods
- Serde serialization/deserialization

#### phase — Phase Transitions
- SubstanceData struct with triple/critical points, latent heats, molar mass
- 8 built-in substances: water, nitrogen, oxygen, CO₂, ammonia, methane, ethanol, mercury
- Clausius-Clapeyron equation (exact slope and integrated pressure form)
- Phase diagram lookup (solid/liquid/gas/supercritical)
- Latent heat helpers: fusion, vaporization, multi-phase heat budget

#### steam — Steam Tables
- 41-entry saturated steam table (273.16 K to 647.096 K, IAPWS-IF97)
- Saturated lookup by temperature and pressure (binary search + hisab lerp)
- 72-entry superheated steam table (9 pressure tiers, bilinear interpolation)
- Quality calculations from volume, enthalpy, entropy
- Wet steam mixture properties

#### cycle — Thermodynamic Cycles
- Otto cycle (spark ignition engine)
- Diesel cycle (compression ignition)
- Brayton cycle (gas turbine) with back-work ratio
- Rankine cycle (steam power plant, saturated and superheated)
- Vapor-compression refrigeration cycle (with hisab bisection for isentropic state)
- Heat pump COP
- T-s and P-v diagram data generation
- Cycle comparison with second-law efficiency

#### numerical — Numerical Methods
- 1D transient conduction: explicit (forward Euler) finite difference
- 1D transient conduction: implicit (Crank-Nicolson) via Thomas algorithm
- Adaptive time stepping for explicit stability
- 2D steady-state conduction (Gauss-Seidel iteration)
- Thermal resistance network solver (via hisab gaussian elimination)
- Boundary conditions: fixed, insulated, convective

#### chem — Chemical Thermodynamics
- Species struct with formation enthalpy, Gibbs energy, heat capacity
- 12 built-in species: H₂, O₂, N₂, H₂O(g), H₂O(l), CO₂, CO, CH₄, C₂H₆, C₃H₈, NH₃, NO
- Hess's law: reaction enthalpy from formation enthalpies
- Gibbs energy of reaction
- Equilibrium constant from ΔG
- Van't Hoff equation (K vs temperature)
- Adiabatic flame temperature (constant-Cp model)

#### stat — Statistical Thermodynamics
- Boltzmann distribution: probability, population ratio
- Entropy from microstates (S = k·ln(W))
- Equipartition theorem: energy per molecule, molar Cv
- Partition functions: translational, rotational, vibrational
- Maxwell-Boltzmann speed distribution: PDF, most probable, mean, RMS speeds
- Einstein model for solid heat capacity
- Debye model for solid heat capacity (via hisab Simpson integration)
- Physical constants: Planck h, Avogadro Nₐ

#### Infrastructure
- Error type: UshmaError with #[non_exhaustive], 16 variants
- 14 feature flags (9 default + steam, ai, logging, full)
- 362 tests (332 unit + 30 integration)
- 57 Criterion benchmarks with CSV history tracking
- 5 examples: basic_thermo, heat_engine, building_thermal, cooling_system, phase_change
- hisab integration: Simpson integration, lerp, bisection, gaussian elimination
- AI integration: daimon/hoosh client (feature-gated)
- Structured logging via USHMA_LOG (feature-gated)
- CI/CD: GitHub Actions (format, clippy, audit, deny, test, coverage, semver, doc)
- Supply chain: cargo-deny, cargo-audit
- Scripts: bench-history.sh, version-bump.sh
