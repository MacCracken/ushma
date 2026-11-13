# Ushma Roadmap

> **Ushma** is the thermodynamics simulation library for AGNOS, in Cyrius. Math foundations come from the Cyrius
> stdlib (ganita) and [hisab](https://github.com/MacCracken/hisab).

This page is what is still to do. What is done is in the [CHANGELOG](../../CHANGELOG.md): 1.0–1.4 (Rust) and 2.0.0,
the Cyrius port. The port's own record is the [port plan](ushma-port-plan.md).

## Scope

Ushma owns the **physics of heat**: how thermal energy transfers, transforms, and equilibrates. It provides the math;
consumers decide what to do with it (simulate engines, model climate, thermal effects in games).

Ushma does NOT own:
- **Physics engine** → impetus (rigid body dynamics, collision)
- **Math primitives** → hisab, ganita (vectors, geometry, calculus, linear algebra, ODE/PDE solvers)
- **Rendering** → soorat/kiran (they consume ushma for thermal visualization)
- **Fluid dynamics** → pravash (CFD, Navier-Stokes)

## How 2.x work lands

- **Rust 1.4.0 stays the oracle.** It is frozen in `rust-old/` and never changes. A change that moves a result away
  from 1.4.0's gets an [ADR 0002](../adr/0002-port-api-divergences.md) row, and the golden rows it moves are mapped
  under it, so the corpus keeps comparing every function against 1.4.0.
- **New physics has no oracle.** It comes with reference tests against the literature and an independent model in
  `scripts/oracle/validate.py`, as the port's own refusals did.
- **The API only grows.** Enums and error codes are append-only, and a new form joins beside the old one rather than
  replacing it. Nothing below needs a 3.0.
- **Versions:** new capability is a minor release (2.x.0); a fix, or catch-up work that changes no result, a patch
  (2.x.y).
- **Demand-gated:** each item waits for a consumer that needs it (kiran, joshua, soorat, or a bridge's crate).

## 2.x.x — planned work

### numerical
- [ ] **2D walls, second order.** A ghost-node treatment for insulated and convective walls. They are first-order at
  the wall now; the 1D grid is second-order already. Moves results: needs an ADR 0002 row.
- [ ] **Gauss–Seidel:** a residual-based stopping criterion, and SOR acceleration, as options beside today's
  largest-change test.

### transfer and bridge
- [ ] **Dittus–Boelter for cooling** (n = 0.3), beside today's heating form (n = 0.4).
- [ ] **`us_velocity_to_convection_coefficient`:** laminar and flat-plate branches. It is a turbulent pipe-flow
  correlation only today.
- [ ] **`us_altitude_to_temperature` above the tropopause:** the isothermal layer past about 11 km. The troposphere's
  lapse rate holds below it.

### steam
- [ ] **A denser superheated grid.** The gap between the 5 and 10 MPa tiers limits h to about 0.2 %.
- [ ] **The compressed-liquid region.**
- [ ] **The saturation dome** in the T-s diagrams of steam cycles. Their processes are straight lines today.

### phase
- [ ] **A sloped melting curve**, for water's negative dT/dP.
- [ ] **Two data points to review:** ETHANOL's triple-point pressure and MERCURY's critical pressure.

### ai
- [ ] **A dedicated error code for network failures.** Today every registration failure is `US_ERR_INVALID_PARAMETER`,
  as 1.4.0's is `InvalidParameter`. A new code, appended, moves the outcome of the golden registrations it covers (a
  refused connection, a timeout, perhaps an HTTP error status): needs an ADR 0002 row.

## Upstream blockers — catch-up when fixes land

Each item is a workaround ushma carries for a defect or a gap upstream. When the fix ships, the catch-up removes the
workaround. Every result stays bit-identical: the corpus, the gates on all four targets and the benchmark record show
it. Each is patch-level unless its row says otherwise. The issues and their repros are in
[`upstream/`](upstream/README.md). An item marked *not filed* waits for the maintainer's word.

### Math, JSON and memory

| Upstream | Status | What ushma does meanwhile | Catch-up when it lands |
|---|---|---|---|
| ganita: correctly rounded exp, expm1, ln, pow, atan2, tanh (`2026-11-09-correctly-rounded-exp-sin-cos-hypot`, ushma's addendum) | filed | vendors kimiya's `src/mathx.cyr` and its tables, plus its own `_us_cr_tanh` (ADR 0004) | call ganita's; delete `mathx.cyr`, `mathx_data.cyr` and `mathx_tables.py`; run `cr_verify.py` against ganita's; every class-C row stays at 0 ulp |
| bayan: a serde_json-compatible shortest f64 printer (`2026-11-11-no-serde-compatible-shortest-f64-printer`) | filed | vendors impetus's `src/fmtf.cyr`, an exact bignum printer: about 3 µs for a full 17-digit double (a 20 × 20 heatmap's document, 1.1–1.3 ms) | write f64 through bayan's; delete `fmtf.cyr` and `vendor_fmtf.py`; `fmtf.tsv`'s 26,088 rows and every JSON row unchanged; re-measure the JSON benches |
| bayan: refuse raw control bytes and invalid UTF-8 in strings (`2026-11-11-json-strings-accept-raw-control-bytes-and-invalid-utf8`) | filed | `json.cyr` pre-scans each document (`_us_jr_raw_ok`) and validates each string (`_us_utf8_valid`) | delete both; the JSON rows unchanged |
| bayan: the integer `-0` keeps its sign (pravash's `2026-11-11-json-integer-negative-zero-loses-its-sign`, ushma's addendum) | filed | the decoder rows that read an integer `-0` are accepted under ADR 5 | those rows become ok; drop the harness's `-0` acceptance; amend ADR 0005 |
| bayan: an opt-in for duplicate keys and deep nesting in ignored members (`2026-11-09-json-duplicate-key-has-no-opt-in-for-ignored-members`, ushma's addendum) | filed | such documents are refused, ADR 5's mapped rows | opt in; those rows become ok; amend ADR 0005 |
| bayan: integers above i64::MAX (`2026-11-09-json-integer-above-i64-max-refuses-the-document`, ushma's addendum) | filed | refused, ADR 5's mapped rows (a usize or an f64 field ≥ 2⁶³) | read them into f64 fields. A usize ≥ 2⁶³ needs a decision first: ushma holds counts as i64 words. Then those rows become ok, under an ADR 0005 amendment |
| ganita: an allocator-taking `mat_lu_solve` (`2026-11-11-mat-lu-solve-allocates-per-call`) | filed | `numerical.cyr` runs the forward and back substitution itself over ganita's factors (ADR 0003, second amendment) | call ganita's with the network's scratch; the 600-network bit-identity test, all 476 solve rows and the 0-byte soak unchanged; re-measure `thermal_network_5` |
| cyrius stdlib (math): Rust-semantics `f64_min` / `f64_max` (`2026-11-11-f64-min-max-drop-nan-only-as-first-operand`, ushma's addendum); `f64_is_nan`, `f64_is_finite`, `F64_INF`, `F64_NAN`, `F64_EPSILON` (*not filed*: add to impetus's filing first) | min/max filed; the rest not filed | `_us_min`, `_us_max`, `_us_is_nan`, `_us_is_finite`, `_us_nan`, `_US_INF`, `_US_EPSILON` | use the stdlib's; delete the private copies |
| cyrius stdlib (string): `memcpy` and `memset` copy a byte at a time (ranga's `2026-11-02` issue) | open; ushma's numbers not yet added | copies words itself (soorat's `_us_sr_copy`; the cycle results store in place instead, M4a) | use `memcpy`; re-measure `soorat/*` |
| cyrius stdlib (fmt): `fmt_float_buf` at 0 decimals writes `123.0`, and rounds in f64 rather than exactly | *not filed* | `docs/examples/print.cyr` cuts the text at the point | drop the cut; `scripts/examples.sh` stays byte-identical |

### Tooling, logging and the network

| Upstream | Status | What ushma does meanwhile | Catch-up when it lands |
|---|---|---|---|
| hisab: drop its `[deps.sakshi]` (`2026-11-11-deps-sakshi-package-warns-every-consumer`) | filed | every resolve warns `refusing to overwrite stdlib leaf 'sakshi'`, filtered by the gates | bump the hisab pin and `cyrius deps --relock`; drop the warning from `gates.sh`'s known list and `state.md` |
| sakshi: a level-taking span entry (`2026-11-11-span-enter-has-no-level-check`) | filed | `_us_span_enter` / `_us_span_exit` test the level first (probe 8: 192 ns into the ring otherwise) | call sakshi's; re-measure the span cost |
| sandhi: `total_ms` bounds a read blocked before the first byte (`2026-11-11-total-ms-does-not-bound-a-blocked-first-read`) | filed | the daimon client sets `read_ms` as well as `total_ms` (ADR 0002 row 7) | `read_ms` may stay (it is harmless) or go: amend row 7, and re-run `tests/ai.tcyr`'s silent-server case with `total_ms` alone |
| cyrius: `#inline` takes a small body with an `if` or a local (pravash's `2026-11-11-inline-refuses-small-numeric-helpers`, ushma's addendum) | filed | the finite checks are written as one expression; the level-gated spans can't be, and cost 7 ns a call | `#inline` the span helpers and the predicates; re-measure |
| cyrius: the `field-store-of-struct-call-copies-a-field-argument` miscompile (`2026-11-09`) | open (not ushma's filing) | CLAUDE.md's workaround: bind a struct-returning call to a local before storing it into a field | relax the rule in CLAUDE.md |
| cyrius stdlib (fdlopen, dynlib): six raw syscalls with no PE route (`2026-11-11-fdlopen-dynlib-raw-syscalls-warn-on-every-pe-build`) | filed | every `[lib.ai]` PE build keeps one warning, filtered as `^warning:lib/` | drop the note from `state.md` |
| sigil: its TLS buffers' `large static data` warning (prakash's record) | open (not ushma's filing) | filtered by the gates | drop it from `gates.sh`'s known list and `state.md` |

### Filed, with nothing to catch up

These are fixed in ushma's design, so a fix upstream changes nothing here:
- **cyrius `mod` scope** (`2026-11-11-mod-scope-runs-past-the-included-file`) and **duplicate functions across
  bundles** (`2026-11-11-duplicate-fn-across-bundles-last-wins`). ushma's names carry the written `us_` prefix (ADR
  0001), and that prefix is the 2.x API.
- **cyrius `#derive(Serialize)`** (`2026-11-11-derive-writer-faults-on-a-zero-str-field` and the derive-json shapes
  addendum). ushma's JSON codecs are hand-written (ADR 0005).
- **ganita `mat_gaussian_elim`'s pivot scale** (`2026-11-11-gaussian-elim-pivot-scale-includes-rhs`). ushma solves
  on LU (ADR 0003).
- **`cyrius port` clobbering project docs** (`2026-11-08`). It was a one-time scaffold.

Notes for the sibling ports (kimiya, pravash, impetus, tanmatra) that are not filed while those ports are in progress
are in the port plan's §8. They need nothing from ushma.

## Consumers

| Consumer | What it uses |
|----------|-------------|
| **kiran** | Thermal simulation in game worlds (heat zones, fire, cooling) |
| **joshua** | Simulation: heat engines, thermal experiments, climate models |
| **soorat** | Thermal visualization (the soorat data types) |
| **badal** | Weather/atmospheric modeling (the badal bridges) |
| **dravya** | Material/substance modeling (the dravya bridges) |

None consumes the Cyrius library yet; they use the Rust crate, 1.4.0. The current list is in [state.md](state.md).

## Boundary with Other Crates

| Feature | ushma | other |
|---------|-------|-------|
| Heat transfer math | Yes | — |
| Rigid body dynamics | — | impetus |
| ODE/PDE solvers | — | hisab |
| Fluid dynamics (CFD) | — | pravash |
| Vector/matrix math | — | hisab, ganita |
| Thermal visualization | — | kiran/soorat |
| Material science (stress/strain) | — | dravya |
