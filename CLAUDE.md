# Ushma — Claude Code Instructions

> **Core rule**: this file is **preferences, process, and procedures** —
> durable rules that change rarely. Volatile state (current version,
> module line counts, port progress, test counts, consumers) lives in
> [`docs/development/state.md`](docs/development/state.md).
> Do not inline state here.

## Project Identity

**Ushma** (Sanskrit: heat/thermal energy) — Thermodynamics simulation — heat transfer, entropy, equations of state, thermal material properties

- **Type**: Flat library (Cyrius), ported from the Rust crate 1.4.0; the port ships as 2.0.0
- **License**: GPL-3.0-only
- **Language**: Cyrius. The toolchain is pinned in `cyrius.cyml` `[package].cyrius`
- **Version**: `VERSION` at the project root is the source of truth. Don't inline the number here
- **Rust oracle**: ushma 1.4.0 (tag `1.4.0`), frozen at `rust-old/`. It is never edited and never
  re-released. It is the parity reference
- **Plan of record**: [`docs/development/ushma-port-plan.md`](docs/development/ushma-port-plan.md),
  with ADRs 0001–0005 in [`docs/adr/`](docs/adr/)
- **Standards**: [First-Party Standards](https://github.com/MacCracken/agnosticos/blob/main/docs/development/applications/first-party-standards.md) · [First-Party Documentation](https://github.com/MacCracken/agnosticos/blob/main/docs/development/applications/first-party-documentation.md)

## Goal

ushma owns the physics of heat for AGNOS: how thermal energy transfers, transforms and
equilibrates. That covers conduction, convection and radiation, equations of state, entropy,
phase change and steam, heat engines and cycles, numerical conduction, chemical and statistical
thermodynamics, the cross-crate bridges, the soorat visualization data, and the daimon/hoosh
client. It provides the math; consumers decide what to do with it.

## Consumers

kiran (thermal simulation), joshua (heat engines, thermal experiments) and soorat (thermal
visualization, its `thermo` feature) use the Rust crate. None consumes the Cyrius port yet. The
current list is in `state.md`.

## Quick Start

```sh
cyrius deps                                   # resolve dependencies (vendors lib/, incl. hisab)
cyrius build src/main.cyr build/ushma         # compile
cyrius test                                   # run tests/*.tcyr (tests/parity.tcyr replays the golden corpus)
cyrius bench                                  # run benchmarks
cyrius distlib --all                          # dist/ushma.cyr and dist/ushma-ai.cyr (from M1)
cargo test --manifest-path rust-old/Cargo.toml --all-features   # the Rust oracle, read-only
```

```sh
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check   # golden corpus + generated .cyr up to date
python3 -I -B scripts/oracle/validate.py .                                         # independent check of the corpus
python3 -I -B scripts/oracle/mathx_tables.py --check                               # src/mathx_data.cyr current
python3 -I -B scripts/vendor_fmtf.py --check                                      # src/fmtf.cyr is impetus's, renamed
OUT=/tmp/ushma-gates bash scripts/gates.sh                                         # every gate, every target
```

## Development Process

### P(-1): Scaffold Hardening (before any new features)

0. Read the roadmap, the port plan, the CHANGELOG and the open issues — know what was intended
   before auditing what was built
1. Test + benchmark sweep of existing code
2. Cleanliness check (below)
3. Get baseline benchmarks
4. Internal deep review — gaps, optimizations, security, logging/errors, docs
5. External research — domain completeness, missing capabilities, best practices, world-class accuracy
6. Cleanliness check — must be clean after review
7. Additional tests/benchmarks from findings
8. Post-review benchmarks — prove the wins
9. Repeat if heavy

### Cleanliness check

```sh
cyrius fmt --check <files>
cyrius lint <files>
cyrius vet src/main.cyr
cyrius doc --check <files>      # every public fn documented
cyrius build src/main.cyr build/ushma
cyrius test                     # must be warning-free, not just passing
cyrius test --bounds
cyrius bench
cyrius coverage --exec          # target ≥ 80 %
cyrius distlib --all && cyrius distlib --check   # once [lib] exists (M1)
cyrius deps --locked
cyrius audit
```

`cyrius test` must be **warning-free**. The f64 type diagnostics are warnings, not errors, and
each one is a bug:
- *integer arithmetic with an f64 right operand*
- *f64 arithmetic with a non-f64 right operand*
- *an integer stored into an f64 slot*
- *comparison mixes f64 and integer*

### Work Loop / Working Loop (continuous)

1. Work phase — one module at a time, in the order the port plan gives
2. Cleanliness check
3. Test + benchmark additions: golden parity rows, single-step rows, reference tests, JSON rows, benches
4. Run benchmarks, recorded against the Rust numbers (`docs/benchmarks-rust-v-cyrius.md`)
5. Internal review — performance, memory (flat heap), security, throughput, correctness, and a
   surface-parity diff against `rust-old/src/`
6. Cleanliness check — must be clean after audit
7. Deeper tests/benchmarks from audit observations
8. Run benchmarks again — prove the wins
9. If audit heavy → return to step 5
10. Documentation — CHANGELOG, `state.md`, roadmap, the port plan's tracking table
11. Version check — `VERSION` and `cyrius.cyml` (`${file:VERSION}`) in sync
12. Return to step 1

### Task Sizing

- **Low/Medium effort**: Batch freely — multiple items per work loop cycle
- **Large effort**: Small bites only — break into sub-tasks, verify each before moving to the next. Never batch large items together
- **If unsure**: Treat it as large. Smaller bites are always safer than overcommitting

### Refactoring

- Refactor when the code tells you to — duplication, unclear boundaries, performance bottlenecks
- Never refactor speculatively. Wait for the third instance before extracting an abstraction
- Refactoring is part of the work loop, not a separate phase. If a review (step 5) reveals structural issues, refactor before moving to step 6
- Every refactor must pass the same cleanliness + benchmark gates as new code

## Key Principles

- **Everything ports.** Every Rust module, feature, test, benchmark and example has a Cyrius
  counterpart. Nothing is dropped or deferred. If something seems to have no Cyrius counterpart,
  grep `lib/` and the Cyrius docs first, then show the evidence and ask.
- **Cyrius-native, not Rust-shaped.** Rust is the oracle for *results*, not the design. Port each capability the
  way a Cyrius library builds it (stdlib and sibling idioms), never Rust's mechanics: no thiserror-style formatted
  error strings, no trait or derive plumbing recreated for its own sake. Decide the Cyrius form, record it in an ADR,
  and compare against Rust only on results.
- **The Rust code is frozen.** `rust-old/` is never edited and never re-released. When the
  oracle is wrong, the port does it right and an ADR 0002 row maps the oracle's rows.
- **Cross-check against `rust-old/`.** The correctness bar is "matches what Rust did", within the
  tolerance class the port plan (§5) assigns. Diverge only with an ADR.
- **Stdlib- and hisab-first.** Use the Cyrius stdlib (math, ganita, sakshi, log, sandhi, bayan, …)
  and hisab for everything they provide. Where they fall short, file the gap upstream with the
  numbers that show it (port plan §8).
- **Accuracy over glibc parity.** Transcendentals are correctly rounded and give the same bits
  on every target (ADR 0004). Never pick a primitive because it matches glibc's bits.
- **Never skip benchmarks.** Numbers don't lie. The benchmark record is the proof. Speed is
  measured, not gated.
- **Tests + benchmarks are the way.** Minimum 80%+ coverage target (`cyrius coverage --exec`).
- **Own the stack.** If an AGNOS library wraps an external capability, depend on the AGNOS library.
- **No magic.** Every operation is measurable, auditable, traceable.
- **Enums are append-only** with explicit values that are never renumbered (the Cyrius form of
  `#[non_exhaustive]`). Variant names are global, so they carry the `US_` prefix. Property
  functions `match` with an `_ =>` arm.
- **`#must_use` on all pure functions**, plus `#pure` where it applies.
- **`#inline` on hot paths** where a benchmark shows it helps.
- **`str_builder` into caller buffers over temporaries** (Rust's "`write!` over `format!`").
- **Borrow pointers over copies** (Rust's "Cow over clone"). Stepped state is passed as `p: *T`
  and called with `&x`, because a by-value struct over 8 B is a copy.
- **Column tables and sorted arrays over `hashmap`** (Rust's "Vec arena over HashMap"). Never
  iterate a `hashmap`: its seed is random per process.
- **Profiles for optional stacks.** `[lib]` is the physics. `[lib.ai]` adds the sandhi/TLS client.
- **Spans on every operation Rust instrumented.** sakshi spans at DEBUG under the Rust span
  names, and `log_*` events.
- **No aborts in library code.** `src/` never uses `+?`, `-?` or `*?` (they exit the process),
  `assert_fatal`, `panic`, `unwrap` or `result_unwrap`. Errors return a value-form `Result` with
  a `US_ERR_*` code.
- **f64 hygiene.**
  - Bind every builtin/call result to a `: f64` local before using it in an operator or comparison.
  - The left operand decides the operator's type: write `0.0 - x` / `x * 2.0`, never `0 - x` / `2 * x`.
  - Pass f64 literals as `3.0`, never `3`, because an integer keeps its bits in an f64 slot.
  - Never put hex in a `const`: it is an integer.
  - No exponent notation in literals; generated files write decimals positionally.
  - Use the canonical `ganita_f64_*` names, not the deprecated aliases.
- **Never test an f64 for truth** (`if (x)`, `!x`): it tests the raw bits (−0.0 is true) and the compiler doesn't
  flag it (probe 12). Compare explicitly.
- **Spans are level-gated.** `sakshi_span_enter` has no level check (192 ns per call even into its ring buffer,
  probe 8), so instrumented functions use `_us_span_enter` / `_us_span_exit`, which test the level first (+1.6 ns
  with DEBUG off).
- **Results through out-pointers.** A wide result fills the caller's `out: *T` and the function returns a value-form
  `Result` (Ok(0) / Err(code)). Never a freshly allocated record in the payload, which leaks 56 B per call (probe 3).
- **No per-call allocation on repeated paths.** Point lists go into a caller buffer or caller vecs that are refilled
  (`vec_truncate`); stepped state keeps its scratch, or resets an arena each step (probes 2, 4). The global heap
  never frees.
- **Toolchain workarounds** (open issues at the pin; port plan §2.8).
  - Bind a struct-returning call to a local before storing it into a field (`s.f = g(t.e)` miscompiles).
  - Never put an attribute (`#must_use`, `#inline`) directly before a `var` or `const`; it
    attaches to the next fn.
  - Wide results go through out-pointers, never a struct `Result` payload.
  - Never `alloc_reset`; per-step scratch lives in an arena.
- **Constants and tables are generated, never typed by hand.** `scripts/oracle/` writes them from
  the oracle. Regenerate them; never edit them.
- **Names.** Public names carry the `us` namespace (`us_*` fns; `US_*` constants, enum members,
  error codes), by the mechanism ADR 0001 records (P6). Private helpers use `_us_`.
- Build with `cyrius build`, not raw `cat file | cc5` — the manifest auto-resolves deps
- Source files only need project includes — stdlib auto-resolves from `cyrius.cyml`
- `var buf[N]` = N **bytes**, not N entries

## DO NOT

- **Do not commit or push** — the user handles all git operations (commit, push, tag)
- **NEVER use `gh` CLI** — use `curl` to GitHub API only
- **Do not modify `rust-old/`** — it is the frozen parity oracle. Tools that drive it live in
  `scripts/oracle/` and depend on it by path, read-only. No Rust fixes, no Rust releases
- **Do not drop or defer Rust functionality** — everything ports
- Do not modify `lib/` files (vendored stdlib / dep symlinks)
- Do not hardcode toolchain versions in CI YAML — `cyrius = "X.Y.Z"` in `cyrius.cyml` is the source of truth
- Do not widen a parity tolerance class to make a row pass
- Do not add unnecessary dependencies — keep it lean
- Do not abort or panic in library code
- Do not skip benchmarks before claiming performance improvements
- Do not skip tests before claiming changes work
- Do not commit `build/` or `rust-old/target/`. `rust-old/Cargo.lock` and `scripts/oracle/Cargo.lock`
  **are** committed: the oracle must rebuild against its exact pins

## Documentation Structure

```
Root files (required):
  README.md          — quick start, features, dependency stack, consumers, license
  CHANGELOG.md       — per-version changes (Added/Changed/Fixed/Removed)
  CLAUDE.md          — this file (development process, principles, DO NOTs)
  CONTRIBUTING.md    — fork, branch, cleanliness check, PR workflow
  SECURITY.md        — supported versions, scope, reporting
  CODE_OF_CONDUCT.md — Contributor Covenant
  LICENSE            — GPL-3.0

docs/:
  development/ushma-port-plan.md — the port plan of record
  development/state.md           — live state
  development/roadmap.md         — completed items, backlog, future features, v1.0 criteria
  development/port-parity.md     — measured parity baselines (M0 step 8)
  benchmarks-rust-v-cyrius.md    — Rust vs Cyrius benchmark record
  adr/NNNN-title.md              — architectural decision records (template.md)
  architecture/                  — module map, data flow, non-obvious constraints
  guides/                        — task-oriented how-tos
  examples/                      — runnable Cyrius examples (M8)
rust-old/                        — the frozen Rust 1.4.0 crate, its docs, benchmarks and CI
```

ADR format:
  # NNNN — Title
  ## Status: Accepted/Superseded
  ## Context: Why this decision was needed
  ## Decision: What we chose
  ## Consequences: Trade-offs, what changes

## CHANGELOG Format

Follow [Keep a Changelog](https://keepachangelog.com/):

```markdown
# Changelog

## [Unreleased]
### Added — new features
### Changed — changes to existing features
### Fixed — bug fixes
### Removed — removed features
### Security — vulnerability fixes
### Performance — benchmark-proven improvements (include numbers)
```

Rules:
- Every change that alters behavior gets a CHANGELOG entry
- Performance claims MUST include benchmark numbers (Rust vs Cyrius where relevant)
- Breaking changes get a **Breaking** section with a migration guide
- Group by module when multiple changes land in one release
- Link the ADR when a change was driven by an architectural decision
