# 0001 — Port everything to Cyrius, with Rust 1.4.0 frozen as the parity oracle

**Status**: Accepted
**Date**: 2026-11-11

## Context

AGNOS libraries are moving from Rust to Cyrius. ushma's math dependency (hisab 3.4.1) and the
services its `ai` module talks to (daimon 2.5.6, hoosh 2.12.3) are already Cyrius-only. ushma
1.4.0 is the result of a full correctness audit: IAPWS-IF97 steam tables, NaN/±∞ rejected by
every guard, and 80 regression tests with IF97, Incropera and NIST reference values. All 479
tests pass on rustc 1.99.1, and fmt, clippy, rustdoc, audit and deny are clean.

A translation needs an oracle it can check value by value, and only a repaired crate can be one.
The Cyrius 6.7.18 stdlib covers every Rust dependency ushma has:
- ganita: pow, tanh, atan2, cbrt, Gaussian elimination, LU;
- sakshi and log: spans, levels;
- sandhi: an HTTP/1.1 + HTTP/2 client;
- bayan: JSON;
- io: `getenv`.

## Decision

ushma 2.0.0 is a Cyrius library, and **everything ports**: all 18 Rust files and every feature,
including `bridge`, `integration/soorat`, `logging` and `ai`. It also carries serde as JSON
codecs, all 449 unit and 30 integration tests, the 63 benchmarks and the 5 examples. Nothing is
dropped or deferred.

- **The oracle is frozen.** The Rust 1.4.0 crate (tag `1.4.0`, commit `2688099`, whose Rust code
  is identical to `9528e1f`) moves to `rust-old/` as it stands. It is never edited and never
  re-released. When the oracle is wrong, the port does it right, and ADR 0002 maps the rows.
- **Initialization.** `cyrius port` on Cyrius 6.7.18, then the repairs in the port plan §4 step 2.
- **Scope, milestones, parity policy.** As in `docs/development/ushma-port-plan.md`, M0–M8, with
  parity classes A–D.
- **Stdlib- and hisab-first.** The port calls the Cyrius stdlib and hisab 3.4.1 for everything
  they provide, and files gaps upstream (plan §8).
- **Bundles (P5).** `[lib]` holds every module except the HTTP client, logging included.
  `[lib.ai]` is `[lib]` plus `src/ai.cyr`, so only a consumer that wants the daimon client links
  sandhi's TLS stack, as with Rust's `ai` feature.
- **Names (P6).** Cyrius has one global namespace. A duplicate non-private fn is a warning, and
  the last definition wins; enum variant names are global too. So public names carry the `us`
  namespace: `us_*` fns, `US_*` constants, enum members and error codes. The mechanism, the
  language's `mod us;` mangling or a written prefix, is chosen by M0 probe 6 and recorded here
  as an amendment.
- **Oracle mechanics.** Generators live in `scripts/oracle/` and depend on `rust-old/` by path,
  read-only. `rust-old/Cargo.lock` is committed, so the oracle always builds against hisab 1.4.0,
  serde_json 1.0.154 and reqwest 0.13.

## Amendment, 2026-11-11 (M0 probe 6): the written prefix

`mod us;` is unusable for a bundle. Its scope runs past the included file to the end of the compilation unit, so
it mangles the consumer's own functions (`main` became `us_main`). No form ends it, globals aren't mangled, and a
bare call between the module's own functions doesn't resolve without a `use`. So public names carry a **written
`us_` prefix**: `us_*` functions, and `US_*` constants, enum members and error codes. Private helpers use `_us_`.
The collision is real: kimiya's bundle defines `ideal_gas_pressure`, and including it beside an unprefixed ushma
is a duplicate-fn error, or a silent last-definition-wins when the signatures match.

## Consequences

- **Positive**
  - Correctness is measured, not asserted: every public function is diffed against the
    repaired Rust crate.
  - The audit's reference-value tests travel with the code.
  - No Rust consumer loses a feature by moving to the port.
- **Negative**
  - Two toolchains are needed until the port completes: cargo for the oracle and its generators,
    cyrius for the library.
  - The oracle's transcendentals come from the host glibc (2.44 at freeze), so the corpus header
    records it.
  - hisab is a production dependency: `steam` calls `calc_lerp` and `stat` calls
    `calc_integral_simpson`. Every consumer therefore vendors hisab and sees the resolver's
    benign sakshi warning until hisab drops its `[deps.sakshi]` (plan §2.7, §8). Probe 10
    measures the cost.
- **Neutral**
  - `rust-old/README.md` is the tagged README, restored there because `Cargo.toml` names it.
  - The Rust CI, scripts, benchmarks and architecture overview are kept under `rust-old/`.

## Alternatives considered

- **Port only the physics and drop or defer `ai`, `logging`, `bridge`, `soorat-compat` and
  serde.** Rejected by the maintainer (2026-11-11): every one has a Cyrius counterpart, and
  consumers would lose working features.
- **Fix oracle defects in Rust 1.4.x first.** Rejected by the maintainer (2026-11-11). The Rust
  code is legacy, and the services it talks to exist only in Cyrius.
- **Port without an oracle.** Rejected: hand-checked values drift, and the siblings' corpora
  caught defects their test suites would not have.
