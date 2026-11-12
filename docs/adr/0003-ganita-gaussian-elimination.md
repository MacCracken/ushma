# 0003 — `ThermalNetwork::solve` on ganita's Gaussian elimination

**Status**: Accepted
**Date**: 2026-11-11

## Context

Rust 1.4.0 solves the thermal network's conductance system with Rust hisab 1.4.0's
`num::gaussian_elimination`. That routine uses partial pivoting and calls a pivot singular when
it is below `EPSILON_F64` in absolute value. ushma 1.4.0 scales each row by its diagonal first,
so that tiny conductances (R ≈ 1e14 K/W) don't trip that absolute threshold. Cyrius hisab 3.4.1
has no dense direct solver, but the stdlib's ganita 1.2.15 has two:
- `ganita_mat_gaussian_elim(aug)`: partial pivoting, in place on [A|b], a **relative**
  singularity threshold over the whole augmented matrix, and refusal of any non-finite entry;
- `ganita_mat_lu` + `ganita_mat_lu_solve`.

## Decision

`ThermalNetwork::solve` uses `ganita_mat_gaussian_elim`, keeping 1.4.0's diagonal row scaling in
front of it. Nothing is transcribed.

- **Parity.** If the pivot choices and operation order match Rust hisab's, the rows are class B.
  If they don't, class D. M0 probe 9 measures this on ushma's network matrices.
- **Singularity.** ganita decides singularity by its relative threshold, not by Rust's absolute
  `ε`. On near-singular inputs the two can disagree. Those rows are accepted here, and the
  harness reports them by name.
- **Non-finite input.** ganita refuses it, where Rust hisab would have propagated NaN. ushma
  1.4.0 validates its inputs before solving, so this should never be reached. A probe-9 row
  confirms it.

## Consequences

- **Positive** — no private solver to own; ushma uses the stdlib's, which is audited for
  non-finite input.
- **Negative** — near-singular networks may get an `Err` where Rust returned a value (or the
  reverse); each such row is named in the parity report.
- **Neutral** — kimiya planned to transcribe the same routine (its ADR 0003), so ushma tells
  kimiya ganita already has one (plan §8).

## Alternatives considered

- **Transcribe Rust hisab 1.4.0's routine** into a private `_us_*` solver for exact class B.
  Rejected: the stdlib already has the solver, and the port is stdlib-first.
- **`ganita_mat_lu` + `ganita_mat_lu_solve`.** Kept as the fallback if probe 9 shows
  `ganita_mat_gaussian_elim` diverges more.
