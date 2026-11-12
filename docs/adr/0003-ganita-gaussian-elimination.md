# 0003 — `ThermalNetwork::solve` on ganita's Gaussian elimination

**Status**: Accepted; amended 2026-11-11 (probe 9): LU, not elimination
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

## Amendment, 2026-11-11 (M0 probe 9): LU

Probe 9 ran 2,178 seeded connected networks (R from 1e-3 to 1e14, systems built and row-scaled exactly as 1.4.0's
`solve` builds them; `port-parity.md` §3):
- **`ganita_mat_gaussian_elim`** refuses 259 of the 2,000 systems 1.4.0 solves. Its tolerance scales with the whole
  [A|b], and the right-hand side holds temperatures. Where it solves, it is bit-identical.
- **`ganita_mat_lu` + `ganita_mat_lu_solve`** is bit-identical on all 2,000 (37,424 unknowns), and refuses exactly the
  178 systems 1.4.0 refuses. Its threshold is 1e-12 · max|aᵢⱼ| over A alone. After 1.4.0's row scaling, every
  diagonal is exactly 1 and every other entry is at most 1 in magnitude, so that is 1.4.0's absolute
  `|pivot| < 1e-12`. The two differ only at a pivot of exactly 1e-12 (ganita's test is `≤`).

So **`ThermalNetwork::solve` uses `ganita_mat_lu` + `ganita_mat_lu_solve`** after the same row scaling, and its rows
are **class B** (0 ulp measured). The singularity divergence this ADR accepted above doesn't occur on the measured
systems. A pivot of exactly 1e-12 is the one remaining case, and it is named if a row ever hits it. The 178 systems
1.4.0 refuses are connected networks whose resistances span enough decades that a pivot falls below 1e-12 after row
scaling. The port keeps that decision. Solving them by a better-conditioned method would be a divergence for M6 to
propose (ADR 0002), not a parity row.

`ganita_mat_lu_solve` allocates its forward-substitution vector on every call (`alloc(n * 8)`, never freed: 320 B at
n = 40). The request is an allocator-taking `ganita_mat_lu_solve_a`, filed upstream (port plan §8). Until it lands,
M6 measures a local substitution over ganita's factors against `_lu_solve`, bit for bit: the same operations in the
same order. If the bits match, M6 uses it as a named stopgap, so that a consumer's repeated `solve` doesn't grow the
heap.

## Amendment, 2026-11-12 (M6c): the stopgap, measured

`us_network_solve` factors with `ganita_mat_lu` and runs the forward and back substitution itself, over ganita's
factors, in `_lu_solve`'s operations and order. The right-hand side, solution and forward vector live in the network's
scratch:
- **The measurement.** On 600 seeded networks (`tests/numerical.tcyr`, drawn as probe 9's), the network's own scaled
  system and right-hand side went through `ganita_mat_lu` + `ganita_mat_lu_solve`. 537 solved, and all 10,221
  unknowns matched the stopgap's bit for bit. Both refused the other 63. The golden corpus agrees: 398 solutions bit
  for bit, 69 refusals alike. So the stopgap holds, and a repeated solve allocates nothing.
- **ganita's checks.** `_lu_solve` refuses non-finite factors or a non-finite right-hand side. The stopgap catches
  both in the back substitution instead, with the same outcome. A NaN or ±∞ there reaches the solution: ∞ times any
  finite value is ±∞ or NaN, and each of U's diagonal entries is tested as it divides. So the solve is refused exactly
  when `_lu_solve` would refuse it, or would return a non-finite temperature.
- **Non-finite input is reachable.** The Decision above expected 1.4.0's validation to keep non-finite values out of
  the solver. It doesn't. 1/R overflows for a subnormal R (finite and > 0, so valid), and g·T overflows for a large T.
  1.4.0 then returns NaN or ±∞ temperatures as a success. The port refuses such a system (ADR 0002 row 16; 9 golden
  rows accepted by name).

When ganita ships an allocator-taking `_lu_solve_a` (port plan §8), the stopgap gives way to it. Its bits are the
same, so no golden row moves.

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
  `ganita_mat_gaussian_elim` diverges more. It did, and the amendment above adopts LU.
