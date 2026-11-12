# `ganita_mat_lu_solve` allocates its forward-substitution vector on every call, so a repeated solve leaks — open

**Status:** 🟡 **OPEN**. Verified 2026-11-11 on ganita 1.2.15 (the cyrius 6.7.18 bundle, x86_64). The repro exits 1.
**Placement:** unpinned.
**Discovered:** 2026-11-11, by ushma's Rust → Cyrius port (M0 probe 9, `ushma/docs/development/port-parity.md` §3).
**Reporter:** ushma 2.0.0 port (`ushma/docs/development/ushma-port-plan.md` §8; ushma ADR 0003 amendment)
**Related:** `2026-11-11-mat-inv-allocates-per-call.md` (pravash), the same family. ganita's own comment at
`ganita_mat_inv` (`lib/ganita.cyr:1017–1019` in the 6.7.18 snapshot) already names this allocation: 1.2.6 stopped
`mat_inv` calling `lu_solve` because of it. cyrius S1 (`2026-11-09-no-scoped-allocator-for-alloc.md`) would also
contain it.
**Severity:** Medium for a repeated solve. Nothing fails, but memory grows by 8·n bytes per call and is never
returned.
**Affects:** ganita 1.2.15, `ganita_mat_lu_solve` (`var y = alloc(n * 8)`, `lib/ganita.cyr:938` in the snapshot).

## Summary

`ganita_mat_lu_solve(lu, piv, b, out_x)` takes the factors, the pivots, the right-hand side and the solution from
the caller, then allocates its own `y` for the forward substitution on every call. The global allocator never frees,
so a program that factors once and solves many times grows the heap linearly. At n = 40 that is **320 B per call**,
with every caller buffer reused.

ushma's `ThermalNetwork::solve` moved to `ganita_mat_lu` + `ganita_mat_lu_solve` at M0. It is bit-identical to Rust
hisab 1.4.0 on 2,000 row-scaled network systems, and refuses the same 178. So `lu_solve` is ushma's network solver,
and a consumer that re-solves a network (a parameter sweep, an optimiser) leaks through it.

## Reproduction

`repros/2026-11-11-mat-lu-solve-allocates-per-call.cyr` factors a 40 × 40 tridiagonal matrix once, then calls
`ganita_mat_lu_solve` 1,000 times. It prints `heap bytes per ganita_mat_lu_solve at n = 40: 320` and exits 1.

## Root cause

`ganita_mat_lu_solve` allocates `y` with `alloc(n * 8)` after its argument checks. Nothing frees it, and there is no
caller-scratch or allocator-taking form.

## Proposed fix

Either form gives the same bits:
- `ganita_mat_lu_solve_a(a, lu, piv, b, out_x)`, the stdlib's `_a` convention, so the caller routes the scratch to an
  arena it resets;
- `ganita_mat_lu_solve_into(lu, piv, b, out_x, scratch)`, with `scratch` holding n words, matching pravash's
  `mat_inv_into` ask.

Alternatively, `lu_solve` needs no `y` when `out_x` is not `b`. The forward substitution can write into `out_x`, and
the back substitution then overwrites it in place: it reads `y[i]` before writing `x[i]`, and reads `x[j]` only for
j > i. That is the same operations in the same order, so the same bits, with no allocation and no new name. The
catch is aliasing. Today `out_x == b` works, because the forward pass finishes reading `b[piv[i]]` before anything
is written. An in-place forward pass would overwrite `b` first, so that case still needs a scratch vector.

## Consumer-side workaround

ushma's plan: until the fix lands, M6 measures a local forward and back substitution over `ganita_mat_lu`'s factors
against `ganita_mat_lu_solve`, bit for bit. If it matches, ushma uses it as a named stopgap.
