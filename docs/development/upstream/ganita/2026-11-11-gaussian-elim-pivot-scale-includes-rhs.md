# `ganita_mat_gaussian_elim` scales its singularity threshold by the right-hand side too, so b's units decide whether A is singular — open

**Status:** 🟡 **OPEN**. Verified 2026-11-11 on ganita 1.2.15 (the cyrius 6.7.18 bundle, x86_64). The repro exits 1.
**Placement:** unpinned.
**Discovered:** 2026-11-11, by ushma's Rust → Cyrius port (M0 probe 9, `ushma/docs/development/port-parity.md` §3).
**Reporter:** ushma 2.0.0 port (`ushma/docs/development/ushma-port-plan.md` §8; ushma ADR 0003 amendment)
**Related:** `archived/2026-09-07-linalg-eps-absolute-tolerance.md`, resolved in 1.2.4. It made the tolerance
relative to "the matrix scale for the pivot tests in `lu` and `gaussian_elim`". For `gaussian_elim` that scale is
taken over the augmented matrix, so the 1.2.4 problem (a verdict that depends on the caller's units) is back
through b.
**Severity:** Medium. A well-posed system is refused, `0` with no solution, depending on the magnitude of its
right-hand side. `ganita_mat_lu` decides the same matrix correctly.
**Affects:** ganita 1.2.15, `ganita_mat_gaussian_elim`:
`var ge_tol = f64_mul(LINALG_EPS, _linalg_mat_scale(aug))` (`lib/ganita.cyr:1375` in the snapshot), against
`ganita_mat_lu`'s `_linalg_mat_scale(out_lu)` over A alone (:851).

## Summary

Whether A is singular is a property of A, but `gaussian_elim` measures its pivots against
1e-12 · max|[A|b]|. Multiply b by 1,000, by writing temperatures in K instead of kK, and the same A becomes singular.

ushma's thermal-network systems show the effect at scale. M0 probe 9 ran 2,178 seeded connected networks, each system
row-scaled exactly as ushma 1.4.0 builds it (unit diagonal, |aᵢⱼ| ≤ 1, b = conductance-weighted temperatures up to
1,000 K):

| | Rust hisab 1.4.0 | `ganita_mat_lu` + `_lu_solve` | `ganita_mat_gaussian_elim` |
|---|---|---|---|
| 2,000 systems Rust solves | solved | 2,000, bit-identical | 1,741 bit-identical, **259 refused** |
| 178 systems Rust refuses | refused | 178 refused | 178 refused |

On these matrices LU's relative test, 1e-12 · max|A| = 1e-12, is exactly Rust hisab's absolute one. Elimination's
threshold is up to 1,000 times larger.

## Reproduction

`repros/2026-11-11-gaussian-elim-pivot-scale-includes-rhs.cyr` uses A = [[1, −(1 − 1e-10)], [−1, 1]]: unit diagonal,
smallest pivot 1e-10, cond ≈ 4e10. Its output:

```
ganita_mat_lu factors A: 1
ganita_mat_gaussian_elim, b = (1, 0) (kK): 1
ganita_mat_gaussian_elim, b = (1000, 0) (K): 0
```

It exits 1.

## Root cause

`ge_tol` is computed from `_linalg_mat_scale(aug)`, every entry of the n × (n + 1) augmented matrix, column n
included.

## Proposed fix

Scale over the first n columns only, as `lu` does: the largest |aᵢⱼ| with j < n. The repro's two verdicts then agree,
and ushma's 259 refusals become solutions bit-identical to LU's. (Probe 9 found elimination bit-identical wherever it
solved.) Add both directions as assertions, as 1.2.4's fix did: a system that must not be called singular whatever
b's scale, and a singular one that must be called singular whatever b's scale.

## Consumer-side workaround

ushma uses `ganita_mat_lu` + `ganita_mat_lu_solve` (its ADR 0003 amendment), which isn't affected.
