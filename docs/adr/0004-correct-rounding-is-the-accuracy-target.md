# 0004 — Correct rounding is the accuracy target

**Status**: Accepted
**Date**: 2026-11-11

## Context

The oracle's transcendentals come from glibc 2.44 through Rust std: `powf` (29 calls), `ln`
(23), `exp` (14), `atan` (11), `tanh` (2) and `exp_m1` (2). glibc is not correctly rounded
everywhere. The Cyrius stdlib offers:
- the `f64_exp`/`f64_ln`/`f64_atan` builtins: x87 on x86 and a polyfill on aarch64, so the bits
  differ by target. The open `exp-polyfill-double-rounds-subnormal-results` issue measures up to
  39 ulp on aarch64 subnormals, and about 5 % misrounding on x87;
- `_f64_exp_polyfill`, `_f64_ln_polyfill` and `_f64_atan_polyfill`: ≤ 1 ulp, the same bits
  everywhere;
- `ganita_f64_pow`, `ganita_f64_tanh` and `ganita_f64_cbrt` (correctly rounded);
- **no expm1**.

The maintainer's standing policy, set in kimiya (its ADRs 0005 and 0006), is to be better than
glibc. Primitives are correctly rounded, and parity with the oracle is functional, not
bit-copying glibc's errors.

## Decision

1. **Every transcendental primitive ushma calls is correctly rounded** and target-independent:
   exp, expm1, ln, pow, atan (as atan2(x, 1), which equals atan exactly) and tanh. They live in
   `src/mathx.cyr` as private `_us_cr_*` fns, vendored from kimiya's `src/mathx.cyr` (Ziv
   strategy: a fast path with a proven bound and a rounding test, then an accurate path).
   `ganita_f64_cbrt` is used directly, since it is already correctly rounded. `ganita_f64_tanh`
   is used only if M0 step 8 shows it correctly rounded; otherwise `_us_cr_tanh`.
2. **The class-C reference is the corrected oracle.** Rust's operation sequence evaluated with
   correctly rounded primitives (MPFR via `mpfr_ref.py`). Rows where glibc misrounds are
   generated into `tests/golden/<module>.cr.tsv` and counted "accepted under ADR 0004" per
   function. A port value equal to glibc's misrounded value on such a row is a failure.
3. **Class C is compared exactly** once a function is on the `_us_cr_*` primitives: 0 ulp
   expected.
4. **Functional parity is unchanged.** The port keeps the oracle's results, `Err` codes, check
   order and branches. A corrected primitive that changes a branch makes the corrected oracle
   follow the other branch, and the port must too.
5. **Range edges and exact constants are swept per milestone (P2)**, not in a library-wide
   audit afterwards. That covers spurious overflow or underflow of products and quotients, and
   exact rationals Rust misrepresents (`1.0/3.0`, `1.0/6.0` exponents). Each family gets an
   ADR 0002 row.

## Amendment, 2026-11-11 (M0 part 2c-3): tanh

Step 8 measured `ganita_f64_tanh` against the corrected oracle. It is off by 1 ulp on 52 of the 860 golden rows,
where glibc is correctly rounded, so it is not used. `_us_cr_tanh` is ushma's own correctly rounded tanh, in
`src/mathx.cyr` after the vendored code. It computes tanh a = E/(E + 2) with E = expm1(2a), on the vendored expm1's
arithmetic. The fast path is a double-double with a proven bound and a rounding test (the map's relative
sensitivity 2/(E + 2) ≤ 1 keeps E's error from growing). The accurate path is a triple-double to 2^-149, then RN3.
- **Below 2^-27:** x itself.
- **Above `TANH_ONE`:** ±1. 19.061547465398494, from MPFR, generated into `src/mathx_data.cyr`.

tanh is transcendental at every nonzero algebraic argument, so there are no exact cases or ties.
`scripts/oracle/cr_verify.py` checks it against MPFR on 1,082,051 arguments per target, the accurate path forced on
100,045; the fast path falls back on 1.3·10^-5 of them.

## Consequences

- **Positive** — results are the best a double can hold, and identical on x86_64 ELF, PE and
  aarch64.
- **Negative** — ushma owns a third private copy of the correctly rounded primitives (after
  impetus `crmath.cyr` and kimiya `mathx.cyr`), until upstream hosts them (plan §8).
- **Neutral** — step 8 measures the stdlib alternatives against MPFR too, recording what this
  ADR buys.

## Alternatives considered

- **The stdlib polyfills and ganita (≤ 1 ulp, same bits).** Rejected: not correctly rounded,
  so the port would be no better than glibc. There is also no expm1.
- **Extract the shared primitives upstream first.** Rejected for now: kimiya's open audit batches
  still change `mathx.cyr`. Filed as a request instead (plan §8).
