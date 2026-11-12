# Benchmarks: Rust 1.4.0 against the Cyrius port

Speed is measured, not gated (CLAUDE.md). Each module's benchmarks join this record at its milestone.

- **Rust:** criterion 0.8 on the frozen crate at tag `1.4.0` (`rust-old/bench-history.csv`, commit `2688099`, run
  2026-11-12T02:59Z): 63 benchmarks.
- **Cyrius:** `cyrius bench tests/ushma.bcyr` on cyrius 6.7.18, x86_64 ELF. Each operation runs in a batch of 10⁶
  (10⁵ for the JSON ones) through an indirect call (`fncall1`), so every Cyrius number includes the harness floor
  (the `noop` row) that criterion's inlined loop doesn't have.
- **Machine:** AMD Ryzen 7 5800H, Linux 7.2, both columns.

## M1

| Benchmark | Rust 1.4.0 (ns) | Cyrius (ns) | Notes |
|---|---|---|---|
| noop (harness floor) | — | 4 | an empty fn through `fncall1` |
| material/diffusivity | 0.25 | 7 | k / (ρ·c_p); Rust inlines it to a constant-time load and divide |
| material/volumetric_cp | 0.25 | 7 | ρ·c_p |
| material/to_json | — | 4,567 | Cyrius only. Five f64 through fmtf's exact bignum printer (ADR 0005) plus the name; printing is off every step path |
| material/from_json | — | 4,351 | Cyrius only. bayan's tree on an arena reset per call, then serde's rules |
| span, DEBUG off | — | 11 | Cyrius only. `_us_span_enter` + `_us_span_exit`: two level tests (M0 probe 8 measured 1.6 ns of span on top of its function) |
| transfer/conduction | 13.1 | 31 | two finite checks and a level-gated span. Rust's 13 ns includes its own span |
| transfer/convection | 1.24 | 8 | h·A·ΔT, no checks |
| transfer/radiation | 9.95 | 50 | three checks, two `powi(4)`, a span |
| transfer/resistance_series | 2.19 | 16 | five terms |
| transfer/resistance_parallel | 9.28 | 47 | three terms: a check pass, then the reciprocals |
| transfer/heat_stored | 0.99 | 9 | m·c·ΔT |
| transfer/diffusivity | 4.49 | 33 | three checks and a division |
| transfer/biot_number | 3.99 | 29 | |
| transfer/lumped_capacitance | 5.59 | 54 | correctly rounded exp (ADR 0004); glibc's exp is not correctly rounded |
| transfer/reynolds | 2.57 | 16 | |
| transfer/nusselt_db | 24.9 | 297 | two correctly rounded `pow` calls (ADR 0004), about 145 ns each, against glibc's ~12 ns |

**Where the time goes (M1c).** Cyrius doesn't inline a fn with control flow, so each finite check is two calls
(`_us_finite_ge` → `_us_is_finite`), about 5 ns each. That accounts for most of the gap in the checked functions. The
transcendental rows pay for correct rounding, which is the accuracy target (ADR 0004). M1d's review measures the
checks folded into one branch-free predicate.
