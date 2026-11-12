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
