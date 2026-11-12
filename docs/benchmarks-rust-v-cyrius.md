# Benchmarks: Rust 1.4.0 against the Cyrius port

Speed is measured, not gated (CLAUDE.md). Each module's benchmarks join this record at its milestone.

- **Rust:** criterion 0.8 on the frozen crate at tag `1.4.0` (`rust-old/bench-history.csv`, commit `2688099`, run
  2026-11-12T02:59Z): 63 benchmarks.
- **Cyrius:** `cyrius bench tests/ushma.bcyr` on cyrius 6.7.18, x86_64 ELF. Each operation runs in a batch of 10⁶
  (10⁵ for the JSON ones) through an indirect call (`fncall1`), so every Cyrius number includes the harness floor
  (the `noop` row) that criterion's inlined loop doesn't have.
- **Machine:** AMD Ryzen 7 5800H, Linux 7.2, both columns.

## M1

The Cyrius column is after M1d's review; "before" is the same benchmark before it (M1c's run for part 1, M1d's
first run for part 2).

| Benchmark | Rust 1.4.0 (ns) | Cyrius (ns) | Before review (ns) | Notes |
|---|---|---|---|---|
| noop (harness floor) | — | 4 | 4 | an empty fn through `fncall1` |
| material/diffusivity | 0.25 | 7 | 7 | k / (ρ·c_p); Rust inlines it to a constant-time load and divide |
| material/volumetric_cp | 0.25 | 7 | 7 | ρ·c_p |
| material/to_json | — | 4,567 | — | Cyrius only. Five f64 through fmtf's exact bignum printer (ADR 0005) plus the name; printing is off every step path |
| material/from_json | — | 4,351 | — | Cyrius only. bayan's tree on an arena reset per call, then serde's rules |
| span, DEBUG off | — | 11 | 11 | Cyrius only. `_us_span_enter` + `_us_span_exit`: two level tests (M0 probe 8 measured 1.6 ns of span on top of its function) |
| transfer/conduction | 13.1 | 23 | 31 | two finite checks and a level-gated span. Rust's 13 ns includes its own span |
| transfer/convection | 1.24 | 8 | 8 | h·A·ΔT, no checks |
| transfer/radiation | 9.95 | 25 | 50 | three checks, T⁴ twice, a span |
| transfer/resistance_series | 2.19 | 15 | 16 | five terms |
| transfer/resistance_parallel | 9.28 | 31 | 47 | three terms: a check pass, then the reciprocals |
| transfer/heat_stored | 0.99 | 9 | 9 | m·c·ΔT |
| transfer/diffusivity | 4.49 | 16 | 33 | four checks and a division |
| transfer/biot_number | 3.99 | 15 | 29 | |
| transfer/lumped_capacitance | 5.59 | 54 | 54 | correctly rounded exp (ADR 0004); glibc's exp is not correctly rounded |
| transfer/reynolds | 2.57 | 10 | 16 | |
| transfer/nusselt_db | 24.9 | 287 | 297 | two correctly rounded `pow` calls (ADR 0004), about 140 ns each, against glibc's ~12 ns |
| transfer/fin_rect_heat | 89.76 | 157 | 180 | four checks, two √, a correctly rounded tanh, a span |
| transfer/lmtd_counter | 16.64 | 96 | 106 | two checks, a correctly rounded ln, a span |
| transfer/eff_counter | 8.55 | 74 | 84 | the range checks and a correctly rounded exp |
| transfer/view_factor_disks | 8.19 | 44 | 65 | three checks, four divisions, a √ |
| transfer/view_factor_perpendicular | 54.41 | 590 | 608 | three correctly rounded atan and three ln, a span |

**M1d's review.** Two changes, each bit-identical on every golden row:

- **The finite checks.** `_us_finite_gt` / `_us_finite_ge` were an `if` around a call to `_us_is_finite`. Cyrius
  inlines no fn with control flow or a `var`, so each check cost two calls. They are now one branch-free expression,
  `(x > bound) & (x < +∞)`, which is the whole test for a finite bound (NaN and −∞ fail the first comparison, NaN and
  +∞ the second), and `#inline` takes it. Three checks dropped from 21 ns to 6 (harness floor subtracted). Every
  checked function gained: diffusivity 33 → 16 ns, Biot 29 → 15, parallel resistance 47 → 31. The view factors
  test their three dimensions inline and call the reporting helper only to refuse (disks 50 → 44).
- **Rust's `powi(2)` and `powi(4)`** are x·x and (x·x)·(x·x), the products `__powidf2` forms, so they are written
  as those products instead of the loop. Radiation went from 40 to 25 ns.

**What is left.** The rows with a transcendental pay for correct rounding, which is the accuracy target (ADR 0004).
The rest of the gap is the harness floor (4 ns a call) and, in the seven instrumented functions, the level-gated span
(7 ns for the pair with DEBUG off, the span row less the floor).

## M2

| Benchmark | Rust 1.4.0 (ns) | Cyrius (ns) | Notes |
|---|---|---|---|
| state/ideal_gas_pressure | 9.58 | 24 | three checks, nRT/V, a level-gated span |
| state/ideal_gas_volume | 3.74 | 16 | three checks, nRT/P |
| state/van_der_waals | 12.54 | 38 | five checks and a refusal test on V − nb, a span |
| state/isothermal_work | 9.42 | 83 | four checks, a correctly rounded ln |
| state/adiabatic_temperature | 17.15 | 155 | four checks, a correctly rounded pow |
| state/compressibility_factor | 4.46 | 21 | five checks, PV/(nRT) |
| state/redlich_kwong | 27.27 | 188 | four checks, a correctly rounded T_c^2.5 (ADR 0004), √T, a span |
| state/peng_robinson | 18.85 | 69 | five checks, √(T/T_c), the attractive denominator's check, a span |
| entropy/carnot_efficiency | 8.69 | 20 | two checks, the T_cold < T_hot test, a span |
| entropy/carnot_cop | 3.30 | 14 | three checks |
| entropy/helmholtz | 0.99 | 9 | U − TS |
| entropy/gibbs | 0.99 | 9 | H − TS |
| entropy/mixing | 25.41 | 233 | three fractions: the range and sum checks, three correctly rounded ln, a span |
| entropy/ideal_gas_ds | 15.97 | 151 | six checks, two correctly rounded ln |

**M2's review.** Every check is the inlined branch-free form from M1d, no state or entropy function allocates (Kay's
rules and mixing read the caller's lists; the JSON readers work on the caller's arena), and the gap to Rust is the
correctly rounded ln and pow (ADR 0004), the harness floor and the spans. Nothing in M2 measured worth changing.

## M3

| Benchmark | Rust 1.4.0 (ns) | Cyrius (ns) | Notes |
|---|---|---|---|
| phase/clausius_clapeyron_slope | 3.27 | 19 | three checks, L / (T·Δv) |
| phase/clausius_clapeyron_pressure | 10.62 | 92 | four checks, a correctly rounded exp, the finite check |
| phase/phase_lookup | 22.57 | 205 | water at 300 K, 1 atm: the vapor curve's correctly rounded ln and exp, and the span with its substance field |
| phase/heat_of_fusion | 1.08 | 9 | m·L_fus |
| phase/heat_of_vaporization | 0.99 | 9 | m·L_vap |
| phase/heat_for_phase_change | 11.79 | 49 | four checks, the crossings, the span with its substance field |
