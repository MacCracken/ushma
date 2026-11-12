
## Also needed: ushma (2026-11-12, cyrius 6.7.18): a fourth refusal, `if` and `var`

ushma's port hits a refusal this file doesn't list: **`#inline ignored: body has control flow or a var
declaration`**. It declines a body of 2 parameters and about a dozen tokens, inside the three limits above, for one
`if` or one local:

```text
warning:<source>:8:1: #inline ignored: body has control flow or a var declaration
warning:<source>:15:1: #inline ignored: body has control flow or a var declaration
```

**Where it costs** (ushma 2.0.0, `ushma/docs/benchmarks-rust-v-cyrius.md`, M1):
- **Range checks.** Every fallible function opens with one "finite and above a bound" check per argument. Written as
  an `if` around a call to an `is_finite` helper, three of them cost 21 ns; as one `if` in one fn, 12 ns; rewritten as
  the single expression `(x > b) & (x < INF)`, which `#inline` takes, 6 ns (the harness floor subtracted). The rewrite
  took ushma's thermal diffusivity from 33 to 16 ns (Rust 4.49), its Biot number from 29 to 15 (3.99) and its
  parallel resistance from 47 to 31 (9.28).
- **Level-gated spans.** `if (sakshi_get_level() >= SK_DEBUG) { sakshi_span_enter(name, len); }` has no single
  expression form, so each instrumented function pays 7 ns for its enter and exit with DEBUG off: 7 of
  `ideal_gas_pressure`'s 24 ns (Rust 9.58).

The expression rewrite works only for a predicate. A guard with an early return, a NaN-aware max (pravash's `_pv_max`
above) or a check that calls something has no branch-free form.

**Repro:** `repros/2026-11-11-inline-refuses-small-numeric-helpers-if-var.sh` builds the same check three ways (an
`if`, a local, one expression), confirms the three agree on 2.0 and NaN, and **exits 1 while the first two are
refused**. Run on cyrius 6.7.18, x86_64, 2026-11-12: exit 1, both warnings above.

**Ask (4):** let `#inline` take a body with `if` / `else`, early `return`s and local `var`s, inside the size limit:
each `return` a jump to the end of the inlined block, each local a slot in the caller's frame. With pravash's 1–3,
ushma's checks and span gates would inline as they are written.
