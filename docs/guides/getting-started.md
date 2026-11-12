# Getting started with ushma

## Build

```sh
cyrius deps                              # resolve dependencies
cyrius build src/main.cyr build/ushma    # compile
cyrius test                              # run tests/*.tcyr
```

## Layout

- `src/main.cyr` — smoke entry point (proves the build links). The library ships as bundles:
  `[lib]` → `dist/ushma.cyr`, `[lib.ai]` → `dist/ushma-ai.cyr` (from M1; `cyrius distlib --all`).
- `src/*.cyr` — one Cyrius module per Rust module, in the `[lib]` order of `cyrius.cyml`.
- `tests/` — `.tcyr` tests (parity + reference), `.bcyr` benchmarks, `.fcyr` fuzz targets.
- `lib/` — vendored stdlib leaves and hisab (`cyrius deps`). Never edit.
- `rust-old/` — the frozen Rust 1.4.0 crate, the parity oracle. Never edit.
- `scripts/oracle/` — the golden-value generator that drives `rust-old/` (M0 part 2).

## Adding a module

Follow the port plan's milestone for it ([`../development/ushma-port-plan.md`](../development/ushma-port-plan.md) §6):

1. Port the module from `rust-old/src/` into `src/<module>.cyr` and add it to `[lib]`.
2. Generate its golden rows from the oracle and check them with the parity harness.
3. Port its Rust tests as reference tests, and its benchmarks as Cyrius twins.
4. Run the cleanliness check (CLAUDE.md). `cyrius test` must be warning-free.
5. Add a CHANGELOG entry and update `state.md`.

See [`../adr/template.md`](../adr/template.md) when a non-trivial design choice deserves an ADR.
