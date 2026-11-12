# Contributing to ushma

Thank you for considering a contribution! This document covers the development
workflow, coding standards, and review process.

ushma is being ported from Rust to Cyrius. The Rust 1.4.0 crate is frozen in `rust-old/`
as the parity oracle and is never edited. See
[docs/development/ushma-port-plan.md](docs/development/ushma-port-plan.md).

## Development Setup

```bash
# Clone the repository
git clone https://github.com/MacCracken/ushma.git
cd ushma

# Install the Cyrius toolchain pinned in cyrius.cyml ([package].cyrius), then:
cyrius deps          # vendor the stdlib leaves and hisab into lib/
cyrius build src/main.cyr build/ushma
```

The Rust oracle and the golden-value generator also need a Rust toolchain
(`rust-old/rust-toolchain.toml`).

## Pull Request Process

1. **Fork and branch** — create a feature branch from `main`.
2. **Keep commits focused** — one logical change per commit.
3. **Write tests** — new code requires golden parity rows and reference tests; bug fixes
   require a regression test.
4. **Run the cleanliness check locally** before pushing:
   ```bash
   cyrius fmt --check <files>
   cyrius lint <files>
   cyrius doc --check <files>
   cyrius test            # must be warning-free
   cyrius test --bounds
   cyrius bench
   cyrius deps --locked
   cyrius audit
   ```
5. **Open a PR** against `main` with a clear description of the change.
6. **Address review feedback** — maintainers may request changes before merging.

## Code Style

- `cyrius fmt` defaults (enforced by CI), and zero `cyrius lint` findings.
- `cyrius test` must be warning-free: every f64 type diagnostic is a bug.
- Public API items must have doc comments (`cyrius doc --check`).
- Public names carry the `us` namespace (`us_*` fns; `US_*` constants, enum members and
  error codes).
- Enum values are explicit and append-only.
- Physics functions must document units in their doc comments.
- Never modify `rust-old/` or `lib/`.

## Testing Requirements

- Every public function has golden parity rows against `rust-old/` and the reference
  tests that travel with it.
- Tests go in `tests/*.tcyr`, benchmarks in `*.bcyr`, fuzz targets in `*.fcyr`.
- Benchmarks are recorded against the Rust numbers in
  `docs/benchmarks-rust-v-cyrius.md`.
- Target: coverage ≥ 80 % (`cyrius coverage --exec`), maintained or improved on every PR.

## Commit Messages

Use clear, imperative-mood commit messages:

```
add Van der Waals equation of state

fix entropy calculation for negative heat transfer
```

## License

By contributing, you agree that your contributions will be licensed under the
same license as the project (see [LICENSE](LICENSE)).
