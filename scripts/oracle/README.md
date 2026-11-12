# scripts/oracle — ushma's golden-value generator

Drives the frozen Rust ushma 1.4.0 in `rust-old/` (by path, read-only; it is never edited)
and writes everything the Cyrius port derives from it:
- the golden corpus, `tests/golden/<module>.tsv`, which the parity harness replays;
- the generated Cyrius sources, `src/constants.cyr` and `src/<module>_data.cyr`.

Adapted from kimiya's `scripts/oracle/` (format and rules unchanged; tanmatra's originally).

```sh
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml             # write
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check  # verify
python3 -I -B scripts/oracle/validate.py .                                          # independent check
```

`--check` regenerates every output in memory and compares it byte for byte. It exits 1 on:
- a difference or a missing file;
- a `.tsv` with no emitter;
- a corpus recorded on another glibc platform.

`Cargo.lock` is committed. It was seeded from `rust-old/Cargo.lock`, so every shared
dependency, hisab 1.4.0 above all, is the version the oracle was tested with.

## What's generated (M0 part 2a)

| Output | Rows / lines | What it is |
|---|---|---|
| `tests/golden/constants.tsv` | 9 | every public scalar const of ushma 1.4.0 (8 `f64`, 1 `usize`), class A |
| `tests/golden/math.tsv` | 15,687 | every std f64 method ushma's library calls, over ushma's ranges (port plan §7.1) |
| `tests/golden/modules.txt` | 2 | the module list the harness reads |
| `src/constants.cyr` | 47 lines | the 9 constants as `US_` consts, under their Rust doc comments |

Added at their milestones (port plan §4 step 6, §6):
- the table data (M0 part 2b): `material_data`, `state_data`, `phase_data`, `chem_data`, and
  `steam_data` (private, extracted with a completeness proof);
- each module's rows (M1 on);
- the JSON rows (each module's milestone);
- the request rows (M7);
- the corrected-oracle files `<module>.cr.tsv` (ADR 0004, M0 part 2c).

## Format 1

One row per call, five tab-separated columns:

```text
fn <TAB> class <TAB> inputs <TAB> outputs <TAB> note
```

- **fn** — the Rust function, e.g. `f64::powf` or `transfer::conduction`; `const` for a constant.
- **class** — the parity class (port plan §5):
  - `A` exact (NaN matches any NaN);
  - `B` arithmetic (≤ 2 ulp, expected 0);
  - `C` transcendental (0 ulp against the corrected oracle, ADR 0004);
  - `D` iterative (relative ≤ 1e-12).
- **inputs**, **outputs** — space-separated typed tokens (`-` when empty):

  | Token | Meaning |
  |---|---|
  | `f:<16 hex>` | an f64 by its IEEE-754 bit pattern, lowercase |
  | `i:<decimal>` | an integer; booleans are 0 / 1 |
  | `e:<Enum>.<Variant>` | a Rust enum variant by name |
  | `s:<percent-encoded UTF-8>` | a string; every byte outside printable ASCII, and `%`, is `%XX` |
  | `n:` | `Option::None` |
  | `E:<Variant>` | `Err(UshmaError::<Variant>)`, the variant name up to its payload |
  | `v:<len>` | a sequence; its elements' tokens follow in order |

- **note** — a readable rendering of the call, for review. The harness ignores it, except for
  the ADR 0002 notes that end `(<tokens>)]` and name the port's value.

Each file starts with five `#` lines: the format and module, the oracle provenance (ushma and
hisab versions; `libm: glibc <version> <arch> <fma|no-fma>`), the regenerate command, the
column legend, and the row count.

## Where the math rows come from

The release-built oracle binary's libm imports (`nm -D target/release/oracle`, rustc 1.99.1,
x86_64): `atan@GLIBC_2.2.5`, `exp@GLIBC_2.29`, `expm1@GLIBC_2.2.5`, `log@GLIBC_2.29`,
`pow@GLIBC_2.29`, `tanh@GLIBC_2.2.5`. `powi` is the statically linked compiler-builtins
`__powidf2`, and `sqrt`, `abs`, `ceil`, `min`, `max`, `clamp`, `is_finite` and `is_nan` are
inline. So the class C rows are glibc 2.44's bits, and `--check` refuses another platform.

## validate.py

It shares no code with the generator. It re-derives from first principles:
- the format of every file, and that `modules.txt` lists exactly the module files;
- every constant, against its Rust literal parsed by Python's correctly rounded `float()`;
- every `src/constants.cyr` literal, against its golden bits;
- every class A/B math row, recomputed in binary64 (`powi` in `__powidf2`'s order, `min`/`max`
  with Rust 1.99's recorded semantics: a NaN operand is ignored, and (±0, ∓0) gives the second
  operand).

Class C rows are checked against MPFR by the corrected-oracle tooling (M0 part 2c). Run it
with `python3 -I -B`. A tampered row fails it, as checked when it was written.
