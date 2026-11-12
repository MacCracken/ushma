# Upstream filings from ushma's M0 (drafts)

The port plan's §8, written up in each target repository's own issue format (cyrius's
`docs/development/issues/README.md`, which the fold repositories follow). They are **drafts**: nothing here has
been written into another repository. Each target is the maintainer's to file.

Each new issue has a repro that **exits 1 while the defect is present**. Every repro was run on cyrius 6.7.18
(x86_64; the sandhi one also, through probe 7, on PE and aarch64) on 2026-11-11, and reproduced. The evidence for
each is in `docs/development/port-parity.md` §3.

```sh
cyrius build docs/development/upstream/repros/<name>.cyr build/r && ./build/r; echo $?   # .cyr repros, from the ushma root
bash docs/development/upstream/repros/<name>.sh; echo $?                                  # .sh repros (own temp project)
```

## New issues (12)

When filed, a draft lands at `<target>/docs/development/issues/<file>`, and its repro `repros/<target>-<slug>.*`
lands at `<target>/docs/development/issues/repros/2026-11-11-<slug>.*`.

| Target | File | Severity | Probe |
|---|---|---|---|
| hisab | [`2026-11-11-deps-sakshi-package-warns-every-consumer.md`](hisab/2026-11-11-deps-sakshi-package-warns-every-consumer.md) | Low | 10 |
| ganita | [`2026-11-11-mat-lu-solve-allocates-per-call.md`](ganita/2026-11-11-mat-lu-solve-allocates-per-call.md) | Medium | 9 |
| ganita | [`2026-11-11-gaussian-elim-pivot-scale-includes-rhs.md`](ganita/2026-11-11-gaussian-elim-pivot-scale-includes-rhs.md) | Medium | 9 |
| sandhi | [`2026-11-11-total-ms-does-not-bound-a-blocked-first-read.md`](sandhi/2026-11-11-total-ms-does-not-bound-a-blocked-first-read.md) | Medium | 7 |
| sakshi | [`2026-11-11-span-enter-has-no-level-check.md`](sakshi/2026-11-11-span-enter-has-no-level-check.md) | Low–Medium (request) | 8 |
| bayan | [`2026-11-11-json-integer-minus-zero-loses-its-sign.md`](bayan/2026-11-11-json-integer-minus-zero-loses-its-sign.md) | Low | 11 |
| bayan | [`2026-11-11-json-strings-accept-raw-control-bytes-and-invalid-utf8.md`](bayan/2026-11-11-json-strings-accept-raw-control-bytes-and-invalid-utf8.md) | Medium | 11 |
| bayan | [`2026-11-11-no-serde-compatible-shortest-f64-printer.md`](bayan/2026-11-11-no-serde-compatible-shortest-f64-printer.md) | Low–Medium (request) | 11 |
| cyrius | [`2026-11-11-mod-scope-runs-past-the-included-file.md`](cyrius/2026-11-11-mod-scope-runs-past-the-included-file.md) | Medium | 6 |
| cyrius | [`2026-11-11-duplicate-fn-across-bundles-last-wins.md`](cyrius/2026-11-11-duplicate-fn-across-bundles-last-wins.md) | High (silent) | 6 |
| cyrius | [`2026-11-11-derive-writer-faults-on-a-zero-str-field.md`](cyrius/2026-11-11-derive-writer-faults-on-a-zero-str-field.md) | High (SIGSEGV) | 11 |
| cyrius | [`2026-11-11-fdlopen-dynlib-raw-syscalls-warn-on-every-pe-build.md`](cyrius/2026-11-11-fdlopen-dynlib-raw-syscalls-warn-on-every-pe-build.md) | Low | 10 |

## Addenda to open issues (6)

Each `.addendum.md` is one section to append to an issue another port filed, in the shape already used there
("Seen again", "A third consumer", "Also measured").

| Target issue | Adds |
|---|---|
| ganita [`2026-11-09-correctly-rounded-exp-sin-cos-hypot.md`](ganita/2026-11-09-correctly-rounded-exp-sin-cos-hypot.addendum.md) | a fourth consumer; **tanh** (`ganita_f64_tanh` off on 52 of 860 rows) and **expm1** (no public one); ushma's verified `_us_cr_tanh` to seed the tier |
| cyrius [`2026-11-11-f64-min-max-drop-nan-only-as-first-operand.md`](cyrius/2026-11-11-f64-min-max-drop-nan-only-as-first-operand.addendum.md) | ushma's call sites and its 59-of-600 measurement; a fourth private copy |
| bayan [`2026-11-09-json-duplicate-key-has-no-opt-in-for-ignored-members.md`](bayan/2026-11-09-json-duplicate-key-has-no-opt-in-for-ignored-members.addendum.md) | 46 serde_json documents refused, and the depth cap inside ignored members |
| bayan [`2026-11-09-json-integer-above-i64-max-refuses-the-document.md`](bayan/2026-11-09-json-integer-above-i64-max-refuses-the-document.addendum.md) | both halves of its fix needed: usize at u64::MAX, and an integer beyond i64 for an f64 field |
| cyrius [`2026-11-08-cyrius-port-clobbers-project-docs-and-omits-bench.md`](cyrius/2026-11-08-cyrius-port-clobbers-project-docs-and-omits-bench.addendum.md) | one more reproduction (items 1–4 and pravash's `Cargo.lock` point) |
| cyrius [`2026-11-11-derive-json-no-vector-array-or-tagged-variant-shapes.md`](cyrius/2026-11-11-derive-json-no-vector-array-or-tagged-variant-shapes.addendum.md) | prefixed variant names, the sequence form, the 31-byte name limit, and the decoder's leniency |

## A note for kimiya

kimiya has no issues directory: [`kimiya/2026-11-11-notes-from-ushma.md`](kimiya/2026-11-11-notes-from-ushma.md)
says that ganita's LU reproduces Rust hisab 1.4.0's elimination (its ADR 0003's transcription target), and that
`KM_GAS_CONSTANT` and ushma's R differ in the 11th digit. Where to put it in kimiya is the maintainer's call.

## Not filed, and why

- **cyrius `2026-11-09-field-store-of-struct-call-copies-a-field-argument.md`**: probe 12 re-ran the issue's own
  repro (the same 7 of 12 forms wrong). ushma's grid code is written to the workaround from the start, so it has no
  new shape to add.
- **The `#derive` scanner's three open defects** (2026-11-07): already filed. ushma doesn't use the derive for serde
  JSON (ADR 0005).
- **`large static data (137,984 bytes)`** from sigil's buffers in sandhi's TLS stack: prakash records it.
