# Upstream filings from ushma's M0

The port plan's §8, written up in each target repository's own issue format (cyrius's
`docs/development/issues/README.md`, which the fold repositories follow).

**Filed 2026-11-11** into each target's working tree, at the maintainer's word ("File them all"). Nothing was
committed there: the files are new and uncommitted, or appended in place, for the maintainer to commit. This
directory keeps the copies as ushma's record. The filed repros differ from these in two ways. Each names its real
path in its run line. The fold ones also include their repo's own bundle (`include "dist/ganita.cyr"`, …), and the
cyrius one its `lib/` leaves, as the repos' other repros do. Every filed repro was rebuilt and run from its own
repository, where it reproduced, and the runs left those trees unchanged.

Each new issue has a repro that **exits 1 while the defect is present**. Every repro was run on cyrius 6.7.18
(x86_64; the sandhi one also, through probe 7, on PE and aarch64) on 2026-11-11, and reproduced. The evidence for
each is in `docs/development/port-parity.md` §3.

```sh
cyrius build docs/development/upstream/repros/<name>.cyr build/r && ./build/r; echo $?   # .cyr repros, from the ushma root
bash docs/development/upstream/repros/<name>.sh; echo $?                                  # .sh repros (own temp project)
```

## New issues (11)

Each is filed at `<target>/docs/development/issues/<file>`, with its repro at
`<target>/docs/development/issues/repros/2026-11-11-<slug>.*` (this directory's `repros/<target>-<slug>.*`). sandhi's
issue also has a row in its `docs/development/issues/README.md` table of open defects.

| Target | File | Severity | Probe |
|---|---|---|---|
| hisab | [`2026-11-11-deps-sakshi-package-warns-every-consumer.md`](hisab/2026-11-11-deps-sakshi-package-warns-every-consumer.md) | Low | 10 |
| ganita | [`2026-11-11-mat-lu-solve-allocates-per-call.md`](ganita/2026-11-11-mat-lu-solve-allocates-per-call.md) | Medium | 9 |
| ganita | [`2026-11-11-gaussian-elim-pivot-scale-includes-rhs.md`](ganita/2026-11-11-gaussian-elim-pivot-scale-includes-rhs.md) | Medium | 9 |
| sandhi | [`2026-11-11-total-ms-does-not-bound-a-blocked-first-read.md`](sandhi/2026-11-11-total-ms-does-not-bound-a-blocked-first-read.md) | Medium | 7 |
| sakshi | [`2026-11-11-span-enter-has-no-level-check.md`](sakshi/2026-11-11-span-enter-has-no-level-check.md) | Low–Medium (request) | 8 |
| bayan | [`2026-11-11-json-strings-accept-raw-control-bytes-and-invalid-utf8.md`](bayan/2026-11-11-json-strings-accept-raw-control-bytes-and-invalid-utf8.md) | Medium | 11 |
| bayan | [`2026-11-11-no-serde-compatible-shortest-f64-printer.md`](bayan/2026-11-11-no-serde-compatible-shortest-f64-printer.md) | Low–Medium (request) | 11 |
| cyrius | [`2026-11-11-mod-scope-runs-past-the-included-file.md`](cyrius/2026-11-11-mod-scope-runs-past-the-included-file.md) | Medium | 6 |
| cyrius | [`2026-11-11-duplicate-fn-across-bundles-last-wins.md`](cyrius/2026-11-11-duplicate-fn-across-bundles-last-wins.md) | High (silent) | 6 |
| cyrius | [`2026-11-11-derive-writer-faults-on-a-zero-str-field.md`](cyrius/2026-11-11-derive-writer-faults-on-a-zero-str-field.md) | High (SIGSEGV) | 11 |
| cyrius | [`2026-11-11-fdlopen-dynlib-raw-syscalls-warn-on-every-pe-build.md`](cyrius/2026-11-11-fdlopen-dynlib-raw-syscalls-warn-on-every-pe-build.md) | Low | 10 |

## Addenda to open issues (7)

Each `.addendum.md` is the section appended to an issue another port filed, in the shape already used there ("Seen
again", "A third consumer", "Also measured").

| Target issue | Adds |
|---|---|
| ganita [`2026-11-09-correctly-rounded-exp-sin-cos-hypot.md`](ganita/2026-11-09-correctly-rounded-exp-sin-cos-hypot.addendum.md) | a fourth consumer; **tanh** (`ganita_f64_tanh` off on 52 of 860 rows) and **expm1** (no public one); ushma's verified `_us_cr_tanh` to seed the tier |
| cyrius [`2026-11-11-f64-min-max-drop-nan-only-as-first-operand.md`](cyrius/2026-11-11-f64-min-max-drop-nan-only-as-first-operand.addendum.md) | ushma's call sites and its 59-of-600 measurement; a fourth private copy |
| bayan [`2026-11-09-json-duplicate-key-has-no-opt-in-for-ignored-members.md`](bayan/2026-11-09-json-duplicate-key-has-no-opt-in-for-ignored-members.addendum.md) | 46 serde_json documents refused, and the depth cap inside ignored members |
| bayan [`2026-11-09-json-integer-above-i64-max-refuses-the-document.md`](bayan/2026-11-09-json-integer-above-i64-max-refuses-the-document.addendum.md) | both halves of its fix needed: usize at u64::MAX, and an integer beyond i64 for an f64 field |
| cyrius [`2026-11-08-cyrius-port-clobbers-project-docs-and-omits-bench.md`](cyrius/2026-11-08-cyrius-port-clobbers-project-docs-and-omits-bench.addendum.md) | one more reproduction (items 1–4 and pravash's `Cargo.lock` point) |
| bayan [`2026-11-11-json-integer-negative-zero-loses-its-sign.md`](bayan/2026-11-11-json-integer-negative-zero-loses-its-sign.addendum.md) | ushma's 20 documents. pravash filed this issue at 22:06, two minutes before ushma's own filing of the same defect, so ushma's was withdrawn and became this addendum |
| cyrius [`2026-11-11-derive-json-no-vector-array-or-tagged-variant-shapes.md`](cyrius/2026-11-11-derive-json-no-vector-array-or-tagged-variant-shapes.addendum.md) | prefixed variant names, the sequence form, the 31-byte name limit, and the decoder's leniency |

## A note for kimiya

kimiya had no issues directory. The note is filed as `kimiya/docs/development/issues/2026-11-11-notes-from-ushma.md`
(a new directory, in the layout the other repositories use). It says that ganita's LU reproduces Rust hisab 1.4.0's
elimination (its ADR 0003's transcription target), and that `KM_GAS_CONSTANT` and ushma's R differ in the 11th
digit ([copy](kimiya/2026-11-11-notes-from-ushma.md)).

## Not filed, and why

- **cyrius `2026-11-09-field-store-of-struct-call-copies-a-field-argument.md`**: probe 12 re-ran the issue's own
  repro (the same 7 of 12 forms wrong). ushma's grid code is written to the workaround from the start, so it has no
  new shape to add.
- **The `#derive` scanner's three open defects** (2026-11-07): already filed. ushma doesn't use the derive for serde
  JSON (ADR 0005).
- **`large static data (137,984 bytes)`** from sigil's buffers in sandhi's TLS stack: prakash records it.
