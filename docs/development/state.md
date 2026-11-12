# ushma — Current State

> Refreshed every milestone. CLAUDE.md is preferences/process/procedures
> (durable); this file is **state** (volatile).

## Version

**1.4.0 → 2.0.0-dev.** `VERSION` stays `1.4.0` until M8 ships 2.0.0. The Rust crate's final
release, **1.4.0** (tag `1.4.0`, commit `2688099`), is frozen at `rust-old/` as the parity oracle.
It is never edited and never re-released.

## Port progress

Plan of record: [`ushma-port-plan.md`](ushma-port-plan.md). Decisions P2–P8 were accepted
2026-11-11 and are recorded in ADRs 0001–0005.

| Milestone | Status |
|---|---|
| M0 Prep & initialize | 🟡 in progress. **Part 1 done (2026-11-11):** oracle verified (Rust code identical to `9528e1f`), Rust bench baseline on the tag, `cyrius port`, repairs, manifest, ADRs 0001–0005. **Next, part 2:** golden-value generator (`scripts/oracle/`), tables, parity harness, `src/mathx.cyr`. **Part 3:** baselines on three targets, probes 1–12, upstream filings |
| M1–M8 | ☐ |

## Toolchain

- **Cyrius pin**: `6.7.18` (in `cyrius.cyml [package].cyrius`)
- **Rust oracle**: rustc 1.99.1, glibc 2.44 (the oracle's transcendental bits depend on it)

## Source

- **Rust reference**: `rust-old/`, 11,046 lines in total (`LINES_OF_SOURCE.txt`); ≈5,400 code
  lines in `src/`, 155 public functions, 449 unit + 30 integration tests, 63 benchmarks,
  5 examples. Frozen.
- **Cyrius port**: scaffold only (`src/main.cyr` smoke entry, `tests/ushma.{tcyr,bcyr,fcyr}`
  stubs). No `[lib]` yet; it starts at M1.

## Rust baseline (tag 1.4.0, 2026-11-11)

- **Gates:** 449 unit + 30 integration tests pass; fmt, clippy `-D warnings`, rustdoc
  `-D warnings`, `cargo audit` and `cargo deny` are clean.
- **Benchmarks:** 63, recorded in `rust-old/bench-history.csv` (commit `2688099`, run
  2026-11-12T02:59Z) and `rust-old/benchmarks.md`. These numbers are the Rust column of
  `docs/benchmarks-rust-v-cyrius.md`.

## Gates (scaffold, x86_64 ELF)

`cyrius build`, `cyrius test` and `test --bounds` (2/2), `cyrius bench`, `fmt --check`, `lint`,
`doc --check` (every stub documented), `vet`, `deps --locked` and `cyrius audit` all pass.

**Known upstream diagnostics** (from vendored `lib/`, never edited):
- `warning: large static data (137,984 bytes)` on every compile. This is sigil's crypto buffers,
  pulled in by sandhi → tls (`[lib.ai]`'s stack). prakash records the same diagnostic.
- `cyrius deps`: `refusing to overwrite stdlib leaf 'sakshi'`. hisab 3.4.1 declares an unused
  `[deps.sakshi]`; the stdlib copy wins, and the two are byte-identical (port plan §2.7, §8).
- PE only (`--win`): six `syscall N … is not routed on CYRIUS_TARGET_WIN=1; it returns -38` warnings
  in `lib/fdlopen.cyr` and `lib/dynlib.cyr` (sandhi's companions). Probe 7 checks whether
  `[lib.ai]`'s TLS path on PE reaches them.

`scripts/gates.sh` filters these, and any diagnostic located in `lib/`. ushma's own code must
add no diagnostic.

## Dependencies

Direct (declared in `cyrius.cyml`):
- stdlib: syscalls, string, alloc, str, fmt, io, vec, slice, math, ganita, tagged, result,
  boxed, overflow, fnptr, freelist, sakshi, log, bayan, sandhi, tls, async, random, fdlopen,
  dynlib;
- `[deps.hisab]` 3.4.1 (`dist/hisab.cyr`);
- test stdlib: assert, bench.

`cyrius.lock` pins 80 deps.

## Consumers

Rust 1.x: kiran, joshua and soorat (its `thermo` feature, `ushma = "1"` with `soorat-compat`).
No Cyrius consumer yet. crates.io reverse dependencies are not yet checked; the check queries
crates.io for "ushma" and waits on the maintainer.

## Next

M0 part 2 (port plan §4 steps 4–7).
