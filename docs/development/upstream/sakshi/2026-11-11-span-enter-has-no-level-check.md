# `sakshi_span_enter` / `_exit` take no level and never test one, so a filtered-out debug span costs a clock read and an emit on every call — open

**Status:** 🟡 **OPEN**. Verified 2026-11-11 against sakshi 2.5.10 (the cyrius 6.7.18 stdlib leaf). The repro exits 1.
This is a request, not a defect.
**Placement:** unpinned.
**Discovered:** 2026-11-11, by ushma's Rust → Cyrius port (M0 probe 8, `ushma/docs/development/port-parity.md` §3).
**Reporter:** ushma 2.0.0 port (`ushma/docs/development/ushma-port-plan.md` §8)
**Severity:** Low–Medium. Performance only, and only for libraries that instrument hot functions. Each such library
writes the same level-gated wrapper today.
**Affects:** sakshi 2.5.10, `sakshi_span_enter` (`lib/sakshi.cyr:2037`) and `sakshi_span_exit` (:2070).

## Summary

Rust libraries instrument their public functions with `#[tracing::instrument(level = "debug")]`. Under the default
INFO filter that costs almost nothing: Rust 1.4.0's `ideal_gas_volume` takes 3.7 ns with its debug span. sakshi's span
API has no level:
- `sakshi_span_enter` reads the clock, pushes the span stack and emits an ENTER event on every call, whatever
  `sakshi_set_level` says;
- `sakshi_span_exit` emits the EXIT.

Under the default output, that is a formatted line to stderr per span: an INFO-level program prints every span.

M0 probe 8 measured n·R·T/P (8.1 ns plain):

| Variant | ns per call |
|---|---|
| no span | 8.1 |
| `if (sakshi_get_level() >= SK_DEBUG)` around enter and exit, DEBUG filtered out | 9.7 |
| sakshi span into its ring buffer (`sakshi_output_buffer`, its cheapest output) | 192 |

So every port wraps the span calls in its own level test (ushma: `_us_span_enter` / `_us_span_exit`).

## Reproduction

`repros/2026-11-11-span-enter-has-no-level-check.cyr` sets the level to `SK_ERROR` and the output to the ring buffer,
then opens and closes 10 spans. It prints `events in the ring after 10 spans at level SK_ERROR: 20` and exits 1.

## Root cause

By design: the span functions predate a level parameter, and `_sk_emit_span` doesn't consult `_sk_log_level`.

## Proposed fix

Additive:
- `sakshi_span_enter_at(level, name, name_len)` and a matching `sakshi_span_exit_at(level)`. Each returns before the
  clock read when `level > _sk_log_level`, or when `level > SAKSHI_LEVEL` at compile time, which makes it free.
- Keep the stack balanced: an `_at` exit skips exactly what its `_at` enter skipped.
- The existing `sakshi_span_enter` / `_exit` keep their meaning, as level-less (always-on) spans.

## Consumer-side workaround

The level-gated wrapper above: about 1.6 ns when DEBUG is off.
