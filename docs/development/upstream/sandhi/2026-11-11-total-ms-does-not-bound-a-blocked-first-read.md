# `total_ms` doesn't bound a request whose server accepts and then stays silent; only `read_ms` does — open

**Status:** 🟡 **OPEN**. Verified 2026-11-11 against sandhi 1.13.7 (the cyrius 6.7.18 bundle) on x86_64, x86_64 PE
(wine) and aarch64 (qemu). The repro exits 1.
**Placement:** unpinned.
**Discovered:** 2026-11-11, by ushma's Rust → Cyrius port (M0 probe 7, `ushma/docs/development/port-parity.md` §3).
**Reporter:** ushma 2.0.0 port (`ushma/docs/development/ushma-port-plan.md` §8; ushma ADR 0002 row 7)
**Severity:** Medium. A client that sets only the overall deadline (the natural port of reqwest's `timeout(30 s)`)
hangs for as long as the server stays silent. That is unbounded against a server that accepts and never answers.
**Affects:** sandhi 1.13.7: `sandhi_http_*_opts` with `total_ms` set and `read_ms` unset; `idle_ms` doesn't help.

## Summary

Against a server that accepts, waits 3 s, then answers 200:

| Options | Result |
|---|---|
| `total_ms` 1000 | 200 after 3,001 ms |
| `total_ms` + `idle_ms` 1000 | 200 after 3,000 ms |
| `total_ms` + `read_ms` 1000 | `err_kind` 4 (timeout) after 1,010 ms |

The same happens on all three targets. The documented `total_ms` is a ceiling on the whole request, so the first two
rows should time out at about 1 s.

## Reproduction

`repros/2026-11-11-total-ms-does-not-bound-a-blocked-first-read.cyr`. Its header has a one-line Python server that
accepts, sleeps 3 s and answers. The repro makes the three requests above and exits 1 while either of the first two
runs past 1.5 s.

## Root cause

`_sandhi_http_do_impl_a` turns `total_ms` into a monotonic `deadline_ms`, and clamps the connect timeout to it
(`_sandhi_http_clamp_ms`). The receive side checks the deadline only **between** `recv` calls
(`sandhi_conn_recv_all_deadline`: "each iteration of the read loop checks the deadline before the next recv"). The
receive itself is bounded only by `SO_RCVTIMEO`, which is armed from `read_ms`. With `read_ms` unset, the first `recv`
blocks until the server sends or closes, and the deadline is never consulted.

## Proposed fix

Treat the receive the way connect is already treated. Before each `recv` (or once per exchange, if setting the
socket option per call costs too much), arm `SO_RCVTIMEO` with `min(read_ms, deadline − now)`, using
`_sandhi_http_clamp_ms`. A `recv` that times out under the deadline then reports `err_kind` 4, as the `read_ms` row
does. The same applies to `SO_SNDTIMEO` and the request write.

## Consumer-side workaround

Set `read_ms` as well as `total_ms`. ushma's daimon client sets both to 30 s (its ADR 0002 row 7).
