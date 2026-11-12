# 0002 — API divergences from Rust 1.4.0

**Status**: Accepted
**Date**: 2026-11-11

## Context

The Cyrius API follows the Rust API, but some Rust shapes have no Cyrius equivalent (String
error payloads, nested `Vec`, async). Others are wrong in the oracle and must not be copied. The
Rust code is frozen (ADR 0001), so every such difference is made in the port only. Each one has
to be named, so the parity harness can still compare the function rather than exempt it.

## Decision

Each row below differs from the oracle on purpose. Its last column says how the golden corpus
still compares it. Rows are appended one per function family as the port finds them; anything
that moves results noticeably goes to the maintainer first, and the decision is quoted here.

| # | Rust 1.4.0 | Cyrius 2.0.0 | Golden comparison |
|---|---|---|---|
| 1 | `UshmaError` variants with `String` and f64 payloads | `US_ERR_*` code enum (append-only) in the value-form `Result`'s one payload word. The detail fields go in an error record (out-pointer or per-call error slot). serde's JSON shape is kept, and `us_error_message` prints thiserror's Display text | code, JSON document and message text: class A |
| 2 | `ThermalGrid2D` stores `Vec<Vec<f64>>`; Crank–Nicolson allocates per step | flat row-major buffer + dims; caller-reusable scratch; the JSON keeps the nested-array shape | state dump flattened the same way, class B; JSON class A |
| 3 | lookups by `&str` name | table index (−1 when absent) plus a name search | same index set, class A |
| 4 | `register_agent` is `async` (tokio + reqwest) | synchronous `sandhi_http_post_opts`, same 30 s timeout. Cyrius async doesn't run on aarch64 | request bytes (method, headers, body) class A; path per #6 |
| 5 | `EnvFilter` directives in `USHMA_LOG` | the level names trace / debug / info / warn / error, and `ushma=<level>`; any other directive is refused with a warning | parsing rows class A; the log line format is sakshi's |
| 6 | `POST /v1/agents/register`, reads `agent_id`. No daimon has ever served that path or sent that field | `POST /v1/agents`, reads `id`, as daimon 2.5.6 serves (prakash made the same correction in 2.5.1) | the oracle's path is mapped to `/v1/agents`; response handling is checked against a stub answering as daimon 2.5.6 |

Exact-constant misses in Rust's formulas (exponents written as rounded decimals, such as
`powf(1.0/3.0)` and `powf(1.0/6.0)`) are fixed under ADR 0004 as each module ports, one row per
family here.

## Consequences

- **Positive** — the port's API is native Cyrius and talks to the real daimon, and no function
  leaves parity: each divergence has a mapping.
- **Negative** — these functions can't reuse the Rust tests one to one; their ledger entries
  say "folded into golden rows with mapping #N".
- **Neutral** — `docs/guides/migration-2.0.md` (M8) lists every row for Rust users.

## Alternatives considered

- **Copy Rust's shapes exactly.** Rejected: nested vecs and String payloads don't exist in
  Cyrius's value-form Results, and the daimon path would never register.
- **Fix the oracle in Rust first.** Rejected by ADR 0001: the Rust code is frozen.
