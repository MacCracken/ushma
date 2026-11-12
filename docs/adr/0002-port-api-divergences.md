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
| 1 | `UshmaError` variants (thiserror) with `String` and f64 payloads, and formatted Display text | `US_ERR_*` codes (append-only) in the value-form `Result`'s payload; `us_error_name(code)` and `us_error_message(code)` give each code a fixed message. The failing values go to a level-gated sakshi event at the failure site, not into the error. Rust's formatted message strings are Rust machinery, not a capability, so they are not reproduced (the port is Cyrius-native: Rust is the oracle for results) | the error **code**, class A; message text is ushma's own and not compared |
| 2 | `ThermalGrid2D` stores `Vec<Vec<f64>>`; Crank–Nicolson allocates per step | flat row-major buffer + dims; caller-reusable scratch; the JSON keeps the nested-array shape | state dump flattened the same way, class B; JSON class A |
| 3 | lookups by `&str` name | table index (−1 when absent) plus a name search | same index set, class A |
| 4 | `register_agent` is `async` (tokio + reqwest) | synchronous `sandhi_http_post_opts`, same 30 s timeout. Cyrius async doesn't run on aarch64 | request bytes (method, headers, body) class A; path per #6 |
| 5 | `EnvFilter` directives in `USHMA_LOG` | the level names trace / debug / info / warn / error, and `ushma=<level>`; any other directive is refused with a warning | parsing rows class A; the log line format is sakshi's |
| 6 | `POST /v1/agents/register`, reads `agent_id`. No daimon has ever served that path or sent that field | `POST /v1/agents`, reads `id`, as daimon 2.5.6 serves (prakash made the same correction in 2.5.1) | the oracle's path is mapped to `/v1/agents`; response handling is checked against a stub answering as daimon 2.5.6 |
| 7 | reqwest's `timeout(30 s)` bounds the whole request | sandhi's `total_ms` alone doesn't bound a server that is silent before its first byte (probe 7), so the client sets `total_ms` **and** `read_ms`, 30 s each | stub rows: a silent server gives the timeout error |
| 8 | reqwest sends exactly host, content-type, content-length, authorization and `accept: */*` | sandhi adds `User-Agent: sandhi/<v>`, `Accept-Encoding: identity`, `TE: trailers` and `Connection: close`. The port sets `Accept: */*` and writes the body in serde's key order (`capabilities` before `name`) | body bytes class A; the oracle's five headers by value, names case-insensitive; sandhi's four transport headers listed as extra |
| 9 | correlation exponents written as Rust doubles: Dittus–Boelter's `powf(0.8)`, `powf(0.4)`; Churchill–Chu's `powf(1.0 / 6.0)`, `powf(9.0 / 16.0)`, `powf(-8.0 / 27.0)`; Pohlhausen's `powf(1.0 / 3.0)` | the same doubles, through the correctly rounded `_us_cr_pow`. Each is a fitted or approximate correlation's exponent, accurate to a few percent, so restoring the "exact" 1/6 or 1/3 would change results by under 1 ulp of a quantity known to ~10 %; the double Rust wrote is the definition. kimiya's #16 differs on purpose: there 2^(1/6) is a physical constant (the Lennard-Jones minimum) | class C against the corrected oracle: no mapped rows |
| 10 | `otto_cycle`, `diesel_cycle` and `brayton_cycle` return a `CycleResult` by value that owns a `Vec<StatePoint>` and a `Vec<ProcessKind>`, built fresh on every analysis; a refused allocation aborts the process | the caller owns a `UsCycleResult` (`us_cycle_result_init`, or `_init_a` on an allocator), which each analysis fills through `r: *UsCycleResult` and returns Ok(0). Its state points (five f64 words each, flat) and processes are two stdlib vecs that grow on the result's own allocator. An analysis sizes them and stores in place, so a result reused across calls allocates nothing after its first (1,000 analyses, 0 heap bytes). `us_cycle_state` / `us_cycle_process` read it, `us_cycle_push_state` / `us_cycle_push_process` build one by hand, and the decoder builds one on its arena. A refused allocation is `Err(US_ERR_OUT_OF_MEMORY)`, a code appended after row 1's, never an abort | each number, each state point's five words and each process code: class B, class C where a row rounds an ln or pow; JSON class A |
| 11 | `rankine_cycle(p_condenser, p_boiler, t_superheat: Option<f64>)`; `refrigeration_cycle` returns a `RefrigerationResult { cycle: CycleResult, cop_refrigeration, cop_heat_pump }` by value | the `Option` is a flag and a value, the form impetus's ports use: `us_rankine_cycle(p_condenser, p_boiler, superheated, t_superheat, r)`, where `superheated` 0 is `None` and `t_superheat` is then unread. A sentinel can't stand for `None`, since every f64 (NaN included) is a `Some` the lookup must refuse. `UsRefrigerationResult` nests its `UsCycleResult` as the first field, and `us_refrigeration_cycle` fills it through `out` as row 10's cycles do. `&rr.cycle` is no Cyrius expression, so `us_refrigeration_result_cycle(rr)` gives the cycle's address for the row 10 accessors (the struct's own address) | `n:` (None) passes 0 and `f:` passes 1 with the value: class C; the result's tokens are the cycle's then the two COPs; JSON class A. No mapped rows |
| 12 | `cycle_ts_diagram` / `cycle_pv_diagram` return a fresh `Vec<DiagramPoint>`, preallocating up to `MAX_DIAGRAM_PREALLOC` (2¹⁶); `compare_cycles(&[(CycleKind, &CycleResult, f64, f64)])` returns a fresh `Vec<CycleComparisonEntry>` | the diagrams write `UsDiagramPoint` records into the caller's buffer of `cap` records and return Ok(the count), or Err(`US_ERR_INVALID_PARAMETER`) when they would not fit (probe 4: 2.0 µs and no heap, against 7 µs and about 8 KB leaked for a fresh list). `us_cycle_diagram_len` sizes the buffer; the check divides, so a huge point count cannot overflow. With no allocation there is nothing to preallocate, and the private constant has no role (it is still generated, `_US_CYCLE_MAX_DIAGRAM_PREALLOC`). `us_compare_cycles(kinds, results, t_highs, t_lows, n, out)` takes the tuples as four parallel lists, as Kay's rules take theirs, and writes n `UsCycleComparison` records into the caller's buffer; `us_compare_cycle` is one entry | a diagram row's `v:count` and every point, class C; a comparison row's `v:m` and every entry, class B; JSON class A. A buffer too small has no oracle row: it is the Cyrius shape's own refusal, tested in `tests/cycle.tcyr`. No mapped rows |
| 13 | chem's reaction sides are `&[(f64, &Species)]`: `reaction_enthalpy(products, reactants)`, `reaction_gibbs(products, reactants)`, `adiabatic_flame_temperature(reactants, products, t_initial)` | a side is a term list: n consecutive word pairs (the coefficient, an f64; the species, a `*UsSpecies`), the slice of pairs laid out flat, with its length beside it: `us_reaction_enthalpy(products, np, reactants, nr)`, `us_adiabatic_flame_temperature(reactants, nr, products, np, t_initial)`. Rust's pair belongs together, so it stays one record, where row 12's four-field tuple became parallel lists. Rust's sides keep their order | each term written as its coefficient and its species by value, so edge species sit beside the presets: class B, class C where exp enters. No mapped rows |

Exact-constant misses in Rust's formulas are reviewed as each module ports, one row per family here. An
exponent that is a physical constant is restored to its exact value; a correlation's fitted exponent is kept as the
double Rust wrote (row 9).

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
