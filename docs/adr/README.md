# Architecture Decision Records

Decisions about ushma — what we chose, the context, and the consequences we accept. Use these when a future reader would reasonably ask *"why did we do it this way?"*

## Conventions

- **Filename**: `NNNN-kebab-case-title.md`, zero-padded to four digits. Never renumber.
- **One decision per ADR.** If a decision supersedes a prior one, add a new ADR and set the old one's status to `Superseded by NNNN`.
- **Status lifecycle**: `Proposed` → `Accepted` → (optionally) `Superseded` or `Deprecated`.
- Use [`template.md`](template.md) as the starting point.

## ADR vs. architecture note vs. guide

| Kind | Lives in | Answers |
|---|---|---|
| ADR | `docs/adr/` | *Why did we choose X over Y?* |
| Architecture note | `docs/architecture/` | *What non-obvious constraint is true about the code?* |
| Guide | `docs/guides/` | *How do I do X?* |

## Index

| ADR | Title | Status |
|---|---|---|
| [0001](0001-cyrius-port-with-rust-oracle.md) | Port everything to Cyrius, with Rust 1.4.0 frozen as the parity oracle | Accepted |
| [0002](0002-port-api-divergences.md) | API divergences from Rust 1.4.0 | Accepted |
| [0003](0003-ganita-gaussian-elimination.md) | `ThermalNetwork::solve` on ganita's Gaussian elimination | Accepted |
| [0004](0004-correct-rounding-is-the-accuracy-target.md) | Correct rounding is the accuracy target | Accepted |
| [0005](0005-serde-compatible-json.md) | serde-compatible JSON codecs on bayan | Accepted |

The Rust 1.x crate had no ADRs.
