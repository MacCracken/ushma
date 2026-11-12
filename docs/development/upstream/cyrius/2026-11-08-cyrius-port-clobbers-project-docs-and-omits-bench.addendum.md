<!-- Draft addendum: appended to cyrius/docs/development/issues/2026-11-08-cyrius-port-clobbers-project-docs-and-omits-bench.md -->

## Seen again: ushma (2026-11-11, cyrius 6.7.18)

ushma's port (a Rust library with CLAUDE.md, roadmap.md, state.md and `VERSION` 1.4.0) hit items 1–4 again:
- the documents and both workflows were overwritten;
- the manifest got `version = "0.1.0"` and `[test] stdlib = ["assert"]`;
- the scaffold's stubs failed `cyrius doc --check` until they were documented.

It also hit pravash's `Cargo.lock` point: the library's `.gitignore` still ignored `Cargo.lock` after the move, so
`rust-old/Cargo.lock`, which the parity oracle must rebuild against, was untracked until the `.gitignore` was
rewritten. Nothing new beyond pravash's list; one more reproduction.
