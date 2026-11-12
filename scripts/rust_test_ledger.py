#!/usr/bin/env python3
"""The Rust-test ledger (port plan §6 and M8): every test of the Rust oracle is ported, with nothing deferred.

    python3 -I -B scripts/rust_test_ledger.py .           # write docs/development/rust-test-ledger.md
    python3 -I -B scripts/rust_test_ledger.py . --check   # exit 1 if the ledger is stale or a due test is unported

It lists each `#[test]` in rust-old/src/**.rs and rust-old/tests/integration.rs. A Cyrius test function names the
Rust tests it ports in its doc comment, one or more lines of

    # rust <file under rust-old/>: <test>, <test>, ...

and a Rust test whose assertions are spread over several groups is named by each. A test is due when its module's
milestone has come (PORTED); an integration test when every module it calls is PORTED. A due test with no port is
an error, as is a marker naming no Rust test. Tests not yet due are counted under their milestone.
"""
import re
import sys
from pathlib import Path

root = Path(sys.argv[1]).resolve()
check = "--check" in sys.argv[2:]
out_path = root / "docs/development/rust-test-ledger.md"

PORTED = ["error", "validate", "logging", "material", "transfer", "state", "entropy", "phase", "steam", "cycle", "chem", "stat"]
MILESTONE = {
    "error": "M1", "validate": "M1", "logging": "M1", "material": "M1", "transfer": "M1",
    "state": "M2", "entropy": "M2", "phase": "M3", "steam": "M3", "cycle": "M4", "chem": "M5", "stat": "M5",
    "numerical": "M6", "bridge": "M7", "integration/soorat": "M7", "ai": "M7",
}
MODULE_RE = re.compile(r"\b(" + "|".join(m for m in MILESTONE if "/" not in m) + r")::")


def rust_tests(rel):
    """(name, body) of each #[test] in rust-old/<rel>."""
    src = (root / "rust-old" / rel).read_text()
    out = []
    for m in re.finditer(r"#\[test\]\s*\n\s*fn (\w+)\(\)[^{]*\{", src):
        depth, i = 1, m.end()
        while depth and i < len(src):
            depth += {"{": 1, "}": -1}.get(src[i], 0)
            i += 1
        out.append((m.group(1), src[m.end():i]))
    return out


def markers():
    """(rust file, test) -> [(cyrius file, fn)], and the markers' own problems."""
    found, problems = {}, []
    for path in sorted((root / "tests").glob("*.tcyr")):
        pending = []
        for n, line in enumerate(path.read_text().splitlines(), 1):
            m = re.match(r"# rust (\S+): (.+)$", line)
            if m:
                pending += [(m.group(1), t.strip(), n) for t in m.group(2).split(",") if t.strip()]
                continue
            fm = re.match(r"fn (\w+)\(", line)
            if fm and pending:
                for rel, t, at in pending:
                    found.setdefault((rel, t), []).append((path.name, fm.group(1)))
                pending = []
            elif not line.startswith("#") and pending:
                problems.append(f"tests/{path.name}:{pending[0][2]}: a `# rust` marker not followed by its fn")
                pending = []
    return found, problems


ported, errors = markers()
files = sorted(str(p.relative_to(root / "rust-old")) for p in (root / "rust-old/src").rglob("*.rs"))
files.append("tests/integration.rs")
rows, all_tests = [], set()
for rel in files:
    tests = rust_tests(rel)
    if not tests:
        continue
    mod = rel[len("src/"):-3] if rel.startswith("src/") else None
    counts = {"ported": 0, "due": 0, "later": {}}
    for name, body in tests:
        all_tests.add((rel, name))
        if mod is not None:
            due, ms = mod in PORTED, MILESTONE[mod]
        else:
            used = set(MODULE_RE.findall(body))
            due = used <= set(PORTED)
            ms = max((MILESTONE[u] for u in used), default="M1")
        if (rel, name) in ported:
            counts["ported"] += 1
        elif due:
            counts["due"] += 1
            errors.append(f"rust-old/{rel}: {name} is due ({ms}) and has no `# rust` marker")
        else:
            counts["later"][ms] = counts["later"].get(ms, 0) + 1
    rows.append((rel, len(tests), counts))
for (rel, name), where in sorted(ported.items()):
    if (rel, name) not in all_tests:
        errors.append(f"tests/{where[0][0]} {where[0][1]}: names rust-old/{rel} {name}, which is no Rust test")

total = sum(r[1] for r in rows)
n_ported = sum(r[2]["ported"] for r in rows)
later = {}
for _, _, c in rows:
    for ms, k in c["later"].items():
        later[ms] = later.get(ms, 0) + k
md = [
    "# Rust-test ledger",
    "",
    "> **Generated** by `python3 -I -B scripts/rust_test_ledger.py .`; `--check` gates it (port plan §6). Every test",
    "> of the frozen Rust oracle (`rust-old/`, ushma 1.4.0) ports to a Cyrius test whose doc comment names it",
    "> (`# rust <file>: <test>`), at its module's milestone. Nothing is deferred: a test not ported yet is one whose",
    "> milestone hasn't come. The golden corpus (`tests/golden/`) checks every ported function on its grid besides.",
    "",
    f"**{total} Rust tests:** {n_ported} ported; "
    + (", ".join(f"{k} at {ms}" for ms, k in sorted(later.items())) or "none left") + ".",
    "",
    "| Rust file | Tests | Ported | Due, unported | Later (milestone) |",
    "|---|---:|---:|---:|---|",
]
for rel, n, c in rows:
    lt = ", ".join(f"{k} ({ms})" for ms, k in sorted(c["later"].items())) or "—"
    md.append(f"| `{rel}` | {n} | {c['ported']} | {c['due']} | {lt} |")
n_due = sum(r[2]["due"] for r in rows)
md += [f"| **Total** | **{total}** | **{n_ported}** | **{n_due}** | |", "", "## Where each is ported", ""]
md += ["| Rust test | Cyrius test |", "|---|---|"]
for (rel, name), where in sorted(ported.items()):
    if (rel, name) in all_tests:
        md.append(f"| `{rel}` `{name}` | " + ", ".join(f"`tests/{f}` `{fn}`" for f, fn in where) + " |")
md.append("")
text = "\n".join(md)

for e in errors:
    print("LEDGER  ", e)
if check:
    old = out_path.read_text() if out_path.exists() else ""
    if old != text:
        print(f"STALE    {out_path.relative_to(root)} (regenerate: python3 -I -B scripts/rust_test_ledger.py .)")
        errors.append("stale")
else:
    out_path.write_text(text)
print(f"rust-test ledger: {total} Rust tests, {n_ported} ported, {sum(later.values())} later; {len(errors)} problem(s)")
sys.exit(1 if errors else 0)
