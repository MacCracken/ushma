#!/usr/bin/env bash
# ushma — the full gate run (port plan §6): oracle, three targets, bundles, hygiene.
# Adapted from kimiya's handoff/gates/gates.sh. Steps whose tooling doesn't exist yet (the
# oracle generator lands in M0 part 2, [lib] in M1, the parity scripts later) are reported
# as SKIP, never silently passed.
#
#   OUT=/tmp/ushma-gates bash scripts/gates.sh > /tmp/ushma-gates/gates.log 2>&1
#
# Rules: run sequentially; never run two wine processes at once; a truncated wine run can
# exit 0, so check that each test log carries its summary line.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${OUT:-/tmp/ushma-gates}"
mkdir -p "$OUT"
cd "$ROOT" || exit 2

# Diagnostics that come from vendored lib/ and are recorded in docs/development/state.md.
# Anything else is ushma's and fails the warning count.
KNOWN='large static data|refusing to overwrite stdlib leaf .sakshi.|^warning:lib/'

# Every step that exits nonzero is named at the end, and the run then exits 1.
FAILED=()
run() { # PAT=<grep -E> [TAILN=n] run <cmd...>
    echo "=== $*"
    "$@" > "$OUT/step.out" 2>&1; rc=$?
    grep -v '^note:' "$OUT/step.out" | grep -iE "${PAT:-.}" | tail -"${TAILN:-4}"
    echo "--- exit $rc"
    [ "$rc" -eq 0 ] || FAILED+=("$* (exit $rc)")
}
skip() { echo "=== SKIP: $*"; }

# 1. The Rust oracle (read-only) and the golden-value generator.
PAT='test result' TAILN=6 run cargo test --manifest-path rust-old/Cargo.toml --all-features
if [ -f scripts/oracle/Cargo.toml ]; then
    PAT='ok|fail|error|warn|pass|check' run cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check
    PAT='warn|error|Finished' run cargo clippy --release --locked --manifest-path scripts/oracle/Cargo.toml --all-targets -- -D warnings
    PAT='self-checks|Error|assert' run python3 -I -B scripts/oracle/mpfr_ref.py
    for py in validate.py probe_tables.py cr_verify.py mathx_tables.py; do
        if [ -f "scripts/oracle/$py" ]; then
            case "$py" in
                mathx_tables.py)
                    flag=--check
                    [ -f ../kimiya/src/mathx_data.cyr ] && flag=--check-kimiya
                    PAT='STALE|DIFFERS|current' run python3 -I -B "scripts/oracle/$py" "$flag" ;;
                probe_tables.py) PAT='FAIL|checked' run python3 -I -B "scripts/oracle/$py" . --all-targets ;;
                cr_verify.py) PAT='MISMATCH|DIFFER|cr_verify:|special values' TAILN=6 run python3 -I -B "scripts/oracle/$py" ;;
                *) PAT='FAIL|rows|ok|checked|mismatch' TAILN=6 run python3 -I -B "scripts/oracle/$py" . ;;
            esac
        fi
    done
else
    skip "scripts/oracle/ (M0 part 2)"
fi

# 2. Build, bundles.
PAT='OK|error|warn|FAIL' run cyrius build src/main.cyr build/ushma
if grep -q '^\[lib\]' cyrius.cyml; then
    PAT='OK|error|warn|lines' run cyrius distlib --all
    PAT='.' run cyrius distlib --check
else
    skip "cyrius distlib (no [lib] until M1)"
fi

# 3. Tests on every target: ELF, ELF --bounds, PE (wine), aarch64 (qemu).
for t in '' --bounds --win --aarch64; do
    log="$OUT/test${t:-_elf}.log"
    echo "=== cyrius test $t"
    cyrius test $t > "$log" 2>&1
    rc=$?
    echo "--- exit $rc"
    grep -E 'parity summary' "$log"
    grep -E '^[0-9]+ passed, [0-9]+ failed' "$log" | tail -1
    warns=$(grep -i warning "$log" | grep -viE "$KNOWN" | wc -l)
    echo "ushma warning lines: $warns"
    [ "$rc" -eq 0 ] || FAILED+=("cyrius test $t (exit $rc)")
    grep -qE 'parity summary: .* 0 failing' "$log" || FAILED+=("cyrius test $t (no clean parity summary)")
    grep -qE '^[1-9][0-9]* passed, 0 failed' "$log" || FAILED+=("cyrius test $t (no clean assert summary)")
    [ "$warns" -eq 0 ] || FAILED+=("cyrius test $t ($warns ushma warning lines)")
done

# 4. Benchmarks, hygiene, lock, audit.
PAT='passed|failed|ns' TAILN=18 run cyrius bench
echo "=== fmt/lint/doc per file"
fails=0
for f in src/*.cyr tests/*.tcyr tests/*.bcyr tests/*.fcyr tests/*.cyr; do
    [ -f "$f" ] || continue
    cyrius fmt --check "$f" > /dev/null 2>&1 || { echo "fmt FAIL $f"; fails=$((fails + 1)); }
    cyrius lint "$f" 2>&1 | grep -E '^[1-9][0-9]* warnings?$|: warning' && { echo "  ^ lint $f"; fails=$((fails + 1)); }
done
for f in src/*.cyr; do
    cyrius doc --check "$f" 2>&1 | grep -E '[1-9][0-9]* undocumented' && { echo "  ^ doc $f"; fails=$((fails + 1)); }
done
echo "--- fmt/lint/doc failures $fails"
[ "$fails" -eq 0 ] || FAILED+=("fmt/lint/doc ($fails)")
PAT='.' run cyrius vet src/main.cyr
PAT='.' TAILN=2 run cyrius deps --locked
PAT='.' TAILN=14 run cyrius audit

# 5. Parity bookkeeping (from M1).
for py in surface_parity.py rust_test_ledger.py; do
    if [ -f "scripts/$py" ]; then PAT='parity|missing|ported|FAIL' run python3 -I -B "scripts/$py" . ; else skip "scripts/$py (from M1)"; fi
done

echo "=== staged (must be empty: the maintainer does all git)"
git diff --cached --name-only
[ -z "$(git diff --cached --name-only)" ] || FAILED+=("staged files (the maintainer does all git)")

if [ "${#FAILED[@]}" -gt 0 ]; then
    echo "=== FAILED: ${#FAILED[@]} step(s)"
    printf '  %s\n' "${FAILED[@]}"
    exit 1
fi
echo "=== all steps exited 0"
