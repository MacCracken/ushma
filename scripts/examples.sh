#!/usr/bin/env bash
# The five examples (port plan M8): each Cyrius example in docs/examples/ is built on the published bundle
# (dist/ushma.cyr, as a consumer builds it) and run, and its output must be docs/examples/<name>.expected byte for
# byte. That file is Rust 1.4.0's example's output (rust-old/examples/<name>.rs).
#
#   bash scripts/examples.sh                  # x86_64 ELF
#   bash scripts/examples.sh --all-targets    # ELF, PE (wine, one at a time) and aarch64 (qemu)
#   bash scripts/examples.sh --rust           # also rerun the Rust examples: .expected must be their output
#   bash scripts/examples.sh --write          # rewrite .expected from the Rust examples (then check the Cyrius ones)
#
# The Rust examples build into rust-old/target (ignored); rust-old/ is never edited. Exits 1 on any difference.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT" || exit 2
EXAMPLES="basic_thermo heat_engine building_thermal cooling_system phase_change"
rust=0
write=0
targets="elf"
for a in "$@"; do
    case "$a" in
        --rust) rust=1 ;;
        --write) rust=1; write=1 ;;
        --all-targets) targets="elf win aarch64" ;;
        *) echo "examples.sh: unknown argument $a"; exit 2 ;;
    esac
done
mkdir -p build/examples
fails=0
[ -f dist/ushma.cyr ] || { echo "examples.sh: dist/ushma.cyr is missing (cyrius distlib --all)"; exit 1; }
for e in $EXAMPLES; do
    want="docs/examples/$e.expected"
    if [ "$rust" -eq 1 ]; then
        got="build/examples/$e.rust.out"
        if ! cargo run --release --locked --quiet --manifest-path rust-old/Cargo.toml --all-features --example "$e" \
            > "$got" 2> build/examples/$e.rust.err; then
            echo "FAIL rust $e: cargo run failed (build/examples/$e.rust.err)"
            fails=$((fails + 1))
            continue
        fi
        if [ "$write" -eq 1 ]; then
            cp "$got" "$want"
            echo "wrote $want ($(wc -l < "$want") lines)"
        elif cmp -s "$got" "$want"; then
            echo "ok   rust $e: its output is $want"
        else
            echo "FAIL rust $e: its output differs from $want"
            diff "$want" "$got" | head -10
            fails=$((fails + 1))
        fi
    fi
    for t in $targets; do
        case "$t" in
            elf) flag=""; bin="build/examples/$e"; runner="" ;;
            win) flag="--win"; bin="build/examples/$e.exe"; runner="wine" ;;
            aarch64) flag="--aarch64"; bin="build/examples/$e.arm"; runner="qemu-aarch64" ;;
        esac
        if ! cyrius build $flag "docs/examples/$e.cyr" "$bin" > build/examples/$e.$t.build 2>&1; then
            echo "FAIL $t $e: build failed (build/examples/$e.$t.build)"
            fails=$((fails + 1))
            continue
        fi
        got="build/examples/$e.$t.out"
        if ! timeout 300 $runner "$bin" > "$got" 2> build/examples/$e.$t.err; then
            echo "FAIL $t $e: exited nonzero"
            fails=$((fails + 1))
            continue
        fi
        if cmp -s "$got" "$want"; then
            echo "ok   $t $e: $(wc -l < "$got") lines, Rust's output"
        else
            echo "FAIL $t $e: output differs from $want"
            diff "$want" "$got" | head -10
            fails=$((fails + 1))
        fi
    done
done
echo "examples: $fails failure(s)"
[ "$fails" -eq 0 ]
