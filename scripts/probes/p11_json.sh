#!/usr/bin/env bash
# Probe 11 (port plan §4 step 9.11): bayan's JSON and the f64 printers against serde_json 1.0.154. It writes the
# oracle's documents (scripts/oracle/src/bin/probe_json.rs -> build/probe-json/), copies impetus's src/fmtf.cyr (the
# printer P7 vendors at M1; a sibling checkout at ../impetus, its sha256 printed) to build/p11/, and builds and runs
# scripts/probes/p11_json.cyr. Nothing tracked changes.
#
#   OUT=build/p11 bash scripts/probes/p11_json.sh
set -u
here=$(cd "$(dirname "$0")/../.." && pwd)
cd "$here" || exit 1
fmtf=${FMTF:-$here/../impetus/src/fmtf.cyr}
if [ ! -f "$fmtf" ]; then
    echo "probe 11: no impetus checkout at $fmtf (set FMTF=…/fmtf.cyr)"
    exit 1
fi
mkdir -p build/p11
cp "$fmtf" build/p11/fmtf.cyr
echo "fmtf.cyr: $(sha256sum build/p11/fmtf.cyr | cut -c1-16)… ($(git -C "$(dirname "$fmtf")" log -1 --format=%h -- fmtf.cyr 2>/dev/null))"
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml --bin probe_json 2>&1 | tail -1
cyrius build scripts/probes/p11_json.cyr build/p11/p11 2>&1 \
    | grep -v -E 'large static data|refusing to overwrite stdlib leaf|^warning:lib/|unreachable fns|dead code eliminated' \
    | grep -E 'warning|error|OK'
timeout 600 ./build/p11/p11
echo "exit $?"
