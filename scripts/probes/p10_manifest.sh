#!/usr/bin/env bash
# Probe 10 (port plan §4 step 9.10): the manifest, measured in a scratch copy so nothing tracked changes.
# 1. sakshi twice: `cyrius deps` with hisab 3.4.1's [deps.sakshi] package beside the sakshi stdlib leaf (log, sandhi),
#    and the byte identity of the two copies at the pin.
# 2. [lib] and [lib.ai]: the real src/ modules plus a stand-in that reaches every production dependency
#    (p10_manifest/core.cyr) and the sandhi client (p10_manifest/ai.cyr); `cyrius distlib --all`, the .deps
#    sidecars, the requires blocks and `distlib --check`.
# 3. Consumers of each bundle, with the sidecar's leaves as their stdlib and with none: built for x86_64 ELF, PE and
#    aarch64 and run (wine and qemu), the ai one against the probe-7 stub. One wine process at a time.
#
#   OUT=build/p10 bash scripts/probes/p10_manifest.sh
set -u
here=$(cd "$(dirname "$0")/../.." && pwd)
OUT=${OUT:-$here/build/p10}
case "$OUT" in /*) ;; *) OUT=$here/$OUT ;; esac
rm -rf "$OUT"
mkdir -p "$OUT"
S=$OUT/ushma
mkdir -p "$S"
cp -r "$here/cyrius.cyml" "$here/cyrius.lock" "$here/VERSION" "$here/src" "$here/lib" "$S/"
cp "$here/scripts/probes/p10_manifest/core.cyr" "$S/src/p10_core.cyr"
cp "$here/scripts/probes/p10_manifest/ai.cyr" "$S/src/p10_ai.cyr"
core='"src/constants.cyr", "src/constants_private.cyr", "src/mathx_data.cyr", "src/mathx.cyr", "src/material_data.cyr",
    "src/state_data.cyr", "src/phase_data.cyr", "src/chem_data.cyr", "src/steam_data.cyr", "src/p10_core.cyr"'
cat >> "$S/cyrius.cyml" <<EOF

[lib]
modules = [
    $core,
]

[lib.ai]
modules = [
    $core,
    "src/p10_ai.cyr",
]
EOF

say() { printf '\n=== %s\n' "$*"; }

say "1. cyrius deps (sakshi: hisab's package and the stdlib leaf)"
(cd "$S" && cyrius deps 2>&1 | grep -v '^\s*$' | grep -i -E 'sakshi|warn|error|refus|hisab' | sort -u)
say "sakshi copies (sha256)"
sha256sum "$S/lib/sakshi.cyr" "$HOME/.cyrius/versions/6.7.18/lib/sakshi.cyr" "$HOME/.cyrius/deps/sakshi/2.5.10/dist/sakshi.cyr" 2>&1 \
    | sed "s|$HOME|~|"

say "2. cyrius distlib --all"
(cd "$S" && cyrius distlib --all 2>&1 | grep -v '^\s*$' | sed "s|$S/||g")
say "dist/"
(cd "$S/dist" && wc -c ./* | sed 's|\./||')
for d in "$S"/dist/*.deps; do
    say "$(basename "$d")"
    grep -v '^#' "$d" | tr '\n' ' '
    echo
done
for b in "$S"/dist/ushma.cyr "$S"/dist/ushma-ai.cyr; do
    say "$(basename "$b"): requires block, and a hisab include?"
    grep -n '^include\|^# Requires' "$b" | head -20
    grep -c 'hisab' "$b" | sed 's/^/lines naming hisab: /'
done
say "cyrius distlib --check"
(cd "$S" && cyrius distlib --check 2>&1 | grep -v '^\s*$' | sed "s|$S/||g")

# A consumer project: $1 name, $2 bundle module, $3 stdlib list (TOML array body, may be empty), $4 extra build flags.
consumer() {
    local c=$OUT/$1
    mkdir -p "$c/src"
    cp "$here/scripts/probes/p10_manifest/consumer.cyr" "$c/src/main.cyr"
    cat > "$c/cyrius.cyml" <<EOF
[package]
name = "$1"
version = "0.0.1"
language = "cyrius"
cyrius = "6.7.18"

[build]
entry = "src/main.cyr"
output = "build/$1"
dce = true

[deps]
stdlib = [$3]

[deps.ushma]
path = "../ushma"
modules = ["$2"]
EOF
    say "3. consumer $1 ($2; stdlib [$3])"
    (cd "$c" && cyrius deps 2>&1 | grep -v '^\s*$' | grep -i -E 'warn|error|refus|sakshi|hisab' | sort -u)
    local t
    for t in elf win aarch64; do
        local flag="" ext=""
        [ "$t" = win ] && flag=--win && ext=.exe
        [ "$t" = aarch64 ] && flag=--aarch64
        # shellcheck disable=SC2086
        local log
        log=$(cd "$c" && cyrius build $flag $4 src/main.cyr "build/$1-$t$ext" 2>&1)
        local rc=$?
        printf '%-8s build rc %d: %s\n' "$t" "$rc" "$(printf '%s\n' "$log" | grep -E '^OK|error|undefined|warning' \
            | grep -v -E 'large static data|refusing to overwrite stdlib leaf|^warning:lib/' | tr '\n' ' ')"
        [ $rc -ne 0 ] && continue
        local run
        case $t in
            elf) run=$("$c/build/$1-$t" 2>&1) ;;
            win) run=$(wine "$c/build/$1-$t$ext" 2>/dev/null | tr -d '\r') ;;
            aarch64) run=$(qemu-aarch64 "$c/build/$1-$t" 2>&1) ;;
        esac
        local rrc=$?
        printf '%-8s run rc %d: %s\n' "$t" $rrc "$(printf '%s' "$run" | tr '\n' ' ' | cut -c1-300)"
    done
}

core_deps=$(grep -v '^#' "$S/dist/ushma.deps" 2>/dev/null | sed 's/.*/"&"/' | paste -sd, -)
ai_deps=$(grep -v '^#' "$S/dist/ushma-ai.deps" 2>/dev/null | sed 's/.*/"&"/' | paste -sd, -)
consumer c_core "dist/ushma.cyr" "$core_deps" ""
consumer c_core_bare "dist/ushma.cyr" "" ""
python3 -I -B "$here/scripts/probes/p7_daimon_stub.py" 18731 "$OUT/stub.log" &
stub=$!
python3 -c 'import time; time.sleep(0.5)'
consumer c_ai "dist/ushma-ai.cyr" "$ai_deps" "-D P10_AI"
consumer c_ai_bare "dist/ushma-ai.cyr" "" "-D P10_AI"
kill $stub 2>/dev/null
wait $stub 2>/dev/null
say "stub saw"
grep -c '^POST /v1/agents HTTP/1.1' "$OUT/stub.log" 2>/dev/null | sed 's/$/ registration POSTs/' || echo "0 registration POSTs"
