#!/usr/bin/env python3
"""Compile every generated constant and table cell through the Cyrius compiler and check each value.

    python3 -I -B scripts/oracle/probe_tables.py .                # x86_64 ELF
    python3 -I -B scripts/oracle/probe_tables.py . --all-targets  # + aarch64 (qemu) and PE (wine): identical output

This proves the generated sources reproduce the expected bits as the compiler reads them
(decimal-literal conversion, u16 offsets, the NUL-separated name blobs), not just as Python
reads them. The expected values come from two independent places:

- **golden rows** (the oracle's API) for everything public: `src/constants.cyr` against
  `constants.tsv`; the four public tables (`src/material_data.cyr`, `src/state_data.cyr`,
  `src/phase_data.cyr`, `src/chem_data.cyr`) against their `<module>::<PRESET>` and
  `<module>::<ALL>` rows: the count, every index constant, every name and every f64 cell;
- **this script's own parse of `rust-old/src`** for the private data, which the API can't reach:
  `src/steam_data.cyr` (R4_N, the Region 4 scalars, the saturated and superheated tables, the
  tier pressures, the counts) and `src/constants_private.cyr`. The literals are re-read with
  Python's correctly rounded `float()`, and the one constant expression (`70.0 * 4184.0 /
  86_400.0`) is folded in binary64, left to right, as Rust's const evaluation folds it. It
  shares no code with the Rust generator.

Every `var` and `const` in each generated file must be probed. It writes `build/oracle-probe/`
(gitignored), builds there with the project's `cyrius.cyml`, and exits 1 on any mismatch, any
unprobed item, or any compiler diagnostic other than the vendored ones recorded in
`docs/development/state.md`. Adapted from kimiya's `probe_tables.py`.
"""
import re
import subprocess
import sys
from pathlib import Path

args = [a for a in sys.argv[1:] if not a.startswith("--")]
all_targets = "--all-targets" in sys.argv[1:]
root = Path(args[0] if args else ".").resolve()
work = root / "build/oracle-probe"
work.mkdir(parents=True, exist_ok=True)

# Diagnostics from vendored lib/, recorded in docs/development/state.md.
KNOWN = re.compile(r"large static data|refusing to overwrite stdlib leaf 'sakshi'|^\s+(kept|skipped) \(|"
                   r"declare it in \[deps\] stdlib|^warning:lib/")

checks = []   # (statement printing one value, kind 'f' | 'i' | 's', expected text, label)
probed = {}   # file -> set of names its checks read
fails = 0


def fail(msg):
    global fails
    fails += 1
    if fails <= 40:
        print("FAIL", msg)


def check(file, names, stmt, kind, want, label):
    probed.setdefault(file, set()).update(names)
    checks.append((stmt, kind, want, label))


def rows(name):
    path = root / f"tests/golden/{name}.tsv"
    for line in path.read_text().splitlines():
        if line and not line.startswith("#"):
            fn, cls, ins, outs, _note = line.split("\t")
            yield fn, cls, ([] if ins == "-" else ins.split(" ")), ([] if outs == "-" else outs.split(" "))


def unpct(t):
    return re.sub(rb"%([0-9A-F]{2})", lambda m: bytes([int(m.group(1), 16)]), t.encode()).decode("utf-8")


def f64_bits_hex(x):
    import struct
    return format(struct.unpack("<Q", struct.pack("<d", x))[0], "x")


def rust_float(lit):
    return float(lit.strip().replace("_", ""))


def strip_comments(text):
    return "\n".join(line.split("//", 1)[0] for line in text.splitlines())


def nontest(src):
    at = src.find("#[cfg(test)]")
    return src if at < 0 else src[:at]


# ── src/constants.cyr against constants.tsv ─────────────────────────────────────────────────
golden_const = {ins[0][2:]: outs[0] for fn, cls, ins, outs in rows("constants")}
cyr = (root / "src/constants.cyr").read_text()
for path, name in re.findall(r"^# Rust: ushma::(\S+)\nconst (US_[A-Z0-9_]+) = ", cyr, re.M):
    want = golden_const.get(path)
    if want is None:
        fail(f"src/constants.cyr {name}: ushma::{path} has no golden row")
        continue
    if want.startswith("f:"):
        check("src/constants.cyr", {name}, f"x = {name}; fmt_hex(x);", "f", want[2:], f"const {path}")
    else:
        check("src/constants.cyr", {name}, f"fmt_int({name});", "i", want[2:], f"const {path}")


# ── the public tables against their golden rows ─────────────────────────────────────────────
def struct_fields(src, ty):
    m = re.search(rf"pub struct {ty} \{{(.*?)\n\}}", src, re.S)
    if not m:
        raise SystemExit(f"no struct {ty}")
    return re.findall(r"^\s*pub (\w+): (\S+),", m.group(1), re.M)


TABLES = [  # (file, module, struct, prefix, ALL slice, prefix the index constant drops)
    ("src/material_data.cyr", "material", "ThermalMaterial", "MATERIAL", "ALL_MATERIALS", ""),
    ("src/state_data.cyr", "state", "GasData", "GAS", "ALL_GASES", "GAS_"),
    ("src/phase_data.cyr", "phase", "SubstanceData", "SUBSTANCE", "ALL_SUBSTANCES", ""),
    ("src/chem_data.cyr", "chem", "Species", "SPECIES", "ALL_SPECIES", ""),
]
table_cells = {}
for file, module, ty, prefix, all_name, strip in TABLES:
    fields = [(n, t) for n, t in struct_fields((root / f"rust-old/src/{module}.rs").read_text(), ty)]
    if fields[0][0] != "name" or any(t != "f64" for _, t in fields[1:]):
        fail(f"{ty}: expected `name` then f64 fields, found {fields}")
        continue
    f64_fields = [n for n, _ in fields[1:]]
    presets, order = {}, None
    for fn, cls, ins, outs in rows(module):
        if fn == f"{module}::{all_name}":
            order = [t[2:] for t in outs[1:]]
            if outs[0] != f"v:{len(order)}":
                fail(f"{fn}: {outs[0]} but {len(order)} names")
        elif fn.startswith(f"{module}::") and not ins:
            presets[fn.split("::", 1)[1]] = outs
    if order is None:
        fail(f"tests/golden/{module}.tsv: no {module}::{all_name} row")
        continue
    check(file, {f"US_{prefix}_COUNT"}, f"fmt_int(US_{prefix}_COUNT);", "i", str(len(order)), f"US_{prefix}_COUNT")
    for r, rust in enumerate(order):
        outs = presets.get(rust)
        if outs is None or len(outs) != 1 + len(f64_fields):
            fail(f"{module}::{rust}: golden row missing or the wrong width")
            continue
        idx = f"US_{prefix}_{rust[len(strip):] if rust.startswith(strip) else rust}"
        check(file, {idx}, f"fmt_int({idx});", "i", str(r), f"{idx} (row of {module}::{rust})")
        blob, off = f"_US_{prefix}_NAME_BYTES", f"_US_{prefix}_NAME_OFF"
        check(file, {blob, off}, f"p = &{blob} + {off}[{idx}]; print(p, strlen(p));", "s", unpct(outs[0][2:]),
              f"{module}::{rust}.name")
        for j, field in enumerate(f64_fields):
            col = f"US_{prefix}_{field.upper()}"
            check(file, {col}, f"x = {col}[{idx}]; fmt_hex(x);", "f", outs[1 + j][2:], f"{module}::{rust}.{field}")
    table_cells[file] = len(order) * (len(f64_fields) + 2) + 1


# ── src/steam_data.cyr against this script's parse of rust-old/src/steam.rs ────────────────
steam = strip_comments(nontest((root / "rust-old/src/steam.rs").read_text()))
SF = "src/steam_data.cyr"


def block(name):
    m = re.search(rf"const {name}: [^=]+= &?\[(.*?)\n\];", steam, re.S)
    if not m:
        raise SystemExit(f"steam.rs: no {name}")
    return m.group(1)


def entries(text, ty):
    out = []
    for body in re.findall(rf"{ty} \{{(.*?)\}}", text, re.S):
        out.append([(f, rust_float(v)) for f, v in re.findall(r"(\w+):\s*([^,]+),", body)])
    return out


r4 = [rust_float(v) for v in block("R4_N").split(",") if v.strip()]
if len(r4) != 10:
    fail(f"steam.rs R4_N: {len(r4)} coefficients")
for k, v in enumerate(r4):
    check(SF, {f"_US_STEAM_R4_N{k + 1}"}, f"x = _US_STEAM_R4_N{k + 1}; fmt_hex(x);", "f", f64_bits_hex(v), f"R4_N[{k}]")
for name in ("R4_T_MIN", "T_CRIT", "P_CRIT", "QUALITY_ROUNDOFF"):
    v = rust_float(re.search(rf"const {name}: f64 = ([^;]+);", steam).group(1))
    check(SF, {f"_US_STEAM_{name}"}, f"x = _US_STEAM_{name}; fmt_hex(x);", "f", f64_bits_hex(v), f"steam::{name}")
per = int(re.search(r"const TEMPS_PER_TIER: usize = ([^;]+);", steam).group(1).replace("_", ""))
check(SF, {"_US_STEAM_TEMPS_PER_TIER"}, "fmt_int(_US_STEAM_TEMPS_PER_TIER);", "i", str(per), "TEMPS_PER_TIER")
for table, ty, prefix, count in (("SATURATED_TABLE", "SaturatedEntry", "SAT", "SATURATED_COUNT"),
                                 ("SUPERHEATED_DATA", "SuperheatedEntry", "SH", "SUPERHEATED_COUNT")):
    ents = entries(block(table), ty)
    check(SF, {f"_US_STEAM_{count}"}, f"fmt_int(_US_STEAM_{count});", "i", str(len(ents)), f"len({table})")
    for k, e in enumerate(ents):
        for field, v in e:
            col = f"_US_STEAM_{prefix}_{field.upper()}"
            check(SF, {col}, f"x = {col}[{k}]; fmt_hex(x);", "f", f64_bits_hex(v), f"{table}[{k}].{field}")
tiers = [rust_float(v) for v in block("SUPERHEAT_PRESSURES").split(",") if v.strip()]
check(SF, {"_US_STEAM_SUPERHEAT_TIER_COUNT"}, "fmt_int(_US_STEAM_SUPERHEAT_TIER_COUNT);", "i", str(len(tiers)),
      "len(SUPERHEAT_PRESSURES)")
for k, v in enumerate(tiers):
    check(SF, {"_US_STEAM_SUPERHEAT_PRESSURES"}, f"x = _US_STEAM_SUPERHEAT_PRESSURES[{k}]; fmt_hex(x);", "f",
          f64_bits_hex(v), f"SUPERHEAT_PRESSURES[{k}]")


# ── src/constants_private.cyr against this script's parse of rust-old/src ───────────────────
PF = "src/constants_private.cyr"


def fold(expr):
    """A literal, or literals joined by * and /, folded left to right in binary64 (Rust's const evaluation)."""
    toks = re.split(r"\s*([*/])\s*", expr.strip())
    acc = rust_float(toks[0])
    for op, lit in zip(toks[1::2], toks[2::2]):
        acc = acc * rust_float(lit) if op == "*" else acc / rust_float(lit)
    return acc


pcyr = (root / PF).read_text()
for module, name in re.findall(r"^# Rust: ushma::(\w+)::(?:\w+::)?(\w+) \(private", pcyr, re.M):
    src = strip_comments(nontest((root / f"rust-old/src/{module}.rs").read_text()))
    m = re.search(rf"\bconst {name}: (\w+) = ([^;]+);", src)
    if not m:
        fail(f"{PF}: no `const {name}` in rust-old/src/{module}.rs")
        continue
    ty, expr = m.groups()
    cname = f"_US_{module.upper()}_{name}"
    if ty == "f64":
        check(PF, {cname}, f"x = {cname}; fmt_hex(x);", "f", f64_bits_hex(fold(expr)), f"{module}::{name}")
    else:
        val = (int(expr.split("<<")[0].replace("_", "")) << int(expr.split("<<")[1].replace("_", ""))
               if "<<" in expr else int(expr.replace("_", "")))
        check(PF, {cname}, f"fmt_int({cname});", "i", str(val), f"{module}::{name}")


# ── coverage: every var and const of every generated file is probed ─────────────────────────
GENERATED = ["src/constants.cyr", PF, SF] + [t[0] for t in TABLES]
for file in GENERATED:
    text = (root / file).read_text()
    declared = set(re.findall(r"^(?:var|const) (_?US_[A-Z0-9_]+)", text, re.M))
    missing = declared - probed.get(file, set())
    gone = probed.get(file, set()) - declared
    if missing or gone:
        fail(f"{file}: not probed {sorted(missing)}; probed but not declared {sorted(gone)}")


# ── the probe program ───────────────────────────────────────────────────────────────────────
src = [f'include "{f}"' for f in GENERATED] + ["", "fn main() {", "    var x: f64 = 0.0;", "    var p = 0;"]
for k, (stmt, kind, want, label) in enumerate(checks):
    src.append(f'    {stmt} println(" {k}");')
src += ["    return 0;", "}", "", "var r = main();", "sys_exit_group(r);", ""]
probe = work / "tables_probe.cyr"
probe.write_text("\n".join(src))

b = subprocess.run(["cyrius", "build", str(probe), str(work / "tables_probe")], capture_output=True, text=True,
                   cwd=root)
diag = [line for line in (b.stdout + b.stderr).splitlines()
        if re.search(r"warn|error|refus", line, re.I) and not KNOWN.search(line)]
if b.returncode != 0 or diag:
    print("BUILD:", b.returncode, *diag, sep="\n")
    sys.exit(2)
r = subprocess.run([str(work / "tables_probe")], capture_output=True, text=True)
got = {}
for line in r.stdout.splitlines():
    if " " in line:
        val, k = line.rsplit(" ", 1)
        if k.isdigit():
            got[int(k)] = val

for k, (stmt, kind, want, label) in enumerate(checks):
    g = got.get(k)
    if kind == "f":
        ok = g is not None and re.fullmatch(r"[0-9a-fA-F]+", g) is not None and int(g, 16) == int(want, 16)
    else:
        ok = g == want
    if not ok:
        fail(f"{label}: {stmt} -> {g!r}, want {want!r}")

# The cells are decimal literals the compiler converts, so every target must print the same bytes.
targets = "x86_64 ELF"
if all_targets:
    for flag, exe, runner in (("--aarch64", "tables_probe_a64", ["qemu-aarch64"]),
                              ("--win", "tables_probe.exe", ["wine"])):
        bt = subprocess.run(["cyrius", "build", flag, str(probe), str(work / exe)], capture_output=True, text=True,
                            cwd=root)
        if bt.returncode != 0:
            fail(f"build {flag}: exit {bt.returncode}")
            continue
        rt = subprocess.run(runner + [str(work / exe)], capture_output=True, text=True)
        if rt.stdout.replace("\r", "") != r.stdout:
            fail(f"{flag}: output differs from x86_64 ELF")
    targets = "x86_64 ELF, aarch64 and PE (identical)"
print("checked " + ", ".join(f"{Path(f).name}" for f in GENERATED) + f": {len(checks)} values through the "
      f"Cyrius compiler on {targets}; failures: {fails}")
sys.exit(1 if fails else 0)
