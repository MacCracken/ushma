#!/usr/bin/env python3
"""Independent check of ushma's golden corpus (port plan §4 step 5).

    python3 -I -B scripts/oracle/validate.py .

It shares no code with the Rust generator. It re-derives what it can from first principles and
compares bit for bit:

- the format: every `tests/golden/*.tsv` has the format-1 header, a correct row count, five
  tab-separated columns, a class letter and well-formed tokens, and `modules.txt` lists exactly
  the module files;
- `constants.tsv`: every row equals the Rust literal in `rust-old/src/<module>.rs`, parsed by
  Python's correctly rounded `float()`;
- `src/constants.cyr`: every `const US_<NAME> = <literal>;` under a `# Rust: ushma::<path>` line
  parses to that constant's golden bits;
- the public table rows (`<module>::<PRESET>` in `chem`, `material`, `phase`, `state`): every
  preset's name and fields equal the struct literal in `rust-old/src/<module>.rs`, and the
  `<module>::<ALL>` row names exactly the slice's entries in order;
- `math.tsv`, classes A and B: sqrt, abs, ceil, min, max, clamp, is_finite, is_nan and powi are
  recomputed in IEEE double arithmetic (Python floats are binary64; powi follows
  `__powidf2`'s order);
- the corrected oracle (ADR 0004), `tests/golden/math.cr.tsv`: every class-C math row (exp, exp_m1,
  ln, powf, atan, tanh) is recomputed correctly rounded with MPFR (`mpfr_ref.py`). A row whose
  oracle output (glibc's) differs gets a correction: the golden row's fn and inputs byte for byte,
  the correctly rounded outputs, and a note naming the oracle's outputs and the call glibc
  misrounded. The file is checked against that recomputation. Each listed call is checked against
  MPFR directly, and that glibc's value differs from it. On the corpus's recorded libm platform it
  is also checked that glibc's value is this machine's libm result, and the whole file text.

    python3 -I -B scripts/oracle/validate.py .                       # check
    python3 -I -B scripts/oracle/validate.py . --write-corrections   # (re)write tests/golden/*.cr.tsv

Exit status 1 on any failure; the last line is `FAILURES: <n>`.
"""

import importlib.util
import math
import pathlib
import platform
import re
import struct
import sys

# mpfr_ref.py sits beside this script; `-I` keeps the script's directory off sys.path, so load it by path.
_mpfr_spec = importlib.util.spec_from_file_location("mpfr_ref", pathlib.Path(__file__).resolve().parent / "mpfr_ref.py")
MPFR = importlib.util.module_from_spec(_mpfr_spec)
_mpfr_spec.loader.exec_module(MPFR)

TOKEN = re.compile(r"^(f:[0-9a-f]{16}|i:-?[0-9]+|e:[A-Za-z0-9_]+\.[A-Za-z0-9_]+|s:\S*|n:|E:[A-Za-z0-9_]+|v:[0-9]+)$")


def pct(text):
    """Format 1's `s:` encoding: every byte outside printable ASCII (0x21-0x7e), and `%`, is `%XX`."""
    return "".join(chr(b) if 0x21 <= b <= 0x7E and b != 0x25 else f"%{b:02X}" for b in text.encode())


def bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(b):
    return struct.unpack("<d", struct.pack("<Q", b))[0]


def tok_f(tok):
    if not tok.startswith("f:"):
        raise ValueError(f"expected an f: token, got {tok}")
    return from_bits(int(tok[2:], 16))


def same(x, y):
    """Class A equality: identical bits, or both NaN (any NaN matches any NaN)."""
    if math.isnan(x) and math.isnan(y):
        return True
    return bits(x) == bits(y)


class Report:
    def __init__(self):
        self.failures = 0
        self.checked = 0

    def fail(self, where, msg):
        self.failures += 1
        if self.failures <= 40:
            print(f"FAIL {where}: {msg}")

    def ok(self, n=1):
        self.checked += n


def read_tsv(path, rep):
    """The rows of one golden file, after checking its header and format."""
    lines = path.read_text().splitlines()
    name = path.name[: -len(".tsv")]
    want = [
        f"# ushma golden corpus, format 1: {name}",
        "# oracle: ",
        "# regenerate: cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml",
        "# columns: fn<TAB>class<TAB>inputs<TAB>outputs<TAB>note (tokens: scripts/oracle/README.md)",
        "# rows: ",
    ]
    for k, w in enumerate(want):
        if k >= len(lines) or not lines[k].startswith(w):
            rep.fail(path.name, f"header line {k + 1} should start {w!r}")
            return []
    if "libm: " not in lines[1]:
        rep.fail(path.name, "the oracle line names no libm platform")
    rows = [ln for ln in lines[5:] if ln]
    if lines[4] != f"# rows: {len(rows)}":
        rep.fail(path.name, f"{lines[4]!r}, but {len(rows)} rows follow")
    out = []
    for n, ln in enumerate(rows, start=6):
        cols = ln.split("\t")
        if len(cols) != 5:
            rep.fail(f"{path.name}:{n}", f"{len(cols)} columns")
            continue
        fn, cls, ins, outs, _note = cols
        if cls not in "ABCD" or len(cls) != 1:
            rep.fail(f"{path.name}:{n}", f"class {cls!r}")
            continue
        toks_in = [] if ins == "-" else ins.split(" ")
        toks_out = [] if outs == "-" else outs.split(" ")
        bad = [t for t in toks_in + toks_out if not TOKEN.match(t)]
        if bad:
            rep.fail(f"{path.name}:{n}", f"malformed tokens {bad}")
            continue
        out.append((n, fn, cls, toks_in, toks_out))
    return out


def rust_literal(text):
    """A Rust numeric literal's value: underscores dropped, then Python's correctly rounded parse."""
    t = text.replace("_", "")
    if re.fullmatch(r"-?[0-9]+", t):
        return int(t)
    return float(t)


def check_constants(root, rows, rep):
    golden = {}
    for n, fn, cls, ins, outs in rows:
        where = f"constants.tsv:{n}"
        if fn != "const" or cls != "A" or len(ins) != 1 or len(outs) != 1 or not ins[0].startswith("s:"):
            rep.fail(where, "not a `const A s:<path> <value>` row")
            continue
        path = ins[0][2:]
        module, name = path.split("::")
        src = (root / "rust-old/src" / f"{module}.rs").read_text()
        m = re.search(rf"^pub const {name}: (\w+) = ([^;]+);", src, re.M)
        if not m:
            rep.fail(where, f"no `pub const {name}` literal in rust-old/src/{module}.rs")
            continue
        value = rust_literal(m.group(2).strip())
        if isinstance(value, int):
            ok = outs[0] == f"i:{value}"
        else:
            ok = outs[0].startswith("f:") and bits(tok_f(outs[0])) == bits(value)
        if ok:
            rep.ok()
        else:
            rep.fail(where, f"{path}: golden {outs[0]}, Rust literal {m.group(2)}")
        golden[path] = outs[0]
    return golden


def check_constants_cyr(root, golden, rep):
    text = (root / "src/constants.cyr").read_text()
    pairs = re.findall(r"^# Rust: ushma::(\S+)\nconst US_(\w+) = ([^;]+);$", text, re.M)
    if len(pairs) != len(golden):
        rep.fail("src/constants.cyr", f"{len(pairs)} constants, the corpus has {len(golden)}")
    for path, name, lit in pairs:
        where = f"src/constants.cyr US_{name}"
        if path.split("::")[1] != name:
            rep.fail(where, f"named for {path}")
            continue
        if path not in golden:
            rep.fail(where, f"{path} has no golden row")
            continue
        if re.search(r"[eExX]", lit):
            rep.fail(where, f"{lit} is not a positional decimal")
            continue
        g = golden[path]
        if g.startswith("i:"):
            ok = lit == g[2:]
        else:
            ok = bits(float(lit)) == int(g[2:], 16)
        if ok:
            rep.ok()
        else:
            rep.fail(where, f"{lit} is not the golden {g}")


# ── the class A / B math functions, recomputed ──

def rust_sqrt(x):
    if math.isnan(x) or x < 0:
        return math.nan
    if x == 0 or math.isinf(x):
        return x  # keeps -0.0 and +inf
    return math.sqrt(x)


def rust_abs(x):
    return from_bits(bits(x) & 0x7FFF_FFFF_FFFF_FFFF)


def rust_ceil(x):
    if math.isnan(x) or math.isinf(x) or x == math.floor(x):
        return x  # integers, ±0, ±inf and NaN are their own ceiling
    r = float(math.ceil(x))
    return -0.0 if r == 0 and x < 0 else r


def rust_max(x, y):
    """Rust 1.99 f64::max as recorded: a NaN operand is ignored; (±0, ∓0) gives the second operand."""
    if math.isnan(x):
        return y
    if math.isnan(y):
        return x
    if x == y:
        return y
    return x if x > y else y


def rust_min(x, y):
    if math.isnan(x):
        return y
    if math.isnan(y):
        return x
    if x == y:
        return y
    return x if x < y else y


def rust_clamp(x, lo, hi):
    if x < lo:
        x = lo
    if x > hi:
        x = hi
    return x


def powidf2(a, b):
    """compiler-builtins `__powidf2`: square-and-multiply, C division, then 1/r for b < 0."""
    recip = b < 0
    r = 1.0
    while True:
        if b & 1:
            r *= a
        b = int(b / 2)
        if b == 0:
            break
        a *= a
    if recip:
        # IEEE 1/±0 is ±inf (Python raises on a float division by zero, so it is spelled out).
        return math.copysign(math.inf, r) if r == 0 else 1.0 / r
    return r


def check_math(rows, rep):
    unary = {
        "f64::sqrt": rust_sqrt,
        "f64::abs": rust_abs,
        "f64::ceil": rust_ceil,
    }
    counts = {}
    for n, fn, cls, ins, outs in rows:
        where = f"math.tsv:{n}"
        counts[cls] = counts.get(cls, 0) + 1
        try:
            if fn in unary:
                got, want = tok_f(outs[0]), unary[fn](tok_f(ins[0]))
            elif fn == "f64::max":
                got, want = tok_f(outs[0]), rust_max(tok_f(ins[0]), tok_f(ins[1]))
            elif fn == "f64::min":
                got, want = tok_f(outs[0]), rust_min(tok_f(ins[0]), tok_f(ins[1]))
            elif fn == "f64::clamp":
                got, want = tok_f(outs[0]), rust_clamp(tok_f(ins[0]), tok_f(ins[1]), tok_f(ins[2]))
            elif fn == "f64::is_finite":
                x = tok_f(ins[0])
                got, want = outs[0], f"i:{int(math.isfinite(x))}"
            elif fn == "f64::is_nan":
                x = tok_f(ins[0])
                got, want = outs[0], f"i:{int(math.isnan(x))}"
            elif fn == "f64::powi":
                got, want = tok_f(outs[0]), powidf2(tok_f(ins[0]), int(ins[1][2:]))
            else:
                if cls in "AB":
                    rep.fail(where, f"{fn} is class {cls} but has no independent model here")
                continue
        except (ValueError, IndexError, OverflowError) as e:
            rep.fail(where, f"{fn}: {e}")
            continue
        if cls not in "AB":
            rep.fail(where, f"{fn} is modelled as exact but marked class {cls}")
            continue
        if isinstance(got, str):
            ok = got == want
        else:
            ok = same(got, want)
        if ok:
            rep.ok()
        else:
            rep.fail(where, f"{fn}{tuple(ins)} golden {outs[0]}, model {want!r}")
    return counts


TABLES = [  # (module, struct, ALL slice)
    ("material", "ThermalMaterial", "ALL_MATERIALS"),
    ("state", "GasData", "ALL_GASES"),
    ("phase", "SubstanceData", "ALL_SUBSTANCES"),
    ("chem", "Species", "ALL_SPECIES"),
]


def check_tables(root, tables, rep):
    """The public presets' golden rows against their struct literals in the Rust source."""
    for module, ty, all_name in TABLES:
        src = (root / "rust-old/src" / f"{module}.rs").read_text()
        src = src[: src.find("#[cfg(test)]")] if "#[cfg(test)]" in src else src
        literals = {}
        for name, body in re.findall(rf"^pub const (\w+): {ty} = {ty} \{{(.*?)\n\}};", src, re.M | re.S):
            body = "\n".join(ln.split("//", 1)[0] for ln in body.splitlines())
            fields = re.findall(r"(\w+):\s*([^,]+),", body)
            if not fields or fields[0][0] != "name":
                rep.fail(f"{module}.rs {name}", "no `name` field first")
                continue
            m = re.fullmatch(r'Cow::Borrowed\("([^"]*)"\)', fields[0][1].strip())
            if not m:
                rep.fail(f"{module}.rs {name}", f"name is {fields[0][1]!r}")
                continue
            literals[name] = (m.group(1), [float(v.strip().replace("_", "")) for _, v in fields[1:]])
        all_m = re.search(rf"^pub const {all_name}: &\[&{ty}\] = &\[(.*?)\];", src, re.M | re.S)
        all_order = re.findall(r"&(\w+)", all_m.group(1)) if all_m else []
        seen = set()
        for n, fn, cls, ins, outs in tables.get(f"{module}.tsv", []):
            where = f"{module}.tsv:{n}"
            if not fn.startswith(f"{module}::") or ins:
                continue
            rust = fn.split("::", 1)[1]
            if rust == all_name:
                names = [t[2:] for t in outs[1:]]
                if outs[0] != f"v:{len(names)}" or names != all_order:
                    rep.fail(where, f"{fn} names {names}, the slice is {all_order}")
                else:
                    rep.ok()
                continue
            if rust not in literals:
                continue  # a function row (M1+), not a preset
            seen.add(rust)
            name, vals = literals[rust]
            want = [f"s:{pct(name)}"] + [f"f:{bits(v):016x}" for v in vals]
            if outs == want:
                rep.ok()
            else:
                rep.fail(where, f"{fn}: golden {outs}, Rust literal {want}")
        if seen != set(literals):
            rep.fail(f"{module}.tsv", f"presets without a golden row: {sorted(set(literals) - seen)}")


# ── the corrected oracle (ADR 0004) ──

CR_FN = {  # golden fn -> (MPFR function, the libm function glibc's value came from)
    "f64::exp": ("exp", math.exp),
    "f64::exp_m1": ("expm1", math.expm1),
    "f64::ln": ("log", math.log),
    "f64::powf": ("pow", math.pow),
    "f64::atan": ("atan", math.atan),
    "f64::tanh": ("tanh", math.tanh),
}
CR_SUFFIX = ".cr.tsv"


def running_platform():
    """`glibc 2.44 x86_64 fma`, as scripts/oracle/src/glibc.rs names it."""
    _lib, ver = platform.libc_ver()
    arch = platform.machine()
    cpu = "generic"
    if arch == "x86_64":
        flags = set()
        try:
            for line in open("/proc/cpuinfo"):
                if line.startswith("flags"):
                    flags = set(line.split(":", 1)[1].split())
                    break
        except OSError:
            pass
        cpu = "fma" if {"fma", "avx2"} <= flags else "no-fma"
    return f"glibc {ver} {arch} {cpu}"


def recorded_platform(text):
    head = text.split("\n", 2)[1]
    return head.split("libm: ", 1)[1].split(". ", 1)[0]


def fbits(x):
    return f"f:{bits(x):016x}"


def toks_text(toks):
    return " ".join(toks) if toks else "-"


def mpfr_call(name, args):
    return MPFR.cr(name, args[0]) if len(args) == 1 else MPFR.cr2(name, args[0], args[1])


def corrections(rows):
    """The rows of math.tsv whose correctly rounded outputs differ from the oracle's: (fn, inputs, outputs, note)."""
    out = []
    for n, fn, cls, ins, outs in rows:
        if fn not in CR_FN or cls != "C":
            continue
        name = CR_FN[fn][0]
        args = [tok_f(t) for t in ins]
        g = tok_f(outs[0])
        c = mpfr_call(name, args)
        if same(c, g):
            continue
        call = f"{name}({', '.join(fbits(a) for a in args)}) glibc {fbits(g)} correctly rounded {fbits(c)}"
        out.append((fn, toks_text(ins), fbits(c), f"oracle {toks_text(outs)}; {call}"))
    return out


def cr_text(mod, golden_text, rows):
    prov = golden_text.split("\n", 2)[1].split(". Generated by", 1)[0]
    head = [
        f"# ushma corrected oracle (ADR 0004), format 1: {mod}",
        f"{prov}. Correctly rounded: MPFR {MPFR.VERSION}. Generated by scripts/oracle/validate.py; do not edit.",
        "# regenerate: python3 -I -B scripts/oracle/validate.py . --write-corrections",
        "# columns: fn<TAB>inputs<TAB>outputs<TAB>note (fn, inputs: the golden row's, byte for byte; outputs: "
        "the oracle's operations with every primitive result correctly rounded; note: the oracle's outputs, then "
        "each call glibc misrounded)",
        f"# rows: {len(rows)}",
    ]
    return "\n".join(head + ["\t".join(r) for r in rows]) + "\n"


CALL = re.compile(r"(\w+)\((f:[0-9a-f]{16}(?:, f:[0-9a-f]{16})?)\) glibc f:([0-9a-f]{16}) correctly rounded "
                  r"f:([0-9a-f]{16})")


def check_corrections(root, tables, rep, write):
    """Every module's corrected oracle (math, and each modelled module): its `.cr.tsv` against the MPFR
    recomputation, each listed call against MPFR and glibc. Returns the corrected rows per function."""
    gold = root / "tests/golden"
    mods = ["math"] + [m for m in MODELS if f"{m}.tsv" in tables]
    for stale in sorted(gold.glob("*" + CR_SUFFIX)):
        if stale.name[: -len(CR_SUFFIX)] not in mods:
            rep.fail(stale.name, "no such modelled module")
    per = {}
    for mod in mods:
        rows = tables.get(f"{mod}.tsv", [])
        want = corrections(rows) if mod == "math" else model_corrections(mod, rows)
        for fn, k in check_module_corrections(root, mod, rows, want, rep, write).items():
            per[fn] = per.get(fn, 0) + k
    return per


def check_module_corrections(root, mod, rows, want, rep, write):
    """tests/golden/<mod>.cr.tsv against `want`; each listed call against MPFR and glibc."""
    gold = root / "tests/golden"
    golden_text = (gold / f"{mod}.tsv").read_text()
    exact = recorded_platform(golden_text) == running_platform()
    path = gold / f"{mod}{CR_SUFFIX}"
    if write:
        if not exact:
            rep.fail("--write-corrections", f"this machine ({running_platform()}) is not the corpus's libm platform "
                     f"({recorded_platform(golden_text)}); nothing written")
        else:
            if want or path.exists():
                path.write_text(cr_text(mod, golden_text, want))
                print(f"wrote tests/golden/{mod}{CR_SUFFIX}: {len(want)} rows")
    if not path.exists():
        if want:
            rep.fail(path.name, f"missing; {len(want)} rows need a correction (run with --write-corrections)")
        return {}
    text = path.read_text()
    lines = text.split("\n")
    if lines[-1] != "":
        rep.fail(path.name, "no final newline")
    body = [ln for ln in lines[5:] if ln and not ln.startswith("#")]
    golden_keys = {(fn, toks_text(ins)): toks_text(outs) for n, fn, cls, ins, outs in rows}
    got = []
    for k, line in enumerate(body, 6):
        cols = line.split("\t")
        if len(cols) != 4:
            rep.fail(f"{path.name}:{k}", f"{len(cols)} fields")
            continue
        fn, ins, outs, note = cols
        got.append((fn, ins, outs))
        if (fn, ins) not in golden_keys:
            rep.fail(f"{path.name}:{k}", f"no golden row {fn} {ins}")
            continue
        parts = note.split("; ")
        if parts[0] != f"oracle {golden_keys[(fn, ins)]}":
            rep.fail(f"{path.name}:{k}", f"the note's oracle outputs {parts[0]!r} are not the golden row's")
        for call in parts[1:]:
            m = CALL.fullmatch(call)
            if not m or m.group(1) not in {v[0] for v in CR_FN.values()}:
                rep.fail(f"{path.name}:{k}", f"malformed call {call!r}")
                continue
            name = m.group(1)
            args = [from_bits(int(a[2:], 16)) for a in m.group(2).split(", ")]
            g, c = from_bits(int(m.group(3), 16)), from_bits(int(m.group(4), 16))
            ok = same(c, mpfr_call(name, args)) and not same(g, c)
            if mod == "math":
                ok = ok and outs == fbits(c)
            if exact:
                lib = {v[0]: v[1] for v in CR_FN.values()}[name]
                try:
                    ok = ok and same(g, lib(*args))
                except (OverflowError, ValueError):
                    pass  # Python raises where C returns ±inf / NaN; MPFR and the golden row still decide
            if ok:
                rep.ok()
            else:
                rep.fail(f"{path.name}:{k}", f"{call}: not MPFR's correctly rounded value, or glibc's isn't different")
    if got != [w[:3] for w in want]:
        rep.fail(path.name, f"{len(got)} rows, the MPFR recomputation gives {len(want)} (or a different set/order)")
    if exact and text != cr_text(mod, golden_text, want):
        rep.fail(path.name, "differs from the corrected oracle (regenerate with --write-corrections)")
    per = {}
    for fn, *_ in want:
        per[fn] = per.get(fn, 0) + 1
    return per


# ── module functions, modelled (M1c onward) ──
# Each model recomputes one ushma 1.4.0 function from its inputs in Rust's operation order, in binary64 (Python
# floats), with its transcendentals from a Prim: glibc's (ctypes libm, C's own special cases) to reproduce the golden
# row exactly, or MPFR's correctly rounded to give the corrected oracle (ADR 0004). A model's result is its output
# tokens.

import ctypes as _ct

_LIBM = _ct.CDLL("libm.so.6")
for _n in ("exp", "pow", "log", "tanh", "atan"):
    getattr(_LIBM, _n).restype = _ct.c_double
    getattr(_LIBM, _n).argtypes = [_ct.c_double, _ct.c_double] if _n == "pow" else [_ct.c_double]


class Prim:
    """exp, pow, ln, tanh and atan from glibc (cr=False) or MPFR (cr=True); with cr=True, records each call glibc
    misrounds."""

    def __init__(self, cr):
        self.cr = cr
        self.calls = []

    def _one(self, name, args, g, c):
        if not self.cr:
            return g
        if not same(g, c):
            self.calls.append(f"{name}({', '.join(fbits(a) for a in args)}) glibc {fbits(g)} correctly rounded {fbits(c)}")
        return c

    def exp(self, x):
        return self._one("exp", [x], _LIBM.exp(x), MPFR.exp(x))

    def pow(self, x, y):
        return self._one("pow", [x, y], _LIBM.pow(x, y), MPFR.pow(x, y))

    def log(self, x):
        return self._one("log", [x], _LIBM.log(x), MPFR.log(x))

    def tanh(self, x):
        return self._one("tanh", [x], _LIBM.tanh(x), MPFR.tanh(x))

    def atan(self, x):
        return self._one("atan", [x], _LIBM.atan(x), MPFR.atan(x))


def fdiv(a, b):
    """IEEE a / b (Python raises on a zero divisor)."""
    if b == 0.0:
        if a == 0.0 or math.isnan(a):
            return math.nan
        return math.copysign(math.inf, a) * math.copysign(1.0, b)
    return a / b


def fsqrt(x):
    return math.nan if x < 0 or math.isnan(x) else math.sqrt(x)


def fge(x, b):
    return math.isfinite(x) and x >= b


def fgt(x, b):
    return math.isfinite(x) and x > b


def F(x):
    return [fbits(x)]


def E(name):
    return [f"E:{name}"]


TINY = 1e-30


def m_conduction(a, p):
    k, area, th, tc, l = a
    if not fge(k, 0.0):
        return E("InvalidConductivity")
    if not fge(abs(l), TINY):
        return E("DivisionByZero")
    return F(fdiv(k * area * (th - tc), l))


def m_radiation(a, p):
    e, area, ts, tr = a
    if not (0.0 <= e <= 1.0):
        return E("InvalidParameter")
    if not fge(ts, 0.0) or not fge(tr, 0.0):
        return E("InvalidTemperature")
    return F(e * 5.670374419e-8 * area * (powidf2(ts, 4) - powidf2(tr, 4)))


def m_r_cond(a, p):
    k, area, l = a
    if not fgt(k, 0.0):
        return E("InvalidConductivity")
    if not fgt(area, 0.0) or not fge(l, 0.0):
        return E("InvalidParameter")
    d = k * area
    if not fge(abs(d), TINY):
        return E("DivisionByZero")
    return F(fdiv(l, d))


def m_r_conv(a, p):
    h, area = a
    if not fge(h, 0.0) or not fge(area, 0.0):
        return E("InvalidParameter")
    d = h * area
    if not fge(abs(d), TINY):
        return E("DivisionByZero")
    return F(fdiv(1.0, d))


def m_series(a, p):
    s = -0.0
    for r in a:
        s = s + r
    return F(s)


def m_parallel(a, p):
    for r in a:
        if not fgt(r, 0.0):
            return E("InvalidParameter")
    s = -0.0
    for r in a:
        s = s + fdiv(1.0, r)
    if not fge(abs(s), TINY):
        return E("DivisionByZero")
    return F(fdiv(1.0, s))


def m_diffusivity(a, p):
    k, rho, cp = a
    if not fge(k, 0.0):
        return E("InvalidConductivity")
    if not fge(rho, 0.0):
        return E("InvalidParameter")
    if not fge(cp, 0.0):
        return E("InvalidSpecificHeat")
    d = rho * cp
    if not fge(abs(d), TINY):
        return E("DivisionByZero")
    return F(fdiv(k, d))


def m_biot(a, p):
    h, l, k = a
    if not fge(h, 0.0) or not fge(l, 0.0):
        return E("InvalidParameter")
    if not fge(k, TINY):
        return E("DivisionByZero")
    return F(fdiv(h * l, k))


def m_reynolds(a, p):
    v, l, nu = a
    if not fgt(nu, 0.0):
        return E("InvalidParameter")
    return F(fdiv(v * l, nu))


def m_prandtl(a, p):
    nu, al = a
    if not fgt(nu, 0.0):
        return E("InvalidParameter")
    if not fgt(al, 0.0):
        return E("DivisionByZero")
    return F(fdiv(nu, al))


def m_nusselt(a, p):
    h, l, k = a
    if not fge(k, TINY):
        return E("DivisionByZero")
    return F(fdiv(h * l, k))


def m_dittus(a, p):
    re_, pr = a
    if not fgt(re_, 0.0) or not fgt(pr, 0.0):
        return E("InvalidParameter")
    return F(0.023 * p.pow(re_, 0.8) * p.pow(pr, 0.4))


def m_natural(a, p):
    ra, pr = a
    if not fge(ra, 0.0) or not fgt(pr, 0.0):
        return E("InvalidParameter")
    psi = p.pow(1.0 + p.pow(fdiv(0.492, pr), 9.0 / 16.0), -8.0 / 27.0)
    return F(powidf2(0.825 + 0.387 * p.pow(ra, 1.0 / 6.0) * psi, 2))


def _blasius(x, re_):
    if not fgt(x, 0.0) or not fgt(re_, 0.0):
        return None
    return fdiv(5.0 * x, fsqrt(re_))


def m_blasius(a, p):
    d = _blasius(*a)
    return E("InvalidParameter") if d is None else F(d)


def m_thermal_bl(a, p):
    x, re_, pr = a
    if not fgt(pr, 0.0):
        return E("InvalidParameter")
    d = _blasius(x, re_)
    if d is None:
        return E("InvalidParameter")
    return F(fdiv(d, p.pow(pr, 1.0 / 3.0)))


def _fin_m(h, per, k, ac):
    """fin_parameter: (error name, None) or (None, m)."""
    d = k * ac
    if not fgt(k, 0.0) or not fgt(ac, 0.0):
        return "DivisionByZero", None
    if not fge(h, 0.0) or not fgt(per, 0.0):
        return "InvalidParameter", None
    return None, fsqrt(fdiv(h * per, d))


def m_fin_parameter(a, p):
    e, m = _fin_m(*a)
    return E(e) if e else F(m)


def m_fin_heat(a, p):
    h, per, k, ac, length, tb, tf = a
    e, m = _fin_m(h, per, k, ac)
    if e:
        return E(e)
    if not fgt(length, 0.0):
        return E("InvalidParameter")
    return F(fsqrt(h * per * k * ac) * (tb - tf) * p.tanh(m * length))


def m_fin_efficiency(a, p):
    m, length = a
    if not fge(m, 0.0) or not fge(length, 0.0):
        return E("InvalidParameter")
    ml = m * length
    if abs(ml) < TINY:
        return F(1.0)
    return F(fdiv(p.tanh(ml), ml))


def m_fin_effectiveness(a, p):
    q, h, ac, tb, tf = a
    q0 = h * ac * (tb - tf)
    if not fge(abs(q0), TINY):
        return E("DivisionByZero")
    return F(fdiv(q, q0))


def _lmtd(dt1, dt2, p):
    if not fgt(dt1, 0.0) or not fgt(dt2, 0.0):
        return E("InvalidParameter")
    ratio = fdiv(dt1, dt2)
    if abs(ratio - 1.0) < 1e-10:
        return F(0.5 * (dt1 + dt2))
    return F(fdiv(dt1 - dt2, p.log(ratio)))


def m_ntu(a, p):
    u, area, c = a
    if not fgt(c, 0.0):
        return E("DivisionByZero")
    return F(fdiv(u * area, c))


def _entu_ok(n, cr):
    return fge(n, 0.0) and 0.0 <= cr <= 1.0


def m_eff_parallel(a, p):
    n, cr = a
    if not _entu_ok(n, cr):
        return E("InvalidParameter")
    if cr < TINY:
        return F(1.0 - p.exp(-n))
    return F(fdiv(1.0 - p.exp(-(1.0 + cr) * n), 1.0 + cr))


def m_eff_counter(a, p):
    n, cr = a
    if not _entu_ok(n, cr):
        return E("InvalidParameter")
    if cr < TINY:
        return F(1.0 - p.exp(-n))
    if abs(cr - 1.0) < 1e-10:
        return F(fdiv(n, 1.0 + n))
    e = p.exp(-(1.0 - cr) * n)
    return F(fdiv(1.0 - e, 1.0 - cr * e))


def _dims_ok(a):
    return all(fgt(x, 0.0) for x in a)


def m_vf_parallel(a, p):
    if not _dims_ok(a):
        return E("InvalidParameter")
    w, h, d = a
    x, y = fdiv(w, d), fdiv(h, d)
    x2, y2 = x * x, y * y
    if rust_max(x, y) < 1e-2:
        return F(fdiv(x * y, math.pi) * (1.0 - fdiv(x2 + y2, 3.0)))
    t1 = p.log(fsqrt(fdiv((1.0 + x2) * (1.0 + y2), 1.0 + x2 + y2)))
    sy, sx = fsqrt(1.0 + y2), fsqrt(1.0 + x2)
    t2 = x * sy * p.atan(fdiv(x, sy))
    t3 = y * sx * p.atan(fdiv(y, sx))
    t4, t5 = x * p.atan(x), y * p.atan(y)
    return F(fdiv(2.0, math.pi * x * y) * (t1 + t2 + t3 - t4 - t5))


def m_vf_perpendicular(a, p):
    if not _dims_ok(a):
        return E("InvalidParameter")
    width, h1, h2 = a
    w, h = fdiv(h1, width), fdiv(h2, width)
    w2, hh = w * w, h * h
    s = hh + w2
    rs = fsqrt(s)
    at = w * p.atan(fdiv(1.0, w)) + h * p.atan(fdiv(1.0, h)) - rs * p.atan(fdiv(1.0, rs))
    lt = (
        p.log(fdiv((1.0 + w2) * (1.0 + hh), 1.0 + s))
        + w2 * p.log(fdiv(w2 * (1.0 + s), (1.0 + w2) * s))
        + hh * p.log(fdiv(hh * (1.0 + s), (1.0 + hh) * s))
    )
    return F(rust_clamp(fdiv(at + 0.25 * lt, math.pi * w), 0.0, 1.0))


def m_vf_disks(a, p):
    if not _dims_ok(a):
        return E("InvalidParameter")
    r1, r2, d = a
    r1d, r2d = fdiv(r1, d), fdiv(r2, d)
    s = 1.0 + fdiv(1.0 + r2d * r2d, r1d * r1d)
    k2 = powidf2(fdiv(r2, r1), 2)
    return F(fdiv(2.0 * k2, s + fsqrt(s * s - 4.0 * k2)))


GAS_R = 8.31446261815324  # state::GAS_CONSTANT


def m_ig_pressure(a, p):
    n, temp, v = a
    if not fge(n, 0.0):
        return E("InvalidParameter")
    if not fge(temp, 0.0):
        return E("InvalidTemperature")
    if not fgt(v, 0.0):
        return E("InvalidVolume")
    return F(fdiv(n * GAS_R * temp, v))


def m_ig_volume(a, p):
    n, temp, pr = a
    if not fge(n, 0.0):
        return E("InvalidParameter")
    if not fge(temp, 0.0):
        return E("InvalidTemperature")
    if not fgt(pr, 0.0):
        return E("InvalidPressure")
    return F(fdiv(n * GAS_R * temp, pr))


def m_ig_temperature(a, p):
    pr, v, n = a
    if not fge(n, 0.0):
        return E("InvalidParameter")
    if not fge(pr, 0.0):
        return E("InvalidPressure")
    if not fge(v, 0.0):
        return E("InvalidVolume")
    if not fge(abs(n), TINY):
        return E("DivisionByZero")
    return F(fdiv(pr * v, n * GAS_R))


def m_vdw(a, p):
    n, temp, v, aa, b = a
    if not fge(n, 0.0):
        return E("InvalidParameter")
    if not fgt(temp, 0.0):
        return E("InvalidTemperature")
    if not fgt(v, 0.0):
        return E("InvalidVolume")
    if not math.isfinite(aa) or not fge(b, 0.0):
        return E("InvalidParameter")
    v_eff = v - n * b
    if not fgt(v_eff, 0.0):
        return E("InvalidParameter")
    d = fdiv(n, v)
    return F(fdiv(n * GAS_R * temp, v_eff) - aa * d * d)


def m_iso_work(a, p):
    n, temp, v1, v2 = a
    if not fge(n, 0.0):
        return E("InvalidParameter")
    if not fge(temp, 0.0):
        return E("InvalidTemperature")
    if not fgt(v1, 0.0) or not fgt(v2, 0.0):
        return E("InvalidVolume")
    return F(n * GAS_R * temp * p.log(fdiv(v2, v1)))


def m_adiabatic(a, p):
    t1, v1, v2, g = a
    if not fgt(t1, 0.0):
        return E("InvalidTemperature")
    if not fgt(g, 1.0):
        return E("InvalidParameter")
    if not fgt(v1, 0.0) or not fgt(v2, 0.0):
        return E("InvalidVolume")
    return F(t1 * p.pow(fdiv(v1, v2), g - 1.0))


def m_z(a, p):
    pr, v, n, temp = a
    if not fgt(pr, 0.0):
        return E("InvalidPressure")
    if not fgt(v, 0.0):
        return E("InvalidVolume")
    if not fgt(temp, 0.0):
        return E("InvalidTemperature")
    if not math.isfinite(n):
        return E("InvalidParameter")
    d = n * GAS_R * temp
    if not fge(abs(d), TINY):
        return E("DivisionByZero")
    return F(fdiv(pr * v, d))


def _critical_ok(tc, pc):
    if not fgt(tc, 0.0):
        return "InvalidTemperature"
    if not fgt(pc, 0.0):
        return "InvalidPressure"
    return None


def rk_params(tc, pc, p):
    return (fdiv(0.42748 * GAS_R * GAS_R * p.pow(tc, 2.5), pc), fdiv(0.08664 * GAS_R * tc, pc))


def pr_params(tc, pc, om):
    return (fdiv(0.45724 * GAS_R * GAS_R * tc * tc, pc), fdiv(0.07780 * GAS_R * tc, pc),
            0.37464 + 1.54226 * om - 0.26992 * om * om)


def m_rk(a, p):
    temp, vm, tc, pc = a
    if not fgt(temp, 0.0):
        return E("InvalidTemperature")
    e = _critical_ok(tc, pc)
    if e:
        return E(e)
    if not fgt(vm, 0.0):
        return E("InvalidVolume")
    aa, b = rk_params(tc, pc, p)
    vm_b = vm - b
    if not fgt(vm_b, 0.0):
        return E("InvalidParameter")
    return F(fdiv(GAS_R * temp, vm_b) - fdiv(aa, fsqrt(temp) * vm * (vm + b)))


def m_pr(a, p):
    temp, vm, tc, pc, om = a
    if not fgt(temp, 0.0):
        return E("InvalidTemperature")
    e = _critical_ok(tc, pc)
    if e:
        return E(e)
    if not fgt(vm, 0.0):
        return E("InvalidVolume")
    if not math.isfinite(om):
        return E("InvalidParameter")
    aa, b, k = pr_params(tc, pc, om)
    vm_b = vm - b
    if not fgt(vm_b, 0.0):
        return E("InvalidParameter")
    al = 1.0 + k * (1.0 - fsqrt(fdiv(temp, tc)))
    d = vm * vm + 2.0 * b * vm - b * b
    if not fge(abs(d), TINY):
        return E("DivisionByZero")
    return F(fdiv(GAS_R * temp, vm_b) - fdiv(aa * (al * al), d))


def _pitzer_b(tr, p):
    return 0.083 - fdiv(0.422, p.pow(tr, 1.6)), 0.139 - fdiv(0.172, p.pow(tr, 4.2))


def m_pitzer_b(a, p):
    temp, tc, pc, om = a
    b0, b1 = _pitzer_b(fdiv(temp, tc), p)
    return F(fdiv((b0 + om * b1) * GAS_R * tc, pc))


def m_virial(a, p):
    temp, vm, b = a[:3]
    if not fgt(temp, 0.0):
        return E("InvalidTemperature")
    if not fgt(vm, 0.0):
        return E("InvalidVolume")
    z = 1.0 + fdiv(b, vm)
    if len(a) == 4:
        z = z + fdiv(a[3], vm * vm)
    return F(fdiv(z * GAS_R * temp, vm))


def m_z_pitzer(a, p):
    pr, tr, om = a
    if not fgt(tr, 0.0) or not fge(pr, 0.0) or not math.isfinite(om):
        return E("InvalidParameter")
    b0, b1 = _pitzer_b(tr, p)
    z = (1.0 + fdiv(b0 * pr, tr)) + om * fdiv(b1 * pr, tr)
    return F(z) if fgt(z, 0.0) else E("InvalidParameter")


def m_kays(a, p):
    ys, xs = a
    if len(ys) != len(xs) or not all(0.0 <= y <= 1.0 for y in ys) or not all(math.isfinite(x) for x in xs):
        return E("InvalidParameter")
    s = -0.0
    for y in ys:
        s = s + y
    if abs(s - 1.0) > 1e-10:
        return E("InvalidParameter")
    r = -0.0
    for y, x in zip(ys, xs):
        r = r + y * x
    return F(r)


def m_ig_ds(a, p):
    n, cv, t1, t2, v1, v2 = a
    if not fgt(n, 0.0):
        return E("InvalidParameter")
    if not fgt(cv, 0.0):
        return E("InvalidSpecificHeat")
    if not fgt(t1, 0.0) or not fgt(t2, 0.0):
        return E("InvalidTemperature")
    if not fgt(v1, 0.0) or not fgt(v2, 0.0):
        return E("InvalidVolume")
    return F(n * (cv * p.log(fdiv(t2, t1)) + GAS_R * p.log(fdiv(v2, v1))))


def m_iso_ds(a, p):
    n, v1, v2 = a
    if not fgt(n, 0.0):
        return E("InvalidParameter")
    if not fgt(v1, 0.0) or not fgt(v2, 0.0):
        return E("InvalidVolume")
    return F(n * GAS_R * p.log(fdiv(v2, v1)))


def m_heat_ds(a, p):
    q, temp = a
    if not math.isfinite(q):
        return E("InvalidParameter")
    if not fgt(temp, 0.0):
        return E("InvalidTemperature")
    return F(fdiv(q, temp))


def m_carnot(a, p):
    th, tc = a
    if not fgt(th, 0.0) or not fgt(tc, 0.0):
        return E("InvalidTemperature")
    if tc >= th:
        return E("InvalidParameter")
    return F(1.0 - fdiv(tc, th))


def m_cop(a, p):
    th, tc = a
    if not fgt(th, 0.0) or not fgt(tc, 0.0):
        return E("InvalidTemperature")
    dt = th - tc
    if not fgt(dt, 0.0):
        return E("InvalidParameter")
    return F(fdiv(tc, dt))


def m_mixing(a, p):
    n, xs = a
    if not fgt(n, 0.0) or not all(0.0 <= x <= 1.0 for x in xs):
        return E("InvalidParameter")
    fs = -0.0
    for x in xs:
        fs = fs + x
    if abs(fs - 1.0) > 1e-10:
        return E("InvalidParameter")
    s = -0.0
    for x in xs:
        if x > 0.0:
            s = s + x * p.log(x)
    return F(-n * GAS_R * s)


ATM = 101325.0  # state::ATM
F64_EPS = 2.220446049250313e-16


class Sub:
    """A SubstanceData from a row's eleven tokens (the name, then the ten f64 in declaration order)."""

    def __init__(self, a):
        (self.tpt, self.tpp, self.ct, self.cp, self.cd, self.lf, self.lv, self.mp, self.bp, self.mm) = a[1:11]

    def sublimes(self):
        return self.tpp > ATM


def _log_interp(t, t1, p1, t2, p2, p):
    span = fdiv(1.0, t2) - fdiv(1.0, t1)
    if abs(span) < F64_EPS:
        return p1
    frac = fdiv(fdiv(1.0, t) - fdiv(1.0, t1), span)
    return p1 * p.exp(frac * p.log(fdiv(p2, p1)))


def _cc_pressure(pr, tr, tt, lm, p):
    if not fgt(pr, 0.0):
        return "InvalidPressure", None
    if not fgt(tr, 0.0) or not fgt(tt, 0.0):
        return "InvalidTemperature", None
    if not fgt(lm, 0.0):
        return "InvalidParameter", None
    ex = -fdiv(lm, GAS_R) * (fdiv(1.0, tt) - fdiv(1.0, tr))
    v = pr * p.exp(ex)
    return ("InvalidParameter", None) if not math.isfinite(v) else (None, v)


def _vapor(x, t, p):
    if not (math.isfinite(t) and t >= x.tpt) or t > x.ct:
        return "UndefinedPhaseTransition", None
    if x.bp > x.tpt and x.bp < x.ct:
        if t <= x.bp:
            return None, _log_interp(t, x.tpt, x.tpp, x.bp, ATM, p)
        return None, _log_interp(t, x.bp, ATM, x.ct, x.cp, p)
    return None, _log_interp(t, x.tpt, x.tpp, x.ct, x.cp, p)


def _sublimation(x, t, p):
    if not fgt(t, 0.0) or t > x.tpt:
        return "UndefinedPhaseTransition", None
    if x.sublimes() and x.bp < x.tpt:
        return None, _log_interp(t, x.tpt, x.tpp, x.bp, ATM, p)
    return _cc_pressure(x.tpp, x.tpt, t, (x.lf + x.lv) * x.mm, p)


def _lh_sub(x, p):
    d = fdiv(1.0, x.bp) - fdiv(1.0, x.tpt)
    if x.sublimes() and x.bp < x.tpt and fgt(d, 0.0) and x.mm > 0.0:
        return fdiv(fdiv(GAS_R * p.log(fdiv(x.tpp, ATM)), d), x.mm)
    return x.lf + x.lv


def _res(r):
    e, v = r
    return E(e) if e else F(v)


def m_phase_at(a, p):
    x = Sub(a)
    t, pr = a[11], a[12]
    if not fgt(t, 0.0):
        return E("InvalidTemperature")
    if not fgt(pr, 0.0):
        return E("InvalidPressure")
    if t > x.ct and pr > x.cp:
        return ["e:Phase.SupercriticalFluid"]
    if t < x.tpt:
        e, ps = _sublimation(x, t, p)
        if e:
            return E(e)
        return ["e:Phase.Solid" if pr >= ps else "e:Phase.Gas"]
    if t <= x.ct:
        e, ps = _vapor(x, t, p)
        if e:
            return E(e)
        if pr >= ps:
            return ["e:Phase.Solid" if t <= x.mp else "e:Phase.Liquid"]
        return ["e:Phase.Gas"]
    return ["e:Phase.Gas"]


def m_heat_phase(a, p):
    mass, cp, ts, te = a[:4]
    x = Sub(a[4:])
    if not fgt(mass, 0.0):
        return E("InvalidParameter")
    if not fgt(cp, 0.0):
        return E("InvalidSpecificHeat")
    if not fgt(ts, 0.0) or not fgt(te, 0.0):
        return E("InvalidTemperature")
    lo, hi, sign = (ts, te, 1.0) if te >= ts else (te, ts, -1.0)
    q = 0.0
    q += mass * cp * (hi - lo)
    if x.sublimes():
        if lo < x.bp and hi >= x.bp:
            q += mass * _lh_sub(x, p)
    else:
        if lo < x.mp and hi >= x.mp:
            q += mass * x.lf
        if lo < x.bp and hi >= x.bp:
            q += mass * x.lv
    return F(sign * q)


def m_cc_slope(a, p):
    lh, temp, dv = a
    if not fgt(temp, 0.0):
        return E("InvalidTemperature")
    if not fge(abs(dv), TINY):
        return E("DivisionByZero")
    if not fgt(lh, 0.0):
        return E("InvalidParameter")
    return F(fdiv(lh, temp * dv))


# steam: the private tables, parsed from rust-old/src/steam.rs by this script (load_steam, from main), independently of
# the generator's parser.
STEAM = {"rows": [], "r4": []}
SAT_FIELDS = ["temperature", "pressure", "v_f", "v_g", "h_f", "h_fg", "h_g", "s_f", "s_g"]
R4_T_MIN, T_CRIT, P_CRIT = 273.15, 647.096, 22_064_000.0


def load_steam(root):
    src = (root / "rust-old/src/steam.rs").read_text()
    r4 = re.search(r"const R4_N: \[f64; 10\] = \[(.*?)\];", src, re.S).group(1)
    STEAM["r4"] = [rust_literal(v.strip()) for v in r4.split(",") if v.strip()]
    table = src[src.index("const SATURATED_TABLE"):]
    table = table[: table.index("\n];")]
    rows = []
    for body in re.findall(r"SaturatedEntry \{(.*?)\}", table, re.S):
        vals = dict((k, float(rust_literal(v))) for k, v in re.findall(r"(\w+):\s*([-0-9._e]+)", body))
        rows.append([vals[k] for k in SAT_FIELDS])
    STEAM["rows"] = rows
    tiers = re.search(r"const SUPERHEAT_PRESSURES: &\[f64\] = &\[(.*?)\];", src, re.S).group(1)
    STEAM["tiers"] = [float(rust_literal(v.strip())) for v in tiers.split(",") if v.strip()]
    sh = src[src.index("const SUPERHEATED_DATA"):]
    sh = sh[: sh.index("\n];")]
    fields = ["temperature", "pressure", "specific_volume", "specific_enthalpy", "specific_entropy"]
    sh_rows = []
    for body in re.findall(r"SuperheatedEntry \{(.*?)\}", sh, re.S):
        vals = dict((k, float(rust_literal(v))) for k, v in re.findall(r"(\w+):\s*([-0-9._e]+)", body))
        sh_rows.append([vals[k] for k in fields])
    STEAM["sh"] = sh_rows


def _r4p(t):
    n = STEAM["r4"]
    th = t + fdiv(n[8], t - n[9])
    a = th * th + n[0] * th + n[1]
    b = n[2] * th * th + n[3] * th + n[4]
    c = n[5] * th * th + n[6] * th + n[7]
    x = fdiv(2.0 * c, -b + fsqrt(b * b - 4.0 * a * c))
    return powidf2(x, 4) * 1e6


def _r4t(p):
    n = STEAM["r4"]
    beta = fsqrt(fsqrt(p * 1e-6))
    e = beta * beta + n[2] * beta + n[5]
    f = n[0] * beta * beta + n[3] * beta + n[6]
    g = n[1] * beta * beta + n[4] * beta + n[7]
    d = fdiv(2.0 * g, -f - fsqrt(f * f - 4.0 * e * g))
    s = n[9] + d
    return 0.5 * (s - fsqrt(s * s - 4.0 * (n[8] + n[9] * d)))


def m_sat_p(a, p):
    t = a[0]
    if not fge(t, R4_T_MIN) or t > T_CRIT:
        return E("SteamTableOutOfRange")
    return F(rust_min(_r4p(t), P_CRIT))


def m_sat_t(a, p):
    pr = a[0]
    if not fge(pr, _r4p(R4_T_MIN)) or pr > P_CRIT:
        return E("SteamTableOutOfRange")
    return F(rust_clamp(_r4t(pr), R4_T_MIN, T_CRIT))


def _lerp(a, b, t):
    return a + (b - a) * t  # hisab 1.4.0's calc::lerp


def _sat_interp(t, pr):
    rows = STEAM["rows"]
    hi = sum(1 for r in rows if r[0] <= t)
    hi = min(max(hi, 1), len(rows) - 1)
    a, b = rows[hi - 1], rows[hi]
    frac = fdiv(t - a[0], b[0] - a[0])
    za, zb = fdiv(a[3] * a[1], a[0]), fdiv(b[3] * b[1], b[0])
    vg = fdiv(_lerp(za, zb, frac) * t, pr)
    return [t, pr, _lerp(a[2], b[2], frac), vg] + [_lerp(a[k], b[k], frac) for k in range(4, 9)]


def m_sat_by_t(a, p):
    t = a[0]
    rows = STEAM["rows"]
    if not fge(t, rows[0][0]) or t > rows[-1][0]:
        return E("SteamTableOutOfRange")
    return [fbits(x) for x in _sat_interp(t, rust_min(_r4p(t), P_CRIT))]


def m_sat_by_p(a, p):
    pr = a[0]
    rows = STEAM["rows"]
    if not fge(pr, rows[0][1]) or pr > rows[-1][1]:
        return E("SteamTableOutOfRange")
    return [fbits(x) for x in _sat_interp(rust_clamp(_r4t(pr), rows[0][0], T_CRIT), pr)]


def _quality(x):
    if not (-1e-9 <= x <= 1.0 + 1e-9):
        return E("InvalidQuality")
    return F(rust_clamp(x, 0.0, 1.0))


def _q_from(a, lo_k, hi_k):
    v, e = a[0], a[1:10]
    d = e[hi_k] - e[lo_k]
    if not fge(abs(d), TINY):
        return E("DivisionByZero")
    return _quality(fdiv(v - e[lo_k], d))


def m_q_h(a, p):
    h, e = a[0], a[1:10]
    if not fge(abs(e[5]), TINY):
        return E("DivisionByZero")
    return _quality(fdiv(h - e[4], e[5]))


def m_wet(a, p):
    x, e = a[0], a[1:10]
    if not (0.0 <= x <= 1.0):
        return E("InvalidQuality")
    return [fbits(v) for v in (x, e[2] + x * (e[3] - e[2]), e[4] + x * e[5], e[7] + x * (e[8] - e[7]))]


def _sat_by_p_vals(pr):
    """saturated_by_pressure's nine values, or None for its refusal."""
    rows = STEAM["rows"]
    if not fge(pr, rows[0][1]) or pr > rows[-1][1]:
        return None
    return _sat_interp(rust_clamp(_r4t(pr), rows[0][0], T_CRIT), pr)


def _sh_lerp(a, b, t):
    frac = fdiv(t - a[0], b[0] - a[0])
    return [t, a[1], _lerp(a[2], b[2], frac), _lerp(a[3], b[3], frac), _lerp(a[4], b[4], frac)]


def _sh_tier(tier, t):
    rows = STEAM["sh"][tier * 8 : (tier + 1) * 8]
    first, last = rows[0], rows[-1]
    if t > last[0]:
        return None
    if t < first[0]:
        sat = _sat_by_p_vals(first[1])
        if sat is None or t < sat[0]:
            return None
        return _sh_lerp([sat[0], first[1], sat[3], sat[6], sat[8]], first, t)
    hi = min(max(sum(1 for r in rows if r[0] <= t), 1), len(rows) - 1)
    return _sh_lerp(rows[hi - 1], rows[hi], t)


def _sh_lookup(t, pr, p):
    """(error name, None) or (None, the five values)."""
    tiers = STEAM["tiers"]
    if not fgt(t, 0.0):
        return "InvalidTemperature", None
    if not fgt(pr, 0.0):
        return "InvalidPressure", None
    if pr < tiers[0] or pr > tiers[-1]:
        return "SuperheatOutOfRange", None
    hi = sum(1 for x in tiers if x < pr)
    if tiers[hi] <= pr:
        e = _sh_tier(hi, t)
        if e is None:
            return "SuperheatOutOfRange", None
        e[1] = pr
        return None, e
    a, b = _sh_tier(hi - 1, t), _sh_tier(hi, t)
    if a is None or b is None:
        return "SuperheatOutOfRange", None
    frac = fdiv(p.log(fdiv(pr, tiers[hi - 1])), p.log(fdiv(tiers[hi], tiers[hi - 1])))
    v = p.exp(_lerp(p.log(a[2]), p.log(b[2]), frac))
    return None, [t, pr, v, _lerp(a[3], b[3], frac), _lerp(a[4], b[4], frac)]


def _sh_out(r):
    e, v = r
    return E(e) if e else [fbits(x) for x in v]


def _sh_by_s(pr, s, p):
    """superheated_by_entropy: (error name, None) or (None, the five values)."""
    tiers, sh = STEAM["tiers"], STEAM["sh"]
    if not fgt(pr, 0.0):
        return "InvalidPressure", None
    if not math.isfinite(s):
        return "InvalidParameter", None
    if pr < tiers[0] or pr > tiers[-1]:
        return "SuperheatOutOfRange", None
    hi = sum(1 for x in tiers if x < pr)
    on_tier = tiers[hi] <= pr
    sat = _sat_by_p_vals(tiers[hi])
    if sat is None:
        return "SteamTableOutOfRange", None
    t_min = sat[0]
    t_max = sh[hi * 8 + 7][0]
    if not on_tier:
        t_max = rust_min(t_max, sh[(hi - 1) * 8 + 7][0])
    bps = [r[0] for r in sh[hi * 8 : hi * 8 + 8] if t_min < r[0] < t_max]
    if not on_tier:
        bps += [r[0] for r in sh[(hi - 1) * 8 : hi * 8] if t_min < r[0] < t_max]
    bps = sorted(bps + [t_min, t_max])

    def s_at(t):
        e, v = _sh_lookup(t, pr, p)
        return e, (v[4] if v else None)

    lo, hi_i = 0, len(bps) - 1
    e1, s_lo = s_at(bps[lo])
    if e1:
        return e1, None
    e2, s_hi = s_at(bps[hi_i])
    if e2:
        return e2, None
    if s < s_lo or s > s_hi:
        return "SuperheatOutOfRange", None
    while hi_i - lo > 1:
        mid = (lo + hi_i) // 2
        em, s_mid = s_at(bps[mid])
        if em:
            return em, None
        if s_mid <= s:
            lo, s_lo = mid, s_mid
        else:
            hi_i, s_hi = mid, s_mid
    ds = s_hi - s_lo
    tt = bps[lo] + fdiv(s - s_lo, ds) * (bps[hi_i] - bps[lo]) if ds > 0.0 else bps[lo]
    return _sh_lookup(rust_clamp(tt, t_min, t_max), pr, p)


def m_sh_by_s(a, p):
    return _sh_out(_sh_by_s(a[0], a[1], p))


def _cycle_toks(kind, states, procs, nums):
    toks = [f"e:CycleKind.{kind}", f"v:{len(states)}"]
    for st in states:
        toks += [fbits(x) for x in st]
    toks += [f"v:{len(procs)}"] + [f"e:ProcessKind.{k}" for k in procs]
    return toks + [fbits(x) for x in nums]



def _q_check(x):
    """validate_quality: (error name, None) or (None, the clamped quality)."""
    if not (-1e-9 <= x <= 1.0 + 1e-9):
        return "InvalidQuality", None
    return None, rust_clamp(x, 0.0, 1.0)


def _wet(x, e):
    """wet_steam_properties(x, e)'s (v, h, s) for a quality already in [0, 1]."""
    return e[2] + x * (e[3] - e[2]), e[4] + x * e[5], e[7] + x * (e[8] - e[7])


def m_rankine(a, p):
    pc, pb, tsh = a
    if not fgt(pc, 0.0) or not fgt(pb, 0.0):
        return E("InvalidPressure")
    if not fgt(pb, pc):
        return E("InvalidCycleParameter")
    sc = _sat_by_p_vals(pc)
    if sc is None:
        return E("SteamTableOutOfRange")
    t1, v1, h1, s1 = sc[0], sc[2], sc[4], sc[7]
    w_pump = v1 * (pb - pc)
    h2 = h1 + w_pump
    if tsh is not None:
        e, sh = _sh_lookup(tsh, pb, p)
        if e:
            return E(e)
        t3, v3, h3, s3 = sh[0], sh[2], sh[3], sh[4]
    else:
        sb = _sat_by_p_vals(pb)
        if sb is None:
            return E("SteamTableOutOfRange")
        t3, v3, h3, s3 = sb[0], sb[3], sb[6], sb[8]
    s4 = s3
    if s4 <= sc[8]:
        ds = sc[8] - sc[7]
        if not fge(abs(ds), TINY):
            return E("DivisionByZero")
        e, x = _q_check(fdiv(s4 - sc[7], ds))
        if e:
            return E(e)
        v4, h4, _ = _wet(x, sc)
        t4 = sc[0]
    else:
        e, sh4 = _sh_by_s(pc, s4, p)
        if e:
            return E(e)
        t4, v4, h4 = sh4[0], sh4[2], sh4[3]
    q_in, w_t, q_out = h3 - h2, h3 - h4, h4 - h1
    w_net = w_t - w_pump
    st = [[t1, pc, v1, s1, h1], [t1, pb, v1, s1, h2], [t3, pb, v3, s3, h3], [t4, pc, v4, s4, h4]]
    procs = ["Isentropic", "Isobaric", "Isentropic", "Isobaric"]
    return _cycle_toks("Rankine", st, procs, [q_in, q_out, w_net, fdiv(w_net, q_in), fdiv(w_pump, w_t), 0.0])


def m_refrig(a, p):
    pe, pc = a
    if not fgt(pe, 0.0) or not fgt(pc, 0.0):
        return E("InvalidPressure")
    if not fgt(pc, pe):
        return E("InvalidCycleParameter")
    se = _sat_by_p_vals(pe)
    if se is None:
        return E("SteamTableOutOfRange")
    sc = _sat_by_p_vals(pc)
    if sc is None:
        return E("SteamTableOutOfRange")
    t1, v1, h1, s1 = se[0], se[3], se[6], se[8]
    e, sh2 = _sh_by_s(pc, s1, p)
    if e:
        return E(e)
    t2, v2, h2 = sh2[0], sh2[2], sh2[3]
    t3, v3, h3, s3 = sc[0], sc[2], sc[4], sc[7]
    h4 = h3
    if not fge(abs(se[5]), TINY):
        return E("DivisionByZero")
    e, x = _q_check(fdiv(h4 - se[4], se[5]))
    if e:
        return E(e)
    v4, _, s4 = _wet(x, se)
    q_evap, q_cond, w = h1 - h4, h2 - h3, h2 - h1
    cop_r, cop_h = fdiv(q_evap, w), fdiv(q_cond, w)
    st = [[t1, pe, v1, s1, h1], [t2, pc, v2, s1, h2], [t3, pc, v3, s3, h3], [se[0], pe, v4, s4, h4]]
    procs = ["Isentropic", "Isobaric", "Isenthalpic", "Isobaric"]
    toks = _cycle_toks("Refrigeration", st, procs, [q_evap, q_cond, -w, cop_r, 1.0, 0.0])
    return toks + [fbits(cop_r), fbits(cop_h)]


def _parse_cycle(toks, at):
    """A CycleResult's tokens from toks[at]: (kind, state points, processes, its six numbers, the next index)."""
    kind = toks[at][len("e:CycleKind."):]
    n = int(toks[at + 1][2:])
    at += 2
    states = [[tok_f(t) for t in toks[at + 5 * k : at + 5 * k + 5]] for k in range(n)]
    at += 5 * n
    m = int(toks[at][2:])
    procs = [t[len("e:ProcessKind."):] for t in toks[at + 1 : at + 1 + m]]
    at += 1 + m
    return kind, states, procs, [tok_f(t) for t in toks[at : at + 6]], at + 6


def m_compare(ins, p):
    m, at, out = int(ins[0][2:]), 1, [ins[0]]
    for _ in range(m):
        label = ins[at]
        kind, _, _, nums, at = _parse_cycle(ins, at + 1)
        th, tl = tok_f(ins[at]), tok_f(ins[at + 1])
        at += 2
        if not fgt(tl, 0.0) or not fgt(th, tl):
            c = 0.0
        elif kind == "Refrigeration":
            c = fdiv(tl, th - tl)
        else:
            c = 1.0 - fdiv(tl, th)
        second = fdiv(nums[3], c) if c > 0.0 else 0.0
        out += [label] + [fbits(x) for x in (nums[3], nums[2], nums[4], c, second)]
    return out


def _diagram(ins, p, point):
    """A diagram's tokens: `v:count`, then each point's x, y from `point(a, b, process, frac, p)`."""
    _, st, procs, nums, at = _parse_cycle(ins, 0)
    pts = max(int(ins[at][2:]), 2)
    n, out = len(st), []
    for i in range(n):
        a, b = st[i], st[(i + 1) % n]
        pk = procs[i] if i < len(procs) else None
        for j in range(pts):
            out += [fbits(x) for x in point(a, b, pk, fdiv(float(j), float(pts - 1)), nums[5], p)]
    return [f"v:{n * pts}"] + out


def _ts_point(a, b, pk, frac, gamma, p):
    s = a[3] + frac * (b[3] - a[3])
    if gamma > 0.0 and pk in ("Isochoric", "Isobaric") and fgt(a[0], 0.0) and fgt(b[0], 0.0):
        return s, a[0] * p.pow(fdiv(b[0], a[0]), frac)
    return s, a[0] + frac * (b[0] - a[0])


def _pv_point(a, b, pk, frac, gamma, p):
    pa, pb, va, vb = a[1], b[1], a[2], b[2]
    positive = fgt(va, 0.0) and fgt(vb, 0.0) and fgt(pa, 0.0) and fgt(pb, 0.0)
    ln_v = p.log(fdiv(va, vb)) if positive else 0.0
    v = va + frac * (vb - va)
    linear = pa + frac * (pb - pa)
    if pk == "Isentropic" and abs(ln_v) > 1e-12:
        return v, pa * p.pow(fdiv(va, v), fdiv(p.log(fdiv(pb, pa)), ln_v))
    if pk == "Isochoric":
        return va, linear
    if pk == "Isobaric":
        return v, pa
    if pk == "Isothermal" and positive:
        return v, fdiv(pa * va, v)
    return v, linear


# Functions whose rows hold whole cycle results: their models read the raw tokens.
RAW_ARGS = {"cycle::compare_cycles", "cycle::cycle_ts_diagram", "cycle::cycle_pv_diagram"}

def _gas_ok(g, n):
    return fgt(g, 1.0) and fgt(n, 0.0)


def m_otto(a, p):
    t1, p1, r, q, g, n = a
    if not fgt(t1, 0.0):
        return E("InvalidTemperature")
    if not fgt(p1, 0.0):
        return E("InvalidPressure")
    if not fgt(r, 1.0):
        return E("InvalidCompressionRatio")
    if not _gas_ok(g, n) or not fgt(q, 0.0):
        return E("InvalidParameter")
    cv = fdiv(GAS_R, g - 1.0)
    cp = g * cv
    v1 = fdiv(n * GAS_R * t1, p1)
    v2 = fdiv(v1, r)
    t2 = t1 * p.pow(r, g - 1.0)
    p2 = fdiv(n * GAS_R * t2, v2)
    v3 = v2
    q_in = q * n
    t3 = t2 + fdiv(q_in, n * cv)
    p3 = fdiv(n * GAS_R * t3, v3)
    v4 = v1
    t4 = t3 * p.pow(fdiv(1.0, r), g - 1.0)
    p4 = fdiv(n * GAS_R * t4, v4)
    q_out = n * cv * (t4 - t1)
    w = q_in - q_out
    s3 = 0.0 + n * cv * p.log(fdiv(t3, t2))
    st = [[t1, p1, v1, 0.0, n * cp * t1], [t2, p2, v2, 0.0, n * cp * t2], [t3, p3, v3, s3, n * cp * t3],
          [t4, p4, v4, s3, n * cp * t4]]
    return _cycle_toks("Otto", st, ["Isentropic", "Isochoric", "Isentropic", "Isochoric"],
                       [q_in, q_out, w, fdiv(w, q_in), 0.0, g])


def m_diesel(a, p):
    t1, p1, r, rc, g, n = a
    if not fgt(t1, 0.0):
        return E("InvalidTemperature")
    if not fgt(p1, 0.0):
        return E("InvalidPressure")
    if not fgt(r, 1.0):
        return E("InvalidCompressionRatio")
    if not fgt(rc, 1.0):
        return E("InvalidCutoffRatio")
    if rc >= r:
        return E("InvalidCycleParameter")
    if not _gas_ok(g, n):
        return E("InvalidParameter")
    cv = fdiv(GAS_R, g - 1.0)
    cp = g * cv
    v1 = fdiv(n * GAS_R * t1, p1)
    v2 = fdiv(v1, r)
    t2 = t1 * p.pow(r, g - 1.0)
    p2 = fdiv(n * GAS_R * t2, v2)
    p3, v3, t3 = p2, v2 * rc, t2 * rc
    v4 = v1
    t4 = t3 * p.pow(fdiv(v3, v4), g - 1.0)
    p4 = fdiv(n * GAS_R * t4, v4)
    q_in = n * cp * (t3 - t2)
    q_out = n * cv * (t4 - t1)
    w = q_in - q_out
    s3 = 0.0 + n * cp * p.log(fdiv(t3, t2))
    st = [[t1, p1, v1, 0.0, n * cp * t1], [t2, p2, v2, 0.0, n * cp * t2], [t3, p3, v3, s3, n * cp * t3],
          [t4, p4, v4, s3, n * cp * t4]]
    return _cycle_toks("Diesel", st, ["Isentropic", "Isobaric", "Isentropic", "Isochoric"],
                       [q_in, q_out, w, fdiv(w, q_in), 0.0, g])


def m_brayton(a, p):
    t1, p1, rp, t3, g, n = a
    if not fgt(t1, 0.0):
        return E("InvalidTemperature")
    if not fgt(p1, 0.0):
        return E("InvalidPressure")
    if not fgt(t3, 0.0):
        return E("InvalidTemperature")
    if not fgt(rp, 1.0):
        return E("InvalidPressureRatio")
    if not _gas_ok(g, n):
        return E("InvalidParameter")
    cv = fdiv(GAS_R, g - 1.0)
    cp = g * cv
    ex = fdiv(g - 1.0, g)
    v1 = fdiv(n * GAS_R * t1, p1)
    p2 = p1 * rp
    t2 = t1 * p.pow(rp, ex)
    v2 = fdiv(n * GAS_R * t2, p2)
    if not fgt(t3, t2):
        return E("InvalidCycleParameter")
    p3 = p2
    v3 = fdiv(n * GAS_R * t3, p3)
    p4 = p1
    t4 = t3 * p.pow(fdiv(1.0, rp), ex)
    v4 = fdiv(n * GAS_R * t4, p4)
    q_in = n * cp * (t3 - t2)
    q_out = n * cp * (t4 - t1)
    w = q_in - q_out
    bwr = fdiv(n * cp * (t2 - t1), n * cp * (t3 - t4))
    s3 = 0.0 + n * cp * p.log(fdiv(t3, t2))
    st = [[t1, p1, v1, 0.0, n * cp * t1], [t2, p2, v2, 0.0, n * cp * t2], [t3, p3, v3, s3, n * cp * t3],
          [t4, p4, v4, s3, n * cp * t4]]
    return _cycle_toks("Brayton", st, ["Isentropic", "Isobaric", "Isentropic", "Isobaric"],
                       [q_in, q_out, w, fdiv(w, q_in), bwr, g])


def m_adiabatic_p(a, p):
    t1, p1, p2, g = a
    if not fgt(t1, 0.0):
        return E("InvalidTemperature")
    if not fgt(p1, 0.0) or not fgt(p2, 0.0):
        return E("InvalidPressure")
    if not fgt(g, 1.0):
        return E("InvalidParameter")
    return F(t1 * p.pow(fdiv(p2, p1), fdiv(g - 1.0, g)))


def unmodelled(fn, ins):
    """Rows of a modelled module that other checks own: the table presets (no inputs; check_tables) and the JSON
    documents (serde_json's text and verdicts, which the harness compares byte for byte; ADR 0005)."""
    return not ins or fn.endswith("::to_json") or fn.endswith("::from_json")


MODELS = {
    "cycle": {
        "cycle::isochoric_heat": lambda a, p: F(a[0] * a[1] * (a[3] - a[2])),
        "cycle::isobaric_heat": lambda a, p: F(a[0] * a[1] * (a[3] - a[2])),
        "cycle::adiabatic_work": lambda a, p: F(a[0] * a[1] * (a[2] - a[3])),
        "cycle::adiabatic_temp_from_pressure": m_adiabatic_p,
        "cycle::otto_cycle": m_otto,
        "cycle::diesel_cycle": m_diesel,
        "cycle::brayton_cycle": m_brayton,
        "cycle::rankine_cycle": m_rankine,
        "cycle::refrigeration_cycle": m_refrig,
        "cycle::heat_pump_cop": lambda a, p: F(a[0] + 1.0),
        "cycle::compare_cycles": m_compare,
        "cycle::cycle_ts_diagram": lambda a, p: _diagram(a, p, _ts_point),
        "cycle::cycle_pv_diagram": lambda a, p: _diagram(a, p, _pv_point),
    },
    "steam": {
        "steam::saturation_pressure": m_sat_p,
        "steam::saturation_temperature": m_sat_t,
        "steam::saturated_by_temperature": m_sat_by_t,
        "steam::saturated_by_pressure": m_sat_by_p,
        "steam::quality_from_volume": lambda a, p: _q_from(a, 2, 3),
        "steam::quality_from_enthalpy": m_q_h,
        "steam::quality_from_entropy": lambda a, p: _q_from(a, 7, 8),
        "steam::wet_steam_properties": m_wet,
        "steam::superheated_lookup": lambda a, p: _sh_out(_sh_lookup(a[0], a[1], p)),
        "steam::superheated_by_entropy": m_sh_by_s,
    },
    "phase": {
        "phase::SubstanceData::sublimes_at_atmospheric": lambda a, p: [f"i:{int(Sub(a).sublimes())}"],
        "phase::SubstanceData::latent_heat_sublimation": lambda a, p: F(_lh_sub(Sub(a), p)),
        "phase::SubstanceData::vapor_pressure": lambda a, p: _res(_vapor(Sub(a), a[11], p)),
        "phase::SubstanceData::sublimation_pressure": lambda a, p: _res(_sublimation(Sub(a), a[11], p)),
        "phase::SubstanceData::phase_at": m_phase_at,
        "phase::heat_of_fusion": lambda a, p: F(a[11] * Sub(a).lf),
        "phase::heat_of_vaporization": lambda a, p: F(a[11] * Sub(a).lv),
        "phase::heat_for_phase_change": m_heat_phase,
        "phase::clausius_clapeyron_slope": m_cc_slope,
        "phase::clausius_clapeyron_pressure": lambda a, p: _res(_cc_pressure(*a, p)),
    },
    "entropy": {
        "entropy::ideal_gas_entropy_change": m_ig_ds,
        "entropy::isothermal_entropy_change": m_iso_ds,
        "entropy::heat_transfer_entropy": m_heat_ds,
        "entropy::carnot_efficiency": m_carnot,
        "entropy::carnot_cop_refrigeration": m_cop,
        "entropy::helmholtz": lambda a, p: F(a[0] - a[1] * a[2]),
        "entropy::gibbs": lambda a, p: F(a[0] - a[1] * a[2]),
        "entropy::is_spontaneous": lambda a, p: [f"i:{int(a[0] >= 0.0)}"],
        "entropy::entropy_of_mixing": m_mixing,
    },
    "state": {
        "state::ideal_gas_pressure": m_ig_pressure,
        "state::ideal_gas_volume": m_ig_volume,
        "state::ideal_gas_temperature": m_ig_temperature,
        "state::van_der_waals_pressure": m_vdw,
        "state::isothermal_work": m_iso_work,
        "state::isobaric_work": lambda a, p: F(a[0] * (a[2] - a[1])),
        "state::adiabatic_temperature": m_adiabatic,
        "state::compressibility_factor": m_z,
        "state::reduced_temperature": lambda a, p: F(fdiv(a[0], a[1])),
        "state::reduced_pressure": lambda a, p: F(fdiv(a[0], a[1])),
        "state::redlich_kwong_params": lambda a, p: [fbits(x) for x in rk_params(a[0], a[1], p)],
        "state::peng_robinson_params": lambda a, p: [fbits(x) for x in pr_params(*a)],
        "state::redlich_kwong_pressure": m_rk,
        "state::peng_robinson_pressure": m_pr,
        "state::pitzer_second_virial": m_pitzer_b,
        "state::virial_pressure_2nd": m_virial,
        "state::virial_pressure_3rd": m_virial,
        "state::compressibility_pitzer": m_z_pitzer,
        "state::kays_rule_tc": m_kays,
        "state::kays_rule_pc": m_kays,
        "state::kays_rule_omega": m_kays,
    },
    "transfer": {
        "transfer::conduction": m_conduction,
        "transfer::convection": lambda a, p: F(a[0] * a[1] * (a[2] - a[3])),
        "transfer::radiation": m_radiation,
        "transfer::thermal_resistance_conduction": m_r_cond,
        "transfer::thermal_resistance_convection": m_r_conv,
        "transfer::thermal_resistance_series": m_series,
        "transfer::thermal_resistance_parallel": m_parallel,
        "transfer::heat_stored": lambda a, p: F(a[0] * a[1] * a[2]),
        "transfer::thermal_diffusivity": m_diffusivity,
        "transfer::biot_number": m_biot,
        "transfer::lumped_capacitance": lambda a, p: F(a[1] + (a[0] - a[1]) * p.exp(fdiv(-a[2], a[3]))),
        "transfer::reynolds_number": m_reynolds,
        "transfer::prandtl_number": m_prandtl,
        "transfer::nusselt_number": m_nusselt,
        "transfer::nusselt_dittus_boelter": m_dittus,
        "transfer::nusselt_natural_vertical": m_natural,
        "transfer::boundary_layer_thickness": m_blasius,
        "transfer::thermal_boundary_layer": m_thermal_bl,
        "transfer::fin_parameter": m_fin_parameter,
        "transfer::fin_rectangular_heat": m_fin_heat,
        "transfer::fin_efficiency_rectangular": m_fin_efficiency,
        "transfer::fin_effectiveness": m_fin_effectiveness,
        "transfer::lmtd_parallel": lambda a, p: _lmtd(a[0] - a[2], a[1] - a[3], p),
        "transfer::lmtd_counter": lambda a, p: _lmtd(a[0] - a[3], a[1] - a[2], p),
        "transfer::heat_exchanger_lmtd": lambda a, p: F(a[0] * a[1] * a[2]),
        "transfer::ntu": m_ntu,
        "transfer::effectiveness_parallel": m_eff_parallel,
        "transfer::effectiveness_counter": m_eff_counter,
        "transfer::heat_exchanger_ntu": lambda a, p: F(a[0] * a[1] * (a[2] - a[3])),
        "transfer::view_factor_parallel_plates": m_vf_parallel,
        "transfer::view_factor_perpendicular_plates": m_vf_perpendicular,
        "transfer::view_factor_coaxial_disks": m_vf_disks,
    },
}


def model_args(ins):
    """A row's inputs as floats; a `v:n` sequence is the list of its n elements. A row that is one sequence is that
    list itself (the resistance lists); one with two is the pair of lists (Kay's rules)."""
    out, i = [], 0
    while i < len(ins):
        if ins[i].startswith("v:"):
            n = int(ins[i][2:])
            out.append([tok_f(t) for t in ins[i + 1 : i + 1 + n]])
            i += 1 + n
        elif ins[i] == "n:":
            out.append(None)  # Option::None (rankine_cycle's superheat)
            i += 1
        elif ins[i].startswith(("s:", "e:")):
            out.append(ins[i])  # a name (a substance's) or an enum value, as its token
            i += 1
        else:
            out.append(tok_f(ins[i]))
            i += 1
    return out[0] if len(out) == 1 and isinstance(out[0], list) else out


def check_models(tables, rep, exact):
    """Each modelled module's rows, recomputed with glibc's primitives: the outputs must be the golden row's, token for
    token (NaN matching any NaN). Off the corpus's libm platform only rows with no transcendental are checked."""
    for mod, models in MODELS.items():
        for n, fn, cls, ins, outs in tables.get(f"{mod}.tsv", []):
            if unmodelled(fn, ins):
                continue
            model = models.get(fn)
            if model is None:
                rep.fail(f"{mod}.tsv:{n}", f"{fn} has no model")
                continue
            if cls == "C" and not exact:
                continue
            got = model(ins if fn in RAW_ARGS else model_args(ins), Prim(False))
            if len(got) == len(outs) and all(
                g == w or (g.startswith("f:") and w.startswith("f:") and same(tok_f(g), tok_f(w))) for g, w in zip(got, outs)
            ):
                rep.ok()
            else:
                rep.fail(f"{mod}.tsv:{n}", f"{fn} {toks_text(ins)}: model {got}, golden {outs}")


def model_corrections(mod, rows):
    """A modelled module's class-C rows whose correctly rounded outputs differ from the oracle's."""
    out = []
    models = MODELS[mod]
    for n, fn, cls, ins, outs in rows:
        if cls != "C":
            continue
        p = Prim(True)
        got = models[fn](ins if fn in RAW_ARGS else model_args(ins), p)
        if [*got] == [*outs] or not p.calls:
            continue
        out.append((fn, toks_text(ins), " ".join(got), "; ".join([f"oracle {toks_text(outs)}"] + p.calls)))
    return out


# ── fmtf: the f64 text, checked by its properties ──

def unpct(text):
    """Decode format 1's `s:` encoding."""
    out = bytearray()
    k = 0
    while k < len(text):
        if text[k] == "%":
            out.append(int(text[k + 1:k + 3], 16))
            k += 3
        else:
            out.append(ord(text[k]))
            k += 1
    return out.decode()


def sig_digits(text):
    """The significant digits of a decimal text, without sign, point, exponent or leading/trailing zeros."""
    mant = text.lstrip("-").split("e")[0].replace(".", "")
    return mant.strip("0") or "0"


def dec_exp(x):
    """The decimal exponent of a finite nonzero x's shortest digits: x = d.ddd × 10^e."""
    r = repr(abs(x))
    if "e" in r:
        m, e = r.split("e")
        return int(e) + (len(m.split(".")[0]) - 1)
    ip, _, fp = r.partition(".")
    if ip != "0":
        return len(ip) - 1
    return -(len(fp) - len(fp.lstrip("0")) + 1)


def check_fmtf(rows, rep):
    """Each fmtf row's text reads back to its bits, has the shortest digit count (Python's repr), is a closest
    shortest decimal (Rust's Display may break an exact tie upward; serde_json's must equal repr's, ties to even), and
    has its form: Display always positional; serde_json positional while the decimal exponent is in -5..=15."""
    from fractions import Fraction
    for n, fn, _cls, ins, outs in rows:
        where = f"fmtf.tsv:{n}"
        if fn not in ("fmtf::display", "fmtf::json") or len(ins) != 1 or len(outs) != 1 or not outs[0].startswith("s:"):
            rep.fail(where, f"unexpected row {fn} {ins} {outs}")
            continue
        x = tok_f(ins[0])
        text = unpct(outs[0][2:])
        if math.isnan(x) or math.isinf(x):
            want = "null" if fn == "fmtf::json" else ("NaN" if math.isnan(x) else ("inf" if x > 0 else "-inf"))
            if text != want:
                rep.fail(where, f"{text!r} for a non-finite value, expected {want!r}")
            else:
                rep.ok()
            continue
        try:
            back = float(text)
        except ValueError:
            rep.fail(where, f"{text!r} is not a number")
            continue
        if bits(back) != bits(x):
            rep.fail(where, f"{text!r} reads back as {back!r}, not {x!r}")
            continue
        if x == 0.0:
            want = ("-" if math.copysign(1.0, x) < 0 else "") + ("0.0" if fn == "fmtf::json" else "0")
            if text != want:
                rep.fail(where, f"{text!r} for a zero, expected {want!r}")
            else:
                rep.ok()
            continue
        if len(sig_digits(text)) != len(sig_digits(repr(x))):
            rep.fail(where, f"{text!r} is not shortest (repr {x!r})")
            continue
        exact = Fraction(x)
        if abs(Fraction(text) - exact) > abs(Fraction(repr(x)) - exact):
            rep.fail(where, f"{text!r} is farther from the value than {x!r}")
            continue
        e = dec_exp(x)
        if fn == "fmtf::display":
            if "e" in text or text.endswith(".0"):
                rep.fail(where, f"{text!r}: Display is positional, with no '.0'")
                continue
        else:
            if sig_digits(text) != sig_digits(repr(x)):
                rep.fail(where, f"{text!r}: serde_json's digits are repr's (ties to even), {x!r}")
                continue
            positional = "e" not in text
            if positional != (-5 <= e <= 15):
                rep.fail(where, f"{text!r}: the form for decimal exponent {e}")
                continue
            if positional and "." not in text:
                rep.fail(where, f"{text!r}: serde_json writes '.0' on an integer")
                continue
            if not positional and not re.fullmatch(r"-?[0-9](\.[0-9]+)?e[+-][0-9]+", text):
                rep.fail(where, f"{text!r}: serde_json's exponent form is d.ddde±x")
                continue
        rep.ok()


def main():
    root = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    rep = Report()
    golden_dir = root / "tests/golden"
    files = sorted(p for p in golden_dir.glob("*.tsv") if not p.name.endswith(".cr.tsv"))
    listed = [
        ln for ln in (golden_dir / "modules.txt").read_text().splitlines() if ln and not ln.startswith("#")
    ]
    if listed != [p.name[: -len(".tsv")] for p in files]:
        rep.fail("modules.txt", f"lists {listed}, the directory has {[p.name for p in files]}")
    tables = {p.name: read_tsv(p, rep) for p in files}
    load_steam(root)
    golden = check_constants(root, tables.get("constants.tsv", []), rep)
    check_constants_cyr(root, golden, rep)
    check_tables(root, tables, rep)
    counts = check_math(tables.get("math.tsv", []), rep)
    check_fmtf(tables.get("fmtf.tsv", []), rep)
    exact = recorded_platform((golden_dir / "math.tsv").read_text()) == running_platform()
    check_models(tables, rep, exact)
    per = check_corrections(root, tables, rep, "--write-corrections" in sys.argv[2:])
    total = sum(len(v) for v in tables.values())
    print(f"files: {len(files)}; rows: {total}; math classes: {dict(sorted(counts.items()))}")
    print(f"ADR 0004 corrected oracle (MPFR {MPFR.VERSION}): rows glibc misrounds: "
          f"{sum(per.values())} {dict(sorted(per.items()))}")
    print(f"independently recomputed: {rep.checked} (constants, src/constants.cyr, table presets, math classes A and B,"
          " fmtf's text by its properties, modelled modules)")
    print(f"FAILURES: {rep.failures}")
    return 1 if rep.failures else 0


if __name__ == "__main__":
    sys.exit(main())
