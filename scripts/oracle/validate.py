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
    """tests/golden/math.cr.tsv against the MPFR recomputation; each listed call against MPFR and glibc."""
    gold = root / "tests/golden"
    golden_text = (gold / "math.tsv").read_text()
    exact = recorded_platform(golden_text) == running_platform()
    want = corrections(tables.get("math.tsv", []))
    path = gold / f"math{CR_SUFFIX}"
    if write:
        if not exact:
            rep.fail("--write-corrections", f"this machine ({running_platform()}) is not the corpus's libm platform "
                     f"({recorded_platform(golden_text)}); nothing written")
        else:
            path.write_text(cr_text("math", golden_text, want))
            print(f"wrote tests/golden/math{CR_SUFFIX}: {len(want)} rows")
    for stale in sorted(gold.glob("*" + CR_SUFFIX)):
        if stale.name != f"math{CR_SUFFIX}":
            rep.fail(stale.name, "no module but math has corrections")
    if not path.exists():
        if want:
            rep.fail(path.name, f"missing; {len(want)} rows need a correction (run with --write-corrections)")
        return {}
    text = path.read_text()
    lines = text.split("\n")
    if lines[-1] != "":
        rep.fail(path.name, "no final newline")
    body = [ln for ln in lines[5:] if ln and not ln.startswith("#")]
    golden_keys = {(fn, toks_text(ins)): toks_text(outs) for n, fn, cls, ins, outs in tables.get("math.tsv", [])}
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
            ok = same(c, mpfr_call(name, args)) and not same(g, c) and outs == fbits(c)
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
    if exact and text != cr_text("math", golden_text, want):
        rep.fail(path.name, "differs from the corrected oracle (regenerate with --write-corrections)")
    per = {}
    for fn, *_ in want:
        per[fn] = per.get(fn, 0) + 1
    return per


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
    golden = check_constants(root, tables.get("constants.tsv", []), rep)
    check_constants_cyr(root, golden, rep)
    check_tables(root, tables, rep)
    counts = check_math(tables.get("math.tsv", []), rep)
    check_fmtf(tables.get("fmtf.tsv", []), rep)
    per = check_corrections(root, tables, rep, "--write-corrections" in sys.argv[2:])
    total = sum(len(v) for v in tables.values())
    print(f"files: {len(files)}; rows: {total}; math classes: {dict(sorted(counts.items()))}")
    print(f"ADR 0004 corrected oracle (MPFR {MPFR.VERSION}): math rows glibc misrounds: "
          f"{sum(per.values())} {dict(sorted(per.items()))}")
    print(f"independently recomputed: {rep.checked} (constants, src/constants.cyr, table presets, math classes A and B,"
          " fmtf's text by its properties)")
    print(f"FAILURES: {rep.failures}")
    return 1 if rep.failures else 0


if __name__ == "__main__":
    sys.exit(main())
