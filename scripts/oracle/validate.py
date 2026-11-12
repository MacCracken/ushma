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
- `math.tsv`, classes A and B: sqrt, abs, ceil, min, max, clamp, is_finite, is_nan and powi are
  recomputed in IEEE double arithmetic (Python floats are binary64; powi follows
  `__powidf2`'s order). Class C rows (glibc's transcendentals) are checked against MPFR by the
  corrected-oracle tooling (ADR 0004, M0 part 2c), not here.

Exit status 1 on any failure; the last line is `FAILURES: <n>`.
"""

import math
import pathlib
import re
import struct
import sys

TOKEN = re.compile(r"^(f:[0-9a-f]{16}|i:-?[0-9]+|e:[A-Za-z0-9_]+\.[A-Za-z0-9_]+|s:\S*|n:|E:[A-Za-z0-9_]+|v:[0-9]+)$")


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
    counts = check_math(tables.get("math.tsv", []), rep)
    total = sum(len(v) for v in tables.values())
    print(f"files: {len(files)}; rows: {total}; math classes: {dict(sorted(counts.items()))}")
    print(f"independently recomputed: {rep.checked} (constants, src/constants.cyr, math classes A and B)")
    print(f"FAILURES: {rep.failures}")
    return 1 if rep.failures else 0


if __name__ == "__main__":
    sys.exit(main())
