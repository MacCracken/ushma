"""Verify src/mathx.cyr's correctly rounded exp, expm1, ln, pow, atan2 and tanh against MPFR (ADR 0004), on every
target.

Vendored from kimiya's scripts/oracle/cr_verify.py (kimiya commit 5b14c85): its argument families for the five
vendored primitives (exp, expm1, ln, pow, atan2), with the published hard cases, unchanged; and a family of ushma's own
for tanh (M0 part 2c-3). For each primitive it generates seeded arguments (deterministic: random.Random with fixed
seeds), builds a probe program that calls the Cyrius function on every argument, runs it on each target, and compares
every result bit for bit with the correctly rounded MPFR value (scripts/oracle/mpfr_ref.py). Each primitive gets:
- about 10^6 arguments: binade-uniform over the domain plus the dense regions where the code changes path or the
  rounding is delicate (thresholds +- thousands of ulps, reduction boundaries, tiny-argument ends, table edges), every
  special and exact case, and the known hard-to-round arguments (Lefevre–Muller for exp, ln and atan; Hubrecht,
  Jeannerod and Zimmermann for pow);
- the accurate path forced (the _us_cr_<f>_acc entry point) on 100,000 of them plus the special/hard section;
- the fast path checked against the full result wherever its test passes.
tanh's family: binade-uniform |x| < 64 (both signs), typical ranges and the fin range mL in [0, 50], 3,000 ulps either
side of the tiny threshold 2^-27 and of TANH_ONE, the expm1 branch point 2a = 2^-6, and 2a at the exp core's table
boundaries K*ln2/512 and (K + 1/2)*ln2/512.
It also checks x^2.5 (ushma's Redlich–Kwong Tc^2.5) against MPFR over every binade, and the special values of pow,
atan2 and tanh against glibc (ctypes, libm.so.6) as well as MPFR.

    python3 -I -B scripts/oracle/cr_verify.py                      # every function, ELF, PE (wine), aarch64 (qemu)
    python3 -I -B scripts/oracle/cr_verify.py --only tanh,exp      # a subset
    python3 -I -B scripts/oracle/cr_verify.py --targets elf        # one target
    python3 -I -B scripts/oracle/cr_verify.py --work DIR           # probe sources and binaries in DIR
                                                                   # (default build/cr-verify, gitignored)
Exit status 1 on any mismatch or cross-target difference, 2 on a build or run failure.
"""
import importlib.util
import math
import random
import struct
import subprocess
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True   # no __pycache__ in scripts/oracle/
_here = Path(__file__).resolve().parent


def _load(name):
    """A sibling script by path (`-I` keeps the script's directory off sys.path)."""
    spec = importlib.util.spec_from_file_location(name, _here / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mpfr_ref = _load("mpfr_ref")

FUNCS = [  # (id, name, full entry, accurate entry, MPFR reference); ids as kimiya's, tanh ushma's own
    (0, "exp", "_us_cr_exp", "_us_cr_exp_acc", mpfr_ref.exp),
    (2, "expm1", "_us_cr_expm1", "_us_cr_expm1_acc", mpfr_ref.expm1),
    (3, "ln", "_us_cr_ln", "_us_cr_ln_acc", mpfr_ref.log),
    (5, "pow", "_us_cr_pow", "_us_cr_pow_acc", mpfr_ref.pow),
    (8, "atan2", "_us_cr_atan2", "_us_cr_atan2_acc", mpfr_ref.atan2),
    (40, "tanh", "_us_cr_tanh", "_us_cr_tanh_acc", mpfr_ref.tanh),
]
ARITY = {5: 2, 8: 2}   # argument words; 1 elsewhere
N_UNIFORM = 700_000
N_ACC = 100_000

fb = mpfr_ref.from_bits
bits = mpfr_ref.bits
INF = float("inf")
NAN = float("nan")


def ulps(x, k):
    """x moved k ulps along the ordered bit line (k may be negative), staying finite."""
    b = bits(x)
    o = b if b < (1 << 63) else -(b & ((1 << 63) - 1))
    o += k
    nb = o if o >= 0 else ((-o) | (1 << 63))
    y = fb(nb)
    return y if math.isfinite(y) else x


def binade_uniform(rng, n, emin, emax, signs=(1, -1)):
    """n doubles with a uniformly random binade 2^e (emin <= e <= emax; the subnormal binades too) and a
    uniformly random significand in it."""
    out = []
    for _ in range(n):
        e = rng.randint(emin, emax)
        if e < -1022:   # a subnormal binade: the leading bit at 2^e
            k = e + 1074
            x = fb((1 << k) | rng.randrange(0, 1 << k))
        else:
            x = fb(((e + 1023) << 52) | rng.randrange(0, 1 << 52))
        out.append(rng.choice(signs) * x)
    return out


def around(x, k, rng=None, n=None):
    """x and k ulps either side (all of them), or n random offsets within +-k."""
    if n is None:
        return [ulps(x, d) for d in range(-k, k + 1)]
    return [ulps(x, rng.randint(-k, k)) for _ in range(n)]


SPECIALS = [0.0, -0.0, INF, -INF, NAN, fb(1), -fb(1), fb(0x0010000000000000), -fb(0x0010000000000000),
            fb(0x7FEFFFFFFFFFFFFF), -fb(0x7FEFFFFFFFFFFFFF), 1.0, -1.0, 2.0, 0.5, 10.0, 0.1, -0.1,
            fb(0x000FFFFFFFFFFFFF), 2.0 ** -54, -2.0 ** -54, 2.0 ** -53, -2.0 ** -53, 2.0 ** -30, -2.0 ** -30,
            2.0 ** -6, -2.0 ** -6, 708.0, -708.0, 307.0, -307.0]


def _b2f(mant, e, neg=False):
    ip, fp = mant.split(".")
    from fractions import Fraction
    v = Fraction(int(ip + fp, 2), 2 ** len(fp)) * Fraction(2) ** e
    x = float(v)
    return -x if neg else x


# Lefèvre–Muller 2003 worst-case database (Tables 4 and 5): exp and ln hardest-to-round arguments.
HARD_EXP = [
    _b2f("1.1110110100110001100011101111101101100010011111101010", -27, True),
    _b2f("1.0100000000000000000000000000000000000000000000110010", -46, True),
    _b2f("1.0000000000000000000000000000000000000000000000000001", -51, True),
    _b2f("1.1111111111111111111111111111111111111111111111111111", -53),
    _b2f("1.1111111111111111111111111111111111111111111110000000", -46),
    _b2f("1.0001111111111111111111111111111111111111111110101111", -45),
    _b2f("110.00001111010100101111001101111010111011001111110100", 0),
]
HARD_LN = [
    _b2f("1.1110101001110001110110000101110011101110000000100000", -509),
    _b2f("1.1001010001110110111000110000010011001101011111000111", -384),
    _b2f("1.0010011011101001110001001101001100100111100101100000", -232),
    _b2f("1.0110000100111001010101011101110010000000001011111000", -35),
    _b2f("1.0110001010101000100001100001001101100010100110110110", 678),
]
LN2 = 0.6931471805599453
LN10 = 2.302585092994046


def args_exp(rng, lim):
    ovf, unf = lim["EXP_OVF"], lim["EXP_UNF"]
    red = []
    for _ in range(30_000):  # reduction boundaries: (K + 1/2)·ln2/512 (K's rounding flips), K·ln2/512 (r ≈ 0)
        k = rng.randint(-550_000, 524_000)
        red.append(ulps((k + 0.5) * LN2 / 512, rng.randint(-3, 3)))
        red.append(ulps(k * LN2 / 512, rng.randint(-3, 3)))
    return [
        ("special+hard", list(SPECIALS) + HARD_EXP + [-h for h in HARD_EXP]),
        ("binade-uniform", binade_uniform(rng, N_UNIFORM, -1074, 9)),
        ("typical [-700,700]", [rng.uniform(-700.0, 700.0) for _ in range(50_000)]),
        ("[-1,1]", [rng.uniform(-1.0, 1.0) for _ in range(50_000)]),
        ("subnormal results", [rng.uniform(unf - 0.01, -708.0) for _ in range(120_000)]),
        ("overflow end", [rng.uniform(700.0, ovf + 0.01) for _ in range(40_000)]),
        ("thresholds+-ulps", around(ovf, 2000) + around(unf, 2000) + around(708.0, 500) + around(-708.0, 500)
         + around(-708.3964185322641, 2000)),                               # the last: results near DBL_MIN
        ("tiny |x|", [math.copysign(2.0 ** rng.uniform(-60, -25), rng.random() - 0.5) for _ in range(60_000)]
         + around(2.0 ** -54, 300) + around(-2.0 ** -54, 300) + around(2.0 ** -30, 300) + around(-2.0 ** -30, 300)),
        ("reduction edges", red),
    ]


def args_expm1(rng, lim):
    ovf, m1 = lim["EXPM1_OVF"], lim["EXPM1_M1"]
    red = []
    for _ in range(40_000):     # ± multiples of ln 2 / 512 and of ln 2 (j = 0, r ≈ 0; E changes)
        k = rng.randint(-30_000, 30_000)
        red.append(ulps(k * LN2 / 512, rng.randint(-3, 3)))
        red.append(ulps(rng.randint(-60, 60) * LN2, rng.randint(-200, 200)))
    return [
        ("special+hard", list(SPECIALS) + HARD_EXP + [-h for h in HARD_EXP]),
        ("binade-uniform", binade_uniform(rng, N_UNIFORM, -1074, 9)),
        ("typical [-40,40]", [rng.uniform(-40.0, 40.0) for _ in range(50_000)]),
        ("[-1,1]", [rng.uniform(-1.0, 1.0) for _ in range(50_000)]),
        ("near 0 (|x|<2^-6)", [math.copysign(2.0 ** rng.uniform(-60, -6), rng.random() - 0.5) for _ in range(150_000)]
         + around(2.0 ** -6, 2000) + around(-2.0 ** -6, 2000) + around(2.0 ** -53, 500) + around(-2.0 ** -53, 500)),
        ("the -1 end", around(m1, 3000) + [rng.uniform(-40.0, -30.0) for _ in range(30_000)]),
        ("overflow end", [rng.uniform(700.0, ovf + 0.01) for _ in range(20_000)] + around(ovf, 1000) + around(708.0, 500)),
        ("k·ln2/512, k·ln2", red),
    ]


def args_log(rng, lim, log10=False):
    edges = []
    for _ in range(40_000):     # interval edges of the 512-entry table, in random binades
        j = rng.randrange(0, 512)
        e = rng.randint(-1021, 1022)
        b = 0x3FE7FC0000000000 + j * (1 << 43) + rng.randint(-4, 4)
        b = (b & 0x000FFFFFFFFFFFFF) | (((b >> 52) + e) << 52)
        if 0 < b < 0x7FF0000000000000:
            edges.append(fb(b))
    return [
        ("special+hard", list(SPECIALS) + HARD_LN + [-1.0, -INF, -2.5]),
        ("powers of 10", [10.0 ** k for k in range(-30, 31)]
         + [ulps(10.0 ** k, d) for k in range(-25, 26) for d in range(-50, 51) if d]),
        ("binade-uniform", binade_uniform(rng, N_UNIFORM, -1074, 1023, signs=(1,))),
        ("negative", binade_uniform(rng, 20_000, -1074, 1023, signs=(-1,))),
        ("near 1", [rng.uniform(1 - 2.0 ** -6, 1 + 2.0 ** -6) for _ in range(150_000)] + around(1.0, 3000)),
        ("[0.5,2]", [rng.uniform(0.5, 2.0) for _ in range(50_000)]),
        ("subnormal x", [fb(rng.randrange(1, 1 << 52)) for _ in range(30_000)]),
        ("table edges", edges),
    ]


# Published hard pow cases (Hubrecht, Jeannerod and Zimmermann, "Towards a correctly-rounded and fast power
# function in binary64 arithmetic", ARITH 2023, slides 2, 5 and 6): libraries that misrounded them.
HARD_POW = [
    (float.fromhex("0x1.0fbcd29b56829p-1"), float.fromhex("0x1.57347b643eac2p-1")),
    (float.fromhex("0x1.30b3e414e3d3bp-1"), float.fromhex("0x1.a6c0a38da8066p-1")),
    (float.fromhex("0x1.470574d68e0afp+1"), float.fromhex("0x1.02e0706205c0ep+1")),
    (float.fromhex("0x1.f80b060553772p-1"), float.fromhex("0x1.99cp+13")),      # 49 identical bits
    (float.fromhex("0x1.524ebae943097p+1"), float.fromhex("0x1.ep-2")),         # 64 identical bits
]


def pow_specials():
    """Every C99 Annex F special pair, over the special values of x and y."""
    xs = [0.0, -0.0, 1.0, -1.0, INF, -INF, NAN, 0.5, -0.5, 2.0, -2.0, fb(1), -fb(1), 1.5, -1.5,
          fb(0x7FEFFFFFFFFFFFFF), -fb(0x7FEFFFFFFFFFFFFF), 10.0, -10.0]
    ys = [0.0, -0.0, 1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 0.5, -0.5, 2.5, -2.5, INF, -INF, NAN,
          fb(0x7FEFFFFFFFFFFFFF), -fb(0x7FEFFFFFFFFFFFFF), 2.0 ** 53, 2.0 ** 53 + 2, 2.0 ** 53 - 1,
          -(2.0 ** 53 - 1), 2.0 ** 64, 2.0 ** 63 + 2 ** 11, fb(1), -fb(1), 1e-300, 1075.0, -1075.0]
    return [(x, y) for x in xs for y in ys]


def args_pow(rng, lim):
    """(x, y) pairs over the regions that matter, by section."""
    def t_over(lo, hi, n, xgen):
        out = []
        for _ in range(n):
            x = xgen()
            l = math.log(x)
            if l == 0.0:
                continue
            out.append((x, rng.uniform(lo, hi) / l))
        return out
    pos = lambda: binade_uniform(rng, 1, -1074, 1023, signs=(1,))[0]
    near1 = lambda: 1.0 + math.copysign(2.0 ** rng.uniform(-52, -1), rng.random() - 0.5)
    exact = []
    for m in (1, 3, 5, 7, 9, 11, 13, 15, 17, 21, 25, 27, 49, 81, 125, 243, 625, 2187, 3**10, 3**20, 3**33, 5**22, 7**18):
        for e in range(-1100, 1050, 37):
            x = math.ldexp(float(m), e) if -1200 < e + m.bit_length() < 1024 else 0.0
            if not (0 < x < INF) or x != x:
                continue
            for n in range(1, 41):
                exact.append((x, float(n)))
                exact.append((x, -float(n)))
    for a in (1, 3, 5, 7, 9, 11, 13, 15, 27, 31, 101, 1001, 3**6, 65535, 2**26 - 1):
        for e in range(-540, 520, 23):
            for k, root in ((1, 2), (2, 4), (3, 8), (1, 2)):
                if a ** root >= 2 ** 53 or not (-1074 < root * e + (a ** root).bit_length() < 1024):
                    continue
                x = math.ldexp(float(a ** root), root * e)
                if not (0 < x < INF):
                    continue
                for n in range(1, 36, 2):
                    exact.append((x, n / 2 ** k))
                    exact.append((x, -n / 2 ** k))
    for e in range(-1080, 1024):                       # powers of two: exact, the tie 2^-1075, out of range
        for y in (0.5, 1.5, 2.0, 3.0, -0.5, -1.0, 0.25, 1.25):
            exact.append((math.ldexp(1.0, max(e, -1074)), y))
    exact += [(3 * 2.0 ** -215, 5.0), (9 * 2.0 ** -430, 2.5), (0.5, 1075.0), (0.5, 1074.0), (2.0, -1075.0),
              (4.0, -537.5), (0.25, 537.5), (10.0, 23.0), (10.0, 22.0), (27.0, float.fromhex("0x1.5555555555555p-2"))]
    sq = []                                            # perfect squares and cubes with y = 0.5, 1.5, 2.5, 1/3
    for _ in range(20_000):
        a = rng.randrange(1, 1 << 26)
        e = rng.randint(-330, 330)
        xs = math.ldexp(float(a * a), 2 * e)
        xc = math.ldexp(float(a ** 3), 3 * e) if a < (1 << 17) else xs
        for y in (0.5, 1.5, 2.5, -0.5, float.fromhex("0x1.5555555555555p-2")):
            sq.append((xs, y))
            sq.append((xc, y))
    negint = []
    for _ in range(60_000):
        x = -binade_uniform(rng, 1, -60, 60, signs=(1,))[0]
        n = rng.choice([rng.randint(-300, 300), rng.randint(-40, 40), 2 ** 53 + 2 * rng.randint(0, 99),
                        2 ** 53 - 1 - 2 * rng.randint(0, 99), rng.randint(-2000, 2000)])
        negint.append((x, float(n)))
    return [
        ("special", pow_specials() + HARD_POW),
        ("x binade-uniform, y in [-1100,1100]", [(pos(), rng.uniform(-1100.0, 1100.0)) for _ in range(250_000)]),
        ("x binade-uniform, |y ln x| <= 745", t_over(-745.0, 709.0, 200_000, pos)),
        ("typical x in (0,100], y in [-10,10]", [(rng.uniform(1e-3, 100.0), rng.uniform(-10.0, 10.0)) for _ in range(50_000)]),
        ("subnormal results", t_over(-745.2, -708.3, 80_000, pos)),
        ("overflow end", t_over(700.0, 709.8, 60_000, pos)),
        ("x near 1, huge y", t_over(-745.0, 709.0, 120_000, near1)),
        ("x subnormal", [(fb(rng.randrange(1, 1 << 52)), rng.uniform(-1.5, 1.5)) for _ in range(30_000)]),
        ("negative x, integer y", negint),
        ("exact and midpoint families", exact),
        ("squares/cubes, y = 1/2, 3/2, 5/2, 1/3", sq),
        ("tiny |y|", [(pos(), math.copysign(2.0 ** rng.uniform(-80, -20), rng.random() - 0.5)) for _ in range(30_000)]),
        ("y = 2.5 (x^2.5)", [(binade_uniform(rng, 1, -500, 409, signs=(1,))[0], 2.5) for _ in range(40_000)]),
    ]


HARD_ATAN = [float.fromhex(h) for h in (
    "0x1.e000000000546p-21", "0x1.22e8d75e2bc7fp-11", "0x1.a933fe176b375p-3", "0x1.7ba49f739829fp-1",
    "0x1.06b8e2c3cc000p+52")]  # the last: 4621447055448553/2^11, arctan's worst case above 2.25e12 (1 0^45)


def atan2_specials():
    vals = [0.0, -0.0, fb(1), -fb(1), 1.0, -1.0, 2.0, -0.5, fb(0x7FEFFFFFFFFFFFFF), -fb(0x7FEFFFFFFFFFFFFF),
            INF, -INF, NAN, 3.0, -3.0, 2.0 ** -1022, -(2.0 ** -1022), 1e300, -1e-300]
    return [(y, x) for y in vals for x in vals]


def _vec(rng, scale):
    return [rng.uniform(-1.0, 1.0) * scale for _ in range(3)]


def _cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def bond_and_dihedral(rng, n):
    """atan2(|ba x bc|, ba.bc) (potential::bond_angle) and atan2(|b2| b1.n2, n1.n2) (dihedral_angle), as the
    oracle computes them in doubles, at Angstrom (1) and metre (1e-10) scales, with near-collinear and
    perpendicular geometries."""
    bond, dih = [], []
    for _ in range(n):
        sc = rng.choice([1.0, 1e-10, 1.5, 1e-9, 3.7e-11])
        ba = _vec(rng, sc)
        bc = _vec(rng, sc)
        shape = rng.random()
        if shape < 0.25:                         # near-collinear (angle near 0 or pi)
            eps = 2.0 ** rng.uniform(-50, -5)
            sgn = rng.choice([1.0, -1.0])
            bc = [sgn * v * rng.uniform(0.5, 2.0) + eps * sc * rng.uniform(-1, 1) for v in ba]
        elif shape < 0.4:                        # near-perpendicular
            c = _cross(ba, _vec(rng, sc))
            nn = math.sqrt(_dot(c, c)) or 1.0
            bc = [v / nn * sc + 2.0 ** rng.uniform(-50, -5) * sc * rng.uniform(-1, 1) for v in c]
        cr = _cross(ba, bc)
        bond.append((math.sqrt(_dot(cr, cr)), _dot(ba, bc)))
        b1, b2, b3 = _vec(rng, sc), _vec(rng, sc), _vec(rng, sc)
        if shape < 0.2:
            b3 = [v + 2.0 ** rng.uniform(-40, -5) * sc * rng.uniform(-1, 1) for v in b1]
        n1 = _cross(b1, b2)
        n2 = _cross(b2, b3)
        dih.append((math.sqrt(_dot(b2, b2)) * _dot(b1, n2), _dot(n1, n2)))
    return bond, dih


def args_atan2(rng, lim):
    sgn = lambda: rng.choice((1.0, -1.0))
    quad = [(sgn() * binade_uniform(rng, 1, -1074, 1023, signs=(1,))[0],
             sgn() * binade_uniform(rng, 1, -1074, 1023, signs=(1,))[0]) for _ in range(300_000)]
    diffs = []
    for _ in range(200_000):                     # every exponent difference, -2097 ... 2097
        d = rng.randint(-2097, 2097)
        lo, hi = max(-1074, -1074 - d), min(1023, 1023 - d)
        e = rng.randint(lo, hi)
        y = binade_uniform(rng, 1, e + d, e + d, signs=(1,))[0]
        x = binade_uniform(rng, 1, e, e, signs=(1,))[0]
        diffs.append((sgn() * y, sgn() * x))
    near_pi = [(sgn() * x * 2.0 ** -rng.uniform(0, 80), -x) for x in
               (binade_uniform(rng, 1, -900, 900, signs=(1,))[0] for _ in range(60_000))]
    near_pio2 = [(sgn() * y, sgn() * y * 2.0 ** -rng.uniform(0, 80)) for y in
                 (binade_uniform(rng, 1, -900, 900, signs=(1,))[0] for _ in range(60_000))]
    sub = []
    for _ in range(60_000):                      # quotients from 2^-1100 to 2^-1015 (x > 0): subnormal, and below
        eq = rng.randint(-1100, -1015)
        ex = rng.randint(-1074 - eq, 1023)
        x = binade_uniform(rng, 1, ex, ex, signs=(1,))[0]
        y = binade_uniform(rng, 1, eq + ex, eq + ex, signs=(1,))[0]
        sub.append((sgn() * y, x))
    for _ in range(20_000):                      # exact quotients, midpoints (2K+1)*2^-1075 among them
        k = rng.randint(1, 1 << rng.randint(1, 52)) | 1
        b = rng.randint(1, 1021)
        a = rng.choice([-1075, -1075, -1076, -1074, -1060, rng.randint(-1100, -1000)])
        y = math.ldexp(float(k), a + b)
        if y != 0.0 and math.isfinite(y):
            sub.append((sgn() * y, math.ldexp(1.0, b)))
            sub.append((sgn() * y * 3.0, math.ldexp(3.0, b)))
    sub += [(3 * 2.0 ** -60, 2.0 ** 1015), (-3 * 2.0 ** -60, 2.0 ** 1015), (5 * 2.0 ** -61, 2.0 ** 1015),
            (2.0 ** -1074, 1.5), (2.0 ** -1074, 2.0 ** 1023), (fb(3), 2.0)]
    jpts = []
    for _ in range(30_000):                      # |y/x| near (j + 1/2)/256, near j/256, near 1
        j = rng.randint(0, 255)
        x = sgn() * binade_uniform(rng, 1, -500, 500, signs=(1,))[0]
        for z in ((j + 0.5) / 256, j / 256 + 2.0 ** -rng.uniform(20, 52), 1.0 - 2.0 ** -rng.uniform(1, 53)):
            jpts.append((sgn() * abs(x) * z, x))
    jpts += [(v, v) for v in binade_uniform(rng, 2000, -1074, 1023)] + \
        [(v, -v) for v in binade_uniform(rng, 2000, -1074, 1023)]
    bond, dih = bond_and_dihedral(rng, 100_000)
    hard = [(h, 1.0) for h in HARD_ATAN] + [(-h, 1.0) for h in HARD_ATAN] + \
        [(h * 2.0 ** k, 2.0 ** k) for h in HARD_ATAN for k in (-1000, -60, -1, 1, 60, 900)] + \
        [(1.0, h) for h in HARD_ATAN] + [(-1.0, -h) for h in HARD_ATAN]
    return [
        ("special", atan2_specials() + hard),
        ("quadrants, binade-uniform", quad),
        ("every exponent difference", diffs),
        ("near +-pi (|y| << |x|, x < 0)", near_pi),
        ("near +-pi/2 (|x| << |y|)", near_pio2),
        ("subnormal and exact quotients", sub),
        ("|y/x| near j/256 and 1, |y| = |x|", jpts),
        ("bond angles", bond),
        ("dihedral angles", dih),
    ]


def args_tanh(rng, lim):
    """ushma's own (M0 part 2c-3): the regions where _us_cr_tanh changes path or the rounding is delicate."""
    one = lim["TANH_ONE"]
    tiny = 2.0 ** -27
    red = []
    for _ in range(30_000):     # 2a at the exp core's reduction points: K·ln2/512 and (K + 1/2)·ln2/512, 2a ≤ 38.2
        kk = rng.randint(1, 28_200)
        red.append(ulps((kk + 0.5) * LN2 / 1024, rng.randint(-3, 3)))
        red.append(ulps(kk * LN2 / 1024, rng.randint(-3, 3)))
    sgn = lambda: rng.choice((1.0, -1.0))
    return [
        ("special", list(SPECIALS) + [one, -one, ulps(one, 1), -ulps(one, 1), tiny, -tiny, ulps(tiny, -1),
                                      -ulps(tiny, -1), 2.0 ** -7, -2.0 ** -7, 19.0, 20.0, 50.0, 709.0]),
        ("binade-uniform |x| < 64", binade_uniform(rng, N_UNIFORM, -1074, 5)),
        ("typical [-20,20]", [rng.uniform(-20.0, 20.0) for _ in range(100_000)]),
        ("fins: mL in [0,50]", [rng.uniform(0.0, 50.0) for _ in range(50_000)]),
        ("[-1,1]", [rng.uniform(-1.0, 1.0) for _ in range(50_000)]),
        ("tiny end 2^-27", around(tiny, 3000) + around(-tiny, 3000)
         + [sgn() * 2.0 ** rng.uniform(-40, -20) for _ in range(60_000)]),
        ("expm1 branch 2a = 2^-6", around(2.0 ** -7, 2000) + around(-2.0 ** -7, 2000)),
        ("the 1 end", around(one, 3000) + around(-one, 3000) + [sgn() * rng.uniform(15.0, 19.1) for _ in range(30_000)]),
        ("2a at K·ln2/512", red),
    ]


def gen_args(name, lim):
    rng = random.Random({"exp": 1001, "expm1": 1003, "ln": 1004, "pow": 1006, "atan2": 1009, "tanh": 1040}[name])
    sections = {"exp": args_exp, "expm1": args_expm1, "ln": args_log, "pow": args_pow, "atan2": args_atan2,
                "tanh": args_tanh}[name](rng, lim)
    xs = [x for _, sec in sections for x in sec]
    spans = []
    i = 0
    for tag, sec in sections:
        spans.append((tag, i, i + len(sec)))
        i += len(sec)
    head = len(sections[0][1])
    acc = xs[:head] + [xs[i] for i in random.Random(77).sample(range(len(xs)), min(N_ACC, len(xs)))]
    return xs, acc, spans


PROBE = r'''# GENERATED by scripts/oracle/cr_verify.py. Reads job.bin (u64 words: function id, mode, n, then n arguments,
# two words each for pow and atan2) and writes res.bin. Mode 0: for each argument the function's result, the fast
# path's result and one byte, 1 if the fast path's test passed. Mode 1: the accurate path's result alone.
# The probe's own names carry `crv_`: the stdlib leaves share one namespace with it (lib/process.cyr defines `run`).
include "src/mathx_data.cyr"
include "src/mathx.cyr"

var WORDS = WMAX;

# The arity of function id fid, in argument words.
#must_use
fn crv_arity(fid): i64 {
    if (fid == 5 || fid == 8) { return 2; }
    return 1;
}

# (full result, fast result, fast flag) in mode 0, (accurate result, accurate result, 0) in mode 1.
#must_use
fn crv_run(fid, mode, x: f64, y: f64): (f64, f64, i64) {
    if (mode == 1) {
        if (fid == 0) { var a: f64 = _us_cr_exp_acc(x); return (a, a, 0); }
        if (fid == 2) { var c: f64 = _us_cr_expm1_acc(x); return (c, c, 0); }
        if (fid == 3) { var d: f64 = _us_cr_ln_acc(x); return (d, d, 0); }
        if (fid == 8) { var ta: f64 = _us_cr_atan2_acc(x, y); return (ta, ta, 0); }
        if (fid == 40) { var th: f64 = _us_cr_tanh_acc(x); return (th, th, 0); }
        var f: f64 = _us_cr_pow_acc(x, y);
        return (f, f, 0);
    }
    if (fid == 0) { var g: f64 = _us_cr_exp(x); var gf, go = _us_cr_exp_fast(x); return (g, gf, go); }
    if (fid == 2) { var i: f64 = _us_cr_expm1(x); var if1, io = _us_cr_expm1_fast(x); return (i, if1, io); }
    if (fid == 3) { var j: f64 = _us_cr_ln(x); var jf, jo = _us_cr_ln_fast(x); return (j, jf, jo); }
    if (fid == 8) { var o: f64 = _us_cr_atan2(x, y); var of, oo = _us_cr_atan2_fast(x, y); return (o, of, oo); }
    if (fid == 40) { var t: f64 = _us_cr_tanh(x); var tf, to = _us_cr_tanh_fast(x); return (t, tf, to); }
    var l: f64 = _us_cr_pow(x, y);
    var lf, lo = _us_cr_pow_fast(x, y);
    return (l, lf, lo);
}

# Reads the job, runs every argument, writes the results.
fn main(): i64 {
    alloc_init();
    var cap = 8 * (3 + WORDS);
    var buf = alloc(cap);
    var got = file_read_all("job.bin", buf, cap);
    if (got < 24) { return 2; }
    var fid = load64(buf);
    var mode = load64(buf + 8);
    var n = load64(buf + 16);
    var ar = crv_arity(fid);
    if (got != 8 * (3 + ar * n)) { return 3; }
    var out = alloc(17 * n + 8);
    var i = 0;
    while (i < n) {
        var at = buf + 24 + 8 * ar * i;
        var x: f64 = load64(at);
        var y: f64 = 0.0;
        if (ar == 2) { y = load64(at + 8); }
        var r, rf, ok = crv_run(fid, mode, x, y);
        store64(out + 8 * i, r);
        store64(out + 8 * n + 8 * i, rf);
        store8(out + 16 * n + i, ok);
        i = i + 1;
    }
    var w = file_write_all("res.bin", out, 17 * n);
    if (w != 17 * n) { return 4; }
    return 0;
}

var r = main();
sys_exit_group(r);
'''

TARGETS = {
    "elf": ([], "probe", []),
    "pe": (["--win"], "probe.exe", ["wine"]),
    "aarch64": (["--aarch64"], "probe-aarch64", ["qemu-aarch64"]),
}
# Diagnostics from vendored lib/ (docs/development/state.md); anything else fails the build.
KNOWN = ("large static data", "refusing to overwrite stdlib leaf", "warning:lib/")


def build(root, work, wmax, targets):
    src = work / "probe.cyr"
    src.write_text(PROBE.replace("WMAX", str(wmax)))
    for t in targets:
        flags, out, _ = TARGETS[t]
        r = subprocess.run(["cyrius", "build", *flags, str(src), str(work / out)], cwd=root,
                           capture_output=True, text=True)
        log = r.stdout + r.stderr
        bad = [ln for ln in log.splitlines() if "warning" in ln and not any(k in ln for k in KNOWN)]
        if r.returncode != 0 or bad:
            print(log)
            print(f"cr_verify: probe build failed for {t}")
            sys.exit(2)


def run_probe(work, target, fid, mode, args):
    """(results, fast results, fast flags) for the arguments (floats, or (x, y) pairs for pow and atan2)."""
    _, out, runner = TARGETS[target]
    n = len(args)
    flat = [bits(v) for a in args for v in (a if isinstance(a, tuple) else (a,))]
    (work / "job.bin").write_bytes(struct.pack(f"<3Q{len(flat)}Q", fid, mode, n, *flat))
    res = work / "res.bin"
    for _attempt in range(2):            # wine has been seen to exit 0 without writing the file, once
        if res.exists():
            res.unlink()
        r = subprocess.run([*runner, str(work / out)], cwd=work, capture_output=True)
        if r.returncode == 0 and res.exists():
            break
    if r.returncode != 0 or not res.exists():
        print(r.stdout.decode(errors="replace")[-2000:], r.stderr.decode(errors="replace")[-2000:])
        print(f"cr_verify: probe failed on {target} (exit {r.returncode})")
        sys.exit(2)
    data = res.read_bytes()
    words = struct.unpack(f"<{2 * n}Q", data[:16 * n])
    return list(words[:n]), list(words[n:]), data[16 * n:17 * n]


def same(a, b):
    """Bit-identical, or both NaN (the results compared against MPFR, whose NaN has no sign)."""
    if a == b:
        return True
    return math.isnan(fb(a)) and math.isnan(fb(b))


def ref_of(ref, a):
    return ref(*a) if isinstance(a, tuple) else ref(a)


def _libm(name, arity):
    import ctypes
    f = getattr(ctypes.CDLL("libm.so.6"), name)
    f.argtypes = [ctypes.c_double] * arity
    f.restype = ctypes.c_double
    return f


def special_checks(work, targets, plan):
    """The special values of pow, atan2 and tanh: ours, MPFR and glibc. True if ours differs from MPFR anywhere."""
    failed = False
    for name, fid, sp, f, ar in (("pow", 5, pow_specials() + HARD_POW if "pow" in plan else [], "pow", 2),
                                 ("atan2", 8, atan2_specials() if "atan2" in plan else [], "atan2", 2),
                                 ("tanh", 40, SPECIALS if "tanh" in plan else [], "tanh", 1)):
        if not sp:
            continue
        g = _libm(f, ar)
        ref = {"pow": mpfr_ref.pow, "atan2": mpfr_ref.atan2, "tanh": mpfr_ref.tanh}[name]
        got, _, _ = run_probe(work, targets[0], fid, 0, sp)
        n_om = n_og = n_mg = 0
        for i, a in enumerate(sp):
            m = bits(ref_of(ref, a))
            gl = bits(g(*a) if isinstance(a, tuple) else g(a))
            if not same(m, gl):
                n_mg += 1
                print(f"  special: MPFR and glibc differ at {name}{a!r}: MPFR {m:016x}, glibc {gl:016x}")
            n_om += not same(got[i], m)
            n_og += not same(got[i], gl)
        print(f"{name} special values: {len(sp)}; ours vs MPFR {n_om} differ, ours vs glibc {n_og}, "
              f"MPFR vs glibc {n_mg}")
        failed |= bool(n_om)
    return failed


def main():
    argv = sys.argv[1:]
    root = Path(__file__).resolve().parents[2]
    work = root / "build/cr-verify"
    targets = ["elf", "pe", "aarch64"]
    if "--work" in argv:
        work = Path(argv[argv.index("--work") + 1]).resolve()
    if "--targets" in argv:
        targets = argv[argv.index("--targets") + 1].split(",")
    only = argv[argv.index("--only") + 1].split(",") if "--only" in argv else None
    work.mkdir(parents=True, exist_ok=True)
    tables = _load("mathx_tables")
    lim = {k: fb(v) for k, v in tables.limits().items()}
    lim.update({k: fb(v) for k, v in tables.tanh_limits().items()})

    funcs = [f for f in FUNCS if only is None or f[1] in only]
    plan = {}
    wmax = 0
    for fid, name, _, _, _ in funcs:
        xs, acc, spans = gen_args(name, lim)
        plan[name] = (xs, acc, spans)
        wmax = max(wmax, ARITY.get(fid, 1) * len(xs))
    wmax = max(wmax, 2 * 400_000)
    build(root, work, wmax, targets)
    print(f"cr_verify: MPFR {mpfr_ref.VERSION}; targets {', '.join(targets)}; probe in {work}")

    failed = False
    for fid, name, full, accname, ref in funcs:
        xs, acc, spans = plan[name]
        t0 = time.time()
        want = [bits(ref_of(ref, a)) for a in xs]
        want_acc = [bits(ref_of(ref, a)) for a in acc]
        first = None
        for t in targets:
            got, gotf, flags = run_probe(work, t, fid, 0, xs)
            got_acc, _, _ = run_probe(work, t, fid, 1, acc)
            bad = [i for i in range(len(xs)) if not same(got[i], want[i])]
            bad_acc = [i for i in range(len(acc)) if not same(got_acc[i], want_acc[i])]
            bad_fast = [i for i in range(len(xs)) if flags[i] == 1 and gotf[i] != got[i]]
            for i in bad[:5]:
                print(f"  MISMATCH {full}{xs[i]!r} on {t}: got {got[i]:016x}, want {want[i]:016x}")
            for i in bad_acc[:5]:
                print(f"  MISMATCH {accname}{acc[i]!r} on {t}: got {got_acc[i]:016x}, want {want_acc[i]:016x}")
            for i in bad_fast[:5]:
                print(f"  FAST/FULL DIFFER {full}_fast{xs[i]!r} on {t}: {gotf[i]:016x} vs {got[i]:016x}")
            slow = len(xs) - sum(flags)
            if first is None:
                first = (got, got_acc, flags)
                xdiff = 0
            else:
                xdiff = sum(1 for a, b in zip(got, first[0]) if a != b) + \
                    sum(1 for a, b in zip(got_acc, first[1]) if a != b) + \
                    sum(1 for a, b in zip(flags, first[2]) if a != b)
            print(f"{name:6} {t:8} {len(xs):9,} args: {len(bad)} mismatches; accurate path forced on "
                  f"{len(acc):,}: {len(bad_acc)} mismatches; fast path: {len(bad_fast)} differ, fell back "
                  f"{slow:,} ({slow / len(xs):.2e}); bits or flags differing from {targets[0]}: {xdiff}")
            failed |= bool(bad or bad_acc or bad_fast or xdiff)
            if t == targets[0]:
                rates = [f"{tag} {(b - a) - sum(flags[a:b])}/{b - a}" for tag, a, b in spans]
                print("         fallbacks by region: " + "; ".join(rates))
        print(f"       [{time.time() - t0:.0f} s]")

    if "pow" in plan:
        # x^2.5: ushma's Redlich–Kwong Tc^2.5 (state::redlich_kwong_pressure), against MPFR over every binade.
        rng = random.Random(2525)
        x25 = binade_uniform(rng, 200_000, -1074, 1023, signs=(1,)) + \
            [rng.uniform(5.0, 1000.0) for _ in range(100_000)] + \
            [9 * 2.0 ** -430, 2.0 ** -387, 2.0 ** -386, ulps(2.0 ** -387, -1), 4.0, 2.0 ** 409.6, 2.0 ** 410]
        c25, _, _ = run_probe(work, targets[0], 5, 0, [(x, 2.5) for x in x25])
        w25 = [bits(mpfr_ref.pow(x, 2.5)) for x in x25]
        cr_bad = sum(1 for i in range(len(x25)) if not same(c25[i], w25[i]))
        print(f"x^2.5 on {len(x25):,} arguments: _us_cr_pow(x, 2.5) vs MPFR {cr_bad} mismatches")
        failed |= bool(cr_bad)
    failed |= special_checks(work, targets, plan)
    print("cr_verify: " + ("FAILED" if failed else "all results correctly rounded, identical on every target"))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
