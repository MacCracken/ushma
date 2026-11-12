"""Correctly rounded binary64 results from the system MPFR library (ADR 0004).

Each function returns the IEEE-754 binary64 round-to-nearest-even value of the exact mathematical
result: MPFR at 53 bits of precision with binary64's exponent range (emin = -1073, emax = 1024) and
`mpfr_subnormalize`, so subnormal results, overflow to infinity and exact ties (pow(3*2^-215, 5))
round exactly as IEEE round-to-nearest does. A decimal reference at fixed precision can't do that:
it misjudges exact ties (it rounds pow(0.5, 1075) = 2^-1075 up instead of to 0). MPFR's special
values follow ISO C99 Annex F, as glibc's do (pow(-1, ±inf) = 1, pow(x, ±0) = 1, atan(±inf) = ±π/2).

The core is kimiya's scripts/oracle/mpfr_ref.py (`cr`, `cr2` and the self-checks that apply to
ushma's primitives: exp, expm1, log, pow, atan, tanh); kimiya's unbounded-exponent machinery is not
needed here. No Python package is needed: MPFR (libmpfr.so.6, MPFR 4.x) is called through ctypes.
Run `python3 -I -B scripts/oracle/mpfr_ref.py` for the self-checks.
"""
import ctypes
import ctypes.util
import math
import struct

_lib = ctypes.util.find_library("mpfr") or "libmpfr.so.6"
_m = ctypes.CDLL(_lib)

# mpfr_t on LP64: _mpfr_prec (long), _mpfr_sign (int, padded), _mpfr_exp (long), _mpfr_d (limb *).
_MPFR_T = ctypes.c_byte * 32
_RNDN = 0
_P = ctypes.c_void_p

_m.mpfr_get_version.restype = ctypes.c_char_p
_m.mpfr_set_emin.argtypes = [ctypes.c_long]
_m.mpfr_set_emax.argtypes = [ctypes.c_long]
_m.mpfr_init2.argtypes = [_P, ctypes.c_long]
_m.mpfr_set_d.argtypes = [_P, ctypes.c_double, ctypes.c_int]
_m.mpfr_get_d.argtypes = [_P, ctypes.c_int]
_m.mpfr_get_d.restype = ctypes.c_double
_m.mpfr_subnormalize.argtypes = [_P, ctypes.c_int, ctypes.c_int]

_UNARY = ("exp", "expm1", "log", "atan", "tanh", "sqrt")
_BINARY = ("pow", "atan2")
for _name in _UNARY:
    getattr(_m, "mpfr_" + _name).argtypes = [_P, _P, ctypes.c_int]
for _name in _BINARY:
    getattr(_m, "mpfr_" + _name).argtypes = [_P, _P, _P, ctypes.c_int]

_m.mpfr_set_emin(-1073)
_m.mpfr_set_emax(1024)
_a = _MPFR_T()
_b = _MPFR_T()
_r = _MPFR_T()
for _v in (_a, _b, _r):
    _m.mpfr_init2(ctypes.byref(_v), 53)

VERSION = _m.mpfr_get_version().decode()


def bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(b):
    return struct.unpack("<d", struct.pack("<Q", b))[0]


def _finish(t):
    t = _m.mpfr_subnormalize(ctypes.byref(_r), t, _RNDN)
    return _m.mpfr_get_d(ctypes.byref(_r), _RNDN)


def cr(name, x):
    """The correctly rounded binary64 value of the unary MPFR function `name` at x (a float)."""
    if name not in _UNARY:
        raise ValueError(name)
    _m.mpfr_set_d(ctypes.byref(_a), x, _RNDN)
    t = getattr(_m, "mpfr_" + name)(ctypes.byref(_r), ctypes.byref(_a), _RNDN)
    return _finish(t)


def cr2(name, x, y):
    """The correctly rounded binary64 value of the binary MPFR function `name` (pow, atan2)."""
    if name not in _BINARY:
        raise ValueError(name)
    _m.mpfr_set_d(ctypes.byref(_a), x, _RNDN)
    _m.mpfr_set_d(ctypes.byref(_b), y, _RNDN)
    t = getattr(_m, "mpfr_" + name)(ctypes.byref(_r), ctypes.byref(_a), ctypes.byref(_b), _RNDN)
    return _finish(t)


def exp(x): return cr("exp", x)
def expm1(x): return cr("expm1", x)
def log(x): return cr("log", x)
def atan(x): return cr("atan", x)
def tanh(x): return cr("tanh", x)
def pow(x, y): return cr2("pow", x, y)
def atan2(y, x): return cr2("atan2", y, x)


if __name__ == "__main__":
    # Exact ties, rounded to even (glibc 2.44's pow returns 121 here, kimiya port-parity.md §8).
    assert bits(pow(3 * 2.0 ** -215, 5)) == 122
    assert bits(pow(9 * 2.0 ** -430, 2.5)) == 122
    assert pow(0.5, 1075) == 0.0 and bits(pow(0.5, 1074)) == 1
    # Subnormal and overflow ends.
    assert bits(exp(-745.1332191019411)) == 1 and exp(-745.2) == 0.0
    assert exp(709.78) < math.inf and exp(709.8) == math.inf
    # Signed zeros, specials.
    assert bits(expm1(-0.0)) == bits(-0.0) and log(1.0) == 0.0 and math.isnan(log(-1.0))
    assert log(0.0) == -math.inf and exp(-math.inf) == 0.0 and expm1(-math.inf) == -1.0
    assert bits(tanh(-0.0)) == bits(-0.0) and tanh(math.inf) == 1.0 and tanh(-30.0) == -1.0
    assert atan(math.inf) == math.pi / 2 and bits(atan(-0.0)) == bits(-0.0)
    # C99 Annex F special cases of pow, as glibc's (ushma's pow(-1, DBL_MAX) row: ganita gives +inf).
    assert pow(-1.0, math.inf) == 1.0 and pow(-1.0, from_bits(0x7FEFFFFFFFFFFFFF)) == 1.0
    assert pow(math.nan, 0.0) == 1.0 and pow(1.0, math.nan) == 1.0 and pow(2.0, 0.5) == math.sqrt(2.0)
    # Agrees with glibc where glibc is correctly rounded.
    assert exp(1.0) == math.e and log(10.0) == math.log(10.0) and atan2(1.0, -1.0) == 3 * math.pi / 4
    print(f"mpfr_ref: MPFR {VERSION}: self-checks pass")
