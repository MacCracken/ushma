"""Vendor impetus's shortest-digit f64 printer as ushma's src/fmtf.cyr (ADR 0005, port plan P7).

The source is impetus's src/fmtf.cyr at a pinned content hash (a sibling checkout at ../impetus). The copy is
mechanical: impetus's header is replaced by ushma's, `_im_` / `_IM_` become `_us_` / `_US_`, the two public
printers `im_fmt_f64_display` / `im_fmt_f64_json` become ushma-private `_us_fmt_f64_display` / `_us_fmt_f64_json`,
and the one helper it took from impetus's mathx.cyr (`_im_bits`, a bit cast) is defined in the file. Nothing else
changes, so tests/golden/fmtf.tsv (both forms over 13,044 doubles) pins it the way impetus's own golden file pins the
original.

    python3 -I -B scripts/vendor_fmtf.py          # write src/fmtf.cyr
    python3 -I -B scripts/vendor_fmtf.py --check  # src/fmtf.cyr is what this script writes (needs ../impetus)
"""
import hashlib
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOURCE = os.path.join(os.path.dirname(ROOT), "impetus", "src", "fmtf.cyr")
SOURCE_SHA256 = "1d925be3b48c274e52f22b20360a929ed1237278e40307f7649031356ba19d57"
SOURCE_COMMIT = "809ea0f"
TARGET = os.path.join(ROOT, "src", "fmtf.cyr")

HEADER = """\
# fmtf.cyr — f64 to text: the shortest round-trip digits, in Rust's Display and serde_json's forms (ADR 0005).
#
# Vendored mechanically from impetus's src/fmtf.cyr (impetus commit {commit}, sha256 {sha}…) by
# scripts/vendor_fmtf.py; regenerate it there, never edit it here. The only changes: this header, the `_us_` names,
# and `_us_bits` (impetus takes it from its mathx.cyr).
#
# Both printers give the SHORTEST digits that parse back to the same double, choosing the closest to the exact value
# among them: Rust's `Display` (core's Grisu with a Dragon fallback) and serde_json 1.0.154's zmij (Schubfach). They
# differ only on an exact tie between two closest candidates: core rounds it up, zmij to the even digit. The digits
# come from Burger & Dybvig's free-format algorithm ("Printing Floating-Point Numbers Quickly and Accurately", 1996)
# on a small bignum: exact, with no tables. tests/golden/fmtf.tsv pins both forms against Rust 1.99.1 and serde_json
# 1.0.154 over every binade, the powers of ten, the form thresholds and a seeded sweep (M0 probe 11 checked 160,378
# doubles: no difference). bayan's own printer differs from serde_json on 6,276 of them (filed upstream).
#
# Forms:
# - `_us_fmt_f64_display`: Rust's `{{}}` — always positional, no ".0" (`100`, `0.1`, `-0`), `NaN`, `inf`, `-inf`.
# - `_us_fmt_f64_json`: serde_json — positional while the decimal exponent is in −5..=15, with ".0" on an integer
#   (`100.0`, `0.00001`); otherwise `d.ddde±x` (`1e+16`, `1.5e-7`); `null` for NaN and ±∞.

# The bits of an f64 as an integer word: an untyped parameter, an i64 return.
#must_use
fn _us_bits(x): i64 {{
    return x;
}}
"""


def vendor(src):
    """ushma's fmtf.cyr from impetus's text."""
    lines = src.split("\n")
    # impetus's header is the leading comment block; the code starts at its first non-comment line.
    i = 0
    while i < len(lines) and (lines[i].startswith("#") or lines[i].strip() == ""):
        i += 1
    body = "\n".join(lines[i:])
    body = re.sub(r"\bim_fmt_f64_(display|json)\b", r"_us_fmt_f64_\1", body)
    body = re.sub(r"\b_im_", "_us_", body)
    body = re.sub(r"\b_IM_", "_US_", body)
    body = body.replace("the snapshot's JSON form, ADR 0013", "a JSON writer")
    if re.search(r"\b_?[iI][mM]_", body) or "impetus" in body:
        raise SystemExit("vendor_fmtf: an impetus name survived the rename")
    head = HEADER.format(commit=SOURCE_COMMIT, sha=SOURCE_SHA256[:16])
    return head + "\n" + body


def main():
    check = "--check" in sys.argv[1:]
    if not os.path.exists(SOURCE):
        if check:
            print(f"vendor_fmtf: SKIP (no impetus checkout at {SOURCE})")
            return 0
        raise SystemExit(f"vendor_fmtf: no impetus checkout at {SOURCE}")
    src = open(SOURCE, encoding="utf-8").read()
    sha = hashlib.sha256(src.encode("utf-8")).hexdigest()
    if sha != SOURCE_SHA256:
        msg = f"vendor_fmtf: impetus's fmtf.cyr is {sha[:16]}…, not the pinned {SOURCE_SHA256[:16]}…"
        if check:
            print(msg + " (re-vendor deliberately: update the pin and re-run the fmtf golden rows)")
            return 1
        raise SystemExit(msg)
    out = vendor(src)
    if check:
        cur = open(TARGET, encoding="utf-8").read() if os.path.exists(TARGET) else ""
        if cur != out:
            print("vendor_fmtf: src/fmtf.cyr is stale; run scripts/vendor_fmtf.py")
            return 1
        print("vendor_fmtf: src/fmtf.cyr current (impetus " + SOURCE_COMMIT + ")")
        return 0
    open(TARGET, "w", encoding="utf-8").write(out)
    print(f"vendor_fmtf: wrote src/fmtf.cyr ({out.count(chr(10))} lines)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
