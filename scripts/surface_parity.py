#!/usr/bin/env python3
"""Surface parity (port plan §6, every milestone's gate): every public item of the Rust oracle has its Cyrius name.

    python3 -I -B scripts/surface_parity.py .

For each Rust module whose milestone has come (`PORTED`), it lists rust-old/src/<module>.rs's public items before
`#[cfg(test)]` and checks that src/*.cyr defines the name ADR 0001's rules (P6, the written `us_` prefix) give:

- `pub fn f`                                -> `fn us_f(` (`fn us_log_f(` in logging)
- `impl T { pub fn m }`                     -> `fn us_<t>_m(`, with <t> from TYPES
- `pub struct T` / `pub enum T`             -> `struct <TYPES[T]>` / `enum <TYPES[T]>`
- `#[derive(Serialize)]` / `Deserialize`    -> `fn us_<t>_to_json(` / `fn us_<t>_from_json(` (ADR 0005)
- `pub const X: f64` (any module)           -> `const US_X =` (src/constants.cyr, generated)
- `pub const X: T`, T a preset type         -> `const US_<T>_X =` (`US_X` when X already starts `<T>_`) and
                                               `fn us_<t>_preset(`
- `UshmaError::Variant`                     -> `US_ERR_<SNAKE>` in `enum UsError`

Where the rules give no single name, MAPPED lists the Cyrius names that carry the item, and they are checked too:
there are no exemptions. Every file under rust-old/src is PORTED, SCHEDULED (its milestone not yet come; counted) or
WIRING (module declarations, which Cyrius's flat namespace and the [lib] / [lib.ai] lists replace). Exits 1 on a
missing name, an unaccounted file or a MAPPED entry that names no Rust item.
"""
import re
import sys
from pathlib import Path

root = Path(sys.argv[1]).resolve()

# The modules whose milestone has come.
PORTED = ["error", "validate", "logging", "material", "transfer"]
# The rest, with their milestone (port plan §3.1).
SCHEDULED = {
    "state": "M2", "entropy": "M2", "phase": "M3", "steam": "M3", "cycle": "M4", "chem": "M5", "stat": "M5",
    "numerical": "M6", "bridge": "M7", "integration/soorat": "M7", "ai": "M7",
}
WIRING = {
    "lib": "module declarations and re-exports: the [lib] / [lib.ai] module lists in cyrius.cyml",
    "integration/mod": "the `soorat` module declaration",
}
# The free-function prefix per module (default `us_`).
PREFIX = {"logging": "us_log_"}
# Rust type -> (Cyrius type, the `us_<t>_` prefix of its functions and preset constants).
TYPES = {
    "ThermalMaterial": ("UsMaterial", "material"),
    "UshmaError": ("UsError", "error"),
    "GasData": ("UsGas", "gas"),
}
# Rust item -> the Cyrius names that carry it ("fn x", "const X", "enum X", "struct X"), where the rules give none.
MAPPED = {
    "material::const ALL_MATERIALS": ["const US_MATERIAL_COUNT", "fn us_material_preset"],
    "state::const ALL_GASES": ["const US_GAS_COUNT", "fn us_gas_preset"],
    "error::type Result": ["fn is_ok", "fn is_err_result"],
    "error::UshmaError: Display": ["fn us_error_message", "fn us_error_name"],
}


def snake(camel):
    return re.sub(r"(?<!^)([A-Z])", r"_\1", camel).lower()


def rust_source(module):
    src = (root / f"rust-old/src/{module}.rs").read_text()
    cut = src.find("#[cfg(test)]")
    return src if cut < 0 else src[:cut]


def rust_items(module):
    """(label, kind, type, name) for each public item."""
    src = rust_source(module)
    items = []
    impl = None
    depth = 0
    derives = []
    for line in src.splitlines():
        s = line.strip()
        dm = re.match(r"#\[derive\(([^)]*)\)\]", s)
        if dm and depth == 0:
            derives = [d.strip().split("::")[-1] for d in dm.group(1).split(",")]
        if depth == 0:
            im = re.match(r"impl(?:<[^>]*>)?\s+(?:(\w+(?:::\w+)*)\s+for\s+)?(\w+)", s)
            if im:
                impl = im.group(2) if im.group(1) is None else None
            tm = re.match(r"pub (struct|enum) (\w+)", s)
            if tm:
                kind, ty = tm.groups()
                items.append((f"{kind} {ty}", kind, ty, ty))
                for d in ("Serialize", "Deserialize"):
                    if d in derives:
                        items.append((f"{ty}: {d}", d, ty, ty))
                if "Error" in derives:
                    items.append((f"{ty}: Display", "Display", ty, ty))
                derives = []
            cm = re.match(r"pub const ([A-Z][A-Z0-9_]*)\s*:\s*(&?\[?&?\w+)", s)
            if cm:
                name, ty = cm.groups()
                items.append((f"const {name}", "const", ty.lstrip("&["), name))
            am = re.match(r"pub type (\w+)", s)
            if am:
                items.append((f"type {am.group(1)}", "type", None, am.group(1)))
        fm = re.match(r"pub (?:const )?fn (\w+)", s)
        if fm:
            if impl and depth >= 1:
                items.append((f"{impl}::{fm.group(1)}", "method", impl, fm.group(1)))
            elif depth == 0:
                items.append((fm.group(1), "fn", None, fm.group(1)))
        depth += line.count("{") - line.count("}")
        if depth == 0 and s.endswith("}"):
            impl = None
    if module == "error":
        body = src[src.index("pub enum UshmaError"):]
        body = body[: body.index("\n}")]
        for v in re.findall(r"^\s{4}([A-Z]\w+)", body, re.M):
            items.append((f"UshmaError::{v}", "variant", "UshmaError", v))
    return items


cyr = "\n".join(p.read_text() for p in sorted((root / "src").glob("*.cyr")))
stdlib = "\n".join(p.read_text() for p in sorted((root / "lib").glob("*.cyr")))


def defined(spec, text=None):
    """`fn x`, `const X`, `struct X` or `enum X`, defined at the start of a line."""
    kind, name = spec.split(" ", 1)
    pat = {"fn": rf"^fn {name}\(", "const": rf"^(?:const|var) {name}\s*(?::\s*\w+\s*)?=",
           "struct": rf"^struct {name}\b", "enum": rf"^enum {name}\b"}[kind]
    return re.search(pat, text if text is not None else cyr, re.M) is not None


def counterpart(module, kind, ty, name):
    """The Cyrius names an item's rule gives, or None when the rule gives none."""
    t = TYPES.get(ty)
    if kind == "fn":
        return [f"fn {PREFIX.get(module, 'us_')}{name}"]
    if kind == "method" and t:
        return [f"fn us_{t[1]}_{name}"]
    if kind in ("struct", "enum") and t:
        return [f"{kind} {t[0]}"]
    if kind == "Serialize" and t:
        return [f"fn us_{t[1]}_to_json"]
    if kind == "Deserialize" and t:
        return [f"fn us_{t[1]}_from_json"]
    if kind == "const" and ty in ("f64", "usize"):
        return [f"const US_{name}"]
    if kind == "const" and t:
        # A preset already named by its type (GAS_NITROGEN) keeps its name; others take it (COPPER: US_MATERIAL_COPPER).
        pre = t[1].upper() + "_"
        return [f"const US_{name if name.startswith(pre) else pre + name}", f"fn us_{t[1]}_preset"]
    if kind == "variant":
        return [f"const US_ERR_{snake(name).upper()}"]
    return None


def variant_defined(spec):
    """An error code is a member of `enum UsError`."""
    name = spec.split(" ", 1)[1]
    m = re.search(r"^enum UsError \{(.*?)^\}", cyr, re.M | re.S)
    return m is not None and re.search(rf"\b{name}\s*=", m.group(1)) is not None


missing = 0
seen_mapped = set()
for module in PORTED:
    items = rust_items(module)
    found = 0
    for label, kind, ty, name in items:
        key = f"{module}::{label}"
        names = MAPPED.get(key)
        if names is not None:
            seen_mapped.add(key)
        else:
            names = counterpart(module, kind, ty, name)
        if names is None:
            missing += 1
            print(f"MISSING  {key}: no rule gives its Cyrius name and MAPPED has no entry")
            continue
        absent = [n for n in names if not (variant_defined(n) if kind == "variant" else
                                          defined(n) or (key in MAPPED and defined(n, stdlib)))]
        if absent:
            missing += 1
            print(f"MISSING  {key}: {', '.join(absent)}")
        else:
            found += 1
    print(f"{module:<20} {found}/{len(items)} Rust public items have their Cyrius names")
for key in MAPPED:
    module, label = key.split("::", 1)
    if key in seen_mapped:
        continue
    if module in SCHEDULED and label in {lab for lab, *_ in rust_items(module)}:
        continue  # mapped ahead of its milestone; checked once the module is PORTED
    missing += 1
    print(f"STALE    MAPPED names {key}, which is no public Rust item")
# Public constants are generated for every module at once (src/constants.cyr), so they are checked everywhere.
consts = 0
for module in SCHEDULED:
    for label, kind, ty, name in rust_items(module):
        if kind == "const" and ty in ("f64", "usize"):
            consts += 1
            if not defined(f"const US_{name}"):
                missing += 1
                print(f"MISSING  {module}::{label}: const US_{name}")
print(f"{'constants':<20} {consts} public constants of the scheduled modules have their Cyrius names")
for module, ms in SCHEDULED.items():
    print(f"{module:<20} {len(rust_items(module))} Rust public items, scheduled for {ms}")
known = set(PORTED) | set(SCHEDULED) | set(WIRING)
for path in sorted((root / "rust-old/src").rglob("*.rs")):
    mod = str(path.relative_to(root / "rust-old/src"))[:-3]
    if mod not in known:
        missing += 1
        print(f"UNKNOWN  rust-old/src/{mod}.rs is not ported, scheduled or wiring")
for mod, why in WIRING.items():
    print(f"{mod:<20} wiring: {why}")
print(f"surface parity: {missing} missing")
sys.exit(1 if missing else 0)
