# `[deps.sakshi]` is declared but never called, so every consumer that also uses log or sandhi resolves with a warning — open

**Status:** 🟡 **OPEN**. Verified 2026-11-11 against hisab 3.4.1 (`d6703c5`) on cyrius 6.7.18. The repro exits 1.
**Placement:** unpinned.
**Discovered:** 2026-11-11, by ushma's Rust → Cyrius port (M0 probe 10, `ushma/docs/development/port-parity.md` §3).
**Reporter:** ushma 2.0.0 port (`ushma/docs/development/ushma-port-plan.md` §2.7, §8)
**Severity:** Low. Nothing is miscompiled: the two copies are byte-identical at the pin, and the stdlib one wins. But
it is a warning on every resolve that only hisab can remove.
**Affects:** hisab 3.4.1, `cyrius.cyml:78–81`.

## Summary

hisab's manifest declares a sakshi package:

```toml
[deps.sakshi]
git = "https://github.com/MacCracken/sakshi.git"
tag = "2.5.10"
modules = ["dist/sakshi.cyr"]
```

`dist/hisab.cyr` makes no `sakshi_` call (`grep -c sakshi_ dist/hisab.cyr` is 0). A consumer that also needs the
sakshi **stdlib leaf** gets this on every resolve:

```
warning: refusing to overwrite stdlib leaf 'sakshi' -- a dependency resolves a package of the same name.
  kept (stdlib snapshot): lib/sakshi.cyr
  skipped (dep artifact): ~/.cyrius/deps/sakshi/2.5.10/dist/sakshi.cyr
  declare it in [deps] stdlib OR depend on the package, not both.
```

That covers anything using `lib/log.cyr` or sandhi, both of which need the leaf. The consumer can't follow the advice,
because it never declared the package: hisab did.

The warning also travels. ushma depends on hisab in production (`debye_cv`'s Simpson, the steam tables' lerp), and
every resolve warned in M0 probe 10:
- ushma's own `cyrius deps`, `cyrius distlib --all` and `cyrius distlib --check`;
- the own `cyrius deps` of each of four consumer projects of ushma's bundles.

The three copies are byte-identical (sha256 `8a819554…`): `lib/sakshi.cyr`, the 6.7.18 snapshot's and the 2.5.10
package's.

prakash avoided the warning by making hisab `scope = "test"`. ushma can't: hisab is a production dependency there.

## Reproduction

`repros/2026-11-11-deps-sakshi-package-warns-every-consumer.sh` makes a consumer with `stdlib = [..., "sakshi",
"log"]` and `[deps.hisab]` at tag 3.4.1, runs `cyrius deps`, and exits 1 while the warning appears.

## Root cause

`cyrius.cyml:78–81` declares a dependency the bundle doesn't use. The resolver's refusal
(`cbt/deps/stdlib.cyr:251–259`) is correct.

## Proposed fix

Delete the `[deps.sakshi]` block. Nothing in hisab includes sakshi. patra did the same at 1.13.0 ("zero `[deps.*]`
blocks — sakshi moves to the stdlib"). If a later hisab module needs sakshi, declare it in `[deps] stdlib` instead.

## Consumer-side workaround

None inside a consumer. ushma's gates filter the line by its exact text (`refusing to overwrite stdlib leaf
.sakshi.`), and `docs/development/state.md` records it as a known upstream diagnostic.
