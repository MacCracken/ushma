//! `constants.tsv` — every public scalar const of ushma 1.4.0, class A.
//!
//! One row per constant: `const A s:<module>::<NAME> f:<bits>` (or `i:<n>` for the one integer,
//! `numerical::MAX_ADAPTIVE_SUBSTEPS`). The list is shared with the `src/constants.cyr`
//! generator (`t_constants::all`), so the two can't drift apart.

use crate::t_constants::{self, Value};
use crate::tsv::{Class, Module, f, i, s};

pub fn emit() -> Module {
    let mut m = Module::new("constants");
    for c in t_constants::all() {
        let out = match c.value {
            Value::F(x) => f(x),
            Value::I(n) => i(n),
        };
        m.row("const", Class::A, &[s(c.path())], &[out]);
    }
    m
}
