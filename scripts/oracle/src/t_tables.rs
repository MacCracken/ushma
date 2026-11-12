//! The public data tables — `src/material_data.cyr`, `src/state_data.cyr`, `src/phase_data.cyr`,
//! `src/chem_data.cyr` (port plan §4 step 6, §7.3) — and their golden rows.
//!
//! Each Rust table is a set of named `pub const` presets of one struct type plus an `ALL_*` slice
//! of references to them (`material::COPPER`, …, `material::ALL_MATERIALS`). In Cyrius it is:
//! - `US_<TABLE>_COUNT`, the `ALL_*` length;
//! - one index constant per preset, `US_<TABLE>_<PRESET>`, its row in `ALL_*` order;
//! - one typed column array per f64 field, `US_<TABLE>_<FIELD>: f64[n]`;
//! - the names as a NUL-separated blob with offsets (`_US_<TABLE>_NAME_OFF` / `_BYTES`).
//!
//! Every value comes from the oracle's API; the Rust source is read for comments and for the
//! completeness proof: every `pub const` of the table's struct type is a listed preset, and
//! `ALL_*` is a permutation of the presets (matched by name and every field's bits). The golden
//! rows (`<module>.tsv`): one `<module>::<PRESET>` row per preset, its name and fields in
//! declaration order, and one `<module>::<ALL>` row naming the presets in `ALL_*` order.

use crate::cyr::{doc_comment, doc_comment_after, f64_column, header, name_blob};
use crate::rustsrc::const_items;
use crate::tsv::{Class, Module, Tok, f, s};
use std::collections::BTreeSet;
use std::path::Path;

/// One named preset: its Rust const name, its `name` string, its f64 fields in declaration order.
pub struct Preset {
    pub rust: &'static str,
    pub name: String,
    pub fields: Vec<f64>,
}

/// One public table.
pub struct Table {
    /// The generated file, `src/<file>.cyr`.
    pub file: &'static str,
    pub title: &'static str,
    pub module: &'static str,
    pub ty: &'static str,
    pub prefix: &'static str,
    pub all: &'static str,
    pub fields: &'static [&'static str],
    /// Presets in source order.
    pub presets: Vec<Preset>,
    /// `ALL_*` in its order, as (name, fields).
    pub all_rows: Vec<(String, Vec<f64>)>,
    /// A prefix the presets' Rust names carry that the index constant drops (`GAS_`).
    pub strip: &'static str,
}

impl Table {
    /// `ALL_*` order as indices into `presets`, after proving it is a permutation of them.
    fn order(&self) -> Result<Vec<usize>, String> {
        let mut out = Vec::new();
        let mut seen = BTreeSet::new();
        for (name, fields) in &self.all_rows {
            let k = self
                .presets
                .iter()
                .position(|p| {
                    &p.name == name
                        && p.fields.len() == fields.len()
                        && p.fields
                            .iter()
                            .zip(fields)
                            .all(|(a, b)| a.to_bits() == b.to_bits())
                })
                .ok_or_else(|| {
                    format!(
                        "{}::{}: entry `{name}` is no listed preset",
                        self.module, self.all
                    )
                })?;
            if !seen.insert(k) {
                return Err(format!("{}::{}: `{name}` twice", self.module, self.all));
            }
            out.push(k);
        }
        if out.len() != self.presets.len() {
            return Err(format!(
                "{}::{}: {} entries, {} presets",
                self.module,
                self.all,
                out.len(),
                self.presets.len()
            ));
        }
        Ok(out)
    }

    fn index_const(&self, p: &Preset) -> String {
        let bare = p.rust.strip_prefix(self.strip).unwrap_or(p.rust);
        format!("US_{}_{bare}", self.prefix)
    }

    /// Every `pub const` of this table's type in the module is a listed preset, and the `ALL_*`
    /// slice exists.
    fn complete(&self, src: &str) -> Result<(), String> {
        let items = const_items(self.module, src)?;
        let found: BTreeSet<&str> = items
            .iter()
            .filter(|c| c.public && c.ty == self.ty)
            .map(|c| c.name.as_str())
            .collect();
        let listed: BTreeSet<&str> = self.presets.iter().map(|p| p.rust).collect();
        if found != listed {
            return Err(format!(
                "{}: presets in rust-old {found:?}, listed {listed:?}",
                self.module
            ));
        }
        let all_ty = format!("&[&{}]", self.ty);
        if !items
            .iter()
            .any(|c| c.public && c.name == self.all && c.ty == all_ty)
        {
            return Err(format!(
                "{}: no `pub const {}: {all_ty}`",
                self.module, self.all
            ));
        }
        Ok(())
    }

    pub fn emit_cyr(&self, root: &Path, provenance: &str) -> Result<String, String> {
        let src_path = root.join(format!("rust-old/src/{}.rs", self.module));
        let src = std::fs::read_to_string(&src_path)
            .map_err(|e| format!("{}: {e}", src_path.display()))?;
        self.complete(&src)?;
        let order = self.order()?;
        let mut out = header(
            self.title,
            provenance,
            &[&format!(
                "Rows are in `{}::{}` order; each preset's row is its index constant.\n\
                 f64 cells are the oracle's values as round-trip decimals, and\n\
                 tests/golden/{}.tsv pins every cell (one `{}::<PRESET>` row per preset).",
                self.module, self.all, self.module, self.module
            )],
        );
        out.push_str(&doc_comment(&src, &format!("pub struct {}", self.ty))?);
        out.push_str(&doc_comment(&src, &format!("pub const {}:", self.all))?);
        out.push_str(&format!(
            "# Rust: ushma::{}::{} (its length).\n",
            self.module, self.all
        ));
        out.push_str(&format!(
            "const US_{}_COUNT = {};\n\n",
            self.prefix,
            order.len()
        ));
        for (row, &k) in order.iter().enumerate() {
            let p = &self.presets[k];
            out.push_str(&doc_comment(&src, &format!("pub const {}:", p.rust))?);
            out.push_str(&format!("# Rust: ushma::{}::{}\n", self.module, p.rust));
            out.push_str(&format!("const {} = {row};\n\n", self.index_const(p)));
        }
        let anchor = format!("pub struct {} {{", self.ty);
        out.push_str(&doc_comment_after(&src, &anchor, "pub name:")?);
        let names: Vec<String> = order
            .iter()
            .map(|&k| self.presets[k].name.clone())
            .collect();
        out.push_str(&name_blob(
            &format!("_US_{}", self.prefix),
            &names,
            &format!("`{}::name`", self.ty),
        )?);
        out.push('\n');
        let labels: Vec<String> = order
            .iter()
            .map(|&k| self.presets[k].rust.to_string())
            .collect();
        for (j, field) in self.fields.iter().enumerate() {
            out.push_str(&doc_comment_after(&src, &anchor, &format!("pub {field}:"))?);
            out.push_str(&format!(
                "# Indexed by {} row; mirrors `{}::{field}`.\n",
                self.all, self.ty
            ));
            let values: Vec<f64> = order.iter().map(|&k| self.presets[k].fields[j]).collect();
            out.push_str(&f64_column(
                &format!("US_{}_{}", self.prefix, field.to_uppercase()),
                &values,
                &labels,
            )?);
            out.push('\n');
        }
        Ok(out.trim_end().to_string() + "\n")
    }

    /// This table's golden rows, into its module's file.
    pub fn rows(&self, m: &mut Module) -> Result<(), String> {
        let order = self.order()?;
        for p in &self.presets {
            let mut outs: Vec<Tok> = vec![s(p.name.clone())];
            outs.extend(p.fields.iter().map(|&x| f(x)));
            m.row(
                &format!("{}::{}", self.module, p.rust),
                Class::A,
                &[],
                &outs,
            );
        }
        let mut outs = vec![Tok::V(order.len())];
        outs.extend(order.iter().map(|&k| s(self.presets[k].rust)));
        m.row(
            &format!("{}::{}", self.module, self.all),
            Class::A,
            &[],
            &outs,
        );
        Ok(())
    }
}

macro_rules! presets {
    ($module:ident, $ty:ident, [$($c:ident),* $(,)?], |$v:ident| $fields:expr) => {{
        fn fields($v: &ushma::$module::$ty) -> Vec<f64> { $fields }
        let presets = vec![$( Preset {
            rust: stringify!($c),
            name: ushma::$module::$c.name.to_string(),
            fields: fields(&ushma::$module::$c),
        } ),*];
        (presets, fields as fn(&ushma::$module::$ty) -> Vec<f64>)
    }};
}

/// The four tables, in `[lib]` order.
pub fn all() -> Vec<Table> {
    let (mat, mf) = presets!(
        material,
        ThermalMaterial,
        [
            COPPER, ALUMINUM, IRON, STEEL, GLASS, WATER, AIR, WOOD_OAK, CONCRETE, DIAMOND
        ],
        |m| vec![
            m.conductivity,
            m.specific_heat,
            m.density,
            m.melting_point,
            m.boiling_point
        ]
    );
    let (gas, gf) = presets!(
        state,
        GasData,
        [
            GAS_NITROGEN,
            GAS_OXYGEN,
            GAS_CO2,
            GAS_METHANE,
            GAS_WATER,
            GAS_AMMONIA,
            GAS_ETHANE,
            GAS_PROPANE
        ],
        |g| vec![g.critical_t, g.critical_p, g.acentric_factor]
    );
    let (sub, sf) = presets!(
        phase,
        SubstanceData,
        [
            WATER_PHASE,
            NITROGEN,
            OXYGEN,
            CARBON_DIOXIDE,
            AMMONIA,
            METHANE,
            ETHANOL,
            MERCURY
        ],
        |d| vec![
            d.triple_point_t,
            d.triple_point_p,
            d.critical_t,
            d.critical_p,
            d.critical_density,
            d.latent_heat_fusion,
            d.latent_heat_vaporization,
            d.melting_point,
            d.boiling_point,
            d.molar_mass
        ]
    );
    let (spc, cf) = presets!(
        chem,
        Species,
        [
            H2, O2, N2, H2O_GAS, H2O_LIQUID, CO2, CO, CH4, C2H6, C3H8, NH3, NO
        ],
        |c| vec![c.delta_hf, c.delta_gf, c.cp]
    );
    vec![
        Table {
            file: "material_data",
            title: "material data: 10 reference materials at ~300 K",
            module: "material",
            ty: "ThermalMaterial",
            prefix: "MATERIAL",
            all: "ALL_MATERIALS",
            fields: &[
                "conductivity",
                "specific_heat",
                "density",
                "melting_point",
                "boiling_point",
            ],
            all_rows: ushma::material::ALL_MATERIALS
                .iter()
                .map(|m| (m.name.to_string(), mf(m)))
                .collect(),
            presets: mat,
            strip: "",
        },
        Table {
            file: "state_data",
            title: "gas data: critical constants and acentric factors of 8 gases",
            module: "state",
            ty: "GasData",
            prefix: "GAS",
            all: "ALL_GASES",
            fields: &["critical_t", "critical_p", "acentric_factor"],
            all_rows: ushma::state::ALL_GASES
                .iter()
                .map(|g| (g.name.to_string(), gf(g)))
                .collect(),
            presets: gas,
            strip: "GAS_",
        },
        Table {
            file: "phase_data",
            title: "phase data: triple and critical points, latent heats of 8 substances",
            module: "phase",
            ty: "SubstanceData",
            prefix: "SUBSTANCE",
            all: "ALL_SUBSTANCES",
            fields: &[
                "triple_point_t",
                "triple_point_p",
                "critical_t",
                "critical_p",
                "critical_density",
                "latent_heat_fusion",
                "latent_heat_vaporization",
                "melting_point",
                "boiling_point",
                "molar_mass",
            ],
            all_rows: ushma::phase::ALL_SUBSTANCES
                .iter()
                .map(|d| (d.name.to_string(), sf(d)))
                .collect(),
            presets: sub,
            strip: "",
        },
        Table {
            file: "chem_data",
            title: "chemical data: standard formation enthalpies and Gibbs energies of 12 species",
            module: "chem",
            ty: "Species",
            prefix: "SPECIES",
            all: "ALL_SPECIES",
            fields: &["delta_hf", "delta_gf", "cp"],
            all_rows: ushma::chem::ALL_SPECIES
                .iter()
                .map(|c| (c.name.to_string(), cf(c)))
                .collect(),
            presets: spc,
            strip: "",
        },
    ]
}

/// The golden modules these tables seed: `chem`, `material`, `phase`, `state`. Each module's
/// function rows are added to the same module at its milestone (M1–M5).
pub fn golden() -> Result<Vec<Module>, String> {
    let mut out = Vec::new();
    for t in all() {
        let mut m = Module::new(t.module);
        t.rows(&mut m)?;
        out.push(m);
    }
    Ok(out)
}
