//! # Ushma
//!
//! > **Ushma** (Sanskrit: ऊष्मा — heat, thermal energy) — Thermodynamics simulation for AGNOS
//!
//! Provides heat transfer simulation, entropy calculations, equations of state,
//! and thermal material properties. Built on [hisab](https://crates.io/crates/hisab)
//! for math foundations.
//!
//! ## Feature flags
//!
//! | Feature | Default | Description |
//! |---------|---------|-------------|
//! | `transfer` | yes | Conduction, convection, radiation, fins, heat exchangers, view factors |
//! | `state` | yes | Ideal gas, van der Waals, Redlich-Kwong, Peng-Robinson, virial EOS |
//! | `entropy` | yes | Entropy, Carnot, Helmholtz/Gibbs free energy, entropy of mixing |
//! | `material` | yes | 10 built-in materials with thermal properties at ~300 K |
//! | `phase` | yes | Phase transitions, Clausius-Clapeyron, 8 substance data sets |
//! | `cycle` | yes | Otto, Diesel, Brayton, Rankine, refrigeration cycles |
//! | `numerical` | yes | 1D/2D finite difference, Crank-Nicolson, thermal networks |
//! | `chem` | yes | Hess's law, equilibrium constants, Van't Hoff, flame temperature |
//! | `stat` | yes | Boltzmann, partition functions, Maxwell-Boltzmann, Debye/Einstein |
//! | `steam` | no | Saturated + superheated steam tables (IAPWS-IF97) |
//! | `ai` | no | Daimon/hoosh agent integration |
//! | `logging` | no | Structured logging via tracing |
//! | `full` | — | Enables all features |
//!
//! ## Modules
//!
//! - [`transfer`] — Conduction, convection, radiation, fins, heat exchangers, dimensionless numbers
//! - [`state`] — Equations of state, real gas models, compressibility, mixture rules
//! - [`entropy`] — Entropy, free energy, thermodynamic potentials
//! - [`material`] — Thermal properties, specific heat, conductivity tables
//! - [`phase`] — Phase transitions, Clausius-Clapeyron, substance data
//! - [`steam`] — Saturated and superheated steam tables
//! - [`cycle`] — Thermodynamic cycles, diagram generation, efficiency comparison
//! - [`numerical`] — Finite difference solvers, thermal networks
//! - [`chem`] — Chemical thermodynamics, Hess's law, equilibrium
//! - [`stat`] — Statistical thermodynamics, Boltzmann, partition functions
//! - [`error`] — Error types

pub mod error;

#[cfg(feature = "transfer")]
pub mod transfer;

#[cfg(feature = "state")]
pub mod state;

#[cfg(feature = "entropy")]
pub mod entropy;

#[cfg(feature = "material")]
pub mod material;

#[cfg(feature = "phase")]
pub mod phase;

#[cfg(feature = "steam")]
pub mod steam;

#[cfg(feature = "cycle")]
pub mod cycle;

#[cfg(feature = "numerical")]
pub mod numerical;

#[cfg(feature = "chem")]
pub mod chem;

#[cfg(feature = "stat")]
pub mod stat;

#[cfg(feature = "logging")]
pub mod logging;

#[cfg(feature = "ai")]
pub mod ai;

pub use error::UshmaError;
