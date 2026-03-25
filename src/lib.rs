//! Ushma — Thermodynamics simulation for AGNOS
//!
//! Sanskrit: ऊष्मा (ushma) — heat, thermal energy
//!
//! Provides heat transfer simulation, entropy calculations, equations of state,
//! and thermal material properties. Built on [hisab](https://crates.io/crates/hisab)
//! for ODE/PDE solvers and [impetus](https://crates.io/crates/impetus) for body contacts.
//!
//! # Modules
//!
//! - [`transfer`] — Conduction, convection, radiation heat transfer
//! - [`state`] — Equations of state, ideal gas, phase diagrams
//! - [`entropy`] — Entropy, free energy, thermodynamic potentials
//! - [`material`] — Thermal properties, specific heat, conductivity tables
//! - [`phase`] — Phase transitions, Clausius-Clapeyron, substance data
//! - [`steam`] — Saturated and superheated steam tables
//! - [`cycle`] — Thermodynamic cycles (Otto, Diesel, Brayton, Rankine)
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

#[cfg(feature = "logging")]
pub mod logging;

#[cfg(feature = "ai")]
pub mod ai;

pub use error::UshmaError;
