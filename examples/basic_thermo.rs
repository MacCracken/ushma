//! Basic thermodynamics demonstration.
//!
//! Run with: `cargo run --example basic_thermo --features transfer,state,material`

use ushma::material::{self, COPPER, WATER, ALL_MATERIALS};
use ushma::state;
use ushma::transfer;

fn main() {
    println!("=== Ushma — Thermodynamics Demo ===\n");

    // --- Heat Conduction ---
    println!("--- Heat Conduction (Fourier's Law) ---");
    let q = transfer::conduction(
        COPPER.conductivity,
        0.01,    // 1 cm² = 0.01 m²
        373.15,  // 100°C
        293.15,  // 20°C
        0.05,    // 5 cm thick
    )
    .unwrap();
    println!(
        "Copper wall (5cm, 1cm²): {:.1} W for ΔT = 80 K",
        q
    );

    // --- Convection ---
    println!("\n--- Convection (Newton's Cooling) ---");
    let q_conv = transfer::convection(25.0, 1.0, 353.15, 293.15);
    println!(
        "Forced convection (h=25, A=1m²): {:.1} W for ΔT = 60 K",
        q_conv
    );

    // --- Radiation ---
    println!("\n--- Radiation (Stefan-Boltzmann) ---");
    let q_rad = transfer::radiation(0.9, 1.0, 473.15, 293.15).unwrap();
    println!(
        "Hot surface (ε=0.9, A=1m², 200°C): {:.1} W radiated",
        q_rad
    );

    // --- Ideal Gas Law ---
    println!("\n--- Ideal Gas Law ---");
    let v = state::ideal_gas_volume(1.0, state::STANDARD_TEMP, state::ATM).unwrap();
    println!(
        "1 mol at STP: V = {:.4} m³ ({:.2} L)",
        v,
        v * 1000.0
    );

    let p = state::ideal_gas_pressure(2.0, 400.0, 0.05).unwrap();
    println!(
        "2 mol at 400K in 50L: P = {:.0} Pa ({:.2} atm)",
        p,
        p / state::ATM
    );

    // --- Thermal Materials ---
    println!("\n--- Material Properties (at ~300 K) ---");
    println!(
        "{:<20} {:>10} {:>12} {:>10}",
        "Material", "k (W/m⋅K)", "c_p (J/kg⋅K)", "α (m²/s)"
    );
    println!("{}", "-".repeat(56));
    for mat in ALL_MATERIALS {
        println!(
            "{:<20} {:>10.3} {:>12.1} {:>10.2e}",
            mat.name,
            mat.conductivity,
            mat.specific_heat,
            mat.diffusivity(),
        );
    }

    // --- Lumped Capacitance Cooling ---
    println!("\n--- Cooling Curve (lumped capacitance) ---");
    let t0 = 373.15;
    let t_env = 293.15;
    let tau = 120.0; // 2 minute time constant
    for &t in &[0.0, 30.0, 60.0, 120.0, 300.0, 600.0] {
        let temp = transfer::lumped_capacitance(t0, t_env, t, tau);
        println!("  t = {:>4.0}s → T = {:.1} K ({:.1} °C)", t, temp, temp - 273.15);
    }

    println!("\n=== Done ===");
}
