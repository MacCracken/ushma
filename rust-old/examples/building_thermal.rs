//! Building thermal model — walls, windows, insulation with 1D transient conduction.
//!
//! Run with: `cargo run --example building_thermal --features transfer,material,numerical`

use ushma::material;
use ushma::numerical::{BoundaryCondition, ThermalGrid1D};
use ushma::transfer;

fn main() {
    println!("=== Ushma — Building Thermal Model ===\n");

    // --- Wall steady-state heat loss ---
    println!("--- Wall Heat Loss (steady state) ---");
    // Concrete wall: 20 cm thick, 10 m² area
    let thickness = 0.20;
    let area = 10.0;
    let t_inside = 293.15; // 20°C
    let t_outside = 268.15; // -5°C

    let q_concrete = transfer::conduction(
        material::CONCRETE.conductivity,
        area,
        t_inside,
        t_outside,
        thickness,
    )
    .unwrap();
    println!(
        "Concrete wall ({:.0} cm): {:.0} W heat loss",
        thickness * 100.0,
        q_concrete
    );

    // Add 5 cm insulation (wood)
    let r_concrete =
        transfer::thermal_resistance_conduction(material::CONCRETE.conductivity, area, thickness)
            .unwrap();
    let r_wood =
        transfer::thermal_resistance_conduction(material::WOOD_OAK.conductivity, area, 0.05)
            .unwrap();
    let r_total = transfer::thermal_resistance_series(&[r_concrete, r_wood]);
    let q_insulated = (t_inside - t_outside) / r_total;
    println!(
        "With 5 cm wood insulation: {:.0} W ({:.0}% reduction)",
        q_insulated,
        (1.0 - q_insulated / q_concrete) * 100.0
    );

    // --- Window convection loss ---
    println!("\n--- Window Convection ---");
    let h_inside = 10.0; // W/(m²·K), natural convection
    let h_outside = 25.0; // W/(m²·K), wind
    let window_area = 2.0;
    let r_conv_in = transfer::thermal_resistance_convection(h_inside, window_area).unwrap();
    let r_glass =
        transfer::thermal_resistance_conduction(material::GLASS.conductivity, window_area, 0.006)
            .unwrap();
    let r_conv_out = transfer::thermal_resistance_convection(h_outside, window_area).unwrap();
    let r_window = transfer::thermal_resistance_series(&[r_conv_in, r_glass, r_conv_out]);
    let q_window = (t_inside - t_outside) / r_window;
    println!("Single-pane window (2 m²): {:.0} W", q_window);

    // --- Transient cooling simulation ---
    println!("\n--- Transient Wall Cooling (1D FDM) ---");
    // Concrete wall suddenly exposed to cold outside
    let alpha = transfer::thermal_diffusivity(
        material::CONCRETE.conductivity,
        material::CONCRETE.density,
        material::CONCRETE.specific_heat,
    )
    .unwrap();
    let mut grid = ThermalGrid1D::new(
        21,
        thickness,
        alpha,
        t_inside,
        BoundaryCondition::Fixed(t_inside),
        BoundaryCondition::Fixed(t_outside),
    )
    .unwrap();

    let dt = 0.4 * grid.dx * grid.dx / grid.alpha;
    let steps_per_hour = (3600.0 / dt).ceil() as usize;

    for hour in 1..=4 {
        for _ in 0..steps_per_hour {
            grid.step_explicit(dt).unwrap();
        }
        let t_mid = grid.nodes[10];
        println!(
            "  After {hour} hr: wall midpoint T = {:.1} °C",
            t_mid - 273.15
        );
    }

    println!("\n=== Done ===");
}
