//! Cooling system design — fins, heat exchangers, convection.
//!
//! Run with: `cargo run --example cooling_system --features transfer,material`

use ushma::material;
use ushma::transfer;

fn main() {
    println!("=== Ushma — Cooling System Design ===\n");

    // --- Finned heat sink for electronics ---
    println!("--- Aluminum Fin Array ---");
    let h = 25.0; // W/(m²·K), forced air
    let k = material::ALUMINUM.conductivity; // 237 W/(m·K)
    let fin_thickness = 0.002; // 2 mm
    let fin_height = 0.03; // 30 mm
    let fin_width = 0.05; // 50 mm
    let perimeter = 2.0 * (fin_thickness + fin_width);
    let cross_area = fin_thickness * fin_width;
    let t_base = 353.15; // 80°C (chip temperature)
    let t_air = 298.15; // 25°C

    let q_fin =
        transfer::fin_rectangular_heat(h, perimeter, k, cross_area, fin_height, t_base, t_air)
            .unwrap();
    let m = transfer::fin_parameter(h, perimeter, k, cross_area).unwrap();
    let eta = transfer::fin_efficiency_rectangular(m, fin_height).unwrap();
    let effectiveness = transfer::fin_effectiveness(q_fin, h, cross_area, t_base, t_air).unwrap();

    println!(
        "  Fin: {:.0}mm × {:.0}mm × {:.1}mm Al",
        fin_width * 1000.0,
        fin_height * 1000.0,
        fin_thickness * 1000.0
    );
    println!("  Heat per fin: {:.2} W", q_fin);
    println!("  Fin efficiency: {:.1}%", eta * 100.0);
    println!("  Fin effectiveness: {:.1}×", effectiveness);

    let num_fins = 10;
    let q_total = q_fin * num_fins as f64;
    println!("  Total ({num_fins} fins): {:.1} W", q_total);

    // --- Heat exchanger sizing (counter-flow) ---
    println!("\n--- Counter-Flow Heat Exchanger ---");
    let t_h_in = 363.15; // 90°C hot coolant
    let t_h_out = 333.15; // 60°C
    let t_c_in = 293.15; // 20°C cooling water
    let t_c_out = 313.15; // 40°C

    let lmtd = transfer::lmtd_counter(t_h_in, t_h_out, t_c_in, t_c_out).unwrap();
    println!("  LMTD = {:.1} K", lmtd);

    // Target 5 kW heat rejection, U = 800 W/(m²·K)
    let q_target = 5000.0;
    let u = 800.0;
    let area_needed = q_target / (u * lmtd);
    println!("  U = {u} W/(m²·K)");
    println!(
        "  Required area for {:.0} W: {:.3} m²",
        q_target, area_needed
    );

    // Verify with ε-NTU
    let c_h = q_target / (t_h_in - t_h_out);
    let c_c = q_target / (t_c_out - t_c_in);
    let c_min = c_h.min(c_c);
    let c_max = c_h.max(c_c);
    let ntu_val = transfer::ntu(u, area_needed, c_min).unwrap();
    let eff = transfer::effectiveness_counter(ntu_val, c_min / c_max).unwrap();
    let q_check = transfer::heat_exchanger_ntu(eff, c_min, t_h_in, t_c_in);
    println!("  NTU = {:.2}, ε = {:.3}", ntu_val, eff);
    println!("  ε-NTU verification: Q = {:.0} W", q_check);

    // --- Dimensionless analysis ---
    println!("\n--- Dimensionless Numbers (air at 300 K) ---");
    let nu_air = 1.6e-5; // kinematic viscosity (m²/s)
    let alpha_air = 2.2e-5; // thermal diffusivity (m²/s)
    let velocity = 5.0; // m/s
    let length = 0.1; // 10 cm plate

    let re = transfer::reynolds_number(velocity, length, nu_air).unwrap();
    let pr = transfer::prandtl_number(nu_air, alpha_air).unwrap();
    let nu = transfer::nusselt_dittus_boelter(re, pr).unwrap();
    let h_calc = nu * 0.026 / length; // k_air ≈ 0.026 W/(m·K)

    println!("  Re = {:.0}", re);
    println!("  Pr = {:.2}", pr);
    println!("  Nu (Dittus-Boelter) = {:.1}", nu);
    println!("  h = {:.1} W/(m²·K)", h_calc);

    println!("\n=== Done ===");
}
