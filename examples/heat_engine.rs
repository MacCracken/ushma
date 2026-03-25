//! Heat engine cycle simulator — compares Otto, Diesel, and Brayton cycles.
//!
//! Run with: `cargo run --example heat_engine --features cycle,state,entropy`

use ushma::cycle;
use ushma::entropy;

fn main() {
    println!("=== Ushma — Heat Engine Cycle Comparison ===\n");

    // --- Otto cycle (gasoline engine) ---
    let otto = cycle::otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0).unwrap();
    println!("--- Otto Cycle (r=8, γ=1.4) ---");
    println!("  Efficiency: {:.1}%", otto.efficiency * 100.0);
    println!("  Net work:   {:.0} J/mol", otto.net_work);
    println!("  Heat in:    {:.0} J/mol", otto.heat_in);
    for (i, sp) in otto.state_points.iter().enumerate() {
        println!(
            "  State {}: T={:.0} K, P={:.0} Pa, V={:.6} m³",
            i + 1,
            sp.temperature,
            sp.pressure,
            sp.volume
        );
    }

    // --- Diesel cycle ---
    let diesel = cycle::diesel_cycle(300.0, 101_325.0, 20.0, 2.0, 1.4, 1.0).unwrap();
    println!("\n--- Diesel Cycle (r=20, rc=2) ---");
    println!("  Efficiency: {:.1}%", diesel.efficiency * 100.0);
    println!("  Net work:   {:.0} J/mol", diesel.net_work);

    // --- Brayton cycle (gas turbine) ---
    let brayton = cycle::brayton_cycle(300.0, 101_325.0, 10.0, 1400.0, 1.4, 1.0).unwrap();
    println!("\n--- Brayton Cycle (rp=10, T3=1400 K) ---");
    println!("  Efficiency:      {:.1}%", brayton.efficiency * 100.0);
    println!("  Back-work ratio: {:.1}%", brayton.back_work_ratio * 100.0);

    // --- Comparison vs Carnot ---
    println!("\n--- Comparison ---");
    println!(
        "{:<10} {:>8} {:>10} {:>8}",
        "Cycle", "η", "η_Carnot", "η_II"
    );
    println!("{}", "-".repeat(40));
    for (name, result) in [("Otto", &otto), ("Diesel", &diesel), ("Brayton", &brayton)] {
        let t_high = result.state_points[2].temperature;
        let t_low = result.state_points[0].temperature;
        let eta_carnot = entropy::carnot_efficiency(t_high, t_low).unwrap();
        println!(
            "{:<10} {:>7.1}% {:>9.1}% {:>7.1}%",
            name,
            result.efficiency * 100.0,
            eta_carnot * 100.0,
            result.efficiency / eta_carnot * 100.0
        );
    }

    println!("\n=== Done ===");
}
