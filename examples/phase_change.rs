//! Phase change material storage — latent heat, phase transitions, Clausius-Clapeyron.
//!
//! Run with: `cargo run --example phase_change --features phase,transfer`

use ushma::phase;
use ushma::transfer;

fn main() {
    println!("=== Ushma — Phase Change Material Storage ===\n");

    // --- Latent heat storage capacity ---
    println!("--- Water as Phase Change Material ---");
    let mass = 100.0; // 100 kg of water
    let q_fusion = phase::heat_of_fusion(&phase::WATER_PHASE, mass);
    let q_vaporization = phase::heat_of_vaporization(&phase::WATER_PHASE, mass);
    println!("  {:.0} kg water:", mass);
    println!("  Latent heat of fusion:       {:.0} kJ", q_fusion / 1000.0);
    println!(
        "  Latent heat of vaporization: {:.0} kJ",
        q_vaporization / 1000.0
    );

    // Compare: sensible heat for 10 K rise in liquid water
    let q_sensible = transfer::heat_stored(mass, 4186.0, 10.0);
    println!(
        "  Sensible heat (10 K rise):   {:.0} kJ",
        q_sensible / 1000.0
    );
    println!(
        "  Fusion stores {:.0}× more than 10 K sensible",
        q_fusion / q_sensible
    );

    // --- Ice to steam energy budget ---
    println!("\n--- Ice to Steam (1 kg, -20°C to 120°C) ---");
    let q_total =
        phase::heat_for_phase_change(1.0, 2500.0, 253.15, 393.15, &phase::WATER_PHASE).unwrap();
    println!("  Total energy: {:.0} kJ", q_total / 1000.0);

    // --- Phase diagram lookup ---
    println!("\n--- Phase Diagram ---");
    let conditions = [
        (200.0, 101_325.0, "Ice at -73°C"),
        (300.0, 101_325.0, "Liquid at 27°C"),
        (400.0, 101_325.0, "Steam at 127°C"),
        (700.0, 25_000_000.0, "Supercritical"),
    ];
    for (t, p, label) in conditions {
        let ph = phase::WATER_PHASE.phase_at(t, p).unwrap();
        println!("  {label}: T={t:.0} K, P={p:.0} Pa → {ph:?}");
    }

    // --- Clausius-Clapeyron: boiling point at altitude ---
    println!("\n--- Boiling Point vs Pressure ---");
    let l_molar = phase::WATER_PHASE.latent_heat_vaporization * phase::WATER_PHASE.molar_mass;
    let pressures = [101_325.0, 80_000.0, 60_000.0, 40_000.0];
    for &p in &pressures {
        // Invert CC: find T where P_sat = p
        // Use bisection between 330 K and 400 K
        let mut t_lo = 330.0;
        let mut t_hi = 400.0;
        for _ in 0..50 {
            let t_mid = 0.5 * (t_lo + t_hi);
            let p_sat =
                phase::clausius_clapeyron_pressure(101_325.0, 373.15, t_mid, l_molar).unwrap();
            if p_sat < p {
                t_lo = t_mid;
            } else {
                t_hi = t_mid;
            }
        }
        let t_boil = 0.5 * (t_lo + t_hi);
        println!(
            "  P={:.0} kPa → T_boil={:.1} °C",
            p / 1000.0,
            t_boil - 273.15
        );
    }

    // --- Compare substances ---
    println!("\n--- Substance Comparison ---");
    println!(
        "{:<18} {:>10} {:>10} {:>10}",
        "Substance", "T_boil(K)", "T_crit(K)", "L_vap(kJ/kg)"
    );
    println!("{}", "-".repeat(52));
    for sub in phase::ALL_SUBSTANCES {
        println!(
            "{:<18} {:>10.1} {:>10.1} {:>10.0}",
            sub.name,
            sub.boiling_point,
            sub.critical_t,
            sub.latent_heat_vaporization / 1000.0
        );
    }

    println!("\n=== Done ===");
}
