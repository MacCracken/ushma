//! Thermodynamic cycles — Otto, Diesel, Brayton, Rankine, refrigeration.
//!
//! All SI units: kelvins, pascals, joules, moles.
//! Ideal-gas cycles use the calorically perfect gas model.
//! Steam cycles (behind `steam` feature) use IAPWS-IF97 table data.

use serde::{Deserialize, Serialize};

use crate::error::{Result, UshmaError};
use crate::state::GAS_CONSTANT;

/// A thermodynamic state point in a cycle.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct StatePoint {
    /// Temperature (K).
    pub temperature: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Volume (m³). Specific volume for steam cycles, total for ideal-gas.
    pub volume: f64,
    /// Entropy (J/K or J/(kg·K) for steam).
    pub entropy: f64,
    /// Enthalpy (J or J/kg for steam).
    pub enthalpy: f64,
}

/// Type of thermodynamic process connecting two state points.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum ProcessKind {
    /// Constant entropy (adiabatic reversible).
    Isentropic,
    /// Constant volume.
    Isochoric,
    /// Constant pressure.
    Isobaric,
    /// Constant temperature.
    Isothermal,
    /// Constant enthalpy (throttling).
    Isenthalpic,
}

/// Kind of thermodynamic cycle.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum CycleKind {
    Otto,
    Diesel,
    Brayton,
    Rankine,
    Refrigeration,
}

/// Result of a thermodynamic cycle analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CycleResult {
    /// Cycle type.
    pub kind: CycleKind,
    /// State points (one per vertex of the cycle).
    pub state_points: Vec<StatePoint>,
    /// Process types between consecutive state points (wraps around).
    pub processes: Vec<ProcessKind>,
    /// Heat added to the working fluid (J).
    pub heat_in: f64,
    /// Heat rejected from the working fluid (J).
    pub heat_out: f64,
    /// Net work output (J). Positive = work produced.
    pub net_work: f64,
    /// Thermal efficiency (dimensionless, 0-1).
    pub efficiency: f64,
    /// Back-work ratio: W_compressor / W_expander (dimensionless).
    pub back_work_ratio: f64,
    /// Heat capacity ratio γ = Cp/Cv (for ideal-gas cycles).
    pub gamma: f64,
}

/// A point for plotting cycle diagrams (T-s or P-v).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct DiagramPoint {
    /// X-axis value.
    pub x: f64,
    /// Y-axis value.
    pub y: f64,
}

// --- Ideal-gas process helpers ---

/// Heat transfer at constant volume: Q = n·Cv·(T₂ - T₁) (J).
#[inline]
#[must_use]
pub fn isochoric_heat(moles: f64, cv: f64, t1: f64, t2: f64) -> f64 {
    moles * cv * (t2 - t1)
}

/// Heat transfer at constant pressure: Q = n·Cp·(T₂ - T₁) (J).
#[inline]
#[must_use]
pub fn isobaric_heat(moles: f64, cp: f64, t1: f64, t2: f64) -> f64 {
    moles * cp * (t2 - t1)
}

/// Work done in an adiabatic process: W = n·Cv·(T₁ - T₂) (J).
///
/// Positive work = expansion (T drops). Negative = compression (T rises).
#[inline]
#[must_use]
pub fn adiabatic_work(moles: f64, cv: f64, t1: f64, t2: f64) -> f64 {
    moles * cv * (t1 - t2)
}

/// Adiabatic temperature from pressure ratio: T₂ = T₁·(P₂/P₁)^((γ-1)/γ).
///
/// For isentropic compression/expansion in Brayton cycles.
pub fn adiabatic_temp_from_pressure(t1: f64, p1: f64, p2: f64, gamma: f64) -> Result<f64> {
    if t1 <= 0.0 {
        return Err(UshmaError::InvalidTemperature { kelvin: t1 });
    }
    if p1 <= 0.0 {
        return Err(UshmaError::InvalidPressure { pascals: p1 });
    }
    if p2 <= 0.0 {
        return Err(UshmaError::InvalidPressure { pascals: p2 });
    }
    if gamma <= 1.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("gamma {gamma} must be > 1"),
        });
    }
    Ok(t1 * (p2 / p1).powf((gamma - 1.0) / gamma))
}

// --- Cycle analyses ---

/// Otto cycle analysis (spark-ignition engine).
///
/// Four processes on an ideal gas:
/// 1→2: Isentropic compression
/// 2→3: Isochoric heat addition
/// 3→4: Isentropic expansion
/// 4→1: Isochoric heat rejection
///
/// - `t1`: initial temperature (K)
/// - `p1`: initial pressure (Pa)
/// - `compression_ratio`: r = V₁/V₂ (must be > 1)
/// - `heat_in`: heat added per mole in process 2→3 (J/mol)
/// - `gamma`: heat capacity ratio Cp/Cv (must be > 1)
/// - `moles`: amount of gas (mol)
pub fn otto_cycle(
    t1: f64,
    p1: f64,
    compression_ratio: f64,
    heat_in: f64,
    gamma: f64,
    moles: f64,
) -> Result<CycleResult> {
    if t1 <= 0.0 {
        return Err(UshmaError::InvalidTemperature { kelvin: t1 });
    }
    if p1 <= 0.0 {
        return Err(UshmaError::InvalidPressure { pascals: p1 });
    }
    if compression_ratio <= 1.0 {
        return Err(UshmaError::InvalidCompressionRatio {
            ratio: compression_ratio,
        });
    }
    if gamma <= 1.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("gamma {gamma} must be > 1"),
        });
    }
    if moles <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("moles {moles} must be positive"),
        });
    }

    if heat_in <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("heat_in {heat_in} J/mol must be positive"),
        });
    }

    let r = compression_ratio;
    let cv = GAS_CONSTANT / (gamma - 1.0);
    let cp = gamma * cv;
    let n = moles;

    // State 1: given
    let v1 = n * GAS_CONSTANT * t1 / p1;

    // State 2: isentropic compression
    let v2 = v1 / r;
    let t2 = t1 * r.powf(gamma - 1.0);
    let p2 = n * GAS_CONSTANT * t2 / v2;

    // State 3: isochoric heat addition
    let v3 = v2;
    let q_in = heat_in * n;
    let t3 = t2 + q_in / (n * cv);
    let p3 = n * GAS_CONSTANT * t3 / v3;

    // State 4: isentropic expansion
    let v4 = v1;
    let t4 = t3 * (1.0 / r).powf(gamma - 1.0);
    let p4 = n * GAS_CONSTANT * t4 / v4;

    // Energy balance
    let q_out = n * cv * (t4 - t1);
    let w_net = q_in - q_out;
    let efficiency = w_net / q_in;

    // Entropy at each state (relative to state 1 = 0)
    let s1 = 0.0;
    let s2 = s1; // isentropic
    let s3 = s2 + n * cv * (t3 / t2).ln(); // isochoric: ds = nCv ln(T3/T2)
    let s4 = s3; // isentropic

    Ok(CycleResult {
        kind: CycleKind::Otto,
        state_points: vec![
            StatePoint {
                temperature: t1,
                pressure: p1,
                volume: v1,
                entropy: s1,
                enthalpy: n * cp * t1,
            },
            StatePoint {
                temperature: t2,
                pressure: p2,
                volume: v2,
                entropy: s2,
                enthalpy: n * cp * t2,
            },
            StatePoint {
                temperature: t3,
                pressure: p3,
                volume: v3,
                entropy: s3,
                enthalpy: n * cp * t3,
            },
            StatePoint {
                temperature: t4,
                pressure: p4,
                volume: v4,
                entropy: s4,
                enthalpy: n * cp * t4,
            },
        ],
        processes: vec![
            ProcessKind::Isentropic,
            ProcessKind::Isochoric,
            ProcessKind::Isentropic,
            ProcessKind::Isochoric,
        ],
        heat_in: q_in,
        heat_out: q_out,
        net_work: w_net,
        efficiency,
        back_work_ratio: 0.0, // no separate compressor/expander in reciprocating engine
        gamma,
    })
}

/// Diesel cycle analysis (compression-ignition engine).
///
/// Four processes on an ideal gas:
/// 1→2: Isentropic compression
/// 2→3: Isobaric heat addition
/// 3→4: Isentropic expansion
/// 4→1: Isochoric heat rejection
///
/// - `t1`: initial temperature (K)
/// - `p1`: initial pressure (Pa)
/// - `compression_ratio`: r = V₁/V₂ (must be > 1)
/// - `cutoff_ratio`: rc = V₃/V₂ (must be > 1 and < r)
/// - `gamma`: heat capacity ratio Cp/Cv (must be > 1)
/// - `moles`: amount of gas (mol)
pub fn diesel_cycle(
    t1: f64,
    p1: f64,
    compression_ratio: f64,
    cutoff_ratio: f64,
    gamma: f64,
    moles: f64,
) -> Result<CycleResult> {
    if t1 <= 0.0 {
        return Err(UshmaError::InvalidTemperature { kelvin: t1 });
    }
    if p1 <= 0.0 {
        return Err(UshmaError::InvalidPressure { pascals: p1 });
    }
    if compression_ratio <= 1.0 {
        return Err(UshmaError::InvalidCompressionRatio {
            ratio: compression_ratio,
        });
    }
    if cutoff_ratio <= 1.0 {
        return Err(UshmaError::InvalidCutoffRatio {
            ratio: cutoff_ratio,
        });
    }
    if cutoff_ratio >= compression_ratio {
        return Err(UshmaError::InvalidCycleParameter {
            reason: format!(
                "cutoff ratio {cutoff_ratio} must be less than compression ratio {compression_ratio}"
            ),
        });
    }
    if gamma <= 1.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("gamma {gamma} must be > 1"),
        });
    }
    if moles <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("moles {moles} must be positive"),
        });
    }

    let r = compression_ratio;
    let rc = cutoff_ratio;
    let cv = GAS_CONSTANT / (gamma - 1.0);
    let cp = gamma * cv;
    let n = moles;

    // State 1
    let v1 = n * GAS_CONSTANT * t1 / p1;

    // State 2: isentropic compression
    let v2 = v1 / r;
    let t2 = t1 * r.powf(gamma - 1.0);
    let p2 = n * GAS_CONSTANT * t2 / v2;

    // State 3: isobaric heat addition
    let p3 = p2;
    let v3 = v2 * rc;
    let t3 = t2 * rc;

    // State 4: isentropic expansion
    let v4 = v1;
    let t4 = t3 * (v3 / v4).powf(gamma - 1.0);
    let p4 = n * GAS_CONSTANT * t4 / v4;

    let q_in = n * cp * (t3 - t2);
    let q_out = n * cv * (t4 - t1);
    let w_net = q_in - q_out;
    let efficiency = w_net / q_in;

    let s1 = 0.0;
    let s2 = s1;
    let s3 = s2 + n * cp * (t3 / t2).ln(); // isobaric: ds = nCp ln(T3/T2)
    let s4 = s3;

    Ok(CycleResult {
        kind: CycleKind::Diesel,
        state_points: vec![
            StatePoint {
                temperature: t1,
                pressure: p1,
                volume: v1,
                entropy: s1,
                enthalpy: n * cp * t1,
            },
            StatePoint {
                temperature: t2,
                pressure: p2,
                volume: v2,
                entropy: s2,
                enthalpy: n * cp * t2,
            },
            StatePoint {
                temperature: t3,
                pressure: p3,
                volume: v3,
                entropy: s3,
                enthalpy: n * cp * t3,
            },
            StatePoint {
                temperature: t4,
                pressure: p4,
                volume: v4,
                entropy: s4,
                enthalpy: n * cp * t4,
            },
        ],
        processes: vec![
            ProcessKind::Isentropic,
            ProcessKind::Isobaric,
            ProcessKind::Isentropic,
            ProcessKind::Isochoric,
        ],
        heat_in: q_in,
        heat_out: q_out,
        net_work: w_net,
        efficiency,
        back_work_ratio: 0.0,
        gamma,
    })
}

/// Brayton cycle analysis (gas turbine).
///
/// Four processes on an ideal gas:
/// 1→2: Isentropic compression
/// 2→3: Isobaric heat addition
/// 3→4: Isentropic expansion
/// 4→1: Isobaric heat rejection
///
/// - `t1`: compressor inlet temperature (K)
/// - `p1`: compressor inlet pressure (Pa)
/// - `pressure_ratio`: rp = P₂/P₁ (must be > 1)
/// - `t3`: turbine inlet temperature (K), must be > T₂
/// - `gamma`: heat capacity ratio Cp/Cv (must be > 1)
/// - `moles`: amount of gas (mol)
pub fn brayton_cycle(
    t1: f64,
    p1: f64,
    pressure_ratio: f64,
    t3: f64,
    gamma: f64,
    moles: f64,
) -> Result<CycleResult> {
    if t1 <= 0.0 {
        return Err(UshmaError::InvalidTemperature { kelvin: t1 });
    }
    if p1 <= 0.0 {
        return Err(UshmaError::InvalidPressure { pascals: p1 });
    }
    if pressure_ratio <= 1.0 {
        return Err(UshmaError::InvalidPressureRatio {
            ratio: pressure_ratio,
        });
    }
    if gamma <= 1.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("gamma {gamma} must be > 1"),
        });
    }
    if moles <= 0.0 {
        return Err(UshmaError::InvalidParameter {
            reason: format!("moles {moles} must be positive"),
        });
    }

    let rp = pressure_ratio;
    let cv = GAS_CONSTANT / (gamma - 1.0);
    let cp = gamma * cv;
    let n = moles;

    // State 1
    let v1 = n * GAS_CONSTANT * t1 / p1;

    // State 2: isentropic compression
    let p2 = p1 * rp;
    let t2 = t1 * rp.powf((gamma - 1.0) / gamma);
    let v2 = n * GAS_CONSTANT * t2 / p2;

    if t3 <= t2 {
        return Err(UshmaError::InvalidCycleParameter {
            reason: format!("turbine inlet T3={t3} K must exceed compressor outlet T2={t2:.1} K"),
        });
    }

    // State 3: isobaric heat addition
    let p3 = p2;
    let v3 = n * GAS_CONSTANT * t3 / p3;

    // State 4: isentropic expansion
    let p4 = p1;
    let t4 = t3 * (1.0 / rp).powf((gamma - 1.0) / gamma);
    let v4 = n * GAS_CONSTANT * t4 / p4;

    let q_in = n * cp * (t3 - t2);
    let q_out = n * cp * (t4 - t1);
    let w_net = q_in - q_out;
    let efficiency = w_net / q_in;

    let w_compressor = n * cp * (t2 - t1);
    let w_turbine = n * cp * (t3 - t4);
    let back_work_ratio = w_compressor / w_turbine;

    let s1 = 0.0;
    let s2 = s1;
    let s3 = s2 + n * cp * (t3 / t2).ln();
    let s4 = s3;

    Ok(CycleResult {
        kind: CycleKind::Brayton,
        state_points: vec![
            StatePoint {
                temperature: t1,
                pressure: p1,
                volume: v1,
                entropy: s1,
                enthalpy: n * cp * t1,
            },
            StatePoint {
                temperature: t2,
                pressure: p2,
                volume: v2,
                entropy: s2,
                enthalpy: n * cp * t2,
            },
            StatePoint {
                temperature: t3,
                pressure: p3,
                volume: v3,
                entropy: s3,
                enthalpy: n * cp * t3,
            },
            StatePoint {
                temperature: t4,
                pressure: p4,
                volume: v4,
                entropy: s4,
                enthalpy: n * cp * t4,
            },
        ],
        processes: vec![
            ProcessKind::Isentropic,
            ProcessKind::Isobaric,
            ProcessKind::Isentropic,
            ProcessKind::Isobaric,
        ],
        heat_in: q_in,
        heat_out: q_out,
        net_work: w_net,
        efficiency,
        back_work_ratio,
        gamma,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_state_point_construction() {
        let sp = StatePoint {
            temperature: 300.0,
            pressure: 101_325.0,
            volume: 0.02241,
            entropy: 150.0,
            enthalpy: 6000.0,
        };
        assert!((sp.temperature - 300.0).abs() < 1e-10);
    }

    #[test]
    fn test_process_kind_traits() {
        let p = ProcessKind::Isentropic;
        assert_eq!(p, ProcessKind::Isentropic);
        assert_ne!(p, ProcessKind::Isochoric);
        let debug = format!("{p:?}");
        assert!(debug.contains("Isentropic"));
    }

    #[test]
    fn test_cycle_kind_traits() {
        assert_eq!(CycleKind::Otto, CycleKind::Otto);
        assert_ne!(CycleKind::Otto, CycleKind::Diesel);
    }

    #[test]
    fn test_state_point_serde_roundtrip() {
        let sp = StatePoint {
            temperature: 500.0,
            pressure: 200_000.0,
            volume: 0.05,
            entropy: 200.0,
            enthalpy: 10_000.0,
        };
        let json = serde_json::to_string(&sp).unwrap();
        let back: StatePoint = serde_json::from_str(&json).unwrap();
        assert!((back.temperature - 500.0).abs() < 1e-10);
    }

    #[test]
    fn test_cycle_kind_serde_roundtrip() {
        let json = serde_json::to_string(&CycleKind::Brayton).unwrap();
        let back: CycleKind = serde_json::from_str(&json).unwrap();
        assert_eq!(back, CycleKind::Brayton);
    }

    #[test]
    fn test_process_kind_serde_roundtrip() {
        let json = serde_json::to_string(&ProcessKind::Isenthalpic).unwrap();
        let back: ProcessKind = serde_json::from_str(&json).unwrap();
        assert_eq!(back, ProcessKind::Isenthalpic);
    }

    // --- Process helper tests ---

    #[test]
    fn test_isochoric_heat() {
        // 1 mol, Cv=20.8 J/(mol·K), 300→400 K → Q = 2080 J
        let q = isochoric_heat(1.0, 20.8, 300.0, 400.0);
        assert!((q - 2080.0).abs() < 1e-10);
    }

    #[test]
    fn test_isobaric_heat() {
        // 1 mol, Cp=29.1 J/(mol·K), 300→400 K → Q = 2910 J
        let q = isobaric_heat(1.0, 29.1, 300.0, 400.0);
        assert!((q - 2910.0).abs() < 1e-10);
    }

    #[test]
    fn test_adiabatic_work_expansion() {
        // Expansion: T drops → positive work
        let w = adiabatic_work(1.0, 20.8, 500.0, 300.0);
        assert!(w > 0.0);
        assert!((w - 20.8 * 200.0).abs() < 1e-10);
    }

    #[test]
    fn test_adiabatic_work_compression() {
        // Compression: T rises → negative work
        let w = adiabatic_work(1.0, 20.8, 300.0, 500.0);
        assert!(w < 0.0);
    }

    #[test]
    fn test_adiabatic_temp_from_pressure_roundtrip() {
        let t2 = adiabatic_temp_from_pressure(300.0, 100_000.0, 1_000_000.0, 1.4).unwrap();
        let t_back = adiabatic_temp_from_pressure(t2, 1_000_000.0, 100_000.0, 1.4).unwrap();
        assert!((t_back - 300.0).abs() < 1e-10);
    }

    #[test]
    fn test_adiabatic_temp_from_pressure_invalid() {
        assert!(adiabatic_temp_from_pressure(0.0, 100_000.0, 200_000.0, 1.4).is_err());
        assert!(adiabatic_temp_from_pressure(300.0, 0.0, 200_000.0, 1.4).is_err());
        assert!(adiabatic_temp_from_pressure(300.0, 100_000.0, 0.0, 1.4).is_err());
        assert!(adiabatic_temp_from_pressure(300.0, 100_000.0, 200_000.0, 1.0).is_err());
    }

    // --- Otto cycle tests ---

    #[test]
    fn test_otto_efficiency_analytical() {
        // η = 1 - 1/r^(γ-1) = 1 - 1/8^0.4 = 1 - 0.4353 = 0.5647
        let result = otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0).unwrap();
        let expected_eta = 1.0 - 1.0 / 8.0_f64.powf(0.4);
        assert!(
            (result.efficiency - expected_eta).abs() < 1e-6,
            "Otto η={}, expected {}",
            result.efficiency,
            expected_eta
        );
    }

    #[test]
    fn test_otto_energy_conservation() {
        let result = otto_cycle(300.0, 101_325.0, 10.0, 40_000.0, 1.4, 1.0).unwrap();
        let balance = (result.heat_in - result.net_work - result.heat_out).abs();
        assert!(balance < 1e-6, "Energy balance error: {balance} J");
    }

    #[test]
    fn test_otto_four_state_points() {
        let result = otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0).unwrap();
        assert_eq!(result.state_points.len(), 4);
        assert_eq!(result.processes.len(), 4);
    }

    #[test]
    fn test_otto_temperature_ordering() {
        let r = otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0).unwrap();
        let sp = &r.state_points;
        // T2 > T1 (compression heats), T3 > T2 (heat added), T4 < T3 (expansion cools)
        assert!(sp[1].temperature > sp[0].temperature);
        assert!(sp[2].temperature > sp[1].temperature);
        assert!(sp[3].temperature < sp[2].temperature);
    }

    #[test]
    fn test_otto_isochoric_volumes() {
        let r = otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0).unwrap();
        let sp = &r.state_points;
        // V2 = V3 (isochoric heat addition)
        assert!((sp[1].volume - sp[2].volume).abs() < 1e-10);
        // V4 = V1 (isochoric heat rejection)
        assert!((sp[3].volume - sp[0].volume).abs() < 1e-10);
    }

    #[test]
    fn test_otto_isentropic_entropy() {
        let r = otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0).unwrap();
        let sp = &r.state_points;
        // S1 = S2 (isentropic compression)
        assert!((sp[0].entropy - sp[1].entropy).abs() < 1e-10);
        // S3 = S4 (isentropic expansion)
        assert!((sp[2].entropy - sp[3].entropy).abs() < 1e-10);
    }

    #[test]
    fn test_otto_higher_compression_higher_efficiency() {
        let r8 = otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0).unwrap();
        let r10 = otto_cycle(300.0, 101_325.0, 10.0, 50_000.0, 1.4, 1.0).unwrap();
        assert!(r10.efficiency > r8.efficiency);
    }

    // --- Diesel cycle tests ---

    #[test]
    fn test_diesel_efficiency_analytical() {
        // η = 1 - (rc^γ - 1) / (γ(rc-1)r^(γ-1))
        let r = 20.0;
        let rc = 2.0;
        let gamma = 1.4;
        let result = diesel_cycle(300.0, 101_325.0, r, rc, gamma, 1.0).unwrap();
        let expected = 1.0 - (rc.powf(gamma) - 1.0) / (gamma * (rc - 1.0) * r.powf(gamma - 1.0));
        assert!(
            (result.efficiency - expected).abs() < 1e-6,
            "Diesel η={}, expected {}",
            result.efficiency,
            expected
        );
    }

    #[test]
    fn test_diesel_energy_conservation() {
        let result = diesel_cycle(300.0, 101_325.0, 18.0, 2.5, 1.4, 1.0).unwrap();
        let balance = (result.heat_in - result.net_work - result.heat_out).abs();
        assert!(balance < 1e-6, "Diesel energy balance error: {balance} J");
    }

    #[test]
    fn test_diesel_isobaric_pressures() {
        let r = diesel_cycle(300.0, 101_325.0, 20.0, 2.0, 1.4, 1.0).unwrap();
        let sp = &r.state_points;
        // P2 = P3 (isobaric heat addition)
        assert!((sp[1].pressure - sp[2].pressure).abs() / sp[1].pressure < 1e-10);
    }

    #[test]
    fn test_diesel_invalid_inputs() {
        assert!(diesel_cycle(300.0, 101_325.0, 1.0, 2.0, 1.4, 1.0).is_err());
        assert!(diesel_cycle(300.0, 101_325.0, 20.0, 1.0, 1.4, 1.0).is_err());
        assert!(diesel_cycle(300.0, 101_325.0, 20.0, 25.0, 1.4, 1.0).is_err()); // rc >= r
    }

    // --- Brayton cycle tests ---

    #[test]
    fn test_brayton_efficiency_analytical() {
        // η = 1 - 1/rp^((γ-1)/γ)
        let rp = 10.0;
        let gamma = 1.4;
        let result = brayton_cycle(300.0, 101_325.0, rp, 1400.0, gamma, 1.0).unwrap();
        let expected = 1.0 - 1.0 / rp.powf((gamma - 1.0) / gamma);
        assert!(
            (result.efficiency - expected).abs() < 1e-6,
            "Brayton η={}, expected {}",
            result.efficiency,
            expected
        );
    }

    #[test]
    fn test_brayton_energy_conservation() {
        let result = brayton_cycle(300.0, 101_325.0, 10.0, 1400.0, 1.4, 1.0).unwrap();
        let balance = (result.heat_in - result.net_work - result.heat_out).abs();
        assert!(balance < 1e-6, "Brayton energy balance error: {balance} J");
    }

    #[test]
    fn test_brayton_back_work_ratio() {
        let result = brayton_cycle(300.0, 101_325.0, 10.0, 1400.0, 1.4, 1.0).unwrap();
        // Gas turbines typically have BWR 40-80%
        assert!(result.back_work_ratio > 0.3 && result.back_work_ratio < 0.9);
    }

    #[test]
    fn test_brayton_t3_must_exceed_t2() {
        // rp=10, γ=1.4 → T2 = 300 * 10^(0.4/1.4) ≈ 579 K
        // T3=500 < T2 → error
        assert!(brayton_cycle(300.0, 101_325.0, 10.0, 500.0, 1.4, 1.0).is_err());
    }

    #[test]
    fn test_brayton_invalid_inputs() {
        assert!(brayton_cycle(300.0, 101_325.0, 1.0, 1400.0, 1.4, 1.0).is_err());
        assert!(brayton_cycle(300.0, 101_325.0, 0.5, 1400.0, 1.4, 1.0).is_err());
    }

    #[test]
    fn test_otto_invalid_inputs() {
        assert!(otto_cycle(0.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0).is_err());
        assert!(otto_cycle(300.0, 0.0, 8.0, 50_000.0, 1.4, 1.0).is_err());
        assert!(otto_cycle(300.0, 101_325.0, 1.0, 50_000.0, 1.4, 1.0).is_err());
        assert!(otto_cycle(300.0, 101_325.0, 0.5, 50_000.0, 1.4, 1.0).is_err());
        assert!(otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.0, 1.0).is_err());
        assert!(otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.4, 0.0).is_err());
        // heat_in must be positive
        assert!(otto_cycle(300.0, 101_325.0, 8.0, 0.0, 1.4, 1.0).is_err());
        assert!(otto_cycle(300.0, 101_325.0, 8.0, -100.0, 1.4, 1.0).is_err());
    }

    #[test]
    fn test_diagram_point() {
        let dp = DiagramPoint { x: 1.0, y: 2.0 };
        assert!((dp.x - 1.0).abs() < 1e-10);
        assert!((dp.y - 2.0).abs() < 1e-10);
    }
}
