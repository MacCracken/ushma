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

// --- Steam-based cycles ---

/// Rankine cycle analysis (steam power plant).
///
/// Four processes:
/// 1→2: Isentropic compression (pump, liquid)
/// 2→3: Isobaric heat addition (boiler)
/// 3→4: Isentropic expansion (turbine)
/// 4→1: Isobaric heat rejection (condenser)
///
/// - `p_condenser`: condenser pressure (Pa)
/// - `p_boiler`: boiler pressure (Pa)
/// - `t_superheat`: superheated steam temperature (K), or `None` for saturated
#[cfg(feature = "steam")]
pub fn rankine_cycle(
    p_condenser: f64,
    p_boiler: f64,
    t_superheat: Option<f64>,
) -> Result<CycleResult> {
    use crate::steam;

    if p_condenser <= 0.0 {
        return Err(UshmaError::InvalidPressure {
            pascals: p_condenser,
        });
    }
    if p_boiler <= 0.0 {
        return Err(UshmaError::InvalidPressure { pascals: p_boiler });
    }
    if p_boiler <= p_condenser {
        return Err(UshmaError::InvalidCycleParameter {
            reason: format!(
                "boiler pressure {p_boiler} Pa must exceed condenser pressure {p_condenser} Pa"
            ),
        });
    }

    // State 1: saturated liquid at condenser pressure
    let sat_cond = steam::saturated_by_pressure(p_condenser)?;
    let h1 = sat_cond.h_f;
    let s1 = sat_cond.s_f;
    let v1 = sat_cond.v_f;
    let t1 = sat_cond.temperature;

    // State 2: compressed liquid (pump, incompressible approximation)
    let w_pump = v1 * (p_boiler - p_condenser);
    let h2 = h1 + w_pump;
    let s2 = s1; // isentropic
    let t2 = t1; // negligible temperature change in pump
    let v2 = v1; // incompressible

    // State 3: boiler exit
    let (h3, s3, t3, v3) = if let Some(t_sh) = t_superheat {
        let sh = steam::superheated_lookup(t_sh, p_boiler)?;
        (
            sh.specific_enthalpy,
            sh.specific_entropy,
            sh.temperature,
            sh.specific_volume,
        )
    } else {
        let sat_boil = steam::saturated_by_pressure(p_boiler)?;
        (
            sat_boil.h_g,
            sat_boil.s_g,
            sat_boil.temperature,
            sat_boil.v_g,
        )
    };

    // State 4: isentropic expansion to condenser pressure
    let s4 = s3;
    let (h4, t4, v4) = if s4 <= sat_cond.s_g {
        // Wet region
        let x4 = steam::quality_from_entropy(s4, &sat_cond)?;
        let props = steam::wet_steam_properties(x4, &sat_cond)?;
        (
            props.specific_enthalpy,
            sat_cond.temperature,
            props.specific_volume,
        )
    } else {
        // Superheated at condenser pressure (unusual)
        (sat_cond.h_g, sat_cond.temperature, sat_cond.v_g)
    };

    let q_in = h3 - h2;
    let w_turbine = h3 - h4;
    let q_out = h4 - h1;
    let w_net = w_turbine - w_pump;
    let efficiency = w_net / q_in;
    let back_work_ratio = w_pump / w_turbine;

    Ok(CycleResult {
        kind: CycleKind::Rankine,
        state_points: vec![
            StatePoint {
                temperature: t1,
                pressure: p_condenser,
                volume: v1,
                entropy: s1,
                enthalpy: h1,
            },
            StatePoint {
                temperature: t2,
                pressure: p_boiler,
                volume: v2,
                entropy: s2,
                enthalpy: h2,
            },
            StatePoint {
                temperature: t3,
                pressure: p_boiler,
                volume: v3,
                entropy: s3,
                enthalpy: h3,
            },
            StatePoint {
                temperature: t4,
                pressure: p_condenser,
                volume: v4,
                entropy: s4,
                enthalpy: h4,
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
        gamma: 0.0,
    })
}

/// Refrigeration cycle result with COP values.
#[cfg(feature = "steam")]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RefrigerationResult {
    /// Cycle state points and energy data.
    pub cycle: CycleResult,
    /// Coefficient of performance for refrigeration: Q_cold / W.
    pub cop_refrigeration: f64,
    /// Coefficient of performance for heat pump: Q_hot / W = COP_ref + 1.
    pub cop_heat_pump: f64,
}

/// Vapor-compression refrigeration cycle analysis.
///
/// Four processes:
/// 1→2: Isentropic compression (compressor)
/// 2→3: Isobaric heat rejection (condenser → saturated liquid)
/// 3→4: Isenthalpic expansion (throttle valve)
/// 4→1: Isobaric heat absorption (evaporator → saturated vapor)
///
/// - `p_evaporator`: evaporator pressure (Pa)
/// - `p_condenser`: condenser pressure (Pa)
#[cfg(feature = "steam")]
pub fn refrigeration_cycle(p_evaporator: f64, p_condenser: f64) -> Result<RefrigerationResult> {
    use crate::steam;

    if p_evaporator <= 0.0 {
        return Err(UshmaError::InvalidPressure {
            pascals: p_evaporator,
        });
    }
    if p_condenser <= 0.0 {
        return Err(UshmaError::InvalidPressure {
            pascals: p_condenser,
        });
    }
    if p_condenser <= p_evaporator {
        return Err(UshmaError::InvalidCycleParameter {
            reason: format!(
                "condenser pressure {p_condenser} Pa must exceed evaporator pressure {p_evaporator} Pa"
            ),
        });
    }

    let sat_evap = steam::saturated_by_pressure(p_evaporator)?;
    let sat_cond = steam::saturated_by_pressure(p_condenser)?;

    // State 1: saturated vapor at evaporator
    let h1 = sat_evap.h_g;
    let s1 = sat_evap.s_g;
    let t1 = sat_evap.temperature;
    let v1 = sat_evap.v_g;

    // State 2: isentropic compression to condenser pressure
    // Find superheated state where s ≈ s1 at p_condenser via bisection
    let s_target = s1;
    let t_lo_search = sat_cond.temperature + 1.0;
    // Upper bound: try progressively higher temperatures
    let mut t_hi_search = sat_cond.temperature + 200.0;
    // Ensure upper bound is within superheated table range
    if steam::superheated_lookup(t_hi_search, p_condenser).is_err() {
        t_hi_search = sat_cond.temperature + 100.0;
    }

    let mut t_lo_b = t_lo_search;
    let mut t_hi_b = t_hi_search;
    for _ in 0..50 {
        let t_mid = 0.5 * (t_lo_b + t_hi_b);
        if let Ok(sh) = steam::superheated_lookup(t_mid, p_condenser) {
            if sh.specific_entropy < s_target {
                t_lo_b = t_mid;
            } else {
                t_hi_b = t_mid;
            }
        } else {
            t_hi_b = t_mid;
        }
    }
    let t2 = 0.5 * (t_lo_b + t_hi_b);
    let sh2 = steam::superheated_lookup(t2, p_condenser)?;
    let h2 = sh2.specific_enthalpy;
    let s2 = s1;
    let v2 = sh2.specific_volume;

    // State 3: saturated liquid at condenser
    let h3 = sat_cond.h_f;
    let s3 = sat_cond.s_f;
    let t3 = sat_cond.temperature;
    let v3 = sat_cond.v_f;

    // State 4: isenthalpic expansion (throttle), h4 = h3
    let h4 = h3;
    let x4 = steam::quality_from_enthalpy(h4, &sat_evap)?;
    let props4 = steam::wet_steam_properties(x4, &sat_evap)?;
    let t4 = sat_evap.temperature;
    let v4 = props4.specific_volume;
    let s4 = props4.specific_entropy;

    let q_evap = h1 - h4;
    let q_cond = h2 - h3;
    let w_comp = h2 - h1;

    let cop_ref = q_evap / w_comp;
    let cop_hp = q_cond / w_comp;

    Ok(RefrigerationResult {
        cycle: CycleResult {
            kind: CycleKind::Refrigeration,
            state_points: vec![
                StatePoint {
                    temperature: t1,
                    pressure: p_evaporator,
                    volume: v1,
                    entropy: s1,
                    enthalpy: h1,
                },
                StatePoint {
                    temperature: t2,
                    pressure: p_condenser,
                    volume: v2,
                    entropy: s2,
                    enthalpy: h2,
                },
                StatePoint {
                    temperature: t3,
                    pressure: p_condenser,
                    volume: v3,
                    entropy: s3,
                    enthalpy: h3,
                },
                StatePoint {
                    temperature: t4,
                    pressure: p_evaporator,
                    volume: v4,
                    entropy: s4,
                    enthalpy: h4,
                },
            ],
            processes: vec![
                ProcessKind::Isentropic,
                ProcessKind::Isobaric,
                ProcessKind::Isenthalpic,
                ProcessKind::Isobaric,
            ],
            heat_in: q_evap,
            heat_out: q_cond,
            net_work: w_comp,
            efficiency: cop_ref,
            back_work_ratio: 1.0,
            gamma: 0.0,
        },
        cop_refrigeration: cop_ref,
        cop_heat_pump: cop_hp,
    })
}

/// Heat pump COP from refrigeration COP.
///
/// COP_hp = COP_ref + 1 (energy conservation: Q_hot = Q_cold + W).
#[inline]
#[must_use]
pub fn heat_pump_cop(cop_refrigeration: f64) -> f64 {
    cop_refrigeration + 1.0
}

/// Entry for comparing multiple cycles.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CycleComparisonEntry {
    /// Cycle type.
    pub kind: CycleKind,
    /// Thermal efficiency (or COP for refrigeration).
    pub efficiency: f64,
    /// Net work output (J).
    pub net_work: f64,
    /// Back-work ratio.
    pub back_work_ratio: f64,
    /// Carnot efficiency at the same temperature limits.
    pub carnot_efficiency: f64,
    /// Second-law efficiency: η / η_carnot.
    pub second_law_efficiency: f64,
}

/// Compare multiple cycles against their Carnot limits.
///
/// Each entry is `(kind, result, t_high, t_low)` where t_high/t_low
/// define the Carnot reference temperatures (K).
pub fn compare_cycles(
    entries: &[(CycleKind, &CycleResult, f64, f64)],
) -> Vec<CycleComparisonEntry> {
    entries
        .iter()
        .map(|(kind, result, t_high, t_low)| {
            let eta_carnot = if *t_high > *t_low && *t_high > 0.0 {
                1.0 - t_low / t_high
            } else {
                0.0
            };
            let eta_second = if eta_carnot > 0.0 {
                result.efficiency / eta_carnot
            } else {
                0.0
            };
            CycleComparisonEntry {
                kind: *kind,
                efficiency: result.efficiency,
                net_work: result.net_work,
                back_work_ratio: result.back_work_ratio,
                carnot_efficiency: eta_carnot,
                second_law_efficiency: eta_second,
            }
        })
        .collect()
}

// --- Diagram data generation ---

/// Generate T-s diagram points for a cycle.
///
/// Returns interpolated points along each process path.
/// - Isentropic: vertical line (constant s, T varies)
/// - Isochoric: curve from (s1,T1) to (s2,T2) via ds = nCv·ln(T/T_start)
/// - Isobaric: curve from (s1,T1) to (s2,T2) via ds = nCp·ln(T/T_start)
///
/// X-axis = entropy (J/K), Y-axis = temperature (K).
pub fn cycle_ts_diagram(result: &CycleResult, points_per_process: usize) -> Vec<DiagramPoint> {
    let n = result.state_points.len();
    let pts = points_per_process.max(2);
    let mut out = Vec::with_capacity(n * pts);

    for i in 0..n {
        let a = &result.state_points[i];
        let b = &result.state_points[(i + 1) % n];

        for j in 0..pts {
            let frac = j as f64 / (pts - 1) as f64;
            let t = a.temperature + frac * (b.temperature - a.temperature);
            let s = a.entropy + frac * (b.entropy - a.entropy);
            out.push(DiagramPoint { x: s, y: t });
        }
    }

    out
}

/// Generate P-v diagram points for a cycle.
///
/// Returns interpolated points along each process path.
/// - Isentropic: curve via PV^γ = const → P = P_a·(V_a/V)^γ
/// - Isochoric: vertical line (constant V, P varies)
/// - Isobaric: horizontal line (constant P, V varies)
///
/// X-axis = volume (m³), Y-axis = pressure (Pa).
pub fn cycle_pv_diagram(result: &CycleResult, points_per_process: usize) -> Vec<DiagramPoint> {
    let n = result.state_points.len();
    let pts = points_per_process.max(2);
    let gamma = result.gamma;
    let mut out = Vec::with_capacity(n * pts);

    for i in 0..n {
        let a = &result.state_points[i];
        let b = &result.state_points[(i + 1) % n];
        let process = result.processes[i];

        for j in 0..pts {
            let frac = j as f64 / (pts - 1) as f64;

            let (v, p) = match process {
                ProcessKind::Isentropic => {
                    // PV^γ = const → P = P_a * (V_a/V)^γ
                    let v = a.volume + frac * (b.volume - a.volume);
                    let p = a.pressure * (a.volume / v).powf(gamma);
                    (v, p)
                }
                ProcessKind::Isochoric => {
                    // V constant, P varies linearly with T
                    let v = a.volume;
                    let p = a.pressure + frac * (b.pressure - a.pressure);
                    (v, p)
                }
                ProcessKind::Isobaric => {
                    let v = a.volume + frac * (b.volume - a.volume);
                    let p = a.pressure;
                    (v, p)
                }
                _ => {
                    // Linear fallback for other process types
                    let v = a.volume + frac * (b.volume - a.volume);
                    let p = a.pressure + frac * (b.pressure - a.pressure);
                    (v, p)
                }
            };

            out.push(DiagramPoint { x: v, y: p });
        }
    }

    out
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

    // --- Diagram tests ---

    #[test]
    fn test_ts_diagram_point_count() {
        let r = otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0).unwrap();
        let pts = cycle_ts_diagram(&r, 10);
        // 4 processes × 10 points each
        assert_eq!(pts.len(), 40);
    }

    #[test]
    fn test_pv_diagram_point_count() {
        let r = brayton_cycle(300.0, 101_325.0, 10.0, 1400.0, 1.4, 1.0).unwrap();
        let pts = cycle_pv_diagram(&r, 20);
        assert_eq!(pts.len(), 80);
    }

    #[test]
    fn test_ts_diagram_isentropic_vertical() {
        let r = otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0).unwrap();
        let pts = cycle_ts_diagram(&r, 10);
        // Process 0 (1→2) is isentropic: all x (entropy) values should be equal
        let s_ref = pts[0].x;
        for p in &pts[0..10] {
            assert!(
                (p.x - s_ref).abs() < 1e-10,
                "Isentropic segment not vertical: s={}, expected {}",
                p.x,
                s_ref
            );
        }
    }

    #[test]
    fn test_pv_diagram_isochoric_vertical() {
        let r = otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0).unwrap();
        let pts = cycle_pv_diagram(&r, 10);
        // Process 1 (2→3) is isochoric: all x (volume) values should be equal
        let v_ref = pts[10].x;
        for p in &pts[10..20] {
            assert!(
                (p.x - v_ref).abs() < 1e-10,
                "Isochoric segment not vertical: v={}, expected {}",
                p.x,
                v_ref
            );
        }
    }

    #[test]
    fn test_pv_diagram_isobaric_horizontal() {
        let r = brayton_cycle(300.0, 101_325.0, 10.0, 1400.0, 1.4, 1.0).unwrap();
        let pts = cycle_pv_diagram(&r, 10);
        // Process 1 (2→3) is isobaric: all y (pressure) values should be equal
        let p_ref = pts[10].y;
        for p in &pts[10..20] {
            assert!(
                (p.y - p_ref).abs() / p_ref < 1e-10,
                "Isobaric segment not horizontal"
            );
        }
    }

    #[test]
    fn test_pv_diagram_closed_loop() {
        let r = otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0).unwrap();
        let pts = cycle_pv_diagram(&r, 10);
        // Last point of last process should equal first point
        let first = &pts[0];
        let last = &pts[pts.len() - 1];
        assert!(
            (first.x - last.x).abs() / first.x < 1e-10,
            "P-v loop not closed (V)"
        );
        assert!(
            (first.y - last.y).abs() / first.y < 1e-10,
            "P-v loop not closed (P)"
        );
    }

    #[test]
    fn test_diagram_point() {
        let dp = DiagramPoint { x: 1.0, y: 2.0 };
        assert!((dp.x - 1.0).abs() < 1e-10);
        assert!((dp.y - 2.0).abs() < 1e-10);
    }

    // --- Rankine tests ---

    #[cfg(feature = "steam")]
    #[test]
    fn test_rankine_basic() {
        // 10 kPa condenser, 2 MPa boiler, superheated to 573.15 K
        let r = rankine_cycle(10_000.0, 2_000_000.0, Some(573.15)).unwrap();
        assert!(r.efficiency > 0.15 && r.efficiency < 0.40);
        assert!(r.net_work > 0.0);
    }

    #[cfg(feature = "steam")]
    #[test]
    fn test_rankine_energy_conservation() {
        let r = rankine_cycle(10_000.0, 1_000_000.0, Some(573.15)).unwrap();
        let balance = (r.heat_in - r.net_work - r.heat_out).abs();
        assert!(
            balance / r.heat_in < 0.01,
            "Rankine energy balance: {balance}"
        );
    }

    #[cfg(feature = "steam")]
    #[test]
    fn test_rankine_superheated_better_than_saturated() {
        let sat = rankine_cycle(10_000.0, 1_000_000.0, None).unwrap();
        let sup = rankine_cycle(10_000.0, 1_000_000.0, Some(573.15)).unwrap();
        assert!(sup.efficiency > sat.efficiency);
    }

    #[cfg(feature = "steam")]
    #[test]
    fn test_rankine_low_back_work_ratio() {
        let r = rankine_cycle(10_000.0, 2_000_000.0, Some(573.15)).unwrap();
        // Rankine BWR typically 1-3%
        assert!(r.back_work_ratio < 0.05);
    }

    #[cfg(feature = "steam")]
    #[test]
    fn test_rankine_invalid_pressures() {
        assert!(rankine_cycle(0.0, 1_000_000.0, None).is_err());
        assert!(rankine_cycle(10_000.0, 0.0, None).is_err());
        assert!(rankine_cycle(1_000_000.0, 10_000.0, None).is_err()); // boiler < condenser
    }

    // --- Refrigeration tests ---

    #[cfg(feature = "steam")]
    #[test]
    fn test_refrigeration_basic() {
        // Evaporator at ~7 kPa (~40°C), condenser at ~47 kPa (~80°C)
        let r = refrigeration_cycle(7_384.0, 47_390.0).unwrap();
        assert!(r.cop_refrigeration > 0.0);
        assert!(r.cop_heat_pump > r.cop_refrigeration);
    }

    #[cfg(feature = "steam")]
    #[test]
    fn test_refrigeration_cop_hp_identity() {
        let r = refrigeration_cycle(7_384.0, 47_390.0).unwrap();
        // COP_hp = COP_ref + 1 (within floating point tolerance)
        assert!((r.cop_heat_pump - (r.cop_refrigeration + 1.0)).abs() < 0.1);
    }

    #[cfg(feature = "steam")]
    #[test]
    fn test_refrigeration_energy_balance() {
        let r = refrigeration_cycle(7_384.0, 47_390.0).unwrap();
        let c = &r.cycle;
        // Q_cond = Q_evap + W (within tolerance for approximate state 2)
        let balance = (c.heat_out - c.heat_in - c.net_work).abs();
        assert!(
            balance / c.heat_out < 0.05,
            "Refrigeration energy balance: {balance}"
        );
    }

    #[cfg(feature = "steam")]
    #[test]
    fn test_refrigeration_invalid() {
        assert!(refrigeration_cycle(0.0, 47_390.0).is_err());
        assert!(refrigeration_cycle(47_390.0, 7_384.0).is_err()); // condenser < evaporator
    }

    // --- Heat pump + comparison tests ---

    #[test]
    fn test_heat_pump_cop_identity() {
        assert!((heat_pump_cop(3.0) - 4.0).abs() < 1e-10);
        assert!((heat_pump_cop(0.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_compare_cycles_second_law() {
        let otto = otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0).unwrap();
        let brayton = brayton_cycle(300.0, 101_325.0, 10.0, 1400.0, 1.4, 1.0).unwrap();

        let sp_otto = &otto.state_points;
        let t_high_otto = sp_otto[2].temperature;
        let t_low_otto = sp_otto[0].temperature;

        let sp_bray = &brayton.state_points;
        let t_high_bray = sp_bray[2].temperature;
        let t_low_bray = sp_bray[0].temperature;

        let comparison = compare_cycles(&[
            (CycleKind::Otto, &otto, t_high_otto, t_low_otto),
            (CycleKind::Brayton, &brayton, t_high_bray, t_low_bray),
        ]);

        assert_eq!(comparison.len(), 2);
        for entry in &comparison {
            assert!(
                entry.second_law_efficiency <= 1.0,
                "η_II > 1: {}",
                entry.second_law_efficiency
            );
            assert!(entry.second_law_efficiency > 0.0);
            assert!(entry.carnot_efficiency > entry.efficiency);
        }
    }
}
