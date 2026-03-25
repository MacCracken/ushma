//! Steam tables — saturated and superheated water/steam properties.
//!
//! All SI units: kelvins, pascals, joules, kilograms.
//! Data source: IAPWS-IF97 / NIST Steam Tables.

use serde::{Deserialize, Serialize};

use crate::error::{Result, UshmaError};

/// Properties at a saturated (liquid-vapor equilibrium) state point.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SaturatedEntry {
    /// Temperature (K).
    pub temperature: f64,
    /// Saturation pressure (Pa).
    pub pressure: f64,
    /// Specific volume of saturated liquid v_f (m³/kg).
    pub v_f: f64,
    /// Specific volume of saturated vapor v_g (m³/kg).
    pub v_g: f64,
    /// Specific enthalpy of saturated liquid h_f (J/kg).
    pub h_f: f64,
    /// Specific enthalpy of vaporization h_fg (J/kg).
    pub h_fg: f64,
    /// Specific enthalpy of saturated vapor h_g (J/kg).
    pub h_g: f64,
    /// Specific entropy of saturated liquid s_f (J/(kg·K)).
    pub s_f: f64,
    /// Specific entropy of saturated vapor s_g (J/(kg·K)).
    pub s_g: f64,
}

// --- Saturated steam table (temperature-indexed) ---
// Source: IAPWS-IF97 / Cengel & Boles Thermodynamics, Appendix A-4
// 41 entries from triple point (273.16 K) to critical point (647.096 K)

const SATURATED_TABLE: &[SaturatedEntry] = &[
    // 273.16 K (0.01°C) — triple point
    SaturatedEntry {
        temperature: 273.16,
        pressure: 611.7,
        v_f: 0.001_000_2,
        v_g: 206.1,
        h_f: 0.0,
        h_fg: 2_501_000.0,
        h_g: 2_501_000.0,
        s_f: 0.0,
        s_g: 9156.0,
    },
    // 278.15 K (5°C)
    SaturatedEntry {
        temperature: 278.15,
        pressure: 872.6,
        v_f: 0.001_000_0,
        v_g: 147.1,
        h_f: 21_000.0,
        h_fg: 2_489_000.0,
        h_g: 2_510_000.0,
        s_f: 76.0,
        s_g: 9025.0,
    },
    // 283.15 K (10°C)
    SaturatedEntry {
        temperature: 283.15,
        pressure: 1228.0,
        v_f: 0.001_000_4,
        v_g: 106.4,
        h_f: 42_000.0,
        h_fg: 2_478_000.0,
        h_g: 2_520_000.0,
        s_f: 151.0,
        s_g: 8901.0,
    },
    // 288.15 K (15°C)
    SaturatedEntry {
        temperature: 288.15,
        pressure: 1706.0,
        v_f: 0.001_001_0,
        v_g: 77.93,
        h_f: 63_000.0,
        h_fg: 2_466_000.0,
        h_g: 2_529_000.0,
        s_f: 224.0,
        s_g: 8781.0,
    },
    // 293.15 K (20°C)
    SaturatedEntry {
        temperature: 293.15,
        pressure: 2339.0,
        v_f: 0.001_002_0,
        v_g: 57.79,
        h_f: 84_000.0,
        h_fg: 2_454_000.0,
        h_g: 2_538_000.0,
        s_f: 296.0,
        s_g: 8667.0,
    },
    // 298.15 K (25°C)
    SaturatedEntry {
        temperature: 298.15,
        pressure: 3169.0,
        v_f: 0.001_003_0,
        v_g: 43.36,
        h_f: 105_000.0,
        h_fg: 2_442_000.0,
        h_g: 2_547_000.0,
        s_f: 367.0,
        s_g: 8558.0,
    },
    // 303.15 K (30°C)
    SaturatedEntry {
        temperature: 303.15,
        pressure: 4246.0,
        v_f: 0.001_004_0,
        v_g: 32.90,
        h_f: 126_000.0,
        h_fg: 2_431_000.0,
        h_g: 2_556_000.0,
        s_f: 437.0,
        s_g: 8453.0,
    },
    // 313.15 K (40°C)
    SaturatedEntry {
        temperature: 313.15,
        pressure: 7384.0,
        v_f: 0.001_008_0,
        v_g: 19.52,
        h_f: 168_000.0,
        h_fg: 2_407_000.0,
        h_g: 2_574_000.0,
        s_f: 572.0,
        s_g: 8257.0,
    },
    // 323.15 K (50°C)
    SaturatedEntry {
        temperature: 323.15,
        pressure: 12_349.0,
        v_f: 0.001_012_0,
        v_g: 12.03,
        h_f: 209_000.0,
        h_fg: 2_383_000.0,
        h_g: 2_592_000.0,
        s_f: 704.0,
        s_g: 8076.0,
    },
    // 333.15 K (60°C)
    SaturatedEntry {
        temperature: 333.15,
        pressure: 19_940.0,
        v_f: 0.001_017_0,
        v_g: 7.671,
        h_f: 251_000.0,
        h_fg: 2_359_000.0,
        h_g: 2_609_000.0,
        s_f: 831.0,
        s_g: 7909.0,
    },
    // 343.15 K (70°C)
    SaturatedEntry {
        temperature: 343.15,
        pressure: 31_190.0,
        v_f: 0.001_023_0,
        v_g: 5.042,
        h_f: 293_000.0,
        h_fg: 2_334_000.0,
        h_g: 2_627_000.0,
        s_f: 955.0,
        s_g: 7755.0,
    },
    // 353.15 K (80°C)
    SaturatedEntry {
        temperature: 353.15,
        pressure: 47_390.0,
        v_f: 0.001_029_0,
        v_g: 3.407,
        h_f: 335_000.0,
        h_fg: 2_309_000.0,
        h_g: 2_644_000.0,
        s_f: 1075.0,
        s_g: 7612.0,
    },
    // 363.15 K (90°C)
    SaturatedEntry {
        temperature: 363.15,
        pressure: 70_140.0,
        v_f: 0.001_036_0,
        v_g: 2.361,
        h_f: 377_000.0,
        h_fg: 2_283_000.0,
        h_g: 2_660_000.0,
        s_f: 1193.0,
        s_g: 7479.0,
    },
    // 373.15 K (100°C) — normal boiling point
    SaturatedEntry {
        temperature: 373.15,
        pressure: 101_350.0,
        v_f: 0.001_044_0,
        v_g: 1.673,
        h_f: 419_000.0,
        h_fg: 2_257_000.0,
        h_g: 2_676_000.0,
        s_f: 1307.0,
        s_g: 7355.0,
    },
    // 383.15 K (110°C)
    SaturatedEntry {
        temperature: 383.15,
        pressure: 143_300.0,
        v_f: 0.001_052_0,
        v_g: 1.210,
        h_f: 461_000.0,
        h_fg: 2_230_000.0,
        h_g: 2_691_000.0,
        s_f: 1419.0,
        s_g: 7239.0,
    },
    // 393.15 K (120°C)
    SaturatedEntry {
        temperature: 393.15,
        pressure: 198_500.0,
        v_f: 0.001_060_0,
        v_g: 0.8919,
        h_f: 504_000.0,
        h_fg: 2_203_000.0,
        h_g: 2_706_000.0,
        s_f: 1528.0,
        s_g: 7130.0,
    },
    // 403.15 K (130°C)
    SaturatedEntry {
        temperature: 403.15,
        pressure: 270_100.0,
        v_f: 0.001_070_0,
        v_g: 0.6685,
        h_f: 546_000.0,
        h_fg: 2_174_000.0,
        h_g: 2_720_000.0,
        s_f: 1635.0,
        s_g: 7027.0,
    },
    // 413.15 K (140°C)
    SaturatedEntry {
        temperature: 413.15,
        pressure: 361_300.0,
        v_f: 0.001_080_0,
        v_g: 0.5089,
        h_f: 589_000.0,
        h_fg: 2_145_000.0,
        h_g: 2_734_000.0,
        s_f: 1739.0,
        s_g: 6930.0,
    },
    // 423.15 K (150°C)
    SaturatedEntry {
        temperature: 423.15,
        pressure: 475_800.0,
        v_f: 0.001_091_0,
        v_g: 0.3928,
        h_f: 632_000.0,
        h_fg: 2_114_000.0,
        h_g: 2_746_000.0,
        s_f: 1842.0,
        s_g: 6838.0,
    },
    // 433.15 K (160°C)
    SaturatedEntry {
        temperature: 433.15,
        pressure: 617_800.0,
        v_f: 0.001_102_0,
        v_g: 0.3071,
        h_f: 675_000.0,
        h_fg: 2_083_000.0,
        h_g: 2_758_000.0,
        s_f: 1943.0,
        s_g: 6750.0,
    },
    // 443.15 K (170°C)
    SaturatedEntry {
        temperature: 443.15,
        pressure: 791_700.0,
        v_f: 0.001_114_0,
        v_g: 0.2428,
        h_f: 719_000.0,
        h_fg: 2_050_000.0,
        h_g: 2_769_000.0,
        s_f: 2042.0,
        s_g: 6666.0,
    },
    // 453.15 K (180°C)
    SaturatedEntry {
        temperature: 453.15,
        pressure: 1_002_000.0,
        v_f: 0.001_127_0,
        v_g: 0.1941,
        h_f: 763_000.0,
        h_fg: 2_015_000.0,
        h_g: 2_778_000.0,
        s_f: 2139.0,
        s_g: 6586.0,
    },
    // 463.15 K (190°C)
    SaturatedEntry {
        temperature: 463.15,
        pressure: 1_254_000.0,
        v_f: 0.001_141_0,
        v_g: 0.1565,
        h_f: 807_000.0,
        h_fg: 1_979_000.0,
        h_g: 2_786_000.0,
        s_f: 2236.0,
        s_g: 6507.0,
    },
    // 473.15 K (200°C)
    SaturatedEntry {
        temperature: 473.15,
        pressure: 1_554_000.0,
        v_f: 0.001_157_0,
        v_g: 0.1274,
        h_f: 852_000.0,
        h_fg: 1_941_000.0,
        h_g: 2_793_000.0,
        s_f: 2331.0,
        s_g: 6431.0,
    },
    // 483.15 K (210°C)
    SaturatedEntry {
        temperature: 483.15,
        pressure: 1_907_000.0,
        v_f: 0.001_173_0,
        v_g: 0.1044,
        h_f: 897_000.0,
        h_fg: 1_900_000.0,
        h_g: 2_798_000.0,
        s_f: 2425.0,
        s_g: 6357.0,
    },
    // 493.15 K (220°C)
    SaturatedEntry {
        temperature: 493.15,
        pressure: 2_320_000.0,
        v_f: 0.001_190_0,
        v_g: 0.08620,
        h_f: 943_000.0,
        h_fg: 1_858_000.0,
        h_g: 2_801_000.0,
        s_f: 2518.0,
        s_g: 6284.0,
    },
    // 503.15 K (230°C)
    SaturatedEntry {
        temperature: 503.15,
        pressure: 2_795_000.0,
        v_f: 0.001_209_0,
        v_g: 0.07159,
        h_f: 990_000.0,
        h_fg: 1_813_000.0,
        h_g: 2_803_000.0,
        s_f: 2610.0,
        s_g: 6212.0,
    },
    // 513.15 K (240°C)
    SaturatedEntry {
        temperature: 513.15,
        pressure: 3_344_000.0,
        v_f: 0.001_229_0,
        v_g: 0.05977,
        h_f: 1_037_000.0,
        h_fg: 1_766_000.0,
        h_g: 2_803_000.0,
        s_f: 2702.0,
        s_g: 6141.0,
    },
    // 523.15 K (250°C)
    SaturatedEntry {
        temperature: 523.15,
        pressure: 3_973_000.0,
        v_f: 0.001_251_0,
        v_g: 0.05013,
        h_f: 1_086_000.0,
        h_fg: 1_716_000.0,
        h_g: 2_801_000.0,
        s_f: 2794.0,
        s_g: 6070.0,
    },
    // 533.15 K (260°C)
    SaturatedEntry {
        temperature: 533.15,
        pressure: 4_688_000.0,
        v_f: 0.001_276_0,
        v_g: 0.04221,
        h_f: 1_135_000.0,
        h_fg: 1_662_000.0,
        h_g: 2_797_000.0,
        s_f: 2885.0,
        s_g: 5998.0,
    },
    // 543.15 K (270°C)
    SaturatedEntry {
        temperature: 543.15,
        pressure: 5_499_000.0,
        v_f: 0.001_302_0,
        v_g: 0.03564,
        h_f: 1_185_000.0,
        h_fg: 1_606_000.0,
        h_g: 2_790_000.0,
        s_f: 2976.0,
        s_g: 5925.0,
    },
    // 553.15 K (280°C)
    SaturatedEntry {
        temperature: 553.15,
        pressure: 6_412_000.0,
        v_f: 0.001_332_0,
        v_g: 0.03017,
        h_f: 1_236_000.0,
        h_fg: 1_545_000.0,
        h_g: 2_780_000.0,
        s_f: 3068.0,
        s_g: 5850.0,
    },
    // 563.15 K (290°C)
    SaturatedEntry {
        temperature: 563.15,
        pressure: 7_436_000.0,
        v_f: 0.001_366_0,
        v_g: 0.02557,
        h_f: 1_290_000.0,
        h_fg: 1_478_000.0,
        h_g: 2_767_000.0,
        s_f: 3160.0,
        s_g: 5773.0,
    },
    // 573.15 K (300°C)
    SaturatedEntry {
        temperature: 573.15,
        pressure: 8_581_000.0,
        v_f: 0.001_404_0,
        v_g: 0.02168,
        h_f: 1_345_000.0,
        h_fg: 1_405_000.0,
        h_g: 2_749_000.0,
        s_f: 3254.0,
        s_g: 5693.0,
    },
    // 583.15 K (310°C)
    SaturatedEntry {
        temperature: 583.15,
        pressure: 9_856_000.0,
        v_f: 0.001_448_0,
        v_g: 0.01835,
        h_f: 1_402_000.0,
        h_fg: 1_326_000.0,
        h_g: 2_728_000.0,
        s_f: 3350.0,
        s_g: 5608.0,
    },
    // 593.15 K (320°C)
    SaturatedEntry {
        temperature: 593.15,
        pressure: 11_270_000.0,
        v_f: 0.001_499_0,
        v_g: 0.01549,
        h_f: 1_462_000.0,
        h_fg: 1_239_000.0,
        h_g: 2_700_000.0,
        s_f: 3449.0,
        s_g: 5518.0,
    },
    // 603.15 K (330°C)
    SaturatedEntry {
        temperature: 603.15,
        pressure: 12_845_000.0,
        v_f: 0.001_561_0,
        v_g: 0.01300,
        h_f: 1_526_000.0,
        h_fg: 1_141_000.0,
        h_g: 2_666_000.0,
        s_f: 3551.0,
        s_g: 5420.0,
    },
    // 613.15 K (340°C)
    SaturatedEntry {
        temperature: 613.15,
        pressure: 14_586_000.0,
        v_f: 0.001_638_0,
        v_g: 0.01080,
        h_f: 1_594_000.0,
        h_fg: 1_028_000.0,
        h_g: 2_622_000.0,
        s_f: 3658.0,
        s_g: 5311.0,
    },
    // 623.15 K (350°C)
    SaturatedEntry {
        temperature: 623.15,
        pressure: 16_513_000.0,
        v_f: 0.001_740_0,
        v_g: 0.008_813,
        h_f: 1_671_000.0,
        h_fg: 894_000.0,
        h_g: 2_564_000.0,
        s_f: 3777.0,
        s_g: 5185.0,
    },
    // 633.15 K (360°C)
    SaturatedEntry {
        temperature: 633.15,
        pressure: 18_651_000.0,
        v_f: 0.001_893_0,
        v_g: 0.006_945,
        h_f: 1_761_000.0,
        h_fg: 720_000.0,
        h_g: 2_481_000.0,
        s_f: 3915.0,
        s_g: 5030.0,
    },
    // 643.15 K (370°C)
    SaturatedEntry {
        temperature: 643.15,
        pressure: 21_030_000.0,
        v_f: 0.002_213_0,
        v_g: 0.004_926,
        h_f: 1_891_000.0,
        h_fg: 442_000.0,
        h_g: 2_333_000.0,
        s_f: 4111.0,
        s_g: 4798.0,
    },
    // 647.096 K (373.946°C) — critical point
    SaturatedEntry {
        temperature: 647.096,
        pressure: 22_064_000.0,
        v_f: 0.003_155_0,
        v_g: 0.003_155,
        h_f: 2_084_000.0,
        h_fg: 0.0,
        h_g: 2_084_000.0,
        s_f: 4410.0,
        s_g: 4410.0,
    },
];

/// Look up saturated steam properties by temperature (K).
///
/// Uses binary search and linear interpolation between table entries.
/// Valid range: 273.16 K (triple point) to 647.096 K (critical point).
pub fn saturated_by_temperature(temperature: f64) -> Result<SaturatedEntry> {
    let table = SATURATED_TABLE;
    let t_min = table[0].temperature;
    let t_max = table[table.len() - 1].temperature;

    if temperature < t_min || temperature > t_max {
        return Err(UshmaError::SteamTableOutOfRange {
            temperature,
            pressure: 0.0,
        });
    }

    interpolate_by(table, temperature, |e| e.temperature)
}

/// Look up saturated steam properties by pressure (Pa).
///
/// Uses binary search and linear interpolation between table entries.
/// Valid range: 611.7 Pa (triple point) to 22,064,000 Pa (critical point).
pub fn saturated_by_pressure(pressure: f64) -> Result<SaturatedEntry> {
    let table = SATURATED_TABLE;
    let p_min = table[0].pressure;
    let p_max = table[table.len() - 1].pressure;

    if pressure < p_min || pressure > p_max {
        return Err(UshmaError::SteamTableOutOfRange {
            temperature: 0.0,
            pressure,
        });
    }

    interpolate_by(table, pressure, |e| e.pressure)
}

/// Binary search + linear interpolation on a sorted table.
fn interpolate_by(
    table: &[SaturatedEntry],
    target: f64,
    key: fn(&SaturatedEntry) -> f64,
) -> Result<SaturatedEntry> {
    // Find the right interval via binary search
    let mut lo = 0;
    let mut hi = table.len() - 1;

    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if key(&table[mid]) <= target {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    let a = &table[lo];
    let b = &table[hi];
    let ka = key(a);
    let kb = key(b);

    // Exact match on boundary
    if (kb - ka).abs() < 1e-30 {
        return Ok(*a);
    }

    let frac = (target - ka) / (kb - ka);

    let l = hisab::calc::lerp;
    Ok(SaturatedEntry {
        temperature: l(a.temperature, b.temperature, frac),
        pressure: l(a.pressure, b.pressure, frac),
        v_f: l(a.v_f, b.v_f, frac),
        v_g: l(a.v_g, b.v_g, frac),
        h_f: l(a.h_f, b.h_f, frac),
        h_fg: l(a.h_fg, b.h_fg, frac),
        h_g: l(a.h_g, b.h_g, frac),
        s_f: l(a.s_f, b.s_f, frac),
        s_g: l(a.s_g, b.s_g, frac),
    })
}

/// Properties of wet steam (two-phase liquid-vapor mixture).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct WetSteamProperties {
    /// Quality (dryness fraction), 0 = saturated liquid, 1 = saturated vapor.
    pub quality: f64,
    /// Specific volume (m³/kg).
    pub specific_volume: f64,
    /// Specific enthalpy (J/kg).
    pub specific_enthalpy: f64,
    /// Specific entropy (J/(kg·K)).
    pub specific_entropy: f64,
}

/// Compute quality (dryness fraction) from specific volume.
///
/// x = (v - v_f) / (v_g - v_f)
pub fn quality_from_volume(v: f64, entry: &SaturatedEntry) -> Result<f64> {
    let dv = entry.v_g - entry.v_f;
    if dv.abs() < 1e-30 {
        return Err(UshmaError::DivisionByZero {
            context: "v_g equals v_f at critical point".into(),
        });
    }
    let x = (v - entry.v_f) / dv;
    validate_quality(x)
}

/// Compute quality (dryness fraction) from specific enthalpy.
///
/// x = (h - h_f) / h_fg
pub fn quality_from_enthalpy(h: f64, entry: &SaturatedEntry) -> Result<f64> {
    if entry.h_fg.abs() < 1e-30 {
        return Err(UshmaError::DivisionByZero {
            context: "h_fg is zero at critical point".into(),
        });
    }
    let x = (h - entry.h_f) / entry.h_fg;
    validate_quality(x)
}

/// Compute quality (dryness fraction) from specific entropy.
///
/// x = (s - s_f) / (s_g - s_f)
pub fn quality_from_entropy(s: f64, entry: &SaturatedEntry) -> Result<f64> {
    let ds = entry.s_g - entry.s_f;
    if ds.abs() < 1e-30 {
        return Err(UshmaError::DivisionByZero {
            context: "s_g equals s_f at critical point".into(),
        });
    }
    let x = (s - entry.s_f) / ds;
    validate_quality(x)
}

/// Compute wet steam mixture properties from quality and saturated state.
///
/// prop = prop_f + x · (prop_g - prop_f)
pub fn wet_steam_properties(quality: f64, entry: &SaturatedEntry) -> Result<WetSteamProperties> {
    if !(0.0..=1.0).contains(&quality) {
        return Err(UshmaError::InvalidQuality { quality });
    }
    Ok(WetSteamProperties {
        quality,
        specific_volume: entry.v_f + quality * (entry.v_g - entry.v_f),
        specific_enthalpy: entry.h_f + quality * entry.h_fg,
        specific_entropy: entry.s_f + quality * (entry.s_g - entry.s_f),
    })
}

fn validate_quality(x: f64) -> Result<f64> {
    if !(0.0..=1.0).contains(&x) {
        return Err(UshmaError::InvalidQuality { quality: x });
    }
    Ok(x)
}

/// Superheated steam properties at a single (T, P) point.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SuperheatedEntry {
    /// Temperature (K).
    pub temperature: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Specific volume (m³/kg).
    pub specific_volume: f64,
    /// Specific enthalpy (J/kg).
    pub specific_enthalpy: f64,
    /// Specific entropy (J/(kg·K)).
    pub specific_entropy: f64,
}

// --- Superheated steam table ---
// Source: IAPWS-IF97 / Cengel & Boles Thermodynamics, Appendix A-6
// Organized by pressure tiers, each with temperature entries above saturation.
// Pressures: 10 kPa, 50 kPa, 100 kPa, 200 kPa, 500 kPa, 1 MPa, 2 MPa, 5 MPa, 10 MPa

/// Pressure tiers used in the superheated table (Pa).
const SUPERHEAT_PRESSURES: &[f64] = &[
    10_000.0,
    50_000.0,
    100_000.0,
    200_000.0,
    500_000.0,
    1_000_000.0,
    2_000_000.0,
    5_000_000.0,
    10_000_000.0,
];

/// Number of temperature entries per pressure tier.
const TEMPS_PER_TIER: usize = 8;

// Flat array: SUPERHEAT_PRESSURES.len() * TEMPS_PER_TIER entries
// Each tier has entries at temperatures above saturation for that pressure.
const SUPERHEATED_DATA: &[SuperheatedEntry] = &[
    // --- 10 kPa (T_sat ≈ 318.96 K / 45.81°C) ---
    SuperheatedEntry {
        temperature: 323.15,
        pressure: 10_000.0,
        specific_volume: 14.87,
        specific_enthalpy: 2_585_000.0,
        specific_entropy: 8175.0,
    },
    SuperheatedEntry {
        temperature: 373.15,
        pressure: 10_000.0,
        specific_volume: 17.20,
        specific_enthalpy: 2_688_000.0,
        specific_entropy: 8449.0,
    },
    SuperheatedEntry {
        temperature: 423.15,
        pressure: 10_000.0,
        specific_volume: 19.51,
        specific_enthalpy: 2_783_000.0,
        specific_entropy: 8689.0,
    },
    SuperheatedEntry {
        temperature: 473.15,
        pressure: 10_000.0,
        specific_volume: 21.83,
        specific_enthalpy: 2_880_000.0,
        specific_entropy: 8904.0,
    },
    SuperheatedEntry {
        temperature: 523.15,
        pressure: 10_000.0,
        specific_volume: 24.14,
        specific_enthalpy: 2_978_000.0,
        specific_entropy: 9101.0,
    },
    SuperheatedEntry {
        temperature: 573.15,
        pressure: 10_000.0,
        specific_volume: 26.45,
        specific_enthalpy: 3_077_000.0,
        specific_entropy: 9283.0,
    },
    SuperheatedEntry {
        temperature: 673.15,
        pressure: 10_000.0,
        specific_volume: 31.06,
        specific_enthalpy: 3_280_000.0,
        specific_entropy: 9615.0,
    },
    SuperheatedEntry {
        temperature: 773.15,
        pressure: 10_000.0,
        specific_volume: 35.68,
        specific_enthalpy: 3_489_000.0,
        specific_entropy: 9916.0,
    },
    // --- 50 kPa (T_sat ≈ 354.48 K / 81.33°C) ---
    SuperheatedEntry {
        temperature: 373.15,
        pressure: 50_000.0,
        specific_volume: 3.418,
        specific_enthalpy: 2_682_000.0,
        specific_entropy: 7695.0,
    },
    SuperheatedEntry {
        temperature: 423.15,
        pressure: 50_000.0,
        specific_volume: 3.889,
        specific_enthalpy: 2_780_000.0,
        specific_entropy: 7940.0,
    },
    SuperheatedEntry {
        temperature: 473.15,
        pressure: 50_000.0,
        specific_volume: 4.356,
        specific_enthalpy: 2_877_000.0,
        specific_entropy: 8158.0,
    },
    SuperheatedEntry {
        temperature: 523.15,
        pressure: 50_000.0,
        specific_volume: 4.821,
        specific_enthalpy: 2_977_000.0,
        specific_entropy: 8356.0,
    },
    SuperheatedEntry {
        temperature: 573.15,
        pressure: 50_000.0,
        specific_volume: 5.284,
        specific_enthalpy: 3_076_000.0,
        specific_entropy: 8538.0,
    },
    SuperheatedEntry {
        temperature: 623.15,
        pressure: 50_000.0,
        specific_volume: 5.747,
        specific_enthalpy: 3_178_000.0,
        specific_entropy: 8707.0,
    },
    SuperheatedEntry {
        temperature: 673.15,
        pressure: 50_000.0,
        specific_volume: 6.209,
        specific_enthalpy: 3_280_000.0,
        specific_entropy: 8866.0,
    },
    SuperheatedEntry {
        temperature: 773.15,
        pressure: 50_000.0,
        specific_volume: 7.134,
        specific_enthalpy: 3_489_000.0,
        specific_entropy: 9167.0,
    },
    // --- 100 kPa (T_sat ≈ 372.76 K / 99.61°C) ---
    SuperheatedEntry {
        temperature: 373.15,
        pressure: 100_000.0,
        specific_volume: 1.694,
        specific_enthalpy: 2_676_000.0,
        specific_entropy: 7361.0,
    },
    SuperheatedEntry {
        temperature: 423.15,
        pressure: 100_000.0,
        specific_volume: 1.936,
        specific_enthalpy: 2_776_000.0,
        specific_entropy: 7614.0,
    },
    SuperheatedEntry {
        temperature: 473.15,
        pressure: 100_000.0,
        specific_volume: 2.172,
        specific_enthalpy: 2_875_000.0,
        specific_entropy: 7834.0,
    },
    SuperheatedEntry {
        temperature: 523.15,
        pressure: 100_000.0,
        specific_volume: 2.406,
        specific_enthalpy: 2_975_000.0,
        specific_entropy: 8034.0,
    },
    SuperheatedEntry {
        temperature: 573.15,
        pressure: 100_000.0,
        specific_volume: 2.639,
        specific_enthalpy: 3_075_000.0,
        specific_entropy: 8216.0,
    },
    SuperheatedEntry {
        temperature: 623.15,
        pressure: 100_000.0,
        specific_volume: 2.871,
        specific_enthalpy: 3_177_000.0,
        specific_entropy: 8386.0,
    },
    SuperheatedEntry {
        temperature: 673.15,
        pressure: 100_000.0,
        specific_volume: 3.103,
        specific_enthalpy: 3_280_000.0,
        specific_entropy: 8544.0,
    },
    SuperheatedEntry {
        temperature: 773.15,
        pressure: 100_000.0,
        specific_volume: 3.565,
        specific_enthalpy: 3_488_000.0,
        specific_entropy: 8845.0,
    },
    // --- 200 kPa (T_sat ≈ 393.36 K / 120.21°C) ---
    SuperheatedEntry {
        temperature: 423.15,
        pressure: 200_000.0,
        specific_volume: 0.9596,
        specific_enthalpy: 2_769_000.0,
        specific_entropy: 7281.0,
    },
    SuperheatedEntry {
        temperature: 473.15,
        pressure: 200_000.0,
        specific_volume: 1.080,
        specific_enthalpy: 2_871_000.0,
        specific_entropy: 7507.0,
    },
    SuperheatedEntry {
        temperature: 523.15,
        pressure: 200_000.0,
        specific_volume: 1.199,
        specific_enthalpy: 2_971_000.0,
        specific_entropy: 7710.0,
    },
    SuperheatedEntry {
        temperature: 573.15,
        pressure: 200_000.0,
        specific_volume: 1.316,
        specific_enthalpy: 3_072_000.0,
        specific_entropy: 7894.0,
    },
    SuperheatedEntry {
        temperature: 623.15,
        pressure: 200_000.0,
        specific_volume: 1.433,
        specific_enthalpy: 3_175_000.0,
        specific_entropy: 8064.0,
    },
    SuperheatedEntry {
        temperature: 673.15,
        pressure: 200_000.0,
        specific_volume: 1.549,
        specific_enthalpy: 3_278_000.0,
        specific_entropy: 8224.0,
    },
    SuperheatedEntry {
        temperature: 723.15,
        pressure: 200_000.0,
        specific_volume: 1.665,
        specific_enthalpy: 3_383_000.0,
        specific_entropy: 8374.0,
    },
    SuperheatedEntry {
        temperature: 773.15,
        pressure: 200_000.0,
        specific_volume: 1.781,
        specific_enthalpy: 3_488_000.0,
        specific_entropy: 8516.0,
    },
    // --- 500 kPa (T_sat ≈ 424.99 K / 151.84°C) ---
    SuperheatedEntry {
        temperature: 473.15,
        pressure: 500_000.0,
        specific_volume: 0.4249,
        specific_enthalpy: 2_855_000.0,
        specific_entropy: 7059.0,
    },
    SuperheatedEntry {
        temperature: 523.15,
        pressure: 500_000.0,
        specific_volume: 0.4744,
        specific_enthalpy: 2_961_000.0,
        specific_entropy: 7271.0,
    },
    SuperheatedEntry {
        temperature: 573.15,
        pressure: 500_000.0,
        specific_volume: 0.5226,
        specific_enthalpy: 3_064_000.0,
        specific_entropy: 7460.0,
    },
    SuperheatedEntry {
        temperature: 623.15,
        pressure: 500_000.0,
        specific_volume: 0.5701,
        specific_enthalpy: 3_168_000.0,
        specific_entropy: 7633.0,
    },
    SuperheatedEntry {
        temperature: 673.15,
        pressure: 500_000.0,
        specific_volume: 0.6173,
        specific_enthalpy: 3_272_000.0,
        specific_entropy: 7794.0,
    },
    SuperheatedEntry {
        temperature: 723.15,
        pressure: 500_000.0,
        specific_volume: 0.6642,
        specific_enthalpy: 3_378_000.0,
        specific_entropy: 7946.0,
    },
    SuperheatedEntry {
        temperature: 773.15,
        pressure: 500_000.0,
        specific_volume: 0.7109,
        specific_enthalpy: 3_484_000.0,
        specific_entropy: 8088.0,
    },
    SuperheatedEntry {
        temperature: 873.15,
        pressure: 500_000.0,
        specific_volume: 0.8041,
        specific_enthalpy: 3_702_000.0,
        specific_entropy: 8353.0,
    },
    // --- 1 MPa (T_sat ≈ 453.03 K / 179.88°C) ---
    SuperheatedEntry {
        temperature: 473.15,
        pressure: 1_000_000.0,
        specific_volume: 0.2060,
        specific_enthalpy: 2_828_000.0,
        specific_entropy: 6694.0,
    },
    SuperheatedEntry {
        temperature: 523.15,
        pressure: 1_000_000.0,
        specific_volume: 0.2327,
        specific_enthalpy: 2_943_000.0,
        specific_entropy: 6926.0,
    },
    SuperheatedEntry {
        temperature: 573.15,
        pressure: 1_000_000.0,
        specific_volume: 0.2579,
        specific_enthalpy: 3_051_000.0,
        specific_entropy: 7124.0,
    },
    SuperheatedEntry {
        temperature: 623.15,
        pressure: 1_000_000.0,
        specific_volume: 0.2825,
        specific_enthalpy: 3_158_000.0,
        specific_entropy: 7301.0,
    },
    SuperheatedEntry {
        temperature: 673.15,
        pressure: 1_000_000.0,
        specific_volume: 0.3066,
        specific_enthalpy: 3_264_000.0,
        specific_entropy: 7465.0,
    },
    SuperheatedEntry {
        temperature: 723.15,
        pressure: 1_000_000.0,
        specific_volume: 0.3304,
        specific_enthalpy: 3_371_000.0,
        specific_entropy: 7619.0,
    },
    SuperheatedEntry {
        temperature: 773.15,
        pressure: 1_000_000.0,
        specific_volume: 0.3541,
        specific_enthalpy: 3_479_000.0,
        specific_entropy: 7763.0,
    },
    SuperheatedEntry {
        temperature: 873.15,
        pressure: 1_000_000.0,
        specific_volume: 0.4011,
        specific_enthalpy: 3_698_000.0,
        specific_entropy: 8029.0,
    },
    // --- 2 MPa (T_sat ≈ 485.57 K / 212.42°C) ---
    SuperheatedEntry {
        temperature: 523.15,
        pressure: 2_000_000.0,
        specific_volume: 0.1115,
        specific_enthalpy: 2_903_000.0,
        specific_entropy: 6546.0,
    },
    SuperheatedEntry {
        temperature: 573.15,
        pressure: 2_000_000.0,
        specific_volume: 0.1255,
        specific_enthalpy: 3_024_000.0,
        specific_entropy: 6768.0,
    },
    SuperheatedEntry {
        temperature: 623.15,
        pressure: 2_000_000.0,
        specific_volume: 0.1386,
        specific_enthalpy: 3_137_000.0,
        specific_entropy: 6957.0,
    },
    SuperheatedEntry {
        temperature: 673.15,
        pressure: 2_000_000.0,
        specific_volume: 0.1512,
        specific_enthalpy: 3_248_000.0,
        specific_entropy: 7127.0,
    },
    SuperheatedEntry {
        temperature: 723.15,
        pressure: 2_000_000.0,
        specific_volume: 0.1635,
        specific_enthalpy: 3_358_000.0,
        specific_entropy: 7284.0,
    },
    SuperheatedEntry {
        temperature: 773.15,
        pressure: 2_000_000.0,
        specific_volume: 0.1757,
        specific_enthalpy: 3_468_000.0,
        specific_entropy: 7432.0,
    },
    SuperheatedEntry {
        temperature: 873.15,
        pressure: 2_000_000.0,
        specific_volume: 0.1996,
        specific_enthalpy: 3_690_000.0,
        specific_entropy: 7702.0,
    },
    SuperheatedEntry {
        temperature: 973.15,
        pressure: 2_000_000.0,
        specific_volume: 0.2232,
        specific_enthalpy: 3_917_000.0,
        specific_entropy: 7946.0,
    },
    // --- 5 MPa (T_sat ≈ 536.67 K / 263.92°C) ---
    SuperheatedEntry {
        temperature: 573.15,
        pressure: 5_000_000.0,
        specific_volume: 0.04532,
        specific_enthalpy: 2_925_000.0,
        specific_entropy: 6209.0,
    },
    SuperheatedEntry {
        temperature: 623.15,
        pressure: 5_000_000.0,
        specific_volume: 0.05195,
        specific_enthalpy: 3_069_000.0,
        specific_entropy: 6451.0,
    },
    SuperheatedEntry {
        temperature: 673.15,
        pressure: 5_000_000.0,
        specific_volume: 0.05781,
        specific_enthalpy: 3_196_000.0,
        specific_entropy: 6647.0,
    },
    SuperheatedEntry {
        temperature: 723.15,
        pressure: 5_000_000.0,
        specific_volume: 0.06330,
        specific_enthalpy: 3_317_000.0,
        specific_entropy: 6819.0,
    },
    SuperheatedEntry {
        temperature: 773.15,
        pressure: 5_000_000.0,
        specific_volume: 0.06857,
        specific_enthalpy: 3_434_000.0,
        specific_entropy: 6977.0,
    },
    SuperheatedEntry {
        temperature: 873.15,
        pressure: 5_000_000.0,
        specific_volume: 0.07869,
        specific_enthalpy: 3_664_000.0,
        specific_entropy: 7260.0,
    },
    SuperheatedEntry {
        temperature: 973.15,
        pressure: 5_000_000.0,
        specific_volume: 0.08849,
        specific_enthalpy: 3_895_000.0,
        specific_entropy: 7512.0,
    },
    SuperheatedEntry {
        temperature: 1073.15,
        pressure: 5_000_000.0,
        specific_volume: 0.09811,
        specific_enthalpy: 4_132_000.0,
        specific_entropy: 7740.0,
    },
    // --- 10 MPa (T_sat ≈ 584.15 K / 311.0°C) ---
    SuperheatedEntry {
        temperature: 623.15,
        pressure: 10_000_000.0,
        specific_volume: 0.02242,
        specific_enthalpy: 2_924_000.0,
        specific_entropy: 5944.0,
    },
    SuperheatedEntry {
        temperature: 673.15,
        pressure: 10_000_000.0,
        specific_volume: 0.02641,
        specific_enthalpy: 3_097_000.0,
        specific_entropy: 6212.0,
    },
    SuperheatedEntry {
        temperature: 723.15,
        pressure: 10_000_000.0,
        specific_volume: 0.02975,
        specific_enthalpy: 3_241_000.0,
        specific_entropy: 6419.0,
    },
    SuperheatedEntry {
        temperature: 773.15,
        pressure: 10_000_000.0,
        specific_volume: 0.03279,
        specific_enthalpy: 3_374_000.0,
        specific_entropy: 6598.0,
    },
    SuperheatedEntry {
        temperature: 823.15,
        pressure: 10_000_000.0,
        specific_volume: 0.03564,
        specific_enthalpy: 3_500_000.0,
        specific_entropy: 6758.0,
    },
    SuperheatedEntry {
        temperature: 873.15,
        pressure: 10_000_000.0,
        specific_volume: 0.03837,
        specific_enthalpy: 3_624_000.0,
        specific_entropy: 6903.0,
    },
    SuperheatedEntry {
        temperature: 973.15,
        pressure: 10_000_000.0,
        specific_volume: 0.04358,
        specific_enthalpy: 3_868_000.0,
        specific_entropy: 7168.0,
    },
    SuperheatedEntry {
        temperature: 1073.15,
        pressure: 10_000_000.0,
        specific_volume: 0.04859,
        specific_enthalpy: 4_114_000.0,
        specific_entropy: 7403.0,
    },
];

/// Look up superheated steam properties by temperature (K) and pressure (Pa).
///
/// Uses bilinear interpolation across the pressure-temperature grid.
/// Temperature must be above saturation temperature at the given pressure.
pub fn superheated_lookup(temperature: f64, pressure: f64) -> Result<SuperheatedEntry> {
    if temperature <= 0.0 {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    if pressure <= 0.0 {
        return Err(UshmaError::InvalidPressure { pascals: pressure });
    }

    let pressures = SUPERHEAT_PRESSURES;
    let n = TEMPS_PER_TIER;

    // Check pressure bounds
    if pressure < pressures[0] || pressure > pressures[pressures.len() - 1] {
        return Err(UshmaError::SuperheatOutOfRange {
            temperature,
            pressure,
        });
    }

    // Find pressure bracket
    let mut p_lo = 0;
    let mut p_hi = pressures.len() - 1;
    while p_hi - p_lo > 1 {
        let mid = (p_lo + p_hi) / 2;
        if pressures[mid] <= pressure {
            p_lo = mid;
        } else {
            p_hi = mid;
        }
    }

    // Interpolate within each pressure tier, then between tiers
    let entry_lo = interp_tier(&SUPERHEATED_DATA[p_lo * n..(p_lo + 1) * n], temperature)?;
    let entry_hi = interp_tier(&SUPERHEATED_DATA[p_hi * n..(p_hi + 1) * n], temperature)?;

    // Pressure interpolation fraction
    if (pressures[p_hi] - pressures[p_lo]).abs() < 1e-30 {
        return Ok(entry_lo);
    }
    let frac = (pressure - pressures[p_lo]) / (pressures[p_hi] - pressures[p_lo]);

    let l = hisab::calc::lerp;
    Ok(SuperheatedEntry {
        temperature,
        pressure,
        specific_volume: l(entry_lo.specific_volume, entry_hi.specific_volume, frac),
        specific_enthalpy: l(entry_lo.specific_enthalpy, entry_hi.specific_enthalpy, frac),
        specific_entropy: l(entry_lo.specific_entropy, entry_hi.specific_entropy, frac),
    })
}

/// Interpolate within a single pressure tier by temperature.
fn interp_tier(tier: &[SuperheatedEntry], temperature: f64) -> Result<SuperheatedEntry> {
    let t_min = tier[0].temperature;
    let t_max = tier[tier.len() - 1].temperature;

    if temperature < t_min || temperature > t_max {
        return Err(UshmaError::SuperheatOutOfRange {
            temperature,
            pressure: tier[0].pressure,
        });
    }

    let mut lo = 0;
    let mut hi = tier.len() - 1;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if tier[mid].temperature <= temperature {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    let a = &tier[lo];
    let b = &tier[hi];
    let dt = b.temperature - a.temperature;
    if dt.abs() < 1e-30 {
        return Ok(*a);
    }
    let frac = (temperature - a.temperature) / dt;

    let l = hisab::calc::lerp;
    Ok(SuperheatedEntry {
        temperature,
        pressure: a.pressure,
        specific_volume: l(a.specific_volume, b.specific_volume, frac),
        specific_enthalpy: l(a.specific_enthalpy, b.specific_enthalpy, frac),
        specific_entropy: l(a.specific_entropy, b.specific_entropy, frac),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_saturated_by_temp_boiling_point() {
        // 373.15 K (100°C): P_sat ≈ 101,325 Pa
        let e = saturated_by_temperature(373.15).unwrap();
        assert!((e.pressure - 101_325.0).abs() / 101_325.0 < 0.001);
    }

    #[test]
    fn test_saturated_by_temp_triple_point() {
        let e = saturated_by_temperature(273.16).unwrap();
        assert!((e.pressure - 611.7).abs() < 1.0);
        assert!(e.h_fg > 2_400_000.0); // latent heat at triple point
    }

    #[test]
    fn test_saturated_by_temp_critical_point() {
        let e = saturated_by_temperature(647.096).unwrap();
        assert!((e.pressure - 22_064_000.0).abs() / 22_064_000.0 < 0.001);
        assert!(e.h_fg.abs() < 1000.0); // h_fg → 0 at critical point
    }

    #[test]
    fn test_saturated_by_temp_interpolation() {
        // 50°C = 323.15 K is in the table; 55°C = 328.15 K is between entries
        let e = saturated_by_temperature(328.15).unwrap();
        // Should be between 50°C and 60°C values
        let e50 = saturated_by_temperature(323.15).unwrap();
        let e60 = saturated_by_temperature(333.15).unwrap();
        assert!(e.pressure > e50.pressure && e.pressure < e60.pressure);
        assert!(e.h_fg < e50.h_fg && e.h_fg > e60.h_fg); // h_fg decreases with T
    }

    #[test]
    fn test_saturated_by_temp_out_of_range() {
        assert!(saturated_by_temperature(200.0).is_err()); // below triple point
        assert!(saturated_by_temperature(700.0).is_err()); // above critical point
    }

    #[test]
    fn test_saturated_by_pressure_1atm() {
        let e = saturated_by_pressure(101_325.0).unwrap();
        // Should give ~373.15 K
        assert!((e.temperature - 373.15).abs() < 1.0);
    }

    #[test]
    fn test_saturated_by_pressure_triple() {
        let e = saturated_by_pressure(611.7).unwrap();
        assert!((e.temperature - 273.16).abs() < 0.1);
    }

    #[test]
    fn test_saturated_by_pressure_out_of_range() {
        assert!(saturated_by_pressure(100.0).is_err()); // below triple point pressure
        assert!(saturated_by_pressure(25_000_000.0).is_err()); // above critical
    }

    #[test]
    fn test_saturated_enthalpy_consistency() {
        // h_g = h_f + h_fg at every table point
        for entry in SATURATED_TABLE {
            let diff = (entry.h_g - (entry.h_f + entry.h_fg)).abs();
            assert!(
                diff < 2000.0, // within 2 kJ/kg tolerance for rounded values
                "h_g != h_f + h_fg at T={} K: diff={}",
                entry.temperature,
                diff
            );
        }
    }

    #[test]
    fn test_saturated_table_monotonic() {
        // Temperature, pressure must be strictly increasing
        for w in SATURATED_TABLE.windows(2) {
            assert!(
                w[1].temperature > w[0].temperature,
                "T not increasing: {} -> {}",
                w[0].temperature,
                w[1].temperature
            );
            assert!(
                w[1].pressure > w[0].pressure,
                "P not increasing: {} -> {}",
                w[0].pressure,
                w[1].pressure
            );
        }
    }

    #[test]
    fn test_saturated_table_h_fg_decreasing() {
        // h_fg decreases as we approach critical point
        for w in SATURATED_TABLE.windows(2) {
            assert!(
                w[1].h_fg <= w[0].h_fg,
                "h_fg not decreasing: {} -> {} at T={}",
                w[0].h_fg,
                w[1].h_fg,
                w[1].temperature
            );
        }
    }

    #[test]
    fn test_saturated_serde_roundtrip() {
        let e = saturated_by_temperature(373.15).unwrap();
        let json = serde_json::to_string(&e).unwrap();
        let back: SaturatedEntry = serde_json::from_str(&json).unwrap();
        assert!((back.temperature - e.temperature).abs() < 1e-10);
        assert!((back.pressure - e.pressure).abs() < 1e-10);
    }

    // --- Quality tests ---

    #[test]
    fn test_quality_x0_is_liquid() {
        let e = saturated_by_temperature(373.15).unwrap();
        let x = quality_from_volume(e.v_f, &e).unwrap();
        assert!(x.abs() < 1e-10);
    }

    #[test]
    fn test_quality_x1_is_vapor() {
        let e = saturated_by_temperature(373.15).unwrap();
        let x = quality_from_volume(e.v_g, &e).unwrap();
        assert!((x - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_quality_midpoint() {
        let e = saturated_by_temperature(373.15).unwrap();
        let v_mid = 0.5 * (e.v_f + e.v_g);
        let x = quality_from_volume(v_mid, &e).unwrap();
        assert!((x - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_quality_from_enthalpy() {
        let e = saturated_by_temperature(373.15).unwrap();
        // x=0 at h_f
        let x0 = quality_from_enthalpy(e.h_f, &e).unwrap();
        assert!(x0.abs() < 1e-10);
        // x=1 at h_g
        let x1 = quality_from_enthalpy(e.h_g, &e).unwrap();
        assert!((x1 - 1.0).abs() < 0.01); // tolerance for rounded h_g
    }

    #[test]
    fn test_quality_from_entropy() {
        let e = saturated_by_temperature(373.15).unwrap();
        let x0 = quality_from_entropy(e.s_f, &e).unwrap();
        assert!(x0.abs() < 1e-10);
        let x1 = quality_from_entropy(e.s_g, &e).unwrap();
        assert!((x1 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_quality_out_of_range() {
        let e = saturated_by_temperature(373.15).unwrap();
        // Volume below v_f → x < 0
        assert!(quality_from_volume(e.v_f - 0.001, &e).is_err());
        // Volume above v_g → x > 1
        assert!(quality_from_volume(e.v_g + 1.0, &e).is_err());
    }

    #[test]
    fn test_wet_steam_properties_x0() {
        let e = saturated_by_temperature(373.15).unwrap();
        let p = wet_steam_properties(0.0, &e).unwrap();
        assert!((p.specific_volume - e.v_f).abs() < 1e-10);
        assert!((p.specific_enthalpy - e.h_f).abs() < 1e-10);
        assert!((p.specific_entropy - e.s_f).abs() < 1e-10);
    }

    #[test]
    fn test_wet_steam_properties_x1() {
        let e = saturated_by_temperature(373.15).unwrap();
        let p = wet_steam_properties(1.0, &e).unwrap();
        assert!((p.specific_volume - e.v_g).abs() < 1e-10);
        assert!((p.specific_enthalpy - e.h_g).abs() < 1e-10);
        assert!((p.specific_entropy - e.s_g).abs() < 1e-10);
    }

    #[test]
    fn test_wet_steam_properties_invalid_quality() {
        let e = saturated_by_temperature(373.15).unwrap();
        assert!(wet_steam_properties(-0.1, &e).is_err());
        assert!(wet_steam_properties(1.1, &e).is_err());
    }

    #[test]
    fn test_quality_enthalpy_roundtrip() {
        let e = saturated_by_temperature(373.15).unwrap();
        let target_x = 0.75;
        let props = wet_steam_properties(target_x, &e).unwrap();
        let recovered_x = quality_from_enthalpy(props.specific_enthalpy, &e).unwrap();
        assert!((recovered_x - target_x).abs() < 0.01);
    }

    // --- Superheated tests ---

    #[test]
    fn test_superheated_200c_100kpa() {
        // 473.15 K at 100 kPa: v ≈ 2.172 m³/kg, h ≈ 2875 kJ/kg
        let e = superheated_lookup(473.15, 100_000.0).unwrap();
        assert!((e.specific_volume - 2.172).abs() / 2.172 < 0.01);
        assert!((e.specific_enthalpy - 2_875_000.0).abs() / 2_875_000.0 < 0.01);
    }

    #[test]
    fn test_superheated_500c_1mpa() {
        // 773.15 K at 1 MPa
        let e = superheated_lookup(773.15, 1_000_000.0).unwrap();
        assert!((e.specific_volume - 0.3541).abs() / 0.3541 < 0.01);
        assert!((e.specific_enthalpy - 3_479_000.0).abs() / 3_479_000.0 < 0.01);
    }

    #[test]
    fn test_superheated_interpolation_between_pressures() {
        // 573.15 K at 150 kPa (between 100 kPa and 200 kPa tiers)
        let e = superheated_lookup(573.15, 150_000.0).unwrap();
        let e100 = superheated_lookup(573.15, 100_000.0).unwrap();
        let e200 = superheated_lookup(573.15, 200_000.0).unwrap();
        // Enthalpy should be between the two tiers
        assert!(e.specific_enthalpy >= e200.specific_enthalpy.min(e100.specific_enthalpy));
        assert!(e.specific_enthalpy <= e200.specific_enthalpy.max(e100.specific_enthalpy));
    }

    #[test]
    fn test_superheated_out_of_range_pressure() {
        assert!(superheated_lookup(573.15, 5_000.0).is_err()); // below 10 kPa
        assert!(superheated_lookup(573.15, 20_000_000.0).is_err()); // above 10 MPa
    }

    #[test]
    fn test_superheated_out_of_range_temperature() {
        // Temperature below the lowest entry in the tier
        assert!(superheated_lookup(300.0, 1_000_000.0).is_err());
    }

    #[test]
    fn test_superheated_invalid_inputs() {
        assert!(superheated_lookup(0.0, 100_000.0).is_err());
        assert!(superheated_lookup(500.0, 0.0).is_err());
    }

    #[test]
    fn test_superheated_serde_roundtrip() {
        let e = superheated_lookup(573.15, 100_000.0).unwrap();
        let json = serde_json::to_string(&e).unwrap();
        let back: SuperheatedEntry = serde_json::from_str(&json).unwrap();
        assert!((back.specific_enthalpy - e.specific_enthalpy).abs() < 1e-10);
    }

    #[test]
    fn test_superheated_table_data_valid() {
        // Every entry should have positive values
        for entry in SUPERHEATED_DATA {
            assert!(entry.temperature > 0.0);
            assert!(entry.pressure > 0.0);
            assert!(entry.specific_volume > 0.0);
            assert!(entry.specific_enthalpy > 0.0);
            assert!(entry.specific_entropy > 0.0);
        }
    }

    #[test]
    fn test_wet_steam_serde_roundtrip() {
        let e = saturated_by_temperature(373.15).unwrap();
        let p = wet_steam_properties(0.5, &e).unwrap();
        let json = serde_json::to_string(&p).unwrap();
        let back: WetSteamProperties = serde_json::from_str(&json).unwrap();
        assert!((back.quality - 0.5).abs() < 1e-10);
    }
}
