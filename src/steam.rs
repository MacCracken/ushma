//! Steam tables — saturated and superheated water/steam properties.
//!
//! All SI units: kelvins, pascals, joules, kilograms.
//! Data source: IAPWS-IF97 / NIST Steam Tables.

use serde::{Deserialize, Serialize};

use crate::error::{Result, UshmaError};
use crate::validate::{finite_ge, finite_gt, invalid_param};

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
// 42 entries from triple point (273.16 K) to critical point (647.096 K).
// Generated from IAPWS-IF97 (Region 4 saturation line, Regions 1/2 properties);
// h_fg = h_g − h_f exactly. Critical row uses IAPWS critical constants
// (ρ_c = 322 kg/m³, h_c ≈ 2087.5 kJ/kg, s_c ≈ 4.412 kJ/(kg·K)).

const SATURATED_TABLE: &[SaturatedEntry] = &[
    // 273.16 K (0.01°C) — triple point
    SaturatedEntry {
        temperature: 273.16,
        pressure: 611.657,
        v_f: 0.001_000_2,
        v_g: 206.00,
        h_f: 0.0,
        h_fg: 2_500_900.0,
        h_g: 2_500_900.0,
        s_f: 0.0,
        s_g: 9155.5,
    },
    // 278.15 K (5°C)
    SaturatedEntry {
        temperature: 278.15,
        pressure: 872.575,
        v_f: 0.001_000_1,
        v_g: 147.02,
        h_f: 21_000.0,
        h_fg: 2_489_100.0,
        h_g: 2_510_100.0,
        s_f: 76.3,
        s_g: 9024.9,
    },
    // 283.15 K (10°C)
    SaturatedEntry {
        temperature: 283.15,
        pressure: 1228.18,
        v_f: 0.001_000_3,
        v_g: 106.31,
        h_f: 42_000.0,
        h_fg: 2_477_200.0,
        h_g: 2_519_200.0,
        s_f: 151.1,
        s_g: 8899.8,
    },
    // 288.15 K (15°C)
    SaturatedEntry {
        temperature: 288.15,
        pressure: 1705.74,
        v_f: 0.001_000_9,
        v_g: 77.881,
        h_f: 63_000.0,
        h_fg: 2_465_400.0,
        h_g: 2_528_400.0,
        s_f: 224.5,
        s_g: 8780.4,
    },
    // 293.15 K (20°C)
    SaturatedEntry {
        temperature: 293.15,
        pressure: 2339.21,
        v_f: 0.001_001_8,
        v_g: 57.761,
        h_f: 83_900.0,
        h_fg: 2_453_600.0,
        h_g: 2_537_500.0,
        s_f: 296.5,
        s_g: 8666.1,
    },
    // 298.15 K (25°C)
    SaturatedEntry {
        temperature: 298.15,
        pressure: 3169.75,
        v_f: 0.001_003_0,
        v_g: 43.341,
        h_f: 104_800.0,
        h_fg: 2_441_700.0,
        h_g: 2_546_500.0,
        s_f: 367.3,
        s_g: 8556.8,
    },
    // 303.15 K (30°C)
    SaturatedEntry {
        temperature: 303.15,
        pressure: 4246.69,
        v_f: 0.001_004_4,
        v_g: 32.882,
        h_f: 125_700.0,
        h_fg: 2_429_900.0,
        h_g: 2_555_600.0,
        s_f: 436.8,
        s_g: 8452.1,
    },
    // 313.15 K (40°C)
    SaturatedEntry {
        temperature: 313.15,
        pressure: 7384.43,
        v_f: 0.001_007_9,
        v_g: 19.517,
        h_f: 167_500.0,
        h_fg: 2_406_000.0,
        h_g: 2_573_500.0,
        s_f: 572.4,
        s_g: 8255.7,
    },
    // 323.15 K (50°C)
    SaturatedEntry {
        temperature: 323.15,
        pressure: 12_351.3,
        v_f: 0.001_012_1,
        v_g: 12.028,
        h_f: 209_300.0,
        h_fg: 2_382_000.0,
        h_g: 2_591_300.0,
        s_f: 703.8,
        s_g: 8074.9,
    },
    // 333.15 K (60°C)
    SaturatedEntry {
        temperature: 333.15,
        pressure: 19_945.8,
        v_f: 0.001_017_1,
        v_g: 7.6677,
        h_f: 251_200.0,
        h_fg: 2_357_600.0,
        h_g: 2_608_800.0,
        s_f: 831.2,
        s_g: 7908.2,
    },
    // 343.15 K (70°C)
    SaturatedEntry {
        temperature: 343.15,
        pressure: 31_200.6,
        v_f: 0.001_022_8,
        v_g: 5.0397,
        h_f: 293_000.0,
        h_fg: 2_333_100.0,
        h_g: 2_626_100.0,
        s_f: 955.0,
        s_g: 7754.0,
    },
    // 353.15 K (80°C)
    SaturatedEntry {
        temperature: 353.15,
        pressure: 47_414.7,
        v_f: 0.001_029_0,
        v_g: 3.4053,
        h_f: 334_900.0,
        h_fg: 2_308_100.0,
        h_g: 2_643_000.0,
        s_f: 1075.4,
        s_g: 7611.0,
    },
    // 363.15 K (90°C)
    SaturatedEntry {
        temperature: 363.15,
        pressure: 70_182.4,
        v_f: 0.001_035_9,
        v_g: 2.3591,
        h_f: 377_000.0,
        h_fg: 2_282_500.0,
        h_g: 2_659_500.0,
        s_f: 1192.7,
        s_g: 7478.1,
    },
    // 373.15 K (100°C) — normal boiling point
    SaturatedEntry {
        temperature: 373.15,
        pressure: 101_418.0,
        v_f: 0.001_043_5,
        v_g: 1.6719,
        h_f: 419_100.0,
        h_fg: 2_256_500.0,
        h_g: 2_675_600.0,
        s_f: 1307.0,
        s_g: 7354.1,
    },
    // 383.15 K (110°C)
    SaturatedEntry {
        temperature: 383.15,
        pressure: 143_376.0,
        v_f: 0.001_051_6,
        v_g: 1.2094,
        h_f: 461_400.0,
        h_fg: 2_229_700.0,
        h_g: 2_691_100.0,
        s_f: 1418.7,
        s_g: 7238.0,
    },
    // 393.15 K (120°C)
    SaturatedEntry {
        temperature: 393.15,
        pressure: 198_665.0,
        v_f: 0.001_060_3,
        v_g: 0.891_30,
        h_f: 503_800.0,
        h_fg: 2_202_100.0,
        h_g: 2_705_900.0,
        s_f: 1527.8,
        s_g: 7129.1,
    },
    // 403.15 K (130°C)
    SaturatedEntry {
        temperature: 403.15,
        pressure: 270_260.0,
        v_f: 0.001_069_7,
        v_g: 0.668_08,
        h_f: 546_400.0,
        h_fg: 2_173_700.0,
        h_g: 2_720_100.0,
        s_f: 1634.6,
        s_g: 7026.4,
    },
    // 413.15 K (140°C)
    SaturatedEntry {
        temperature: 413.15,
        pressure: 361_501.0,
        v_f: 0.001_079_8,
        v_g: 0.508_52,
        h_f: 589_200.0,
        h_fg: 2_144_200.0,
        h_g: 2_733_400.0,
        s_f: 1739.3,
        s_g: 6929.3,
    },
    // 423.15 K (150°C)
    SaturatedEntry {
        temperature: 423.15,
        pressure: 476_101.0,
        v_f: 0.001_090_5,
        v_g: 0.392_50,
        h_f: 632_300.0,
        h_fg: 2_113_600.0,
        h_g: 2_745_900.0,
        s_f: 1842.0,
        s_g: 6837.0,
    },
    // 433.15 K (160°C)
    SaturatedEntry {
        temperature: 433.15,
        pressure: 618_139.0,
        v_f: 0.001_102_0,
        v_g: 0.306_82,
        h_f: 675_600.0,
        h_fg: 2_081_800.0,
        h_g: 2_757_400.0,
        s_f: 1942.8,
        s_g: 6749.1,
    },
    // 443.15 K (170°C)
    SaturatedEntry {
        temperature: 443.15,
        pressure: 792_053.0,
        v_f: 0.001_114_3,
        v_g: 0.242_62,
        h_f: 719_200.0,
        h_fg: 2_048_700.0,
        h_g: 2_767_900.0,
        s_f: 2041.9,
        s_g: 6664.9,
    },
    // 453.15 K (180°C)
    SaturatedEntry {
        temperature: 453.15,
        pressure: 1_002_635.0,
        v_f: 0.001_127_4,
        v_g: 0.193_86,
        h_f: 763_200.0,
        h_fg: 2_014_000.0,
        h_g: 2_777_200.0,
        s_f: 2139.5,
        s_g: 6584.1,
    },
    // 463.15 K (190°C)
    SaturatedEntry {
        temperature: 463.15,
        pressure: 1_255_018.0,
        v_f: 0.001_141_4,
        v_g: 0.156_38,
        h_f: 807_600.0,
        h_fg: 1_977_700.0,
        h_g: 2_785_300.0,
        s_f: 2235.8,
        s_g: 6506.0,
    },
    // 473.15 K (200°C)
    SaturatedEntry {
        temperature: 473.15,
        pressure: 1_554_672.0,
        v_f: 0.001_156_5,
        v_g: 0.127_22,
        h_f: 852_400.0,
        h_fg: 1_939_700.0,
        h_g: 2_792_100.0,
        s_f: 2330.8,
        s_g: 6430.3,
    },
    // 483.15 K (210°C)
    SaturatedEntry {
        temperature: 483.15,
        pressure: 1_907_391.0,
        v_f: 0.001_172_7,
        v_g: 0.104_30,
        h_f: 897_700.0,
        h_fg: 1_899_700.0,
        h_g: 2_797_400.0,
        s_f: 2424.8,
        s_g: 6356.5,
    },
    // 493.15 K (220°C)
    SaturatedEntry {
        temperature: 493.15,
        pressure: 2_319_288.0,
        v_f: 0.001_190_2,
        v_g: 0.086_101,
        h_f: 943_600.0,
        h_fg: 1_857_500.0,
        h_g: 2_801_100.0,
        s_f: 2517.8,
        s_g: 6284.2,
    },
    // 503.15 K (230°C)
    SaturatedEntry {
        temperature: 503.15,
        pressure: 2_796_792.0,
        v_f: 0.001_209_0,
        v_g: 0.071_510,
        h_f: 990_200.0,
        h_fg: 1_812_800.0,
        h_g: 2_803_000.0,
        s_f: 2610.2,
        s_g: 6213.1,
    },
    // 513.15 K (240°C)
    SaturatedEntry {
        temperature: 513.15,
        pressure: 3_346_652.0,
        v_f: 0.001_229_5,
        v_g: 0.059_710,
        h_f: 1_037_500.0,
        h_fg: 1_765_600.0,
        h_g: 2_803_100.0,
        s_f: 2701.9,
        s_g: 6142.5,
    },
    // 523.15 K (250°C)
    SaturatedEntry {
        temperature: 523.15,
        pressure: 3_975_939.0,
        v_f: 0.001_251_7,
        v_g: 0.050_087,
        h_f: 1_085_700.0,
        h_fg: 1_715_300.0,
        h_g: 2_801_000.0,
        s_f: 2793.4,
        s_g: 6072.2,
    },
    // 533.15 K (260°C)
    SaturatedEntry {
        temperature: 533.15,
        pressure: 4_692_071.0,
        v_f: 0.001_276_1,
        v_g: 0.042_175,
        h_f: 1_134_800.0,
        h_fg: 1_661_800.0,
        h_g: 2_796_600.0,
        s_f: 2884.7,
        s_g: 6001.7,
    },
    // 543.15 K (270°C)
    SaturatedEntry {
        temperature: 543.15,
        pressure: 5_502_839.0,
        v_f: 0.001_303_0,
        v_g: 0.035_622,
        h_f: 1_185_100.0,
        h_fg: 1_604_600.0,
        h_g: 2_789_700.0,
        s_f: 2976.2,
        s_g: 5930.4,
    },
    // 553.15 K (280°C)
    SaturatedEntry {
        temperature: 553.15,
        pressure: 6_416_459.0,
        v_f: 0.001_332_8,
        v_g: 0.030_154,
        h_f: 1_236_700.0,
        h_fg: 1_543_100.0,
        h_g: 2_779_800.0,
        s_f: 3068.1,
        s_g: 5857.8,
    },
    // 563.15 K (290°C)
    SaturatedEntry {
        temperature: 563.15,
        pressure: 7_441_643.0,
        v_f: 0.001_366_3,
        v_g: 0.025_557,
        h_f: 1_289_800.0,
        h_fg: 1_476_800.0,
        h_g: 2_766_600.0,
        s_f: 3160.8,
        s_g: 5783.2,
    },
    // 573.15 K (300°C)
    SaturatedEntry {
        temperature: 573.15,
        pressure: 8_587_708.0,
        v_f: 0.001_404_2,
        v_g: 0.021_663,
        h_f: 1_344_800.0,
        h_fg: 1_404_800.0,
        h_g: 2_749_600.0,
        s_f: 3254.7,
        s_g: 5705.8,
    },
    // 583.15 K (310°C)
    SaturatedEntry {
        temperature: 583.15,
        pressure: 9_864_746.0,
        v_f: 0.001_447_9,
        v_g: 0.018_339,
        h_f: 1_402_000.0,
        h_fg: 1_325_900.0,
        h_g: 2_727_900.0,
        s_f: 3350.6,
        s_g: 5624.3,
    },
    // 593.15 K (320°C)
    SaturatedEntry {
        temperature: 593.15,
        pressure: 11_283_856.0,
        v_f: 0.001_499_1,
        v_g: 0.015_476,
        h_f: 1_462_100.0,
        h_fg: 1_238_600.0,
        h_g: 2_700_700.0,
        s_f: 3449.1,
        s_g: 5537.3,
    },
    // 603.15 K (330°C)
    SaturatedEntry {
        temperature: 603.15,
        pressure: 12_857_522.0,
        v_f: 0.001_560_6,
        v_g: 0.012_984,
        h_f: 1_525_700.0,
        h_fg: 1_140_500.0,
        h_g: 2_666_200.0,
        s_f: 3551.6,
        s_g: 5442.5,
    },
    // 613.15 K (340°C)
    SaturatedEntry {
        temperature: 613.15,
        pressure: 14_600_181.0,
        v_f: 0.001_637_5,
        v_g: 0.010_784,
        h_f: 1_594_400.0,
        h_fg: 1_027_700.0,
        h_g: 2_622_100.0,
        s_f: 3659.9,
        s_g: 5335.9,
    },
    // 623.15 K (350°C)
    SaturatedEntry {
        temperature: 623.15,
        pressure: 16_529_164.0,
        v_f: 0.001_740_1,
        v_g: 0.008_800_9,
        h_f: 1_670_900.0,
        h_fg: 892_700.0,
        h_g: 2_563_600.0,
        s_f: 3778.3,
        s_g: 5210.9,
    },
    // 633.15 K (360°C)
    SaturatedEntry {
        temperature: 633.15,
        pressure: 18_666_371.0,
        v_f: 0.001_894_5,
        v_g: 0.006_945_0,
        h_f: 1_761_500.0,
        h_fg: 719_500.0,
        h_g: 2_481_000.0,
        s_f: 3916.4,
        s_g: 5052.7,
    },
    // 643.15 K (370°C)
    SaturatedEntry {
        temperature: 643.15,
        pressure: 21_043_367.0,
        v_f: 0.002_222_1,
        v_g: 0.004_946_2,
        h_f: 1_892_600.0,
        h_fg: 440_900.0,
        h_g: 2_333_500.0,
        s_f: 4114.2,
        s_g: 4799.6,
    },
    // 647.096 K (373.946°C) — critical point
    SaturatedEntry {
        temperature: 647.096,
        pressure: 22_064_000.0,
        v_f: 0.003_105_6,
        v_g: 0.003_105_6,
        h_f: 2_087_500.0,
        h_fg: 0.0,
        h_g: 2_087_500.0,
        s_f: 4412.0,
        s_g: 4412.0,
    },
];

// --- IAPWS-IF97 Region 4: the saturation line ---

/// IAPWS-IF97 Region 4 coefficients n₁ … n₁₀.
const R4_N: [f64; 10] = [
    0.116_705_214_527_67e4,
    -0.724_213_165_981_38e6,
    -0.170_738_469_400_92e2,
    0.120_208_247_024_70e5,
    -0.323_255_503_223_33e7,
    0.149_151_086_135_30e2,
    -0.482_326_573_615_91e4,
    0.405_113_405_420_57e6,
    -0.238_555_575_678_49,
    0.650_175_348_447_98e3,
];

/// Lower validity limit of IAPWS-IF97 Region 4 (K).
const R4_T_MIN: f64 = 273.15;
/// Critical temperature of water (K).
const T_CRIT: f64 = 647.096;
/// Critical pressure of water (Pa).
const P_CRIT: f64 = 22_064_000.0;

/// IF97 Region 4 saturation-pressure equation (no range check), in Pa.
#[inline]
fn r4_saturation_pressure(t: f64) -> f64 {
    let n = &R4_N;
    let th = t + n[8] / (t - n[9]);
    let a = th * th + n[0] * th + n[1];
    let b = n[2] * th * th + n[3] * th + n[4];
    let c = n[5] * th * th + n[6] * th + n[7];
    let x = 2.0 * c / (-b + (b * b - 4.0 * a * c).sqrt());
    x.powi(4) * 1e6
}

/// IF97 Region 4 saturation-temperature backward equation (no range check), in K.
#[inline]
fn r4_saturation_temperature(p: f64) -> f64 {
    let n = &R4_N;
    let beta = (p * 1e-6).sqrt().sqrt();
    let e = beta * beta + n[2] * beta + n[5];
    let f = n[0] * beta * beta + n[3] * beta + n[6];
    let g = n[1] * beta * beta + n[4] * beta + n[7];
    let d = 2.0 * g / (-f - (f * f - 4.0 * e * g).sqrt());
    let s = n[9] + d;
    0.5 * (s - (s * s - 4.0 * (n[8] + n[9] * d)).sqrt())
}

/// Saturation pressure of water (Pa) from the IAPWS-IF97 Region 4 equation.
///
/// Valid from 273.15 K to the critical point (647.096 K). Agrees with the
/// IF97 verification values to better than 1e-8 relative.
#[inline]
pub fn saturation_pressure(temperature: f64) -> Result<f64> {
    if !finite_ge(temperature, R4_T_MIN) || temperature > T_CRIT {
        return Err(UshmaError::SteamTableOutOfRange {
            temperature,
            pressure: f64::NAN,
        });
    }
    Ok(r4_saturation_pressure(temperature).min(P_CRIT))
}

/// Saturation temperature of water (K) from the IAPWS-IF97 Region 4 backward equation.
///
/// Valid from 611.213 Pa (saturation at 273.15 K) to the critical pressure (22.064 MPa).
/// This is the exact algebraic inverse of [`saturation_pressure`].
#[inline]
pub fn saturation_temperature(pressure: f64) -> Result<f64> {
    let p_min = r4_saturation_pressure(R4_T_MIN);
    if !finite_ge(pressure, p_min) || pressure > P_CRIT {
        return Err(UshmaError::SteamTableOutOfRange {
            temperature: f64::NAN,
            pressure,
        });
    }
    Ok(r4_saturation_temperature(pressure).clamp(R4_T_MIN, T_CRIT))
}

/// Look up saturated steam properties by temperature (K).
///
/// The saturation pressure comes from the IAPWS-IF97 Region 4 equation; the
/// other properties are interpolated linearly in T between table entries,
/// except v_g, which varies by orders of magnitude: its compressibility-like
/// ratio P·v_g/T is interpolated instead (≈ 4× lower mean error than
/// log-linear interpolation, and no transcendental functions).
/// Valid range: 273.16 K (triple point) to 647.096 K (critical point).
#[inline]
pub fn saturated_by_temperature(temperature: f64) -> Result<SaturatedEntry> {
    let table = SATURATED_TABLE;
    let t_min = table[0].temperature;
    let t_max = table[table.len() - 1].temperature;

    if !finite_ge(temperature, t_min) || temperature > t_max {
        return Err(UshmaError::SteamTableOutOfRange {
            temperature,
            pressure: f64::NAN,
        });
    }

    let pressure = r4_saturation_pressure(temperature).min(P_CRIT);
    Ok(interpolate_saturated(temperature, pressure))
}

/// Look up saturated steam properties by pressure (Pa).
///
/// The saturation temperature comes from the IAPWS-IF97 Region 4 backward
/// equation; the other properties are then interpolated as in
/// [`saturated_by_temperature`].
/// Valid range: 611.657 Pa (triple point) to 22,064,000 Pa (critical point).
#[inline]
pub fn saturated_by_pressure(pressure: f64) -> Result<SaturatedEntry> {
    let table = SATURATED_TABLE;
    let p_min = table[0].pressure;
    let p_max = table[table.len() - 1].pressure;

    if !finite_ge(pressure, p_min) || pressure > p_max {
        return Err(UshmaError::SteamTableOutOfRange {
            temperature: f64::NAN,
            pressure,
        });
    }

    let t = r4_saturation_temperature(pressure).clamp(table[0].temperature, T_CRIT);
    Ok(interpolate_saturated(t, pressure))
}

/// Interpolate the saturated table at temperature `t`, where the saturation
/// pressure is `p` (caller checks range).
#[inline]
fn interpolate_saturated(t: f64, p: f64) -> SaturatedEntry {
    let table = SATURATED_TABLE;
    let hi = table
        .partition_point(|e| e.temperature <= t)
        .clamp(1, table.len() - 1);
    let a = &table[hi - 1];
    let b = &table[hi];
    let frac = (t - a.temperature) / (b.temperature - a.temperature);

    let l = hisab::calc::lerp;
    // P·v_g/T is smooth (≈ R/M·Z), unlike v_g itself.
    let za = a.v_g * a.pressure / a.temperature;
    let zb = b.v_g * b.pressure / b.temperature;
    SaturatedEntry {
        temperature: t,
        pressure: p,
        v_f: l(a.v_f, b.v_f, frac),
        v_g: l(za, zb, frac) * t / p,
        h_f: l(a.h_f, b.h_f, frac),
        h_fg: l(a.h_fg, b.h_fg, frac),
        h_g: l(a.h_g, b.h_g, frac),
        s_f: l(a.s_f, b.s_f, frac),
        s_g: l(a.s_g, b.s_g, frac),
    }
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
#[inline]
pub fn quality_from_volume(v: f64, entry: &SaturatedEntry) -> Result<f64> {
    let dv = entry.v_g - entry.v_f;
    if !finite_ge(dv.abs(), 1e-30) {
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
#[inline]
pub fn quality_from_enthalpy(h: f64, entry: &SaturatedEntry) -> Result<f64> {
    if !finite_ge(entry.h_fg.abs(), 1e-30) {
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
#[inline]
pub fn quality_from_entropy(s: f64, entry: &SaturatedEntry) -> Result<f64> {
    let ds = entry.s_g - entry.s_f;
    if !finite_ge(ds.abs(), 1e-30) {
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
#[inline]
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

/// Tolerance for round-off when a quality is computed from saturated properties
/// (e.g. x(h_g) = 1.000_000_000_000_000_2).
const QUALITY_ROUNDOFF: f64 = 1e-9;

/// Accept x ∈ [0, 1] (with round-off tolerance, clamped); reject NaN and anything else.
#[inline]
fn validate_quality(x: f64) -> Result<f64> {
    if !(-QUALITY_ROUNDOFF..=1.0 + QUALITY_ROUNDOFF).contains(&x) {
        return Err(UshmaError::InvalidQuality { quality: x });
    }
    Ok(x.clamp(0.0, 1.0))
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
// Source: generated from IAPWS-IF97 Region 2 at each (T, P) grid point.
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
        specific_volume: 14.867,
        specific_enthalpy: 2_592_000.0,
        specific_entropy: 8174.1,
    },
    SuperheatedEntry {
        temperature: 373.15,
        pressure: 10_000.0,
        specific_volume: 17.197,
        specific_enthalpy: 2_687_400.0,
        specific_entropy: 8448.8,
    },
    SuperheatedEntry {
        temperature: 423.15,
        pressure: 10_000.0,
        specific_volume: 19.514,
        specific_enthalpy: 2_783_000.0,
        specific_entropy: 8689.2,
    },
    SuperheatedEntry {
        temperature: 473.15,
        pressure: 10_000.0,
        specific_volume: 21.826,
        specific_enthalpy: 2_879_600.0,
        specific_entropy: 8904.8,
    },
    SuperheatedEntry {
        temperature: 523.15,
        pressure: 10_000.0,
        specific_volume: 24.136,
        specific_enthalpy: 2_977_400.0,
        specific_entropy: 9101.4,
    },
    SuperheatedEntry {
        temperature: 573.15,
        pressure: 10_000.0,
        specific_volume: 26.446,
        specific_enthalpy: 3_076_700.0,
        specific_entropy: 9282.7,
    },
    SuperheatedEntry {
        temperature: 673.15,
        pressure: 10_000.0,
        specific_volume: 31.064,
        specific_enthalpy: 3_279_900.0,
        specific_entropy: 9609.3,
    },
    SuperheatedEntry {
        temperature: 773.15,
        pressure: 10_000.0,
        specific_volume: 35.680,
        specific_enthalpy: 3_489_700.0,
        specific_entropy: 9899.7,
    },
    // --- 50 kPa (T_sat ≈ 354.47 K / 81.32°C) ---
    SuperheatedEntry {
        temperature: 373.15,
        pressure: 50_000.0,
        specific_volume: 3.4188,
        specific_enthalpy: 2_682_400.0,
        specific_entropy: 7695.2,
    },
    SuperheatedEntry {
        temperature: 423.15,
        pressure: 50_000.0,
        specific_volume: 3.8899,
        specific_enthalpy: 2_780_200.0,
        specific_entropy: 7941.2,
    },
    SuperheatedEntry {
        temperature: 473.15,
        pressure: 50_000.0,
        specific_volume: 4.3563,
        specific_enthalpy: 2_877_800.0,
        specific_entropy: 8159.1,
    },
    SuperheatedEntry {
        temperature: 523.15,
        pressure: 50_000.0,
        specific_volume: 4.8207,
        specific_enthalpy: 2_976_200.0,
        specific_entropy: 8356.8,
    },
    SuperheatedEntry {
        temperature: 573.15,
        pressure: 50_000.0,
        specific_volume: 5.2841,
        specific_enthalpy: 3_075_800.0,
        specific_entropy: 8538.6,
    },
    SuperheatedEntry {
        temperature: 623.15,
        pressure: 50_000.0,
        specific_volume: 5.7470,
        specific_enthalpy: 3_176_800.0,
        specific_entropy: 8707.6,
    },
    SuperheatedEntry {
        temperature: 673.15,
        pressure: 50_000.0,
        specific_volume: 6.2095,
        specific_enthalpy: 3_279_300.0,
        specific_entropy: 8865.8,
    },
    SuperheatedEntry {
        temperature: 773.15,
        pressure: 50_000.0,
        specific_volume: 7.1339,
        specific_enthalpy: 3_489_200.0,
        specific_entropy: 9156.5,
    },
    // --- 100 kPa (T_sat ≈ 372.76 K / 99.61°C) ---
    SuperheatedEntry {
        temperature: 373.15,
        pressure: 100_000.0,
        specific_volume: 1.6960,
        specific_enthalpy: 2_675_800.0,
        specific_entropy: 7361.0,
    },
    SuperheatedEntry {
        temperature: 423.15,
        pressure: 100_000.0,
        specific_volume: 1.9367,
        specific_enthalpy: 2_776_600.0,
        specific_entropy: 7614.7,
    },
    SuperheatedEntry {
        temperature: 473.15,
        pressure: 100_000.0,
        specific_volume: 2.1725,
        specific_enthalpy: 2_875_500.0,
        specific_entropy: 7835.6,
    },
    SuperheatedEntry {
        temperature: 523.15,
        pressure: 100_000.0,
        specific_volume: 2.4062,
        specific_enthalpy: 2_974_500.0,
        specific_entropy: 8034.6,
    },
    SuperheatedEntry {
        temperature: 573.15,
        pressure: 100_000.0,
        specific_volume: 2.6389,
        specific_enthalpy: 3_074_500.0,
        specific_entropy: 8217.1,
    },
    SuperheatedEntry {
        temperature: 623.15,
        pressure: 100_000.0,
        specific_volume: 2.8710,
        specific_enthalpy: 3_175_800.0,
        specific_entropy: 8386.5,
    },
    SuperheatedEntry {
        temperature: 673.15,
        pressure: 100_000.0,
        specific_volume: 3.1027,
        specific_enthalpy: 3_278_500.0,
        specific_entropy: 8545.1,
    },
    SuperheatedEntry {
        temperature: 773.15,
        pressure: 100_000.0,
        specific_volume: 3.5656,
        specific_enthalpy: 3_488_700.0,
        specific_entropy: 8836.1,
    },
    // --- 200 kPa (T_sat ≈ 393.36 K / 120.21°C) ---
    SuperheatedEntry {
        temperature: 423.15,
        pressure: 200_000.0,
        specific_volume: 0.959_89,
        specific_enthalpy: 2_769_100.0,
        specific_entropy: 7280.9,
    },
    SuperheatedEntry {
        temperature: 473.15,
        pressure: 200_000.0,
        specific_volume: 1.0805,
        specific_enthalpy: 2_870_800.0,
        specific_entropy: 7508.1,
    },
    SuperheatedEntry {
        temperature: 523.15,
        pressure: 200_000.0,
        specific_volume: 1.1989,
        specific_enthalpy: 2_971_300.0,
        specific_entropy: 7710.0,
    },
    SuperheatedEntry {
        temperature: 573.15,
        pressure: 200_000.0,
        specific_volume: 1.3162,
        specific_enthalpy: 3_072_100.0,
        specific_entropy: 7894.0,
    },
    SuperheatedEntry {
        temperature: 623.15,
        pressure: 200_000.0,
        specific_volume: 1.4330,
        specific_enthalpy: 3_173_900.0,
        specific_entropy: 8064.3,
    },
    SuperheatedEntry {
        temperature: 673.15,
        pressure: 200_000.0,
        specific_volume: 1.5493,
        specific_enthalpy: 3_277_000.0,
        specific_entropy: 8223.5,
    },
    SuperheatedEntry {
        temperature: 723.15,
        pressure: 200_000.0,
        specific_volume: 1.6655,
        specific_enthalpy: 3_381_500.0,
        specific_entropy: 8373.3,
    },
    SuperheatedEntry {
        temperature: 773.15,
        pressure: 200_000.0,
        specific_volume: 1.7814,
        specific_enthalpy: 3_487_600.0,
        specific_entropy: 8515.1,
    },
    // --- 500 kPa (T_sat ≈ 424.99 K / 151.84°C) ---
    SuperheatedEntry {
        temperature: 473.15,
        pressure: 500_000.0,
        specific_volume: 0.425_03,
        specific_enthalpy: 2_855_900.0,
        specific_entropy: 7061.1,
    },
    SuperheatedEntry {
        temperature: 523.15,
        pressure: 500_000.0,
        specific_volume: 0.474_43,
        specific_enthalpy: 2_961_100.0,
        specific_entropy: 7272.6,
    },
    SuperheatedEntry {
        temperature: 573.15,
        pressure: 500_000.0,
        specific_volume: 0.522_60,
        specific_enthalpy: 3_064_600.0,
        specific_entropy: 7461.4,
    },
    SuperheatedEntry {
        temperature: 623.15,
        pressure: 500_000.0,
        specific_volume: 0.570_14,
        specific_enthalpy: 3_168_100.0,
        specific_entropy: 7634.5,
    },
    SuperheatedEntry {
        temperature: 673.15,
        pressure: 500_000.0,
        specific_volume: 0.617_29,
        specific_enthalpy: 3_272_300.0,
        specific_entropy: 7795.4,
    },
    SuperheatedEntry {
        temperature: 723.15,
        pressure: 500_000.0,
        specific_volume: 0.664_21,
        specific_enthalpy: 3_377_700.0,
        specific_entropy: 7946.4,
    },
    SuperheatedEntry {
        temperature: 773.15,
        pressure: 500_000.0,
        specific_volume: 0.710_95,
        specific_enthalpy: 3_484_400.0,
        specific_entropy: 8089.1,
    },
    SuperheatedEntry {
        temperature: 873.15,
        pressure: 500_000.0,
        specific_volume: 0.804_10,
        specific_enthalpy: 3_702_500.0,
        specific_entropy: 8354.3,
    },
    // --- 1 MPa (T_sat ≈ 453.04 K / 179.89°C) ---
    SuperheatedEntry {
        temperature: 473.15,
        pressure: 1_000_000.0,
        specific_volume: 0.206_00,
        specific_enthalpy: 2_828_300.0,
        specific_entropy: 6695.5,
    },
    SuperheatedEntry {
        temperature: 523.15,
        pressure: 1_000_000.0,
        specific_volume: 0.232_74,
        specific_enthalpy: 2_943_200.0,
        specific_entropy: 6926.6,
    },
    SuperheatedEntry {
        temperature: 573.15,
        pressure: 1_000_000.0,
        specific_volume: 0.257_98,
        specific_enthalpy: 3_051_700.0,
        specific_entropy: 7124.7,
    },
    SuperheatedEntry {
        temperature: 623.15,
        pressure: 1_000_000.0,
        specific_volume: 0.282_49,
        specific_enthalpy: 3_158_200.0,
        specific_entropy: 7302.8,
    },
    SuperheatedEntry {
        temperature: 673.15,
        pressure: 1_000_000.0,
        specific_volume: 0.306_59,
        specific_enthalpy: 3_264_400.0,
        specific_entropy: 7466.8,
    },
    SuperheatedEntry {
        temperature: 723.15,
        pressure: 1_000_000.0,
        specific_volume: 0.330_44,
        specific_enthalpy: 3_371_200.0,
        specific_entropy: 7619.8,
    },
    SuperheatedEntry {
        temperature: 773.15,
        pressure: 1_000_000.0,
        specific_volume: 0.354_11,
        specific_enthalpy: 3_479_000.0,
        specific_entropy: 7764.0,
    },
    SuperheatedEntry {
        temperature: 873.15,
        pressure: 1_000_000.0,
        specific_volume: 0.401_11,
        specific_enthalpy: 3_698_600.0,
        specific_entropy: 8030.9,
    },
    // --- 2 MPa (T_sat ≈ 485.53 K / 212.38°C) ---
    SuperheatedEntry {
        temperature: 523.15,
        pressure: 2_000_000.0,
        specific_volume: 0.111_48,
        specific_enthalpy: 2_903_200.0,
        specific_entropy: 6547.4,
    },
    SuperheatedEntry {
        temperature: 573.15,
        pressure: 2_000_000.0,
        specific_volume: 0.125_50,
        specific_enthalpy: 3_024_300.0,
        specific_entropy: 6768.5,
    },
    SuperheatedEntry {
        temperature: 623.15,
        pressure: 2_000_000.0,
        specific_volume: 0.138_59,
        specific_enthalpy: 3_137_600.0,
        specific_entropy: 6958.2,
    },
    SuperheatedEntry {
        temperature: 673.15,
        pressure: 2_000_000.0,
        specific_volume: 0.151_21,
        specific_enthalpy: 3_248_200.0,
        specific_entropy: 7129.0,
    },
    SuperheatedEntry {
        temperature: 723.15,
        pressure: 2_000_000.0,
        specific_volume: 0.163_54,
        specific_enthalpy: 3_358_100.0,
        specific_entropy: 7286.3,
    },
    SuperheatedEntry {
        temperature: 773.15,
        pressure: 2_000_000.0,
        specific_volume: 0.175_68,
        specific_enthalpy: 3_468_100.0,
        specific_entropy: 7433.5,
    },
    SuperheatedEntry {
        temperature: 873.15,
        pressure: 2_000_000.0,
        specific_volume: 0.199_61,
        specific_enthalpy: 3_690_700.0,
        specific_entropy: 7704.2,
    },
    SuperheatedEntry {
        temperature: 973.15,
        pressure: 2_000_000.0,
        specific_volume: 0.223_26,
        specific_enthalpy: 3_918_200.0,
        specific_entropy: 7950.9,
    },
    // --- 5 MPa (T_sat ≈ 537.09 K / 263.94°C) ---
    SuperheatedEntry {
        temperature: 573.15,
        pressure: 5_000_000.0,
        specific_volume: 0.045_347,
        specific_enthalpy: 2_925_600.0,
        specific_entropy: 6210.9,
    },
    SuperheatedEntry {
        temperature: 623.15,
        pressure: 5_000_000.0,
        specific_volume: 0.051_971,
        specific_enthalpy: 3_069_300.0,
        specific_entropy: 6451.5,
    },
    SuperheatedEntry {
        temperature: 673.15,
        pressure: 5_000_000.0,
        specific_volume: 0.057_840,
        specific_enthalpy: 3_196_600.0,
        specific_entropy: 6648.1,
    },
    SuperheatedEntry {
        temperature: 723.15,
        pressure: 5_000_000.0,
        specific_volume: 0.063_325,
        specific_enthalpy: 3_317_000.0,
        specific_entropy: 6820.8,
    },
    SuperheatedEntry {
        temperature: 773.15,
        pressure: 5_000_000.0,
        specific_volume: 0.068_583,
        specific_enthalpy: 3_434_500.0,
        specific_entropy: 6977.8,
    },
    SuperheatedEntry {
        temperature: 873.15,
        pressure: 5_000_000.0,
        specific_volume: 0.078_703,
        specific_enthalpy: 3_666_800.0,
        specific_entropy: 7260.4,
    },
    SuperheatedEntry {
        temperature: 973.15,
        pressure: 5_000_000.0,
        specific_volume: 0.088_515,
        specific_enthalpy: 3_900_500.0,
        specific_entropy: 7513.7,
    },
    SuperheatedEntry {
        temperature: 1073.15,
        pressure: 5_000_000.0,
        specific_volume: 0.098_151,
        specific_enthalpy: 4_137_900.0,
        specific_entropy: 7745.9,
    },
    // --- 10 MPa (T_sat ≈ 584.15 K / 311.00°C) ---
    SuperheatedEntry {
        temperature: 623.15,
        pressure: 10_000_000.0,
        specific_volume: 0.022_442,
        specific_enthalpy: 2_924_000.0,
        specific_entropy: 5945.8,
    },
    SuperheatedEntry {
        temperature: 673.15,
        pressure: 10_000_000.0,
        specific_volume: 0.026_439,
        specific_enthalpy: 3_097_400.0,
        specific_entropy: 6213.9,
    },
    SuperheatedEntry {
        temperature: 723.15,
        pressure: 10_000_000.0,
        specific_volume: 0.029_785,
        specific_enthalpy: 3_242_300.0,
        specific_entropy: 6421.7,
    },
    SuperheatedEntry {
        temperature: 773.15,
        pressure: 10_000_000.0,
        specific_volume: 0.032_813,
        specific_enthalpy: 3_375_100.0,
        specific_entropy: 6599.3,
    },
    SuperheatedEntry {
        temperature: 823.15,
        pressure: 10_000_000.0,
        specific_volume: 0.035_655,
        specific_enthalpy: 3_501_900.0,
        specific_entropy: 6758.4,
    },
    SuperheatedEntry {
        temperature: 873.15,
        pressure: 10_000_000.0,
        specific_volume: 0.038_377,
        specific_enthalpy: 3_625_800.0,
        specific_entropy: 6904.5,
    },
    SuperheatedEntry {
        temperature: 973.15,
        pressure: 10_000_000.0,
        specific_volume: 0.043_594,
        specific_enthalpy: 3_870_300.0,
        specific_entropy: 7169.6,
    },
    SuperheatedEntry {
        temperature: 1073.15,
        pressure: 10_000_000.0,
        specific_volume: 0.048_624,
        specific_enthalpy: 4_114_700.0,
        specific_entropy: 7408.7,
    },
];

// The flat table must hold exactly one block of rows per pressure tier.
const _: () = assert!(SUPERHEATED_DATA.len() == SUPERHEAT_PRESSURES.len() * TEMPS_PER_TIER);

/// Look up superheated steam properties by temperature (K) and pressure (Pa).
///
/// Within a pressure tier, properties are interpolated linearly in T, using the
/// saturated-vapor state at that pressure as the lowest point so the whole
/// range from T_sat upward is covered. Between tiers, h and s are interpolated
/// linearly in ln P and v log-linearly in ln P (v ∝ 1/P for a near-ideal gas,
/// s ∝ −R·ln P), which is far more accurate than interpolating linearly in P.
///
/// At a tier pressure the temperature must be at or above saturation; between
/// tiers it must be at or above the saturation temperature of the higher tier.
pub fn superheated_lookup(temperature: f64, pressure: f64) -> Result<SuperheatedEntry> {
    if !finite_gt(temperature, 0.0) {
        return Err(UshmaError::InvalidTemperature {
            kelvin: temperature,
        });
    }
    if !finite_gt(pressure, 0.0) {
        return Err(UshmaError::InvalidPressure { pascals: pressure });
    }

    let out_of_range = || UshmaError::SuperheatOutOfRange {
        temperature,
        pressure,
    };
    let pressures = SUPERHEAT_PRESSURES;
    if pressure < pressures[0] || pressure > pressures[pressures.len() - 1] {
        return Err(out_of_range());
    }

    // First tier with p ≥ pressure (in range because of the bounds check above).
    let hi = pressures.partition_point(|&p| p < pressure);
    let p_hi = pressures[hi];
    if p_hi <= pressure {
        // Exactly on a tier: use that tier alone.
        let mut e = interp_tier(hi, temperature).ok_or_else(out_of_range)?;
        e.pressure = pressure;
        return Ok(e);
    }

    let lo = hi - 1;
    let p_lo = pressures[lo];
    let a = interp_tier(lo, temperature).ok_or_else(out_of_range)?;
    let b = interp_tier(hi, temperature).ok_or_else(out_of_range)?;
    let frac = (pressure / p_lo).ln() / (p_hi / p_lo).ln();

    let l = hisab::calc::lerp;
    Ok(SuperheatedEntry {
        temperature,
        pressure,
        specific_volume: l(a.specific_volume.ln(), b.specific_volume.ln(), frac).exp(),
        specific_enthalpy: l(a.specific_enthalpy, b.specific_enthalpy, frac),
        specific_entropy: l(a.specific_entropy, b.specific_entropy, frac),
    })
}

/// Find the superheated state at `pressure` (Pa) whose specific entropy equals
/// `entropy` (J/(kg·K)), e.g. the exit of an isentropic compressor or turbine.
///
/// Solves s(T, P) = s exactly for the table's interpolation model (s is
/// piecewise linear in T) over the temperature range that
/// [`superheated_lookup`] covers at this pressure, using ~7 table lookups. Returns
/// [`UshmaError::SuperheatOutOfRange`] (reporting the violated temperature
/// bound) when the requested entropy lies outside that range — in particular
/// when it is below s_g, i.e. the state is actually wet.
#[tracing::instrument(level = "debug")]
pub fn superheated_by_entropy(pressure: f64, entropy: f64) -> Result<SuperheatedEntry> {
    if !finite_gt(pressure, 0.0) {
        return Err(UshmaError::InvalidPressure { pascals: pressure });
    }
    if !entropy.is_finite() {
        return Err(invalid_param(format_args!(
            "entropy must be finite, got {entropy}"
        )));
    }
    let pressures = SUPERHEAT_PRESSURES;
    if pressure < pressures[0] || pressure > pressures[pressures.len() - 1] {
        return Err(UshmaError::SuperheatOutOfRange {
            temperature: f64::NAN,
            pressure,
        });
    }

    // s(T) at constant P is piecewise linear in T (linear in T within each tier,
    // and a fixed-weight blend of two tiers between tier pressures), with
    // breakpoints at the tabulated temperatures and the saturated-vapor anchor.
    // So the inverse is exact: bracket by binary search over the breakpoints,
    // then solve the linear segment.
    let hi = pressures.partition_point(|&p| p < pressure);
    let on_tier = pressures[hi] <= pressure;
    let t_min = saturated_by_pressure(pressures[hi])?.temperature;
    let rows = |tier: usize| &SUPERHEATED_DATA[tier * TEMPS_PER_TIER..(tier + 1) * TEMPS_PER_TIER];

    let mut t_max = rows(hi)[TEMPS_PER_TIER - 1].temperature;
    if !on_tier {
        t_max = t_max.min(rows(hi - 1)[TEMPS_PER_TIER - 1].temperature);
    }
    let mut bps = [0.0_f64; 2 * TEMPS_PER_TIER + 2];
    let mut nb = 0;
    let mut push = |t: f64| {
        if t > t_min && t < t_max {
            bps[nb] = t;
            nb += 1;
        }
    };
    rows(hi).iter().for_each(|e| push(e.temperature));
    if !on_tier {
        rows(hi - 1).iter().for_each(|e| push(e.temperature));
    }
    bps[nb] = t_min;
    bps[nb + 1] = t_max;
    nb += 2;
    let bps = &mut bps[..nb];
    bps.sort_unstable_by(f64::total_cmp);

    let s_at = |t: f64| superheated_lookup(t, pressure).map(|e| e.specific_entropy);
    let (mut lo, mut hi_i) = (0, bps.len() - 1);
    let (mut s_lo, mut s_hi) = (s_at(bps[lo])?, s_at(bps[hi_i])?);
    if entropy < s_lo || entropy > s_hi {
        return Err(UshmaError::SuperheatOutOfRange {
            temperature: if entropy < s_lo { t_min } else { t_max },
            pressure,
        });
    }
    while hi_i - lo > 1 {
        let mid = (lo + hi_i) / 2;
        let s_mid = s_at(bps[mid])?;
        if s_mid <= entropy {
            lo = mid;
            s_lo = s_mid;
        } else {
            hi_i = mid;
            s_hi = s_mid;
        }
    }
    let ds = s_hi - s_lo;
    let t = if ds > 0.0 {
        bps[lo] + (entropy - s_lo) / ds * (bps[hi_i] - bps[lo])
    } else {
        bps[lo]
    };
    superheated_lookup(t.clamp(t_min, t_max), pressure)
}

/// Interpolate within pressure tier `tier` at `temperature`.
///
/// The saturated-vapor state at the tier pressure serves as the lowest point.
/// Returns `None` outside [T_sat, T_max] for the tier.
fn interp_tier(tier: usize, temperature: f64) -> Option<SuperheatedEntry> {
    let rows = &SUPERHEATED_DATA[tier * TEMPS_PER_TIER..(tier + 1) * TEMPS_PER_TIER];
    let first = &rows[0];
    let last = &rows[rows.len() - 1];

    if temperature > last.temperature {
        return None;
    }
    if temperature < first.temperature {
        let sat = saturated_by_pressure(first.pressure).ok()?;
        if temperature < sat.temperature {
            return None;
        }
        let anchor = SuperheatedEntry {
            temperature: sat.temperature,
            pressure: first.pressure,
            specific_volume: sat.v_g,
            specific_enthalpy: sat.h_g,
            specific_entropy: sat.s_g,
        };
        return Some(lerp_superheated(&anchor, first, temperature));
    }

    let hi = rows
        .partition_point(|e| e.temperature <= temperature)
        .clamp(1, rows.len() - 1);
    Some(lerp_superheated(&rows[hi - 1], &rows[hi], temperature))
}

/// Linear interpolation in temperature between two rows of the same tier.
#[inline]
fn lerp_superheated(a: &SuperheatedEntry, b: &SuperheatedEntry, t: f64) -> SuperheatedEntry {
    let frac = (t - a.temperature) / (b.temperature - a.temperature);
    let l = hisab::calc::lerp;
    SuperheatedEntry {
        temperature: t,
        pressure: a.pressure,
        specific_volume: l(a.specific_volume, b.specific_volume, frac),
        specific_enthalpy: l(a.specific_enthalpy, b.specific_enthalpy, frac),
        specific_entropy: l(a.specific_entropy, b.specific_entropy, frac),
    }
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

    // --- Regression tests from the 1.4.0 audit ---

    #[test]
    fn test_saturation_pressure_if97_verification_values() {
        // IAPWS-IF97 Table 35 verification values.
        for (t, p_mpa) in [
            (300.0, 0.353_658_941e-2),
            (500.0, 2.638_897_76),
            (600.0, 12.344_314_6),
        ] {
            let p = saturation_pressure(t).unwrap();
            assert!((p / (p_mpa * 1e6) - 1.0).abs() < 1e-8, "T={t}: {p}");
        }
    }

    #[test]
    fn test_saturation_temperature_if97_verification_values() {
        // IAPWS-IF97 Table 36 verification values.
        for (p_mpa, t_ref) in [(0.1, 372.755_919), (1.0, 453.035_632), (10.0, 584.149_488)] {
            let t = saturation_temperature(p_mpa * 1e6).unwrap();
            assert!((t - t_ref).abs() < 1e-5, "P={p_mpa} MPa: {t}");
        }
    }

    #[test]
    fn test_saturation_line_roundtrip_and_range() {
        for t in [273.16, 300.0, 400.0, 500.0, 600.0, 647.0] {
            let p = saturation_pressure(t).unwrap();
            assert!((saturation_temperature(p).unwrap() - t).abs() < 1e-6);
        }
        assert!(saturation_pressure(273.0).is_err());
        assert!(saturation_pressure(650.0).is_err());
        assert!(saturation_pressure(f64::NAN).is_err());
        assert!(saturation_temperature(f64::NAN).is_err());
        assert!(saturation_temperature(3e7).is_err());
    }

    #[test]
    fn test_saturated_pressure_between_entries_matches_if97() {
        // Linear-in-T interpolation used to give 5815 Pa here (+3.3 %).
        let e = saturated_by_temperature(308.15).unwrap();
        assert!((e.pressure / 5628.62 - 1.0).abs() < 1e-5, "{}", e.pressure);
        let e = saturated_by_temperature(368.15).unwrap();
        assert!((e.pressure / 84_608.9 - 1.0).abs() < 1e-5, "{}", e.pressure);
    }

    #[test]
    fn test_saturated_by_pressure_roundtrip() {
        for t in [280.0, 333.0, 420.0, 555.0, 640.0] {
            let p = saturated_by_temperature(t).unwrap().pressure;
            let back = saturated_by_pressure(p).unwrap();
            assert!((back.temperature - t).abs() < 1e-6);
            assert!((back.pressure - p).abs() < 1e-9 * p);
        }
    }

    #[test]
    fn test_saturated_table_enthalpy_identity_exact() {
        for e in SATURATED_TABLE {
            assert!((e.h_g - e.h_f - e.h_fg).abs() < 1e-6, "T={}", e.temperature);
        }
        // s_g − s_f = h_fg / T (Clapeyron), to table rounding.
        for e in &SATURATED_TABLE[..SATURATED_TABLE.len() - 1] {
            let lhs = e.s_g - e.s_f;
            let rhs = e.h_fg / e.temperature;
            assert!(
                (lhs - rhs).abs() < 0.5,
                "T={}: {lhs} vs {rhs}",
                e.temperature
            );
        }
    }

    #[test]
    fn test_quality_of_saturated_vapor_is_one_between_entries() {
        // Inconsistent h_fg used to make x(h_g) = 1.0005 → InvalidQuality.
        for t in [478.0, 488.0, 528.0, 598.0] {
            let e = saturated_by_temperature(t).unwrap();
            let x = quality_from_enthalpy(e.h_g, &e).unwrap();
            assert!((x - 1.0).abs() < 1e-9);
        }
    }

    #[test]
    fn test_critical_point_specific_volume() {
        let e = saturated_by_temperature(647.096).unwrap();
        assert!((e.v_f - 1.0 / 322.0).abs() < 1e-7);
        assert!((e.pressure - 22_064_000.0).abs() < 1.0);
    }

    #[test]
    fn test_superheated_exact_tier_entries() {
        // These are table rows; the lookup used to reject them because it also
        // required the next tier up to cover T.
        for (t, p, v) in [
            (423.15, 200_000.0, 0.959_64),
            (373.15, 100_000.0, 1.695_9),
            (473.15, 1_000_000.0, 0.206_00),
            (523.15, 2_000_000.0, 0.111_50),
            (573.15, 5_000_000.0, 0.045_35),
        ] {
            let e = superheated_lookup(t, p).unwrap();
            assert!(
                (e.specific_volume / v - 1.0).abs() < 2e-3,
                "({t}, {p}): {}",
                e.specific_volume
            );
            assert!((e.pressure - p).abs() < 1e-9);
        }
    }

    #[test]
    fn test_superheated_between_tiers_matches_if97() {
        // (T, P, v, h, s) from IAPWS-IF97 Region 2.
        for (t, p, v, h, s) in [
            (573.15, 30_000.0, 8.8111, 3_076_248.0, 8775.0),
            (573.15, 150_000.0, 1.757_11, 3_073_315.0, 8028.4),
            (773.15, 7_500_000.0, 0.044_75, 3_405_334.0, 6762.0),
            (673.15, 300_000.0, 1.031_54, 3_275_424.0, 8034.6),
        ] {
            let e = superheated_lookup(t, p).unwrap();
            assert!(
                (e.specific_volume / v - 1.0).abs() < 0.01,
                "v({t},{p}) = {}",
                e.specific_volume
            );
            // 0.2 %: limited by the coarse 5 → 10 MPa tier spacing.
            assert!(
                (e.specific_enthalpy / h - 1.0).abs() < 2e-3,
                "h({t},{p}) = {}",
                e.specific_enthalpy
            );
            assert!(
                (e.specific_entropy - s).abs() < 15.0,
                "s({t},{p}) = {}",
                e.specific_entropy
            );
        }
    }

    #[test]
    fn test_superheated_near_saturation_uses_saturated_anchor() {
        // 10 kPa: T_sat = 318.96 K, first tabulated row at 323.15 K.
        let e = superheated_lookup(320.0, 10_000.0).unwrap();
        assert!((e.specific_volume - 14.72).abs() < 0.05);
        assert!((e.specific_enthalpy - 2_585_908.0).abs() < 1_000.0);
        assert!((e.specific_entropy - 8155.2).abs() < 2.0);
        // Below saturation is still rejected.
        assert!(superheated_lookup(318.0, 10_000.0).is_err());
    }

    #[test]
    fn test_superheated_out_of_range_reports_requested_pressure() {
        match superheated_lookup(300.0, 150_000.0) {
            Err(UshmaError::SuperheatOutOfRange { pressure, .. }) => {
                assert!((pressure - 150_000.0).abs() < 1e-9);
            }
            other => panic!("expected SuperheatOutOfRange, got {other:?}"),
        }
    }

    #[test]
    fn test_superheated_by_entropy_roundtrip() {
        for (t, p) in [(514.0, 100_000.0), (700.0, 1_000_000.0), (600.0, 30_000.0)] {
            let s = superheated_lookup(t, p).unwrap().specific_entropy;
            let e = superheated_by_entropy(p, s).unwrap();
            assert!((e.temperature - t).abs() < 1e-6, "{} vs {t}", e.temperature);
        }
        // IF97: s = 8.0 kJ/(kg·K) at 0.1 MPa → T = 514.13 K.
        let e = superheated_by_entropy(100_000.0, 8000.0).unwrap();
        assert!((e.temperature - 514.13).abs() < 0.5, "{}", e.temperature);
    }

    #[test]
    fn test_superheated_by_entropy_out_of_range() {
        // Below s_g at 100 kPa → wet, not superheated.
        assert!(matches!(
            superheated_by_entropy(100_000.0, 7000.0),
            Err(UshmaError::SuperheatOutOfRange { .. })
        ));
        assert!(superheated_by_entropy(100_000.0, 20_000.0).is_err());
        assert!(superheated_by_entropy(100_000.0, f64::NAN).is_err());
        assert!(superheated_by_entropy(f64::NAN, 8000.0).is_err());
        assert!(superheated_by_entropy(1.0, 8000.0).is_err());
    }

    #[test]
    fn test_steam_lookups_reject_nan() {
        assert!(saturated_by_temperature(f64::NAN).is_err());
        assert!(saturated_by_pressure(f64::NAN).is_err());
        assert!(superheated_lookup(f64::NAN, 100_000.0).is_err());
        assert!(superheated_lookup(500.0, f64::NAN).is_err());
        assert!(superheated_lookup(f64::INFINITY, 100_000.0).is_err());
    }
}
