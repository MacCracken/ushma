//! `transfer.tsv` (M1c onward): every public function of Rust 1.4.0's `transfer`, over seeded physical grids and the
//! edges every argument can take.
//!
//! For each function: a seeded sweep over its arguments' physical ranges, then, for each argument in turn, a physical
//! base case with that argument replaced by each edge value (0, −0, −1, NaN, ±∞, the 1e-30 divisor threshold and its
//! neighbours, 1e-300, 1e300), then each branch switch the function has. Outputs are the f64 result or
//! `E:<Variant>`. Classes: B for + − × ÷ √ and `powi`; C where `exp`, `ln`, `powf`, `tanh` or `atan` enters, compared
//! against the corrected oracle (`transfer.cr.tsv`, ADR 0004).

use ushma::transfer as t;

use crate::grids::Rng;
use crate::sweep::{around, edges, res, row, tuples};
use crate::tsv::{Class, Module, Tok, f};

fn part1(m: &mut Module) {
    // conduction(k, A, T_hot, T_cold, L)
    for a in tuples(
        &[
            vec![401.0, 1.0, 400.0, 300.0, 0.1],
            vec![0.0263, 2.5, 293.15, 313.15, 0.005],
        ],
        400,
        1,
        &|r| {
            vec![
                r.log_uniform(-3.0, 3.5),
                r.log_uniform(-4.0, 2.0),
                r.uniform(1.0, 2000.0),
                r.uniform(1.0, 2000.0),
                r.signed_log_uniform(-6.0, 1.0),
            ]
        },
    ) {
        row(
            m,
            "transfer::conduction",
            Class::B,
            &a,
            res(t::conduction(a[0], a[1], a[2], a[3], a[4])),
        );
    }
    // convection(h, A, T_s, T_f)
    for a in tuples(&[vec![25.0, 2.0, 350.0, 300.0]], 400, 2, &|r| {
        vec![
            r.log_uniform(0.0, 5.0),
            r.log_uniform(-4.0, 2.0),
            r.uniform(1.0, 2000.0),
            r.uniform(1.0, 2000.0),
        ]
    }) {
        row(
            m,
            "transfer::convection",
            Class::B,
            &a,
            vec![f(t::convection(a[0], a[1], a[2], a[3]))],
        );
    }
    // radiation(ε, A, T_s, T_surr)
    let mut rad = tuples(&[vec![0.9, 1.0, 500.0, 300.0]], 400, 3, &|r| {
        vec![
            r.unit(),
            r.log_uniform(-4.0, 2.0),
            r.uniform(0.0, 3000.0),
            r.uniform(0.0, 3000.0),
        ]
    });
    rad.extend([
        vec![1.0, 1.0, 300.0, 300.0],
        vec![1.0 + f64::EPSILON, 1.0, 300.0, 0.0],
        vec![-0.0, 1.0, 1e77, 0.0],
    ]);
    for a in rad {
        row(
            m,
            "transfer::radiation",
            Class::B,
            &a,
            res(t::radiation(a[0], a[1], a[2], a[3])),
        );
    }
    // thermal_resistance_conduction(k, A, L)
    for a in tuples(&[vec![401.0, 1.0, 0.1]], 300, 4, &|r| {
        vec![
            r.log_uniform(-3.0, 3.5),
            r.log_uniform(-4.0, 2.0),
            r.log_uniform(-6.0, 1.0),
        ]
    }) {
        row(
            m,
            "transfer::thermal_resistance_conduction",
            Class::B,
            &a,
            res(t::thermal_resistance_conduction(a[0], a[1], a[2])),
        );
    }
    // thermal_resistance_convection(h, A)
    for a in tuples(&[vec![25.0, 2.0]], 300, 5, &|r| {
        vec![r.log_uniform(-1.0, 5.0), r.log_uniform(-4.0, 2.0)]
    }) {
        row(
            m,
            "transfer::thermal_resistance_convection",
            Class::B,
            &a,
            res(t::thermal_resistance_convection(a[0], a[1])),
        );
    }
    // thermal_resistance_series / parallel over lists (`v:n` then the elements)
    let mut lists: Vec<Vec<f64>> = vec![
        vec![],
        vec![1.0],
        vec![0.5, 0.25, 0.125],
        vec![1e300, 1e300],
        vec![-0.0],
        vec![0.0, 0.0],
    ];
    for e in edges() {
        lists.push(vec![2.0, e, 3.0]);
    }
    let mut r = Rng::new(6);
    for _ in 0..200 {
        let n = r.int(1, 8) as usize;
        lists.push((0..n).map(|_| r.log_uniform(-4.0, 4.0)).collect());
    }
    for l in &lists {
        let mut ins = vec![Tok::V(l.len())];
        ins.extend(l.iter().map(|&x| f(x)));
        m.row(
            "transfer::thermal_resistance_series",
            Class::B,
            &ins,
            &[f(t::thermal_resistance_series(l))],
        );
        m.row(
            "transfer::thermal_resistance_parallel",
            Class::B,
            &ins,
            &res(t::thermal_resistance_parallel(l)),
        );
    }
    // heat_stored(m, c, ΔT)
    for a in tuples(&[vec![2.0, 4186.0, 25.0]], 300, 7, &|r| {
        vec![
            r.log_uniform(-3.0, 4.0),
            r.log_uniform(2.0, 4.0),
            r.signed_log_uniform(-2.0, 3.0),
        ]
    }) {
        row(
            m,
            "transfer::heat_stored",
            Class::B,
            &a,
            vec![f(t::heat_stored(a[0], a[1], a[2]))],
        );
    }
    // thermal_diffusivity(k, ρ, c_p)
    for a in tuples(&[vec![401.0, 8960.0, 385.0]], 300, 8, &|r| {
        vec![
            r.log_uniform(-2.0, 3.5),
            r.log_uniform(0.0, 4.5),
            r.log_uniform(2.0, 4.0),
        ]
    }) {
        row(
            m,
            "transfer::thermal_diffusivity",
            Class::B,
            &a,
            res(t::thermal_diffusivity(a[0], a[1], a[2])),
        );
    }
    // biot_number(h, L, k)
    for a in tuples(&[vec![25.0, 0.01, 50.0]], 300, 9, &|r| {
        vec![
            r.log_uniform(0.0, 4.0),
            r.log_uniform(-4.0, 0.0),
            r.log_uniform(-2.0, 3.5),
        ]
    }) {
        row(
            m,
            "transfer::biot_number",
            Class::B,
            &a,
            res(t::biot_number(a[0], a[1], a[2])),
        );
    }
    // lumped_capacitance(T0, T_env, t, τ)
    for a in tuples(&[vec![500.0, 300.0, 60.0, 120.0]], 400, 10, &|r| {
        vec![
            r.uniform(1.0, 2000.0),
            r.uniform(1.0, 2000.0),
            r.log_uniform(-2.0, 5.0),
            r.log_uniform(-1.0, 4.0),
        ]
    }) {
        row(
            m,
            "transfer::lumped_capacitance",
            Class::C,
            &a,
            vec![f(t::lumped_capacitance(a[0], a[1], a[2], a[3]))],
        );
    }
    // reynolds_number(V, L, ν)
    for a in tuples(&[vec![2.0, 0.05, 1.5e-5]], 300, 11, &|r| {
        vec![
            r.log_uniform(-3.0, 2.0),
            r.log_uniform(-4.0, 1.0),
            r.log_uniform(-7.0, -3.0),
        ]
    }) {
        row(
            m,
            "transfer::reynolds_number",
            Class::B,
            &a,
            res(t::reynolds_number(a[0], a[1], a[2])),
        );
    }
    // prandtl_number(ν, α)
    for a in tuples(&[vec![1.5e-5, 2.2e-5]], 300, 12, &|r| {
        vec![r.log_uniform(-7.0, -3.0), r.log_uniform(-8.0, -3.0)]
    }) {
        row(
            m,
            "transfer::prandtl_number",
            Class::B,
            &a,
            res(t::prandtl_number(a[0], a[1])),
        );
    }
    // nusselt_number(h, L, k)
    for a in tuples(&[vec![25.0, 0.1, 0.026]], 300, 13, &|r| {
        vec![
            r.log_uniform(0.0, 4.0),
            r.log_uniform(-4.0, 1.0),
            r.log_uniform(-2.0, 3.0),
        ]
    }) {
        row(
            m,
            "transfer::nusselt_number",
            Class::B,
            &a,
            res(t::nusselt_number(a[0], a[1], a[2])),
        );
    }
    // nusselt_dittus_boelter(Re, Pr)
    for a in tuples(&[vec![50_000.0, 0.7]], 600, 14, &|r| {
        vec![r.log_uniform(3.0, 7.0), r.log_uniform(-1.0, 3.0)]
    }) {
        row(
            m,
            "transfer::nusselt_dittus_boelter",
            Class::C,
            &a,
            res(t::nusselt_dittus_boelter(a[0], a[1])),
        );
    }
    // nusselt_natural_vertical(Ra, Pr)
    for a in tuples(&[vec![1e9, 0.71], vec![0.0, 7.0]], 600, 15, &|r| {
        vec![r.log_uniform(-2.0, 14.0), r.log_uniform(-2.0, 3.0)]
    }) {
        row(
            m,
            "transfer::nusselt_natural_vertical",
            Class::C,
            &a,
            res(t::nusselt_natural_vertical(a[0], a[1])),
        );
    }
    // boundary_layer_thickness(x, Re_x)
    for a in tuples(&[vec![0.5, 1e5]], 300, 16, &|r| {
        vec![r.log_uniform(-3.0, 1.0), r.log_uniform(1.0, 6.0)]
    }) {
        row(
            m,
            "transfer::boundary_layer_thickness",
            Class::B,
            &a,
            res(t::boundary_layer_thickness(a[0], a[1])),
        );
    }
    // thermal_boundary_layer(x, Re_x, Pr)
    for a in tuples(&[vec![0.5, 1e5, 0.7]], 600, 17, &|r| {
        vec![
            r.log_uniform(-3.0, 1.0),
            r.log_uniform(1.0, 6.0),
            r.log_uniform(-2.0, 3.0),
        ]
    }) {
        row(
            m,
            "transfer::thermal_boundary_layer",
            Class::C,
            &a,
            res(t::thermal_boundary_layer(a[0], a[1], a[2])),
        );
    }
}

/// The ε-NTU grid: NTU 0…100 × C_r 0…1, with C_r at both branch switches (1e-30, and 1 ± 1e-10) and their
/// neighbours.
fn ntu_grid() -> Vec<Vec<f64>> {
    let ntus = [
        0.0, -0.0, 1e-300, 1e-12, 0.001, 0.01, 0.1, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0,
        50.0, 100.0,
    ];
    let mut crs = vec![0.0, -0.0, 1e-31, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0];
    crs.extend(around(1e-30, 2));
    crs.extend(around(1.0 - 1e-10, 2));
    crs.extend(around(1.0 - 2e-10, 1));
    crs.push(f64::from_bits(1.0_f64.to_bits() - 1));
    let mut v = Vec::new();
    for &n in &ntus {
        for &c in &crs {
            v.push(vec![n, c]);
        }
    }
    v
}

/// LMTD's inputs (ΔT₁, ΔT₂, 0, 0), so the differences are ΔT₁ and ΔT₂ exactly: equal, 1 ulp apart, and across the
/// |ΔT₁/ΔT₂ − 1| < 1e-10 switch (with ΔT₂ = 1 the ratio is ΔT₁ exactly and its distance from 1 is exact).
fn lmtd_switch() -> Vec<Vec<f64>> {
    let mut v = Vec::new();
    for dt2 in [1.0, 50.0, 1e-3, 1e6] {
        for dt1 in around(dt2, 2) {
            v.push(vec![dt1, dt2, 0.0, 0.0]);
        }
    }
    for dt1 in around(1.0 + 1e-10, 3)
        .into_iter()
        .chain(around(1.0 - 1e-10, 3))
    {
        v.push(vec![dt1, 1.0, 0.0, 0.0]);
    }
    for k in -12_i32..=12 {
        v.push(vec![50.0 * (1.0 + f64::from(k) * 1e-11), 50.0, 0.0, 0.0]);
    }
    v
}

fn part2(m: &mut Module) {
    // fin_parameter(h, P, k, A_c)
    for a in tuples(&[vec![25.0, 0.04, 200.0, 1e-4]], 300, 18, &|r| {
        vec![
            r.log_uniform(0.0, 4.0),
            r.log_uniform(-3.0, 0.0),
            r.log_uniform(-1.0, 3.5),
            r.log_uniform(-7.0, -2.0),
        ]
    }) {
        row(
            m,
            "transfer::fin_parameter",
            Class::B,
            &a,
            res(t::fin_parameter(a[0], a[1], a[2], a[3])),
        );
    }
    // fin_rectangular_heat(h, P, k, A_c, L, T_base, T_fluid)
    for a in tuples(
        &[vec![25.0, 0.04, 200.0, 1e-4, 0.05, 373.15, 293.15]],
        600,
        19,
        &|r| {
            vec![
                r.log_uniform(0.0, 4.0),
                r.log_uniform(-3.0, 0.0),
                r.log_uniform(-1.0, 3.5),
                r.log_uniform(-7.0, -2.0),
                r.log_uniform(-3.0, 0.5),
                r.uniform(1.0, 1500.0),
                r.uniform(1.0, 1500.0),
            ]
        },
    ) {
        row(
            m,
            "transfer::fin_rectangular_heat",
            Class::C,
            &a,
            res(t::fin_rectangular_heat(
                a[0], a[1], a[2], a[3], a[4], a[5], a[6],
            )),
        );
    }
    // fin_efficiency_rectangular(m, L), with m·L at the 1e-30 switch and where tanh saturates
    let mut eff = tuples(&[vec![10.0, 0.05]], 400, 20, &|r| {
        vec![r.log_uniform(-2.0, 3.0), r.log_uniform(-4.0, 1.0)]
    });
    for ml in around(1e-30, 2) {
        eff.push(vec![ml, 1.0]);
    }
    eff.extend([
        vec![0.0, 0.0],
        vec![0.0, 5.0],
        vec![1e-20, 1e-20],
        vec![1e3, 10.0],
        vec![1e200, 1e200],
        vec![1e-160, 1e-160],
    ]);
    for a in eff {
        row(
            m,
            "transfer::fin_efficiency_rectangular",
            Class::C,
            &a,
            res(t::fin_efficiency_rectangular(a[0], a[1])),
        );
    }
    // fin_effectiveness(Q_fin, h, A_c, T_base, T_fluid)
    let mut fe = tuples(&[vec![5.0, 25.0, 1e-4, 373.15, 293.15]], 300, 21, &|r| {
        vec![
            r.log_uniform(-2.0, 3.0),
            r.log_uniform(0.0, 4.0),
            r.log_uniform(-7.0, -2.0),
            r.uniform(1.0, 1500.0),
            r.uniform(1.0, 1500.0),
        ]
    });
    fe.extend([
        vec![5.0, 25.0, 1e-4, 300.0, 300.0],
        vec![5.0, 1e-20, 1e-11, 301.0, 300.0],
        vec![5.0, 25.0, 1e-4, 293.15, 373.15],
    ]);
    for a in fe {
        row(
            m,
            "transfer::fin_effectiveness",
            Class::B,
            &a,
            res(t::fin_effectiveness(a[0], a[1], a[2], a[3], a[4])),
        );
    }
    // lmtd_parallel / lmtd_counter(T_h,in, T_h,out, T_c,in, T_c,out)
    let lm_draw = |r: &mut Rng| {
        let tc_in = r.uniform(250.0, 600.0);
        let tc_out = tc_in + r.uniform(0.0, 200.0);
        let th_out = tc_out + r.log_uniform(-1.0, 2.5);
        let th_in = th_out + r.uniform(0.0, 300.0);
        vec![th_in, th_out, tc_in, tc_out]
    };
    let mut lp = tuples(&[vec![400.0, 350.0, 300.0, 320.0]], 400, 22, &lm_draw);
    lp.extend(lmtd_switch());
    for a in lp {
        row(
            m,
            "transfer::lmtd_parallel",
            Class::C,
            &a,
            res(t::lmtd_parallel(a[0], a[1], a[2], a[3])),
        );
    }
    let mut lc = tuples(&[vec![400.0, 350.0, 300.0, 320.0]], 400, 23, &lm_draw);
    lc.extend(lmtd_switch());
    for a in lc {
        row(
            m,
            "transfer::lmtd_counter",
            Class::C,
            &a,
            res(t::lmtd_counter(a[0], a[1], a[2], a[3])),
        );
    }
    // heat_exchanger_lmtd(U, A, LMTD)
    for a in tuples(&[vec![500.0, 10.0, 25.0]], 300, 24, &|r| {
        vec![
            r.log_uniform(0.0, 4.0),
            r.log_uniform(-2.0, 3.0),
            r.log_uniform(-1.0, 2.5),
        ]
    }) {
        row(
            m,
            "transfer::heat_exchanger_lmtd",
            Class::B,
            &a,
            vec![f(t::heat_exchanger_lmtd(a[0], a[1], a[2]))],
        );
    }
    // ntu(U, A, C_min)
    for a in tuples(&[vec![500.0, 10.0, 2000.0]], 300, 25, &|r| {
        vec![
            r.log_uniform(0.0, 4.0),
            r.log_uniform(-2.0, 3.0),
            r.log_uniform(1.0, 5.0),
        ]
    }) {
        row(
            m,
            "transfer::ntu",
            Class::B,
            &a,
            res(t::ntu(a[0], a[1], a[2])),
        );
    }
    // effectiveness_parallel / effectiveness_counter(NTU, C_r)
    let ntu_draw = |r: &mut Rng| vec![r.log_uniform(-3.0, 2.0), r.unit()];
    let mut ep = tuples(&[vec![2.0, 0.5]], 600, 26, &ntu_draw);
    ep.extend(ntu_grid());
    for a in ep {
        row(
            m,
            "transfer::effectiveness_parallel",
            Class::C,
            &a,
            res(t::effectiveness_parallel(a[0], a[1])),
        );
    }
    let mut ec = tuples(&[vec![2.0, 0.5]], 600, 27, &ntu_draw);
    ec.extend(ntu_grid());
    for a in ec {
        row(
            m,
            "transfer::effectiveness_counter",
            Class::C,
            &a,
            res(t::effectiveness_counter(a[0], a[1])),
        );
    }
    // heat_exchanger_ntu(ε, C_min, T_h,in, T_c,in)
    for a in tuples(&[vec![0.7, 2000.0, 400.0, 300.0]], 300, 28, &|r| {
        vec![
            r.unit(),
            r.log_uniform(1.0, 5.0),
            r.uniform(250.0, 1500.0),
            r.uniform(250.0, 1500.0),
        ]
    }) {
        row(
            m,
            "transfer::heat_exchanger_ntu",
            Class::B,
            &a,
            vec![f(t::heat_exchanger_ntu(a[0], a[1], a[2], a[3]))],
        );
    }
    // view_factor_parallel_plates(W, H, D), with D = 1 so X = W and Y = H across the max(X, Y) < 0.01 switch
    let mut vp = tuples(&[vec![1.0, 1.0, 1.0], vec![2.0, 0.5, 1.0]], 600, 29, &|r| {
        vec![
            r.log_uniform(-4.0, 2.0),
            r.log_uniform(-4.0, 2.0),
            r.log_uniform(-2.0, 1.0),
        ]
    });
    let mut sw = around(0.01, 2);
    sw.extend([0.005, 0.0099, 0.0101, 0.02, 1e-6]);
    for &w in &sw {
        for &h in &sw {
            vp.push(vec![w, h, 1.0]);
        }
    }
    vp.extend([
        vec![1e-150, 1e-150, 1.0],
        vec![1e150, 1e150, 1.0],
        vec![1e6, 1e-6, 1.0],
    ]);
    for a in vp {
        row(
            m,
            "transfer::view_factor_parallel_plates",
            Class::C,
            &a,
            res(t::view_factor_parallel_plates(a[0], a[1], a[2])),
        );
    }
    // view_factor_perpendicular_plates(W, H₁, H₂), including aspect ratios where the clamp acts
    let mut vq = tuples(&[vec![1.0, 1.0, 1.0], vec![1.0, 2.0, 0.5]], 600, 30, &|r| {
        vec![
            r.log_uniform(-2.0, 2.0),
            r.log_uniform(-3.0, 3.0),
            r.log_uniform(-3.0, 3.0),
        ]
    });
    for w in [1e-8, 1e-4, 1e4, 1e8, 1e150] {
        vq.push(vec![1.0, w, 1.0]);
        vq.push(vec![1.0, 1.0, w]);
        vq.push(vec![1.0, w, w]);
    }
    for a in vq {
        row(
            m,
            "transfer::view_factor_perpendicular_plates",
            Class::C,
            &a,
            res(t::view_factor_perpendicular_plates(a[0], a[1], a[2])),
        );
    }
    // view_factor_coaxial_disks(r₁, r₂, d), including small disks far apart
    let mut vd = tuples(
        &[vec![1.0, 1.0, 1.0], vec![0.5, 1.0, 0.25]],
        400,
        31,
        &|r| {
            vec![
                r.log_uniform(-3.0, 1.0),
                r.log_uniform(-3.0, 1.0),
                r.log_uniform(-3.0, 1.0),
            ]
        },
    );
    vd.extend([
        vec![1e-6, 1e-6, 1e3],
        vec![1e-6, 2e-6, 1e3],
        vec![1e3, 1e-3, 1e-3],
        vec![1e150, 1e150, 1e-150],
    ]);
    for a in vd {
        row(
            m,
            "transfer::view_factor_coaxial_disks",
            Class::B,
            &a,
            res(t::view_factor_coaxial_disks(a[0], a[1], a[2])),
        );
    }
}

pub fn emit() -> Module {
    let mut m = Module::new("transfer");
    part1(&mut m);
    part2(&mut m);
    m
}
