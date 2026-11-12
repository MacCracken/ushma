# Rust-test ledger

> **Generated** by `python3 -I -B scripts/rust_test_ledger.py .`; `--check` gates it (port plan §6). Every test
> of the frozen Rust oracle (`rust-old/`, ushma 1.4.0) ports to a Cyrius test whose doc comment names it
> (`# rust <file>: <test>`), at its module's milestone. Nothing is deferred: a test not ported yet is one whose
> milestone hasn't come. The golden corpus (`tests/golden/`) checks every ported function on its grid besides.

**479 Rust tests:** 102 ported; 78 at M2, 95 at M3, 64 at M4, 62 at M5, 27 at M6, 51 at M7.

| Rust file | Tests | Ported | Due, unported | Later (milestone) |
|---|---:|---:|---:|---|
| `src/ai.rs` | 7 | 0 | 0 | 7 (M7) |
| `src/bridge.rs` | 31 | 0 | 0 | 31 (M7) |
| `src/chem.rs` | 23 | 0 | 0 | 23 (M5) |
| `src/cycle.rs` | 59 | 0 | 0 | 59 (M4) |
| `src/entropy.rs` | 23 | 0 | 0 | 23 (M2) |
| `src/error.rs` | 18 | 18 | 0 | — |
| `src/integration/soorat.rs` | 13 | 0 | 0 | 13 (M7) |
| `src/logging.rs` | 1 | 1 | 0 | — |
| `src/material.rs` | 8 | 8 | 0 | — |
| `src/numerical.rs` | 27 | 0 | 0 | 27 (M6) |
| `src/phase.rs` | 44 | 0 | 0 | 44 (M3) |
| `src/stat.rs` | 33 | 0 | 0 | 33 (M5) |
| `src/state.rs` | 47 | 0 | 0 | 47 (M2) |
| `src/steam.rs` | 46 | 0 | 0 | 46 (M3) |
| `src/transfer.rs` | 67 | 67 | 0 | — |
| `src/validate.rs` | 2 | 2 | 0 | — |
| `tests/integration.rs` | 30 | 6 | 0 | 8 (M2), 5 (M3), 5 (M4), 6 (M5) |
| **Total** | **479** | **102** | **0** | |

## Where each is ported

| Rust test | Cyrius test |
|---|---|
| `src/error.rs` `test_division_by_zero` | `tests/foundation.tcyr` `test_rust_error_variants` |
| `src/error.rs` `test_error_is_debug` | `tests/foundation.tcyr` `test_rust_error_variants` |
| `src/error.rs` `test_invalid_compression_ratio` | `tests/foundation.tcyr` `test_rust_error_variants` |
| `src/error.rs` `test_invalid_conductivity` | `tests/foundation.tcyr` `test_rust_error_variants` |
| `src/error.rs` `test_invalid_cutoff_ratio` | `tests/foundation.tcyr` `test_rust_error_variants` |
| `src/error.rs` `test_invalid_cycle_parameter` | `tests/foundation.tcyr` `test_rust_error_variants` |
| `src/error.rs` `test_invalid_parameter` | `tests/foundation.tcyr` `test_rust_error_variants` |
| `src/error.rs` `test_invalid_pressure` | `tests/foundation.tcyr` `test_rust_error_variants` |
| `src/error.rs` `test_invalid_pressure_ratio` | `tests/foundation.tcyr` `test_rust_error_variants` |
| `src/error.rs` `test_invalid_quality` | `tests/foundation.tcyr` `test_rust_error_variants` |
| `src/error.rs` `test_invalid_specific_heat` | `tests/foundation.tcyr` `test_rust_error_variants` |
| `src/error.rs` `test_invalid_temperature` | `tests/foundation.tcyr` `test_rust_error_variants` |
| `src/error.rs` `test_invalid_volume` | `tests/foundation.tcyr` `test_rust_error_variants` |
| `src/error.rs` `test_phase_transition` | `tests/foundation.tcyr` `test_rust_error_variants` |
| `src/error.rs` `test_result_alias` | `tests/foundation.tcyr` `test_rust_error_variants` |
| `src/error.rs` `test_steam_table_out_of_range` | `tests/foundation.tcyr` `test_rust_error_variants` |
| `src/error.rs` `test_substance_not_found` | `tests/foundation.tcyr` `test_rust_error_variants` |
| `src/error.rs` `test_superheat_out_of_range` | `tests/foundation.tcyr` `test_rust_error_variants` |
| `src/logging.rs` `init_does_not_panic` | `tests/foundation.tcyr` `test_log_init` |
| `src/material.rs` `test_air_low_conductivity` | `tests/material.tcyr` `test_rust_material` |
| `src/material.rs` `test_all_materials_count` | `tests/material.tcyr` `test_rust_material` |
| `src/material.rs` `test_copper_diffusivity` | `tests/material.tcyr` `test_rust_material` |
| `src/material.rs` `test_diamond_highest_conductivity` | `tests/material.tcyr` `test_rust_material` |
| `src/material.rs` `test_material_serde_roundtrip` | `tests/material.tcyr` `test_rust_material` |
| `src/material.rs` `test_volumetric_heat_capacity` | `tests/material.tcyr` `test_rust_material` |
| `src/material.rs` `test_water_phase_transitions` | `tests/material.tcyr` `test_rust_material` |
| `src/material.rs` `test_water_specific_heat` | `tests/material.tcyr` `test_rust_material` |
| `src/transfer.rs` `test_biot_number` | `tests/transfer.tcyr` `test_storage_diffusivity_biot` |
| `src/transfer.rs` `test_biot_number_zero_conductivity` | `tests/transfer.tcyr` `test_storage_diffusivity_biot` |
| `src/transfer.rs` `test_boundary_layer_invalid` | `tests/transfer.tcyr` `test_boundary_layers` |
| `src/transfer.rs` `test_boundary_layer_thickness` | `tests/transfer.tcyr` `test_boundary_layers` |
| `src/transfer.rs` `test_conduction_positive_flux` | `tests/transfer.tcyr` `test_conduction` |
| `src/transfer.rs` `test_conduction_rejects_nan_and_inf` | `tests/transfer.tcyr` `test_conduction` |
| `src/transfer.rs` `test_conduction_zero_gradient` | `tests/transfer.tcyr` `test_conduction` |
| `src/transfer.rs` `test_conduction_zero_thickness` | `tests/transfer.tcyr` `test_conduction` |
| `src/transfer.rs` `test_convection` | `tests/transfer.tcyr` `test_convection_radiation` |
| `src/transfer.rs` `test_dittus_boelter` | `tests/transfer.tcyr` `test_dimensionless` |
| `src/transfer.rs` `test_dittus_boelter_invalid` | `tests/transfer.tcyr` `test_dimensionless` |
| `src/transfer.rs` `test_effectiveness_counter_balanced` | `tests/transfer.tcyr` `test_entu` |
| `src/transfer.rs` `test_effectiveness_counter_gte_parallel` | `tests/transfer.tcyr` `test_entu` |
| `src/transfer.rs` `test_effectiveness_cr_zero` | `tests/transfer.tcyr` `test_entu` |
| `src/transfer.rs` `test_effectiveness_invalid` | `tests/transfer.tcyr` `test_entu` |
| `src/transfer.rs` `test_effectiveness_parallel_increases_with_ntu` | `tests/transfer.tcyr` `test_entu` |
| `src/transfer.rs` `test_fin_effectiveness_greater_than_one` | `tests/transfer.tcyr` `test_fins` |
| `src/transfer.rs` `test_fin_efficiency_long_fin` | `tests/transfer.tcyr` `test_fins` |
| `src/transfer.rs` `test_fin_efficiency_short_fin` | `tests/transfer.tcyr` `test_fins` |
| `src/transfer.rs` `test_fin_efficiency_zero_length` | `tests/transfer.tcyr` `test_fins` |
| `src/transfer.rs` `test_fin_parameter` | `tests/transfer.tcyr` `test_fins` |
| `src/transfer.rs` `test_fin_parameter_invalid` | `tests/transfer.tcyr` `test_fins` |
| `src/transfer.rs` `test_fin_rectangular_heat` | `tests/transfer.tcyr` `test_fins` |
| `src/transfer.rs` `test_heat_exchanger_lmtd_basic` | `tests/transfer.tcyr` `test_lmtd` |
| `src/transfer.rs` `test_heat_exchanger_ntu_basic` | `tests/transfer.tcyr` `test_entu` |
| `src/transfer.rs` `test_heat_stored` | `tests/transfer.tcyr` `test_storage_diffusivity_biot` |
| `src/transfer.rs` `test_heat_stored_negative_delta` | `tests/transfer.tcyr` `test_storage_diffusivity_biot` |
| `src/transfer.rs` `test_invalid_conductivity` | `tests/transfer.tcyr` `test_conduction` |
| `src/transfer.rs` `test_lmtd_counter` | `tests/transfer.tcyr` `test_lmtd` |
| `src/transfer.rs` `test_lmtd_counter_gte_parallel` | `tests/transfer.tcyr` `test_lmtd` |
| `src/transfer.rs` `test_lmtd_invalid` | `tests/transfer.tcyr` `test_lmtd` |
| `src/transfer.rs` `test_lmtd_parallel` | `tests/transfer.tcyr` `test_lmtd` |
| `src/transfer.rs` `test_lumped_capacitance` | `tests/transfer.tcyr` `test_storage_diffusivity_biot` |
| `src/transfer.rs` `test_ntu_basic` | `tests/transfer.tcyr` `test_entu` |
| `src/transfer.rs` `test_nusselt_natural_vertical` | `tests/transfer.tcyr` `test_dimensionless` |
| `src/transfer.rs` `test_nusselt_natural_vertical_zero_ra` | `tests/transfer.tcyr` `test_dimensionless` |
| `src/transfer.rs` `test_nusselt_number` | `tests/transfer.tcyr` `test_dimensionless` |
| `src/transfer.rs` `test_prandtl_number` | `tests/transfer.tcyr` `test_dimensionless` |
| `src/transfer.rs` `test_radiation` | `tests/transfer.tcyr` `test_convection_radiation` |
| `src/transfer.rs` `test_radiation_invalid_emissivity` | `tests/transfer.tcyr` `test_convection_radiation` |
| `src/transfer.rs` `test_radiation_nan_emissivity` | `tests/transfer.tcyr` `test_convection_radiation` |
| `src/transfer.rs` `test_radiation_negative_temp` | `tests/transfer.tcyr` `test_convection_radiation` |
| `src/transfer.rs` `test_reynolds_invalid` | `tests/transfer.tcyr` `test_dimensionless` |
| `src/transfer.rs` `test_reynolds_number` | `tests/transfer.tcyr` `test_dimensionless` |
| `src/transfer.rs` `test_thermal_boundary_layer` | `tests/transfer.tcyr` `test_boundary_layers` |
| `src/transfer.rs` `test_thermal_diffusivity` | `tests/transfer.tcyr` `test_storage_diffusivity_biot` |
| `src/transfer.rs` `test_thermal_diffusivity_zero_density` | `tests/transfer.tcyr` `test_storage_diffusivity_biot` |
| `src/transfer.rs` `test_thermal_resistance_convection_zero` | `tests/transfer.tcyr` `test_resistances` |
| `src/transfer.rs` `test_thermal_resistance_parallel` | `tests/transfer.tcyr` `test_resistances` |
| `src/transfer.rs` `test_thermal_resistance_parallel_single` | `tests/transfer.tcyr` `test_resistances` |
| `src/transfer.rs` `test_thermal_resistance_series` | `tests/transfer.tcyr` `test_resistances` |
| `src/transfer.rs` `test_transfer_rejects_nan` | `tests/transfer.tcyr` `test_convection_radiation`, `tests/transfer.tcyr` `test_dimensionless`, `tests/transfer.tcyr` `test_lmtd`, `tests/transfer.tcyr` `test_entu`, `tests/transfer.tcyr` `test_view_factors` |
| `src/transfer.rs` `test_transfer_validation_consistency` | `tests/transfer.tcyr` `test_resistances`, `tests/transfer.tcyr` `test_storage_diffusivity_biot`, `tests/transfer.tcyr` `test_dimensionless`, `tests/transfer.tcyr` `test_fins` |
| `src/transfer.rs` `test_view_factor_coaxial_disks_equal` | `tests/transfer.tcyr` `test_view_factors` |
| `src/transfer.rs` `test_view_factor_coaxial_disks_far` | `tests/transfer.tcyr` `test_view_factors` |
| `src/transfer.rs` `test_view_factor_coaxial_disks_range` | `tests/transfer.tcyr` `test_view_factors` |
| `src/transfer.rs` `test_view_factor_coaxial_small_disks` | `tests/transfer.tcyr` `test_view_factor_regressions` |
| `src/transfer.rs` `test_view_factor_invalid` | `tests/transfer.tcyr` `test_view_factors` |
| `src/transfer.rs` `test_view_factor_parallel_plates_close` | `tests/transfer.tcyr` `test_view_factors` |
| `src/transfer.rs` `test_view_factor_parallel_plates_far` | `tests/transfer.tcyr` `test_view_factors` |
| `src/transfer.rs` `test_view_factor_parallel_plates_invalid` | `tests/transfer.tcyr` `test_view_factors` |
| `src/transfer.rs` `test_view_factor_parallel_plates_range` | `tests/transfer.tcyr` `test_view_factors` |
| `src/transfer.rs` `test_view_factor_parallel_small_plates` | `tests/transfer.tcyr` `test_view_factor_regressions` |
| `src/transfer.rs` `test_view_factor_perpendicular_invalid` | `tests/transfer.tcyr` `test_view_factors` |
| `src/transfer.rs` `test_view_factor_perpendicular_range` | `tests/transfer.tcyr` `test_view_factors` |
| `src/transfer.rs` `test_view_factor_perpendicular_reciprocity` | `tests/transfer.tcyr` `test_view_factor_regressions` |
| `src/transfer.rs` `test_view_factor_perpendicular_reference_values` | `tests/transfer.tcyr` `test_view_factor_regressions` |
| `src/validate.rs` `finite_ge_accepts_bound` | `tests/foundation.tcyr` `test_validate_predicates` |
| `src/validate.rs` `finite_gt_rejects_nan_and_inf` | `tests/foundation.tcyr` `test_validate_predicates` |
| `tests/integration.rs` `heat_conduction_through_copper_wall` | `tests/integration.tcyr` `heat_conduction_through_copper_wall` |
| `tests/integration.rs` `lmtd_vs_ntu_same_heat_exchanger` | `tests/integration.tcyr` `lmtd_vs_ntu_same_heat_exchanger` |
| `tests/integration.rs` `material_thermal_diffusivity_matches_manual` | `tests/integration.tcyr` `material_thermal_diffusivity_matches_manual` |
| `tests/integration.rs` `nusselt_to_h_to_convection_chain` | `tests/integration.tcyr` `nusselt_to_h_to_convection_chain` |
| `tests/integration.rs` `thermal_equilibrium_lumped_capacitance` | `tests/integration.tcyr` `thermal_equilibrium_lumped_capacitance` |
| `tests/integration.rs` `view_factor_reciprocity` | `tests/integration.tcyr` `view_factor_reciprocity` |
