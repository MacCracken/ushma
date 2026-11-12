# Rust-test ledger

> **Generated** by `python3 -I -B scripts/rust_test_ledger.py .`; `--check` gates it (port plan §6). Every test
> of the frozen Rust oracle (`rust-old/`, ushma 1.4.0) ports to a Cyrius test whose doc comment names it
> (`# rust <file>: <test>`), at its module's milestone. Nothing is deferred: a test not ported yet is one whose
> milestone hasn't come. The golden corpus (`tests/golden/`) checks every ported function on its grid besides.

**479 Rust tests:** 260 ported; 15 at M3, 64 at M4, 62 at M5, 27 at M6, 51 at M7.

| Rust file | Tests | Ported | Due, unported | Later (milestone) |
|---|---:|---:|---:|---|
| `src/ai.rs` | 7 | 0 | 0 | 7 (M7) |
| `src/bridge.rs` | 31 | 0 | 0 | 31 (M7) |
| `src/chem.rs` | 23 | 0 | 0 | 23 (M5) |
| `src/cycle.rs` | 59 | 0 | 0 | 59 (M4) |
| `src/entropy.rs` | 23 | 23 | 0 | — |
| `src/error.rs` | 18 | 18 | 0 | — |
| `src/integration/soorat.rs` | 13 | 0 | 0 | 13 (M7) |
| `src/logging.rs` | 1 | 1 | 0 | — |
| `src/material.rs` | 8 | 8 | 0 | — |
| `src/numerical.rs` | 27 | 0 | 0 | 27 (M6) |
| `src/phase.rs` | 44 | 44 | 0 | — |
| `src/stat.rs` | 33 | 0 | 0 | 33 (M5) |
| `src/state.rs` | 47 | 47 | 0 | — |
| `src/steam.rs` | 46 | 32 | 0 | 14 (M3) |
| `src/transfer.rs` | 67 | 67 | 0 | — |
| `src/validate.rs` | 2 | 2 | 0 | — |
| `tests/integration.rs` | 30 | 18 | 0 | 1 (M3), 5 (M4), 6 (M5) |
| **Total** | **479** | **260** | **0** | |

## Where each is ported

| Rust test | Cyrius test |
|---|---|
| `src/entropy.rs` `test_carnot_cop` | `tests/entropy.tcyr` `test_carnot` |
| `src/entropy.rs` `test_carnot_cop_invalid` | `tests/entropy.tcyr` `test_carnot` |
| `src/entropy.rs` `test_carnot_efficiency` | `tests/entropy.tcyr` `test_carnot` |
| `src/entropy.rs` `test_carnot_invalid_temps` | `tests/entropy.tcyr` `test_carnot` |
| `src/entropy.rs` `test_entropy_functions_reject_nan` | `tests/entropy.tcyr` `test_heat_transfer`, `tests/entropy.tcyr` `test_carnot`, `tests/entropy.tcyr` `test_mixing` |
| `src/entropy.rs` `test_entropy_functions_validate_moles` | `tests/entropy.tcyr` `test_entropy_changes`, `tests/entropy.tcyr` `test_mixing` |
| `src/entropy.rs` `test_entropy_of_mixing` | `tests/entropy.tcyr` `test_mixing` |
| `src/entropy.rs` `test_entropy_of_mixing_bad_fractions` | `tests/entropy.tcyr` `test_mixing` |
| `src/entropy.rs` `test_entropy_of_mixing_pure_substance` | `tests/entropy.tcyr` `test_mixing` |
| `src/entropy.rs` `test_entropy_of_mixing_three_components` | `tests/entropy.tcyr` `test_mixing` |
| `src/entropy.rs` `test_entropy_of_mixing_zero_fraction` | `tests/entropy.tcyr` `test_mixing` |
| `src/entropy.rs` `test_gibbs` | `tests/entropy.tcyr` `test_potentials` |
| `src/entropy.rs` `test_heat_transfer_entropy` | `tests/entropy.tcyr` `test_heat_transfer` |
| `src/entropy.rs` `test_heat_transfer_entropy_negative_heat` | `tests/entropy.tcyr` `test_heat_transfer` |
| `src/entropy.rs` `test_heat_transfer_entropy_zero_temp` | `tests/entropy.tcyr` `test_heat_transfer` |
| `src/entropy.rs` `test_helmholtz` | `tests/entropy.tcyr` `test_potentials` |
| `src/entropy.rs` `test_ideal_gas_entropy_change` | `tests/entropy.tcyr` `test_entropy_changes` |
| `src/entropy.rs` `test_ideal_gas_entropy_zero_temp` | `tests/entropy.tcyr` `test_entropy_changes` |
| `src/entropy.rs` `test_ideal_gas_entropy_zero_volume` | `tests/entropy.tcyr` `test_entropy_changes` |
| `src/entropy.rs` `test_is_spontaneous_boundary` | `tests/entropy.tcyr` `test_potentials` |
| `src/entropy.rs` `test_isothermal_entropy_compression` | `tests/entropy.tcyr` `test_entropy_changes` |
| `src/entropy.rs` `test_isothermal_entropy_expansion` | `tests/entropy.tcyr` `test_entropy_changes` |
| `src/entropy.rs` `test_spontaneous` | `tests/entropy.tcyr` `test_potentials` |
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
| `src/phase.rs` `test_all_substances_count` | `tests/phase.tcyr` `test_table` |
| `src/phase.rs` `test_all_substances_have_valid_data` | `tests/phase.tcyr` `test_table` |
| `src/phase.rs` `test_ammonia_boiling_point` | `tests/phase.tcyr` `test_table` |
| `src/phase.rs` `test_clausius_clapeyron_pressure_invalid` | `tests/phase.tcyr` `test_clausius_clapeyron` |
| `src/phase.rs` `test_clausius_clapeyron_pressure_overflow_is_error` | `tests/phase.tcyr` `test_clausius_clapeyron` |
| `src/phase.rs` `test_clausius_clapeyron_pressure_roundtrip` | `tests/phase.tcyr` `test_clausius_clapeyron` |
| `src/phase.rs` `test_clausius_clapeyron_pressure_water` | `tests/phase.tcyr` `test_clausius_clapeyron` |
| `src/phase.rs` `test_clausius_clapeyron_slope_invalid` | `tests/phase.tcyr` `test_clausius_clapeyron` |
| `src/phase.rs` `test_clausius_clapeyron_slope_water` | `tests/phase.tcyr` `test_clausius_clapeyron` |
| `src/phase.rs` `test_co2_critical_point` | `tests/phase.tcyr` `test_table` |
| `src/phase.rs` `test_co2_heat_for_phase_change_uses_sublimation` | `tests/phase.tcyr` `test_latent_heats` |
| `src/phase.rs` `test_co2_phase_regions` | `tests/phase.tcyr` `test_phase_at` |
| `src/phase.rs` `test_co2_sublimation_curve` | `tests/phase.tcyr` `test_curves` |
| `src/phase.rs` `test_co2_supercritical` | `tests/phase.tcyr` `test_phase_at` |
| `src/phase.rs` `test_co2_vapor_pressure` | `tests/phase.tcyr` `test_curves` |
| `src/phase.rs` `test_compressed_liquid_above_critical_pressure` | `tests/phase.tcyr` `test_phase_at` |
| `src/phase.rs` `test_heat_for_phase_change_cooling` | `tests/phase.tcyr` `test_latent_heats` |
| `src/phase.rs` `test_heat_for_phase_change_ice_to_steam` | `tests/phase.tcyr` `test_latent_heats` |
| `src/phase.rs` `test_heat_for_phase_change_invalid` | `tests/phase.tcyr` `test_latent_heats` |
| `src/phase.rs` `test_heat_for_phase_change_no_transition` | `tests/phase.tcyr` `test_latent_heats` |
| `src/phase.rs` `test_heat_of_fusion_scales_with_mass` | `tests/phase.tcyr` `test_latent_heats` |
| `src/phase.rs` `test_heat_of_fusion_water` | `tests/phase.tcyr` `test_latent_heats` |
| `src/phase.rs` `test_heat_of_vaporization_water` | `tests/phase.tcyr` `test_latent_heats` |
| `src/phase.rs` `test_mercury_melting_point` | `tests/phase.tcyr` `test_table` |
| `src/phase.rs` `test_nitrogen_boiling_point` | `tests/phase.tcyr` `test_table` |
| `src/phase.rs` `test_nitrogen_gas_at_room_temp` | `tests/phase.tcyr` `test_phase_at` |
| `src/phase.rs` `test_nitrogen_liquid_at_cryogenic` | `tests/phase.tcyr` `test_phase_at` |
| `src/phase.rs` `test_phase_at_invalid_inputs` | `tests/phase.tcyr` `test_phase_at` |
| `src/phase.rs` `test_phase_enum_traits` | `tests/phase.tcyr` `test_table` |
| `src/phase.rs` `test_phase_functions_reject_nan` | `tests/phase.tcyr` `test_latent_heats`, `tests/phase.tcyr` `test_clausius_clapeyron`, `tests/phase.tcyr` `test_phase_at` |
| `src/phase.rs` `test_phase_serde_roundtrip` | `tests/phase.tcyr` `test_json` |
| `src/phase.rs` `test_substance_serde_roundtrip` | `tests/phase.tcyr` `test_json` |
| `src/phase.rs` `test_vapor_pressure_meets_anchor_points` | `tests/phase.tcyr` `test_curves` |
| `src/phase.rs` `test_water_below_triple_point_low_pressure_is_vapor` | `tests/phase.tcyr` `test_phase_at` |
| `src/phase.rs` `test_water_critical_point` | `tests/phase.tcyr` `test_table` |
| `src/phase.rs` `test_water_gas_above_boiling` | `tests/phase.tcyr` `test_phase_at` |
| `src/phase.rs` `test_water_high_pressure_liquid` | `tests/phase.tcyr` `test_phase_at` |
| `src/phase.rs` `test_water_high_temp_low_pressure_is_gas` | `tests/phase.tcyr` `test_phase_at` |
| `src/phase.rs` `test_water_latent_heats` | `tests/phase.tcyr` `test_table` |
| `src/phase.rs` `test_water_liquid_at_room_temp` | `tests/phase.tcyr` `test_phase_at` |
| `src/phase.rs` `test_water_solid_below_freezing` | `tests/phase.tcyr` `test_phase_at` |
| `src/phase.rs` `test_water_supercritical` | `tests/phase.tcyr` `test_phase_at` |
| `src/phase.rs` `test_water_triple_point` | `tests/phase.tcyr` `test_table` |
| `src/phase.rs` `test_water_vapor_pressure_never_exceeds_critical` | `tests/phase.tcyr` `test_curves` |
| `src/state.rs` `test_adiabatic_invalid_gamma` | `tests/state.tcyr` `test_adiabatic` |
| `src/state.rs` `test_adiabatic_roundtrip` | `tests/state.tcyr` `test_adiabatic` |
| `src/state.rs` `test_adiabatic_temperature` | `tests/state.tcyr` `test_adiabatic` |
| `src/state.rs` `test_adiabatic_temperature_validates_t1` | `tests/state.tcyr` `test_adiabatic` |
| `src/state.rs` `test_all_gases_count` | `tests/state.tcyr` `test_gas_data` |
| `src/state.rs` `test_compressibility_factor_validates_inputs` | `tests/state.tcyr` `test_compressibility` |
| `src/state.rs` `test_compressibility_factor_zero_temp` | `tests/state.tcyr` `test_compressibility` |
| `src/state.rs` `test_compressibility_ideal` | `tests/state.tcyr` `test_compressibility` |
| `src/state.rs` `test_compressibility_pitzer_ideal` | `tests/state.tcyr` `test_pitzer` |
| `src/state.rs` `test_compressibility_pitzer_invalid` | `tests/state.tcyr` `test_pitzer` |
| `src/state.rs` `test_compressibility_pitzer_rejects_nonphysical_z` | `tests/state.tcyr` `test_pitzer` |
| `src/state.rs` `test_cubic_eos_validate_critical_constants` | `tests/state.tcyr` `test_redlich_kwong`, `tests/state.tcyr` `test_peng_robinson` |
| `src/state.rs` `test_gas_constant_is_exact` | `tests/state.tcyr` `test_gas_data` |
| `src/state.rs` `test_gas_data_serde_roundtrip` | `tests/state.tcyr` `test_gas_data` |
| `src/state.rs` `test_gas_data_valid` | `tests/state.tcyr` `test_gas_data` |
| `src/state.rs` `test_ideal_gas_pressure` | `tests/state.tcyr` `test_ideal_gas` |
| `src/state.rs` `test_ideal_gas_rejects_negative_moles` | `tests/state.tcyr` `test_ideal_gas`, `tests/state.tcyr` `test_van_der_waals`, `tests/state.tcyr` `test_work` |
| `src/state.rs` `test_ideal_gas_roundtrip` | `tests/state.tcyr` `test_ideal_gas` |
| `src/state.rs` `test_ideal_gas_stp` | `tests/state.tcyr` `test_ideal_gas` |
| `src/state.rs` `test_ideal_gas_temperature_zero_moles` | `tests/state.tcyr` `test_ideal_gas` |
| `src/state.rs` `test_ideal_gas_zero_temp` | `tests/state.tcyr` `test_ideal_gas` |
| `src/state.rs` `test_isobaric_work` | `tests/state.tcyr` `test_work` |
| `src/state.rs` `test_isothermal_work_compression` | `tests/state.tcyr` `test_work` |
| `src/state.rs` `test_isothermal_work_expansion` | `tests/state.tcyr` `test_work` |
| `src/state.rs` `test_isothermal_work_negative_temp` | `tests/state.tcyr` `test_work` |
| `src/state.rs` `test_isothermal_work_zero_temp` | `tests/state.tcyr` `test_work` |
| `src/state.rs` `test_kays_rule_binary` | `tests/state.tcyr` `test_kays` |
| `src/state.rs` `test_kays_rule_invalid` | `tests/state.tcyr` `test_kays` |
| `src/state.rs` `test_kays_rule_pure` | `tests/state.tcyr` `test_kays` |
| `src/state.rs` `test_mixture_rules_reject_bad_fractions` | `tests/state.tcyr` `test_kays` |
| `src/state.rs` `test_negative_temperature` | `tests/state.tcyr` `test_ideal_gas` |
| `src/state.rs` `test_pitzer_second_virial` | `tests/state.tcyr` `test_virial` |
| `src/state.rs` `test_pr_invalid` | `tests/state.tcyr` `test_peng_robinson` |
| `src/state.rs` `test_pr_params` | `tests/state.tcyr` `test_peng_robinson` |
| `src/state.rs` `test_pr_pressure_co2_stp` | `tests/state.tcyr` `test_peng_robinson` |
| `src/state.rs` `test_pr_vs_rk_differ` | `tests/state.tcyr` `test_peng_robinson` |
| `src/state.rs` `test_reduced_coordinates` | `tests/state.tcyr` `test_gas_data` |
| `src/state.rs` `test_rk_invalid` | `tests/state.tcyr` `test_redlich_kwong` |
| `src/state.rs` `test_rk_params` | `tests/state.tcyr` `test_redlich_kwong` |
| `src/state.rs` `test_rk_pressure_co2_stp` | `tests/state.tcyr` `test_redlich_kwong` |
| `src/state.rs` `test_van_der_waals` | `tests/state.tcyr` `test_van_der_waals` |
| `src/state.rs` `test_van_der_waals_negative_temp` | `tests/state.tcyr` `test_van_der_waals` |
| `src/state.rs` `test_van_der_waals_volume_too_small` | `tests/state.tcyr` `test_van_der_waals` |
| `src/state.rs` `test_virial_3rd_order` | `tests/state.tcyr` `test_virial` |
| `src/state.rs` `test_virial_ideal_at_large_vm` | `tests/state.tcyr` `test_virial` |
| `src/state.rs` `test_virial_invalid` | `tests/state.tcyr` `test_virial` |
| `src/state.rs` `test_zero_volume` | `tests/state.tcyr` `test_ideal_gas` |
| `src/steam.rs` `test_critical_point_specific_volume` | `tests/steam.tcyr` `test_lookups` |
| `src/steam.rs` `test_quality_enthalpy_roundtrip` | `tests/steam.tcyr` `test_wet_steam` |
| `src/steam.rs` `test_quality_from_enthalpy` | `tests/steam.tcyr` `test_quality` |
| `src/steam.rs` `test_quality_from_entropy` | `tests/steam.tcyr` `test_quality` |
| `src/steam.rs` `test_quality_midpoint` | `tests/steam.tcyr` `test_quality` |
| `src/steam.rs` `test_quality_of_saturated_vapor_is_one_between_entries` | `tests/steam.tcyr` `test_quality` |
| `src/steam.rs` `test_quality_out_of_range` | `tests/steam.tcyr` `test_quality` |
| `src/steam.rs` `test_quality_x0_is_liquid` | `tests/steam.tcyr` `test_quality` |
| `src/steam.rs` `test_quality_x1_is_vapor` | `tests/steam.tcyr` `test_quality` |
| `src/steam.rs` `test_saturated_by_pressure_1atm` | `tests/steam.tcyr` `test_lookups` |
| `src/steam.rs` `test_saturated_by_pressure_out_of_range` | `tests/steam.tcyr` `test_lookups` |
| `src/steam.rs` `test_saturated_by_pressure_roundtrip` | `tests/steam.tcyr` `test_lookups` |
| `src/steam.rs` `test_saturated_by_pressure_triple` | `tests/steam.tcyr` `test_lookups` |
| `src/steam.rs` `test_saturated_by_temp_boiling_point` | `tests/steam.tcyr` `test_lookups` |
| `src/steam.rs` `test_saturated_by_temp_critical_point` | `tests/steam.tcyr` `test_lookups` |
| `src/steam.rs` `test_saturated_by_temp_interpolation` | `tests/steam.tcyr` `test_lookups` |
| `src/steam.rs` `test_saturated_by_temp_out_of_range` | `tests/steam.tcyr` `test_lookups` |
| `src/steam.rs` `test_saturated_by_temp_triple_point` | `tests/steam.tcyr` `test_lookups` |
| `src/steam.rs` `test_saturated_enthalpy_consistency` | `tests/steam.tcyr` `test_table` |
| `src/steam.rs` `test_saturated_pressure_between_entries_matches_if97` | `tests/steam.tcyr` `test_lookups` |
| `src/steam.rs` `test_saturated_serde_roundtrip` | `tests/steam.tcyr` `test_wet_steam` |
| `src/steam.rs` `test_saturated_table_enthalpy_identity_exact` | `tests/steam.tcyr` `test_table` |
| `src/steam.rs` `test_saturated_table_h_fg_decreasing` | `tests/steam.tcyr` `test_table` |
| `src/steam.rs` `test_saturated_table_monotonic` | `tests/steam.tcyr` `test_table` |
| `src/steam.rs` `test_saturation_line_roundtrip_and_range` | `tests/steam.tcyr` `test_region4` |
| `src/steam.rs` `test_saturation_pressure_if97_verification_values` | `tests/steam.tcyr` `test_region4` |
| `src/steam.rs` `test_saturation_temperature_if97_verification_values` | `tests/steam.tcyr` `test_region4` |
| `src/steam.rs` `test_steam_lookups_reject_nan` | `tests/steam.tcyr` `test_lookups` |
| `src/steam.rs` `test_wet_steam_properties_invalid_quality` | `tests/steam.tcyr` `test_wet_steam` |
| `src/steam.rs` `test_wet_steam_properties_x0` | `tests/steam.tcyr` `test_wet_steam` |
| `src/steam.rs` `test_wet_steam_properties_x1` | `tests/steam.tcyr` `test_wet_steam` |
| `src/steam.rs` `test_wet_steam_serde_roundtrip` | `tests/steam.tcyr` `test_wet_steam` |
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
| `tests/integration.rs` `carnot_engine_efficiency` | `tests/integration.tcyr` `carnot_engine_efficiency` |
| `tests/integration.rs` `clausius_clapeyron_vs_steam_table` | `tests/integration.tcyr` `clausius_clapeyron_vs_steam_table` |
| `tests/integration.rs` `energy_conservation_adiabatic_process` | `tests/integration.tcyr` `energy_conservation_adiabatic_process` |
| `tests/integration.rs` `entropy_increase_irreversible_heat_transfer` | `tests/integration.tcyr` `entropy_increase_irreversible_heat_transfer` |
| `tests/integration.rs` `eos_comparison_co2_stp` | `tests/integration.tcyr` `eos_comparison_co2_stp` |
| `tests/integration.rs` `heat_conduction_through_copper_wall` | `tests/integration.tcyr` `heat_conduction_through_copper_wall` |
| `tests/integration.rs` `ideal_gas_law_roundtrip` | `tests/integration.tcyr` `ideal_gas_law_roundtrip` |
| `tests/integration.rs` `lmtd_vs_ntu_same_heat_exchanger` | `tests/integration.tcyr` `lmtd_vs_ntu_same_heat_exchanger` |
| `tests/integration.rs` `material_thermal_diffusivity_matches_manual` | `tests/integration.tcyr` `material_thermal_diffusivity_matches_manual` |
| `tests/integration.rs` `mixture_z_factor_kays_rule` | `tests/integration.tcyr` `mixture_z_factor_kays_rule` |
| `tests/integration.rs` `nusselt_to_h_to_convection_chain` | `tests/integration.tcyr` `nusselt_to_h_to_convection_chain` |
| `tests/integration.rs` `phase_lookup_agrees_with_steam_table` | `tests/integration.tcyr` `phase_lookup_agrees_with_steam_table` |
| `tests/integration.rs` `quality_consistency_across_methods` | `tests/integration.tcyr` `quality_consistency_across_methods` |
| `tests/integration.rs` `thermal_equilibrium_lumped_capacitance` | `tests/integration.tcyr` `thermal_equilibrium_lumped_capacitance` |
| `tests/integration.rs` `van_der_waals_vs_ideal_gas` | `tests/integration.tcyr` `van_der_waals_vs_ideal_gas` |
| `tests/integration.rs` `view_factor_reciprocity` | `tests/integration.tcyr` `view_factor_reciprocity` |
| `tests/integration.rs` `virial_approaches_ideal_at_low_density` | `tests/integration.tcyr` `virial_approaches_ideal_at_low_density` |
| `tests/integration.rs` `water_boiling_point_consistency` | `tests/integration.tcyr` `water_boiling_point_consistency` |
