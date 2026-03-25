# Ushma

> **Ushma** (Sanskrit: ऊष्मा — heat, thermal energy) — Thermodynamics simulation for AGNOS

Heat transfer, entropy, equations of state, phase transitions, thermodynamic cycles, numerical solvers, chemical and statistical thermodynamics.

Built on [hisab](https://crates.io/crates/hisab) for math foundations and [impetus](https://crates.io/crates/impetus) for body contacts.

**MSRV:** 1.89 · **Edition:** 2024 · **License:** GPL-3.0-only

## Installation

```toml
[dependencies]
ushma = "1"
```

## Features

| Feature | Default | Description |
|---------|---------|-------------|
| `transfer` | yes | Conduction, convection, radiation, fins, heat exchangers, view factors, dimensionless numbers |
| `state` | yes | Ideal gas, van der Waals, Redlich-Kwong, Peng-Robinson, virial EOS, mixture rules |
| `entropy` | yes | Entropy change, Carnot efficiency/COP, Helmholtz/Gibbs free energy, entropy of mixing |
| `material` | yes | 10 built-in materials with thermal properties at ~300 K |
| `phase` | yes | Phase transitions, Clausius-Clapeyron, 8 substances with critical/triple points |
| `cycle` | yes | Otto, Diesel, Brayton, Rankine, refrigeration cycles, T-s/P-v diagrams |
| `numerical` | yes | 1D/2D finite difference (explicit, Crank-Nicolson), thermal networks, adaptive stepping |
| `chem` | yes | Hess's law, Gibbs energy of reaction, equilibrium constants, Van't Hoff, adiabatic flame temperature |
| `stat` | yes | Boltzmann distribution, partition functions, Maxwell-Boltzmann speeds, Einstein/Debye solid Cv |
| `steam` | no | Saturated + superheated steam tables (IAPWS-IF97), quality calculations |
| `ai` | no | Daimon/hoosh agent integration (network dependencies) |
| `logging` | no | Structured logging via `USHMA_LOG` env var |
| `full` | — | Enables all features |

## Architecture

```
ushma (this crate)
  ├── hisab (math: vectors, geometry, calculus)
  └── impetus (physics: rigid body dynamics, collision)

Consumers:
  kiran       ──→ ushma (thermal simulation in game worlds)
  joshua      ──→ ushma (heat engines, thermal experiments, climate models)
```

## Quick Start

```rust
use ushma::transfer;
use ushma::state;
use ushma::entropy;
use ushma::material::COPPER;

// Fourier's law: heat conduction through a copper wall
let q = transfer::conduction(COPPER.conductivity, 0.01, 373.15, 293.15, 0.1).unwrap();

// Ideal gas law at STP
let v = state::ideal_gas_volume(1.0, state::STANDARD_TEMP, state::ATM).unwrap();

// Carnot efficiency
let eta = entropy::carnot_efficiency(500.0, 300.0).unwrap(); // 40%

// Otto cycle
use ushma::cycle;
let otto = cycle::otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0).unwrap();
assert!(otto.efficiency > 0.5);

// Peng-Robinson equation of state
let p = state::peng_robinson_pressure(300.0, 0.02241, 304.13, 7_375_000.0, 0.224).unwrap();

// Adiabatic flame temperature
use ushma::chem;
let t_flame = chem::adiabatic_flame_temperature(
    &[(1.0, &chem::CH4), (2.0, &chem::O2), (7.52, &chem::N2)],
    &[(1.0, &chem::CO2), (2.0, &chem::H2O_GAS), (7.52, &chem::N2)],
    298.15,
).unwrap();
```

## Modules

| Module | Key items | Description |
|--------|-----------|-------------|
| `transfer` | `conduction`, `convection`, `radiation`, `fin_*`, `lmtd_*`, `effectiveness_*`, `reynolds_number`, `nusselt_*`, `view_factor_*` | Heat transfer, fins, exchangers, dimensionless numbers |
| `state` | `ideal_gas_*`, `van_der_waals_*`, `redlich_kwong_*`, `peng_robinson_*`, `virial_*`, `kays_rule_*`, `GasData` | Equations of state, real gas models, mixtures |
| `entropy` | `carnot_efficiency`, `carnot_cop_*`, `helmholtz`, `gibbs`, `entropy_of_mixing` | Thermodynamic potentials, second law |
| `material` | `ThermalMaterial`, `COPPER`...`DIAMOND`, `ALL_MATERIALS` | 10 reference materials at ~300 K |
| `phase` | `SubstanceData`, `Phase`, `clausius_clapeyron_*`, `heat_of_*`, `phase_at` | Phase transitions, 8 substances |
| `steam` | `saturated_by_*`, `superheated_lookup`, `quality_from_*`, `wet_steam_properties` | IAPWS-IF97 steam tables |
| `cycle` | `otto_cycle`, `diesel_cycle`, `brayton_cycle`, `rankine_cycle`, `refrigeration_cycle`, `cycle_*_diagram` | Thermodynamic cycles, diagrams |
| `numerical` | `ThermalGrid1D`, `ThermalGrid2D`, `ThermalNetwork`, `BoundaryCondition` | FDM solvers, thermal networks |
| `chem` | `reaction_enthalpy`, `reaction_gibbs`, `equilibrium_constant`, `vant_hoff_k`, `adiabatic_flame_temperature` | Chemical thermodynamics |
| `stat` | `boltzmann_*`, `partition_*`, `maxwell_boltzmann_*`, `einstein_cv`, `debye_cv`, `equipartition_*` | Statistical mechanics |

## Examples

```sh
cargo run --example basic_thermo --features transfer,state,material
cargo run --example heat_engine --features cycle,state,entropy
cargo run --example building_thermal --features transfer,material,numerical
cargo run --example cooling_system --features transfer,material
cargo run --example phase_change --features phase,transfer
```

## Building

```sh
cargo build --all-features
cargo test --all-features
cargo bench --all-features
```

## Documentation

- [Roadmap](docs/development/roadmap.md)
- [Contributing](CONTRIBUTING.md)
- [Security Policy](SECURITY.md)
- [Changelog](CHANGELOG.md)

## License

GPL-3.0-only — see [LICENSE](LICENSE).
