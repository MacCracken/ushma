# Ushma

> **Ushma** (Sanskrit: ऊष्मा — heat, thermal energy) — Thermodynamics simulation for AGNOS

Heat transfer, entropy, equations of state, and thermal material properties. Built on [hisab](https://crates.io/crates/hisab) for math foundations and [impetus](https://crates.io/crates/impetus) for body contacts.

## Features

| Feature | Default | Description |
|---------|---------|-------------|
| `transfer` | yes | Conduction (Fourier), convection (Newton), radiation (Stefan-Boltzmann), thermal resistance, lumped capacitance |
| `state` | yes | Ideal gas law, van der Waals, isothermal/isobaric/adiabatic processes, compressibility |
| `entropy` | yes | Entropy change, Carnot efficiency/COP, Helmholtz/Gibbs free energy, entropy of mixing |
| `material` | yes | 10 built-in materials (copper, aluminum, iron, steel, glass, water, air, wood, concrete, diamond) |
| `ai` | no | Daimon/hoosh integration (network deps) |
| `logging` | no | Structured logging via `USHMA_LOG` env var |
| `full` | — | Enables all features |

## Architecture

```
ushma (this crate)
  ├── hisab (math: vectors, geometry, calculus, numerical methods)
  └── impetus (physics: rigid body dynamics, collision detection)

Consumers:
  kiran       ──→ ushma (thermal simulation in game worlds)
  joshua      ──→ ushma (simulation: heat engines, thermal experiments)
```

## Quick Start

```rust
use ushma::transfer::{conduction, convection, radiation, lumped_capacitance};
use ushma::state::{ideal_gas_pressure, ideal_gas_volume, van_der_waals_pressure, ATM, STANDARD_TEMP};
use ushma::entropy::{carnot_efficiency, helmholtz, gibbs};
use ushma::material::{COPPER, WATER, ALL_MATERIALS};

// Fourier's law: heat conduction through a copper wall
let q = conduction(COPPER.conductivity, 0.01, 373.15, 293.15, 0.1).unwrap();

// Ideal gas law: volume of 1 mol at STP
let v = ideal_gas_volume(1.0, STANDARD_TEMP, ATM).unwrap(); // ~0.02241 m³

// Carnot efficiency: steam engine at 500K/300K
let eta = carnot_efficiency(500.0, 300.0).unwrap(); // 0.4 (40%)

// Material properties
let alpha = COPPER.diffusivity(); // thermal diffusivity (m²/s)
```

## Modules

| Module | Functions | Description |
|--------|-----------|-------------|
| `transfer` | `conduction`, `convection`, `radiation`, `thermal_resistance_*`, `heat_stored`, `thermal_diffusivity`, `biot_number`, `lumped_capacitance` | SI units (W, m, K, s) |
| `state` | `ideal_gas_*`, `van_der_waals_pressure`, `isothermal_work`, `isobaric_work`, `adiabatic_temperature`, `compressibility_factor` | Gas constants (R, ATM, standard temp) |
| `entropy` | `ideal_gas_entropy_change`, `isothermal_entropy_change`, `heat_transfer_entropy`, `carnot_efficiency`, `carnot_cop_refrigeration`, `helmholtz`, `gibbs`, `entropy_of_mixing` | Second law, thermodynamic potentials |
| `material` | `ThermalMaterial`, `COPPER`, `ALUMINUM`, `IRON`, `STEEL`, `GLASS`, `WATER`, `AIR`, `WOOD_OAK`, `CONCRETE`, `DIAMOND` | Reference data at ~300 K |

## Building

```sh
cargo build --all-features
cargo test --all-features
```

## Roadmap

See [docs/development/roadmap.md](docs/development/roadmap.md).

## License

GPL-3.0 — see [LICENSE](LICENSE).
