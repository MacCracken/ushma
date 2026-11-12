# Architecture Overview

## Module map

| Module | Feature | Depends on | Role |
|--------|---------|-----------|------|
| `error` | always | — | `UshmaError` (`#[non_exhaustive]`), `Result<T>` |
| `validate` | always (crate-private) | — | `finite_gt` / `finite_ge`: positivity/range checks that also reject NaN and ±∞; `invalid_param`: `#[cold]` out-of-line error construction so validated hot paths stay inlinable |
| `bridge` | always | — | Primitive-value conversions for bijli, kimiya, dravya, badal, pravash, sharira (infallible, return 0 on invalid input) |
| `transfer` | `transfer` | — | Conduction, convection, radiation, resistances, dimensionless numbers, fins, heat exchangers, view factors |
| `state` | `state` | — | Ideal gas, vdW, Redlich-Kwong, Peng-Robinson, virial, Pitzer, Kay's rule; `GAS_CONSTANT`, `ATM` |
| `entropy` | `entropy` | `state` | Entropy changes, Carnot limits, free energies, mixing |
| `material` | `material` | — | `ThermalMaterial` reference data (~300 K) |
| `phase` | `phase` | `state` | `SubstanceData`, vapor-pressure / sublimation curves, `phase_at`, latent heat budgets |
| `steam` | `steam` | `phase` | IAPWS-IF97 Region 4 saturation line; saturated and superheated tables; quality |
| `cycle` | `cycle` | `state`, `entropy` (+ `steam` for Rankine/refrigeration) | Otto, Diesel, Brayton, Rankine, vapor-compression; T-s / P-v diagram data |
| `numerical` | `numerical` | `transfer` | 1D explicit / Crank-Nicolson / adaptive FD, 2D Gauss-Seidel, thermal resistance networks |
| `chem` | `chem` | `state`, `entropy` | Hess's law, ΔG°, K and ln K, Van't Hoff, adiabatic flame temperature |
| `stat` | `stat` | `state`, `transfer` | Boltzmann, partition functions, Maxwell-Boltzmann, Einstein/Debye Cv |
| `integration::soorat` | `soorat-compat` (+ `numerical` / `cycle`) | — | Render-ready data for soorat: heatmaps, profiles, cycle diagrams, network graphs, flux fields |
| `ai` | `ai` | reqwest, serde_json | daimon agent registration client |
| `logging` | `logging` | tracing-subscriber | `USHMA_LOG` env-filter subscriber |

## Data flow

```
material / phase / steam / state   (reference data, equations of state)
        │
        ▼
transfer ──► numerical            (rates → discretized fields and networks)
        │
        ▼
entropy ──► cycle ──► integration::soorat   (cycle analysis → plot data)
chem, stat                        (independent: reaction and molecular thermodynamics)
```

All quantities are SI (K, Pa, m³, J, W, mol, kg). Fallible functions return
`Result<T, UshmaError>`; pure closed-form helpers return `f64` and are
`#[must_use]`. Input guards reject NaN/±∞ as well as out-of-range values.

## Dependency stack

- **hisab** (AGNOS math): `calc::lerp`, `calc::integral_simpson`, `num::bisection`, `num::gaussian_elimination`
- **serde** (derive), **thiserror**, **tracing**
- Optional: **reqwest** + **serde_json** (`ai`), **tracing-subscriber** (`logging`)

## Consumers

- **kiran** — thermal simulation in game worlds
- **joshua** — heat engines, thermal experiments, climate models
- **badal**, **dravya**, **bijli**, **kimiya**, **pravash**, **sharira** — via `bridge`
- **soorat** — via `integration::soorat`
