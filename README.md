# Ushma

> **Ushma** (Sanskrit: ऊष्मा — heat, thermal energy) — Thermodynamics simulation for AGNOS

Heat transfer, entropy, equations of state, phase transitions, steam tables, thermodynamic
cycles, numerical solvers, chemical and statistical thermodynamics, cross-crate bridges,
visualization data for soorat, and a daimon/hoosh client.

**Language:** [Cyrius](https://github.com/MacCracken/cyrius) (toolchain pinned in `cyrius.cyml`) · **License:** GPL-3.0-only

> 🚧 **Porting to Cyrius.** The Rust crate's final release, **1.4.0**, is frozen in
> [`rust-old/`](rust-old/) as the parity oracle. The Cyrius library ships as **2.0.0**, with
> every module, feature, test, benchmark and example ported. Plan of record:
> [docs/development/ushma-port-plan.md](docs/development/ushma-port-plan.md). Live progress:
> [docs/development/state.md](docs/development/state.md).
>
> Rust users: `ushma = "1"` on crates.io stays at 1.4.0. Its README is
> [`rust-old/README.md`](rust-old/README.md).

## Bundles

| Bundle | Contents |
|--------|----------|
| `dist/ushma.cyr` (`[lib]`) | every module: transfer, state, entropy, material, phase, steam, cycle, numerical, chem, stat, bridge, soorat data, logging (sakshi) |
| `dist/ushma-ai.cyr` (`[lib.ai]`) | the same, plus the daimon/hoosh client on sandhi (pulls the TLS stack) |

A consumer includes exactly one bundle. Each Rust cargo feature maps onto one of them (port plan §3.2).

## Modules

| Module | Covers |
|--------|--------|
| transfer | conduction, convection, radiation, resistances, fins, LMTD and ε-NTU exchangers, Re/Pr/Nu correlations, boundary layers, view factors |
| state | ideal gas, Van der Waals, Redlich–Kwong, Peng–Robinson, virial, Pitzer correlation, Kay's rule, gas data |
| entropy | entropy changes, Carnot efficiency and COP, Helmholtz and Gibbs energies, Clausius inequality, mixing |
| material | 10 reference materials at ~300 K |
| phase | phase diagrams for 8 substances, vapor and sublimation curves, latent heats, Clausius–Clapeyron |
| steam | IAPWS-IF97 saturated and superheated tables, exact Region 4 saturation line, quality, isentropic lookup |
| cycle | Otto, Diesel, Brayton, Rankine, refrigeration, heat pump, cycle comparison, T-s and P-v diagrams |
| numerical | 1D explicit / Crank–Nicolson / adaptive conduction, 2D Gauss–Seidel steady state, thermal networks |
| chem | Hess's law, reaction Gibbs energy, equilibrium constants, van 't Hoff, adiabatic flame temperature |
| stat | Boltzmann, partition functions, Maxwell–Boltzmann speeds, equipartition, Einstein and Debye heat capacities |
| bridge | conversions for bijli, kimiya, dravya, badal, pravash and sharira |
| soorat | heatmap, temperature-profile, cycle-diagram, network-graph and heat-flux data for soorat |
| logging | sakshi spans and levels, `USHMA_LOG` |
| ai | daimon agent registration and hoosh configuration (`[lib.ai]`) |

## Dependency stack

```
ushma
  ├── Cyrius stdlib: math, ganita (pow, tanh, cbrt, Gaussian elimination), sakshi + log,
  │                  bayan (JSON), sandhi (HTTP client, [lib.ai] only), io, result, vec, alloc
  └── hisab 3.4.1 (interpolation, Simpson integration)

Consumers:
  kiran   ──→ ushma (thermal simulation in game worlds)
  joshua  ──→ ushma (heat engines, thermal experiments, climate models)
  soorat  ──→ ushma (thermal visualization)
```

## Building

```sh
cyrius deps                              # vendor the stdlib leaves and hisab into lib/
cyrius build src/main.cyr build/ushma    # compile
cyrius test                              # parity + reference tests
cyrius bench                             # benchmarks
```

## Documentation

- [Port plan](docs/development/ushma-port-plan.md) and [state](docs/development/state.md)
- [Roadmap](docs/development/roadmap.md)
- [Architecture decisions](docs/adr/)
- [Contributing](CONTRIBUTING.md)
- [Security Policy](SECURITY.md)
- [Changelog](CHANGELOG.md)

## License

GPL-3.0-only — see [LICENSE](LICENSE).
