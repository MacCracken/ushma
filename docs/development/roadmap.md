# Ushma Roadmap

> **Ushma** is the thermodynamics simulation crate. Math foundations come from [hisab](https://github.com/MacCracken/hisab). Body contacts come from [impetus](https://github.com/MacCracken/impetus).

## Scope

Ushma owns the **physics of heat**: how thermal energy transfers, transforms, and equilibrates. It provides the math; consumers decide what to do with it (simulate engines, model climate, thermal effects in games).

Ushma does NOT own:
- **Physics engine** → impetus (rigid body dynamics, collision)
- **Math primitives** → hisab (vectors, geometry, calculus, ODE/PDE solvers)
- **Rendering** → soorat/kiran (they consume ushma for thermal visualization)
- **Fluid dynamics** → future crate (CFD, Navier-Stokes)

## V0.1 — Foundation (done)

### transfer
- [x] Fourier's law (conduction through solids)
- [x] Newton's law of cooling (convection)
- [x] Stefan-Boltzmann law (radiation)
- [x] Thermal resistance: conduction, convection, series, parallel
- [x] Heat stored (Q = mcΔT)
- [x] Thermal diffusivity (α = k/ρc_p)
- [x] Biot number
- [x] Lumped capacitance transient model

### state
- [x] Ideal gas law (P, V, T)
- [x] Van der Waals equation
- [x] Isothermal work
- [x] Isobaric work
- [x] Adiabatic temperature change
- [x] Compressibility factor

### entropy
- [x] Ideal gas entropy change
- [x] Isothermal entropy change
- [x] Heat transfer entropy (ΔS = Q/T)
- [x] Carnot efficiency
- [x] Carnot COP (refrigeration)
- [x] Helmholtz free energy
- [x] Gibbs free energy
- [x] Clausius inequality (spontaneity check)
- [x] Entropy of mixing

### material
- [x] ThermalMaterial struct (k, c_p, ρ, melting/boiling points)
- [x] 10 built-in materials (copper, aluminum, iron, steel, glass, water, air, wood, concrete, diamond)
- [x] Diffusivity and volumetric heat capacity methods

## V0.2 — Phase Transitions (done)

- [x] Clausius-Clapeyron equation (phase boundary slope dP/dT)
- [x] Latent heat of fusion and vaporization
- [x] Phase diagram lookup (solid/liquid/gas regions)
- [x] Triple point and critical point data for common substances
- [x] Saturated steam tables (temperature and pressure indexed)
- [x] Quality (dryness fraction) for wet steam
- [x] Superheated steam properties

## V0.3 — Heat Engines & Cycles (done)

- [x] Otto cycle (spark ignition engine)
- [x] Diesel cycle (compression ignition)
- [x] Rankine cycle (steam power plant)
- [x] Brayton cycle (gas turbine)
- [x] Refrigeration cycle (vapor compression)
- [x] Heat pump COP
- [x] Cycle T-s and P-v diagram data generation
- [x] Thermal efficiency comparison across cycles

## V0.4 — Extended Heat Transfer (done)

- [x] Fin heat transfer (rectangular, insulated tip)
- [x] Fin efficiency and effectiveness
- [x] Heat exchangers: LMTD method (parallel/counter flow)
- [x] Heat exchangers: ε-NTU method
- [x] Nusselt number correlations (Dittus-Boelter, Churchill-Chu)
- [x] Reynolds and Prandtl number calculations
- [x] Thermal boundary layer thickness (Blasius)
- [x] View factors for radiation (parallel plates, perpendicular plates, coaxial disks)

## V0.5 — Numerical Methods (done)

- [x] 1D transient conduction (explicit finite difference)
- [x] 1D transient conduction (implicit/Crank-Nicolson)
- [x] 2D steady-state conduction (Gauss-Seidel iteration)
- [x] Thermal network solver (conductance matrix via hisab gaussian elimination)
- [x] Integration with hisab (Simpson integration, lerp, bisection, gaussian elimination)
- [x] Adaptive time stepping for stability

## V0.6 — Real Gas Models (done)

- [x] Redlich-Kwong equation of state
- [x] Peng-Robinson equation of state
- [x] Virial equation of state (B, C coefficients)
- [x] Generalized compressibility charts (Pitzer correlation)
- [x] Acentric factor database (8 gases)
- [x] Mixture rules (Kay's rule for Tc, Pc, ω)

## V0.7 — Chemical Thermodynamics (done)

- [x] Standard enthalpy of formation (12 species)
- [x] Hess's law (reaction enthalpy from formation enthalpies)
- [x] Gibbs free energy of reaction
- [x] Equilibrium constant from ΔG
- [x] Van't Hoff equation (K vs temperature)
- [x] Adiabatic flame temperature

## V0.8 — Statistical Thermodynamics (done)

- [x] Boltzmann distribution
- [x] Partition functions (translational, rotational, vibrational)
- [x] Maxwell-Boltzmann speed distribution
- [x] Einstein and Debye models for solid heat capacity
- [x] Equipartition theorem
- [x] Entropy from microstates (S = k⋅ln(W))

## V0.9 — Examples & Documentation (done)

- [x] Example: heat engine cycle simulator
- [x] Example: building thermal model (walls, windows, insulation)
- [x] Example: cooling system design
- [x] Example: phase change material storage
- [x] Module-level documentation with feature tables
- [x] README with full module/feature matrix

## V1.0 — Stable Release (done)

- [x] API review: naming consistency, parameter ordering
- [x] `#[must_use]` on all pure functions and accessors
- [x] Feature gate audit (all 10 features compile independently)
- [x] Documentation coverage check (rustdoc -D warnings clean)
- [x] README update with full feature matrix
- [x] hisab integration (Simpson, lerp, bisection, gaussian elimination)
- [x] Performance: 57 benchmarked functions with CSV history

## V1.1 — P(-1) Hardening (done)

- [x] Remove unused `impetus` dependency
- [x] Fix `unwrap()` in library code (numerical.rs, ai.rs)
- [x] Standardize temperature validation (T=0 handling)
- [x] Add input validation: `ideal_gas_temperature`, `heat_for_phase_change`
- [x] Fix `entropy_of_mixing` to accept zero mole fractions
- [x] Add `#[tracing::instrument]` to 23 core functions
- [x] Baseline benchmarks (57 functions, 14 runs in history)

## Cross-Crate Bridges

- [ ] **`bridge.rs` module** — primitive-value conversions for cross-crate thermodynamics
- [ ] **bijli bridge**: Joule heating power (W) → temperature rise rate (K/s); EM absorption → volumetric heat source
- [ ] **kimiya bridge**: reaction enthalpy (J/mol) → heat release rate; equilibrium constant → temperature-dependent yield
- [ ] **dravya bridge**: temperature (K) → thermal expansion strain; thermal gradient (K/m) → thermal stress (Pa)
- [ ] **badal bridge**: altitude (m), lapse rate (K/m) → atmospheric temperature; humidity ratio → wet-bulb temperature
- [ ] **pravash bridge**: fluid velocity → convective heat transfer coefficient; turbulent kinetic energy → eddy thermal diffusivity

---

## Soorat Integration (rendering visualization)

- [ ] **`integration/soorat.rs` module** — feature-gated `soorat-compat`, visualization data structures for soorat to consume
- [ ] **Thermal grid heatmap**: Expose `ThermalGridVisualization` from `ThermalGrid2D` (temperature values, dimensions, spacing, min/max bounds) for 2D heatmap rendering
- [ ] **1D temperature profile**: Expose `TemperatureProfile` from `ThermalGrid1D` (node temperatures, spacing, boundary conditions) for line/ribbon rendering
- [ ] **Cycle diagrams**: Expose `CycleDiagramData` wrapping `cycle_ts_diagram()` / `cycle_pv_diagram()` output (Vec<DiagramPoint>, CycleKind, state points) for T-s and P-v plot rendering
- [ ] **Thermal network graph**: Expose `ThermalNetworkVisualization` from `ThermalNetwork` (node positions/temperatures, resistance edges) for node-link rendering
- [ ] **Heat flux vectors**: Expose `HeatFluxField` — grid of flux vectors from `ThermalGrid2D` gradients for arrow/streamline rendering

---

## Consumers

| Consumer | What it uses |
|----------|-------------|
| **kiran** | Thermal simulation in game worlds (heat zones, fire, cooling) |
| **joshua** | Simulation: heat engines, thermal experiments, climate models |
| **badal** | Weather/atmospheric modeling (thermo feature) |
| **dravya** | Material/substance modeling |

## Boundary with Other Crates

| Feature | ushma | other |
|---------|-------|-------|
| Heat transfer math | Yes | — |
| Rigid body dynamics | — | impetus |
| ODE/PDE solvers | — | hisab |
| Fluid dynamics (CFD) | — | future crate |
| Vector/matrix math | — | hisab |
| Thermal visualization | — | kiran/soorat |
| Material science (stress/strain) | — | impetus |
