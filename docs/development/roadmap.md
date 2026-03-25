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

## V0.2 — Phase Transitions

- [ ] Clausius-Clapeyron equation (phase boundary slope dP/dT)
- [ ] Latent heat of fusion and vaporization
- [ ] Phase diagram lookup (solid/liquid/gas regions)
- [ ] Triple point and critical point data for common substances
- [ ] Saturated steam tables (temperature and pressure indexed)
- [ ] Quality (dryness fraction) for wet steam
- [ ] Superheated steam properties

## V0.3 — Heat Engines & Cycles

- [ ] Otto cycle (spark ignition engine)
- [ ] Diesel cycle (compression ignition)
- [ ] Rankine cycle (steam power plant)
- [ ] Brayton cycle (gas turbine)
- [ ] Refrigeration cycle (vapor compression)
- [ ] Heat pump COP
- [ ] Cycle T-s and P-v diagram data generation
- [ ] Thermal efficiency comparison across cycles

## V0.4 — Extended Heat Transfer

- [ ] Fin heat transfer (rectangular, triangular, pin fins)
- [ ] Fin efficiency and effectiveness
- [ ] Heat exchangers: LMTD method (parallel/counter flow)
- [ ] Heat exchangers: ε-NTU method
- [ ] Nusselt number correlations (forced/natural convection)
- [ ] Reynolds and Prandtl number calculations
- [ ] Thermal boundary layer thickness
- [ ] View factors for radiation between surfaces

## V0.5 — Numerical Methods

- [ ] 1D transient conduction (explicit finite difference)
- [ ] 1D transient conduction (implicit/Crank-Nicolson)
- [ ] 2D steady-state conduction (Gauss-Seidel iteration)
- [ ] Thermal network solver (resistance matrix)
- [ ] Integration with hisab ODE solvers for transient problems
- [ ] Adaptive time stepping for stability

## V0.6 — Real Gas Models

- [ ] Redlich-Kwong equation of state
- [ ] Peng-Robinson equation of state
- [ ] Virial equation of state (B, C coefficients)
- [ ] Generalized compressibility charts (reduced coordinates)
- [ ] Acentric factor database
- [ ] Mixture rules (Kay's rule, mixing parameters)

## V0.7 — Chemical Thermodynamics

- [ ] Standard enthalpy of formation
- [ ] Hess's law (reaction enthalpy from formation enthalpies)
- [ ] Gibbs free energy of reaction
- [ ] Equilibrium constant from ΔG
- [ ] Van't Hoff equation (K vs temperature)
- [ ] Adiabatic flame temperature

## V0.8 — Statistical Thermodynamics

- [ ] Boltzmann distribution
- [ ] Partition functions (translational, rotational, vibrational)
- [ ] Maxwell-Boltzmann speed distribution
- [ ] Einstein and Debye models for solid heat capacity
- [ ] Equipartition theorem
- [ ] Entropy from microstates (S = k⋅ln(W))

## V0.9 — Examples & Documentation

- [ ] Example: heat engine cycle simulator
- [ ] Example: building thermal model (walls, windows, insulation)
- [ ] Example: cooling system design
- [ ] Example: phase change material storage
- [ ] Physics explanations in module-level documentation
- [ ] Worked examples in doc comments

## V1.0 — Stable Release

- [ ] API review: naming consistency, parameter ordering
- [ ] `#[must_use]` on all pure functions and accessors
- [ ] Feature gate audit (minimize compile time for partial usage)
- [ ] Documentation coverage check (all public items documented)
- [ ] README update with full feature matrix
- [ ] Performance audit and optimization pass

## Consumers

| Consumer | What it uses |
|----------|-------------|
| **kiran** | Thermal simulation in game worlds (heat zones, fire, cooling) |
| **joshua** | Simulation: heat engines, thermal experiments, climate models |
| **impetus** | Body contact thermal transfer (friction heating) |

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
