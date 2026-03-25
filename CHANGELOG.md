# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-03-25

### Added

#### transfer — Heat Transfer
- Fourier's law: conduction through solid materials
- Newton's law of cooling: convective heat transfer
- Stefan-Boltzmann law: radiative heat transfer
- Thermal resistance: conduction, convection, series, parallel networks
- Heat stored (Q = mcΔT), thermal diffusivity, Biot number
- Lumped capacitance transient cooling model
- Reynolds, Prandtl, Nusselt number calculations
- Dittus-Boelter (turbulent pipe) and Churchill-Chu (natural convection) correlations
- Laminar and thermal boundary layer thickness (Blasius)
- Fin heat transfer: rectangular fins, efficiency, effectiveness
- Heat exchangers: LMTD method (parallel/counter flow)
- Heat exchangers: ε-NTU method with special cases (Cr=0, Cr=1)
- Radiation view factors: parallel plates, perpendicular plates, coaxial disks

#### state — Equations of State
- Ideal gas law: pressure, volume, temperature (PV = nRT)
- Van der Waals equation for real gases
- Redlich-Kwong equation of state
- Peng-Robinson equation of state (with acentric factor)
- Virial equation of state (2nd and 3rd order)
- Pitzer correlation for second virial coefficient
- Generalized compressibility (Pitzer Z⁰ + ωZ¹)
- Isothermal, isobaric, and adiabatic work
- Compressibility factor
- GasData struct with 8 built-in gases (N₂, O₂, CO₂, CH₄, H₂O, NH₃, C₂H₆, C₃H₈)
- Kay's rule mixture rules for Tc, Pc, ω

#### entropy — Thermodynamic Potentials
- Ideal gas entropy change, isothermal entropy change
- Heat transfer entropy (ΔS = Q/T)
- Carnot efficiency and COP for refrigeration
- Helmholtz and Gibbs free energy
- Clausius inequality check (spontaneity)
- Entropy of mixing for ideal gases

#### material — Thermal Properties
- ThermalMaterial struct with Cow<'static, str> name
- 10 built-in materials: copper, aluminum, iron, steel, glass, water, air, wood (oak), concrete, diamond
- Thermal diffusivity and volumetric heat capacity methods
- Serde serialization/deserialization

#### phase — Phase Transitions
- SubstanceData struct with triple/critical points, latent heats, molar mass
- 8 built-in substances: water, nitrogen, oxygen, CO₂, ammonia, methane, ethanol, mercury
- Clausius-Clapeyron equation (exact slope and integrated pressure form)
- Phase diagram lookup (solid/liquid/gas/supercritical)
- Latent heat helpers: fusion, vaporization, multi-phase heat budget

#### steam — Steam Tables
- 41-entry saturated steam table (273.16 K to 647.096 K, IAPWS-IF97)
- Saturated lookup by temperature and pressure (binary search + hisab lerp)
- 72-entry superheated steam table (9 pressure tiers, bilinear interpolation)
- Quality calculations from volume, enthalpy, entropy
- Wet steam mixture properties

#### cycle — Thermodynamic Cycles
- Otto cycle (spark ignition engine)
- Diesel cycle (compression ignition)
- Brayton cycle (gas turbine) with back-work ratio
- Rankine cycle (steam power plant, saturated and superheated)
- Vapor-compression refrigeration cycle (with hisab bisection for isentropic state)
- Heat pump COP
- T-s and P-v diagram data generation
- Cycle comparison with second-law efficiency

#### numerical — Numerical Methods
- 1D transient conduction: explicit (forward Euler) finite difference
- 1D transient conduction: implicit (Crank-Nicolson) via Thomas algorithm
- Adaptive time stepping for explicit stability
- 2D steady-state conduction (Gauss-Seidel iteration)
- Thermal resistance network solver (via hisab gaussian elimination)
- Boundary conditions: fixed, insulated, convective

#### chem — Chemical Thermodynamics
- Species struct with formation enthalpy, Gibbs energy, heat capacity
- 12 built-in species: H₂, O₂, N₂, H₂O(g), H₂O(l), CO₂, CO, CH₄, C₂H₆, C₃H₈, NH₃, NO
- Hess's law: reaction enthalpy from formation enthalpies
- Gibbs energy of reaction
- Equilibrium constant from ΔG
- Van't Hoff equation (K vs temperature)
- Adiabatic flame temperature (constant-Cp model)

#### stat — Statistical Thermodynamics
- Boltzmann distribution: probability, population ratio
- Entropy from microstates (S = k·ln(W))
- Equipartition theorem: energy per molecule, molar Cv
- Partition functions: translational, rotational, vibrational
- Maxwell-Boltzmann speed distribution: PDF, most probable, mean, RMS speeds
- Einstein model for solid heat capacity
- Debye model for solid heat capacity (via hisab Simpson integration)
- Physical constants: Planck h, Avogadro Nₐ

#### Infrastructure
- Error type: UshmaError with #[non_exhaustive], 16 variants
- 14 feature flags (9 default + steam, ai, logging, full)
- 362 tests (332 unit + 30 integration)
- 57 Criterion benchmarks with CSV history tracking
- 5 examples: basic_thermo, heat_engine, building_thermal, cooling_system, phase_change
- hisab integration: Simpson integration, lerp, bisection, gaussian elimination
- AI integration: daimon/hoosh client (feature-gated)
- Structured logging via USHMA_LOG (feature-gated)
- CI/CD: GitHub Actions (format, clippy, audit, deny, test, coverage, semver, doc)
- Supply chain: cargo-deny, cargo-audit
- Scripts: bench-history.sh, version-bump.sh
