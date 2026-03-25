# Changelog

## [0.1.0] - 2026-03-24

### Added

#### transfer — Heat Transfer
- Fourier's law: conduction through solid materials (q = kAΔT/L)
- Newton's law of cooling: convective heat transfer (q = hAΔT)
- Stefan-Boltzmann law: radiative heat transfer (q = εσA(T⁴ - T_surr⁴))
- Thermal resistance: conduction (L/kA), convection (1/hA), series, parallel
- Heat stored in a body (Q = mcΔT)
- Thermal diffusivity (α = k/ρc_p)
- Biot number for lumped capacitance validity
- Lumped capacitance transient cooling model
- Physical constants: Stefan-Boltzmann σ, Boltzmann k_B

#### state — Equations of State
- Ideal gas law: pressure, volume, temperature (PV = nRT)
- Van der Waals equation for real gases
- Isothermal work (W = nRT⋅ln(V₂/V₁))
- Isobaric work (W = PΔV)
- Adiabatic temperature change (PV^γ = const)
- Compressibility factor (Z = PV/nRT)
- Constants: R, ATM, standard temperature

#### entropy — Thermodynamic Potentials
- Ideal gas entropy change (ΔS = nCv⋅ln(T₂/T₁) + nR⋅ln(V₂/V₁))
- Isothermal entropy change
- Heat transfer entropy (ΔS = Q/T)
- Carnot efficiency (η = 1 - T_cold/T_hot)
- Carnot COP for refrigeration
- Helmholtz free energy (A = U - TS)
- Gibbs free energy (G = H - TS)
- Clausius inequality check (spontaneity)
- Entropy of mixing for ideal gases

#### material — Thermal Properties
- ThermalMaterial struct with conductivity, specific heat, density, phase transitions
- 10 built-in materials: copper, aluminum, iron, steel, glass, water, air, wood (oak), concrete, diamond
- Thermal diffusivity and volumetric heat capacity methods
- Serde serialization support

#### Infrastructure
- Error type: UshmaError with #[non_exhaustive], 8 variants
- AI integration: daimon/hoosh client (feature-gated)
- Structured logging via USHMA_LOG (feature-gated)
- CI/CD: GitHub Actions workflows (ci.yml, release.yml)
- Criterion benchmarks: 23 benchmarked functions
- Integration tests: 8 cross-module tests
- SECURITY.md, CONTRIBUTING.md, CODE_OF_CONDUCT.md
- deny.toml, codecov.yml, Makefile
- scripts/bench-history.sh, scripts/version-bump.sh
