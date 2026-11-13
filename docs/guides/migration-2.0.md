# Migrating from ushma 1.4 (Rust) to ushma 2.0 (Cyrius)

ushma 2.0.0 is the Rust crate ported to Cyrius. Every module, function, type, test, benchmark and example of 1.4.0 has
its counterpart, and the port is checked against 1.4.0 itself: 154,448 golden rows replayed on every target, the
Rust tests ported one for one, and the five examples printing 1.4.0's output byte for byte. The Rust crate stays at
1.4.0, frozen at the tag `1.4.0` (and in `rust-old/` until that is retired); it gets no fixes and no releases.

This guide is for code that used the Rust crate. It covers getting the library, the naming rule, how calls and
results look, and every place where the Cyrius API differs from Rust's on purpose, each with its record in
[ADR 0002](../adr/0002-port-api-divergences.md).

## 1. Getting the library

ushma ships as two bundles, built by `cyrius distlib --all`:

| Bundle | Profile | Holds | Rust features |
|---|---|---|---|
| `dist/ushma.cyr` | `[lib]` | every module but the ai client | `transfer`, `state`, `entropy`, `material`, `phase`, `steam`, `cycle`, `numerical`, `chem`, `stat`, `soorat-compat`, `logging` |
| `dist/ushma-ai.cyr` | `[lib.ai]` | the same, plus the daimon/hoosh client on sandhi (and its TLS stack) | the above and `ai` (`full`) |

Rust's features were compile-time gates. In Cyrius every physics module is in `[lib]`, since DCE drops what a program
doesn't call; only the ai client, which pulls the HTTP and TLS stack, is a separate profile. Depend on one bundle:

```toml
[deps.ushma]
git = "https://github.com/MacCracken/ushma.git"
tag = "2.0.0"
modules = ["dist/ushma.cyr"]      # or "dist/ushma-ai.cyr" for the daimon client
```

and `include "dist/ushma.cyr"` in the program (the examples in [`../examples/`](../examples/) do exactly this).
ushma brings hisab with it; its stdlib leaves are in the bundle's sidecar (`dist/ushma.deps`).

## 2. Names

Cyrius has one namespace, so every public name carries ushma's prefix (ADR 0001, the P6 amendment):

| Rust | Cyrius |
|---|---|
| `transfer::conduction(..)` (a free function) | `us_conduction(..)` |
| `ThermalGrid1D::step_explicit(&mut self, ..)` (a method) | `us_grid1d_step_explicit(g, ..)`, the value first |
| `ThermalGrid1D` (a type) | `UsGrid1D` |
| `material::COPPER` (a preset) | `US_MATERIAL_COPPER`, an index: `us_material_preset(US_MATERIAL_COPPER, &m)` |
| `state::GAS_CONSTANT` (a constant) | `US_GAS_CONSTANT` |
| `Phase::Liquid` (an enum variant) | `US_PHASE_LIQUID`, a member of `enum UsPhase` |
| `UshmaError::InvalidTemperature { .. }` | `US_ERR_INVALID_TEMPERATURE`, a member of `enum UsError` |
| `serde_json::to_string(&x)` / `from_str(..)` | `us_<t>_to_json(sb, x)` / `us_<t>_from_json(buf, n, out)` (`_from_json_a` on an allocator) |
| `impl Default for T` | `us_<t>_default(out)` |

The types, and the prefix of their functions:

| Rust | Cyrius | Prefix |
|---|---|---|
| `ThermalMaterial` | `UsMaterial` | `us_material_` |
| `GasData` | `UsGas` | `us_gas_` |
| `SubstanceData`, `Phase` | `UsSubstance`, `UsPhase` | `us_substance_`, `us_phase_` |
| `SaturatedEntry`, `WetSteamProperties`, `SuperheatedEntry` | `UsSaturated`, `UsWetSteam`, `UsSuperheated` | `us_saturated_`, `us_wet_steam_`, `us_superheated_` |
| `StatePoint`, `ProcessKind`, `CycleKind`, `CycleResult` | `UsStatePoint`, `UsProcess`, `UsCycle`, `UsCycleResult` | `us_state_point_`, `us_process_`, `us_cycle_`, `us_cycle_result_` |
| `RefrigerationResult`, `CycleComparisonEntry`, `DiagramPoint` | `UsRefrigerationResult`, `UsCycleComparison`, `UsDiagramPoint` | `us_refrigeration_result_`, `us_cycle_comparison_`, `us_diagram_point_` |
| `Species` | `UsSpecies` | `us_species_` |
| `BoundaryCondition`, `Side` | `UsBoundary` (the kind) and `UsBoundaryCondition` (the value), `UsSide` | `us_boundary_`, `us_side_` |
| `ThermalGrid1D`, `ThermalGrid2D`, `ThermalNetwork` | `UsGrid1D`, `UsGrid2D`, `UsNetwork` | `us_grid1d_`, `us_grid2d_`, `us_network_` |
| `ThermalGridVisualization`, `TemperatureProfile` | `UsGridVisualization`, `UsTemperatureProfile` | `us_grid_visualization_`, `us_temperature_profile_` |
| `CycleDiagramData`, `CycleStatePoint` | `UsCycleDiagramData`, `UsCycleStatePoint` | `us_cycle_diagram_data_`, `us_cycle_state_point_` |
| `ThermalNetworkVisualization`, `HeatFluxField` | `UsNetworkVisualization`, `UsHeatFluxField` | `us_network_visualization_`, `us_heat_flux_field_` |
| `DaimonConfig`, `HooshConfig`, `DaimonClient` | `UsDaimonConfig`, `UsHooshConfig`, `UsDaimonClient` | `us_daimon_config_`, `us_hoosh_config_`, `us_daimon_client_` |

`scripts/surface_parity.py` checks that every public Rust item has its Cyrius name. Names starting `_us_` are
private.

## 3. Calls and results

**Results** (row 1). A fallible function returns a value-form `Result`, read with `is_ok` / `is_err_result`:

```rust
let q = transfer::conduction(401.0, 0.01, 373.15, 293.15, 0.05)?;
```
```cyrius
var t, v = us_conduction(401.0, 0.01, 373.15, 293.15, 0.05);
if (is_err_result(t) == 1) { return Err(v); }     # v is the US_ERR_* code
var q: f64 = v;                                     # Ok's payload is the f64's bits
```

The error is a code, not an enum with a payload: `us_error_name(code)` gives Rust's variant name
(`"InvalidTemperature"`) and `us_error_message(code)` a fixed message. The failing value Rust put in the payload goes
to a DEBUG event at the failure site instead (§5). Two codes have no Rust variant: `US_ERR_DESERIALIZE` (serde's
error) and `US_ERR_OUT_OF_MEMORY`, for an allocation Rust would have aborted the process on. Functions that cannot
fail in Rust return their value directly (`us_convection`, every bridge conversion).

**f64 in Cyrius.** Bind a call's f64 result to a `: f64` local before using it, and write float literals as floats
(`300.0`, never `300`): an integer keeps its bits in an f64 slot.

**Values the caller owns** (rows 2, 10, 15–17). Where Rust returns a fresh struct that owns `Vec`s, the Cyrius
function fills one you own, in place, and returns `Ok(0)`. Make it once with `us_<t>_init` (or `_init_a` on an
allocator of your choosing), and reuse it: after the first call, a reused value allocates nothing.

```rust
let r = cycle::otto_cycle(300.0, 101_325.0, 8.0, 50_000.0, 1.4, 1.0)?;
println!("{}", r.efficiency);
```
```cyrius
var r: UsCycleResult;
var ti, vi = us_cycle_result_init(&r);
var t, v = us_otto_cycle(300.0, 101325.0, 8.0, 50000.0, 1.4, 1.0, &r);
var eff: f64 = r.efficiency;
```

**Lists** go in as a pointer and a count (`us_thermal_resistance_series(rs, n)`), and come out into a buffer you
pass with its capacity, the function returning the count (`us_cycle_ts_diagram(r, pts, out, cap)`). **Strings** are
bytes and a length (`name`, `name_len`). An **`Option`** is a flag and a value (`us_rankine_cycle(.., superheated,
t_superheat, ..)`), or a 0 pointer for a string (`api_key`). A fixed **array** argument (`[f64; 3]`) is that many f64
words at a pointer.

## 4. Every divergence, by module

Each row of ADR 0002 says how the golden corpus still compares the function; none is exempt.

### error (row 1)
`UshmaError`'s formatted messages are not reproduced: they are Rust's machinery, not a capability. Match on the code.

### material, state, phase, chem: the reference tables (row 3)
1.4.0's presets are `pub const` values and `ALL_*` slices. The port's are table indices: `US_MATERIAL_COPPER` …,
`US_MATERIAL_COUNT`, `us_material_preset(i, &m)` to fill a `UsMaterial`, and `us_material_name(i)`. Iterate `0 …
US_MATERIAL_COUNT − 1` for `ALL_MATERIALS`, in the same order. The same holds for gases (`US_GAS_*`), substances
(`US_SUBSTANCE_*`, water is `US_SUBSTANCE_WATER_PHASE`) and species (`US_SPECIES_*`).

### logging (row 5)
`USHMA_LOG` takes a level name (`trace` … `error`, `off`) or `ushma=<level>`; other `EnvFilter` directives are
refused with one warning. Call `us_log_init()` (default `info`) or `us_log_init_with_level("debug")`.

### transfer (row 9)
The correlations' fitted exponents (Dittus–Boelter's 0.8 and 0.4, Churchill–Chu's 1/6, 9/16 and −8/27, Pohlhausen's
1/3) are the doubles Rust wrote, through a correctly rounded pow.

### cycle (rows 10–12)
- The ideal-gas cycles fill a `UsCycleResult` you own (row 10): read its state points with
  `us_cycle_state(r, i, &sp)` and its processes with `us_cycle_process(r, i)`.
- `rankine_cycle`'s `Option<f64>` superheat is a flag and a value (row 11).
- `RefrigerationResult` nests its cycle first: `us_refrigeration_result_cycle(rr)` gives its address (row 11).
- The diagrams write `UsDiagramPoint` records into your buffer; size it with `us_cycle_diagram_len` (row 12).
- `compare_cycles` takes its tuples as four parallel lists (row 12).

### chem, stat (rows 13, 14)
- A reaction side is a term list: n pairs of (coefficient, `*UsSpecies`) laid out flat, with n beside it (row 13).
- `partition_rotational`'s symmetry number is an integer, and a negative one is refused (row 14).

### numerical (rows 2, 15, 16)
- **A boundary condition** is a kind (`US_BOUNDARY_FIXED`, `_INSULATED`, `_CONVECTIVE`) and a value,
  `UsBoundaryCondition`, filled by `us_boundary_fixed(&bc, t)` and the others (row 15). A 1D grid keeps its
  Crank–Nicolson scratch, so no step allocates after the first.
- **A 2D grid** stores its nodes flat and row-major. Rust's `nodes[j][i]` is `us_grid2d_node(g, i, j)`, and the
  setter is `us_grid2d_set_node`. A JSON document whose rows aren't ny × nx is refused at decode, where 1.4.0 refuses
  it at the first solve (row 2).
- **A network** keeps its lists and its solve scratch, and `us_network_solve(net, out, cap)` writes the temperatures
  into your buffer (row 16). **One result changes:** a system whose numbers overflow f64 (a resistance below about
  1e-305 K/W, or a temperature above about 1e305 K) is refused, where 1.4.0 returned NaN or ±∞ temperatures as a
  success. Some of those were wrong where the true answer is finite.

### bridge
No divergence. Every conversion keeps Rust's signature, f64s in and one f64 out, and its fallbacks.

### soorat (row 17)
The six visualization types are values you own, filled in place by their constructors
(`us_grid_visualization_from_grid_2d(&v, &grid)` …). Their `Vec`s are flat vecs: a point or flux vector two words,
an edge two, a state point four. Their fixed arrays are fields (`nx`, `ny`, `origin_x` …). The diagram's kind is
`kind`/`kind_len`.

### ai (rows 4, 6, 7, 8, 18)
- **Blocking:** `us_daimon_client_register_agent(&client)` is synchronous, on sandhi, where Rust's was `async` (row 4).
- **The path:** it posts to `/v1/agents`, the path daimon serves, not 1.4.0's `/v1/agents/register`. It returns
  `Ok(id)`, daimon's integer `id`, where 1.4.0 read a string `agent_id` that no daimon sends (row 6).
- **Timeouts:** it bounds the request with `total_ms` and `read_ms`, 30 s each, which are fields of the client
  (row 7).
- **On the wire:** the request is 1.4.0's: its body bytes, `Accept: */*` and the bearer key. sandhi adds its four
  transport headers (row 8).
- **Refusals:** an endpoint or key holding a control byte is refused before anything is sent, and no redirect is
  followed (row 18).
- **The text:** Rust's redacting `Debug` is `us_daimon_config_to_text`, `endpoint=… api_key=<redacted>` (row 18).

## 5. What else to know

**Accuracy** (ADR 0004). `exp`, `expm1`, `ln`, `pow`, `atan` and `tanh` are correctly rounded and give the same bits on
every target. Where glibc misrounds them, a result can differ from 1.4.0's in its last bit: the corpus has 68 such
rows among 154,448. Everything else is 1.4.0's arithmetic in 1.4.0's order, bit for bit.

**JSON** (ADR 0005). Every serde type writes and reads serde_json's documents, byte for byte, so Rust programs' JSON
passes through unchanged. A decoder refuses what serde refuses. It also refuses four things serde reads:
- a duplicate key in an unknown field;
- nesting past 128 levels;
- an integer ≥ 2⁶³;
- an integer `-0`, which reads as +0, not −0.0.

**Logging.** The functions Rust instrumented open DEBUG spans under Rust's span names, and a refusal reports its
value in a DEBUG event (`us_conduction: … conductivity=-1`). With DEBUG off they cost a level test. The daimon key is
never logged.

**Memory.** The global heap never frees, so a value meant to be reused (a grid, a network, a cycle result, a
visualization) is made once and refilled. The `_a` variants take an allocator, such as an arena you reset.

**Examples.** [`../examples/`](../examples/) has the five Rust examples as Cyrius programs on the bundle, each
printing 1.4.0's output (`scripts/examples.sh` checks it on x86_64, PE and aarch64).
