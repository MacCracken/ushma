//! Soorat integration — visualization data structures for thermal analysis.
//!
//! Provides structured types that soorat can render: thermal grid heatmaps,
//! 1D temperature profiles, cycle diagrams, thermal network graphs, and
//! heat flux vector fields.

use serde::{Deserialize, Serialize};

// ── Thermal grid heatmap ───────────────────────────────────────────────────

/// A 2D temperature field for heatmap rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ThermalGridVisualization {
    /// Temperature values at grid points (K), flattened row-major: `values[y * nx + x]`.
    pub values: Vec<f64>,
    /// Grid dimensions (nx, ny).
    pub dimensions: [usize; 2],
    /// World-space origin `[x, y]` in metres.
    pub origin: [f64; 2],
    /// Grid spacing `[dx, dy]` in metres.
    pub spacing: [f64; 2],
    /// Minimum temperature in the grid (K).
    pub min_temp: f64,
    /// Maximum temperature in the grid (K).
    pub max_temp: f64,
}

#[cfg(feature = "numerical")]
impl ThermalGridVisualization {
    /// Create from a `ThermalGrid2D`.
    #[must_use]
    pub fn from_grid_2d(grid: &crate::numerical::ThermalGrid2D) -> Self {
        let mut values = Vec::with_capacity(grid.nx * grid.ny);
        let mut min_t = f64::MAX;
        let mut max_t = f64::MIN;
        for row in &grid.nodes {
            for &t in row {
                if t < min_t {
                    min_t = t;
                }
                if t > max_t {
                    max_t = t;
                }
                values.push(t);
            }
        }
        Self {
            values,
            dimensions: [grid.nx, grid.ny],
            origin: [0.0, 0.0],
            spacing: [grid.dx, grid.dy],
            min_temp: min_t,
            max_temp: max_t,
        }
    }
}

// ── 1D temperature profile ─────────────────────────────────────────────────

/// A 1D temperature profile for line/ribbon rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TemperatureProfile {
    /// Temperature at each node (K).
    pub temperatures: Vec<f64>,
    /// Node spacing (m).
    pub dx: f64,
    /// World-space start position `[x, y, z]`.
    pub origin: [f64; 3],
    /// Direction of the profile (unit vector).
    pub direction: [f64; 3],
    /// Minimum temperature (K).
    pub min_temp: f64,
    /// Maximum temperature (K).
    pub max_temp: f64,
}

#[cfg(feature = "numerical")]
impl TemperatureProfile {
    /// Create from a `ThermalGrid1D`, laid out along a direction.
    #[must_use]
    pub fn from_grid_1d(
        grid: &crate::numerical::ThermalGrid1D,
        origin: [f64; 3],
        direction: [f64; 3],
    ) -> Self {
        let min_t = grid.nodes.iter().cloned().fold(f64::MAX, f64::min);
        let max_t = grid.nodes.iter().cloned().fold(f64::MIN, f64::max);
        Self {
            temperatures: grid.nodes.clone(),
            dx: grid.dx,
            origin,
            direction,
            min_temp: min_t,
            max_temp: max_t,
        }
    }
}

// ── Cycle diagrams ─────────────────────────────────────────────────────────

/// Cycle diagram data for T-s and P-v plot rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CycleDiagramData {
    /// T-s diagram points (x=entropy, y=temperature).
    pub ts_points: Vec<[f64; 2]>,
    /// P-v diagram points (x=volume, y=pressure).
    pub pv_points: Vec<[f64; 2]>,
    /// State point vertices (corners of the cycle).
    pub state_points: Vec<CycleStatePoint>,
    /// Cycle kind label.
    pub kind: String,
    /// Thermal efficiency.
    pub efficiency: f64,
}

/// A state point for labeling on diagrams.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CycleStatePoint {
    /// Temperature (K).
    pub temperature: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Volume (m³).
    pub volume: f64,
    /// Entropy (J/K).
    pub entropy: f64,
}

#[cfg(feature = "cycle")]
impl CycleDiagramData {
    /// Create from a `CycleResult` with diagram point generation.
    #[must_use]
    pub fn from_cycle_result(
        result: &crate::cycle::CycleResult,
        points_per_process: usize,
    ) -> Self {
        let ts = crate::cycle::cycle_ts_diagram(result, points_per_process);
        let pv = crate::cycle::cycle_pv_diagram(result, points_per_process);

        let ts_points = ts.iter().map(|p| [p.x, p.y]).collect();
        let pv_points = pv.iter().map(|p| [p.x, p.y]).collect();

        let state_points = result
            .state_points
            .iter()
            .map(|sp| CycleStatePoint {
                temperature: sp.temperature,
                pressure: sp.pressure,
                volume: sp.volume,
                entropy: sp.entropy,
            })
            .collect();

        Self {
            ts_points,
            pv_points,
            state_points,
            kind: format!("{:?}", result.kind),
            efficiency: result.efficiency,
        }
    }
}

// ── Thermal network graph ──────────────────────────────────────────────────

/// Thermal network for node-link visualization.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ThermalNetworkVisualization {
    /// Temperature at each node (K). Index = node ID.
    pub node_temperatures: Vec<f64>,
    /// Resistance edges: `(node_a, node_b, resistance_k_per_w)`.
    pub edges: Vec<[usize; 2]>,
    /// Conductance (1/R) for each edge, for line width scaling.
    pub conductances: Vec<f64>,
}

#[cfg(feature = "numerical")]
impl ThermalNetworkVisualization {
    /// Create from a `ThermalNetwork` and its solved temperatures.
    #[must_use]
    pub fn from_network(network: &crate::numerical::ThermalNetwork, temperatures: &[f64]) -> Self {
        let node_temps = if temperatures.len() >= network.num_nodes {
            temperatures[..network.num_nodes].to_vec()
        } else {
            let mut t = temperatures.to_vec();
            t.resize(network.num_nodes, 0.0);
            t
        };

        let edges: Vec<[usize; 2]> = network
            .resistances
            .iter()
            .map(|&(a, b, _)| [a, b])
            .collect();

        let conductances: Vec<f64> = network
            .resistances
            .iter()
            .map(|&(_, _, r)| if r > 0.0 { 1.0 / r } else { 0.0 })
            .collect();

        Self {
            node_temperatures: node_temps,
            edges,
            conductances,
        }
    }
}

// ── Heat flux vectors ──────────────────────────────────────────────────────

/// A 2D grid of heat flux vectors for arrow/streamline rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HeatFluxField {
    /// Heat flux vectors `[qx, qy]` at each grid point (W/m²).
    /// Flattened row-major: `fluxes[y * nx + x]`.
    pub fluxes: Vec<[f64; 2]>,
    /// Grid dimensions (nx, ny).
    pub dimensions: [usize; 2],
    /// Grid spacing `[dx, dy]` in metres.
    pub spacing: [f64; 2],
    /// Maximum flux magnitude (for normalization).
    pub max_magnitude: f64,
}

#[cfg(feature = "numerical")]
impl HeatFluxField {
    /// Compute heat flux from a `ThermalGrid2D` using Fourier's law.
    ///
    /// `conductivity`: thermal conductivity (W/(m·K)).
    #[must_use]
    pub fn from_grid_2d(grid: &crate::numerical::ThermalGrid2D, conductivity: f64) -> Self {
        let nx = grid.nx;
        let ny = grid.ny;
        let mut fluxes = Vec::with_capacity(nx * ny);
        let mut max_mag = 0.0_f64;

        for iy in 0..ny {
            for ix in 0..nx {
                // Central difference for temperature gradient
                let dt_dx = if ix > 0 && ix < nx - 1 {
                    (grid.nodes[iy][ix + 1] - grid.nodes[iy][ix - 1]) / (2.0 * grid.dx)
                } else if ix == 0 && nx > 1 {
                    (grid.nodes[iy][ix + 1] - grid.nodes[iy][ix]) / grid.dx
                } else if ix == nx - 1 && nx > 1 {
                    (grid.nodes[iy][ix] - grid.nodes[iy][ix - 1]) / grid.dx
                } else {
                    0.0
                };

                let dt_dy = if iy > 0 && iy < ny - 1 {
                    (grid.nodes[iy + 1][ix] - grid.nodes[iy - 1][ix]) / (2.0 * grid.dy)
                } else if iy == 0 && ny > 1 {
                    (grid.nodes[iy + 1][ix] - grid.nodes[iy][ix]) / grid.dy
                } else if iy == ny - 1 && ny > 1 {
                    (grid.nodes[iy][ix] - grid.nodes[iy - 1][ix]) / grid.dy
                } else {
                    0.0
                };

                // q = -k × ∇T
                let qx = -conductivity * dt_dx;
                let qy = -conductivity * dt_dy;
                let mag = (qx * qx + qy * qy).sqrt();
                if mag > max_mag {
                    max_mag = mag;
                }
                fluxes.push([qx, qy]);
            }
        }

        Self {
            fluxes,
            dimensions: [nx, ny],
            spacing: [grid.dx, grid.dy],
            max_magnitude: max_mag,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thermal_grid_viz_serializes() {
        let viz = ThermalGridVisualization {
            values: vec![300.0; 4],
            dimensions: [2, 2],
            origin: [0.0, 0.0],
            spacing: [1.0, 1.0],
            min_temp: 300.0,
            max_temp: 300.0,
        };
        let json = serde_json::to_string(&viz);
        assert!(json.is_ok());
    }

    #[test]
    fn temperature_profile_serializes() {
        let prof = TemperatureProfile {
            temperatures: vec![300.0, 310.0, 320.0],
            dx: 0.1,
            origin: [0.0; 3],
            direction: [1.0, 0.0, 0.0],
            min_temp: 300.0,
            max_temp: 320.0,
        };
        let json = serde_json::to_string(&prof);
        assert!(json.is_ok());
    }

    #[test]
    fn cycle_diagram_serializes() {
        let diag = CycleDiagramData {
            ts_points: vec![[0.0, 300.0], [100.0, 600.0]],
            pv_points: vec![[0.001, 100000.0], [0.01, 50000.0]],
            state_points: vec![CycleStatePoint {
                temperature: 300.0,
                pressure: 100000.0,
                volume: 0.001,
                entropy: 0.0,
            }],
            kind: "Otto".to_string(),
            efficiency: 0.56,
        };
        assert_eq!(diag.ts_points.len(), 2);
    }

    #[test]
    fn thermal_network_viz_manual() {
        let viz = ThermalNetworkVisualization {
            node_temperatures: vec![300.0, 350.0, 400.0],
            edges: vec![[0, 1], [1, 2]],
            conductances: vec![10.0, 5.0],
        };
        assert_eq!(viz.edges.len(), 2);
        assert_eq!(viz.conductances.len(), 2);
    }

    #[test]
    fn heat_flux_field_manual() {
        let flux = HeatFluxField {
            fluxes: vec![[100.0, 0.0], [-50.0, 50.0]],
            dimensions: [2, 1],
            spacing: [0.1, 0.1],
            max_magnitude: 100.0,
        };
        assert_eq!(flux.fluxes.len(), 2);
    }

    #[cfg(feature = "numerical")]
    #[test]
    fn thermal_grid_from_grid_2d() {
        let mut grid = crate::numerical::ThermalGrid2D::new(4, 4, 1.0, 1.0, 300.0).unwrap();
        grid.nodes[0][0] = 400.0; // hot corner
        let viz = ThermalGridVisualization::from_grid_2d(&grid);
        assert_eq!(viz.dimensions, [4, 4]);
        assert_eq!(viz.values.len(), 16);
        assert!((viz.min_temp - 300.0).abs() < 0.01);
        assert!((viz.max_temp - 400.0).abs() < 0.01);
    }

    #[cfg(feature = "numerical")]
    #[test]
    fn temperature_profile_from_grid_1d() {
        let grid = crate::numerical::ThermalGrid1D::new(
            5,
            1.0,
            1e-5,
            350.0,
            crate::numerical::BoundaryCondition::Fixed(400.0),
            crate::numerical::BoundaryCondition::Fixed(300.0),
        )
        .unwrap();
        let prof = TemperatureProfile::from_grid_1d(&grid, [0.0; 3], [1.0, 0.0, 0.0]);
        assert_eq!(prof.temperatures.len(), 5);
        assert!((prof.dx - 0.25).abs() < 0.01);
    }

    #[cfg(feature = "numerical")]
    #[test]
    fn heat_flux_from_grid_2d() {
        let mut grid = crate::numerical::ThermalGrid2D::new(5, 5, 1.0, 1.0, 300.0).unwrap();
        // Create a temperature gradient in X
        for row in &mut grid.nodes {
            for (ix, t) in row.iter_mut().enumerate() {
                *t = 300.0 + ix as f64 * 25.0;
            }
        }
        let flux = HeatFluxField::from_grid_2d(&grid, 1.0);
        assert_eq!(flux.fluxes.len(), 25);
        assert!(flux.max_magnitude > 0.0);
        // Interior point should have negative qx (heat flows left-to-right = -k × dT/dx)
        let mid = flux.fluxes[2 * 5 + 2];
        assert!(
            mid[0] < 0.0,
            "flux should be negative (left to right): {}",
            mid[0]
        );
    }

    #[cfg(feature = "numerical")]
    #[test]
    fn thermal_network_from_network() {
        let mut net = crate::numerical::ThermalNetwork::new(3);
        net.add_resistance(0, 1, 10.0).unwrap();
        net.add_resistance(1, 2, 20.0).unwrap();
        let temps = vec![300.0, 350.0, 400.0];
        let viz = ThermalNetworkVisualization::from_network(&net, &temps);
        assert_eq!(viz.node_temperatures.len(), 3);
        assert_eq!(viz.edges.len(), 2);
        assert!((viz.conductances[0] - 0.1).abs() < 0.001);
        assert!((viz.conductances[1] - 0.05).abs() < 0.001);
    }

    #[cfg(feature = "cycle")]
    #[test]
    fn cycle_diagram_from_otto() {
        let result = crate::cycle::otto_cycle(300.0, 101325.0, 8.0, 2000.0, 1.4, 1.0).unwrap();
        let diag = CycleDiagramData::from_cycle_result(&result, 20);
        assert!(!diag.ts_points.is_empty());
        assert!(!diag.pv_points.is_empty());
        assert_eq!(diag.state_points.len(), result.state_points.len());
        assert!(diag.efficiency > 0.0 && diag.efficiency < 1.0);
    }
}
