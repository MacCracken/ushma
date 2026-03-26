//! Numerical methods — finite difference, thermal networks, grid solvers.
//!
//! All SI units: kelvins, meters, seconds, W/(m·K).

use serde::{Deserialize, Serialize};

use crate::error::{Result, UshmaError};

/// Boundary condition for a thermal grid edge.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[non_exhaustive]
pub enum BoundaryCondition {
    /// Fixed temperature (Dirichlet).
    Fixed(f64),
    /// Insulated / adiabatic (Neumann, zero flux).
    Insulated,
    /// Convective (Robin): -k·dT/dx = h·(T - T_fluid).
    Convective { h: f64, t_fluid: f64, k: f64 },
}

/// Side of a 2D grid.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum Side {
    Left,
    Right,
    Top,
    Bottom,
}

/// 1D transient conduction grid.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThermalGrid1D {
    /// Temperature at each node (K).
    pub nodes: Vec<f64>,
    /// Node spacing (m).
    pub dx: f64,
    /// Thermal diffusivity α (m²/s).
    pub alpha: f64,
    /// Left boundary condition.
    pub bc_left: BoundaryCondition,
    /// Right boundary condition.
    pub bc_right: BoundaryCondition,
}

impl ThermalGrid1D {
    /// Create a new 1D grid with uniform initial temperature.
    pub fn new(
        num_nodes: usize,
        length: f64,
        alpha: f64,
        t_initial: f64,
        bc_left: BoundaryCondition,
        bc_right: BoundaryCondition,
    ) -> Result<Self> {
        if num_nodes < 3 {
            return Err(UshmaError::InvalidParameter {
                reason: format!("need at least 3 nodes, got {num_nodes}"),
            });
        }
        if length <= 0.0 {
            return Err(UshmaError::InvalidParameter {
                reason: format!("length {length} m must be positive"),
            });
        }
        if alpha <= 0.0 {
            return Err(UshmaError::InvalidParameter {
                reason: format!("thermal diffusivity {alpha} must be positive"),
            });
        }
        let dx = length / (num_nodes - 1) as f64;
        let mut nodes = vec![t_initial; num_nodes];

        // Apply fixed BCs to endpoints
        if let BoundaryCondition::Fixed(t) = bc_left {
            nodes[0] = t;
        }
        if let BoundaryCondition::Fixed(t) = bc_right {
            nodes[num_nodes - 1] = t;
        }

        Ok(Self {
            nodes,
            dx,
            alpha,
            bc_left,
            bc_right,
        })
    }

    /// Fourier number for a given time step.
    #[inline]
    #[must_use]
    pub fn fourier_number(&self, dt: f64) -> f64 {
        self.alpha * dt / (self.dx * self.dx)
    }

    /// Advance one time step using explicit (forward Euler) finite difference.
    ///
    /// Stability requires Fo = α·dt/dx² ≤ 0.5.
    #[tracing::instrument(level = "debug", skip(self))]
    pub fn step_explicit(&mut self, dt: f64) -> Result<()> {
        let fo = self.fourier_number(dt);
        if fo > 0.5 {
            return Err(UshmaError::InvalidParameter {
                reason: format!("Fourier number {fo:.4} exceeds stability limit 0.5; reduce dt"),
            });
        }
        let n = self.nodes.len();
        let mut new = self.nodes.clone();

        // Interior nodes: stencil needs i-1, i, i+1 from old values
        #[allow(clippy::needless_range_loop)]
        for i in 1..n - 1 {
            new[i] =
                self.nodes[i] + fo * (self.nodes[i - 1] - 2.0 * self.nodes[i] + self.nodes[i + 1]);
        }

        // Boundary conditions
        self.apply_bc_left(&mut new, fo);
        self.apply_bc_right(&mut new, fo);

        self.nodes = new;
        Ok(())
    }

    /// Advance one time step using Crank-Nicolson (implicit, unconditionally stable).
    ///
    /// Solves the tridiagonal system via the Thomas algorithm — O(n).
    #[tracing::instrument(level = "debug", skip(self))]
    pub fn step_crank_nicolson(&mut self, dt: f64) -> Result<()> {
        let fo = self.fourier_number(dt);
        let n = self.nodes.len();

        // Tridiagonal coefficients: a_i·T_{i-1} + b_i·T_i + c_i·T_{i+1} = d_i
        let mut a = vec![0.0; n];
        let mut b = vec![0.0; n];
        let mut c = vec![0.0; n];
        let mut d = vec![0.0; n];

        // Interior nodes
        for i in 1..n - 1 {
            a[i] = -fo / 2.0;
            b[i] = 1.0 + fo;
            c[i] = -fo / 2.0;
            d[i] = (fo / 2.0) * self.nodes[i - 1]
                + (1.0 - fo) * self.nodes[i]
                + (fo / 2.0) * self.nodes[i + 1];
        }

        // Boundary conditions
        match self.bc_left {
            BoundaryCondition::Fixed(t) => {
                b[0] = 1.0;
                d[0] = t;
            }
            BoundaryCondition::Insulated => {
                // Ghost node: T_{-1} = T_1
                b[0] = 1.0 + fo;
                c[0] = -fo;
                d[0] = (1.0 - fo) * self.nodes[0] + fo * self.nodes[1];
            }
            BoundaryCondition::Convective { h, t_fluid, k } => {
                let bi = h * self.dx / k;
                b[0] = 1.0 + fo * (1.0 + bi);
                c[0] = -fo;
                d[0] = (1.0 - fo * (1.0 + bi)) * self.nodes[0]
                    + fo * self.nodes[1]
                    + 2.0 * fo * bi * t_fluid;
            }
        }

        match self.bc_right {
            BoundaryCondition::Fixed(t) => {
                b[n - 1] = 1.0;
                d[n - 1] = t;
            }
            BoundaryCondition::Insulated => {
                a[n - 1] = -fo;
                b[n - 1] = 1.0 + fo;
                d[n - 1] = fo * self.nodes[n - 2] + (1.0 - fo) * self.nodes[n - 1];
            }
            BoundaryCondition::Convective { h, t_fluid, k } => {
                let bi = h * self.dx / k;
                a[n - 1] = -fo;
                b[n - 1] = 1.0 + fo * (1.0 + bi);
                d[n - 1] = fo * self.nodes[n - 2]
                    + (1.0 - fo * (1.0 + bi)) * self.nodes[n - 1]
                    + 2.0 * fo * bi * t_fluid;
            }
        }

        // Thomas algorithm (forward sweep + back substitution)
        let mut c_star = vec![0.0; n];
        let mut d_star = vec![0.0; n];

        c_star[0] = c[0] / b[0];
        d_star[0] = d[0] / b[0];

        for i in 1..n {
            let m = b[i] - a[i] * c_star[i - 1];
            if m.abs() < 1e-30 {
                return Err(UshmaError::DivisionByZero {
                    context: "Thomas algorithm pivot is zero".into(),
                });
            }
            c_star[i] = c[i] / m;
            d_star[i] = (d[i] - a[i] * d_star[i - 1]) / m;
        }

        self.nodes[n - 1] = d_star[n - 1];
        for i in (0..n - 1).rev() {
            self.nodes[i] = d_star[i] - c_star[i] * self.nodes[i + 1];
        }

        Ok(())
    }

    /// Advance with adaptive time stepping.
    ///
    /// Subdivides dt if Fo > 0.5 for explicit stability. Returns actual total time advanced.
    pub fn step_adaptive(&mut self, dt_target: f64) -> Result<f64> {
        let fo_max = 0.5;
        let dt_stable = fo_max * self.dx * self.dx / self.alpha;

        if dt_target <= dt_stable {
            self.step_explicit(dt_target)?;
            return Ok(dt_target);
        }

        // Subdivide into stable sub-steps
        let num_steps = (dt_target / dt_stable).ceil() as usize;
        let dt_sub = dt_target / num_steps as f64;

        for _ in 0..num_steps {
            self.step_explicit(dt_sub)?;
        }

        Ok(dt_target)
    }

    fn apply_bc_left(&self, new: &mut [f64], fo: f64) {
        match self.bc_left {
            BoundaryCondition::Fixed(t) => new[0] = t,
            BoundaryCondition::Insulated => {
                // Ghost node symmetry: T_{-1} = T_1
                new[0] = self.nodes[0] + fo * (2.0 * self.nodes[1] - 2.0 * self.nodes[0]);
            }
            BoundaryCondition::Convective { h, t_fluid, k } => {
                let bi = h * self.dx / k;
                new[0] = self.nodes[0]
                    + fo * (2.0 * self.nodes[1] - 2.0 * (1.0 + bi) * self.nodes[0]
                        + 2.0 * bi * t_fluid);
            }
        }
    }

    fn apply_bc_right(&self, new: &mut [f64], fo: f64) {
        let n = self.nodes.len();
        match self.bc_right {
            BoundaryCondition::Fixed(t) => new[n - 1] = t,
            BoundaryCondition::Insulated => {
                new[n - 1] =
                    self.nodes[n - 1] + fo * (2.0 * self.nodes[n - 2] - 2.0 * self.nodes[n - 1]);
            }
            BoundaryCondition::Convective { h, t_fluid, k } => {
                let bi = h * self.dx / k;
                new[n - 1] = self.nodes[n - 1]
                    + fo * (2.0 * self.nodes[n - 2] - 2.0 * (1.0 + bi) * self.nodes[n - 1]
                        + 2.0 * bi * t_fluid);
            }
        }
    }
}

/// 2D steady-state conduction grid.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThermalGrid2D {
    /// Temperature at each node (K). Row-major: `nodes[row][col]`.
    pub nodes: Vec<Vec<f64>>,
    /// Number of columns (x-direction).
    pub nx: usize,
    /// Number of rows (y-direction).
    pub ny: usize,
    /// Node spacing in x (m).
    pub dx: f64,
    /// Node spacing in y (m).
    pub dy: f64,
    /// Boundary conditions: left, right, bottom, top.
    pub bc_left: BoundaryCondition,
    pub bc_right: BoundaryCondition,
    pub bc_bottom: BoundaryCondition,
    pub bc_top: BoundaryCondition,
}

impl ThermalGrid2D {
    /// Create a 2D grid with uniform initial temperature.
    pub fn new(nx: usize, ny: usize, lx: f64, ly: f64, t_initial: f64) -> Result<Self> {
        if nx < 3 || ny < 3 {
            return Err(UshmaError::InvalidParameter {
                reason: format!("need at least 3×3 grid, got {nx}×{ny}"),
            });
        }
        if lx <= 0.0 || ly <= 0.0 {
            return Err(UshmaError::InvalidParameter {
                reason: format!("dimensions must be positive: lx={lx}, ly={ly}"),
            });
        }
        Ok(Self {
            nodes: vec![vec![t_initial; nx]; ny],
            nx,
            ny,
            dx: lx / (nx - 1) as f64,
            dy: ly / (ny - 1) as f64,
            bc_left: BoundaryCondition::Insulated,
            bc_right: BoundaryCondition::Insulated,
            bc_bottom: BoundaryCondition::Insulated,
            bc_top: BoundaryCondition::Insulated,
        })
    }

    /// Set a boundary condition on one side.
    pub fn set_boundary(&mut self, side: Side, bc: BoundaryCondition) {
        match side {
            Side::Left => self.bc_left = bc,
            Side::Right => self.bc_right = bc,
            Side::Bottom => self.bc_bottom = bc,
            Side::Top => self.bc_top = bc,
        }
        self.apply_boundaries();
    }

    /// Solve steady-state via Gauss-Seidel iteration.
    ///
    /// Returns the number of iterations to converge.
    /// For square grids (dx = dy): T_ij = (T_{i-1,j} + T_{i+1,j} + T_{i,j-1} + T_{i,j+1}) / 4.
    #[tracing::instrument(level = "debug", skip(self))]
    pub fn solve_steady_state(&mut self, tolerance: f64, max_iterations: usize) -> Result<usize> {
        let rx = 1.0 / (self.dx * self.dx);
        let ry = 1.0 / (self.dy * self.dy);
        let denom = 2.0 * (rx + ry);

        for iter in 0..max_iterations {
            let mut max_change = 0.0_f64;

            for j in 1..self.ny - 1 {
                for i in 1..self.nx - 1 {
                    let old = self.nodes[j][i];
                    let new_val = (rx * (self.nodes[j][i - 1] + self.nodes[j][i + 1])
                        + ry * (self.nodes[j - 1][i] + self.nodes[j + 1][i]))
                        / denom;
                    self.nodes[j][i] = new_val;
                    max_change = max_change.max((new_val - old).abs());
                }
            }

            // Update non-fixed boundary nodes (insulated = zero gradient)
            self.apply_insulated_boundaries();

            if max_change < tolerance {
                return Ok(iter + 1);
            }
        }

        Err(UshmaError::InvalidParameter {
            reason: format!("Gauss-Seidel did not converge in {max_iterations} iterations"),
        })
    }

    fn apply_insulated_boundaries(&mut self) {
        // Zero gradient: boundary node = nearest interior node
        if matches!(self.bc_left, BoundaryCondition::Insulated) {
            for j in 1..self.ny - 1 {
                self.nodes[j][0] = self.nodes[j][1];
            }
        }
        if matches!(self.bc_right, BoundaryCondition::Insulated) {
            for j in 1..self.ny - 1 {
                self.nodes[j][self.nx - 1] = self.nodes[j][self.nx - 2];
            }
        }
        if matches!(self.bc_bottom, BoundaryCondition::Insulated) {
            for i in 1..self.nx - 1 {
                self.nodes[0][i] = self.nodes[1][i];
            }
        }
        if matches!(self.bc_top, BoundaryCondition::Insulated) {
            for i in 1..self.nx - 1 {
                self.nodes[self.ny - 1][i] = self.nodes[self.ny - 2][i];
            }
        }
    }

    fn apply_boundaries(&mut self) {
        // Left (column 0)
        if let BoundaryCondition::Fixed(t) = self.bc_left {
            for j in 0..self.ny {
                self.nodes[j][0] = t;
            }
        }
        // Right (column nx-1)
        if let BoundaryCondition::Fixed(t) = self.bc_right {
            for j in 0..self.ny {
                self.nodes[j][self.nx - 1] = t;
            }
        }
        // Bottom (row 0)
        if let BoundaryCondition::Fixed(t) = self.bc_bottom {
            for i in 0..self.nx {
                self.nodes[0][i] = t;
            }
        }
        // Top (row ny-1)
        if let BoundaryCondition::Fixed(t) = self.bc_top {
            for i in 0..self.nx {
                self.nodes[self.ny - 1][i] = t;
            }
        }
    }
}

/// Thermal resistance network solver.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThermalNetwork {
    /// Number of nodes.
    pub num_nodes: usize,
    /// Resistance connections: (node_a, node_b, resistance in K/W).
    pub resistances: Vec<(usize, usize, f64)>,
    /// Fixed temperature nodes: (node_index, temperature).
    pub fixed: Vec<(usize, f64)>,
}

impl ThermalNetwork {
    /// Create an empty network with n nodes.
    pub fn new(num_nodes: usize) -> Self {
        Self {
            num_nodes,
            resistances: Vec::new(),
            fixed: Vec::new(),
        }
    }

    /// Add a thermal resistance between two nodes.
    pub fn add_resistance(&mut self, node_a: usize, node_b: usize, resistance: f64) -> Result<()> {
        if node_a >= self.num_nodes || node_b >= self.num_nodes {
            return Err(UshmaError::InvalidParameter {
                reason: format!(
                    "node index out of range: a={node_a}, b={node_b}, max={}",
                    self.num_nodes - 1
                ),
            });
        }
        if resistance <= 0.0 {
            return Err(UshmaError::InvalidParameter {
                reason: format!("resistance {resistance} K/W must be positive"),
            });
        }
        self.resistances.push((node_a, node_b, resistance));
        Ok(())
    }

    /// Set a node to a fixed temperature.
    pub fn set_fixed_temperature(&mut self, node: usize, temperature: f64) -> Result<()> {
        if node >= self.num_nodes {
            return Err(UshmaError::InvalidParameter {
                reason: format!("node {node} out of range, max={}", self.num_nodes - 1),
            });
        }
        self.fixed.push((node, temperature));
        Ok(())
    }

    /// Solve for steady-state temperatures at all nodes.
    ///
    /// Builds a conductance matrix and solves via [`hisab::num::gaussian_elimination`].
    #[tracing::instrument(level = "debug", skip(self))]
    pub fn solve(&self) -> Result<Vec<f64>> {
        let n = self.num_nodes;

        let mut is_fixed = vec![false; n];
        let mut fixed_t = vec![0.0; n];
        for &(node, t) in &self.fixed {
            is_fixed[node] = true;
            fixed_t[node] = t;
        }

        // Count free (unknown) nodes
        let free_indices: Vec<usize> = (0..n).filter(|i| !is_fixed[*i]).collect();
        let nf = free_indices.len();

        if nf == 0 {
            return Ok(fixed_t);
        }

        // Map node index → free index
        let mut free_map = vec![usize::MAX; n];
        for (fi, &ni) in free_indices.iter().enumerate() {
            free_map[ni] = fi;
        }

        // Build augmented matrix [G | b] for free nodes
        let mut matrix = vec![vec![0.0; nf + 1]; nf];

        for &(a, b, r) in &self.resistances {
            let g = 1.0 / r;

            if !is_fixed[a] && !is_fixed[b] {
                let fa = free_map[a];
                let fb = free_map[b];
                matrix[fa][fa] += g;
                matrix[fa][fb] -= g;
                matrix[fb][fb] += g;
                matrix[fb][fa] -= g;
            } else if !is_fixed[a] && is_fixed[b] {
                let fa = free_map[a];
                matrix[fa][fa] += g;
                matrix[fa][nf] += g * fixed_t[b]; // RHS
            } else if is_fixed[a] && !is_fixed[b] {
                let fb = free_map[b];
                matrix[fb][fb] += g;
                matrix[fb][nf] += g * fixed_t[a]; // RHS
            }
        }

        let solution = hisab::num::gaussian_elimination(&mut matrix).map_err(|e| {
            UshmaError::InvalidParameter {
                reason: format!("network solve failed: {e}"),
            }
        })?;

        // Assemble full temperature vector
        let mut temps = fixed_t;
        for (fi, &ni) in free_indices.iter().enumerate() {
            temps[ni] = solution[fi];
        }

        Ok(temps)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- Grid 1D construction ---

    #[test]
    fn test_grid1d_construction() {
        let g = ThermalGrid1D::new(
            10,
            1.0,
            1e-4,
            300.0,
            BoundaryCondition::Fixed(400.0),
            BoundaryCondition::Fixed(300.0),
        )
        .unwrap();
        assert_eq!(g.nodes.len(), 10);
        assert_eq!(g.nodes[0], 400.0); // fixed left
        assert_eq!(g.nodes[9], 300.0); // fixed right
        assert_eq!(g.nodes[5], 300.0); // interior = initial
    }

    #[test]
    fn test_grid1d_invalid() {
        assert!(
            ThermalGrid1D::new(
                2,
                1.0,
                1e-4,
                300.0,
                BoundaryCondition::Insulated,
                BoundaryCondition::Insulated
            )
            .is_err()
        );
        assert!(
            ThermalGrid1D::new(
                10,
                0.0,
                1e-4,
                300.0,
                BoundaryCondition::Insulated,
                BoundaryCondition::Insulated
            )
            .is_err()
        );
        assert!(
            ThermalGrid1D::new(
                10,
                1.0,
                0.0,
                300.0,
                BoundaryCondition::Insulated,
                BoundaryCondition::Insulated
            )
            .is_err()
        );
    }

    // --- 1D explicit ---

    #[test]
    fn test_explicit_uniform_stays_uniform() {
        let mut g = ThermalGrid1D::new(
            20,
            1.0,
            1e-4,
            300.0,
            BoundaryCondition::Fixed(300.0),
            BoundaryCondition::Fixed(300.0),
        )
        .unwrap();
        for _ in 0..100 {
            g.step_explicit(0.1).unwrap();
        }
        for &t in &g.nodes {
            assert!((t - 300.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_explicit_converges_to_linear() {
        // Fixed 400K left, 300K right → steady state is linear
        let mut g = ThermalGrid1D::new(
            21,
            1.0,
            1e-4,
            350.0,
            BoundaryCondition::Fixed(400.0),
            BoundaryCondition::Fixed(300.0),
        )
        .unwrap();
        let dt = 0.4 * g.dx * g.dx / g.alpha; // Fo < 0.5
        for _ in 0..100_000 {
            g.step_explicit(dt).unwrap();
        }
        // Check linearity at midpoint
        let t_mid = g.nodes[10];
        assert!(
            (t_mid - 350.0).abs() < 1.0,
            "midpoint T={t_mid}, expected ~350"
        );
    }

    #[test]
    fn test_explicit_stability_rejected() {
        let g = ThermalGrid1D::new(
            10,
            1.0,
            1e-4,
            300.0,
            BoundaryCondition::Insulated,
            BoundaryCondition::Insulated,
        )
        .unwrap();
        // dx = 1/9 ≈ 0.111, α=1e-4, Fo = α·dt/dx² > 0.5 → dt > 0.5*dx²/α ≈ 61.7
        let dt_unstable = 100.0;
        let mut g2 = g;
        assert!(g2.step_explicit(dt_unstable).is_err());
    }

    // --- Crank-Nicolson ---

    #[test]
    fn test_crank_nicolson_converges_to_linear() {
        let mut g = ThermalGrid1D::new(
            21,
            1.0,
            1e-4,
            350.0,
            BoundaryCondition::Fixed(400.0),
            BoundaryCondition::Fixed(300.0),
        )
        .unwrap();
        // Large dt — CN is unconditionally stable
        for _ in 0..10_000 {
            g.step_crank_nicolson(1.0).unwrap();
        }
        let t_mid = g.nodes[10];
        assert!(
            (t_mid - 350.0).abs() < 1.0,
            "CN midpoint T={t_mid}, expected ~350"
        );
    }

    #[test]
    fn test_crank_nicolson_insulated() {
        // Insulated both ends, initial hot spot → should equilibrate
        let mut g = ThermalGrid1D::new(
            11,
            1.0,
            1e-4,
            300.0,
            BoundaryCondition::Insulated,
            BoundaryCondition::Insulated,
        )
        .unwrap();
        g.nodes[5] = 500.0; // hot spot in middle
        let total_energy_before: f64 = g.nodes.iter().sum();
        for _ in 0..10_000 {
            g.step_crank_nicolson(0.1).unwrap();
        }
        let total_energy_after: f64 = g.nodes.iter().sum();
        // Energy conservation within tolerance for insulated BCs
        let rel_err = (total_energy_before - total_energy_after).abs() / total_energy_before;
        assert!(rel_err < 0.01, "energy drift {rel_err:.4}");
        // Should equilibrate to uniform
        let t_avg = total_energy_after / g.nodes.len() as f64;
        for &t in &g.nodes {
            assert!(
                (t - t_avg).abs() < 1.0,
                "not equilibrated: T={t}, avg={t_avg}"
            );
        }
    }

    // --- Adaptive ---

    #[test]
    fn test_adaptive_small_dt_passes_through() {
        let mut g = ThermalGrid1D::new(
            10,
            1.0,
            1e-4,
            300.0,
            BoundaryCondition::Fixed(400.0),
            BoundaryCondition::Fixed(300.0),
        )
        .unwrap();
        let dt_small = 0.1 * g.dx * g.dx / g.alpha;
        let actual = g.step_adaptive(dt_small).unwrap();
        assert!((actual - dt_small).abs() < 1e-10);
    }

    #[test]
    fn test_adaptive_large_dt_subdivides() {
        let mut g = ThermalGrid1D::new(
            10,
            1.0,
            1e-4,
            300.0,
            BoundaryCondition::Fixed(400.0),
            BoundaryCondition::Fixed(300.0),
        )
        .unwrap();
        let dt_large = 10.0 * g.dx * g.dx / g.alpha; // Fo >> 0.5
        let actual = g.step_adaptive(dt_large).unwrap();
        assert!((actual - dt_large).abs() < 1e-10);
    }

    // --- 2D grid ---

    #[test]
    fn test_grid2d_construction() {
        let g = ThermalGrid2D::new(10, 10, 1.0, 1.0, 300.0).unwrap();
        assert_eq!(g.nodes.len(), 10);
        assert_eq!(g.nodes[0].len(), 10);
    }

    #[test]
    fn test_grid2d_uniform_boundaries() {
        let mut g = ThermalGrid2D::new(11, 11, 1.0, 1.0, 300.0).unwrap();
        g.set_boundary(Side::Left, BoundaryCondition::Fixed(400.0));
        g.set_boundary(Side::Right, BoundaryCondition::Fixed(400.0));
        g.set_boundary(Side::Top, BoundaryCondition::Fixed(400.0));
        g.set_boundary(Side::Bottom, BoundaryCondition::Fixed(400.0));
        let iters = g.solve_steady_state(1e-6, 10_000).unwrap();
        assert!(iters > 0);
        // All interior should be 400 K
        for j in 1..10 {
            for i in 1..10 {
                assert!(
                    (g.nodes[j][i] - 400.0).abs() < 0.01,
                    "T[{j}][{i}]={}",
                    g.nodes[j][i]
                );
            }
        }
    }

    #[test]
    fn test_grid2d_two_fixed_sides() {
        let mut g = ThermalGrid2D::new(21, 21, 1.0, 1.0, 350.0).unwrap();
        g.set_boundary(Side::Left, BoundaryCondition::Fixed(400.0));
        g.set_boundary(Side::Right, BoundaryCondition::Fixed(300.0));
        g.set_boundary(Side::Top, BoundaryCondition::Insulated);
        g.set_boundary(Side::Bottom, BoundaryCondition::Insulated);
        let iters = g.solve_steady_state(1e-6, 50_000).unwrap();
        assert!(iters > 0);
        // Midpoint should be ~350 K (linear profile in x)
        let t_mid = g.nodes[10][10];
        assert!((t_mid - 350.0).abs() < 2.0, "2D midpoint T={t_mid}");
    }

    #[test]
    fn test_grid2d_invalid() {
        assert!(ThermalGrid2D::new(2, 10, 1.0, 1.0, 300.0).is_err());
        assert!(ThermalGrid2D::new(10, 10, 0.0, 1.0, 300.0).is_err());
    }

    // --- Thermal network ---

    #[test]
    fn test_network_series() {
        // Three nodes: 0—R₁—1—R₂—2, fixed T at 0 and 2
        let mut net = ThermalNetwork::new(3);
        net.add_resistance(0, 1, 1.0).unwrap();
        net.add_resistance(1, 2, 1.0).unwrap();
        net.set_fixed_temperature(0, 400.0).unwrap();
        net.set_fixed_temperature(2, 300.0).unwrap();
        let temps = net.solve().unwrap();
        assert!((temps[0] - 400.0).abs() < 1e-6);
        assert!((temps[1] - 350.0).abs() < 1e-6);
        assert!((temps[2] - 300.0).abs() < 1e-6);
    }

    #[test]
    fn test_network_parallel() {
        // Two paths from 0 to 2: 0—R₁—1—R₂—2 and 0—R₃—2
        let mut net = ThermalNetwork::new(3);
        net.add_resistance(0, 1, 2.0).unwrap();
        net.add_resistance(1, 2, 2.0).unwrap();
        net.add_resistance(0, 2, 2.0).unwrap();
        net.set_fixed_temperature(0, 400.0).unwrap();
        net.set_fixed_temperature(2, 300.0).unwrap();
        let temps = net.solve().unwrap();
        // Node 1: connected to 0 via R=2 and 2 via R=2 → T1 = (400/2 + 300/2)/(1/2+1/2) = 350
        assert!((temps[1] - 350.0).abs() < 1e-6);
    }

    #[test]
    fn test_network_invalid() {
        let mut net = ThermalNetwork::new(3);
        assert!(net.add_resistance(0, 5, 1.0).is_err()); // out of range
        assert!(net.add_resistance(0, 1, 0.0).is_err()); // zero R
        assert!(net.set_fixed_temperature(5, 300.0).is_err()); // out of range
    }

    // --- BC serde ---

    #[test]
    fn test_boundary_condition_serde() {
        let bc = BoundaryCondition::Convective {
            h: 25.0,
            t_fluid: 293.15,
            k: 401.0,
        };
        let json = serde_json::to_string(&bc).unwrap();
        let back: BoundaryCondition = serde_json::from_str(&json).unwrap();
        if let BoundaryCondition::Convective { h, .. } = back {
            assert!((h - 25.0).abs() < 1e-10);
        } else {
            panic!("wrong BC variant");
        }
    }
}
