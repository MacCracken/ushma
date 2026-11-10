//! Numerical methods — finite difference, thermal networks, grid solvers.
//!
//! All SI units: kelvins, meters, seconds, W/(m·K).

use serde::{Deserialize, Serialize};

use crate::error::{Result, UshmaError};
use crate::validate::{finite_ge, finite_gt, invalid_param};

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

impl BoundaryCondition {
    /// Check that the condition's parameters are physical and finite.
    fn validate(&self) -> Result<()> {
        match *self {
            BoundaryCondition::Fixed(t) if !t.is_finite() => Err(invalid_param(format_args!(
                "fixed boundary temperature {t} K must be finite"
            ))),
            BoundaryCondition::Convective { h, t_fluid, k }
                if !finite_ge(h, 0.0) || !finite_gt(k, 0.0) || !t_fluid.is_finite() =>
            {
                Err(invalid_param(format_args!(
                    "convective boundary needs finite h ≥ 0, k > 0 and T_fluid: h={h}, k={k}, T_fluid={t_fluid}"
                )))
            }
            _ => Ok(()),
        }
    }

    /// Cell Biot number h·Δ/k for a convective boundary, 0 otherwise.
    #[inline]
    fn biot(&self, spacing: f64) -> f64 {
        match *self {
            BoundaryCondition::Convective { h, k, .. } => h * spacing / k,
            _ => 0.0,
        }
    }
}

/// Upper limit on explicit sub-steps taken by [`ThermalGrid1D::step_adaptive`].
pub const MAX_ADAPTIVE_SUBSTEPS: usize = 1_000_000;

/// Relative slack on the explicit stability limit, so a step computed as
/// exactly `max_stable_dt()` is not rejected because of round-off.
const STABILITY_ROUNDOFF: f64 = 1e-12;

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
            return Err(invalid_param(format_args!(
                "need at least 3 nodes, got {num_nodes}"
            )));
        }
        if !finite_gt(length, 0.0) {
            return Err(invalid_param(format_args!(
                "length {length} m must be positive"
            )));
        }
        if !finite_gt(alpha, 0.0) {
            return Err(invalid_param(format_args!(
                "thermal diffusivity {alpha} must be positive"
            )));
        }
        if !t_initial.is_finite() {
            return Err(invalid_param(format_args!(
                "initial temperature {t_initial} K must be finite"
            )));
        }
        bc_left.validate()?;
        bc_right.validate()?;
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

    /// Largest Fourier number for which the explicit scheme is stable:
    /// 0.5 / (1 + Bi_max), where Bi = h·dx/k of a convective boundary node.
    #[inline]
    fn max_stable_fourier(&self) -> f64 {
        let bi = self.bc_left.biot(self.dx).max(self.bc_right.biot(self.dx));
        0.5 / (1.0 + bi)
    }

    /// Largest stable explicit time step (s): Fo·(1 + Bi_max) ≤ 0.5.
    ///
    /// Interior nodes need Fo ≤ 0.5; a convective boundary node needs the
    /// stricter Fo·(1 + Bi) ≤ 0.5 (Incropera, Table 5.3).
    #[inline]
    #[must_use]
    pub fn max_stable_dt(&self) -> f64 {
        self.max_stable_fourier() * self.dx * self.dx / self.alpha
    }

    /// Validate grid shape and parameters (fields are public and the type is
    /// deserializable, so the constructor's checks may have been bypassed).
    fn check(&self) -> Result<()> {
        if self.nodes.len() < 3 {
            return Err(invalid_param(format_args!(
                "need at least 3 nodes, got {}",
                self.nodes.len()
            )));
        }
        if !finite_gt(self.dx, 0.0) || !finite_gt(self.alpha, 0.0) {
            return Err(invalid_param(format_args!(
                "dx={} and alpha={} must be positive and finite",
                self.dx, self.alpha
            )));
        }
        self.bc_left.validate()?;
        self.bc_right.validate()
    }

    fn check_dt(dt: f64) -> Result<()> {
        if !finite_gt(dt, 0.0) {
            return Err(invalid_param(format_args!(
                "time step {dt} s must be positive and finite"
            )));
        }
        Ok(())
    }

    /// Advance one time step using explicit (forward Euler) finite difference.
    ///
    /// Stability requires Fo = α·dt/dx² ≤ 0.5, or Fo·(1 + Bi) ≤ 0.5 with a
    /// convective boundary — see [`max_stable_dt`](Self::max_stable_dt).
    #[tracing::instrument(level = "debug", skip(self))]
    pub fn step_explicit(&mut self, dt: f64) -> Result<()> {
        self.check()?;
        Self::check_dt(dt)?;
        let fo = self.fourier_number(dt);
        let fo_max = self.max_stable_fourier();
        if fo > fo_max * (1.0 + STABILITY_ROUNDOFF) {
            return Err(invalid_param(format_args!(
                "Fourier number {fo:.4} exceeds stability limit {fo_max:.4}; reduce dt to ≤ {:.6e} s",
                self.max_stable_dt()
            )));
        }
        self.explicit_update(fo);
        Ok(())
    }

    /// One forward-Euler update in place (no allocation). Caller checks stability.
    fn explicit_update(&mut self, fo: f64) {
        let n = self.nodes.len();
        let nodes = &mut self.nodes;
        // Old boundary-adjacent values needed by the boundary updates.
        let (old_0, old_1) = (nodes[0], nodes[1]);
        let (old_nm2, old_nm1) = (nodes[n - 2], nodes[n - 1]);

        // Interior nodes; `prev` carries the old value of node i-1.
        let mut prev = old_0;
        for i in 1..n - 1 {
            let cur = nodes[i];
            nodes[i] = cur + fo * (prev - 2.0 * cur + nodes[i + 1]);
            prev = cur;
        }

        nodes[0] = Self::explicit_boundary(self.bc_left, old_0, old_1, fo, self.dx);
        nodes[n - 1] = Self::explicit_boundary(self.bc_right, old_nm1, old_nm2, fo, self.dx);
    }

    /// Explicit update of a boundary node from its old value and its neighbor's
    /// old value (ghost-node / half-cell energy balance).
    #[inline]
    fn explicit_boundary(bc: BoundaryCondition, t_b: f64, t_in: f64, fo: f64, dx: f64) -> f64 {
        match bc {
            BoundaryCondition::Fixed(t) => t,
            // Ghost node symmetry: T_ghost = T_in
            BoundaryCondition::Insulated => t_b + fo * (2.0 * t_in - 2.0 * t_b),
            BoundaryCondition::Convective { t_fluid, .. } => {
                let bi = bc.biot(dx);
                t_b + fo * (2.0 * t_in - 2.0 * (1.0 + bi) * t_b + 2.0 * bi * t_fluid)
            }
        }
    }

    /// Crank-Nicolson row i: (a_i, b_i, c_i, d_i) of a·T_{i-1} + b·T_i + c·T_{i+1} = d.
    #[inline]
    fn cn_row(&self, i: usize, fo: f64) -> (f64, f64, f64, f64) {
        let t = &self.nodes;
        let n = t.len();
        let boundary = |bc: BoundaryCondition, t_b: f64, t_in: f64| -> (f64, f64, f64) {
            // Returns (diag, off-diagonal toward the interior, rhs).
            match bc {
                BoundaryCondition::Fixed(tf) => (1.0, 0.0, tf),
                // Ghost node: T_ghost = T_in
                BoundaryCondition::Insulated => (1.0 + fo, -fo, (1.0 - fo) * t_b + fo * t_in),
                BoundaryCondition::Convective { t_fluid, .. } => {
                    let bi = bc.biot(self.dx);
                    (
                        1.0 + fo * (1.0 + bi),
                        -fo,
                        (1.0 - fo * (1.0 + bi)) * t_b + fo * t_in + 2.0 * fo * bi * t_fluid,
                    )
                }
            }
        };
        if i == 0 {
            let (b, c, d) = boundary(self.bc_left, t[0], t[1]);
            (0.0, b, c, d)
        } else if i == n - 1 {
            let (b, a, d) = boundary(self.bc_right, t[n - 1], t[n - 2]);
            (a, b, 0.0, d)
        } else {
            let half = 0.5 * fo;
            (
                -half,
                1.0 + fo,
                -half,
                half * t[i - 1] + (1.0 - fo) * t[i] + half * t[i + 1],
            )
        }
    }

    /// Advance one time step using Crank-Nicolson (implicit, unconditionally stable).
    ///
    /// Solves the tridiagonal system via the Thomas algorithm — O(n), two
    /// scratch vectors.
    #[tracing::instrument(level = "debug", skip(self))]
    pub fn step_crank_nicolson(&mut self, dt: f64) -> Result<()> {
        self.check()?;
        Self::check_dt(dt)?;
        let fo = self.fourier_number(dt);
        let n = self.nodes.len();

        // Forward sweep. Rows are built from the old temperatures, which are not
        // overwritten until back substitution.
        let mut c_star = vec![0.0; n];
        let mut d_star = vec![0.0; n];
        let mut c_prev = 0.0;
        let mut d_prev = 0.0;
        for i in 0..n {
            let (a, b, c, d) = self.cn_row(i, fo);
            let m = b - a * c_prev;
            if !finite_ge(m.abs(), 1e-30) {
                return Err(UshmaError::DivisionByZero {
                    context: "Thomas algorithm pivot is zero".into(),
                });
            }
            c_prev = c / m;
            d_prev = (d - a * d_prev) / m;
            c_star[i] = c_prev;
            d_star[i] = d_prev;
        }

        self.nodes[n - 1] = d_star[n - 1];
        for i in (0..n - 1).rev() {
            self.nodes[i] = d_star[i] - c_star[i] * self.nodes[i + 1];
        }

        Ok(())
    }

    /// Advance with adaptive time stepping.
    ///
    /// Subdivides `dt_target` into equal explicit sub-steps no longer than
    /// [`max_stable_dt`](Self::max_stable_dt). Returns the total time advanced.
    /// Errors if more than [`MAX_ADAPTIVE_SUBSTEPS`] sub-steps would be needed —
    /// use [`step_crank_nicolson`](Self::step_crank_nicolson) for such stiff cases.
    pub fn step_adaptive(&mut self, dt_target: f64) -> Result<f64> {
        self.check()?;
        Self::check_dt(dt_target)?;
        let dt_stable = self.max_stable_dt();

        // The (1 − ε) factor keeps exact multiples of dt_stable from gaining an
        // extra sub-step through round-off.
        let ratio = dt_target / dt_stable;
        let steps = (ratio * (1.0 - STABILITY_ROUNDOFF)).ceil().max(1.0);
        if !finite_ge(steps, 1.0) || steps > MAX_ADAPTIVE_SUBSTEPS as f64 {
            return Err(invalid_param(format_args!(
                "adaptive step needs {steps:.3e} explicit sub-steps (limit {MAX_ADAPTIVE_SUBSTEPS}); use step_crank_nicolson"
            )));
        }
        let num_steps = steps as usize;
        let fo_sub = self.fourier_number(dt_target / steps);

        for _ in 0..num_steps {
            self.explicit_update(fo_sub);
        }

        Ok(dt_target)
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
            return Err(invalid_param(format_args!(
                "need at least 3×3 grid, got {nx}×{ny}"
            )));
        }
        if !finite_gt(lx, 0.0) || !finite_gt(ly, 0.0) {
            return Err(invalid_param(format_args!(
                "dimensions must be positive: lx={lx}, ly={ly}"
            )));
        }
        if !t_initial.is_finite() {
            return Err(invalid_param(format_args!(
                "initial temperature {t_initial} K must be finite"
            )));
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
    ///
    /// Fixed temperatures are written to the boundary immediately; all
    /// conditions are validated by [`solve_steady_state`](Self::solve_steady_state).
    pub fn set_boundary(&mut self, side: Side, bc: BoundaryCondition) {
        match side {
            Side::Left => self.bc_left = bc,
            Side::Right => self.bc_right = bc,
            Side::Bottom => self.bc_bottom = bc,
            Side::Top => self.bc_top = bc,
        }
        if self.check_shape().is_ok() {
            self.apply_boundaries();
        }
    }

    /// Validate grid shape (fields are public and the type is deserializable).
    fn check_shape(&self) -> Result<()> {
        if self.nx < 3 || self.ny < 3 {
            return Err(invalid_param(format_args!(
                "need at least 3×3 grid, got {}×{}",
                self.nx, self.ny
            )));
        }
        if self.nodes.len() != self.ny || self.nodes.iter().any(|row| row.len() != self.nx) {
            return Err(invalid_param(format_args!(
                "node array shape does not match nx={} × ny={}",
                self.nx, self.ny
            )));
        }
        if !finite_gt(self.dx, 0.0) || !finite_gt(self.dy, 0.0) {
            return Err(invalid_param(format_args!(
                "dx={} and dy={} must be positive and finite",
                self.dx, self.dy
            )));
        }
        Ok(())
    }

    /// Solve steady-state via Gauss-Seidel iteration.
    ///
    /// Returns the number of iterations to converge (largest per-sweep change
    /// below `tolerance`). For square grids (dx = dy):
    /// T_ij = (T_{i-1,j} + T_{i+1,j} + T_{i,j-1} + T_{i,j+1}) / 4.
    ///
    /// Boundaries: `Fixed` (Dirichlet), `Insulated` (T_b = T_interior) and
    /// `Convective` (T_b = (T_interior + Bi·T_fluid)/(1 + Bi), Bi = h·Δ/k), all
    /// first-order accurate at the wall. Errors if the field becomes non-finite.
    #[tracing::instrument(level = "debug", skip(self))]
    pub fn solve_steady_state(&mut self, tolerance: f64, max_iterations: usize) -> Result<usize> {
        self.check_shape()?;
        for bc in [self.bc_left, self.bc_right, self.bc_bottom, self.bc_top] {
            bc.validate()?;
        }
        if !finite_gt(tolerance, 0.0) {
            return Err(invalid_param(format_args!(
                "tolerance {tolerance} must be positive and finite"
            )));
        }
        self.apply_boundaries();

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
                    // NaN-propagating max (f64::max would drop a NaN change).
                    let d = (new_val - old).abs();
                    if d > max_change || d.is_nan() {
                        max_change = d;
                    }
                }
            }

            // Update non-fixed boundary nodes (insulated / convective)
            let d = self.apply_flux_boundaries();
            if d > max_change || d.is_nan() {
                max_change = d;
            }

            if !max_change.is_finite() {
                return Err(invalid_param(format_args!(
                    "Gauss-Seidel produced non-finite temperatures at iteration {}",
                    iter + 1
                )));
            }
            if max_change < tolerance {
                return Ok(iter + 1);
            }
        }

        Err(invalid_param(format_args!(
            "Gauss-Seidel did not converge in {max_iterations} iterations"
        )))
    }

    /// New value of a non-fixed boundary node from its interior neighbor.
    #[inline]
    fn flux_boundary_value(bc: BoundaryCondition, t_in: f64, spacing: f64) -> Option<f64> {
        match bc {
            BoundaryCondition::Fixed(_) => None,
            BoundaryCondition::Insulated => Some(t_in),
            BoundaryCondition::Convective { t_fluid, .. } => {
                let bi = bc.biot(spacing);
                Some((t_in + bi * t_fluid) / (1.0 + bi))
            }
        }
    }

    /// Apply insulated / convective boundaries; returns the largest change.
    fn apply_flux_boundaries(&mut self) -> f64 {
        let (nx, ny) = (self.nx, self.ny);
        let mut max_change = 0.0_f64;
        let mut set = |cell: &mut f64, v: f64| {
            let d = (v - *cell).abs();
            if d > max_change || d.is_nan() {
                max_change = d;
            }
            *cell = v;
        };
        for j in 1..ny - 1 {
            if let Some(v) = Self::flux_boundary_value(self.bc_left, self.nodes[j][1], self.dx) {
                set(&mut self.nodes[j][0], v);
            }
            if let Some(v) =
                Self::flux_boundary_value(self.bc_right, self.nodes[j][nx - 2], self.dx)
            {
                set(&mut self.nodes[j][nx - 1], v);
            }
        }
        for i in 1..nx - 1 {
            if let Some(v) = Self::flux_boundary_value(self.bc_bottom, self.nodes[1][i], self.dy) {
                set(&mut self.nodes[0][i], v);
            }
            if let Some(v) = Self::flux_boundary_value(self.bc_top, self.nodes[ny - 2][i], self.dy)
            {
                set(&mut self.nodes[ny - 1][i], v);
            }
        }
        max_change
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
    #[must_use]
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
            return Err(invalid_param(format_args!(
                "node index out of range: a={node_a}, b={node_b}, num_nodes={}",
                self.num_nodes
            )));
        }
        if !finite_gt(resistance, 0.0) {
            return Err(invalid_param(format_args!(
                "resistance {resistance} K/W must be positive"
            )));
        }
        self.resistances.push((node_a, node_b, resistance));
        Ok(())
    }

    /// Set a node to a fixed temperature.
    pub fn set_fixed_temperature(&mut self, node: usize, temperature: f64) -> Result<()> {
        if node >= self.num_nodes {
            return Err(invalid_param(format_args!(
                "node {node} out of range, num_nodes={}",
                self.num_nodes
            )));
        }
        if !temperature.is_finite() {
            return Err(invalid_param(format_args!(
                "fixed temperature {temperature} K must be finite"
            )));
        }
        self.fixed.push((node, temperature));
        Ok(())
    }

    /// Validate indices and values (fields are public and deserializable).
    fn check(&self) -> Result<()> {
        let n = self.num_nodes;
        for &(a, b, r) in &self.resistances {
            if a >= n || b >= n || !finite_gt(r, 0.0) {
                return Err(invalid_param(format_args!(
                    "invalid resistance ({a}, {b}, {r}) for num_nodes={n}"
                )));
            }
        }
        for &(node, t) in &self.fixed {
            if node >= n || !t.is_finite() {
                return Err(invalid_param(format_args!(
                    "invalid fixed temperature ({node}, {t}) for num_nodes={n}"
                )));
            }
        }
        Ok(())
    }

    /// Every free node must be connected (through resistances) to a fixed
    /// node; otherwise its temperature is undetermined and the conductance
    /// matrix is singular.
    fn check_connectivity(&self, is_fixed: &[bool]) -> Result<()> {
        // Union-find over resistance edges.
        fn find(parent: &mut [usize], mut x: usize) -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        }
        let n = self.num_nodes;
        let mut parent: Vec<usize> = (0..n).collect();
        for &(a, b, _) in &self.resistances {
            let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
            if ra != rb {
                parent[ra] = rb;
            }
        }
        let mut anchored = vec![false; n];
        for (node, &fixed) in is_fixed.iter().enumerate() {
            if fixed {
                let root = find(&mut parent, node);
                anchored[root] = true;
            }
        }
        for (node, &fixed) in is_fixed.iter().enumerate() {
            if !fixed {
                let root = find(&mut parent, node);
                if !anchored[root] {
                    return Err(invalid_param(format_args!(
                        "node {node} is not connected to any fixed-temperature node; its temperature is undetermined"
                    )));
                }
            }
        }
        Ok(())
    }

    /// Solve for steady-state temperatures at all nodes.
    ///
    /// Builds a conductance matrix and solves via [`hisab::num::gaussian_elimination`].
    /// Each row is scaled by its diagonal so the pivot test is relative to the
    /// conductance scale. Errors if any free node is not connected to a fixed
    /// node (the system would be singular).
    #[tracing::instrument(level = "debug", skip(self))]
    pub fn solve(&self) -> Result<Vec<f64>> {
        self.check()?;
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
        self.check_connectivity(&is_fixed)?;

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

        // Row scaling: diagonal → 1 (positive because every free node is connected).
        for (fi, row) in matrix.iter_mut().enumerate() {
            let diag = row[fi];
            if finite_gt(diag, 0.0) {
                row.iter_mut().for_each(|v| *v /= diag);
            }
        }

        let solution = hisab::num::gaussian_elimination(&mut matrix)
            .map_err(|e| invalid_param(format_args!("network solve failed: {e}")))?;

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

    // --- Regression tests from the 1.4.0 audit ---

    fn convective(h: f64, k: f64, t_fluid: f64) -> BoundaryCondition {
        BoundaryCondition::Convective { h, t_fluid, k }
    }

    #[test]
    fn test_explicit_stability_accounts_for_biot() {
        // n=11, L=0.1 → dx=0.01; Bi = 1000·0.01/1 = 10 → Fo_max = 0.5/11.
        let bc = convective(1000.0, 1.0, 300.0);
        let mut g = ThermalGrid1D::new(11, 0.1, 1e-5, 400.0, bc, bc).unwrap();
        // Fo = 0.4 passed the old Fo ≤ 0.5 check and blew up (T0 → −400 K).
        assert!(g.step_explicit(4.0).is_err());
        let dt = g.max_stable_dt();
        assert!((g.fourier_number(dt) - 0.5 / 11.0).abs() < 1e-12);
        for _ in 0..200 {
            g.step_explicit(dt).unwrap();
        }
        assert!(
            g.nodes.iter().all(|&t| (300.0..=400.0).contains(&t)),
            "{:?}",
            g.nodes
        );
    }

    #[test]
    fn test_step_adaptive_exact_multiple_of_stable_dt() {
        // Used to fail with Fo = 0.5000000000000001 from round-off.
        let bc = BoundaryCondition::Fixed(300.0);
        let mut g = ThermalGrid1D::new(51, 0.1, 1e-5, 400.0, bc, bc).unwrap();
        assert_eq!(g.step_adaptive(1.0).unwrap(), 1.0);
        let dt = g.max_stable_dt();
        assert!(g.step_adaptive(dt).is_ok());
        assert!(g.step_adaptive(3.0 * dt).is_ok());
    }

    #[test]
    fn test_step_adaptive_substep_cap_and_bad_dt() {
        let bc = BoundaryCondition::Insulated;
        let mut g = ThermalGrid1D::new(3, 2e-6, 1e-4, 300.0, bc, bc).unwrap();
        assert!(g.step_adaptive(3600.0).is_err()); // ~7e11 sub-steps
        assert!(g.step_adaptive(f64::NAN).is_err());
        assert!(g.step_adaptive(-1.0).is_err());
        assert!(g.step_explicit(f64::NAN).is_err());
        assert!(g.step_crank_nicolson(-1.0).is_err());
        assert!(g.step_crank_nicolson(f64::INFINITY).is_err());
    }

    #[test]
    fn test_explicit_matches_crank_nicolson_for_small_steps() {
        let bc_l = BoundaryCondition::Fixed(400.0);
        let bc_r = convective(50.0, 15.0, 290.0);
        let mut ex = ThermalGrid1D::new(21, 0.2, 4e-6, 300.0, bc_l, bc_r).unwrap();
        let mut cn = ex.clone();
        let dt = ex.max_stable_dt() * 0.5;
        for _ in 0..2000 {
            ex.step_explicit(dt).unwrap();
            cn.step_crank_nicolson(dt).unwrap();
        }
        for (a, b) in ex.nodes.iter().zip(&cn.nodes) {
            assert!((a - b).abs() < 0.05, "{a} vs {b}");
        }
    }

    #[test]
    fn test_grid_1d_bad_state_is_error_not_panic() {
        let bc = BoundaryCondition::Insulated;
        let mut g = ThermalGrid1D::new(5, 1.0, 1e-5, 300.0, bc, bc).unwrap();
        g.nodes.truncate(1);
        assert!(g.step_explicit(1.0).is_err());
        assert!(g.step_crank_nicolson(1.0).is_err());
        assert!(g.step_adaptive(1.0).is_err());
        g.nodes.clear();
        assert!(g.step_explicit(1.0).is_err());
        assert!(ThermalGrid1D::new(5, 1.0, 1e-5, f64::NAN, bc, bc).is_err());
        assert!(ThermalGrid1D::new(5, 1.0, 1e-5, 300.0, convective(10.0, 0.0, 300.0), bc).is_err());
    }

    #[test]
    fn test_grid_2d_convective_boundary_matches_analytic() {
        // Left fixed 400 K, right convective to 300 K, top/bottom insulated:
        // a 1D slab with q = (T1 − T∞)/(L/k + 1/h). The discrete solution is exact
        // for this linear profile.
        let (l, k, h) = (0.5, 2.0, 10.0);
        let mut g = ThermalGrid2D::new(11, 5, l, 0.2, 350.0).unwrap();
        g.set_boundary(Side::Left, BoundaryCondition::Fixed(400.0));
        g.set_boundary(Side::Right, convective(h, k, 300.0));
        g.solve_steady_state(1e-12, 200_000).unwrap();
        let t_wall = 300.0 + 100.0 * (1.0 / h) / (l / k + 1.0 / h);
        let got = g.nodes[2][10];
        assert!((got - t_wall).abs() < 1e-6, "{got} vs {t_wall}");
        // The convective side is no longer pinned at the initial temperature.
        assert!((got - 350.0).abs() > 1.0);
    }

    #[test]
    fn test_grid_2d_nan_and_bad_shape() {
        assert!(ThermalGrid2D::new(5, 5, 1.0, 1.0, f64::NAN).is_err());
        let mut g = ThermalGrid2D::new(5, 5, 1.0, 1.0, 300.0).unwrap();
        g.bc_left = BoundaryCondition::Fixed(f64::NAN);
        assert!(g.solve_steady_state(1e-6, 100).is_err());
        let mut g = ThermalGrid2D::new(5, 5, 1.0, 1.0, 300.0).unwrap();
        g.nodes[2][2] = f64::NAN;
        assert!(g.solve_steady_state(1e-6, 100).is_err());
        let mut g = ThermalGrid2D::new(5, 5, 1.0, 1.0, 300.0).unwrap();
        g.nodes.pop();
        assert!(g.solve_steady_state(1e-6, 100).is_err());
        g.set_boundary(Side::Left, BoundaryCondition::Fixed(400.0)); // must not panic
        assert!(
            ThermalGrid2D::new(5, 5, 1.0, 1.0, 300.0)
                .unwrap()
                .solve_steady_state(f64::NAN, 10)
                .is_err()
        );
    }

    #[test]
    fn test_network_empty_and_out_of_range() {
        let mut net = ThermalNetwork::new(0);
        assert!(net.add_resistance(0, 0, 1.0).is_err()); // used to underflow num_nodes − 1
        assert!(net.set_fixed_temperature(0, 300.0).is_err());
        let mut net = ThermalNetwork::new(2);
        assert!(net.set_fixed_temperature(0, f64::NAN).is_err());
        // Deserialized garbage must not panic.
        net.resistances.push((0, 5, 1.0));
        assert!(net.solve().is_err());
    }

    #[test]
    fn test_network_unanchored_nodes_are_errors() {
        // Triangle of free nodes with large conductances and no fixed node:
        // used to "solve" to all zeros.
        let mut net = ThermalNetwork::new(3);
        net.add_resistance(0, 1, 1e-3).unwrap();
        net.add_resistance(1, 2, 1e-3).unwrap();
        net.add_resistance(0, 2, 3e-4).unwrap();
        assert!(net.solve().is_err());
        // A floating island next to an anchored part.
        let mut net = ThermalNetwork::new(4);
        net.add_resistance(0, 1, 1.0).unwrap();
        net.add_resistance(2, 3, 1.0).unwrap();
        net.set_fixed_temperature(0, 300.0).unwrap();
        assert!(net.solve().is_err());
    }

    #[test]
    fn test_network_tiny_conductances_solve() {
        // R = 1e14 K/W → g = 1e-14, below hisab's absolute pivot threshold
        // without row scaling.
        let mut net = ThermalNetwork::new(3);
        net.add_resistance(0, 1, 1e14).unwrap();
        net.add_resistance(1, 2, 1e14).unwrap();
        net.set_fixed_temperature(0, 400.0).unwrap();
        net.set_fixed_temperature(2, 300.0).unwrap();
        let t = net.solve().unwrap();
        assert!((t[1] - 350.0).abs() < 1e-9);
    }
}
