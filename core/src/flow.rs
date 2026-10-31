//! Cone-preserving discrete (vertex-scaling) Ricci flow: port of
//! `embeddings/true_flow.py::DiscreteConeFlow` (seeds smart / gluing / polytope,
//! damped Newton with the exact sparse Jacobian + Sherman-Morrison rank-one rho
//! term, adaptive variable-step BDF2 with a backward-Euler first step, area
//! gauge, diagnostics, frames).  rho is an OUTPUT (Gauss-Bonnet).
use crate::mesh::{area_jacobian, cot_laplacian, defects_and_areas, Mesh, VertexGeom};
use crate::sparse::{Coo, SparseLu};
use std::f64::consts::PI;

const TWO_PI: f64 = 2.0 * PI;

/// Starting metric of the flow.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Seed {
    /// round metric * prod |z - p_i|^{2(beta_i - 1)} (flat-cone far field)
    Smart,
    /// round sphere with cone caps glued in
    Gluing,
    /// flat faces, curvature in cones + 4 extra corners
    Polytope,
}

impl Seed {
    pub fn parse(s: &str) -> Option<Seed> {
        match s.to_ascii_lowercase().as_str() {
            "smart" => Some(Seed::Smart),
            "gluing" => Some(Seed::Gluing),
            "polytope" => Some(Seed::Polytope),
            _ => None,
        }
    }
}

/// Solver settings (defaults = Python defaults).
#[derive(Clone, Debug)]
pub struct FlowOptions {
    pub times: Vec<f64>,
    pub newton_tol: f64,
    pub max_newton: usize,
    pub lte_tol: f64,
    pub dt0: f64,
    pub dt_max: f64,
    pub seed: Seed,
}

impl Default for FlowOptions {
    fn default() -> Self {
        FlowOptions {
            times: default_times(16.0, 12, 8, 0.05, 1e-3),
            newton_tol: 1e-11,
            max_newton: 30,
            lte_tol: 2e-3,
            dt0: 1e-4,
            dt_max: 2.0,
            seed: Seed::Smart,
        }
    }
}

fn geomspace(a: f64, b: f64, n: usize) -> Vec<f64> {
    // numpy.geomspace: 10**linspace(log10 a, log10 b), endpoints exact
    if n == 1 {
        return vec![a];
    }
    let (la, lb) = (a.log10(), b.log10());
    let mut v: Vec<f64> = (0..n).map(|i| 10f64.powf(la + (lb - la) * i as f64 / (n - 1) as f64)).collect();
    v[0] = a;
    v[n - 1] = b;
    v
}

/// `true_flow.default_times`.
pub fn default_times(t_end: f64, n_early: usize, n_late: usize, t_split: f64, t_first: f64) -> Vec<f64> {
    let mut t = vec![0.0];
    t.extend(geomspace(t_first, t_split, n_early));
    t.extend(geomspace(t_split * 1.4, t_end, n_late));
    t
}

/// Solver counters (same keys as Python `stats`).
#[derive(Clone, Debug, Default)]
pub struct FlowStats {
    pub newton: usize,
    pub steps: usize,
    pub rejected: usize,
    pub newton_fail: usize,
}

/// Per-frame diagnostics (`DiscreteConeFlow.diagnostics`).
#[derive(Clone, Debug)]
pub struct Diagnostics {
    pub valid: bool,
    pub gb_resid: f64,
    pub cone_defect_err: f64,
    pub rho: f64,
    pub area: f64,
    pub k_smooth_min: f64,
    pub n_k_nonpos: usize,
    pub n_defect_nonpos_all: usize,
    pub kdens_rel_dev_max: f64,
    pub kdens_rel_dev_rms: f64,
}

/// One intrinsic metric frame.
#[derive(Clone, Debug)]
pub struct Frame {
    pub k: usize,
    pub t: f64,
    /// log conformal factor per vertex
    pub u: Vec<f64>,
    /// edge lengths (mesh.e order)
    pub l: Vec<f64>,
    /// angle defects per vertex
    pub k_defect: Vec<f64>,
    /// barycentric vertex areas (face area / 3)
    pub a_vertex: Vec<f64>,
    pub diag: Diagnostics,
}

/// BDF2 history: (u_prev, dt_prev, dt_next).
#[derive(Clone, Debug)]
pub struct History {
    pub u_prev: Vec<f64>,
    pub dt_prev: f64,
    pub dt_next: f64,
}

pub type FlowResult<T> = Result<T, String>;

struct Eval {
    r: Vec<f64>,
    l: Vec<f64>,
    g: VertexGeom,
    rho: f64,
}

/// Cone-preserving discrete Ricci flow.
pub struct DiscreteConeFlow {
    pub mesh: Mesh,
    pub betas: Vec<f64>,
    pub opts: FlowOptions,
    pub l0: Vec<f64>,
    pub smooth: Vec<bool>,
    pub cone_k: Vec<f64>,
    /// total smooth curvature 4 pi - sum cone defects
    pub c: f64,
    pub stats: FlowStats,
    /// polytope seed: extra corner vertex ids
    pub seed_corners: Vec<usize>,
    /// gluing seed: cap radius
    pub gluing_radius: Option<f64>,
}

fn maxabs(v: &[f64]) -> f64 {
    v.iter().fold(0.0f64, |m, x| m.max(x.abs()))
}
fn vdot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let w = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt()
}
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

impl DiscreteConeFlow {
    pub fn new(mesh: Mesh, betas: Vec<f64>, opts: FlowOptions) -> Self {
        let nv = mesh.nv();
        let l0 = mesh.chord_lengths();
        let mut smooth = vec![true; nv];
        for &c in &mesh.cones {
            smooth[c] = false;
        }
        let cone_k: Vec<f64> = betas.iter().map(|b| TWO_PI * (1.0 - b)).collect();
        let c = 4.0 * PI - cone_k.iter().sum::<f64>();
        DiscreteConeFlow {
            mesh,
            betas,
            opts,
            l0,
            smooth,
            cone_k,
            c,
            stats: FlowStats::default(),
            seed_corners: vec![],
            gluing_radius: None,
        }
    }

    pub fn nv(&self) -> usize {
        self.mesh.nv()
    }

    /// l_ij = exp((u_i + u_j)/2) l0_ij
    pub fn lengths(&self, u: &[f64]) -> Vec<f64> {
        self.mesh.e.iter().zip(&self.l0).map(|(&[i, j], l0)| (0.5 * (u[i] + u[j])).exp() * l0).collect()
    }

    pub fn geom(&self, u: &[f64]) -> (Vec<f64>, VertexGeom) {
        let l = self.lengths(u);
        let g = defects_and_areas(&self.mesh.f, &self.mesh.opp, &l, self.nv());
        (l, g)
    }

    /// Shift u so the total area is 4 pi (K and rho A_v unchanged).
    pub fn gauge(&self, u: &[f64]) -> Vec<f64> {
        let (_, g) = self.geom(u);
        let s = 0.5 * (4.0 * PI / g.av.iter().sum::<f64>()).ln();
        u.iter().map(|x| x + s).collect()
    }

    fn smooth_area(&self, av: &[f64]) -> f64 {
        av.iter().zip(&self.smooth).filter(|(_, &s)| s).map(|(a, _)| a).sum()
    }

    fn residual(&self, u: &[f64], hist: &[f64], mn: &[f64], dt: f64, c0: f64) -> Eval {
        let (l, g) = self.geom(u);
        let rho = self.c / self.smooth_area(&g.av);
        let mut r: Vec<f64> =
            (0..u.len()).map(|v| mn[v] * (c0 * u[v] - hist[v]) - dt * (rho * g.av[v] - g.k[v])).collect();
        for (ci, &c) in self.mesh.cones.iter().enumerate() {
            r[c] = g.k[c] - self.cone_k[ci];
        }
        Eval { r, l, g, rho }
    }

    /// Damped Newton for  Mn (c0 u - hist) - dt (rho A - K) = 0 (smooth rows),
    /// K_c = target (cone rows).  Returns (u, ok, iterations).
    pub fn newton(&mut self, u0: &[f64], hist: &[f64], mn: &[f64], dt: f64, c0: f64) -> (Vec<f64>, bool, usize) {
        let nv = self.nv();
        let mut u = u0.to_vec();
        let mut ev = self.residual(&u, hist, mn, dt, c0);
        if !ev.g.face.valid {
            return (u, false, 0);
        }
        let s: Vec<f64> = self.smooth.iter().map(|&b| if b { 1.0 } else { 0.0 }).collect();
        for it in 0..self.opts.max_newton {
            let rn = maxabs(&ev.r);
            if rn < self.opts.newton_tol {
                return (u, true, it);
            }
            let rho = ev.rho;
            let a_s = self.smooth_area(&ev.g.av);
            let mut h = Coo::new(nv);
            cot_laplacian(&self.mesh.f, &ev.g.face.cot, &mut h);
            let mut da = Coo::new(nv);
            area_jacobian(&self.mesh.f, &self.mesh.opp, &ev.l, &ev.g.face.cot, &mut da);
            let q: Vec<f64> = da.matvec_t(&s).iter().map(|x| -(rho / a_s) * x).collect();
            let mut m = Coo::new(nv);
            for v in 0..nv {
                if self.smooth[v] {
                    m.push(v, v, c0 * mn[v]);
                }
            }
            for k in 0..h.vals.len() {
                let i = h.rows[k];
                m.push(i, h.cols[k], if self.smooth[i] { dt * h.vals[k] } else { h.vals[k] });
            }
            for k in 0..da.vals.len() {
                let i = da.rows[k];
                if self.smooth[i] {
                    m.push(i, da.cols[k], -dt * rho * da.vals[k]);
                }
            }
            let p: Vec<f64> = (0..nv).map(|v| -dt * s[v] * ev.g.av[v]).collect();
            let lu = match SparseLu::factor(&m) {
                Ok(lu) => lu,
                Err(_) => return (u, false, it),
            };
            let negr: Vec<f64> = ev.r.iter().map(|x| -x).collect();
            let sol = lu.solve_many(&[&negr, &p]);
            let (x, y) = (&sol[0], &sol[1]);
            let fac = vdot(&q, x) / (1.0 + vdot(&q, y));
            let du: Vec<f64> = x.iter().zip(y).map(|(a, b)| a - b * fac).collect();
            let mut step = 1.0;
            let mut accepted = None;
            while step > 1e-6 {
                let un: Vec<f64> = u.iter().zip(&du).map(|(a, d)| a + step * d).collect();
                let evn = self.residual(&un, hist, mn, dt, c0);
                if evn.g.face.valid && maxabs(&evn.r) < (1.0 - 1e-4 * step) * rn {
                    accepted = Some((un, evn));
                    break;
                }
                step *= 0.5;
            }
            match accepted {
                None => return (u, false, it),
                Some((un, evn)) => {
                    u = un;
                    ev = evn;
                }
            }
            self.stats.newton += 1;
        }
        let ok = maxabs(&ev.r) < 100.0 * self.opts.newton_tol;
        (u, ok, self.opts.max_newton)
    }

    /// Cone vertex starts at the mean of its smooth neighbours (edge-list order).
    fn cone_neighbour_mean(&self, u: &mut [f64]) {
        for &c in &self.mesh.cones {
            let mut nb: Vec<usize> = self.mesh.e.iter().filter(|e| e[0] == c).map(|e| e[1]).collect();
            nb.extend(self.mesh.e.iter().filter(|e| e[1] == c).map(|e| e[0]));
            let nb: Vec<usize> = nb.into_iter().filter(|&j| self.smooth[j]).collect();
            u[c] = nb.iter().map(|&j| u[j]).sum::<f64>() / nb.len() as f64;
        }
    }

    /// `far_field_u`: sum_c (beta_c - 1) log|Y - Y_c| at smooth vertices; cones at neighbour mean.
    pub fn far_field_u(&self) -> Vec<f64> {
        let y = &self.mesh.y;
        let nv = self.nv();
        let mut u = vec![0.0; nv];
        for (ci, &c) in self.mesh.cones.iter().enumerate() {
            let b = self.betas[ci];
            for v in 0..nv {
                let d = if self.smooth[v] { dist3(y[v], y[c]) } else { 1.0 };
                u[v] += (b - 1.0) * d.ln();
            }
        }
        self.cone_neighbour_mean(&mut u);
        u
    }

    /// Starting log-scale u(0) for the chosen seed.
    pub fn initial_u(&mut self) -> FlowResult<Vec<f64>> {
        match self.opts.seed {
            Seed::Gluing => self.gluing_u(),
            Seed::Polytope => self.polytope_u(4),
            Seed::Smart => self.smart_u(),
        }
    }

    fn algebraic_cone_solve(&mut self, u: Vec<f64>, msg: &str) -> FlowResult<Vec<f64>> {
        let ones = vec![1.0; self.nv()];
        let hist = u.clone();
        let (u, ok, _) = self.newton(&u, &hist, &ones, 0.0, 1.0);
        if !ok {
            return Err(msg.to_string());
        }
        Ok(self.gauge(&u))
    }

    /// Far field + algebraic solve of the cone constraints (smooth u fixed).
    pub fn smart_u(&mut self) -> FlowResult<Vec<f64>> {
        let mut u = self.far_field_u();
        let (sum, n) = u
            .iter()
            .zip(&self.smooth)
            .filter(|(_, &s)| s)
            .fold((0.0, 0usize), |(a, n), (x, _)| (a + x, n + 1));
        let mean = sum / n as f64;
        for x in u.iter_mut() {
            *x -= mean;
        }
        self.algebraic_cone_solve(u, "cone constraints at t=0 failed (invalid far-field triangles?)")
    }

    /// Round sphere with cone caps of radius (min cone distance)/4 glued in.
    pub fn gluing_u(&mut self) -> FlowResult<Vec<f64>> {
        let y = self.mesh.y.clone();
        let cones = self.mesh.cones.clone();
        let r = if cones.len() > 1 {
            let mut m = f64::INFINITY;
            for a in 0..cones.len() {
                for b in a + 1..cones.len() {
                    m = m.min(dot3(y[cones[a]], y[cones[b]]).clamp(-1.0, 1.0).acos());
                }
            }
            m / 4.0
        } else {
            PI / 4.0
        };
        self.gluing_radius = Some(r);
        let nv = self.nv();
        let mut u = vec![0.0; nv];
        for (ci, &c) in cones.iter().enumerate() {
            let b = self.betas[ci];
            for v in 0..nv {
                let d = dot3(y[v], y[c]).clamp(-1.0, 1.0).acos();
                if d < r && self.smooth[v] {
                    u[v] += (b - 1.0) * (d.max(1e-12) / r).ln();
                }
            }
        }
        self.cone_neighbour_mean(&mut u);
        self.algebraic_cone_solve(u, "cone gluing seed failed (cone fan triangles became invalid)")
    }

    /// Flat metric with all curvature in corners (cones + `ncorner` farthest-point corners).
    pub fn polytope_u(&mut self, ncorner: usize) -> FlowResult<Vec<f64>> {
        let y = self.mesh.y.clone();
        let nv = self.nv();
        let mut d: Vec<f64> = (0..nv)
            .map(|v| self.mesh.cones.iter().map(|&c| dist3(y[v], y[c])).fold(f64::INFINITY, f64::min))
            .collect();
        let mut corners = Vec::new();
        for _ in 0..ncorner {
            let mut j = 0;
            for v in 1..nv {
                if d[v] > d[j] {
                    j = v;
                }
            }
            corners.push(j);
            for v in 0..nv {
                d[v] = d[v].min(dist3(y[v], y[j]));
            }
        }
        let mut tgt = vec![0.0; nv];
        for (ci, &c) in self.mesh.cones.iter().enumerate() {
            tgt[c] = self.cone_k[ci];
        }
        for &j in &corners {
            tgt[j] = self.c / ncorner as f64;
        }
        let mut u = self.smart_u()?;
        let mut converged = false;
        for _ in 0..80 {
            let (_, g) = self.geom(&u);
            let r: Vec<f64> = (0..nv).map(|v| g.k[v] - tgt[v]).collect();
            let rn = maxabs(&r);
            if rn < 1e-10 {
                converged = true;
                break;
            }
            let mut h = Coo::new(nv);
            cot_laplacian(&self.mesh.f, &g.face.cot, &mut h);
            for v in 0..nv {
                h.push(v, v, 1e-10); // K is shift-invariant
            }
            let lu = SparseLu::factor(&h)?;
            let negr: Vec<f64> = r.iter().map(|x| -x).collect();
            let du = lu.solve(&negr);
            let mut step = 1.0;
            let mut acc = None;
            while step > 1e-8 {
                let un: Vec<f64> = u.iter().zip(&du).map(|(a, b)| a + step * b).collect();
                let (_, gn) = self.geom(&un);
                let rnn = (0..nv).map(|v| (gn.k[v] - tgt[v]).abs()).fold(0.0, f64::max);
                if gn.face.valid && rnn < (1.0 - 1e-4 * step) * rn {
                    acc = Some(un);
                    break;
                }
                step *= 0.5;
            }
            match acc {
                Some(un) => u = un,
                None => return Err(format!("polytope seed failed (Newton stalled at residual {rn:.2e})")),
            }
        }
        if !converged {
            return Err("polytope seed failed (no convergence)".into());
        }
        self.seed_corners = corners;
        Ok(self.gauge(&u))
    }

    /// One backward-Euler step (comparison / tests).
    pub fn step(&mut self, u: &[f64], dt: f64) -> FlowResult<Vec<f64>> {
        let (_, g) = self.geom(u);
        let (un, ok, _) = self.newton(u, u, &g.av, dt, 1.0);
        if !ok {
            return Err("implicit step failed".into());
        }
        Ok(self.gauge(&un))
    }

    fn centred(&self, u: &[f64]) -> Vec<f64> {
        let v: Vec<f64> = u.iter().zip(&self.smooth).filter(|(_, &s)| s).map(|(x, _)| *x).collect();
        let m = v.iter().sum::<f64>() / v.len() as f64;
        v.into_iter().map(|x| x - m).collect()
    }

    /// Adaptive variable-step BDF2 from (t0, u) to exactly t1 (first step backward Euler
    /// when `hist` is None).  Returns (u, hist).
    pub fn integrate(
        &mut self,
        u0: &[f64],
        t0: f64,
        t1: f64,
        hist: Option<History>,
    ) -> FlowResult<(Vec<f64>, Option<History>)> {
        let mut u = u0.to_vec();
        let mut hist = hist;
        let mut t = t0;
        let mut dt = match &hist {
            None => self.opts.dt0,
            Some(h) => h.dt_next,
        };
        while t < t1 - 1e-14 * t1.max(1.0) {
            let h = dt.min(t1 - t).min(self.opts.dt_max);
            let (_, g) = self.geom(&u);
            let mut err = 0.0;
            let (un, ok) = match &hist {
                None => {
                    let uc = u.clone();
                    let (un, ok, _) = self.newton(&u, &uc, &g.av, h, 1.0);
                    (un, ok)
                }
                Some(hs) => {
                    let (up, hp) = (hs.u_prev.clone(), hs.dt_prev);
                    let w = h / hp;
                    let (c0, c1, c2) = ((1.0 + 2.0 * w) / (1.0 + w), 1.0 + w, w * w / (1.0 + w));
                    let hv: Vec<f64> = u.iter().zip(&up).map(|(a, b)| c1 * a - c2 * b).collect();
                    let (un, ok, _) = self.newton(&u, &hv, &g.av, h, c0);
                    if ok {
                        let pred: Vec<f64> = u.iter().zip(&up).map(|(a, b)| a + w * (a - b)).collect();
                        let (cu, cp) = (self.centred(&un), self.centred(&pred));
                        let dmax = cu.iter().zip(&cp).fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
                        err = dmax * w / (1.0 + 2.0 * w);
                    }
                    (un, ok)
                }
            };
            if !ok {
                self.stats.newton_fail += 1;
                dt = h / 4.0;
                if dt < 1e-12 {
                    return Err("implicit step failed".into());
                }
                continue;
            }
            if err > 2.0 * self.opts.lte_tol {
                self.stats.rejected += 1;
                dt = h * 0.25f64.max(0.9 * (self.opts.lte_tol / err).sqrt());
                continue;
            }
            let un = self.gauge(&un);
            let fac = if err == 0.0 { 2.0 } else { 2.0f64.min(0.9 * (self.opts.lte_tol / err).sqrt()) };
            let mut dt_next = (h * fac).max(1e-12);
            if h < dt {
                dt_next = dt_next.max(dt); // step clipped to an output time: keep dt
            }
            hist = Some(History { u_prev: u, dt_prev: h, dt_next });
            u = un;
            t += h;
            dt = dt_next;
            self.stats.steps += 1;
        }
        Ok((u, hist))
    }

    pub fn diagnostics_of(&self, g: &VertexGeom) -> Diagnostics {
        let rho = self.c / self.smooth_area(&g.av);
        let sm: Vec<usize> = (0..self.nv()).filter(|&v| self.smooth[v]).collect();
        let kg: Vec<f64> = sm.iter().map(|&v| g.k[v] / g.av[v] / rho - 1.0).collect();
        Diagnostics {
            valid: g.face.valid,
            gb_resid: g.k.iter().sum::<f64>() - 4.0 * PI,
            cone_defect_err: self
                .mesh
                .cones
                .iter()
                .enumerate()
                .map(|(ci, &c)| (g.k[c] - self.cone_k[ci]).abs())
                .fold(0.0, f64::max),
            rho,
            area: g.av.iter().sum(),
            k_smooth_min: sm.iter().map(|&v| g.k[v]).fold(f64::INFINITY, f64::min),
            n_k_nonpos: sm.iter().filter(|&&v| g.k[v] <= 0.0).count(),
            n_defect_nonpos_all: g.k.iter().filter(|&&k| k <= 0.0).count(),
            kdens_rel_dev_max: kg.iter().fold(0.0, |m: f64, x| m.max(x.abs())),
            kdens_rel_dev_rms: (kg.iter().map(|x| x * x).sum::<f64>() / kg.len() as f64).sqrt(),
        }
    }

    pub fn diagnostics(&self, u: &[f64]) -> Diagnostics {
        let (_, g) = self.geom(u);
        self.diagnostics_of(&g)
    }

    pub fn make_frame(&self, k: usize, t: f64, u: &[f64]) -> Frame {
        let (l, g) = self.geom(u);
        let diag = self.diagnostics_of(&g);
        Frame { k, t, u: u.to_vec(), l, k_defect: g.k, a_vertex: g.av, diag }
    }

    /// Metric at an intermediate time integrated from `frame` (`DiscreteConeFlow.between`).
    pub fn between(&mut self, frame: &Frame, t_mid: f64) -> FlowResult<Frame> {
        let (u, _) = self.integrate(&frame.u, frame.t, t_mid, None)?;
        Ok(self.make_frame(frame.k, t_mid, &u))
    }

    /// Start a streaming run (computes the seed); call `FlowRun::next_frame` repeatedly.
    pub fn start(mut self) -> FlowResult<FlowRun> {
        let u = self.initial_u()?;
        Ok(FlowRun { flow: self, u, tcur: 0.0, hist: None, k: 0 })
    }

    /// All frames at `opts.times` (`DiscreteConeFlow.frames`); returns the solver too
    /// (stats, seed_corners).
    pub fn frames(self) -> FlowResult<(Vec<Frame>, DiscreteConeFlow)> {
        let mut run = self.start()?;
        let mut out = Vec::new();
        while let Some(fr) = run.next_frame()? {
            out.push(fr);
        }
        Ok((out, run.flow))
    }
}

/// Streaming frame iterator (one frame per call, for incremental UI updates).
pub struct FlowRun {
    pub flow: DiscreteConeFlow,
    pub u: Vec<f64>,
    pub tcur: f64,
    pub hist: Option<History>,
    pub k: usize,
}

impl FlowRun {
    pub fn n_frames(&self) -> usize {
        self.flow.opts.times.len()
    }
    pub fn next_frame(&mut self) -> FlowResult<Option<Frame>> {
        if self.k >= self.flow.opts.times.len() {
            return Ok(None);
        }
        let tk = self.flow.opts.times[self.k];
        if tk > self.tcur {
            let (u, h) = self.flow.integrate(&self.u, self.tcur, tk, self.hist.take())?;
            self.u = u;
            self.hist = h;
            self.tcur = tk;
        }
        let fr = self.flow.make_frame(self.k, tk, &self.u);
        self.k += 1;
        Ok(Some(fr))
    }
}

// ------------------------------------------------------------------ wasm binding
use wasm_bindgen::prelude::*;

/// JS-facing streaming flow.  `y` flat (nv*3) unit-sphere points, `f` flat (nf*3)
/// outward faces, `cones` vertex ids, `betas` per cone, `seed` "smart"|"gluing"|"polytope".
/// Uses the default output times.  Each `next_frame` returns false when done; the
/// current frame is then read through the getters.
#[wasm_bindgen]
pub struct WasmFlow {
    run: FlowRun,
    cur: Option<Frame>,
}

#[wasm_bindgen]
impl WasmFlow {
    #[wasm_bindgen(constructor)]
    pub fn new(y: &[f64], f: &[u32], cones: &[u32], betas: &[f64], seed: &str) -> Result<WasmFlow, JsValue> {
        let yv: Vec<[f64; 3]> = y.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let fv: Vec<[usize; 3]> = f.chunks(3).map(|c| [c[0] as usize, c[1] as usize, c[2] as usize]).collect();
        let mesh = Mesh::from_parts(yv, fv, cones.iter().map(|&c| c as usize).collect());
        let seed = Seed::parse(seed).ok_or_else(|| JsValue::from_str("unknown seed"))?;
        let flow = DiscreteConeFlow::new(mesh, betas.to_vec(), FlowOptions { seed, ..Default::default() });
        let run = flow.start().map_err(|e| JsValue::from_str(&e))?;
        Ok(WasmFlow { run, cur: None })
    }
    pub fn n_frames(&self) -> usize {
        self.run.n_frames()
    }
    pub fn next_frame(&mut self) -> Result<bool, JsValue> {
        self.cur = self.run.next_frame().map_err(|e| JsValue::from_str(&e))?;
        Ok(self.cur.is_some())
    }
    pub fn t(&self) -> f64 {
        self.cur.as_ref().map_or(f64::NAN, |f| f.t)
    }
    pub fn rho(&self) -> f64 {
        self.cur.as_ref().map_or(f64::NAN, |f| f.diag.rho)
    }
    pub fn gb_resid(&self) -> f64 {
        self.cur.as_ref().map_or(f64::NAN, |f| f.diag.gb_resid)
    }
    pub fn u(&self) -> Vec<f64> {
        self.cur.as_ref().map_or(vec![], |f| f.u.clone())
    }
    pub fn lengths(&self) -> Vec<f64> {
        self.cur.as_ref().map_or(vec![], |f| f.l.clone())
    }
    pub fn k_defect(&self) -> Vec<f64> {
        self.cur.as_ref().map_or(vec![], |f| f.k_defect.clone())
    }
    pub fn a_vertex(&self) -> Vec<f64> {
        self.cur.as_ref().map_or(vec![], |f| f.a_vertex.clone())
    }
    /// flat (ne*2) edge list matching `lengths`
    pub fn edges(&self) -> Vec<u32> {
        self.run.flow.mesh.e.iter().flat_map(|e| [e[0] as u32, e[1] as u32]).collect()
    }
    pub fn seed_corners(&self) -> Vec<u32> {
        self.run.flow.seed_corners.iter().map(|&c| c as u32).collect()
    }
}
