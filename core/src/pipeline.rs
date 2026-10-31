//! End-to-end /flow "True Ricci Flow" pipeline (rho > 0): port of
//! `viewer/backend/flow_true.run_flow_true` (Fast realizer, backward realization) as a steppable
//! [`Session`]: mesh from cone positions -> intrinsic flow frames -> clean start -> backward 3D
//! realization -> Polytope frame 0 + seed-fold exemption + quality note -> per-frame rows with the
//! conformal-factor raster.  rho <= 0 is refused with an error (stage 2).

use crate::cfgrid::{CfGrid, GRID_LAT, GRID_LON};
use crate::flow::{default_times, DiscreteConeFlow, FlowOptions, FlowRun, FlowStats, Frame, Seed};
use crate::mesh::Mesh;
use crate::quality::{flow_quality, mark_seed_folds, FrameRow, Quality};
use crate::realize::{self, Backward, FastRealizer, FlowSource, MetricFrame, RealizeMesh, Realized, SeedKind};
use crate::realize::geom::V3;
use crate::seedmesh::{self, FlowMesh};
use std::cell::RefCell;
use std::f64::consts::PI;

pub const NEG_REFUSAL: &str = "True Ricci Flow needs rho > 0 (for rho < 0 intermediate embeddings generally do not exist; the rho < 0 3D path is not ported yet)";

/// Session inputs.
#[derive(Clone, Debug)]
pub struct SessionConfig {
    /// cone positions in the stereographic plane (x, y) = (Re z, Im z)
    pub positions: Vec<[f64; 2]>,
    pub betas: Vec<f64>,
    pub seed: Seed,
    /// (N Fibonacci points, m spokes, K rings); see [`seedmesh::page_level`]
    pub level: (usize, usize, usize),
}

/// The flow as a [`FlowSource`] for the realizer (between = `solver.between`, blend of u(0), u(T)).
pub struct FlowSrc {
    pub flow: RefCell<DiscreteConeFlow>,
    pub frames: Vec<Frame>,
    u_t: RefCell<Option<Vec<f64>>>,
}

impl FlowSrc {
    pub fn new(flow: DiscreteConeFlow, frames: Vec<Frame>) -> Self {
        FlowSrc { flow: RefCell::new(flow), frames, u_t: RefCell::new(None) }
    }
    /// u(T) integrated from u(0) in one call, T = max(times, default end) (`clean_start`).
    fn u_end(&self) -> Vec<f64> {
        if let Some(u) = self.u_t.borrow().as_ref() {
            return u.clone();
        }
        let mut fl = self.flow.borrow_mut();
        let dt = default_times(16.0, 12, 8, 0.05, 1e-3);
        let t_end = fl.opts.times.iter().chain(&dt).cloned().fold(f64::NEG_INFINITY, f64::max);
        let u0 = self.frames[0].u.clone();
        let u = match fl.integrate(&u0, 0.0, t_end, None) {
            Ok((u, _)) => u,
            Err(_) => self.frames.last().unwrap().u.clone(),
        };
        *self.u_t.borrow_mut() = Some(u.clone());
        u
    }
}

impl FlowSource for FlowSrc {
    fn between(&self, frame: &MetricFrame, t_mid: f64) -> Vec<f64> {
        let fr = &self.frames[frame.k];
        match self.flow.borrow_mut().between(fr, t_mid) {
            Ok(f) => f.l,
            Err(_) => frame.lengths.clone(),
        }
    }
    fn blend(&self, s: f64) -> Vec<f64> {
        let ut = self.u_end();
        let u0 = &self.frames[0].u;
        let u: Vec<f64> = u0.iter().zip(&ut).map(|(a, b)| s * a + (1.0 - s) * b).collect();
        self.flow.borrow().lengths(&u)
    }
}

/// One output frame (a `run_flow_true` row).
#[derive(Clone, Debug)]
pub struct FrameOut {
    pub k: usize,
    pub t: f64,
    pub u: Vec<f64>,
    /// intrinsic angle defect and barycentric vertex area (K = k_defect / a_vertex)
    pub k_defect: Vec<f64>,
    pub a_vertex: Vec<f64>,
    pub lengths: Vec<f64>,
    pub vertices: Vec<V3>,
    /// Some only when the triangulation differs from the mesh (Polytope frame 0)
    pub faces: Option<Vec<[usize; 3]>>,
    pub len_rel_max: f64,
    pub len_rel_rms: f64,
    pub n_deep: usize,
    pub n_reflex: usize,
    pub min_dihedral: f64,
    pub gb_resid: f64,
    pub cone_defect_err: f64,
    pub rho: f64,
    pub k_dev_max: f64,
    pub seed_folds: bool,
    pub polytope: bool,
    pub start: String,
    pub warning: Option<String>,
    /// conformal-factor raster GRID_LAT x GRID_LON (empty until filled)
    pub cf_w: Vec<f64>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Phase {
    Flow,
    CleanStart,
    BackwardInit,
    Backward,
    Finalize,
    Emit,
    Done,
}

impl Phase {
    pub fn name(&self) -> &'static str {
        match self {
            Phase::Flow => "flowing (implicit BDF2 metric)",
            Phase::CleanStart => "loading the flow in 3D: starting shape",
            Phase::BackwardInit | Phase::Backward => "loading the flow in 3D (backward from the uniformized end)",
            Phase::Finalize => "polytope frame / quality",
            Phase::Emit => "streaming frames",
            Phase::Done => "done",
        }
    }
}

pub struct Session {
    pub cfg: SessionConfig,
    pub mesh: Mesh,
    pub flow_mesh: FlowMesh,
    /// C = 4 pi - sum of cone defects (> 0 here)
    pub c_curv: f64,
    phase: Phase,
    n_frames: usize,
    run: Option<FlowRun>,
    frames: Vec<Frame>,
    src: Option<FlowSrc>,
    rlz: Option<FastRealizer>,
    start: Option<Realized>,
    bw: Option<Backward>,
    mframes: Vec<MetricFrame>,
    pub stats_flow: Option<FlowStats>,
    pub seed_corners: Vec<usize>,
    seed_preview: Option<FrameOut>,
    rows: Vec<FrameOut>,
    emitted: usize,
    pub quality: Option<Quality>,
    cf_main: Option<CfGrid>,
    cf_poly: Option<CfGrid>,
    pub warning: Option<String>,
}

fn rms(v: &[f64]) -> f64 {
    (v.iter().map(|x| x * x).sum::<f64>() / v.len().max(1) as f64).sqrt()
}

impl Session {
    pub fn new(cfg: SessionConfig) -> Result<Session, String> {
        let n = cfg.positions.len();
        if n != cfg.betas.len() || n == 0 {
            return Err("need one beta per cone position".into());
        }
        let c_curv = 4.0 * PI - cfg.betas.iter().map(|b| 2.0 * PI * (1.0 - b)).sum::<f64>();
        if !(c_curv > 0.0) {
            return Err(NEG_REFUSAL.into());
        }
        let (mesh, flow_mesh) = seedmesh::build_mesh(&cfg.positions, &cfg.betas, cfg.level).map_err(|e| format!("mesh failed: {e}"))?;
        let flow = DiscreteConeFlow::new(mesh.clone(), cfg.betas.clone(), FlowOptions { seed: cfg.seed, ..Default::default() });
        let n_frames = flow.opts.times.len();
        let run = flow.start().map_err(|e| format!("seed metric failed: {e}"))?;
        let seed_corners = run.flow.seed_corners.clone();
        Ok(Session {
            cfg,
            mesh,
            flow_mesh,
            c_curv,
            phase: Phase::Flow,
            n_frames,
            run: Some(run),
            frames: vec![],
            src: None,
            rlz: None,
            start: None,
            bw: None,
            mframes: vec![],
            stats_flow: None,
            seed_corners,
            seed_preview: None,
            rows: vec![],
            emitted: 0,
            quality: None,
            cf_main: None,
            cf_poly: None,
            warning: None,
        })
    }

    pub fn phase(&self) -> Phase {
        self.phase
    }
    pub fn n_frames(&self) -> usize {
        self.n_frames
    }
    pub fn is_done(&self) -> bool {
        self.phase == Phase::Done
    }
    /// rho of the uniformized metric (known once the flow frames are computed; NaN before).
    pub fn rho(&self) -> f64 {
        self.frames.last().map_or(f64::NAN, |f| f.diag.rho)
    }
    pub fn frames_intrinsic(&self) -> &[Frame] {
        &self.frames
    }
    pub fn seed_preview(&self) -> Option<&FrameOut> {
        self.seed_preview.as_ref()
    }
    pub fn rows(&self) -> &[FrameOut] {
        &self.rows
    }
    /// Flow statistics after the whole run (frames + realizer in-between metrics).
    pub fn stats_total(&self) -> Option<FlowStats> {
        self.src.as_ref().map(|s| s.flow.borrow().stats.clone())
    }

    /// Fraction of the work done, in [0, 1].
    pub fn progress(&self) -> f64 {
        let n = self.n_frames.max(1) as f64;
        match self.phase {
            Phase::Flow => 0.25 * self.run.as_ref().map_or(0.0, |r| r.k as f64) / n,
            Phase::CleanStart => 0.25,
            Phase::BackwardInit => 0.33,
            Phase::Backward => 0.33 + 0.55 * self.bw.as_ref().map_or(0.0, |b| b.done() as f64) / n,
            Phase::Finalize => 0.9,
            Phase::Emit => 0.92 + 0.08 * self.emitted as f64 / n,
            Phase::Done => 1.0,
        }
    }

    fn row_of(&self, fr: &Frame, r: &Realized) -> FrameOut {
        FrameOut {
            k: fr.k,
            t: fr.t,
            u: fr.u.clone(),
            k_defect: fr.k_defect.clone(),
            a_vertex: fr.a_vertex.clone(),
            lengths: fr.l.clone(),
            vertices: r.x.clone(),
            faces: None,
            len_rel_max: r.len_rel_max,
            len_rel_rms: r.len_rel_rms,
            n_deep: realize::n_deep_folds(&r.x, &self.mesh.f, 0.15),
            n_reflex: r.n_reflex,
            min_dihedral: r.min_dihedral,
            gb_resid: fr.diag.gb_resid,
            cone_defect_err: fr.diag.cone_defect_err,
            rho: fr.diag.rho,
            k_dev_max: fr.diag.kdens_rel_dev_max,
            seed_folds: false,
            polytope: false,
            start: r.start.clone(),
            warning: r.warning.clone(),
            cf_w: vec![],
        }
    }

    /// `polytope_vertices` applied to a row (Polytope seed, t = 0); false if the fit failed.
    fn apply_polytope(&self, row: &mut FrameOut) -> bool {
        let rlz = self.rlz.as_ref().unwrap();
        let corners: Vec<usize> = self.mesh.cones.iter().chain(&self.seed_corners).cloned().collect();
        match realize::polytope_vertices(&rlz.mesh, &row.vertices, &corners, &row.lengths) {
            Some(p) => {
                row.vertices = p.x;
                row.faces = Some(p.faces);
                row.len_rel_max = p.len_rel_max;
                row.len_rel_rms = p.len_rel_rms;
                row.n_reflex = p.n_reflex;
                row.min_dihedral = p.min_dihedral;
                row.polytope = true;
                true
            }
            None => false,
        }
    }

    fn cf_for(&mut self, k: usize, t: f64) -> Option<&CfGrid> {
        let poly = !self.seed_corners.is_empty() && k == 0 && t == 0.0;
        let slot = if poly { &mut self.cf_poly } else { &mut self.cf_main };
        if slot.is_none() {
            let mut sing = self.mesh.cones.clone();
            let mut betas = self.cfg.betas.clone();
            if poly {
                let nc = self.seed_corners.len();
                sing.extend(&self.seed_corners);
                betas.extend(std::iter::repeat(1.0 - self.c_curv / nc as f64 / (2.0 * PI)).take(nc));
            }
            *slot = CfGrid::new(&self.mesh.y, &sing, &betas, GRID_LAT, GRID_LON).ok();
        }
        slot.as_ref()
    }

    fn fill_cf(&mut self, row: &mut FrameOut) {
        if row.cf_w.is_empty() {
            if let Some(g) = self.cf_for(row.k, row.t) {
                row.cf_w = g.grid(&row.u);
            }
        }
    }

    /// Do one chunk of work.  Returns false once everything is computed (frames may remain to
    /// be taken with [`Session::next_frame`]).
    pub fn step(&mut self) -> Result<bool, String> {
        match self.phase {
            Phase::Flow => {
                let run = self.run.as_mut().unwrap();
                match run.next_frame()? {
                    Some(fr) => self.frames.push(fr),
                    None => {
                        let run = self.run.take().unwrap();
                        self.stats_flow = Some(run.flow.stats.clone());
                        let flow = run.flow;
                        let mesh = &self.mesh;
                        let rm = RealizeMesh::new(mesh.f.clone(), mesh.e.clone(), mesh.cones.clone(), mesh.nv());
                        let (x0, kind) = realize::fast_seed(&mesh.y, &mesh.f, &mesh.cones, &self.cfg.betas);
                        if kind != SeedKind::Caps {
                            self.warning = Some("cone-cap seed unavailable: started from the bare sphere".into());
                        }
                        self.rlz = Some(FastRealizer::new(rm, x0, kind));
                        self.mframes =
                            self.frames.iter().map(|f| MetricFrame { k: f.k, t: f.t, lengths: f.l.clone() }).collect();
                        self.src = Some(FlowSrc::new(flow, std::mem::take(&mut self.frames)));
                        self.frames = self.src.as_ref().unwrap().frames.clone();
                        self.phase = Phase::CleanStart;
                    }
                }
            }
            Phase::CleanStart => {
                let (rlz, src) = (self.rlz.as_ref().unwrap(), self.src.as_ref().unwrap());
                self.start = realize::clean_start(rlz, src, &self.frames[0].l, 6);
                if let Some(s) = self.start.clone() {
                    let mut sd = self.row_of(&self.frames[0], &s);
                    if self.cfg.seed == Seed::Polytope {
                        self.apply_polytope(&mut sd);
                    }
                    self.fill_cf(&mut sd);
                    self.seed_preview = Some(sd);
                }
                self.phase = Phase::BackwardInit;
            }
            Phase::BackwardInit => {
                let (rlz, src) = (self.rlz.as_ref().unwrap(), self.src.as_ref().unwrap());
                self.bw = Some(Backward::new(rlz, src, &self.mframes, self.start.take()));
                self.phase = Phase::Backward;
            }
            Phase::Backward => {
                let (rlz, src) = (self.rlz.as_ref().unwrap(), self.src.as_ref().unwrap());
                if !self.bw.as_mut().unwrap().step(rlz, src, &self.mframes) {
                    self.phase = Phase::Finalize;
                }
            }
            Phase::Finalize => {
                let res = self.bw.take().unwrap().finish();
                let mut rows: Vec<FrameOut> = self.frames.iter().zip(&res).map(|(f, r)| self.row_of(f, r)).collect();
                let mut qrows: Vec<FrameRow> = rows
                    .iter()
                    .map(|r| FrameRow { len_rel_max: Some(r.len_rel_max), len_rel_rms: None, n_deep: r.n_deep, ..Default::default() })
                    .collect();
                if self.cfg.seed == Seed::Polytope {
                    // quality reads the smooth fit's errors; the display shows the polytope's own
                    qrows[0].q_len_rel_max = qrows[0].len_rel_max;
                    qrows[0].q_len_rel_rms = None;
                    let mut r0 = rows[0].clone();
                    r0.polytope = true;
                    self.apply_polytope(&mut r0);
                    r0.polytope = true;
                    qrows[0].len_rel_max = Some(r0.len_rel_max);
                    // Python: the q_len_rel_rms key is present (None), so the polytope's rms is not read
                    qrows[0].len_rel_rms = None;
                    rows[0] = r0;
                }
                mark_seed_folds(&mut qrows, self.cfg.seed == Seed::Polytope);
                for (r, q) in rows.iter_mut().zip(&qrows) {
                    r.seed_folds = q.seed_folds;
                }
                self.quality = flow_quality(&qrows);
                self.rows = rows;
                self.phase = Phase::Emit;
            }
            Phase::Emit => {
                if self.emitted >= self.rows.len() {
                    self.phase = Phase::Done;
                }
                return Ok(self.phase != Phase::Done);
            }
            Phase::Done => return Ok(false),
        }
        Ok(true)
    }

    /// Advance: returns the next output frame (with its cf raster) once the realization is done,
    /// otherwise performs one work chunk and returns None.  Check [`Session::is_done`].
    pub fn next_frame(&mut self) -> Result<Option<FrameOut>, String> {
        if self.phase == Phase::Emit {
            if self.emitted < self.rows.len() {
                let mut r = self.rows[self.emitted].clone();
                self.fill_cf(&mut r);
                self.rows[self.emitted].cf_w = r.cf_w.clone();
                self.emitted += 1;
                if self.emitted == self.rows.len() {
                    self.phase = Phase::Done;
                }
                return Ok(Some(r));
            }
            self.phase = Phase::Done;
            return Ok(None);
        }
        if self.phase != Phase::Done {
            self.step()?;
        }
        Ok(None)
    }

    /// Run to completion; returns all frames.
    pub fn run_all(&mut self) -> Result<Vec<FrameOut>, String> {
        let mut out = vec![];
        while !self.is_done() {
            if let Some(f) = self.next_frame()? {
                out.push(f);
            }
        }
        Ok(out)
    }

    /// Mean |len_rel| rms helper for callers that want the rms of an arbitrary shape.
    pub fn len_rel_rms(&self, x: &[V3], lengths: &[f64]) -> f64 {
        rms(&realize::geom::len_rel(x, &self.mesh.e, lengths))
    }
}
