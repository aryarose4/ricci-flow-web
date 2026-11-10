//! End-to-end /flow "True Ricci Flow" pipeline (rho > 0): port of
//! `viewer/backend/flow_true.run_flow_true` (Fast realizer, backward realization) as a steppable
//! [`Session`]: mesh from cone positions -> intrinsic flow frames -> clean start -> backward 3D
//! realization -> Polytope frame 0 + seed-fold exemption + quality note -> per-frame rows with the
//! conformal-factor raster.
//!
//! rho < 0 (port of `run_flow_true_negative`): the intrinsic flow drives and each frame is realized
//! FORWARD by the Fast realizer warm-started from the previous frame; [`crate::tracker::Tracker`]
//! latches the first frame the 3D shape cannot follow.  Frames are produced one per work chunk.
//! rho = 0 (flat) is refused.

use crate::cfgrid::{CfGrid, GRID_LAT, GRID_LON};
use crate::flow::{default_times, DiscreteConeFlow, FlowOptions, FlowRun, FlowStats, Frame, Seed};
use crate::mesh::Mesh;
use crate::quality::{flow_quality, mark_seed_folds, FrameRow, Quality};
use crate::realize::{self, Backward, FastRealizer, FlowSource, MetricFrame, RealizeMesh, Realized, SeedKind};
use crate::realize::geom::V3;
use crate::seedmesh::{self, FlowMesh};
use crate::tracker::{self, Mismatch, State, Tracker};
use std::cell::RefCell;
use std::f64::consts::PI;

pub const NEG_REFUSAL: &str = "True Ricci Flow needs rho != 0 (these cone angles give a flat metric: there is no curvature to flow to)";

/// Session inputs.
#[derive(Clone, Debug)]
pub struct SessionConfig {
    /// cone positions in the stereographic plane (x, y) = (Re z, Im z)
    pub positions: Vec<[f64; 2]>,
    pub betas: Vec<f64>,
    pub seed: Seed,
    /// (N Fibonacci points, m spokes, K rings); see [`seedmesh::page_level`]
    pub level: (usize, usize, usize),
    /// opt-in anti-crumple options for the 3D realization (default: off = the Python port)
    pub anti: AntiCrumple,
}

/// Opt-in anti-crumple options (lane anti-crumple).  `Default` = off: bit-identical to the port.
#[derive(Clone, Debug, Default)]
pub struct AntiCrumple {
    /// (A) adaptive hidden intermediate metrics between output frames
    pub adaptive: Option<realize::Adaptive>,
    /// (B) fading fold-angle prior toward the reference (final / last good) frame's dihedrals
    pub prior: Option<realize::FoldPrior>,
    /// (C) rho > 0: rescue a folded uniformized end from the last clean cold fit among the
    /// `rescue_final` frames before it (0 = off)
    pub rescue_final: usize,
    /// (F, lane forward-flow) rho > 0: realize FORWARD from the t = 0 seed instead of backward
    /// from the uniformized end (A/B/C are then unused)
    pub forward: Option<realize::ForwardOpts>,
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
    /// rho < 0 tracking state ("searching" | "following" | "suspect" | "failed" | "lost"); "" for rho > 0
    pub embed_state: String,
    /// rho < 0: the vertices are not a passing fit of this frame (searching / failed / lost)
    pub embed_lost: bool,
    /// rho < 0: mean / rms ball-averaged curvature error of the 3D attempt, in units of |rho| (NaN for rho > 0)
    pub k_err_signed: f64,
    pub k_err_rms: f64,
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
    /// rho < 0: flow one frame and realize it forward (continuation + tracker)
    Track,
    /// Fast Uniformizer: direct implicit solve for the uniformized metric
    FastSolve,
    /// Fast Uniformizer: 3D realization of that one metric
    FastRealize,
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
            Phase::Track => "flowing (true Ricci flow, rho < 0: implicit BDF2 metric, 3D shape by continuation)",
            Phase::FastSolve => "Fast Uniformizer: solving for the uniformized metric (implicit Newton from the far field)",
            Phase::FastRealize => "Fast Uniformizer: building the 3D surface",
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
    fw: Option<realize::Forward>,
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
    /// rho < 0 state
    pub neg: bool,
    pub tracker: Option<Tracker>,
    prev: Option<Realized>,
    good: Option<Realized>,
    notice: Option<String>,
    /// True Polytope seed: the polytope (mesh, 3D positions, creases)
    pub true_poly: Option<crate::truepoly::TruePoly>,
    /// hidden in-between metrics fitted by the realization (diagnostic)
    pub n_hidden: usize,
    /// (B) was requested but switched off because the uniformized end folded
    pub prior_gated: bool,
    /// (C) the uniformized end was re-fitted
    pub rescued: bool,
    /// (F hybrid) frames taken from the backward realization
    pub n_from_backward: usize,
    /// rho < 0: previous intrinsic frame (for the opt-in hidden metrics)
    prev_fr: Option<Frame>,
    /// Fast Uniformizer mode (see [`Session::set_fast`])
    pub fast: bool,
    /// Fast Uniformizer: intrinsic solve statistics
    pub unif: Option<crate::fastunif::UnifStats>,
}

fn rms(v: &[f64]) -> f64 {
    (v.iter().map(|x| x * x).sum::<f64>() / v.len().max(1) as f64).sqrt()
}

impl Session {
    pub fn new(cfg: SessionConfig) -> Result<Session, String> {
        Self::new_with(cfg, false)
    }

    /// Cheap session for slider / drag previews: True Polytope skips its conformal coordinates.
    pub fn new_preview(cfg: SessionConfig) -> Result<Session, String> {
        Self::new_with(cfg, true)
    }

    fn new_with(cfg: SessionConfig, preview: bool) -> Result<Session, String> {
        let n = cfg.positions.len();
        if n != cfg.betas.len() || n == 0 {
            return Err("need one beta per cone position".into());
        }
        let c_curv = 4.0 * PI - cfg.betas.iter().map(|b| 2.0 * PI * (1.0 - b)).sum::<f64>();
        if !(c_curv.abs() > 1e-12) {
            return Err(NEG_REFUSAL.into());
        }
        crate::fastunif::admissibility(&cfg.betas)?;
        let neg = c_curv < 0.0;
        let mut tp_fallback: Option<String> = None;
        let true_poly = if cfg.seed == Seed::TruePolytope {
            let b = if preview { crate::truepoly::build_preview } else { crate::truepoly::build };
            match b(&cfg.positions, &cfg.betas, cfg.level) {
                Ok(tp) => Some(tp),
                // rho < 0: the star-shaped polyhedron is not available for every configuration;
                // fall back to the Polytope seed and say so
                Err(e) if neg => {
                    tp_fallback = Some(format!("True Polytope unavailable for these cones ({e}); started from the Polytope seed instead"));
                    None
                }
                Err(e) => return Err(format!("mesh failed: {e}")),
            }
        } else {
            None
        };
        let mut cfg = cfg;
        if tp_fallback.is_some() {
            cfg.seed = Seed::Polytope;
        }
        let (mesh, flow_mesh) = match &true_poly {
            Some(tp) => (tp.mesh.clone(), tp.flow_mesh.clone()),
            None => seedmesh::build_mesh(&cfg.positions, &cfg.betas, cfg.level).map_err(|e| format!("mesh failed: {e}"))?,
        };
        let mut flow = DiscreteConeFlow::new(mesh.clone(), cfg.betas.clone(), FlowOptions { seed: cfg.seed, ..Default::default() });
        if let Some(tp) = &true_poly {
            flow.l0 = tp.l0.clone();
        }
        let n_frames = flow.opts.times.len();
        let run = flow.start().map_err(|e| format!("seed metric failed: {e}"))?;
        let seed_corners = match &true_poly {
            Some(tp) => tp.corners.clone(),
            None => run.flow.seed_corners.clone(),
        };
        Ok(Session {
            cfg,
            mesh,
            flow_mesh,
            c_curv,
            phase: if neg { Phase::Track } else { Phase::Flow },
            n_frames,
            run: Some(run),
            frames: vec![],
            src: None,
            rlz: None,
            start: None,
            bw: None,
            fw: None,
            mframes: vec![],
            stats_flow: None,
            seed_corners,
            seed_preview: None,
            rows: vec![],
            emitted: 0,
            quality: None,
            cf_main: None,
            cf_poly: None,
            warning: tp_fallback.or_else(|| true_poly.as_ref().and_then(|t| t.warning.clone())),
            neg,
            tracker: None,
            prev: None,
            good: None,
            notice: None,
            true_poly,
            n_hidden: 0,
            prior_gated: false,
            rescued: false,
            n_from_backward: 0,
            prev_fr: None,
            fast: false,
            unif: None,
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
        if self.neg && !self.fast {
            return if self.phase == Phase::Track { f64::NAN } else { self.rows.last().map_or(f64::NAN, |r| r.rho) };
        }
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
            Phase::Backward => 0.33 + 0.55 * self.bw.as_ref().map_or_else(|| self.fw.as_ref().map_or(0.0, |f| f.done() as f64), |b| b.done() as f64) / n,
            Phase::Finalize => 0.9,
            Phase::Emit if self.neg || self.fast => 0.99,
            Phase::Emit => 0.92 + 0.08 * self.emitted as f64 / n,
            Phase::Track => 0.98 * self.rows.len() as f64 / n,
            Phase::FastSolve => 0.05,
            Phase::FastRealize => 0.5,
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
            embed_state: String::new(),
            embed_lost: false,
            k_err_signed: f64::NAN,
            k_err_rms: f64::NAN,
            cf_w: vec![],
        }
    }

    /// rho < 0 banner text once the 3D tracking has failed (`embed_notice`); after the run it is the
    /// final summary's (an unconfirmed bad last frame or a never-started run also count).
    pub fn embed_notice(&self) -> Option<String> {
        self.notice.clone()
    }

    /// rho < 0: one intrinsic frame + its forward 3D realization; false when the flow is done.
    fn track_step(&mut self) -> Result<bool, String> {
        if self.rlz.is_none() {
            let mesh = &self.mesh;
            let rm = RealizeMesh::new(mesh.f.clone(), mesh.e.clone(), mesh.cones.clone(), mesh.nv());
            let (x0, kind) = realize::fast_seed(&mesh.y, &mesh.f, &mesh.cones, &self.cfg.betas);
            self.rlz = Some(FastRealizer::new(rm, x0, kind));
            self.tracker = Some(Tracker::new(&mesh.f, &mesh.cones));
        }
        let fr = match self.run.as_mut().unwrap().next_frame()? {
            Some(fr) => fr,
            None => {
                let run = self.run.take().unwrap();
                self.stats_flow = Some(run.flow.stats.clone());
                let fail = self.tracker.as_mut().unwrap().finish();
                self.notice = tracker::embed_notice(&fail);
                let qrows: Vec<FrameRow> = self
                    .rows
                    .iter()
                    .map(|r| FrameRow {
                        len_rel_max: Some(r.len_rel_max),
                        len_rel_rms: Some(r.len_rel_rms),
                        n_deep: r.n_deep,
                        embed_lost: r.embed_lost,
                        ..Default::default()
                    })
                    .collect();
                self.quality = flow_quality(&qrows);
                return Ok(false);
            }
        };
        if self.cfg.anti.forward.is_some() {
            return self.neg_forward_step(fr);
        }
        let tp0 = if self.prev.is_none() { self.true_poly_start(&fr.l) } else { None };
        let (rlz, tr) = (self.rlz.as_ref().unwrap(), self.tracker.as_mut().unwrap());
        let rho = fr.diag.rho;
        let (rf, m): (Realized, Mismatch) = match &self.prev {
            None => match tp0 {
                Some(rf) => {
                    let m = tracker::mismatch(tr, &rf.x, &self.mesh.e, &fr.l, &fr.k_defect, &fr.a_vertex, rho);
                    (rf, m)
                }
                None => tracker::cold_start(rlz, tr, &fr.l, &fr.k_defect, &fr.a_vertex, rho),
            },
            Some(p) => {
                // opt-in (A) hidden metrics from the previous intrinsic frame, (B) prior toward the previous
                // accepted shape (the final frame is not realizable here; see PROGRESS_anti-crumple.md)
                let anti = self.cfg.anti.clone();
                let t_end = *self.run.as_ref().unwrap().flow.opts.times.last().unwrap();
                if let Some(pr) = &anti.prior {
                    rlz.set_prior_ref(Some(&p.x), pr.hinge);
                }
                let mut c = p.clone();
                if let (Some(a), Some(pf)) = (&anti.adaptive, &self.prev_fr) {
                    let dl = pf.l.iter().zip(&fr.l).map(|(a, b)| (b / a).ln().abs()).fold(0.0, f64::max);
                    let nsub = ((dl / a.step_rel).ceil() as usize).min(a.max_sub).max(1) - 1;
                    for j in 1..=nsub {
                        let t = pf.t + (fr.t - pf.t) * j as f64 / (nsub + 1) as f64;
                        if let Some(pr) = &anti.prior {
                            rlz.set_prior_w(pr.weight(t, t_end));
                        }
                        let lt = match self.run.as_mut().unwrap().flow.between(pf, t) {
                            Ok(f) => f.l,
                            Err(_) => continue,
                        };
                        let rm = if a.cheap { rlz.realize_cheap(&lt, &c) } else { rlz.realize(&lt, Some(&c)) };
                        self.n_hidden += 1;
                        if rm.ok {
                            c = rm;
                        }
                    }
                }
                if let Some(pr) = &anti.prior {
                    rlz.set_prior_w(pr.weight(fr.t, t_end));
                }
                let rf = rlz.realize(&fr.l, Some(&c));
                if anti.prior.is_some() {
                    rlz.set_prior_ref(None, false);
                    rlz.set_prior_w(0.0);
                }
                let m = tracker::mismatch(tr, &rf.x, &self.mesh.e, &fr.l, &fr.k_defect, &fr.a_vertex, rho);
                (rf, m)
            }
        };
        let state = tr.update(fr.t, &m, rf.ok);
        if state == State::Following {
            self.good = Some(rf.clone());
        }
        if state != State::Searching && rf.ok {
            self.prev = Some(rf.clone());
        }
        if state == State::Failed {
            self.notice = tracker::embed_notice(&tr.summary());
        }
        let mut row = self.row_of(&fr, &rf);
        if self.cfg.anti.adaptive.is_some() {
            self.prev_fr = Some(fr.clone());
        }
        if let Some(g) = &self.good {
            row.vertices = g.x.clone();
        }
        row.n_deep = self.n_deep_off_crease(&row.vertices);
        row.len_rel_max = m.len_rel_max;
        row.len_rel_rms = m.len_rel_rms;
        row.k_err_signed = m.k_err_signed;
        row.k_err_rms = m.k_err_rms;
        row.embed_state = state.name().to_string();
        row.embed_lost = state.lost();
        self.rows.push(row);
        Ok(true)
    }

    /// Fast Uniformizer: switch this (not yet stepped) session to the direct solve.  The seed
    /// metric is ignored: the solve always starts from the flat-cone far field `v_asymp`.
    pub fn set_fast(&mut self) {
        self.fast = true;
        self.phase = Phase::FastSolve;
        self.n_frames = 1;
        self.seed_corners.clear();
    }

    /// Fast Uniformizer: realize the uniformized metric (rho > 0: cold Fast fits from the cone-cap
    /// seed, a few seed scales, cleanest kept; rho < 0: the tracker's multi-start cold fit and its
    /// curvature mismatch, best uniform-curvature attempt, NOT certified as an embedding).
    fn fast_realize(&mut self) {
        let mesh = &self.mesh;
        let rm = RealizeMesh::new(mesh.f.clone(), mesh.e.clone(), mesh.cones.clone(), mesh.nv());
        let (x0, kind) = realize::fast_seed(&mesh.y, &mesh.f, &mesh.cones, &self.cfg.betas);
        if kind != SeedKind::Caps {
            self.warning = Some("cone-cap seed unavailable: started from the bare sphere".into());
        }
        let rlz = FastRealizer::new(rm, x0, kind);
        let fr = self.frames[0].clone();
        let mut row;
        let mut lost = false;
        if self.neg {
            let tr = Tracker::new(&mesh.f, &mesh.cones);
            let (rf, m) = tracker::cold_start(&rlz, &tr, &fr.l, &fr.k_defect, &fr.a_vertex, fr.diag.rho);
            let reasons = tr.reasons(&m, rf.ok);
            lost = !reasons.is_empty();
            row = self.row_of(&fr, &rf);
            row.len_rel_max = m.len_rel_max;
            row.len_rel_rms = m.len_rel_rms;
            row.k_err_signed = m.k_err_signed;
            row.k_err_rms = m.k_err_rms;
            row.embed_state = if lost { "failed" } else { "following" }.into();
            row.embed_lost = lost;
            if lost {
                self.notice = Some(format!(
                    "Fast Uniformizer (rho < 0): no faithful 3D shape found ({}); showing the best attempt",
                    reasons.join("; ")
                ));
            }
        } else {
            let mut best = rlz.realize(&fr.l, None);
            for s in [1.0, 0.9, 1.1] {
                if best.ok && realize::n_deep_folds(&best.x, &mesh.f, 0.15) == 0 {
                    break;
                }
                let xs: Vec<V3> = mesh.y.iter().map(|p| realize::geom::scale(*p, s)).collect();
                let b = rlz.realize_cold_from(&fr.l, &xs);
                if realize::cleaner_picks_b(&rlz, &best, &b) {
                    best = b;
                }
            }
            row = self.row_of(&fr, &best);
        }
        let q = FrameRow {
            len_rel_max: Some(row.len_rel_max),
            len_rel_rms: if self.neg { Some(row.len_rel_rms) } else { None },
            n_deep: row.n_deep,
            embed_lost: lost,
            ..Default::default()
        };
        self.quality = flow_quality(&[q]).map(|mut q| {
            if q.ok {
                q.label = "Certified Target".into();
                q.detail = "The 3D shape matches the uniformized metric (edge lengths within 1 % rms / 5 % max, no folds); NOT a certified embedding".into();
            }
            q
        });
        self.rlz = Some(rlz);
        self.rows = vec![row];
    }

    /// rho < 0 with `anti.forward` (lane neg-forward): guarded forward march.  Every frame is fitted
    /// WARM from the previously SHOWN shape (also past tracking failure: no cold restarts, no frozen
    /// "last good" shape) through hidden flow metrics; a fit that adds deep folds is re-fitted with a
    /// hinge fold prior toward the previous shape at escalating weights, and the candidate with the
    /// fewest deep folds wins (ties: lower length rms).  Fold-free beats a smaller length error: where
    /// R^3 cannot hold the metric the error grows and is reported, instead of hidden in folds.
    fn neg_forward_step(&mut self, fr: Frame) -> Result<bool, String> {
        let fo = self.cfg.anti.forward.clone().unwrap();
        let rho = fr.diag.rho;
        let t_frame = fr.t;
        // the march starts only once some frame has passed (followed); while still searching, keep the
        // tracker's multi-seed cold start (marching warm from a folded never-passing start keeps its folds)
        let started = self.tracker.as_ref().map_or(false, |t| t.last_good_t.is_some());
        let rf: Realized = match (&self.prev, &self.prev_fr) {
            (Some(p), Some(pf)) if started => {
                let rlz = self.rlz.as_ref().unwrap();
                let flow = &mut self.run.as_mut().unwrap().flow;
                let dl = pf.l.iter().zip(&fr.l).map(|(a, b)| (b / a).ln().abs()).fold(0.0, f64::max);
                let nsub = ((dl / fo.step_rel).ceil() as usize).min(fo.max_sub).max(1) - 1;
                let mut lts: Vec<(usize, Vec<Vec<f64>>)> = vec![];
                let mut hidden = |n: usize, flow: &mut crate::flow::DiscreteConeFlow| -> Vec<Vec<f64>> {
                    if let Some((_, v)) = lts.iter().find(|(m, _)| *m == n) {
                        return v.clone();
                    }
                    let v: Vec<Vec<f64>> = (1..=n)
                        .filter_map(|j| flow.between(pf, pf.t + (t_frame - pf.t) * j as f64 / (n + 1) as f64).ok().map(|f| f.l))
                        .collect();
                    lts.push((n, v.clone()));
                    v
                };
                let d0 = realize::n_deep_folds(&p.x, &self.mesh.f, 0.15);
                let mut nh = 0;
                let mut chain = |ls: &[Vec<f64>], w: f64| -> Realized {
                    if w > 0.0 {
                        rlz.set_prior_ref(Some(&p.x), fo.hinge);
                    }
                    rlz.set_prior_w(w);
                    let mut c = Realized { ok: true, ..p.clone() };
                    for lt in ls {
                        let mut rm = rlz.realize_cheap(lt, &c);
                        rm.ok = true; // keep marching warm even where R^3 cannot hold the metric
                        c = rm;
                        nh += 1;
                    }
                    let mut r = rlz.realize(&fr.l, Some(&c));
                    rlz.set_prior_w(0.0);
                    rlz.set_prior_ref(None, false);
                    r.x = realize::geom::kabsch_align(&r.x, &p.x);
                    r
                };
                let key = |r: &Realized| (realize::n_deep_folds(&r.x, &self.mesh.f, 0.15), r.len_rel_rms);
                let ls = hidden(nsub, flow);
                let mut best = chain(&ls, 0.0);
                let mut kb = key(&best);
                if kb.0 > d0 {
                    let ls2 = hidden((2 * nsub + 1).min(4 * fo.max_sub), flow);
                    for (i, w) in [fo.w_guard, 10.0 * fo.w_guard, 100.0 * fo.w_guard].into_iter().enumerate() {
                        let g = chain(if i == 0 { &ls } else { &ls2 }, w);
                        let kg = key(&g);
                        if kg.0 < kb.0 || (kg.0 == kb.0 && kg.1 < kb.1) {
                            best = g;
                            kb = kg;
                        }
                        if kb.0 <= d0 {
                            break;
                        }
                    }
                    best.start += "+guard";
                }
                self.n_hidden += nh;
                best
            }
            _ => match self.true_poly_start(&fr.l).filter(|_| self.prev.is_none()) {
                // True Polytope (rho < 0): frame 0 IS the star-shaped polyhedron
                Some(rf) => rf,
                None => {
                    let (rlz, tr) = (self.rlz.as_ref().unwrap(), self.tracker.as_ref().unwrap());
                    tracker::cold_start(rlz, tr, &fr.l, &fr.k_defect, &fr.a_vertex, rho).0
                }
            },
        };
        let tr = self.tracker.as_mut().unwrap();
        let m = tracker::mismatch(tr, &rf.x, &self.mesh.e, &fr.l, &fr.k_defect, &fr.a_vertex, rho);
        let state = tr.update(fr.t, &m, rf.ok);
        if state == State::Failed || (state == State::Lost && self.notice.is_some()) {
            let base = tracker::embed_notice(&tr.summary()).unwrap_or_default();
            self.notice = Some(format!(
                "{} Shown instead: the fold-free shape closest to the flow, now {}% off in edge lengths (max {}%) and {:+.0}% off in curvature at t = {}.",
                base,
                tracker::fmt_g(100.0 * m.len_rel_rms, 2),
                tracker::fmt_g(100.0 * m.len_rel_max, 2),
                100.0 * m.k_err_signed,
                tracker::fmt_g(fr.t, 3)
            ));
        }
        let mut row = self.row_of(&fr, &rf);
        row.n_deep = self.n_deep_off_crease(&row.vertices);
        row.len_rel_max = m.len_rel_max;
        row.len_rel_rms = m.len_rel_rms;
        row.k_err_signed = m.k_err_signed;
        row.k_err_rms = m.k_err_rms;
        row.embed_state = state.name().to_string();
        row.embed_lost = state.lost();
        self.prev = Some(rf);
        self.prev_fr = Some(fr);
        self.rows.push(row);
        Ok(true)
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

    /// True Polytope: the polyhedron scaled to the gauge of lengths `l` (exact edge lengths).
    fn true_poly_start(&self, l: &[f64]) -> Option<Realized> {
        let tp = self.true_poly.as_ref()?;
        let e = &self.mesh.e;
        let s = l.iter().zip(&tp.l0).map(|(a, b)| a / b).sum::<f64>() / e.len() as f64;
        let x: Vec<V3> = tp.x.iter().map(|q| realize::geom::scale(*q, s)).collect();
        let rel = realize::geom::len_rel(&x, e, l);
        let dh = realize::geom::dihedrals(&x, &self.mesh.f);
        Some(Realized {
            ok: true,
            len_rel_max: rel.iter().fold(0.0f64, |m, v| m.max(v.abs())),
            len_rel_rms: rms(&rel),
            n_reflex: dh.iter().filter(|&&d| d < -1e-9).count(),
            min_dihedral: dh.iter().cloned().fold(f64::INFINITY, f64::min),
            lm_iters: 0,
            start: "true_polytope".into(),
            warning: None,
            x,
        })
    }

    /// deep folds that are NOT concave crease edges of the True Polytope (rho < 0 polyhedra have
    /// saddle corners, so some of their creases are concave by construction)
    fn n_deep_off_crease(&self, x: &[V3]) -> usize {
        let Some(tp) = &self.true_poly else { return realize::n_deep_folds(x, &self.mesh.f, 0.15) };
        let cr: std::collections::HashSet<(usize, usize)> = tp.crease_edges.iter().map(|&k| (self.mesh.e[k][0], self.mesh.e[k][1])).collect();
        let dh = realize::geom::dihedrals(x, &self.mesh.f);
        realize::geom::face_edge_table(&self.mesh.f)
            .iter()
            .zip(&dh)
            .filter(|(([a, b], _), &d)| d < -0.15 && !cr.contains(&(*a.min(b), *a.max(b))))
            .count()
    }

    /// True Polytope frame 0: the polytope itself (scaled to the gauge), rigidly aligned to the fit.
    fn apply_true_polytope(&self, row: &mut FrameOut) {
        let Some(tp) = &self.true_poly else { return };
        let e = &self.mesh.e;
        let s = row.lengths.iter().zip(&tp.l0).map(|(l, l0)| l / l0).sum::<f64>() / e.len() as f64;
        let xs: Vec<V3> = tp.x.iter().map(|q| realize::geom::scale(*q, s)).collect();
        let x = realize::geom::kabsch_align(&xs, &row.vertices);
        let rel = realize::geom::len_rel(&x, e, &row.lengths);
        let dh = realize::geom::dihedrals(&x, &self.mesh.f);
        row.len_rel_max = rel.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        row.len_rel_rms = rms(&rel);
        row.n_reflex = dh.iter().filter(|&&d| d < -1e-9).count();
        row.min_dihedral = dh.iter().cloned().fold(f64::INFINITY, f64::min);
        row.n_deep = self.n_deep_off_crease(&x);
        row.vertices = x;
        row.polytope = true;
    }

    fn cf_for(&mut self, k: usize, t: f64) -> Option<&CfGrid> {
        let poly = !self.seed_corners.is_empty() && k == 0 && t == 0.0 && self.cfg.seed == Seed::Polytope;
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
            // True Polytope: u is relative to the polytope; add its fitted offset vs the round sphere
            let u = match &self.true_poly {
                Some(tp) => row.u.iter().zip(&tp.u_off).map(|(a, b)| a + b).collect(),
                None => row.u.clone(),
            };
            if let Some(g) = self.cf_for(row.k, row.t) {
                row.cf_w = g.grid(&u);
            }
        }
    }

    /// Do one chunk of work.  Returns false once everything is computed (frames may remain to
    /// be taken with [`Session::next_frame`]).
    pub fn step(&mut self) -> Result<bool, String> {
        match self.phase {
            Phase::FastSolve => {
                let mut flow = self.run.take().unwrap().flow;
                let (u, st) = crate::fastunif::uniformize(&mut flow)?;
                self.stats_flow = Some(flow.stats.clone());
                self.unif = Some(st);
                let fr = flow.make_frame(0, 0.0, &u);
                self.frames = vec![fr];
                self.phase = Phase::FastRealize;
            }
            Phase::FastRealize => {
                self.fast_realize();
                self.phase = Phase::Emit;
            }
            Phase::Track => {
                if !self.track_step()? {
                    self.phase = if self.emitted >= self.rows.len() { Phase::Done } else { Phase::Emit };
                }
            }
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
                    self.apply_true_polytope(&mut sd);
                    self.fill_cf(&mut sd);
                    self.seed_preview = Some(sd);
                }
                self.phase = Phase::BackwardInit;
            }
            Phase::BackwardInit => {
                let (rlz, src) = (self.rlz.as_ref().unwrap(), self.src.as_ref().unwrap());
                let a = &self.cfg.anti;
                if let Some(fo) = a.forward.clone() {
                    self.fw = Some(realize::Forward::new(rlz, &self.mframes, self.start.take(), fo));
                    self.phase = Phase::Backward;
                    return Ok(true);
                }
                let mut bw = Backward::new(rlz, src, &self.mframes, self.start.take());
                if a.rescue_final > 0 {
                    bw.rescue_final(rlz, &self.mframes, a.rescue_final);
                }
                self.bw = Some(bw.with_options(rlz, a.adaptive.clone(), a.prior.clone()));
                self.phase = Phase::Backward;
            }
            Phase::Backward => {
                let (rlz, src) = (self.rlz.as_ref().unwrap(), self.src.as_ref().unwrap());
                if let Some(fw) = self.fw.as_mut() {
                    if !fw.step(rlz, src, &self.mframes) {
                        self.phase = Phase::Finalize;
                    }
                    return Ok(true);
                }
                if !self.bw.as_mut().unwrap().step(rlz, src, &self.mframes) {
                    self.phase = Phase::Finalize;
                }
            }
            Phase::Finalize => {
                let res = if let Some(fw) = self.fw.take() {
                    self.n_hidden = fw.n_hidden;
                    let hy = fw.opts.hybrid_rescue;
                    let res = fw.finish();
                    let rlz = self.rlz.as_ref().unwrap();
                    if hy > 0 && realize::any_deep(rlz, &res) {
                        let src = self.src.as_ref().unwrap();
                        let mut bw = Backward::new(rlz, src, &self.mframes, Some(res[0].clone()));
                        bw.rescue_final(rlz, &self.mframes, hy);
                        while bw.step(rlz, src, &self.mframes) {}
                        self.rescued = bw.rescued;
                        let (m, nb) = realize::merge_cleaner(rlz, res, bw.finish());
                        self.n_from_backward = nb;
                        m
                    } else {
                        res
                    }
                } else {
                let bw = self.bw.take().unwrap();
                self.n_hidden = bw.n_hidden;
                self.prior_gated = bw.prior_gated;
                self.rescued = bw.rescued;
                if self.cfg.anti.prior.is_some() {
                    let rlz = self.rlz.as_ref().unwrap();
                    rlz.set_prior_ref(None, false);
                    rlz.set_prior_w(0.0);
                }
                bw.finish()
                };
                let mut rows: Vec<FrameOut> = self.frames.iter().zip(&res).map(|(f, r)| self.row_of(f, r)).collect();
                let mut qrows: Vec<FrameRow> = rows
                    .iter()
                    // Python make_row: len_rel_rms = rf.diag.get("len_rel_rms") is None for the Fast realizer
                    // (true_flow.FastRealizer's diag has no len_rel_rms key; fixtures store NaN), so
                    // frame_quality checks only len_rel_max there.
                    .map(|r| FrameRow { len_rel_max: Some(r.len_rel_max), len_rel_rms: None, n_deep: r.n_deep, ..Default::default() })
                    .collect();
                if self.cfg.seed == Seed::Polytope {
                    // quality reads the smooth fit's errors; the display shows the polytope's own
                    qrows[0].q_len_rel_max = qrows[0].len_rel_max;
                    qrows[0].q_len_rel_rms = qrows[0].len_rel_rms; // the fit's: None (see above)
                    qrows[0].has_q = true;
                    let mut r0 = rows[0].clone();
                    r0.polytope = true;
                    self.apply_polytope(&mut r0);
                    r0.polytope = true;
                    qrows[0].len_rel_max = Some(r0.len_rel_max);
                    // display values; quality reads the present q_ keys (has_q), never these
                    qrows[0].len_rel_rms = Some(r0.len_rel_rms);
                    rows[0] = r0;
                }
                if self.true_poly.is_some() {
                    // frame 0 IS the polytope: exact lengths, no folds
                    self.apply_true_polytope(&mut rows[0]);
                    qrows[0].len_rel_max = Some(rows[0].len_rel_max);
                    qrows[0].n_deep = rows[0].n_deep;
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
        if self.phase == Phase::Emit || (self.phase == Phase::Track && self.emitted < self.rows.len()) {
            if self.emitted < self.rows.len() {
                let mut r = self.rows[self.emitted].clone();
                self.fill_cf(&mut r);
                self.rows[self.emitted].cf_w = r.cf_w.clone();
                self.emitted += 1;
                if self.emitted == self.rows.len() && self.phase == Phase::Emit {
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
