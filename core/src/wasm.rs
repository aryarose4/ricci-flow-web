//! JS-facing streaming session (see `core/API.md`).
use crate::flow::Seed;
use crate::pipeline::{FrameOut, Session, SessionConfig};
use crate::seedmesh::page_level;
use wasm_bindgen::prelude::*;

fn flat3(v: &[[f64; 3]]) -> Vec<f64> {
    v.iter().flatten().cloned().collect()
}
fn flatf(v: &[[usize; 3]]) -> Vec<u32> {
    v.iter().flatten().map(|&x| x as u32).collect()
}

/// `FlowSession(positions, betas, seed, level)`: positions flat `[re0, im0, re1, im1, ...]` (cone
/// points in the complex plane), `betas` = cone angle / 2 pi in [0.05, 0.99], seed
/// "gluing" | "smart" | "polytope", level "smoke" | "coarse" | "medium" | "fine" | "extra".
/// Throws (a string) for rho <= 0, bad betas, or a mesh failure.
#[wasm_bindgen]
pub struct FlowSession {
    s: Session,
}

#[wasm_bindgen]
impl FlowSession {
    #[wasm_bindgen(constructor)]
    pub fn new(positions: &[f64], betas: &[f64], seed: &str, level: &str) -> Result<FlowSession, JsValue> {
        let seed = Seed::parse(seed).ok_or_else(|| JsValue::from_str("seed must be gluing | smart | polytope"))?;
        let level = page_level(level).ok_or_else(|| JsValue::from_str("unknown mesh level"))?;
        if positions.len() != 2 * betas.len() {
            return Err(JsValue::from_str("positions must hold 2 numbers per beta"));
        }
        let positions = positions.chunks(2).map(|c| [c[0], c[1]]).collect();
        let s = Session::new(SessionConfig { positions, betas: betas.to_vec(), seed, level }).map_err(|e| JsValue::from_str(&e))?;
        Ok(FlowSession { s })
    }
    /// One work chunk; returns the next frame once realization is done, else undefined.
    pub fn next_frame(&mut self) -> Result<Option<FlowFrame>, JsValue> {
        Ok(self.s.next_frame().map_err(|e| JsValue::from_str(&e))?.map(|f| FlowFrame { f }))
    }
    /// One work chunk without taking frames; false once everything is computed.
    pub fn step(&mut self) -> Result<bool, JsValue> {
        self.s.step().map_err(|e| JsValue::from_str(&e))
    }
    pub fn done(&self) -> bool {
        self.s.is_done()
    }
    pub fn progress(&self) -> f64 {
        self.s.progress()
    }
    pub fn phase(&self) -> String {
        self.s.phase().name().to_string()
    }
    pub fn n_frames(&self) -> usize {
        self.s.n_frames()
    }
    /// rho of the uniformized metric (NaN until the intrinsic flow is done).
    pub fn rho(&self) -> f64 {
        self.s.rho()
    }
    /// Gauss-Bonnet constant C = 4 pi - sum of cone defects (> 0).
    pub fn curvature_total(&self) -> f64 {
        self.s.c_curv
    }
    pub fn n_vertices(&self) -> usize {
        self.s.mesh.nv()
    }
    /// unit-sphere mesh vertices, flat nv*3
    pub fn sphere_vertices(&self) -> Vec<f64> {
        flat3(&self.s.mesh.y)
    }
    /// preview polyhedron (sphere with cone caps), flat nv*3
    pub fn preview_vertices(&self) -> Vec<f64> {
        flat3(&self.s.flow_mesh.vertices)
    }
    /// outward faces, flat nf*3
    pub fn faces(&self) -> Vec<u32> {
        flatf(&self.s.mesh.f)
    }
    /// edges (a < b), flat ne*2, order of `FlowFrame.lengths`
    pub fn edges(&self) -> Vec<u32> {
        self.s.mesh.e.iter().flatten().map(|&x| x as u32).collect()
    }
    /// cone vertex ids (cone i = positions[i])
    pub fn cones(&self) -> Vec<u32> {
        self.s.mesh.cones.iter().map(|&x| x as u32).collect()
    }
    /// Polytope seed's 4 extra corner vertex ids (empty otherwise)
    pub fn seed_corners(&self) -> Vec<u32> {
        self.s.seed_corners.iter().map(|&x| x as u32).collect()
    }
    /// seed betas actually used for the cap mesh (differ from targets if caps were clamped)
    pub fn seed_betas(&self) -> Vec<f64> {
        self.s.flow_mesh.seed.seed_betas.clone()
    }
    /// true if overlapping caps were shrunk
    pub fn clamped(&self) -> bool {
        !self.s.flow_mesh.seed.overlaps.is_empty()
    }
    /// The t = 0 start shape ("seed" event), available after the clean-start phase.
    pub fn seed_frame(&self) -> Option<FlowFrame> {
        self.s.seed_preview().cloned().map(|f| FlowFrame { f })
    }
    pub fn quality_label(&self) -> Option<String> {
        self.s.quality.as_ref().map(|q| q.label.clone())
    }
    pub fn quality_detail(&self) -> Option<String> {
        self.s.quality.as_ref().map(|q| q.detail.clone())
    }
    pub fn quality_ok(&self) -> Option<bool> {
        self.s.quality.as_ref().map(|q| q.ok)
    }
    pub fn warning(&self) -> Option<String> {
        self.s.warning.clone()
    }
    /// [steps, newton, rejected, newton_fail] of the intrinsic flow (frames only)
    pub fn flow_stats(&self) -> Vec<u32> {
        self.s.stats_flow.as_ref().map_or(vec![], |s| vec![s.steps as u32, s.newton as u32, s.rejected as u32, s.newton_fail as u32])
    }
}

/// One output frame.
#[wasm_bindgen]
pub struct FlowFrame {
    f: FrameOut,
}

#[wasm_bindgen]
impl FlowFrame {
    pub fn k(&self) -> usize {
        self.f.k
    }
    pub fn t(&self) -> f64 {
        self.f.t
    }
    /// log conformal factor per vertex
    pub fn u(&self) -> Vec<f64> {
        self.f.u.clone()
    }
    /// angle defect per vertex (K = k_defect / a_vertex)
    pub fn k_defect(&self) -> Vec<f64> {
        self.f.k_defect.clone()
    }
    pub fn a_vertex(&self) -> Vec<f64> {
        self.f.a_vertex.clone()
    }
    /// intrinsic edge lengths (session edge order)
    pub fn lengths(&self) -> Vec<f64> {
        self.f.lengths.clone()
    }
    /// 3D vertices, flat nv*3
    pub fn vertices(&self) -> Vec<f64> {
        flat3(&self.f.vertices)
    }
    /// faces when they differ from the session mesh (Polytope frame 0), else undefined
    pub fn faces(&self) -> Option<Vec<u32>> {
        self.f.faces.as_ref().map(|f| flatf(f))
    }
    pub fn len_rel_max(&self) -> f64 {
        self.f.len_rel_max
    }
    pub fn len_rel_rms(&self) -> f64 {
        self.f.len_rel_rms
    }
    pub fn n_deep(&self) -> usize {
        self.f.n_deep
    }
    pub fn n_reflex(&self) -> usize {
        self.f.n_reflex
    }
    pub fn min_dihedral(&self) -> f64 {
        self.f.min_dihedral
    }
    pub fn gb_resid(&self) -> f64 {
        self.f.gb_resid
    }
    pub fn cone_defect_err(&self) -> f64 {
        self.f.cone_defect_err
    }
    pub fn rho(&self) -> f64 {
        self.f.rho
    }
    pub fn k_dev_max(&self) -> f64 {
        self.f.k_dev_max
    }
    pub fn seed_folds(&self) -> bool {
        self.f.seed_folds
    }
    pub fn polytope(&self) -> bool {
        self.f.polytope
    }
    pub fn warning(&self) -> Option<String> {
        self.f.warning.clone()
    }
    /// conformal-factor raster, 64 x 128 row-major (rows lat 0..pi from north, cols lon -pi..pi)
    pub fn cf_w(&self) -> Vec<f64> {
        self.f.cf_w.clone()
    }
    pub fn cf_lat(&self) -> usize {
        crate::cfgrid::GRID_LAT
    }
    pub fn cf_lon(&self) -> usize {
        crate::cfgrid::GRID_LON
    }
}
