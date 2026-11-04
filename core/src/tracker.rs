//! rho < 0 True Ricci Flow: "can the 3D embedding still follow the flow?"  Port of
//! `viewer/backend/trf_track.py` (`Tracker`, `mismatch`, `cold_start`) and of
//! `flow_true.embed_notice`.
//!
//! The intrinsic flow always runs; each frame is realized in 3D warm from the previous one
//! (continuation).  [`mismatch`] measures the realized surface against the flow metric:
//! relative edge-length error (rms / max) and the ball-averaged (3-ring) angle-defect curvature of
//! the 3D surface minus the same estimator on the flow's lengths, divided by |rho| (signed mean and
//! rms over non-cone vertices outside the cones' 1-rings).  Self-intersections are not computed (they
//! are off the failure criterion by default in Python).

use crate::realize::geom::{self, V3};
use crate::realize::{FastRealizer, Realized};
use std::f64::consts::PI;

pub const LEN_RMS_FAIL: f64 = 5e-3;
pub const K_SIGNED_FAIL: f64 = 0.10;
pub const CONFIRM: usize = 2;
pub const SEED_SCALES: [f64; 4] = [1.0, 0.8, 0.9, 0.7];
pub const K_RING: usize = 3;

#[derive(Clone, Copy, Debug, Default)]
pub struct Mismatch {
    pub len_rel_rms: f64,
    pub len_rel_max: f64,
    pub k_err_signed: f64,
    pub k_err_rms: f64,
}

/// "searching" | "following" | "suspect" | "failed" | "lost"
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum State {
    Searching,
    Following,
    Suspect,
    Failed,
    Lost,
}

impl State {
    pub fn name(&self) -> &'static str {
        match self {
            State::Searching => "searching",
            State::Following => "following",
            State::Suspect => "suspect",
            State::Failed => "failed",
            State::Lost => "lost",
        }
    }
    /// Python `embed_lost`: the shown vertices are not a passing fit of this frame.
    pub fn lost(&self) -> bool {
        matches!(self, State::Searching | State::Failed | State::Lost)
    }
}

#[derive(Clone, Debug, Default)]
pub struct Summary {
    pub failed_t: Option<f64>,
    pub reason: Option<String>,
    pub mismatch: Option<Mismatch>,
    pub last_good_t: Option<f64>,
    pub first_good_t: Option<f64>,
}

pub struct Tracker {
    pub faces: Vec<[usize; 3]>,
    pub nv: usize,
    /// k-ring balls (sorted vertex lists, self included)
    balls: Vec<Vec<usize>>,
    /// 1 off the cones and their 1-ring
    w: Vec<f64>,
    pub failed_t: Option<f64>,
    pub fail_reason: Option<String>,
    pub fail_mismatch: Option<Mismatch>,
    pub last_good_t: Option<f64>,
    pub first_good_t: Option<f64>,
    best_search: Option<Mismatch>,
    pub last_reasons: Vec<&'static str>,
    pending: Option<(f64, Mismatch, Vec<&'static str>, usize)>,
}

fn adjacency(faces: &[[usize; 3]], nv: usize) -> Vec<Vec<usize>> {
    let mut a = vec![vec![]; nv];
    for f in faces {
        for i in 0..3 {
            let (p, q) = (f[i], f[(i + 1) % 3]);
            a[p].push(q);
            a[q].push(p);
        }
    }
    for l in &mut a {
        l.sort_unstable();
        l.dedup();
    }
    a
}

/// `flow_curvature.defect_and_area`: (2 pi - sum of corner angles, barycentric area) per vertex.
pub fn defect_and_area(x: &[V3], faces: &[[usize; 3]], nv: usize) -> (Vec<f64>, Vec<f64>) {
    let mut d = vec![2.0 * PI; nv];
    let mut a = vec![0.0; nv];
    for f in faces {
        for i in 0..3 {
            let (p, q, r) = (f[i], f[(i + 1) % 3], f[(i + 2) % 3]);
            let (u, v) = (geom::sub(x[q], x[p]), geom::sub(x[r], x[p]));
            d[p] -= geom::norm(geom::cross(u, v)).atan2(geom::dot(u, v));
        }
        let ar = 0.5 * geom::norm(geom::cross(geom::sub(x[f[1]], x[f[0]]), geom::sub(x[f[2]], x[f[0]])));
        for &p in f {
            a[p] += ar / 3.0;
        }
    }
    (d, a)
}

impl Tracker {
    pub fn new(faces: &[[usize; 3]], cones: &[usize]) -> Tracker {
        let nv = faces.iter().flatten().max().map_or(0, |m| m + 1);
        let adj = adjacency(faces, nv);
        let mut balls = Vec::with_capacity(nv);
        let mut mark = vec![usize::MAX; nv];
        for v in 0..nv {
            let mut ball = vec![v];
            mark[v] = v;
            let mut front = vec![v];
            for _ in 0..K_RING {
                let mut nf = vec![];
                for &p in &front {
                    for &q in &adj[p] {
                        if mark[q] != v {
                            mark[q] = v;
                            ball.push(q);
                            nf.push(q);
                        }
                    }
                }
                front = nf;
            }
            ball.sort_unstable();
            balls.push(ball);
        }
        let mut w = vec![1.0; nv];
        for &c in cones {
            w[c] = 0.0;
            for &q in &adj[c] {
                w[q] = 0.0;
            }
        }
        Tracker {
            faces: faces.to_vec(),
            nv,
            balls,
            w,
            failed_t: None,
            fail_reason: None,
            fail_mismatch: None,
            last_good_t: None,
            first_good_t: None,
            best_search: None,
            last_reasons: vec![],
            pending: None,
        }
    }

    /// Ball-averaged K (NaN on excluded vertices / empty balls).
    pub fn ball_k(&self, defect: &[f64], area: &[f64]) -> Vec<f64> {
        (0..self.nv)
            .map(|v| {
                let (mut num, mut den) = (0.0, 0.0);
                for &u in &self.balls[v] {
                    num += defect[u] * self.w[u];
                    den += area[u] * self.w[u];
                }
                if den > 0.0 && self.w[v] > 0.0 {
                    num / den
                } else {
                    f64::NAN
                }
            })
            .collect()
    }

    pub fn reasons(&self, m: &Mismatch, ok: bool) -> Vec<&'static str> {
        let mut r = vec![];
        if !ok {
            r.push("solver did not converge");
        }
        if m.len_rel_rms > LEN_RMS_FAIL {
            r.push("edge-length error above 0.5 %");
        }
        if m.k_err_signed.abs() > K_SIGNED_FAIL {
            r.push("curvature error above 10 %");
        }
        r
    }

    pub fn update(&mut self, t: f64, m: &Mismatch, ok: bool) -> State {
        if self.failed_t.is_some() {
            return State::Lost;
        }
        let r = self.reasons(m, ok);
        self.last_reasons = r.clone();
        if !r.is_empty() && self.last_good_t.is_none() {
            if self.best_search.map_or(true, |b| m.len_rel_rms < b.len_rel_rms) {
                self.best_search = Some(*m);
            }
            return State::Searching;
        }
        if !r.is_empty() {
            let p = self.pending.get_or_insert((t, *m, r, 0));
            p.3 += 1;
            if p.3 >= CONFIRM {
                self.latch();
                return State::Failed;
            }
            return State::Suspect;
        }
        self.pending = None;
        if self.first_good_t.is_none() {
            self.first_good_t = Some(t);
        }
        self.last_good_t = Some(t);
        State::Following
    }

    fn latch(&mut self) {
        let (t, m, r, _) = self.pending.take().unwrap();
        self.failed_t = Some(t);
        self.fail_reason = Some(r.join("; "));
        self.fail_mismatch = Some(m);
    }

    /// End of the run: an unconfirmed bad last frame counts as the failure; if nothing ever passed,
    /// the embedding failed from t = 0.
    pub fn finish(&mut self) -> Summary {
        if self.failed_t.is_none() && self.pending.is_some() {
            self.latch();
        }
        if self.failed_t.is_none() && self.last_good_t.is_none() && self.best_search.is_some() {
            self.failed_t = Some(0.0);
            self.fail_mismatch = self.best_search;
            self.fail_reason = Some(format!("never started: {}", self.last_reasons.join("; ")));
        }
        self.summary()
    }

    pub fn summary(&self) -> Summary {
        Summary {
            failed_t: self.failed_t,
            reason: self.fail_reason.clone(),
            mismatch: self.fail_mismatch,
            last_good_t: self.last_good_t,
            first_good_t: self.first_good_t,
        }
    }
}

/// `trf_track.mismatch(X, frame, tr)`; `k_defect` / `a_vertex` are the flow frame's intrinsic data.
pub fn mismatch(tr: &Tracker, x: &[V3], e: &[[usize; 2]], l: &[f64], k_defect: &[f64], a_vertex: &[f64], rho: f64) -> Mismatch {
    let rel = geom::len_rel(x, e, l);
    let (d3, a3) = defect_and_area(x, &tr.faces, tr.nv);
    let (k3, kf) = (tr.ball_k(&d3, &a3), tr.ball_k(k_defect, a_vertex));
    let rho = if rho.abs() == 0.0 { 1.0 } else { rho.abs() };
    let dk: Vec<f64> = k3.iter().zip(&kf).map(|(a, b)| (a - b) / rho).filter(|v| v.is_finite()).collect();
    let n = dk.len() as f64;
    Mismatch {
        len_rel_rms: (rel.iter().map(|v| v * v).sum::<f64>() / rel.len() as f64).sqrt(),
        len_rel_max: rel.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        k_err_signed: if dk.is_empty() { 0.0 } else { dk.iter().sum::<f64>() / n },
        k_err_rms: if dk.is_empty() { 0.0 } else { (dk.iter().map(|v| v * v).sum::<f64>() / n).sqrt() },
    }
}

/// `trf_track.cold_start`: Fast realization from the cone-cap seed scaled by each of
/// [`SEED_SCALES`]; the first passing attempt, else the lowest-rms one.
pub fn cold_start(rlz: &FastRealizer, tr: &Tracker, l: &[f64], kd: &[f64], av: &[f64], rho: f64) -> (Realized, Mismatch) {
    let mut best: Option<(Realized, Mismatch)> = None;
    for s in SEED_SCALES {
        let x0: Vec<V3> = rlz.x0.iter().map(|p| geom::scale(*p, s)).collect();
        let rf = rlz.realize_cold_from(l, &x0);
        let m = mismatch(tr, &rf.x, &rlz.mesh.edges, l, kd, av, rho);
        if tr.reasons(&m, rf.ok).is_empty() {
            return (rf, m);
        }
        if best.as_ref().map_or(true, |b| m.len_rel_rms < b.1.len_rel_rms) {
            best = Some((rf, m));
        }
    }
    best.unwrap()
}

/// Python `"%.<p>g" % x`.
pub fn fmt_g(x: f64, p: usize) -> String {
    if !x.is_finite() {
        return if x.is_nan() { "nan".into() } else if x > 0.0 { "inf".into() } else { "-inf".into() };
    }
    if x == 0.0 {
        return "0".into();
    }
    let p = p.max(1);
    let e = format!("{:.*e}", p - 1, x);
    let exp: i32 = e[e.find('e').unwrap() + 1..].parse().unwrap();
    let strip = |s: String| -> String {
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s
        }
    };
    if exp < -4 || exp >= p as i32 {
        let (mant, _) = e.split_at(e.find('e').unwrap());
        format!("{}e{}{:02}", strip(mant.to_string()), if exp < 0 { '-' } else { '+' }, exp.abs())
    } else {
        strip(format!("{:.*}", (p as i32 - 1 - exp).max(0) as usize, x))
    }
}

/// `flow_true.embed_notice(fail)`: banner text once tracking has failed (None while following).
pub fn embed_notice(fail: &Summary) -> Option<String> {
    let ft = fail.failed_t?;
    let m = fail.mismatch.unwrap_or(Mismatch { len_rel_rms: f64::NAN, len_rel_max: f64::NAN, k_err_signed: f64::NAN, k_err_rms: f64::NAN });
    let reason = fail.reason.clone().unwrap_or_else(|| "None".into());
    let ks = 100.0 * m.k_err_signed;
    let kstr = if ks.is_nan() { "+nan".to_string() } else { format!("{:+.0}", ks) };
    let kstr = if kstr == "-0" { "-0".to_string() } else { kstr };
    if fail.last_good_t.is_none() {
        return Some(format!(
            "The 3D shape could not follow the flow at any time (t = 0 onward); the best attempt is {}% off in edge lengths and {}% off in curvature ({}).",
            fmt_g(100.0 * m.len_rel_rms, 2), kstr, reason
        ));
    }
    Some(format!(
        "The 3D shape can no longer follow the flow from t = {}; it is {}% off in edge lengths and {}% off in curvature ({}).",
        fmt_g(ft, 3), fmt_g(100.0 * m.len_rel_rms, 2), kstr, reason
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g_format() {
        assert_eq!(fmt_g(0.330448899097056, 3), "0.33");
        assert_eq!(fmt_g(0.6771189071109298, 2), "0.68");
        assert_eq!(fmt_g(16.0, 3), "16");
        assert_eq!(fmt_g(1.25e-5, 2), "1.3e-05");
        assert_eq!(fmt_g(0.125, 2), "0.12");
        assert_eq!(fmt_g(2.675, 3), "2.67");
        assert_eq!(fmt_g(123.4, 2), "1.2e+02");
        assert_eq!(format!("{:+.0}", 4.49), "+4");
    }
}
