//! 3D realization of the /flow True Ricci Flow frames (rho > 0 path).
//!
//! Ports of `singular_ricci/dev/embeddings/true_flow.py`: `fast_seed`, `FastRealizer`,
//! `n_deep_folds`, `clean_start`, `realize_backward`, `polytope_snap`, `polytope_fit`, and of
//! `viewer/backend/flow_true.py`: `polytope_vertices`, `hull_faces`; plus the helpers they use
//! (`seed_on_mesh.build_seed_on_mesh`, `convex_alexandrov._dihedrals`, `convex_bi.kabsch_align`).
//!
//! Intrinsic inputs are plain: faces `Vec<[usize;3]>`, edges `Vec<[usize;2]>` (a < b), per-frame
//! edge lengths `Vec<f64>`.  The flow itself is reached through [`FlowSource`] (in-between metrics
//! for hidden substeps and the clean-start blend), so this module does not depend on the flow crate
//! module.  The Bobenko-Izmestiev convex anchor (`convex_anchor`) is also behind [`FlowSource`]
//! (`convex_anchor` returns None by default: BI is not ported; it is only consulted when the
//! cold fit of the uniformized end folds).

use crate::hull;
use crate::lm::{self, FadeSetup, LsqProblem};
use std::cell::RefCell;

pub mod geom {
    pub type V3 = [f64; 3];
    #[inline]
    pub fn sub(a: V3, b: V3) -> V3 {
        [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    }
    #[inline]
    pub fn add(a: V3, b: V3) -> V3 {
        [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
    }
    #[inline]
    pub fn scale(a: V3, s: f64) -> V3 {
        [a[0] * s, a[1] * s, a[2] * s]
    }
    #[inline]
    pub fn dot(a: V3, b: V3) -> f64 {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }
    #[inline]
    pub fn cross(a: V3, b: V3) -> V3 {
        [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
    }
    #[inline]
    pub fn norm(a: V3) -> f64 {
        dot(a, a).sqrt()
    }

    /// `embed4._signed_volume`.
    pub fn signed_volume(x: &[V3], f: &[[usize; 3]]) -> f64 {
        f.iter().map(|t| dot(x[t[0]], cross(x[t[1]], x[t[2]]))).sum::<f64>() / 6.0
    }

    /// |  |X_a - X_b| / L - 1 | per edge.
    pub fn len_rel(x: &[V3], e: &[[usize; 2]], l: &[f64]) -> Vec<f64> {
        e.iter().zip(l).map(|(e, &l)| (norm(sub(x[e[0]], x[e[1]])) / l - 1.0).abs()).collect()
    }

    /// `thick_part.edges_of`: sorted unique (a < b) edges.
    pub fn edges_of(f: &[[usize; 3]]) -> Vec<[usize; 2]> {
        let mut e: Vec<[usize; 2]> = f
            .iter()
            .flat_map(|t| [[t[0], t[1]], [t[1], t[2]], [t[2], t[0]]])
            .map(|[a, b]| if a < b { [a, b] } else { [b, a] })
            .collect();
        e.sort_unstable();
        e.dedup();
        e
    }

    /// `embed4.face_edge_table`: interior edges (a < b) with their two faces.
    pub fn face_edge_table(f: &[[usize; 3]]) -> Vec<([usize; 2], [usize; 2])> {
        let mut inc: std::collections::BTreeMap<(usize, usize), Vec<usize>> = Default::default();
        for (fi, t) in f.iter().enumerate() {
            for k in 0..3 {
                let (u, v) = (t[k], t[(k + 1) % 3]);
                inc.entry((u.min(v), u.max(v))).or_default().push(fi);
            }
        }
        inc.into_iter().filter(|(_, fl)| fl.len() == 2).map(|((u, v), fl)| ([u, v], [fl[0], fl[1]])).collect()
    }

    fn traverses(t: [usize; 3], a: usize, b: usize) -> bool {
        (0..3).any(|k| t[k] == a && t[(k + 1) % 3] == b)
    }

    /// `bend_gn.edge_quads`: (edge (a<b), quad (a', b', c, d)) with face (a', b', c) traversing a'->b'.
    pub fn edge_quads(f: &[[usize; 3]]) -> Vec<([usize; 2], [usize; 4])> {
        face_edge_table(f)
            .into_iter()
            .map(|([a, b], [f0, f1])| {
                let t0 = f[f0];
                let t1 = f[f1];
                let c = *t0.iter().find(|&&v| v != a && v != b).unwrap();
                let d = *t1.iter().find(|&&v| v != a && v != b).unwrap();
                let q = if traverses(t0, a, b) { [a, b, c, d] } else { [b, a, c, d] };
                ([a, b], q)
            })
            .collect()
    }

    /// `convex_alexandrov._dihedrals(X, F)[2]` (the `embed5._dihed_engine` value, convex > 0 for
    /// outward orientation, sign fixed by the signed volume).  Note the engine's value is
    /// atan2(|e| sin(theta), cos(theta)) (unit normals, unnormalised edge vector); it is what the
    /// fold thresholds (`n_deep_folds`, `n_reflex`) are applied to, so it is reproduced as is.
    pub fn dihedrals(x: &[V3], f: &[[usize; 3]]) -> Vec<f64> {
        let s = if signed_volume(x, f) > 0.0 { 1.0 } else { -1.0 };
        face_edge_table(f)
            .into_iter()
            .map(|([a, b], [e0, e1])| {
                let (f1, f2) = if traverses(f[e0], a, b) { (e0, e1) } else { (e1, e0) };
                let (t1, t2) = (f[f1], f[f2]);
                let mut n1 = cross(sub(x[t1[1]], x[t1[0]]), sub(x[t1[2]], x[t1[0]]));
                let mut n2 = cross(sub(x[t2[1]], x[t2[0]]), sub(x[t2[2]], x[t2[0]]));
                n1 = scale(n1, 1.0 / norm(n1));
                n2 = scale(n2, 1.0 / norm(n2));
                let dv = sub(x[b], x[a]);
                s * dot(cross(n1, n2), dv).atan2(dot(n1, n2))
            })
            .collect()
    }

    /// `convex_bi.kabsch_align(X, Xref)`: proper rigid motion of X best matching Xref
    /// (Horn's quaternion method = Kabsch with the reflection correction).
    pub fn kabsch_align(x: &[V3], xr: &[V3]) -> Vec<V3> {
        let n = x.len() as f64;
        let mut cx = [0.0; 3];
        let mut cr = [0.0; 3];
        for (a, b) in x.iter().zip(xr) {
            cx = add(cx, *a);
            cr = add(cr, *b);
        }
        cx = scale(cx, 1.0 / n);
        cr = scale(cr, 1.0 / n);
        let mut h = [[0.0; 3]; 3];
        for (a, b) in x.iter().zip(xr) {
            let p = sub(*a, cx);
            let q = sub(*b, cr);
            for i in 0..3 {
                for j in 0..3 {
                    h[i][j] += p[i] * q[j];
                }
            }
        }
        let (sxx, sxy, sxz) = (h[0][0], h[0][1], h[0][2]);
        let (syx, syy, syz) = (h[1][0], h[1][1], h[1][2]);
        let (szx, szy, szz) = (h[2][0], h[2][1], h[2][2]);
        let nm = [
            [sxx + syy + szz, syz - szy, szx - sxz, sxy - syx],
            [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
            [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
            [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz],
        ];
        let q = top_eigvec4(nm);
        let (w, qx, qy, qz) = (q[0], q[1], q[2], q[3]);
        let r = [
            [w * w + qx * qx - qy * qy - qz * qz, 2.0 * (qx * qy - w * qz), 2.0 * (qx * qz + w * qy)],
            [2.0 * (qx * qy + w * qz), w * w - qx * qx + qy * qy - qz * qz, 2.0 * (qy * qz - w * qx)],
            [2.0 * (qx * qz - w * qy), 2.0 * (qy * qz + w * qx), w * w - qx * qx - qy * qy + qz * qz],
        ];
        x.iter()
            .map(|a| {
                let p = sub(*a, cx);
                add([dot(r[0], p), dot(r[1], p), dot(r[2], p)], cr)
            })
            .collect()
    }

    /// Eigenvector of the largest eigenvalue of a symmetric 4x4 matrix (cyclic Jacobi).
    fn top_eigvec4(mut a: [[f64; 4]; 4]) -> [f64; 4] {
        let mut v = [[0.0; 4]; 4];
        for i in 0..4 {
            v[i][i] = 1.0;
        }
        for _ in 0..100 {
            let off: f64 = (0..4).flat_map(|i| (0..4).filter(move |&j| j != i).map(move |j| (i, j))).map(|(i, j)| a[i][j] * a[i][j]).sum();
            let tot: f64 = (0..4).map(|i| a[i][i] * a[i][i]).sum::<f64>() + off;
            if off <= 1e-30 * tot.max(1e-300) {
                break;
            }
            for p in 0..4 {
                for q in p + 1..4 {
                    if a[p][q].abs() < 1e-300 {
                        continue;
                    }
                    let th = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                    let t = th.signum() / (th.abs() + (th * th + 1.0).sqrt());
                    let t = if th == 0.0 { 1.0 } else { t };
                    let c = 1.0 / (t * t + 1.0).sqrt();
                    let s = t * c;
                    for k in 0..4 {
                        let (akp, akq) = (a[k][p], a[k][q]);
                        a[k][p] = c * akp - s * akq;
                        a[k][q] = s * akp + c * akq;
                    }
                    for k in 0..4 {
                        let (apk, aqk) = (a[p][k], a[q][k]);
                        a[p][k] = c * apk - s * aqk;
                        a[q][k] = s * apk + c * aqk;
                    }
                    for k in 0..4 {
                        let (vkp, vkq) = (v[k][p], v[k][q]);
                        v[k][p] = c * vkp - s * vkq;
                        v[k][q] = s * vkp + c * vkq;
                    }
                }
            }
        }
        let i = (0..4).max_by(|&i, &j| a[i][i].partial_cmp(&a[j][j]).unwrap()).unwrap();
        [v[0][i], v[1][i], v[2][i], v[3][i]]
    }
}

use geom::V3;

// ------------------------------------------------------------------ seed
/// Where the cold-start seed came from (`_Seed.seed_kind`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SeedKind {
    Caps,
    Sphere,
}

fn angle(p: V3, q: V3) -> f64 {
    geom::norm(geom::cross(p, q)).atan2(geom::dot(p, q))
}

fn tip_angle_sum(v: &[V3], ring: &[[usize; 3]], t: usize) -> f64 {
    let mut s = 0.0;
    for f in ring {
        let k = (0..3).find(|&k| f[k] == t).unwrap();
        let a = v[f[k]];
        let u = geom::sub(v[f[(k + 1) % 3]], a);
        let w = geom::sub(v[f[(k + 2) % 3]], a);
        s += geom::norm(geom::cross(u, w)).atan2(geom::dot(u, w));
    }
    s
}

/// Root of f on [a, b] with f(a) > 0 >= f(b) (bisection to machine precision; stands in for brentq
/// with xtol = rtol = 1e-15).
fn bisect(f: &mut dyn FnMut(f64) -> f64, mut a: f64, mut b: f64) -> f64 {
    for _ in 0..200 {
        let m = 0.5 * (a + b);
        if m <= a || m >= b {
            break;
        }
        if f(m) > 0.0 {
            a = m;
        } else {
            b = m;
        }
        if (b - a) <= 1e-15 * (1.0 + b.abs()) {
            break;
        }
    }
    0.5 * (a + b)
}

fn solve_tip(v: &[V3], ring: &[[usize; 3]], t: usize, p: V3, target: f64, h0: Option<f64>) -> Option<f64> {
    let mut w = v.to_vec();
    let mut f = |h: f64| {
        w[t] = geom::scale(p, h);
        tip_angle_sum(&w, ring, t) - target
    };
    if let Some(h0) = h0 {
        let mut d = 1e-3;
        for _ in 0..12 {
            let (a, b) = ((h0 - d).max(0.05), h0 + d);
            if f(a) > 0.0 && 0.0 >= f(b) {
                return Some(bisect(&mut f, a, b));
            }
            d *= 2.0;
        }
    }
    let mut h_ring = f64::MIN;
    for fc in ring {
        for &x in fc {
            if x != t {
                h_ring = h_ring.max(geom::dot(v[x], p));
            }
        }
    }
    let top = h_ring.max(1.0) + 4.0;
    let s0 = h_ring.max(0.05) + 1e-9;
    let scan = |hs: &mut dyn Iterator<Item = f64>, f: &mut dyn FnMut(f64) -> f64| {
        let mut br = None;
        let mut prev: Option<(f64, f64)> = None;
        for h in hs {
            let fv = f(h);
            if let Some((ph, pv)) = prev {
                if pv > 0.0 && 0.0 >= fv {
                    br = Some((ph, h));
                }
            }
            prev = Some((h, fv));
        }
        br
    };
    let lin = |a: f64, b: f64, n: usize| (0..n).map(move |i| a + (b - a) * i as f64 / (n - 1) as f64);
    let mut br = scan(&mut lin(s0, top, 400), &mut f);
    if br.is_none() {
        br = scan(&mut lin(0.05, h_ring.max(0.05), 100), &mut f);
    }
    let (lo, hi) = br?;
    if f(hi) == 0.0 {
        return Some(hi);
    }
    Some(bisect(&mut f, lo, hi))
}

/// `seed_on_mesh.build_seed_on_mesh(faces, xyz, cones, betas)["vertices"]`: the sphere mesh with
/// convex cone caps (tangent cones of angular radius arccos(beta), clamped if caps overlap) whose
/// tips are solved for the exact cone angle.
pub fn build_seed_on_mesh(faces: &[[usize; 3]], y: &[V3], cones: &[usize], betas: &[f64]) -> Option<Vec<V3>> {
    if betas.iter().any(|&b| b >= 1.0 || b <= 0.0) {
        return None;
    }
    let v0: Vec<V3> = y.iter().map(|p| geom::scale(*p, 1.0 / geom::norm(*p))).collect();
    let pp: Vec<V3> = cones.iter().map(|&c| v0[c]).collect();
    let a_t: Vec<f64> = betas.iter().map(|b| b.acos()).collect();
    let nc = cones.len();
    let over = (0..nc).any(|i| (i + 1..nc).any(|j| angle(pp[i], pp[j]) <= a_t[i] + a_t[j]));
    let mut a_u = a_t.clone();
    if over {
        let gap = 0.08;
        for _ in 0..200 {
            let mut changed = false;
            for i in 0..nc {
                for j in i + 1..nc {
                    let d = angle(pp[i], pp[j]);
                    if a_u[i] + a_u[j] > d - gap + 1e-12 {
                        let s = (d - gap).max(1e-3) / (a_u[i] + a_u[j]);
                        a_u[i] *= s;
                        a_u[j] *= s;
                        changed = true;
                    }
                }
            }
            if !changed {
                break;
            }
        }
        for a in a_u.iter_mut() {
            *a = a.max(1e-3);
        }
    }
    let mut v = v0.clone();
    let mut is_tip = vec![false; v.len()];
    for &c in cones {
        is_tip[c] = true;
    }
    for i in 0..nc {
        for k in 0..v.len() {
            if is_tip[k] {
                continue;
            }
            let th = angle(v0[k], pp[i]);
            if th < a_u[i] {
                v[k] = geom::scale(v0[k], 1.0 / (a_u[i] - th).cos());
            }
        }
    }
    let rings: Vec<Vec<[usize; 3]>> = cones.iter().map(|&t| faces.iter().filter(|f| f.contains(&t)).copied().collect()).collect();
    let targets: Vec<f64> = a_u.iter().map(|a| 2.0 * std::f64::consts::PI * a.cos()).collect();
    let mut heights: Vec<f64> = cones.iter().map(|&t| geom::norm(v[t])).collect();
    for sweep in 1..=100 {
        for i in 0..nc {
            let h = solve_tip(&v, &rings[i], cones[i], pp[i], targets[i], if sweep > 1 { Some(heights[i]) } else { None })?;
            v[cones[i]] = geom::scale(pp[i], h);
            heights[i] = h;
        }
        let err = (0..nc).map(|i| (tip_angle_sum(&v, &rings[i], cones[i]) - targets[i]).abs()).fold(0.0, f64::max);
        if err <= 0.1 * 1e-12 {
            break;
        }
    }
    Some(v)
}

/// `true_flow.fast_seed`: cap seed, or the bare sphere (`SeedKind::Sphere`) if impossible.
pub fn fast_seed(y: &[V3], faces: &[[usize; 3]], cones: &[usize], betas: &[f64]) -> (Vec<V3>, SeedKind) {
    match build_seed_on_mesh(faces, y, cones, betas) {
        Some(x) if x.len() == y.len() && x.iter().all(|p| p.iter().all(|c| c.is_finite())) => (x, SeedKind::Caps),
        _ => (y.to_vec(), SeedKind::Sphere),
    }
}

// ------------------------------------------------------------------ realizer
/// Triangulation shared by every frame (faces, edges a < b, cone vertex ids).
#[derive(Clone, Debug)]
pub struct RealizeMesh {
    pub faces: Vec<[usize; 3]>,
    pub edges: Vec<[usize; 2]>,
    pub cones: Vec<usize>,
    pub nv: usize,
}

impl RealizeMesh {
    pub fn new(faces: Vec<[usize; 3]>, edges: Vec<[usize; 2]>, cones: Vec<usize>, nv: usize) -> Self {
        RealizeMesh { faces, edges, cones, nv }
    }
}

/// One realized frame (`RealizedFrame` with the FastRealizer diag fields).
#[derive(Clone, Debug)]
pub struct Realized {
    pub x: Vec<V3>,
    pub ok: bool,
    pub len_rel_max: f64,
    pub len_rel_rms: f64,
    pub n_reflex: usize,
    pub min_dihedral: f64,
    pub lm_iters: usize,
    /// "cold" | "warm" | "warm_full" (+ "+polish") | "clean_start"
    pub start: String,
    pub warning: Option<String>,
}

/// One intrinsic frame: index, time, edge lengths.
#[derive(Clone, Debug)]
pub struct MetricFrame {
    pub k: usize,
    pub t: f64,
    pub lengths: Vec<f64>,
}

/// The flow as seen by the realization driver.
pub trait FlowSource {
    /// Edge lengths of the flow metric at time `t_mid`, integrated from frame `frame` (`solver.between`).
    fn between(&self, frame: &MetricFrame, t_mid: f64) -> Vec<f64>;
    /// Edge lengths of u_s = s u(0) + (1 - s) u(T) (`clean_start`'s blend).
    fn blend(&self, s: f64) -> Vec<f64>;
    /// Convex Alexandrov polytope of a convex metric (BIRealizer), if available.
    fn convex_anchor(&self, _lengths: &[f64]) -> Option<Vec<V3>> {
        None
    }
}

pub const COLD: [(f64, usize); 4] = [(1e-1, 6), (1e-2, 5), (1e-3, 4), (1e-4, 4)];
pub const WARM: [(f64, usize); 2] = [(1e-3, 3), (1e-4, 3)];
pub const POLISH: [(f64, usize); 2] = [(1e-5, 6), (1e-6, 6)];

/// `true_flow.FastRealizer`: fading-bending LM for every frame, warm from the previous one.
pub struct FastRealizer {
    pub mesh: RealizeMesh,
    pub x0: Vec<V3>,
    pub seed_kind: SeedKind,
    pub len_tol: f64,
    setup: RefCell<FadeSetup>,
}

impl FastRealizer {
    pub fn new(mesh: RealizeMesh, x0: Vec<V3>, seed_kind: SeedKind) -> Self {
        let setup = FadeSetup::new(&mesh.faces, &mesh.edges, &mesh.cones, mesh.nv);
        FastRealizer { mesh, x0, seed_kind, len_tol: 3e-3, setup: RefCell::new(setup) }
    }

    fn realize_x(&self, l: &[f64], prev: Option<&Realized>) -> Realized {
        let warm = prev.map_or(false, |p| p.ok);
        let x0 = if warm { prev.unwrap().x.clone() } else { self.x0.clone() };
        let mut runs: Vec<(&str, &[(f64, usize)])> = Vec::new();
        if warm {
            runs.push(("warm", &WARM));
        }
        runs.push((if warm { "warm_full" } else { "cold" }, &COLD));
        let mut s = self.setup.borrow_mut();
        let mut x = x0.clone();
        let mut start = String::new();
        let mut iters = 0;
        for (name, st) in runs {
            start = name.to_string();
            let (mut xr, mut info) = lm::fade_lm(&mut s, l, &x0, st);
            iters = info.n_iter;
            if info.len_rel_max > self.len_tol {
                let (x2, i2) = lm::fade_lm(&mut s, l, &xr, &POLISH);
                xr = x2;
                info = i2;
                iters = info.n_iter;
                start += "+polish";
            }
            x = xr;
            if info.len_rel_max <= self.len_tol {
                break;
            }
        }
        drop(s);
        if warm {
            x = geom::kabsch_align(&x, &prev.unwrap().x);
        }
        let mut r = self.diag(x, start, iters);
        if self.seed_kind != SeedKind::Caps {
            r.warning = Some(
                "cone-cap seed unavailable (overlapping caps or beta >= 1): started from the bare sphere, the 3D shape may be on a wrong (dented) branch"
                    .into(),
            );
        }
        r
    }

    fn diag(&self, x: Vec<V3>, start: String, iters: usize) -> Realized {
        diag_of(&self.mesh, x, start, iters, self.len_tol)
    }
}

fn diag_of(mesh: &RealizeMesh, x: Vec<V3>, start: String, iters: usize, len_tol: f64) -> Realized {
    // `len_rel_*` are filled by the caller's lengths: computed lazily in `with_lengths`
    let d = geom::dihedrals(&x, &mesh.faces);
    let n_reflex = d.iter().filter(|&&v| v < -1e-9).count();
    let min_d = d.iter().cloned().fold(f64::INFINITY, f64::min);
    let _ = len_tol;
    Realized { x, ok: false, len_rel_max: f64::NAN, len_rel_rms: f64::NAN, n_reflex, min_dihedral: min_d, lm_iters: iters, start, warning: None }
}

impl FastRealizer {
    /// `FastRealizer.realize(frame, prev)` (diag included; ok = len_rel_max < 5 len_tol).
    pub fn realize(&self, l: &[f64], prev: Option<&Realized>) -> Realized {
        let mut r = self.realize_x(l, prev);
        set_len(&mut r, &self.mesh.edges, l, self.len_tol);
        r
    }
}

fn set_len(r: &mut Realized, e: &[[usize; 2]], l: &[f64], len_tol: f64) {
    let rel = geom::len_rel(&r.x, e, l);
    r.len_rel_max = rel.iter().cloned().fold(0.0, f64::max);
    r.len_rel_rms = (rel.iter().map(|v| v * v).sum::<f64>() / rel.len() as f64).sqrt();
    r.ok = r.len_rel_max < 5.0 * len_tol;
}

/// `true_flow.n_deep_folds`: edges folded inward by more than `thresh` rad (default 0.15).
pub fn n_deep_folds(x: &[V3], faces: &[[usize; 3]], thresh: f64) -> usize {
    geom::dihedrals(x, faces).iter().filter(|&&d| d < -thresh).count()
}

fn deep(rlz: &FastRealizer, r: &Realized) -> usize {
    n_deep_folds(&r.x, &rlz.mesh.faces, 0.15)
}

fn convex_anchor(rlz: &FastRealizer, src: &dyn FlowSource, l: &[f64]) -> Option<Realized> {
    let xb = src.convex_anchor(l)?;
    let bi = Realized { x: xb, ok: true, ..diag_of(&rlz.mesh, vec![], String::new(), 0, rlz.len_tol) };
    Some(rlz.realize(l, Some(&bi)))
}

/// `true_flow.clean_start(solver, realizer)`: a t = 0 shape connected to the uniformized end.
/// `f0` = lengths of u(0), `ft` = lengths of u(T); `src.blend(s)` = lengths of s u(0) + (1-s) u(T).
pub fn clean_start(rlz: &FastRealizer, src: &dyn FlowSource, f0: &[f64], ft: &[f64], nstep: usize) -> Option<Realized> {
    let cold = rlz.realize(f0, None);
    if cold.ok && deep(rlz, &cold) == 0 {
        return Some(cold);
    }
    let mut cur = rlz.realize(ft, None);
    if !cur.ok || deep(rlz, &cur) > 0 {
        if let Some(ca) = convex_anchor(rlz, src, ft) {
            cur = ca;
        }
    }
    for j in 1..nstep {
        let s = j as f64 / nstep as f64;
        let rm = rlz.realize(&src.blend(s), Some(&cur));
        if rm.ok {
            cur = rm;
        }
    }
    let mut r0 = rlz.realize(f0, Some(&cur));
    if !r0.ok || (cold.ok && deep(rlz, &r0) >= deep(rlz, &cold)) {
        return Some(cold);
    }
    r0.start = "clean_start".into();
    Some(r0)
}

/// `_cleaner(a, b)` is b?  (the less folded, then the clearly more accurate; a wins ties)
fn cleaner_picks_b(rlz: &FastRealizer, a: &Realized, b: &Realized) -> bool {
    if !b.ok {
        return false;
    }
    if !a.ok {
        return true;
    }
    let (da, db) = (deep(rlz, a), deep(rlz, b));
    if db != da {
        return db < da;
    }
    b.len_rel_max < 0.5 * a.len_rel_max.max(rlz.len_tol)
}

/// `true_flow.realize_backward`: all frames realized backwards from the uniformized end, with
/// hidden in-between metrics where the lengths jump by more than `step_rel` (0.08, `max_sub` 4),
/// and a lazily computed forward chain from `start` wherever the backward fit folds.
pub fn realize_backward(
    rlz: &FastRealizer,
    src: &dyn FlowSource,
    frs: &[MetricFrame],
    start: Option<Realized>,
    progress: &mut dyn FnMut(usize, usize),
) -> Vec<Realized> {
    let (step_rel, max_sub) = (0.08, 4usize);
    let n = frs.len();
    let mut fwd: Vec<Option<Realized>> = vec![None; n];
    fwd[0] = Some(match start {
        Some(s) => s,
        None => rlz.realize(&frs[0].lengths, None),
    });
    let mut out: Vec<Option<Realized>> = vec![None; n];
    let mut cur = rlz.realize(&frs[n - 1].lengths, None);
    if !cur.ok || deep(rlz, &cur) > 0 {
        if let Some(ca) = convex_anchor(rlz, src, &frs[n - 1].lengths) {
            if cleaner_picks_b(rlz, &cur, &ca) {
                cur = ca;
            }
        }
    }
    out[n - 1] = Some(cur);
    progress(1, n);
    for k in (0..n.saturating_sub(1)).rev() {
        let nf = &frs[k + 1];
        let nx = out[k + 1].clone().unwrap();
        let dl = frs[k].lengths.iter().zip(&nf.lengths).map(|(a, b)| (b / a).ln().abs()).fold(0.0, f64::max);
        let nsub = ((dl / step_rel).ceil() as usize).min(max_sub).max(1) - 1;
        let mut c = nx.clone();
        for j in (1..=nsub).rev() {
            let t = frs[k].t + (nf.t - frs[k].t) * j as f64 / (nsub + 1) as f64;
            let rm = rlz.realize(&src.between(&frs[k], t), Some(&c));
            if rm.ok {
                c = rm;
            }
        }
        let rb = rlz.realize(&frs[k].lengths, Some(&c));
        let mut keep = rb.clone();
        if !rb.ok || deep(rlz, &rb) > 0 {
            // forward(k)
            let mut j = (0..n).filter(|&i| fwd[i].is_some()).max().unwrap();
            while j < k {
                let prev = fwd[j].as_ref().filter(|p| p.ok).cloned();
                fwd[j + 1] = Some(rlz.realize(&frs[j + 1].lengths, prev.as_ref()));
                j += 1;
            }
            let rf = fwd[k].clone().unwrap();
            if cleaner_picks_b(rlz, &rb, &rf) {
                let mut r = rf;
                r.x = geom::kabsch_align(&r.x, &nx.x);
                keep = r;
            }
        }
        out[k] = Some(keep);
        progress(n - k, n);
    }
    out.into_iter().map(|r| r.unwrap()).collect()
}

// ------------------------------------------------------------------ polytope
/// `true_flow.polytope_snap`: every vertex moved radially (from the corners' centroid) onto the
/// convex hull of the corner vertices.  X unchanged if the hull is degenerate.
pub fn polytope_snap(x: &[V3], corner_ids: &[usize]) -> Vec<V3> {
    let p: Vec<V3> = corner_ids.iter().map(|&i| x[i]).collect();
    let Some(h) = hull::convex_hull(&p) else { return x.to_vec() };
    let mut c = [0.0; 3];
    for q in &p {
        c = geom::add(c, *q);
    }
    c = geom::scale(c, 1.0 / p.len() as f64);
    x.iter()
        .map(|xv| {
            let d = geom::sub(*xv, c);
            let mut s = f64::INFINITY;
            for (nn, b) in h.normals.iter().zip(&h.offsets) {
                let nd = geom::dot(d, *nn);
                if nd > 1e-15 {
                    s = s.min((b - geom::dot(*nn, c)) / nd);
                }
            }
            if s.is_finite() { geom::add(c, geom::scale(d, s)) } else { *xv }
        })
        .collect()
}

struct PolyLsq<'a> {
    nv: usize,
    ids: &'a [usize],
    cidx: Vec<Option<usize>>,
    free: Vec<usize>,
    fidx: Vec<usize>, // vertex -> index in free
    tri: Vec<[usize; 3]>, // per free vertex: corner indices (into ids)
    e: &'a [[usize; 2]],
    l: &'a [f64],
}

impl PolyLsq<'_> {
    fn pos(&self, p: &[f64]) -> Vec<V3> {
        let nc = self.ids.len();
        let cc = |i: usize| [p[3 * i], p[3 * i + 1], p[3 * i + 2]];
        let mut out = vec![[0.0; 3]; self.nv];
        for (i, &v) in self.ids.iter().enumerate() {
            out[v] = cc(i);
        }
        for (k, &v) in self.free.iter().enumerate() {
            let t = self.tri[k];
            let (q0, q1) = (p[3 * nc + 2 * k], p[3 * nc + 2 * k + 1]);
            let a = cc(t[0]);
            out[v] = geom::add(a, geom::add(geom::scale(geom::sub(cc(t[1]), a), q0), geom::scale(geom::sub(cc(t[2]), a), q1)));
        }
        out
    }
    /// d P_v / d params as (param index, coefficient on a 3-vector or the vector itself)
    fn vertex_jac(&self, p: &[f64], v: usize, sgn: f64, gvec: V3, out: &mut Vec<(usize, f64)>) {
        let nc = self.ids.len();
        if let Some(i) = self.cidx[v] {
            for c in 0..3 {
                out.push((3 * i + c, sgn * gvec[c]));
            }
        } else {
            let k = self.fidx[v];
            let t = self.tri[k];
            let (q0, q1) = (p[3 * nc + 2 * k], p[3 * nc + 2 * k + 1]);
            let w = [1.0 - q0 - q1, q0, q1];
            for m in 0..3 {
                for c in 0..3 {
                    out.push((3 * t[m] + c, sgn * w[m] * gvec[c]));
                }
            }
            let cc = |i: usize| [p[3 * i], p[3 * i + 1], p[3 * i + 2]];
            let u = geom::sub(cc(t[1]), cc(t[0]));
            let vv = geom::sub(cc(t[2]), cc(t[0]));
            out.push((3 * nc + 2 * k, sgn * geom::dot(u, gvec)));
            out.push((3 * nc + 2 * k + 1, sgn * geom::dot(vv, gvec)));
        }
    }
    fn support(&self, v: usize) -> Vec<usize> {
        let nc = self.ids.len();
        if let Some(i) = self.cidx[v] {
            (3 * i..3 * i + 3).collect()
        } else {
            let k = self.fidx[v];
            // Python's jac_sparsity: all corner columns + the vertex's 2 columns
            let mut s: Vec<usize> = (0..3 * nc).collect();
            s.push(3 * nc + 2 * k);
            s.push(3 * nc + 2 * k + 1);
            s
        }
    }
}

impl LsqProblem for PolyLsq<'_> {
    fn n_params(&self) -> usize {
        3 * self.ids.len() + 2 * self.free.len()
    }
    fn n_resid(&self) -> usize {
        self.e.len()
    }
    fn row_support(&self, i: usize) -> Vec<usize> {
        let mut s = self.support(self.e[i][0]);
        s.extend(self.support(self.e[i][1]));
        s.sort_unstable();
        s.dedup();
        s
    }
    fn residuals(&self, p: &[f64]) -> Vec<f64> {
        let x = self.pos(p);
        self.e.iter().zip(self.l).map(|(e, l)| geom::norm(geom::sub(x[e[0]], x[e[1]])) / l - 1.0).collect()
    }
    fn jac_rows(&self, p: &[f64]) -> Vec<(Vec<usize>, Vec<f64>)> {
        let x = self.pos(p);
        let mut buf = Vec::new();
        self.e
            .iter()
            .zip(self.l)
            .map(|(e, l)| {
                let d = geom::sub(x[e[0]], x[e[1]]);
                let g = geom::scale(d, 1.0 / (geom::norm(d) * l));
                buf.clear();
                self.vertex_jac(p, e[0], 1.0, g, &mut buf);
                self.vertex_jac(p, e[1], -1.0, g, &mut buf);
                buf.sort_by_key(|a| a.0);
                let mut idx: Vec<usize> = Vec::new();
                let mut val: Vec<f64> = Vec::new();
                for &(j, v) in buf.iter() {
                    if idx.last() == Some(&j) {
                        *val.last_mut().unwrap() += v;
                    } else {
                        idx.push(j);
                        val.push(v);
                    }
                }
                (idx, val)
            })
            .collect()
    }
}

/// `true_flow.polytope_fit(X, corner_ids, E, L, rounds=2, max_nfev=30)`: corners free in 3D, every
/// other vertex on a face of their convex hull (2 barycentric coordinates), lengths fit by
/// least squares.  Returns (X, rel_len_err) or (X, None) if the hull fails.
pub fn polytope_fit(x: &[V3], corner_ids: &[usize], e: &[[usize; 2]], l: &[f64]) -> (Vec<V3>, Option<Vec<f64>>) {
    let (rounds, max_nfev) = (2, 30);
    let nv = x.len();
    let ids = corner_ids;
    let mut cidx = vec![None; nv];
    for (i, &v) in ids.iter().enumerate() {
        cidx[v] = Some(i);
    }
    let free: Vec<usize> = (0..nv).filter(|&v| cidx[v].is_none()).collect();
    let mut fidx = vec![usize::MAX; nv];
    for (k, &v) in free.iter().enumerate() {
        fidx[v] = k;
    }
    let mut x = x.to_vec();
    let mut c: Vec<V3> = ids.iter().map(|&i| x[i]).collect();
    let mut rel = None;
    for rnd in 0..rounds {
        let Some(h) = hull::convex_hull(&c) else { return (x, rel) };
        let t = h.simplices;
        let z = if rnd == 0 { polytope_snap(&x, ids) } else { x.clone() };
        let mut tri = Vec::with_capacity(free.len());
        let mut st = Vec::with_capacity(free.len());
        for &v in &free {
            let mut best = (f64::INFINITY, 0usize, [0.0; 2]);
            for (fi, tt) in t.iter().enumerate() {
                let a = c[tt[0]];
                let u = geom::sub(c[tt[1]], a);
                let w = geom::sub(c[tt[2]], a);
                // pinv([u w]) (z - a): 2x2 normal equations
                let r = geom::sub(z[v], a);
                let (uu, uw, ww) = (geom::dot(u, u), geom::dot(u, w), geom::dot(w, w));
                let det = uu * ww - uw * uw;
                let (bu, bw) = (geom::dot(u, r), geom::dot(w, r));
                let q = if det.abs() > 1e-300 { [(ww * bu - uw * bw) / det, (uu * bw - uw * bu) / det] } else { [0.0, 0.0] };
                let pnt = geom::add(a, geom::add(geom::scale(u, q[0]), geom::scale(w, q[1])));
                let pen = (-q[0]).max(0.0) + (-q[1]).max(0.0) + (q[0] + q[1] - 1.0).max(0.0);
                let d = geom::norm(geom::sub(pnt, z[v])) + pen;
                if d < best.0 {
                    best = (d, fi, q);
                }
            }
            tri.push(t[best.1]);
            st.push(best.2);
        }
        let prob = PolyLsq { nv, ids, cidx: cidx.clone(), free: free.clone(), fidx: fidx.clone(), tri, e, l };
        let mut p0: Vec<f64> = c.iter().flat_map(|q| q.iter().copied()).collect();
        for q in &st {
            p0.extend_from_slice(q);
        }
        let res = lm::least_squares_jac_scaled(&prob, &p0, max_nfev);
        x = prob.pos(&res.x);
        c = ids.iter().map(|&i| x[i]).collect();
        rel = Some(prob.residuals(&res.x).into_iter().map(f64::abs).collect());
    }
    (x, rel)
}

/// `flow_true.hull_faces`: re-triangulate points on a convex polytope by their (joggled) hull,
/// outward oriented; None if a vertex would be dropped.
pub fn hull_faces(x: &[V3]) -> Option<Vec<[usize; 3]>> {
    let n = x.len();
    let mut c = [0.0; 3];
    for p in x {
        c = geom::add(c, *p);
    }
    c = geom::scale(c, 1.0 / n as f64);
    let r2: Vec<f64> = x.iter().map(|p| geom::dot(geom::sub(*p, c), geom::sub(*p, c))).collect();
    let rmax = r2.iter().cloned().fold(0.0, f64::max);
    let xb: Vec<V3> = x.iter().zip(&r2).map(|(p, r)| geom::add(c, geom::scale(geom::sub(*p, c), 1.0 + 1e-4 * (1.0 - r / rmax)))).collect();
    let h = hull::convex_hull(&xb)?;
    let mut used = vec![false; n];
    for t in &h.simplices {
        for &v in t {
            used[v] = true;
        }
    }
    if used.iter().any(|u| !u) {
        return None;
    }
    Some(
        h.simplices
            .into_iter()
            .map(|t| {
                let nrm = geom::cross(geom::sub(x[t[1]], x[t[0]]), geom::sub(x[t[2]], x[t[0]]));
                let m = geom::scale(geom::add(x[t[0]], geom::add(x[t[1]], x[t[2]])), 1.0 / 3.0);
                if geom::dot(nrm, geom::sub(m, c)) < 0.0 { [t[0], t[2], t[1]] } else { t }
            })
            .collect(),
    )
}

/// Result of [`polytope_vertices`].
#[derive(Clone, Debug)]
pub struct PolytopeFrame {
    pub x: Vec<V3>,
    pub faces: Vec<[usize; 3]>,
    pub len_rel_max: f64,
    pub len_rel_rms: f64,
    pub n_reflex: usize,
    pub min_dihedral: f64,
}

/// `flow_true.polytope_vertices`: the honest polytope at t = 0 (Polytope seed) on the same mesh,
/// aligned to the smooth fit `x`, drawn on its hull triangulation.  None = keep the fit.
pub fn polytope_vertices(mesh: &RealizeMesh, x: &[V3], corner_ids: &[usize], l: &[f64]) -> Option<PolytopeFrame> {
    let (xp, rel) = polytope_fit(x, corner_ids, &mesh.edges, l);
    let rel = rel?;
    if !xp.iter().all(|p| p.iter().all(|c| c.is_finite())) {
        return None;
    }
    let xp = geom::kabsch_align(&xp, x);
    let fp = hull_faces(&xp).unwrap_or_else(|| mesh.faces.clone());
    let dh = geom::dihedrals(&xp, &fp);
    Some(PolytopeFrame {
        len_rel_max: rel.iter().cloned().fold(0.0, f64::max),
        len_rel_rms: (rel.iter().map(|v| v * v).sum::<f64>() / rel.len() as f64).sqrt(),
        n_reflex: dh.iter().filter(|&&d| d < -1e-9).count(),
        min_dihedral: dh.iter().cloned().fold(f64::INFINITY, f64::min),
        x: xp,
        faces: fp,
    })
}
