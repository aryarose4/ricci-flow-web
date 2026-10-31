//! /flow mesh construction from cone positions + betas: port of
//! `viz/seed_surface.build_seed` (Fibonacci points + tangent cone caps, incremental hull),
//! `viz/cone_mesh.build_flow_mesh` (concentric rings inside each cap) and
//! `embeddings/true_flow.build_mesh` (unit sphere, outward faces), plus the page mesh levels of
//! `viewer/backend/flow_stream.PAGE_MESH` / `_saddle_level`.
//!
//! The hull is the seed builder's own incremental hull (`seed_surface._hull`), not `hull.rs`
//! (Qhull stand-in): the face ORDER of `_hull` depends on CPython's set iteration order over the
//! horizon edges, which is emulated here ([`PySet`]) so the faces come out in the identical order.

use crate::mesh::Mesh;
use std::f64::consts::PI;

type V3 = [f64; 3];

pub const BETA_MIN: f64 = 0.05;
pub const BETA_MAX: f64 = 0.99;
pub const GAP: f64 = 0.08;
const EXCL_FACTORS: [f64; 4] = [1.0, 1.6, 2.4, 3.6];

/// Page mesh levels (N Fibonacci points, m spokes, K rings) — `flow_stream.PAGE_MESH`.
pub fn page_level(name: &str) -> Option<(usize, usize, usize)> {
    match name.trim().to_ascii_lowercase().as_str() {
        "smoke" => Some((150, 12, 2)),
        "coarse" => Some((400, 16, 4)),
        "medium" => Some((1200, 24, 8)),
        "fine" => Some((2500, 32, 12)),
        "extra" => Some((4000, 36, 14)),
        _ => None,
    }
}

#[inline]
fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
fn cross(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
#[inline]
fn norm(a: V3) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// `m9.stereographic_to_sphere` (plane (x, y) -> unit sphere, infinity = north pole).
pub fn stereographic_to_sphere(p: &[[f64; 2]]) -> Vec<V3> {
    p.iter()
        .map(|q| {
            let r2 = q[0] * q[0] + q[1] * q[1];
            let den = 1.0 + r2;
            [2.0 * q[0] / den, 2.0 * q[1] / den, (r2 - 1.0) / den]
        })
        .collect()
}

fn angle(p: V3, q: V3) -> f64 {
    norm(cross(p, q)).atan2(dot(p, q))
}

/// A cap overlap (pair, angle between the cone points, sum of target cap radii).
#[derive(Clone, Debug)]
pub struct Overlap {
    pub i: usize,
    pub j: usize,
    pub angle: f64,
    pub alpha_sum: f64,
}

pub fn cap_overlaps(p: &[V3], alphas: &[f64]) -> Vec<Overlap> {
    let mut out = vec![];
    for i in 0..p.len() {
        for j in i + 1..p.len() {
            let d = angle(p[i], p[j]);
            if d <= alphas[i] + alphas[j] {
                out.push(Overlap { i, j, angle: d, alpha_sum: alphas[i] + alphas[j] });
            }
        }
    }
    out
}

pub fn clamp_alphas(p: &[V3], alphas: &[f64], gap: f64) -> Vec<f64> {
    let mut a = alphas.to_vec();
    for _ in 0..200 {
        let mut changed = false;
        for i in 0..p.len() {
            for j in i + 1..p.len() {
                let d = angle(p[i], p[j]);
                if a[i] + a[j] > d - gap + 1e-12 {
                    let s = (d - gap).max(1e-3) / (a[i] + a[j]);
                    a[i] *= s;
                    a[j] *= s;
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }
    a.iter().map(|&x| x.max(1e-3)).collect()
}

fn fibonacci(n: usize) -> Vec<V3> {
    let g = 3.0 - 5f64.sqrt();
    (0..n)
        .map(|i| {
            let k = i as f64 + 0.5;
            let z = 1.0 - 2.0 * k / n as f64;
            let r = (1.0 - z * z).max(0.0).sqrt();
            let ph = k * PI * g;
            [r * ph.cos(), r * ph.sin(), z]
        })
        .collect()
}

fn tip_height(alpha: f64, beta_s: f64, m: usize) -> f64 {
    let r = alpha.sin();
    let l = r * (PI / m as f64).sin() / (PI * beta_s / m as f64).sin();
    alpha.cos() + (l * l - r * r).max(0.0).sqrt()
}

fn basis(p: V3) -> (V3, V3) {
    let mut ax = 0;
    for c in 1..3 {
        if p[c].abs() < p[ax].abs() {
            ax = c;
        }
    }
    let mut a = [0.0; 3];
    a[ax] = 1.0;
    let e1 = cross(p, a);
    let n = norm(e1);
    let e1 = [e1[0] / n, e1[1] / n, e1[2] / n];
    (e1, cross(p, e1))
}

fn points(p: &[V3], a_u: &[f64], m: usize, nres: usize, factor: f64) -> (Vec<V3>, usize) {
    let n = p.len();
    let mut pts: Vec<V3> = (0..n)
        .map(|i| {
            let h = tip_height(a_u[i], a_u[i].cos(), m);
            [p[i][0] * h, p[i][1] * h, p[i][2] * h]
        })
        .collect();
    let s0 = 0.9 * (4.0 * PI / nres as f64).sqrt();
    let mut n_sph = 0;
    for s in fibonacci(nres) {
        if (0..n).all(|i| norm(cross(s, p[i])).atan2(dot(s, p[i])) > a_u[i] + factor * s0) {
            pts.push(s);
            n_sph += 1;
        }
    }
    for i in 0..n {
        let (e1, e2) = basis(p[i]);
        let (ca, sa) = (a_u[i].cos(), a_u[i].sin());
        for k in 0..m {
            let t = 2.0 * PI * k as f64 / m as f64 + 0.3;
            let (ct, st) = (t.cos(), t.sin());
            let mut q = [0.0; 3];
            for c in 0..3 {
                q[c] = ca * p[i][c] + sa * (ct * e1[c] + st * e2[c]);
            }
            pts.push(q);
        }
    }
    (pts, n_sph)
}

// ------------------------------------------------------------------ CPython set emulation
const XXPRIME_1: u64 = 11400714785074694791;
const XXPRIME_2: u64 = 14029467366897019727;
const XXPRIME_5: u64 = 2870177450012600261;

/// `hash((a, b))` for small non-negative ints (CPython >= 3.8, 64-bit).
fn tuple_hash(a: usize, b: usize) -> u64 {
    let mut acc = XXPRIME_5;
    for lane in [a as u64, b as u64] {
        acc = acc.wrapping_add(lane.wrapping_mul(XXPRIME_2));
        acc = acc.rotate_left(31);
        acc = acc.wrapping_mul(XXPRIME_1);
    }
    acc = acc.wrapping_add(2u64 ^ (XXPRIME_5 ^ 3527539));
    if acc == u64::MAX {
        return 1546275796;
    }
    acc
}

/// Insertion-only CPython `set` of int pairs; iterates in CPython's table order.
struct PySet {
    table: Vec<Option<(u64, (usize, usize))>>,
    fill: usize,
}

impl PySet {
    fn new() -> Self {
        PySet { table: vec![None; 8], fill: 0 }
    }
    fn probe(table: &[Option<(u64, (usize, usize))>], h: u64, key: Option<(usize, usize)>) -> (usize, bool) {
        let mask = table.len() - 1;
        let mut i = (h as usize) & mask;
        let mut perturb = h;
        loop {
            let probes = if i + 9 <= mask { 9 } else { 0 };
            for j in 0..=probes {
                match table[i + j] {
                    None => return (i + j, false),
                    Some((hh, k)) => {
                        if Some(k) == key && hh == h {
                            return (i + j, true);
                        }
                    }
                }
            }
            perturb >>= 5;
            i = (i.wrapping_mul(5).wrapping_add(1).wrapping_add(perturb as usize)) & mask;
        }
    }
    fn contains(&self, k: (usize, usize)) -> bool {
        Self::probe(&self.table, tuple_hash(k.0, k.1), Some(k)).1
    }
    fn add(&mut self, k: (usize, usize)) {
        let h = tuple_hash(k.0, k.1);
        let (slot, found) = Self::probe(&self.table, h, Some(k));
        if found {
            return;
        }
        self.table[slot] = Some((h, k));
        self.fill += 1;
        let mask = self.table.len() - 1;
        if self.fill * 5 < mask * 3 {
            return;
        }
        let used = self.fill;
        let minused = if used > 50000 { used * 2 } else { used * 4 };
        let mut newsize = 8;
        while newsize <= minused {
            newsize <<= 1;
        }
        let old = std::mem::replace(&mut self.table, vec![None; newsize]);
        for e in old.into_iter().flatten() {
            let (s, _) = Self::probe(&self.table, e.0, None);
            self.table[s] = Some(e);
        }
    }
    fn iter(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.table.iter().flatten().map(|e| e.1)
    }
}

/// `seed_surface._hull`: incremental hull in point order, outward triangles, identical face order.
pub fn seed_hull(pts: &[V3]) -> Vec<[usize; 3]> {
    let n = pts.len();
    let i0 = 0;
    let i1 = (0..n).find(|&k| norm(sub(pts[k], pts[i0])) > 1e-6).unwrap();
    let i2 = (0..n).find(|&k| norm(cross(sub(pts[k], pts[i0]), sub(pts[i1], pts[i0]))) > 1e-6).unwrap();
    let n0 = cross(sub(pts[i1], pts[i0]), sub(pts[i2], pts[i0]));
    let i3 = (0..n).find(|&k| dot(sub(pts[k], pts[i0]), n0).abs() > 1e-6).unwrap();
    let mut cen = [0.0; 3];
    for c in 0..3 {
        cen[c] = (pts[i0][c] + pts[i1][c] + pts[i2][c] + pts[i3][c]) / 4.0;
    }
    let mut f: Vec<[usize; 3]> = vec![];
    for [a, b, c] in [[i0, i1, i2], [i0, i1, i3], [i0, i2, i3], [i1, i2, i3]] {
        let nrm = cross(sub(pts[b], pts[a]), sub(pts[c], pts[a]));
        f.push(if dot(nrm, sub(pts[a], cen)) < 0.0 { [a, c, b] } else { [a, b, c] });
    }
    let used = [i0, i1, i2, i3];
    for k in 0..n {
        if used.contains(&k) {
            continue;
        }
        let vis: Vec<bool> = f
            .iter()
            .map(|t| {
                let a = pts[t[0]];
                let nn = cross(sub(pts[t[1]], a), sub(pts[t[2]], a));
                let l = norm(nn);
                let nn = [nn[0] / l, nn[1] / l, nn[2] / l];
                dot(nn, sub(pts[k], a)) > 1e-11
            })
            .collect();
        if !vis.iter().any(|&v| v) {
            continue;
        }
        let mut edges = PySet::new();
        for (t, _) in f.iter().zip(&vis).filter(|(_, &v)| v) {
            for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                edges.add((a, b));
            }
        }
        let mut keep: Vec<[usize; 3]> = f.iter().zip(&vis).filter(|(_, &v)| !v).map(|(t, _)| *t).collect();
        for (a, b) in edges.iter() {
            if !edges.contains((b, a)) {
                keep.push([a, b, k]);
            }
        }
        f = keep;
    }
    f
}

fn tips_valid(f: &[[usize; 3]], n: usize, n_sph: usize, m: usize) -> bool {
    for i in 0..n {
        let rows: Vec<&[usize; 3]> = f.iter().filter(|t| t.contains(&i)).collect();
        if rows.len() != m {
            return false;
        }
        let mut nb: Vec<usize> = rows.iter().flat_map(|t| t.iter().cloned()).filter(|&x| x != i).collect();
        nb.sort_unstable();
        nb.dedup();
        let base = n + n_sph + i * m;
        if nb != (base..base + m).collect::<Vec<_>>() {
            return false;
        }
    }
    true
}

/// `seed_surface.build_seed` output.
#[derive(Clone, Debug)]
pub struct Seed {
    pub vertices: Vec<V3>,
    pub faces: Vec<[usize; 3]>,
    pub tip_ids: Vec<usize>,
    pub ring_ids: Vec<Vec<usize>>,
    pub seed_betas: Vec<f64>,
    pub seed_alphas: Vec<f64>,
    pub positions_xyz: Vec<V3>,
    pub overlaps: Vec<Overlap>,
    pub exclusion_factor: f64,
}

pub fn build_seed(positions: &[[f64; 2]], betas: &[f64], resolution: usize, m: usize) -> Result<Seed, String> {
    if positions.len() != betas.len() {
        return Err("positions and betas differ in length".into());
    }
    if betas.iter().any(|&b| !(b >= BETA_MIN - 1e-12 && b <= BETA_MAX + 1e-12)) {
        return Err(format!("beta outside [{BETA_MIN}, {BETA_MAX}]: beta >= 1 has no convex realization"));
    }
    let p = stereographic_to_sphere(positions);
    let n = p.len();
    let a_t: Vec<f64> = betas.iter().map(|b| b.acos()).collect();
    let over = cap_overlaps(&p, &a_t);
    let mut a_u = if over.is_empty() { a_t.clone() } else { clamp_alphas(&p, &a_t, GAP) };
    let mut found = None;
    'outer: for _ in 0..8 {
        for &factor in &EXCL_FACTORS {
            let (pts, n_sph) = points(&p, &a_u, m, resolution, factor);
            let f = seed_hull(&pts);
            if tips_valid(&f, n, n_sph, m) {
                found = Some((pts, n_sph, f, factor));
                break 'outer;
            }
        }
        a_u.iter_mut().for_each(|a| *a *= 0.9);
    }
    let (pts, n_sph, f, factor) = found.ok_or("could not build a valid seed")?;
    let mut usedv = vec![false; pts.len()];
    for t in &f {
        for &v in t {
            usedv[v] = true;
        }
    }
    let mut remap = vec![usize::MAX; pts.len()];
    let mut vertices = vec![];
    for (i, &u) in usedv.iter().enumerate() {
        if u {
            remap[i] = vertices.len();
            vertices.push(pts[i]);
        }
    }
    let faces: Vec<[usize; 3]> = f.iter().map(|t| [remap[t[0]], remap[t[1]], remap[t[2]]]).collect();
    let ring_ids = (0..n).map(|i| (0..m).map(|k| remap[n + n_sph + i * m + k]).collect()).collect();
    let tip_ids = (0..n).map(|i| remap[i]).collect();
    Ok(Seed {
        vertices,
        faces,
        tip_ids,
        ring_ids,
        seed_betas: a_u.iter().map(|a| a.cos()).collect(),
        seed_alphas: a_u,
        positions_xyz: p,
        overlaps: over,
        exclusion_factor: factor,
    })
}

/// `cone_mesh.ring_fractions` ("power" grading).
pub fn ring_fractions(k: usize, beta: f64) -> Vec<f64> {
    let b = beta.clamp(1e-3, 1.0);
    (1..=k).map(|j| (j as f64 / (k + 1) as f64).powf(1.0 / b)).collect()
}

/// `cone_mesh.build_flow_mesh` output (subset the solver and the page use).
#[derive(Clone, Debug)]
pub struct FlowMesh {
    /// preview polyhedron (tips at the closed-form apex)
    pub vertices: Vec<V3>,
    pub faces: Vec<[usize; 3]>,
    /// radial projection to the unit sphere, tips exactly at the cone points
    pub sphere_xyz: Vec<V3>,
    pub tip_ids: Vec<usize>,
    pub seed: Seed,
}

pub fn build_flow_mesh(positions: &[[f64; 2]], betas: &[f64], nres: usize, m: usize, kk: usize) -> Result<FlowMesh, String> {
    if m < 3 || nres < 20 {
        return Err("need m >= 3, K >= 0, N >= 20".into());
    }
    let seed = build_seed(positions, betas, nres, m)?;
    let mut v = seed.vertices.clone();
    let mut f = seed.faces.clone();
    let tips = seed.tip_ids.clone();
    if kk > 0 {
        let mut nf: Vec<[usize; 3]> = f.iter().filter(|t| !t.iter().any(|x| tips.contains(x))).cloned().collect();
        let mut nv = v.len();
        for (i, (&t_id, ring)) in tips.iter().zip(&seed.ring_ids).enumerate() {
            let tf = ring_fractions(kk, seed.seed_betas[i]);
            let tp = seed.vertices[t_id];
            let mut ids = vec![vec![0usize; m]; kk + 1];
            for jj in 0..kk {
                for (k, &r) in ring.iter().enumerate() {
                    let q = seed.vertices[r];
                    v.push([tp[0] + tf[jj] * (q[0] - tp[0]), tp[1] + tf[jj] * (q[1] - tp[1]), tp[2] + tf[jj] * (q[2] - tp[2])]);
                    ids[jj][k] = nv + k;
                }
                nv += m;
            }
            ids[kk] = ring.clone();
            for k in 0..m {
                let k1 = (k + 1) % m;
                nf.push([t_id, ids[0][k], ids[0][k1]]);
                for jj in 0..kk {
                    let (a, b, c, d) = (ids[jj][k], ids[jj][k1], ids[jj + 1][k], ids[jj + 1][k1]);
                    if (k + jj) % 2 == 0 {
                        nf.push([a, c, d]);
                        nf.push([a, d, b]);
                    } else {
                        nf.push([a, c, b]);
                        nf.push([b, c, d]);
                    }
                }
            }
        }
        f = nf;
    }
    // _orient_out
    for t in f.iter_mut() {
        let (a, b, c) = (v[t[0]], v[t[1]], v[t[2]]);
        let nn = cross(sub(b, a), sub(c, a));
        let s = [a[0] + b[0] + c[0], a[1] + b[1] + c[1], a[2] + b[2] + c[2]];
        if dot(nn, s) < 0.0 {
            t.swap(1, 2);
        }
    }
    let mut s: Vec<V3> = v
        .iter()
        .map(|q| {
            let l = norm(*q);
            [q[0] / l, q[1] / l, q[2] / l]
        })
        .collect();
    for (i, &t) in tips.iter().enumerate() {
        s[t] = seed.positions_xyz[i];
    }
    Ok(FlowMesh { vertices: v, faces: f, sphere_xyz: s, tip_ids: tips, seed })
}

/// `true_flow.build_mesh`: the solver mesh (unit sphere, outward faces, edges, cone ids).
pub fn build_mesh(positions: &[[f64; 2]], betas: &[f64], level: (usize, usize, usize)) -> Result<(Mesh, FlowMesh), String> {
    let fm = build_flow_mesh(positions, betas, level.0, level.1, level.2)?;
    let mesh = Mesh::from_sphere(fm.sphere_xyz.clone(), fm.faces.clone(), fm.tip_ids.clone());
    Ok((mesh, fm))
}
