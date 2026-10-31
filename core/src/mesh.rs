//! Sphere mesh combinatorics + intrinsic (edge-length) geometry.
//! Port of `embeddings/discrete_conformal.py` (face_geometry, defects_and_areas,
//! _cot_laplacian, _area_jacobian, _edge_index), `thick_part.edges_of` and the
//! orientation part of `true_flow.build_mesh`.
use crate::sparse::Coo;
use std::collections::HashMap;
use std::f64::consts::PI;

/// Triangulated unit sphere with cone vertices (fixed combinatorics).
#[derive(Clone, Debug)]
pub struct Mesh {
    /// unit-sphere vertex positions
    pub y: Vec<[f64; 3]>,
    /// faces, outward oriented
    pub f: Vec<[usize; 3]>,
    /// sorted unique edges (i < j), lexicographic order (= np.unique(axis=0))
    pub e: Vec<[usize; 2]>,
    /// for each face corner k, index of the edge opposite to it
    pub opp: Vec<[usize; 3]>,
    /// cone vertex ids (same order as betas)
    pub cones: Vec<usize>,
}

/// `thick_part.edges_of`: sorted unique undirected edges.
pub fn edges_of(f: &[[usize; 3]]) -> Vec<[usize; 2]> {
    let mut e: Vec<[usize; 2]> = Vec::with_capacity(3 * f.len());
    for t in f {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            e.push([a.min(b), a.max(b)]);
        }
    }
    e.sort();
    e.dedup();
    e
}

/// `discrete_conformal._edge_index`.
pub fn edge_index(f: &[[usize; 3]], e: &[[usize; 2]]) -> Vec<[usize; 3]> {
    let key: HashMap<(usize, usize), usize> = e.iter().enumerate().map(|(k, x)| ((x[0], x[1]), k)).collect();
    f.iter()
        .map(|&[a, b, c]| {
            let mut o = [0usize; 3];
            for (k, (x, y)) in [(b, c), (c, a), (a, b)].into_iter().enumerate() {
                o[k] = key[&(x.min(y), x.max(y))];
            }
            o
        })
        .collect()
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// `embed4._signed_volume`.
pub fn signed_volume(x: &[[f64; 3]], f: &[[usize; 3]]) -> f64 {
    f.iter().map(|t| dot(x[t[0]], cross(x[t[1]], x[t[2]]))).sum::<f64>() / 6.0
}

impl Mesh {
    /// Build from sphere points + faces (as `build_mesh` after the hull):
    /// normalise to the unit sphere, flip all faces if the signed volume is
    /// negative, compute edges and the opposite-edge table.
    pub fn from_sphere(y: Vec<[f64; 3]>, mut f: Vec<[usize; 3]>, cones: Vec<usize>) -> Self {
        let y: Vec<[f64; 3]> = y
            .into_iter()
            .map(|p| {
                let n = dot(p, p).sqrt();
                [p[0] / n, p[1] / n, p[2] / n]
            })
            .collect();
        if signed_volume(&y, &f) < 0.0 {
            for t in f.iter_mut() {
                t.swap(1, 2);
            }
        }
        Self::from_parts(y, f, cones)
    }
    /// Use the given arrays verbatim (no normalisation or reorientation).
    pub fn from_parts(y: Vec<[f64; 3]>, f: Vec<[usize; 3]>, cones: Vec<usize>) -> Self {
        let e = edges_of(&f);
        let opp = edge_index(&f, &e);
        Mesh { y, f, e, opp, cones }
    }
    pub fn nv(&self) -> usize {
        self.y.len()
    }
    /// Round-sphere chord lengths L0.
    pub fn chord_lengths(&self) -> Vec<f64> {
        self.e
            .iter()
            .map(|&[i, j]| {
                let d = [self.y[i][0] - self.y[j][0], self.y[i][1] - self.y[j][1], self.y[i][2] - self.y[j][2]];
                dot(d, d).sqrt()
            })
            .collect()
    }
}

/// Per-face geometry: corner cotangents, face areas, validity.
#[derive(Clone, Debug)]
pub struct FaceGeom {
    pub ang: Vec<[f64; 3]>,
    pub cot: Vec<[f64; 3]>,
    pub af: Vec<f64>,
    pub valid: bool,
}

/// `discrete_conformal.face_geometry`.
pub fn face_geometry(opp: &[[usize; 3]], l: &[f64]) -> FaceGeom {
    let nf = opp.len();
    let mut g = FaceGeom { ang: vec![[0.0; 3]; nf], cot: vec![[0.0; 3]; nf], af: vec![0.0; nf], valid: true };
    for t in 0..nf {
        let lf = [l[opp[t][0]], l[opp[t][1]], l[opp[t][2]]];
        let (a, b, c) = (lf[0], lf[1], lf[2]);
        let s = 0.5 * (a + b + c);
        let q = s * (s - a) * (s - b) * (s - c);
        if !(q > 0.0) {
            g.valid = false;
        }
        let af = q.max(0.0).sqrt();
        g.af[t] = af;
        for k in 0..3 {
            let (x, y, z) = (lf[k], lf[(k + 1) % 3], lf[(k + 2) % 3]);
            let num = y * y + z * z - x * x;
            let cs = (num / (2.0 * y * z)).clamp(-1.0, 1.0);
            g.ang[t][k] = cs.acos();
            g.cot[t][k] = num / (4.0 * af.max(1e-300));
        }
    }
    g
}

/// Angle defects K_v, barycentric vertex areas A_v (face area / 3), face geometry.
#[derive(Clone, Debug)]
pub struct VertexGeom {
    pub k: Vec<f64>,
    pub av: Vec<f64>,
    pub face: FaceGeom,
}

/// `discrete_conformal.defects_and_areas`.
pub fn defects_and_areas(f: &[[usize; 3]], opp: &[[usize; 3]], l: &[f64], nv: usize) -> VertexGeom {
    let face = face_geometry(opp, l);
    let mut angsum = vec![0.0; nv];
    let mut asum = vec![0.0; nv];
    for (t, tri) in f.iter().enumerate() {
        for k in 0..3 {
            angsum[tri[k]] += face.ang[t][k];
            asum[tri[k]] += face.af[t];
        }
    }
    let k = angsum.iter().map(|s| 2.0 * PI - s).collect();
    let av = asum.iter().map(|s| s / 3.0).collect();
    VertexGeom { k, av, face }
}

/// `discrete_conformal._cot_laplacian` (= dK/du), as COO.
pub fn cot_laplacian(f: &[[usize; 3]], cot: &[[f64; 3]], out: &mut Coo) {
    for (t, tri) in f.iter().enumerate() {
        for k in 0..3 {
            let (i, j) = (tri[(k + 1) % 3], tri[(k + 2) % 3]);
            let w = 0.5 * cot[t][k];
            out.push(i, j, -w);
            out.push(j, i, -w);
            out.push(i, i, w);
            out.push(j, j, w);
        }
    }
}

/// `discrete_conformal._area_jacobian`: d A_v / d u_w, as COO (scaled).
pub fn area_jacobian(f: &[[usize; 3]], opp: &[[usize; 3]], l: &[f64], cot: &[[f64; 3]], out: &mut Coo) {
    for (t, tri) in f.iter().enumerate() {
        let mut lf2c = [0.0; 3];
        for k in 0..3 {
            let le = l[opp[t][k]];
            lf2c[k] = le * le * cot[t][k] / 4.0;
        }
        for m in 0..3 {
            let g = (lf2c[(m + 1) % 3] + lf2c[(m + 2) % 3]) / 3.0;
            for a in 0..3 {
                out.push(tri[a], tri[m], g);
            }
        }
    }
}
