//! 3D convex hull (incremental, triangulated facets, outward orientation): the parts of
//! `scipy.spatial.ConvexHull` used by `true_flow.polytope_snap` / `polytope_fit` and
//! `flow_true.hull_faces` (`simplices` and `equations`).
//!
//! Coplanar facets come out triangulated (like Qhull's `simplices`), but the triangulation of a
//! flat facet may differ from Qhull's.  Points within `eps` of the current hull are not added.

use crate::realize::geom::{self, V3};
use std::collections::HashMap;

/// Hull triangles (outward, counter-clockwise seen from outside) and plane equations
/// `(n, b)` with unit outward normal n, `n . x <= b` inside (Qhull: `equations = [n, -b]`).
#[derive(Clone, Debug)]
pub struct Hull {
    pub simplices: Vec<[usize; 3]>,
    pub normals: Vec<V3>,
    pub offsets: Vec<f64>,
}

fn plane(p: &[V3], t: [usize; 3]) -> (V3, f64) {
    let n = geom::cross(geom::sub(p[t[1]], p[t[0]]), geom::sub(p[t[2]], p[t[0]]));
    let l = geom::norm(n);
    let n = if l > 0.0 { geom::scale(n, 1.0 / l) } else { n };
    (n, geom::dot(n, p[t[0]]))
}

/// Convex hull of `p`; None if the points are (nearly) coplanar or fewer than 4.
pub fn convex_hull(p: &[V3]) -> Option<Hull> {
    let n = p.len();
    if n < 4 {
        return None;
    }
    let mut lo = p[0];
    let mut hi = p[0];
    for q in p {
        for c in 0..3 {
            lo[c] = lo[c].min(q[c]);
            hi[c] = hi[c].max(q[c]);
        }
    }
    let diam = geom::norm(geom::sub(hi, lo));
    if diam == 0.0 {
        return None;
    }
    let eps = 1e-13 * diam;
    // initial tetrahedron: extreme points
    let i0 = (0..n).min_by(|&a, &b| p[a][0].partial_cmp(&p[b][0]).unwrap())?;
    let i1 = (0..n).max_by(|&a, &b| {
        geom::norm(geom::sub(p[a], p[i0])).partial_cmp(&geom::norm(geom::sub(p[b], p[i0]))).unwrap()
    })?;
    let e01 = geom::sub(p[i1], p[i0]);
    let i2 = (0..n).max_by(|&a, &b| {
        let fa = geom::norm(geom::cross(e01, geom::sub(p[a], p[i0])));
        let fb = geom::norm(geom::cross(e01, geom::sub(p[b], p[i0])));
        fa.partial_cmp(&fb).unwrap()
    })?;
    let nrm = geom::cross(e01, geom::sub(p[i2], p[i0]));
    if geom::norm(nrm) <= eps * diam {
        return None;
    }
    let i3 = (0..n).max_by(|&a, &b| {
        geom::dot(nrm, geom::sub(p[a], p[i0])).abs().partial_cmp(&geom::dot(nrm, geom::sub(p[b], p[i0])).abs()).unwrap()
    })?;
    let vol = geom::dot(nrm, geom::sub(p[i3], p[i0]));
    if vol.abs() <= eps * geom::norm(nrm) {
        return None;
    }
    let mut faces: Vec<Option<[usize; 3]>> = Vec::new();
    let mut planes: Vec<(V3, f64)> = Vec::new();
    for (t, opp) in [([i0, i1, i2], i3), ([i0, i1, i3], i2), ([i0, i2, i3], i1), ([i1, i2, i3], i0)] {
        let (nn, b) = plane(p, t);
        let t = if geom::dot(nn, p[opp]) - b > 0.0 { [t[0], t[2], t[1]] } else { t };
        faces.push(Some(t));
        planes.push(plane(p, t));
    }
    // insertion order: farthest-first is unnecessary for correctness at this size; insert in index order
    let used = [i0, i1, i2, i3];
    for v in 0..n {
        if used.contains(&v) {
            continue;
        }
        let vis: Vec<usize> = (0..faces.len())
            .filter(|&f| faces[f].is_some() && geom::dot(planes[f].0, p[v]) - planes[f].1 > eps)
            .collect();
        if vis.is_empty() {
            continue;
        }
        let mut dir: HashMap<(usize, usize), usize> = HashMap::new();
        for &f in &vis {
            let t = faces[f].unwrap();
            for k in 0..3 {
                dir.insert((t[k], t[(k + 1) % 3]), f);
            }
        }
        let mut horizon = Vec::new();
        for &f in &vis {
            let t = faces[f].unwrap();
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                if !dir.contains_key(&(b, a)) {
                    horizon.push((a, b));
                }
            }
        }
        for &f in &vis {
            faces[f] = None;
        }
        for (a, b) in horizon {
            let t = [a, b, v];
            faces.push(Some(t));
            planes.push(plane(p, t));
        }
    }
    let mut h = Hull { simplices: vec![], normals: vec![], offsets: vec![] };
    for (f, t) in faces.iter().enumerate() {
        if let Some(t) = t {
            h.simplices.push(*t);
            h.normals.push(planes[f].0);
            h.offsets.push(planes[f].1);
        }
    }
    Some(h)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cube_hull() {
        let mut p = vec![];
        for i in 0..8 {
            p.push([(i & 1) as f64, ((i >> 1) & 1) as f64, ((i >> 2) & 1) as f64]);
        }
        p.push([0.5, 0.5, 0.5]);
        let h = convex_hull(&p).unwrap();
        assert_eq!(h.simplices.len(), 12);
        for (nn, b) in h.normals.iter().zip(&h.offsets) {
            for q in &p {
                assert!(geom::dot(*nn, *q) <= b + 1e-12);
            }
        }
    }
}
