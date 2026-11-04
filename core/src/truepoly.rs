//! "True Polytope" seed: an actual convex polytope whose corners carry the cone angles, meshed so
//! that its crease edges ARE mesh edges (each flat face triangulated inside).  The t = 0 metric
//! is exactly the polytope's flat-faced metric and its 3D realization is the polytope itself.
//!
//! Polytope: vertices r_i p_i at the cone directions p_i (stereographic image of the cone
//! positions) plus 4 farthest-point "extra corners" at radius 1.  The cone radii r_i are solved
//! (damped Newton, finite-difference Jacobian on the convex hull) so that the curvature of each
//! cone vertex is exactly 2 pi (1 - beta_i); the extra corners take up the rest
//! (4 pi - sum, Gauss-Bonnet).  The extra corners are ordinary (smooth) flow vertices: their
//! curvature is not preserved, it flows to rho A like every other non-cone vertex.
//!
//! Mesh: crease edges subdivided at spacing h; each (merged coplanar) facet gets a hexagonal
//! lattice of interior points and is triangulated by 2D Delaunay (lower hull of the lifted
//! points).  The solver's unit-sphere positions are the radial projection of the polytope
//! (origin inside), so the cone vertices sit exactly at the cone points.

use crate::hull::convex_hull;
use crate::mesh::Mesh;
use crate::realize::geom::{self, V3};
use crate::seedmesh::{self, FlowMesh};
use crate::sparse::{Coo, SparseLu};
use std::collections::HashMap;
use std::f64::consts::PI;

#[derive(Clone, Debug)]
pub struct TruePoly {
    pub mesh: Mesh,
    pub flow_mesh: FlowMesh,
    /// 3D polytope position of every mesh vertex
    pub x: Vec<V3>,
    /// extra-corner mesh vertex ids (not cones)
    pub corners: Vec<usize>,
    /// mesh edge indices (into mesh.e) lying on polytope crease edges
    pub crease_edges: Vec<usize>,
    /// polytope edge lengths (mesh.e order) = the flow's l0
    pub l0: Vec<f64>,
    /// least-squares vertex-scaling log factor of the polytope metric relative to the round
    /// sphere chords (display offset for the conformal-factor raster)
    pub u_off: Vec<f64>,
    /// solved cone radii and the max cone curvature error of the polytope
    pub radii: Vec<f64>,
    pub defect_err: f64,
}

fn hull_defects(pts: &[V3]) -> Option<(Vec<f64>, crate::hull::Hull)> {
    let h = convex_hull(pts)?;
    let mut ang = vec![0.0; pts.len()];
    let mut on = vec![false; pts.len()];
    for t in &h.simplices {
        for k in 0..3 {
            let (a, b, c) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
            let u = geom::sub(pts[b], pts[a]);
            let v = geom::sub(pts[c], pts[a]);
            ang[a] += geom::norm(geom::cross(u, v)).atan2(geom::dot(u, v));
            on[a] = true;
        }
    }
    let d = (0..pts.len()).map(|i| if on[i] { 2.0 * PI - ang[i] } else { 0.0 }).collect();
    Some((d, h))
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

/// 4 farthest-point corners on the sphere, away from the cones.
fn extra_corners(p: &[V3], n: usize) -> Vec<V3> {
    let cand = fibonacci(4000);
    let mut d: Vec<f64> = cand.iter().map(|c| p.iter().map(|q| geom::norm(geom::sub(*c, *q))).fold(f64::INFINITY, f64::min)).collect();
    let mut out = vec![];
    for _ in 0..n {
        let j = (0..cand.len()).max_by(|&a, &b| d[a].partial_cmp(&d[b]).unwrap()).unwrap();
        out.push(cand[j]);
        for k in 0..cand.len() {
            d[k] = d[k].min(geom::norm(geom::sub(cand[k], cand[j])));
        }
    }
    out
}

fn solve_lin(a: &mut Vec<Vec<f64>>, b: &mut Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for c in 0..n {
        let p = (c..n).max_by(|&i, &j| a[i][c].abs().partial_cmp(&a[j][c].abs()).unwrap())?;
        if a[p][c].abs() < 1e-14 {
            return None;
        }
        a.swap(c, p);
        b.swap(c, p);
        for r in c + 1..n {
            let f = a[r][c] / a[c][c];
            for k in c..n {
                a[r][k] -= f * a[c][k];
            }
            b[r] -= f * b[c];
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        let s: f64 = (r + 1..n).map(|k| a[r][k] * x[k]).sum();
        x[r] = (b[r] - s) / a[r][r];
    }
    Some(x)
}

/// Solve the cone radii; returns (polytope vertices: cones then corners, radii, max error).
pub fn solve_polytope(p: &[V3], betas: &[f64], corners: &[V3]) -> Result<(Vec<V3>, Vec<f64>, f64), String> {
    let n = p.len();
    let tgt: Vec<f64> = betas.iter().map(|b| 2.0 * PI * (1.0 - b)).collect();
    let pts_of = |r: &[f64]| -> Vec<V3> {
        let mut v: Vec<V3> = (0..n).map(|i| geom::scale(p[i], r[i])).collect();
        v.extend_from_slice(corners);
        v
    };
    let res = |r: &[f64]| -> Option<Vec<f64>> {
        let (d, _) = hull_defects(&pts_of(r))?;
        Some((0..n).map(|i| d[i] - tgt[i]).collect())
    };
    let nrm = |v: &[f64]| v.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    // start: tangent-cone apex height 1/beta, capped
    let mut r: Vec<f64> = betas.iter().map(|b| (1.0 / b).min(4.0)).collect();
    let mut f = res(&r).ok_or("polytope hull failed")?;
    for _ in 0..100 {
        let fn0 = nrm(&f);
        if fn0 < 1e-13 {
            break;
        }
        let mut jac = vec![vec![0.0; n]; n];
        for j in 0..n {
            let hstep = 1e-7 * r[j];
            let mut rp = r.clone();
            rp[j] += hstep;
            let mut rm = r.clone();
            rm[j] -= hstep;
            let (fp, fm) = (res(&rp).ok_or("hull failed")?, res(&rm).ok_or("hull failed")?);
            for i in 0..n {
                jac[i][j] = (fp[i] - fm[i]) / (2.0 * hstep);
            }
        }
        let mut b: Vec<f64> = f.iter().map(|x| -x).collect();
        let dr = solve_lin(&mut jac, &mut b).ok_or("true polytope: singular Jacobian (a cone left the hull)")?;
        let mut s = 1.0;
        let mut acc = false;
        while s > 1e-6 {
            let rn: Vec<f64> = r.iter().zip(&dr).map(|(a, d)| (a + s * d).clamp(0.2, 60.0)).collect();
            if let Some(fnew) = res(&rn) {
                if nrm(&fnew) < (1.0 - 1e-4 * s) * fn0 {
                    r = rn;
                    f = fnew;
                    acc = true;
                    break;
                }
            }
            s *= 0.5;
        }
        if !acc {
            break;
        }
    }
    let err = nrm(&f);
    if err > 1e-9 {
        return Err(format!("true polytope: cone curvatures not reached (error {err:.1e})"));
    }
    let pts = pts_of(&r);
    let (d, _) = hull_defects(&pts).unwrap();
    if d.iter().any(|&x| x < 1e-6) {
        return Err("true polytope: an extra corner was swallowed by the hull".into());
    }
    Ok((pts, r, err))
}

/// 2D Delaunay of points (lower hull of the lift); triangles CCW.
fn delaunay2(q: &[[f64; 2]]) -> Option<Vec<[usize; 3]>> {
    let lift: Vec<V3> = q.iter().map(|a| [a[0], a[1], a[0] * a[0] + a[1] * a[1]]).collect();
    let h = convex_hull(&lift)?;
    let mut out = vec![];
    for (t, nn) in h.simplices.iter().zip(&h.normals) {
        if nn[2] < -1e-12 {
            // outward normal pointing down = lower hull; seen from above it is CW -> flip to CCW
            out.push([t[0], t[2], t[1]]);
        }
    }
    Some(out)
}

fn hash01(i: usize, s: u64) -> f64 {
    let mut x = (i as u64).wrapping_mul(0x9E3779B97F4A7C15) ^ s;
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51afd7ed558ccd);
    x ^= x >> 33;
    (x >> 11) as f64 / (1u64 << 53) as f64 - 0.5
}

pub fn build(positions: &[[f64; 2]], betas: &[f64], level: (usize, usize, usize)) -> Result<TruePoly, String> {
    let n = positions.len();
    let p = seedmesh::stereographic_to_sphere(positions);
    let cor = extra_corners(&p, 4);
    let (pv, radii, defect_err) = solve_polytope(&p, betas, &cor)?;
    let nvp = pv.len();
    let hull = convex_hull(&pv).ok_or("polytope hull failed")?;
    // merge coplanar hull triangles into facets
    let nt = hull.simplices.len();
    let mut gid = vec![usize::MAX; nt];
    let mut groups: Vec<Vec<usize>> = vec![];
    for t in 0..nt {
        if gid[t] != usize::MAX {
            continue;
        }
        let g = groups.len();
        let mut mem = vec![];
        for s in t..nt {
            if gid[s] == usize::MAX
                && geom::dot(hull.normals[t], hull.normals[s]) > 1.0 - 1e-12
                && (hull.offsets[t] - hull.offsets[s]).abs() < 1e-10
            {
                gid[s] = g;
                mem.push(s);
            }
        }
        groups.push(mem);
    }
    // facet boundary (directed, CCW seen from outside)
    let mut fb: Vec<Vec<(usize, usize)>> = vec![];
    for mem in &groups {
        let mut cnt: HashMap<(usize, usize), usize> = HashMap::new();
        for &t in mem {
            let s = hull.simplices[t];
            for k in 0..3 {
                cnt.insert((s[k], s[(k + 1) % 3]), 1);
            }
        }
        let mut b = vec![];
        for (&(a, c), _) in &cnt {
            if !cnt.contains_key(&(c, a)) {
                b.push((a, c));
            }
        }
        // order into a loop
        let mut nxt: HashMap<usize, usize> = b.iter().cloned().collect();
        let start = b.iter().map(|e| e.0).min().unwrap();
        let mut lp = vec![];
        let mut a = start;
        loop {
            let c = nxt.remove(&a).ok_or("facet boundary not a loop")?;
            lp.push((a, c));
            a = c;
            if a == start {
                break;
            }
        }
        if !nxt.is_empty() {
            return Err("facet boundary not a single loop".into());
        }
        fb.push(lp);
    }
    // spacing h from the page level's vertex count
    let area: f64 = hull
        .simplices
        .iter()
        .map(|t| 0.5 * geom::norm(geom::cross(geom::sub(pv[t[1]], pv[t[0]]), geom::sub(pv[t[2]], pv[t[0]]))))
        .sum();
    // same vertex budget as the sphere solver mesh of this page level
    let target_nv = match seedmesh::build_mesh(positions, betas, level) {
        Ok((m, _)) => m.nv() as f64,
        Err(_) => (level.0 + n * level.1 * level.2.max(1)) as f64,
    };
    let h = (2.0 * area / (3f64.sqrt() * target_nv)).sqrt();
    // vertices: polytope vertices first, then crease points, then facet interiors
    let mut x: Vec<V3> = pv.clone();
    let mut edge_pts: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for lp in &fb {
        for &(a, c) in lp {
            let key = (a.min(c), a.max(c));
            if edge_pts.contains_key(&key) {
                continue;
            }
            let (a0, c0) = key;
            let len = geom::norm(geom::sub(pv[c0], pv[a0]));
            let k = ((len / h).round() as usize).max(1);
            let mut ids = vec![a0];
            for j in 1..k {
                let s = j as f64 / k as f64;
                ids.push(x.len());
                x.push(geom::add(pv[a0], geom::scale(geom::sub(pv[c0], pv[a0]), s)));
            }
            ids.push(c0);
            edge_pts.insert(key, ids);
        }
    }
    let mut faces: Vec<[usize; 3]> = vec![];
    for (g, lp) in fb.iter().enumerate() {
        let nn = hull.normals[groups[g][0]];
        let o = pv[lp[0].0];
        let e1 = {
            let d = geom::sub(pv[lp[0].1], o);
            geom::scale(d, 1.0 / geom::norm(d))
        };
        let e2 = geom::cross(nn, e1);
        let to2 = |q: V3| -> [f64; 2] {
            let d = geom::sub(q, o);
            [geom::dot(d, e1) / h, geom::dot(d, e2) / h]
        };
        // boundary ids in loop order
        let mut ids: Vec<usize> = vec![];
        for &(a, c) in lp {
            let key = (a.min(c), a.max(c));
            let mut ch = edge_pts[&key].clone();
            if ch[0] != a {
                ch.reverse();
            }
            ids.extend_from_slice(&ch[..ch.len() - 1]);
        }
        let mut q: Vec<[f64; 2]> = ids.iter().map(|&i| to2(x[i])).collect();
        let poly2: Vec<[f64; 2]> = lp.iter().map(|&(a, _)| to2(pv[a])).collect();
        // inward half-planes (CCW polygon: inside is to the left)
        let lines: Vec<([f64; 2], [f64; 2])> = (0..poly2.len())
            .map(|k| {
                let (a, b) = (poly2[k], poly2[(k + 1) % poly2.len()]);
                let d = [b[0] - a[0], b[1] - a[1]];
                let l = (d[0] * d[0] + d[1] * d[1]).sqrt();
                (a, [-d[1] / l, d[0] / l])
            })
            .collect();
        let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
        for a in &poly2 {
            for c in 0..2 {
                lo[c] = lo[c].min(a[c]);
                hi[c] = hi[c].max(a[c]);
            }
        }
        let nb = ids.len();
        let rows = ((hi[1] - lo[1]) / (0.5 * 3f64.sqrt())).ceil() as i64 + 1;
        let cols = (hi[0] - lo[0]).ceil() as i64 + 2;
        for jr in 0..rows {
            for jc in -1..cols {
                let yy = lo[1] + jr as f64 * 0.5 * 3f64.sqrt();
                let xx = lo[0] + jc as f64 + if jr % 2 == 1 { 0.5 } else { 0.0 };
                let dmin = lines.iter().map(|(a, nrm)| (xx - a[0]) * nrm[0] + (yy - a[1]) * nrm[1]).fold(f64::INFINITY, f64::min);
                if dmin > 0.55 {
                    let k = q.len();
                    q.push([xx + 1e-6 * hash01(k, g as u64), yy + 1e-6 * hash01(k, 7 + g as u64)]);
                }
            }
        }
        let mut gidx = ids.clone();
        for k in nb..q.len() {
            let [u, v] = q[k];
            gidx.push(x.len());
            x.push(geom::add(o, geom::add(geom::scale(e1, u * h), geom::scale(e2, v * h))));
        }
        let tri = delaunay2(&q).ok_or("facet triangulation failed")?;
        // check coverage
        let a_poly: f64 = (0..poly2.len())
            .map(|k| {
                let (a, b) = (poly2[k], poly2[(k + 1) % poly2.len()]);
                a[0] * b[1] - a[1] * b[0]
            })
            .sum::<f64>()
            * 0.5;
        let a_tri: f64 = tri
            .iter()
            .map(|t| {
                let (a, b, c) = (q[t[0]], q[t[1]], q[t[2]]);
                0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
            })
            .sum();
        if (a_tri - a_poly).abs() > 1e-6 * a_poly.abs().max(1.0) {
            return Err(format!("facet triangulation does not cover the facet ({a_tri} vs {a_poly})"));
        }
        for t in tri {
            faces.push([gidx[t[0]], gidx[t[1]], gidx[t[2]]]);
        }
    }
    // unit-sphere positions: radial projection (origin is inside the polytope)
    let y: Vec<V3> = x.iter().map(|q| geom::scale(*q, 1.0 / geom::norm(*q))).collect();
    let cones: Vec<usize> = (0..n).collect();
    let corners: Vec<usize> = (n..nvp).collect();
    let mesh = Mesh::from_sphere(y.clone(), faces.clone(), cones.clone());
    // sanity: closed 2-manifold sphere
    let ne = mesh.e.len();
    if mesh.nv() as i64 - ne as i64 + mesh.f.len() as i64 != 2 {
        return Err("true polytope mesh is not a sphere".into());
    }
    let l0: Vec<f64> = mesh.e.iter().map(|&[i, j]| geom::norm(geom::sub(x[i], x[j]))).collect();
    let crease: std::collections::HashSet<(usize, usize)> = edge_pts
        .values()
        .flat_map(|ids| ids.windows(2).map(|w| (w[0].min(w[1]), w[0].max(w[1]))).collect::<Vec<_>>())
        .collect();
    let crease_edges: Vec<usize> = (0..ne).filter(|&k| crease.contains(&(mesh.e[k][0], mesh.e[k][1]))).collect();
    // display offset: least squares (u_i + u_j)/2 = log(l0 / chord)
    let chord = mesh.chord_lengths();
    let nv = mesh.nv();
    let mut a = Coo::new(nv);
    let mut rhs = vec![0.0; nv];
    for (k, &[i, j]) in mesh.e.iter().enumerate() {
        let g = (l0[k] / chord[k]).ln();
        for (r, c) in [(i, i), (i, j), (j, i), (j, j)] {
            a.push(r, c, 0.25);
        }
        rhs[i] += 0.5 * g;
        rhs[j] += 0.5 * g;
    }
    for v in 0..nv {
        a.push(v, v, 1e-12);
    }
    let u_off = SparseLu::factor(&a)?.solve(&rhs);
    let seed = seedmesh::Seed {
        vertices: pv.clone(),
        faces: hull.simplices.clone(),
        tip_ids: cones.clone(),
        ring_ids: vec![],
        seed_betas: betas.to_vec(),
        seed_alphas: betas.iter().map(|b| b.acos()).collect(),
        positions_xyz: p,
        overlaps: vec![],
        exclusion_factor: 0.0,
    };
    let flow_mesh = FlowMesh { vertices: x.clone(), faces: mesh.f.clone(), sphere_xyz: mesh.y.clone(), tip_ids: cones, seed };
    Ok(TruePoly { mesh, flow_mesh, x, corners, crease_edges, l0, u_off, radii, defect_err })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::defects_and_areas;
    #[test]
    fn page_default_polytope() {
        let pos = [[0.7, 0.0], [-0.7, 0.0], [0.0, 0.7], [0.0, -0.7]];
        let betas = [0.6, 0.7, 0.5, 0.8];
        let tp = build(&pos, &betas, seedmesh::page_level("coarse").unwrap()).unwrap();
        let g = defects_and_areas(&tp.mesh.f, &tp.mesh.opp, &tp.l0, tp.mesh.nv());
        for (i, b) in betas.iter().enumerate() {
            assert!((g.k[i] - 2.0 * PI * (1.0 - b)).abs() < 1e-10);
        }
        let tot: f64 = g.k.iter().sum();
        assert!((tot - 4.0 * PI).abs() < 1e-10);
        for v in 4 + tp.corners.len()..tp.mesh.nv() {
            assert!(g.k[v].abs() < 1e-10, "flat vertex {v} has defect {}", g.k[v]);
        }
        // creases are exactly the folded edges
        let dh = geom::dihedrals(&tp.x, &tp.mesh.f);
        let te = geom::face_edge_table(&tp.mesh.f);
        let mut folded = vec![];
        for (k, ([a, b], _)) in te.iter().enumerate() {
            if dh[k].abs() > 1e-9 {
                folded.push((*a.min(b), *a.max(b)));
            }
            assert!(dh[k] > -1e-9);
        }
        folded.sort();
        let mut cr: Vec<(usize, usize)> = tp.crease_edges.iter().map(|&k| (tp.mesh.e[k][0], tp.mesh.e[k][1])).collect();
        cr.sort();
        assert_eq!(folded, cr);
    }
}
