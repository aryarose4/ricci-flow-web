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
//! points).  The solver's unit-sphere positions `mesh.y` are discrete-conformal coordinates of
//! the polytope metric (see `conformalize`), and the polytope's cone directions are iterated so
//! the cones' conformal positions equal the user's cone points (see `build_opts`).

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
    /// vertex-scaling log factor of the polytope metric relative to the round-sphere chords of
    /// its conformal coordinates `mesh.y` (display offset for the conformal-factor raster)
    pub u_off: Vec<f64>,
    /// solved cone radii and the max cone curvature error of the polytope
    pub radii: Vec<f64>,
    pub defect_err: f64,
    /// rms per-edge residual of the `u_off` fit (log(l0/chord) - (u_i+u_j)/2)
    pub fit_resid: f64,
    /// max chordal distance between the cones' conformal positions and the user's positions
    pub match_err: f64,
    pub match_iters: usize,
    /// set when the conformal coordinates failed and radial positions are used instead
    pub warning: Option<String>,
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

/// Angle defect at `apex` of the fan over the closed `link` polygon.
fn apex_defect(apex: V3, link: &[V3]) -> f64 {
    let m = link.len();
    let mut ang = 0.0;
    for k in 0..m {
        let u = geom::sub(link[k], apex);
        let v = geom::sub(link[(k + 1) % m], apex);
        ang += geom::norm(geom::cross(u, v)).atan2(geom::dot(u, v));
    }
    2.0 * PI - ang
}

/// rho < 0 "true polyhedron": a STAR-SHAPED (radial-graph) closed polyhedron, hence embedded.
/// Directions = cone directions + 4 farthest-point corners, plus the midpoint of every cone-cone
/// edge of their spherical Delaunay triangulation (convex hull of the unit vectors), repeated until
/// no two cones are adjacent.  Corners sit at radius 1; with every cone link made of corners only,
/// each cone radius is an independent 1D equation (curvature of a pyramid apex over a fixed link,
/// -> 2 pi as r -> inf), solved by bracketing + bisection.  The corners take the rest of the
/// curvature (4 pi - sum < 0: saddles).  Returns (vertices: cones then corners, faces CCW from
/// outside, cone radii, max cone curvature error).
pub fn solve_star(p: &[V3], betas: &[f64]) -> Result<(Vec<V3>, Vec<[usize; 3]>, Vec<f64>, f64), String> {
    let n = p.len();
    let mut dirs: Vec<V3> = p.to_vec();
    dirs.extend(extra_corners(p, 4));
    let mut tri;
    let mut rounds = 0;
    loop {
        tri = convex_hull(&dirs).ok_or("true polyhedron: direction hull failed")?.simplices;
        let mut cc: Vec<(usize, usize)> = vec![];
        for t in &tri {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                if a < n && b < n {
                    cc.push((a.min(b), a.max(b)));
                }
            }
        }
        if cc.is_empty() {
            break;
        }
        rounds += 1;
        if rounds > 8 {
            return Err("true polyhedron: cones stay adjacent".into());
        }
        cc.sort();
        cc.dedup();
        for (a, b) in cc {
            let m = geom::add(dirs[a], dirs[b]);
            let l = geom::norm(m);
            if l < 1e-9 {
                return Err("true polyhedron: antipodal cones".into());
            }
            dirs.push(geom::scale(m, 1.0 / l));
        }
    }
    let mut radii = vec![1.0; n];
    let mut err = 0.0f64;
    for i in 0..n {
        // link of cone i, ordered around p_i
        let mut nxt: HashMap<usize, usize> = HashMap::new();
        for t in &tri {
            for k in 0..3 {
                if t[k] == i {
                    nxt.insert(t[(k + 1) % 3], t[(k + 2) % 3]);
                }
            }
        }
        let s0 = *nxt.keys().min().ok_or("true polyhedron: isolated cone")?;
        let mut link = vec![s0];
        let mut a = nxt[&s0];
        while a != s0 {
            link.push(a);
            a = *nxt.get(&a).ok_or("true polyhedron: open link")?;
            if link.len() > nxt.len() {
                return Err("true polyhedron: link not a cycle".into());
            }
        }
        let lp: Vec<V3> = link.iter().map(|&j| dirs[j]).collect();
        let tgt = 2.0 * PI * (1.0 - betas[i]);
        let f = |r: f64| apex_defect(geom::scale(p[i], r), &lp) - tgt;
        // bracket: grow r from 1 until the defect exceeds the target, or shrink it (shallow cone)
        let (mut lo, mut hi);
        if f(1.0) < 0.0 {
            lo = 1.0;
            hi = 1.25;
            while f(hi) < 0.0 {
                lo = hi;
                hi *= 1.25;
                if hi > 1e5 {
                    return Err(format!("true polyhedron: cone {i} curvature not reached"));
                }
            }
        } else {
            hi = 1.0;
            lo = 0.9;
            while f(lo) > 0.0 {
                hi = lo;
                lo *= 0.9;
                if lo < 0.3 {
                    return Err(format!(
                        "true polyhedron: cone {} (beta {:.2}) is too flat for its star-shaped link (smallest reachable curvature {:.3} > {:.3})",
                        i + 1,
                        betas[i],
                        f(1.0) + tgt,
                        tgt
                    ));
                }
            }
        }
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if mid <= lo || mid >= hi {
                break;
            }
            if f(mid) < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let r = if f(lo).abs() < f(hi).abs() { lo } else { hi };
        radii[i] = r;
        err = err.max(f(r).abs());
    }
    if err > 1e-9 {
        return Err(format!("true polyhedron: cone curvatures not reached (error {err:.1e})"));
    }
    let pv: Vec<V3> = (0..dirs.len()).map(|j| if j < n { geom::scale(dirs[j], radii[j]) } else { dirs[j] }).collect();
    for t in &tri {
        let nn = geom::cross(geom::sub(pv[t[1]], pv[t[0]]), geom::sub(pv[t[2]], pv[t[0]]));
        if geom::dot(nn, pv[t[0]]) <= 1e-12 * geom::norm(nn) {
            return Err("true polyhedron: face not visible from the centre".into());
        }
    }
    Ok((pv, tri, radii, err))
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

/// Build the True Polytope seed.  The polytope's cone directions are iterated (fixed point on
/// the sphere, Möbius-normalised) so that the CONFORMAL positions of its cones -- not their
/// radial projections -- equal the user's positions; `mesh.y` holds those conformal coordinates.
pub fn build(positions: &[[f64; 2]], betas: &[f64], level: (usize, usize, usize)) -> Result<TruePoly, String> {
    build_opts(positions, betas, level, MATCH_ITERS, true)
}

/// Cheap preview build (slider / drag): radial sphere positions, no conformal map, no match loop.
pub fn build_preview(positions: &[[f64; 2]], betas: &[f64], level: (usize, usize, usize)) -> Result<TruePoly, String> {
    build_opts(positions, betas, level, 0, false)
}

/// max fixed-point iterations of the cone-position match (0 = conformal coordinates only)
pub const MATCH_ITERS: usize = 8;
const MATCH_TOL: f64 = 1e-6;
/// fixed-point relaxation (contracts at ~0.7 per step; 2.0 was measured to diverge)
const MATCH_RELAX: f64 = 1.0;

/// `conformal = false`: radial positions + least-squares `u_off` (the pre-conformal behaviour).
/// If the conformal map fails or is degenerate, falls back to that with `warning` set.
pub fn build_opts(
    positions: &[[f64; 2]],
    betas: &[f64],
    level: (usize, usize, usize),
    iters: usize,
    conformal: bool,
) -> Result<TruePoly, String> {
    let n = positions.len();
    let target = seedmesh::stereographic_to_sphere(positions);
    // same vertex budget as the sphere solver mesh of this page level
    let target_nv = match seedmesh::build_mesh(positions, betas, level) {
        Ok((m, _)) => m.nv() as f64,
        Err(_) => (level.0 + n * level.1 * level.2.max(1)) as f64,
    };
    let mut p = target.clone();
    if !conformal {
        return build_raw(&p, betas, target_nv);
    }
    let mut best: Option<(f64, TruePoly)> = None;
    for it in 0..=iters {
        let mut tp = build_raw(&p, betas, target_nv)?;
        if let Err(e) = conformalize(&mut tp, &target) {
            match best {
                // keep the best conformal iterate found so far
                Some(_) => break,
                None => {
                    // radial p = target here, so the radial build is the pre-conformal seed
                    tp.warning = Some(format!("Polytope: conformal coordinates unavailable ({e}); using radial positions (cone conformal positions approximate)"));
                    return Ok(tp);
                }
            }
        }
        let q: Vec<V3> = tp.mesh.cones.iter().map(|&c| tp.mesh.y[c]).collect();
        let err = q.iter().zip(&target).map(|(a, b)| geom::norm(geom::sub(*a, *b))).fold(0.0, f64::max);
        tp.match_err = err;
        tp.match_iters = it;
        // stop at the tolerance, the iteration cap, or once the error stalls at the mesh
        // re-generation noise floor (the facet meshes change with p)
        let stalled = best.as_ref().map_or(false, |(e, _)| err > 0.85 * *e && err < 1e-2);
        let done = err < MATCH_TOL || it == iters || stalled;
        if best.as_ref().map_or(true, |(e, _)| err < *e) {
            best = Some((err, tp));
        }
        if done {
            break;
        }
        // fixed point: move each cone direction by the remaining mismatch
        for i in 0..n {
            let np = geom::add(p[i], geom::scale(geom::sub(target[i], q[i]), MATCH_RELAX));
            p[i] = geom::scale(np, 1.0 / geom::norm(np));
        }
    }
    Ok(best.unwrap().1)
}

/// Replace the radial sphere positions by discrete-conformal coordinates of the polytope metric:
/// uniformize the polytope metric to constant curvature with NO cones (every vertex smooth),
/// lay the resulting near-round metric out on the unit sphere (Gauss-Newton on chord lengths),
/// Möbius-normalise so the cones best match `target`, and refit the display offset `u_off`.
fn conformalize(tp: &mut TruePoly, target: &[V3]) -> Result<(), String> {
    use crate::flow::{DiscreteConeFlow, FlowOptions};
    let mut m0 = tp.mesh.clone();
    m0.cones = vec![];
    let mut flow = DiscreteConeFlow::new(m0, vec![], FlowOptions::default());
    flow.l0 = tp.l0.clone();
    let (u, _) = crate::fastunif::uniformize(&mut flow).map_err(|e| format!("polytope conformal map: {e}"))?;
    let l = flow.lengths(&u);
    let y = sphere_layout(&tp.mesh.y, &tp.mesh.e, &l)?;
    let cones = tp.mesh.cones.clone();
    let src: Vec<V3> = cones.iter().map(|&c| y[c]).collect();
    let mob = fit_mobius(&src, target);
    let det = cadd(cmul(mob[0], mob[3]), cmul([-mob[1][0], -mob[1][1]], mob[2]));
    if !(det[0].hypot(det[1]) > 1e-8) || mob.iter().any(|c| !c[0].is_finite() || !c[1].is_finite()) {
        return Err("degenerate Möbius normalisation".into());
    }
    let y: Vec<V3> = y.iter().map(|q| mobius_apply(&mob, *q)).collect();
    if y.iter().any(|q| !q.iter().all(|c| c.is_finite())) {
        return Err("non-finite sphere layout".into());
    }
    // every face must stay positively oriented (outward) on the sphere
    let flipped = tp.mesh.f.iter().filter(|t| {
        let (a, b, c) = (y[t[0]], y[t[1]], y[t[2]]);
        geom::dot(geom::add(geom::add(a, b), c), geom::cross(geom::sub(b, a), geom::sub(c, a))) <= 0.0
    }).count();
    if flipped > 0 {
        return Err(format!("sphere layout flipped {flipped} faces"));
    }
    let mesh = Mesh::from_parts(y, tp.mesh.f.clone(), cones);
    let u_off = fit_offset(&mesh, &tp.l0)?;
    if u_off.iter().any(|x| !x.is_finite()) {
        return Err("non-finite conformal offset".into());
    }
    tp.u_off = u_off;
    let ch = mesh.chord_lengths();
    let res: Vec<f64> =
        mesh.e.iter().enumerate().map(|(k, &[i, j])| (tp.l0[k] / ch[k]).ln() - 0.5 * (tp.u_off[i] + tp.u_off[j])).collect();
    tp.fit_resid = (res.iter().map(|x| x * x).sum::<f64>() / res.len() as f64).sqrt();
    tp.flow_mesh.sphere_xyz = mesh.y.clone();
    tp.mesh = mesh;
    Ok(())
}

/// least squares (u_i + u_j)/2 = log(l0 / chord)
fn fit_offset(mesh: &Mesh, l0: &[f64]) -> Result<Vec<f64>, String> {
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
    Ok(SparseLu::factor(&a)?.solve(&rhs))
}

/// Unit-sphere points whose chords best match `l` (Gauss-Newton, sphere constraint as a stiff
/// penalty, renormalised each step), started from `y0`.
fn sphere_layout(y0: &[V3], e: &[[usize; 2]], l: &[f64]) -> Result<Vec<V3>, String> {
    let nv = y0.len();
    let mut y = y0.to_vec();
    const W: f64 = 10.0;
    let cost = |y: &[V3]| -> f64 {
        let mut ss = 0.0;
        for (k, &[i, j]) in e.iter().enumerate() {
            let r = geom::norm(geom::sub(y[i], y[j])) - l[k];
            ss += r * r;
        }
        ss
    };
    // Levenberg-Marquardt damping (mu = 0 is the plain Gauss-Newton of the rho > 0 build): a far
    // start (radial projection of long rho < 0 spikes) needs damped steps to avoid folding
    let mut mu = 0.0;
    let mut last = cost(&y);
    for _ in 0..200 {
        let mut a = Coo::new(3 * nv);
        let mut g = vec![0.0; 3 * nv];
        for (k, &[i, j]) in e.iter().enumerate() {
            let d = geom::sub(y[i], y[j]);
            let dn = geom::norm(d);
            let r = dn - l[k];
            let jv = geom::scale(d, 1.0 / dn);
            for a1 in 0..3 {
                for b1 in 0..3 {
                    let v = jv[a1] * jv[b1];
                    a.push(3 * i + a1, 3 * i + b1, v);
                    a.push(3 * j + a1, 3 * j + b1, v);
                    a.push(3 * i + a1, 3 * j + b1, -v);
                    a.push(3 * j + a1, 3 * i + b1, -v);
                }
                g[3 * i + a1] += jv[a1] * r;
                g[3 * j + a1] -= jv[a1] * r;
            }
        }
        for v in 0..nv {
            let r = W * (geom::dot(y[v], y[v]) - 1.0);
            let jv = geom::scale(y[v], 2.0 * W);
            for a1 in 0..3 {
                for b1 in 0..3 {
                    a.push(3 * v + a1, 3 * v + b1, jv[a1] * jv[b1]);
                }
                a.push(3 * v + a1, 3 * v + a1, 1e-9);
                g[3 * v + a1] += jv[a1] * r;
            }
        }
        let lu_base = a;
        let mut accepted = false;
        let mut step = 0.0;
        for _ in 0..30 {
            let mut am = lu_base.clone();
            if mu > 0.0 {
                for k in 0..3 * nv {
                    am.push(k, k, mu);
                }
            }
            let dx = SparseLu::factor(&am)?.solve(&g);
            let yn: Vec<V3> = (0..nv)
                .map(|v| {
                    let q = [y[v][0] - dx[3 * v], y[v][1] - dx[3 * v + 1], y[v][2] - dx[3 * v + 2]];
                    geom::scale(q, 1.0 / geom::norm(q))
                })
                .collect();
            let c = cost(&yn);
            if c.is_finite() && c <= last * (1.0 + 1e-12) {
                step = dx.iter().fold(0.0f64, |m, x| m.max(x.abs()));
                let dec = last - c;
                y = yn;
                mu *= 0.3;
                if mu < 1e-10 {
                    mu = 0.0;
                }
                accepted = true;
                if step < 1e-11 || dec.abs() < 1e-14 * last.max(1e-30) {
                    return Ok(y);
                }
                last = c;
                break;
            }
            mu = if mu == 0.0 { 1e-4 } else { mu * 10.0 };
        }
        if !accepted || step < 1e-11 {
            break;
        }
    }
    Ok(y)
}

type C = [f64; 2];
fn cmul(a: C, b: C) -> C {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}
fn cadd(a: C, b: C) -> C {
    [a[0] + b[0], a[1] + b[1]]
}
/// Möbius map w -> (a w + b)/(c w + d) acting on the unit sphere (stereographic convention of
/// `seedmesh::stereographic_to_sphere`), via homogeneous coordinates (no special case at infinity).
pub fn mobius_apply(m: &[C; 4], y: V3) -> V3 {
    let (z1, z2) = if y[2] < 0.0 { ([y[0], y[1]], [1.0 - y[2], 0.0]) } else { ([1.0 + y[2], 0.0], [y[0], -y[1]]) };
    let p = cadd(cmul(m[0], z1), cmul(m[1], z2));
    let q = cadd(cmul(m[2], z1), cmul(m[3], z2));
    let (pp, qq) = (p[0] * p[0] + p[1] * p[1], q[0] * q[0] + q[1] * q[1]);
    let s = pp + qq;
    let pq = cmul(p, [q[0], -q[1]]);
    [2.0 * pq[0] / s, 2.0 * pq[1] / s, (pp - qq) / s]
}

/// Least-squares Möbius map (chordal residuals) taking `src` to `dst`; d = 1 (near identity).
pub fn fit_mobius(src: &[V3], dst: &[V3]) -> [C; 4] {
    let mk = |x: &[f64; 6]| -> [C; 4] { [[x[0], x[1]], [x[2], x[3]], [x[4], x[5]], [1.0, 0.0]] };
    let res = |x: &[f64; 6]| -> Vec<f64> {
        let m = mk(x);
        src.iter()
            .zip(dst)
            .flat_map(|(s, d)| {
                let q = mobius_apply(&m, *s);
                [q[0] - d[0], q[1] - d[1], q[2] - d[2]]
            })
            .collect()
    };
    let mut x = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    if src.len() < 2 {
        return mk(&x);
    }
    let mut r = res(&x);
    let mut f = r.iter().map(|v| v * v).sum::<f64>();
    let mut mu = 1e-6;
    for _ in 0..200 {
        let h = 1e-7;
        let jac: Vec<Vec<f64>> = (0..6)
            .map(|k| {
                let mut xp = x;
                xp[k] += h;
                res(&xp).iter().zip(&r).map(|(a, b)| (a - b) / h).collect()
            })
            .collect();
        let mut a = vec![vec![0.0; 6]; 6];
        let mut b = vec![0.0; 6];
        for i in 0..6 {
            for j in 0..6 {
                a[i][j] = jac[i].iter().zip(&jac[j]).map(|(p, q)| p * q).sum();
            }
            a[i][i] += mu * (1.0 + a[i][i]);
            b[i] = -jac[i].iter().zip(&r).map(|(p, q)| p * q).sum::<f64>();
        }
        let Some(dx) = solve_lin(&mut a, &mut b) else { break };
        let mut xn = x;
        for k in 0..6 {
            xn[k] += dx[k];
        }
        let rn = res(&xn);
        let fnew = rn.iter().map(|v| v * v).sum::<f64>();
        if fnew < f {
            let small = f - fnew < 1e-14 * f.max(1e-30);
            x = xn;
            r = rn;
            f = fnew;
            mu = (mu * 0.3).max(1e-12);
            if small || f < 1e-28 {
                break;
            }
        } else {
            mu *= 10.0;
            if mu > 1e8 {
                break;
            }
        }
    }
    mk(&x)
}

fn build_raw(p: &[V3], betas: &[f64], target_nv: f64) -> Result<TruePoly, String> {
    let n = p.len();
    let p = p.to_vec();
    let neg = betas.iter().map(|b| 2.0 * PI * (1.0 - b)).sum::<f64>() > 4.0 * PI;
    let (pv, radii, defect_err, hull, groups) = if neg {
        let (pv, simplices, radii, err) = solve_star(&p, betas)?;
        let normals: Vec<V3> = simplices
            .iter()
            .map(|t| {
                let nn = geom::cross(geom::sub(pv[t[1]], pv[t[0]]), geom::sub(pv[t[2]], pv[t[0]]));
                geom::scale(nn, 1.0 / geom::norm(nn))
            })
            .collect();
        let offsets = simplices.iter().zip(&normals).map(|(t, nn)| geom::dot(*nn, pv[t[0]])).collect();
        // every triangle is its own facet (no coplanar merge: non-adjacent faces may share a plane)
        let groups: Vec<Vec<usize>> = (0..simplices.len()).map(|t| vec![t]).collect();
        (pv, radii, err, crate::hull::Hull { simplices, normals, offsets }, groups)
    } else {
        let cor = extra_corners(&p, 4);
        let (pv, radii, defect_err) = solve_polytope(&p, betas, &cor)?;
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
        (pv, radii, defect_err, hull, groups)
    };
    let nvp = pv.len();
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
    let mut y: Vec<V3> = x.iter().map(|q| geom::scale(*q, 1.0 / geom::norm(*q))).collect();
    if neg {
        // star polyhedron: map each face barycentrically onto the flat triangle of its unit
        // directions, then project (long spikes would otherwise crowd into tiny caps)
        let dir: Vec<V3> = pv.iter().map(|q| geom::scale(*q, 1.0 / geom::norm(*q))).collect();
        for (v, q) in x.iter().enumerate() {
            let mut best = (f64::NEG_INFINITY, [0.0; 3], 0);
            for (ti, t) in hull.simplices.iter().enumerate() {
                let (a, b, c) = (pv[t[0]], pv[t[1]], pv[t[2]]);
                let nn = geom::cross(geom::sub(b, a), geom::sub(c, a));
                let n2 = geom::dot(nn, nn);
                let w = [
                    geom::dot(nn, geom::cross(geom::sub(c, b), geom::sub(*q, b))) / n2,
                    geom::dot(nn, geom::cross(geom::sub(a, c), geom::sub(*q, c))) / n2,
                    geom::dot(nn, geom::cross(geom::sub(b, a), geom::sub(*q, a))) / n2,
                ];
                let mn = w[0].min(w[1]).min(w[2]) - (geom::dot(nn, geom::sub(*q, a)) / n2.sqrt()).abs();
                if mn > best.0 {
                    best = (mn, w, ti);
                }
            }
            let t = hull.simplices[best.2];
            let w = best.1;
            let m = geom::add(geom::add(geom::scale(dir[t[0]], w[0]), geom::scale(dir[t[1]], w[1])), geom::scale(dir[t[2]], w[2]));
            y[v] = geom::scale(m, 1.0 / geom::norm(m));
        }
    }
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
    let u_off = fit_offset(&mesh, &l0)?;
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
    Ok(TruePoly { mesh, flow_mesh, x, corners, crease_edges, l0, u_off, radii, defect_err, fit_resid: f64::NAN, match_err: f64::NAN, match_iters: 0, warning: None })
}

/// Number of intersecting pairs among faces that share no vertex (segment-vs-triangle test with
/// a uniform-grid broad phase).  0 for an embedded closed surface.
pub fn self_intersections(x: &[V3], f: &[[usize; 3]]) -> usize {
    let bb: Vec<(V3, V3)> = f
        .iter()
        .map(|t| {
            let mut lo = x[t[0]];
            let mut hi = x[t[0]];
            for &v in &t[1..] {
                for c in 0..3 {
                    lo[c] = lo[c].min(x[v][c]);
                    hi[c] = hi[c].max(x[v][c]);
                }
            }
            (lo, hi)
        })
        .collect();
    let mut glo = bb[0].0;
    let mut ghi = bb[0].1;
    for (lo, hi) in &bb {
        for c in 0..3 {
            glo[c] = glo[c].min(lo[c]);
            ghi[c] = ghi[c].max(hi[c]);
        }
    }
    let cell = bb.iter().map(|(lo, hi)| geom::norm(geom::sub(*hi, *lo))).sum::<f64>() / bb.len() as f64 * 1.5;
    let key = |q: V3| -> [i64; 3] { [0, 1, 2].map(|c| ((q[c] - glo[c]) / cell).floor() as i64) };
    let mut grid: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
    for (i, (lo, hi)) in bb.iter().enumerate() {
        let (a, b) = (key(*lo), key(*hi));
        for i0 in a[0]..=b[0] {
            for i1 in a[1]..=b[1] {
                for i2 in a[2]..=b[2] {
                    grid.entry([i0, i1, i2]).or_default().push(i);
                }
            }
        }
    }
    let seg_tri = |p0: V3, p1: V3, t: [usize; 3]| -> bool {
        let (a, b, c) = (x[t[0]], x[t[1]], x[t[2]]);
        let n = geom::cross(geom::sub(b, a), geom::sub(c, a));
        let (d0, d1) = (geom::dot(n, geom::sub(p0, a)), geom::dot(n, geom::sub(p1, a)));
        let eps = 1e-12 * geom::dot(n, n).sqrt() * cell;
        if (d0 > eps && d1 > eps) || (d0 < -eps && d1 < -eps) || (d0 - d1).abs() < 1e-300 {
            return false;
        }
        let s = d0 / (d0 - d1);
        if !(-1e-9..=1.0 + 1e-9).contains(&s) {
            return false;
        }
        let q = geom::add(p0, geom::scale(geom::sub(p1, p0), s));
        let inside = |u: V3, v: V3| geom::dot(n, geom::cross(geom::sub(v, u), geom::sub(q, u))) > 1e-12 * geom::dot(n, n);
        inside(a, b) && inside(b, c) && inside(c, a)
    };
    let mut seen = std::collections::HashSet::new();
    let mut cnt = 0;
    for v in grid.values() {
        for (k, &i) in v.iter().enumerate() {
            for &j in &v[k + 1..] {
                let (ti, tj) = (f[i], f[j]);
                if ti.iter().any(|a| tj.contains(a)) || !seen.insert((i.min(j), i.max(j))) {
                    continue;
                }
                let (bi, bj) = (bb[i], bb[j]);
                if (0..3).any(|c| bi.1[c] < bj.0[c] || bj.1[c] < bi.0[c]) {
                    continue;
                }
                let hit = (0..3).any(|e| seg_tri(x[ti[e]], x[ti[(e + 1) % 3]], tj)) || (0..3).any(|e| seg_tri(x[tj[e]], x[tj[(e + 1) % 3]], ti));
                if hit {
                    cnt += 1;
                }
            }
        }
    }
    cnt
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::defects_and_areas;
    #[test]
    fn conformal_vs_preview() {
        let pos = [[0.7, 0.0], [-0.7, 0.0], [0.0, 0.7], [0.0, -0.7]];
        let betas = [0.6, 0.7, 0.5, 0.8];
        let lvl = seedmesh::page_level("coarse").unwrap();
        let pv = build_preview(&pos, &betas, lvl).unwrap();
        let tp = build(&pos, &betas, lvl).unwrap();
        assert!(tp.warning.is_none());
        assert!(tp.match_err < 2e-3, "match err {}", tp.match_err);
        assert!(tp.fit_resid < 5e-3, "fit resid {}", tp.fit_resid);
        // preview: radial positions, cones exactly at the user's points, old LS offset
        let target = seedmesh::stereographic_to_sphere(&pos);
        for (i, t) in target.iter().enumerate() {
            assert!(geom::norm(geom::sub(pv.mesh.y[i], *t)) < 1e-12);
        }
        assert!(pv.u_off.iter().all(|x| x.is_finite()));
    }
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
    #[test]
    fn negative_star_polyhedra() {
        let cfgs: Vec<(Vec<f64>, Vec<[f64; 2]>)> = vec![
            (vec![0.2, 0.3, 0.4], vec![[0.7, 0.0], [-0.35, 0.606218], [-0.35, -0.606218]]),
            (vec![0.05, 0.05, 0.05], vec![[0.7, 0.0], [-0.35, 0.606218], [-0.35, -0.606218]]),
            (vec![0.3; 4], vec![[0.7, 0.0], [0.0, 0.7], [-0.7, 0.0], [0.0, -0.7]]),
            (vec![0.25, 0.45, 0.35, 0.3], vec![[0.6, 0.1], [-0.2, 0.7], [-0.6, -0.3], [0.3, -0.5]]),
            (vec![0.1, 0.1, 0.1, 0.95], vec![[0.7, 0.0], [0.0, 0.7], [-0.7, 0.0], [0.0, 0.1]]),
            (vec![0.2; 5], vec![[0.7, 0.0], [0.216, 0.666], [-0.566, 0.411], [-0.566, -0.411], [0.216, -0.666]]),
            (vec![0.1, 0.15, 0.2, 0.25, 0.3, 0.35], vec![[0.7, 0.0], [0.35, 0.606], [-0.35, 0.606], [-0.7, 0.0], [-0.35, -0.606], [0.35, -0.606]]),
        ];
        for (betas, pos) in cfgs {
            let tp = build_preview(&pos, &betas, seedmesh::page_level("coarse").unwrap()).unwrap();
            let g = defects_and_areas(&tp.mesh.f, &tp.mesh.opp, &tp.l0, tp.mesh.nv());
            let n = betas.len();
            for (i, b) in betas.iter().enumerate() {
                assert!((g.k[i] - 2.0 * PI * (1.0 - b)).abs() < 1e-12, "cone {i} of {betas:?}");
            }
            assert!((g.k.iter().sum::<f64>() - 4.0 * PI).abs() < 1e-11);
            assert!(tp.corners.iter().map(|&c| g.k[c]).sum::<f64>() < 0.0);
            for v in n + tp.corners.len()..tp.mesh.nv() {
                assert!(g.k[v].abs() < 1e-11);
            }
            assert_eq!(self_intersections(&tp.x, &tp.mesh.f), 0, "{betas:?}");
        }
    }
    #[test]
    fn self_intersection_detects_overlap() {
        let x = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.2, 0.2, -0.5], [0.3, 0.2, 0.5], [0.2, 0.3, 0.5]];
        assert_eq!(self_intersections(&x, &[[0, 1, 2], [3, 4, 5]]), 1);
        let y: Vec<V3> = x.iter().map(|q| [q[0], q[1], if q[2] == 0.0 { 0.0 } else { q[2] + 1.0 }]).collect();
        assert_eq!(self_intersections(&y, &[[0, 1, 2], [3, 4, 5]]), 0);
    }
}
