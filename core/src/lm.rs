//! Levenberg-Marquardt solvers (port of `embeddings/fast_realizer.fade_lm` and of the
//! `scipy.optimize.least_squares(method="trf", x_scale="jac", jac_sparsity=..)` call used by
//! `true_flow.polytope_fit`), plus the sparse SPD solver they share.
//!
//! The sparse solver is an envelope (skyline) Cholesky with a reverse Cuthill-McKee ordering:
//! for these sphere meshes (a few hundred to a few thousand vertices) the envelope is small and
//! the code has no dependencies, so it builds for `wasm32-unknown-unknown` without threads.
//! Python solves the same SPD systems with SuperLU (MMD ordering); the solution is the same.

use crate::realize::geom::{self, V3};

// ------------------------------------------------------------------ sparse SPD
/// Symmetric positive (semi)definite matrix in envelope storage, permuted by RCM.
#[derive(Clone, Debug)]
pub struct Skyline {
    pub n: usize,
    iperm: Vec<usize>, // original -> permuted
    perm: Vec<usize>,  // permuted -> original
    first: Vec<usize>, // first stored column of permuted row i
    off: Vec<usize>,   // start of row i in vals
    pub vals: Vec<f64>,
}

impl Skyline {
    /// Pattern from residual-row supports: every pair of indices in one row is coupled.
    pub fn from_rows(n: usize, rows: &[Vec<usize>]) -> Skyline {
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        for r in rows {
            for &a in r {
                for &b in r {
                    if a != b {
                        adj[a].push(b);
                    }
                }
            }
        }
        for a in adj.iter_mut() {
            a.sort_unstable();
            a.dedup();
        }
        let perm = rcm(&adj);
        let mut iperm = vec![0; n];
        for (i, &p) in perm.iter().enumerate() {
            iperm[p] = i;
        }
        let mut first: Vec<usize> = (0..n).collect();
        for a in 0..n {
            let ia = iperm[a];
            for &b in &adj[a] {
                let ib = iperm[b];
                if ib < first[ia] {
                    first[ia] = ib;
                }
            }
        }
        let mut off = vec![0; n + 1];
        for i in 0..n {
            off[i + 1] = off[i] + (i - first[i] + 1);
        }
        let nn = off[n];
        Skyline { n, iperm, perm, first, off, vals: vec![0.0; nn] }
    }

    pub fn clear(&mut self) {
        self.vals.iter_mut().for_each(|v| *v = 0.0);
    }

    #[inline]
    fn pos(&self, a: usize, b: usize) -> usize {
        let (i, j) = {
            let (x, y) = (self.iperm[a], self.iperm[b]);
            if x >= y { (x, y) } else { (y, x) }
        };
        debug_assert!(j >= self.first[i]);
        self.off[i] + j - self.first[i]
    }

    /// A += v v^T restricted to the support `idx` (duplicates in idx must be merged by the caller).
    #[inline]
    pub fn add_outer(&mut self, idx: &[usize], v: &[f64]) {
        for p in 0..idx.len() {
            for q in 0..=p {
                let k = self.pos(idx[p], idx[q]);
                self.vals[k] += v[p] * v[q];
            }
        }
    }

    pub fn diagonal(&self) -> Vec<f64> {
        (0..self.n).map(|a| self.vals[self.pos(a, a)]).collect()
    }

    pub fn add_diag(&mut self, d: &[f64]) {
        for a in 0..self.n {
            let k = self.pos(a, a);
            self.vals[k] += d[a];
        }
    }

    /// In-place Cholesky; false if a pivot is not positive.
    pub fn factor(&mut self) -> bool {
        let n = self.n;
        for i in 0..n {
            let fi = self.first[i];
            for j in fi..=i {
                let fj = self.first[j];
                let k0 = fi.max(fj);
                let mut s = self.vals[self.off[i] + j - fi];
                let (ri, rj) = (self.off[i] - fi, self.off[j] - fj);
                for k in k0..j {
                    s -= self.vals[ri + k] * self.vals[rj + k];
                }
                if j < i {
                    s /= self.vals[self.off[j] + j - fj];
                    self.vals[ri + j] = s;
                } else {
                    if !(s > 0.0) || !s.is_finite() {
                        return false;
                    }
                    self.vals[ri + i] = s.sqrt();
                }
            }
        }
        true
    }

    /// Solve with the factor (call after `factor`).
    pub fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = self.n;
        let mut y: Vec<f64> = (0..n).map(|i| b[self.perm[i]]).collect();
        for i in 0..n {
            let fi = self.first[i];
            let r = self.off[i] - fi;
            let mut s = y[i];
            for k in fi..i {
                s -= self.vals[r + k] * y[k];
            }
            y[i] = s / self.vals[r + i];
        }
        for i in (0..n).rev() {
            let fi = self.first[i];
            let r = self.off[i] - fi;
            y[i] /= self.vals[r + i];
            let yi = y[i];
            for k in fi..i {
                y[k] -= self.vals[r + k] * yi;
            }
        }
        let mut x = vec![0.0; n];
        for i in 0..n {
            x[self.perm[i]] = y[i];
        }
        x
    }
}

/// Reverse Cuthill-McKee ordering (per connected component, min-degree start).
fn rcm(adj: &[Vec<usize>]) -> Vec<usize> {
    let n = adj.len();
    let mut seen = vec![false; n];
    let mut order = Vec::with_capacity(n);
    let mut byd: Vec<usize> = (0..n).collect();
    byd.sort_by_key(|&v| (adj[v].len(), v));
    for &s in &byd {
        if seen[s] {
            continue;
        }
        // pseudo-peripheral start: BFS twice
        let mut start = s;
        for _ in 0..2 {
            let lv = bfs_levels(adj, start, &seen);
            let maxl = lv.iter().filter_map(|x| *x).max().unwrap_or(0);
            let far = (0..n)
                .filter(|&v| lv[v] == Some(maxl))
                .min_by_key(|&v| (adj[v].len(), v))
                .unwrap_or(start);
            start = far;
        }
        seen[start] = true;
        let mut q = std::collections::VecDeque::new();
        q.push_back(start);
        while let Some(v) = q.pop_front() {
            order.push(v);
            let mut nb: Vec<usize> = adj[v].iter().copied().filter(|&w| !seen[w]).collect();
            nb.sort_by_key(|&w| (adj[w].len(), w));
            for w in nb {
                seen[w] = true;
                q.push_back(w);
            }
        }
    }
    order.reverse();
    order
}

fn bfs_levels(adj: &[Vec<usize>], s: usize, blocked: &[bool]) -> Vec<Option<usize>> {
    let mut lv = vec![None; adj.len()];
    lv[s] = Some(0);
    let mut q = std::collections::VecDeque::new();
    q.push_back(s);
    while let Some(v) = q.pop_front() {
        let l = lv[v].unwrap();
        for &w in &adj[v] {
            if lv[w].is_none() && !blocked[w] {
                lv[w] = Some(l + 1);
                q.push_back(w);
            }
        }
    }
    lv
}

// ------------------------------------------------------------------ fade_lm
/// Precomputed combinatorics for [`fade_lm`] on one triangulation.
#[derive(Clone, Debug)]
pub struct FadeSetup {
    pub nv: usize,
    pub faces: Vec<[usize; 3]>,
    pub edges: Vec<[usize; 2]>,
    /// bending quads (a, b, c, d) of edges not touching a cone (`bend_gn.edge_quads` order)
    pub quads: Vec<[usize; 4]>,
    sky: Skyline,
}

impl FadeSetup {
    pub fn new(faces: &[[usize; 3]], edges: &[[usize; 2]], cones: &[usize], nv: usize) -> FadeSetup {
        let mut cs = vec![false; nv];
        for &c in cones {
            cs[c] = true;
        }
        let quads: Vec<[usize; 4]> = geom::edge_quads(faces)
            .into_iter()
            .filter(|q| !(cs[q.0[0]] || cs[q.0[1]]))
            .map(|q| q.1)
            .collect();
        let mut rows: Vec<Vec<usize>> = Vec::new();
        for e in edges {
            rows.push(vec3(&[e[0], e[1]]));
        }
        for q in &quads {
            rows.push(vec3(q));
        }
        let sky = Skyline::from_rows(3 * nv, &rows);
        FadeSetup { nv, faces: faces.to_vec(), edges: edges.to_vec(), quads, sky }
    }
}

fn vec3(vs: &[usize]) -> Vec<usize> {
    let mut out = Vec::with_capacity(3 * vs.len());
    for &v in vs {
        for c in 0..3 {
            out.push(3 * v + c);
        }
    }
    out
}

/// Result info of [`fade_lm`].
#[derive(Clone, Debug)]
pub struct FadeInfo {
    pub n_iter: usize,
    pub len_rel_max: f64,
}

#[inline]
fn theta(a: V3, b: V3, c: V3, d: V3) -> f64 {
    let n1 = geom::cross(geom::sub(b, a), geom::sub(c, a));
    let n2 = geom::cross(geom::sub(a, b), geom::sub(d, b));
    let e = geom::sub(b, a);
    let en = geom::norm(e) + 1e-300;
    let e = geom::scale(e, 1.0 / en);
    let s = geom::dot(geom::cross(n1, n2), e);
    let c_ = geom::dot(n1, n2);
    s.atan2(c_)
}

fn phi_of(s: &FadeSetup, x: &[V3], l: &[f64], lam: f64, sgn: f64) -> f64 {
    let mut p = 0.0;
    for (e, &le) in s.edges.iter().zip(l) {
        let r = (geom::norm(geom::sub(x[e[0]], x[e[1]])) - le) / le;
        p += r * r;
    }
    if lam > 0.0 {
        let mut b = 0.0;
        for q in &s.quads {
            let t = sgn * theta(x[q[0]], x[q[1]], x[q[2]], x[q[3]]);
            b += t * t;
        }
        p += lam * b;
    }
    p
}

/// Assemble J^T J into `s.sky` and return g = J^T r (bg._len_jac + bg._dihed_jac, central FD dihedrals).
fn system(s: &mut FadeSetup, x: &[V3], l: &[f64], lam: f64, sgn: f64) -> Vec<f64> {
    let mut g = vec![0.0; 3 * s.nv];
    s.sky.clear();
    let mut idx = [0usize; 12];
    let mut val = [0.0f64; 12];
    for (e, &le) in s.edges.iter().zip(l) {
        let d = geom::sub(x[e[0]], x[e[1]]);
        let ell = geom::norm(d);
        let r = (ell - le) / le;
        for c in 0..3 {
            let gc = d[c] / ell / le;
            idx[c] = 3 * e[0] + c;
            val[c] = gc;
            idx[3 + c] = 3 * e[1] + c;
            val[3 + c] = -gc;
        }
        s.sky.add_outer(&idx[..6], &val[..6]);
        for k in 0..6 {
            g[idx[k]] += val[k] * r;
        }
    }
    if lam > 0.0 {
        let sl = lam.sqrt();
        for qi in 0..s.quads.len() {
            let q = s.quads[qi];
            let p0 = [x[q[0]], x[q[1]], x[q[2]], x[q[3]]];
            let h = 1e-6 * (geom::norm(geom::sub(p0[0], p0[1])) + 1e-12);
            for k in 0..4 {
                for c in 0..3 {
                    let mut pp = p0;
                    let mut pm = p0;
                    pp[k][c] += h;
                    pm[k][c] -= h;
                    let v = sgn * (theta(pp[0], pp[1], pp[2], pp[3]) - theta(pm[0], pm[1], pm[2], pm[3])) / (2.0 * h);
                    idx[3 * k + c] = 3 * q[k] + c;
                    val[3 * k + c] = sl * v;
                }
            }
            s.sky.add_outer(&idx, &val);
            let r = sl * sgn * theta(p0[0], p0[1], p0[2], p0[3]);
            for k in 0..12 {
                g[idx[k]] += val[k] * r;
            }
        }
    }
    g
}

/// `fast_realizer.fade_lm(F, E, L, cones, X0, stages=.., reg="dihedral", tol=1e-6, mu0=1e-3)`:
/// LM on sum ((l - L)/L)^2 + lam * sum(dihedral^2), lam continued over `stages` = [(lam, max_iters)].
pub fn fade_lm(s: &mut FadeSetup, l: &[f64], x0: &[V3], stages: &[(f64, usize)]) -> (Vec<V3>, FadeInfo) {
    let tol = 1e-6;
    let mu0 = 1e-3;
    let mut x = x0.to_vec();
    let sgn = if geom::signed_volume(&x, &s.faces) > 0.0 { 1.0 } else { -1.0 };
    let mut it = 0;
    for &(lam, nit) in stages {
        let mut mu = mu0;
        let mut phi = phi_of(s, &x, l, lam, sgn);
        for _ in 0..nit {
            let g = system(s, &x, l, lam, sgn);
            let dg = s.sky.diagonal();
            let dmax = dg.iter().cloned().fold(f64::MIN, f64::max);
            let base = s.sky.vals.clone();
            let mut accepted: Option<(Vec<V3>, f64)> = None;
            for _ in 0..10 {
                s.sky.vals.copy_from_slice(&base);
                let add: Vec<f64> = dg.iter().map(|d| mu * (d + 1e-9 * dmax)).collect();
                s.sky.add_diag(&add);
                let dx = if s.sky.factor() { Some(s.sky.solve(&g)) } else { None };
                let dx = match dx {
                    Some(d) if d.iter().all(|v| v.is_finite()) => d,
                    _ => {
                        mu *= 8.0;
                        continue;
                    }
                };
                let xn: Vec<V3> = x
                    .iter()
                    .enumerate()
                    .map(|(v, p)| [p[0] - dx[3 * v], p[1] - dx[3 * v + 1], p[2] - dx[3 * v + 2]])
                    .collect();
                if geom::signed_volume(&xn, &s.faces) * sgn <= 0.0 {
                    mu *= 8.0;
                    continue;
                }
                let pn = phi_of(s, &xn, l, lam, sgn);
                if pn < phi {
                    accepted = Some((xn, pn));
                    break;
                }
                mu *= 8.0;
            }
            s.sky.vals.copy_from_slice(&base);
            let Some((xn, pn)) = accepted else { break };
            let dphi = phi - pn;
            x = xn;
            phi = pn;
            mu = (mu / 4.0).max(1e-7);
            it += 1;
            if dphi < tol * phi.max(1e-30) {
                break;
            }
        }
    }
    let rel = geom::len_rel(&x, &s.edges, l);
    let m = rel.iter().cloned().fold(0.0, f64::max);
    (x, FadeInfo { n_iter: it, len_rel_max: m })
}

// ------------------------------------------------------------------ least squares
/// Sparse nonlinear least-squares problem: residuals with an analytic sparse Jacobian.
pub trait LsqProblem {
    fn n_params(&self) -> usize;
    fn n_resid(&self) -> usize;
    /// Support (parameter indices) of each residual row; fixed for the problem.
    fn row_support(&self, i: usize) -> Vec<usize>;
    fn residuals(&self, p: &[f64]) -> Vec<f64>;
    /// Jacobian rows: (indices, values) per residual (indices merged, no duplicates).
    fn jac_rows(&self, p: &[f64]) -> Vec<(Vec<usize>, Vec<f64>)>;
}

#[derive(Clone, Debug)]
pub struct LsqResult {
    pub x: Vec<f64>,
    pub cost: f64,
    pub nfev: usize,
    pub status: i32,
}

/// Replacement for `least_squares(fun, x0, jac_sparsity=S, method="trf", max_nfev, x_scale="jac")`
/// (loss linear, ftol = xtol = gtol = 1e-8): a trust-region Levenberg-Marquardt with the
/// Jacobian column-norm scaling of `x_scale="jac"` (scale only grows, as in scipy), counting
/// residual evaluations against `max_nfev`.  Exact normal-equation steps (sparse Cholesky)
/// instead of scipy's LSMR subproblem: same minimiser family, not bit-identical iterates.
pub fn least_squares_jac_scaled<P: LsqProblem>(prob: &P, x0: &[f64], max_nfev: usize) -> LsqResult {
    let (ftol, xtol, gtol) = (1e-8, 1e-8, 1e-8);
    let n = prob.n_params();
    let rows: Vec<Vec<usize>> = (0..prob.n_resid()).map(|i| prob.row_support(i)).collect();
    let mut sky = Skyline::from_rows(n, &rows);
    let mut x = x0.to_vec();
    let mut r = prob.residuals(&x);
    let mut nfev = 1;
    let mut cost = 0.5 * r.iter().map(|v| v * v).sum::<f64>();
    let mut scale2 = vec![0.0f64; n]; // (1 / x_scale)^2 = max running column norm^2
    let mut lam = -1.0;
    let mut status = 0;
    'outer: while nfev < max_nfev {
        let jr = prob.jac_rows(&x);
        sky.clear();
        let mut g = vec![0.0; n];
        let mut col2 = vec![0.0; n];
        for (i, (idx, val)) in jr.iter().enumerate() {
            sky.add_outer(idx, val);
            for (k, &j) in idx.iter().enumerate() {
                g[j] += val[k] * r[i];
                col2[j] += val[k] * val[k];
            }
        }
        for j in 0..n {
            scale2[j] = scale2[j].max(col2[j]);
            if scale2[j] == 0.0 {
                scale2[j] = 1.0;
            }
        }
        let gn = (0..n).map(|j| (g[j] / scale2[j].sqrt()).abs()).fold(0.0, f64::max);
        if gn < gtol {
            status = 1;
            break;
        }
        if lam < 0.0 {
            lam = 1e-3;
        }
        let base = sky.vals.clone();
        loop {
            if nfev >= max_nfev {
                break 'outer;
            }
            sky.vals.copy_from_slice(&base);
            let add: Vec<f64> = scale2.iter().map(|s| lam * s).collect();
            sky.add_diag(&add);
            if !sky.factor() {
                lam *= 10.0;
                continue;
            }
            let dx: Vec<f64> = sky.solve(&g).iter().map(|v| -v).collect();
            // predicted reduction of the quadratic model
            sky.vals.copy_from_slice(&base);
            let mut jdx = vec![0.0; r.len()];
            for (i, (idx, val)) in jr.iter().enumerate() {
                jdx[i] = idx.iter().zip(val).map(|(&j, v)| v * dx[j]).sum();
            }
            let pred = -(0..r.len()).map(|i| jdx[i] * r[i] + 0.5 * jdx[i] * jdx[i]).sum::<f64>();
            let xn: Vec<f64> = x.iter().zip(&dx).map(|(a, b)| a + b).collect();
            let rn = prob.residuals(&xn);
            nfev += 1;
            let cn = 0.5 * rn.iter().map(|v| v * v).sum::<f64>();
            let act = cost - cn;
            let rho = if pred > 0.0 { act / pred } else { -1.0 };
            if rho > 0.75 {
                lam = (lam / 3.0).max(1e-12);
            } else if rho < 0.25 {
                lam *= 4.0;
            }
            if act > 0.0 && cn.is_finite() {
                let dxn = dx.iter().map(|v| v * v).sum::<f64>().sqrt();
                let xnorm = x.iter().map(|v| v * v).sum::<f64>().sqrt();
                x = xn;
                r = rn;
                let old = cost;
                cost = cn;
                if act < ftol * old {
                    status = 2;
                    break 'outer;
                }
                if dxn < xtol * (xtol + xnorm) {
                    status = 3;
                    break 'outer;
                }
                break;
            }
            if lam > 1e16 {
                break 'outer;
            }
        }
    }
    LsqResult { x, cost, nfev, status }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn skyline_solves_spd() {
        // 1D Laplacian + identity
        let n = 7;
        let rows: Vec<Vec<usize>> = (0..n - 1).map(|i| vec![i, i + 1]).collect();
        let mut s = Skyline::from_rows(n, &rows);
        for i in 0..n - 1 {
            s.add_outer(&[i, i + 1], &[1.0, -1.0]);
        }
        s.add_diag(&vec![1.0; n]);
        let b: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let a = s.clone();
        assert!(s.factor());
        let x = s.solve(&b);
        // residual check: A x
        let mut ax = vec![0.0; n];
        for i in 0..n {
            ax[i] += x[i];
        }
        for i in 0..n - 1 {
            let d = x[i] - x[i + 1];
            ax[i] += d;
            ax[i + 1] -= d;
        }
        let _ = a;
        for i in 0..n {
            assert!((ax[i] - b[i]).abs() < 1e-12);
        }
    }
}
