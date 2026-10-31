//! Smooth raster of the conformal factor u for the page's north/south discs: port of
//! `viz/w_interp.smooth_w_grid` (cone model + cubic RBF of the smooth remainder, linear
//! polynomial tail, = `scipy.interpolate.RBFInterpolator(kernel="cubic", degree=1)`) and of
//! `viewer/backend/flow_true.conformal_grid` (which vertices and singularities are used).
//!
//! The RBF matrix depends only on the point set, so it is LU-factored once ([`CfGrid::new`]) and
//! reused for every frame on the same mesh ([`CfGrid::grid`]).  Output: `n_lat * n_lon` row-major,
//! rows lat in [0, pi] (north first), cols lon in [-pi, pi).

use std::f64::consts::PI;

type V3 = [f64; 3];

pub const GRID_LAT: usize = 64;
pub const GRID_LON: usize = 128;

fn dist(a: V3, b: V3) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

fn unit(p: V3) -> V3 {
    let n = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
    [p[0] / n, p[1] / n, p[2] / n]
}

/// `w_interp.raster_dirs`.
pub fn raster_dirs(n_lat: usize, n_lon: usize) -> Vec<V3> {
    let mut out = Vec::with_capacity(n_lat * n_lon);
    for i in 0..n_lat {
        let lat = (i as f64 + 0.5) * PI / n_lat as f64;
        for j in 0..n_lon {
            let lon = -PI + (j as f64 + 0.5) * 2.0 * PI / n_lon as f64;
            out.push([lat.sin() * lon.cos(), lat.sin() * lon.sin(), lat.cos()]);
        }
    }
    out
}

/// `w_interp.cone_model`: sum_i (beta_i - 1) log max(|X - P_i|, floor).
fn cone_model(x: V3, cu: &[V3], betas: &[f64], floor: f64) -> f64 {
    cu.iter().zip(betas).map(|(p, b)| (b - 1.0) * dist(x, *p).max(floor.max(1e-300)).ln()).sum()
}

/// Dense LU with partial pivoting.
struct Lu {
    n: usize,
    a: Vec<f64>,
    piv: Vec<usize>,
}

impl Lu {
    fn factor(n: usize, mut a: Vec<f64>) -> Option<Lu> {
        let mut piv = vec![0; n];
        for k in 0..n {
            let mut p = k;
            for i in k + 1..n {
                if a[i * n + k].abs() > a[p * n + k].abs() {
                    p = i;
                }
            }
            if a[p * n + k] == 0.0 {
                return None;
            }
            piv[k] = p;
            if p != k {
                for j in 0..n {
                    a.swap(k * n + j, p * n + j);
                }
            }
            let d = a[k * n + k];
            let (top, bot) = a.split_at_mut((k + 1) * n);
            let rk = &top[k * n..];
            for i in 0..n - k - 1 {
                let row = &mut bot[i * n..(i + 1) * n];
                let f = row[k] / d;
                row[k] = f;
                if f != 0.0 {
                    for j in k + 1..n {
                        row[j] -= f * rk[j];
                    }
                }
            }
        }
        Some(Lu { n, a, piv })
    }
    fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = self.n;
        let mut x = b.to_vec();
        for k in 0..n {
            x.swap(k, self.piv[k]);
        }
        for i in 0..n {
            let mut s = x[i];
            for j in 0..i {
                s -= self.a[i * n + j] * x[j];
            }
            x[i] = s;
        }
        for i in (0..n).rev() {
            let mut s = x[i];
            for j in i + 1..n {
                s -= self.a[i * n + j] * x[j];
            }
            x[i] = s / self.a[i * n + i];
        }
        x
    }
}

/// Factored `smooth_w_grid` for one point set / singularity set.
pub struct CfGrid {
    /// mesh vertex ids used as RBF centres (non-cone vertices of the `keep` subset)
    centres: Vec<usize>,
    pts: Vec<V3>,
    model_v: Vec<f64>,
    cu: Vec<V3>,
    betas: Vec<f64>,
    shift: V3,
    scale: V3,
    lu: Lu,
    n_lat: usize,
    n_lon: usize,
    /// per raster cell: cone model (with floor) — constant across frames
    model_dirs: Vec<f64>,
    dirs: Vec<V3>,
}

impl CfGrid {
    /// `conformal_grid` setup: `y` unit-sphere mesh vertices, `sing` singular vertex ids (cones,
    /// plus the polytope seed corners at t = 0), `betas` per singularity.
    pub fn new(y: &[V3], sing: &[usize], betas: &[f64], n_lat: usize, n_lon: usize) -> Result<CfGrid, String> {
        let nv = y.len();
        let mut keep: Vec<usize> = (0..nv).collect();
        if nv > 1500 {
            let rest: Vec<usize> = (0..nv).filter(|v| !sing.contains(v)).collect();
            let stride = (rest.len() as f64 / (1500 - sing.len()) as f64).ceil() as usize;
            keep = sing.to_vec();
            keep.extend(rest.iter().step_by(stride.max(1)));
        }
        let unit_all: Vec<V3> = keep.iter().map(|&v| unit(y[v])).collect();
        let cu: Vec<V3> = sing.iter().map(|&v| unit(y[v])).collect();
        // smooth_w_grid(..., cone_ids=range(len(sing))) indexes the keep subset
        let mut is_cone = vec![false; keep.len()];
        for i in 0..sing.len().min(keep.len()) {
            is_cone[i] = true;
        }
        for (i, x) in unit_all.iter().enumerate() {
            if cu.iter().any(|p| dist(*x, *p) < 1e-9) {
                is_cone[i] = true;
            }
        }
        let idx: Vec<usize> = (0..keep.len()).filter(|&i| !is_cone[i]).collect();
        if idx.len() < 5 {
            return Err("not enough non-cone vertices for smooth interpolation".into());
        }
        let centres: Vec<usize> = idx.iter().map(|&i| keep[i]).collect();
        let pts: Vec<V3> = idx.iter().map(|&i| unit_all[i]).collect();
        let model_v: Vec<f64> = pts.iter().map(|&x| cone_model(x, &cu, betas, 0.0)).collect();
        // scipy: polynomial on (x - shift) / scale, shift/scale from the data bounding box
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for p in &pts {
            for c in 0..3 {
                lo[c] = lo[c].min(p[c]);
                hi[c] = hi[c].max(p[c]);
            }
        }
        let mut shift = [0.0; 3];
        let mut scale = [1.0; 3];
        for c in 0..3 {
            shift[c] = 0.5 * (lo[c] + hi[c]);
            let s = 0.5 * (hi[c] - lo[c]);
            scale[c] = if s == 0.0 { 1.0 } else { s };
        }
        let n = pts.len();
        let m = n + 4;
        let mut a = vec![0.0; m * m];
        for i in 0..n {
            for j in 0..n {
                a[i * m + j] = dist(pts[i], pts[j]).powi(3);
            }
            let q = poly(pts[i], shift, scale);
            for c in 0..4 {
                a[i * m + n + c] = q[c];
                a[(n + c) * m + i] = q[c];
            }
        }
        let lu = Lu::factor(m, a).ok_or("singular RBF system")?;
        let floor = 0.5 * PI / n_lat as f64;
        let dirs = raster_dirs(n_lat, n_lon);
        let model_dirs = dirs.iter().map(|&x| cone_model(x, &cu, betas, floor)).collect();
        Ok(CfGrid { centres, pts, model_v, cu, betas: betas.to_vec(), shift, scale, lu, n_lat, n_lon, model_dirs, dirs })
    }

    /// Raster of the vertex field `u` (len nv), `n_lat * n_lon` row-major.
    pub fn grid(&self, u: &[f64]) -> Vec<f64> {
        let n = self.pts.len();
        let mut rhs: Vec<f64> = self.centres.iter().zip(&self.model_v).map(|(&v, m)| u[v] - m).collect();
        rhs.extend([0.0; 4]);
        let coef = self.lu.solve(&rhs);
        self.dirs
            .iter()
            .zip(&self.model_dirs)
            .map(|(&x, md)| {
                let mut s = 0.0;
                for j in 0..n {
                    s += coef[j] * dist(x, self.pts[j]).powi(3);
                }
                let q = poly(x, self.shift, self.scale);
                for c in 0..4 {
                    s += coef[n + c] * q[c];
                }
                s + md
            })
            .collect()
    }

    pub fn shape(&self) -> (usize, usize) {
        let _ = (&self.cu, &self.betas);
        (self.n_lat, self.n_lon)
    }
}

fn poly(x: V3, shift: V3, scale: V3) -> [f64; 4] {
    [1.0, (x[0] - shift[0]) / scale[0], (x[1] - shift[1]) / scale[1], (x[2] - shift[2]) / scale[2]]
}
