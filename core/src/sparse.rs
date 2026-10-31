//! Minimal sparse assembly (COO, duplicates summed) + sparse LU via `faer`
//! (sequential; compiles to wasm32-unknown-unknown without threads).
use faer::prelude::*;
use faer::sparse::{SparseColMat, Triplet};

/// COO triplet accumulator for an n x n matrix; duplicates are summed.
#[derive(Clone, Debug, Default)]
pub struct Coo {
    pub n: usize,
    pub rows: Vec<usize>,
    pub cols: Vec<usize>,
    pub vals: Vec<f64>,
}

impl Coo {
    pub fn new(n: usize) -> Self {
        Coo { n, rows: Vec::new(), cols: Vec::new(), vals: Vec::new() }
    }
    #[inline]
    pub fn push(&mut self, i: usize, j: usize, v: f64) {
        self.rows.push(i);
        self.cols.push(j);
        self.vals.push(v);
    }
    /// Compressed (sorted by (col,row), duplicates summed) entries.
    pub fn compress(&self) -> Vec<(usize, usize, f64)> {
        let mut idx: Vec<usize> = (0..self.vals.len()).collect();
        idx.sort_by_key(|&k| (self.cols[k], self.rows[k]));
        let mut out: Vec<(usize, usize, f64)> = Vec::with_capacity(idx.len());
        for k in idx {
            let (i, j, v) = (self.rows[k], self.cols[k], self.vals[k]);
            match out.last_mut() {
                Some(last) if last.0 == i && last.1 == j => last.2 += v,
                _ => out.push((i, j, v)),
            }
        }
        out
    }
    /// y = A x
    pub fn matvec(&self, x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0; self.n];
        for k in 0..self.vals.len() {
            y[self.rows[k]] += self.vals[k] * x[self.cols[k]];
        }
        y
    }
    /// y = A^T x
    pub fn matvec_t(&self, x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0; self.n];
        for k in 0..self.vals.len() {
            y[self.cols[k]] += self.vals[k] * x[self.rows[k]];
        }
        y
    }
}

/// Factorized sparse matrix (LU with partial pivoting).
pub struct SparseLu {
    n: usize,
    lu: faer::sparse::linalg::solvers::Lu<usize, f64>,
}

impl SparseLu {
    pub fn factor(a: &Coo) -> Result<Self, String> {
        let trip: Vec<Triplet<usize, usize, f64>> =
            a.compress().into_iter().map(|(i, j, v)| Triplet::new(i, j, v)).collect();
        let m = SparseColMat::<usize, f64>::try_new_from_triplets(a.n, a.n, &trip)
            .map_err(|e| format!("sparse assembly failed: {e:?}"))?;
        let lu = m.sp_lu().map_err(|e| format!("sparse LU failed: {e:?}"))?;
        Ok(SparseLu { n: a.n, lu })
    }
    /// Solve A X = B for several right-hand sides at once.
    pub fn solve_many(&self, rhs: &[&[f64]]) -> Vec<Vec<f64>> {
        let mut b = Mat::<f64>::zeros(self.n, rhs.len());
        for (c, r) in rhs.iter().enumerate() {
            for i in 0..self.n {
                b[(i, c)] = r[i];
            }
        }
        self.lu.solve_in_place(b.as_mut());
        (0..rhs.len()).map(|c| (0..self.n).map(|i| b[(i, c)]).collect()).collect()
    }
    pub fn solve(&self, rhs: &[f64]) -> Vec<f64> {
        self.solve_many(&[rhs]).pop().unwrap()
    }
}
