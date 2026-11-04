//! Fast Uniformizer: solve directly for the uniformized cone metric (K = rho on smooth vertices,
//! prescribed cone defects, rho an OUTPUT via Gauss-Bonnet) without streaming the flow, then
//! realize that one metric in 3D.
//!
//! The intrinsic solve is pseudo-transient continuation on the flow equation: implicit (backward
//! Euler) steps of the cone-preserving discrete Ricci flow with geometrically growing dt, started
//! from the flat-cone far field `v_asymp` (`DiscreteConeFlow::smart_u`).  Each step is the exact
//! damped Newton of `flow.rs` (sparse Jacobian + Sherman-Morrison rank-one rho term), so the
//! stationary Newton is reached without ever taking an explicit step and without fixing rho (the
//! mass term regularises the dilation / constant-shift null direction).  As dt -> infinity the
//! step equation becomes the stationary equation rho A_v = K_v; we stop once that residual is at
//! machine level.  Gauss-Bonnet (sum of defects = 4 pi) is exact by construction of the defects
//! and is reported as the certificate together with the cone-defect error and max |K/rho - 1|.
//!
//! Admissibility (port of `core/model.py::check_admissibility`, sphere, all betas in (0, 1)):
//! one cone: impossible; two cones: Troyanov, equal betas required (football); kappa = sum beta -
//! (n - 2) = 0: flat (no curvature to uniformize to here); kappa < 0: always admissible
//! (McOwen/Troyanov); kappa > 0, n >= 3: Luo-Tian, kappa < 2 min beta.
use crate::flow::DiscreteConeFlow;

/// Statistics of one uniformization solve.
#[derive(Clone, Debug, Default)]
pub struct UnifStats {
    /// accepted implicit continuation steps
    pub steps: usize,
    /// rejected (Newton failed) attempts
    pub rejected: usize,
    /// total Newton iterations
    pub newton: usize,
    /// final stationary residual max_v |rho A_v - K_v| over smooth vertices
    pub resid: f64,
    /// final continuation dt
    pub dt: f64,
}

/// `check_admissibility` for betas in (0, 1) on the sphere.  Err(message) when no
/// constant-curvature cone metric exists (or kappa = 0: flat, refused like True Ricci Flow).
pub fn admissibility(betas: &[f64]) -> Result<(), String> {
    let n = betas.len();
    let kappa = betas.iter().sum::<f64>() - (n as f64 - 2.0);
    if n == 0 {
        return Ok(());
    }
    if n == 1 {
        return Err("inadmissible: a single cone point on the sphere is impossible (any curvature sign)".into());
    }
    if n == 2 {
        if (betas[0] - betas[1]).abs() > 1e-9 {
            return Err(format!(
                "inadmissible: two cones need equal angles (Troyanov 1989: a football); got beta = {:.3} and {:.3}",
                betas[0], betas[1]
            ));
        }
        return Ok(());
    }
    if kappa.abs() <= 1e-12 {
        return Err(crate::pipeline::NEG_REFUSAL.into());
    }
    if kappa < 0.0 {
        return Ok(()); // McOwen / Troyanov: a unique hyperbolic cone metric
    }
    let bmin = betas.iter().cloned().fold(f64::INFINITY, f64::min);
    let margin = 2.0 * bmin - kappa;
    if margin > 1e-9 {
        Ok(())
    } else {
        Err(format!(
            "inadmissible: no spherical metric with these angles (Luo-Tian: need sum beta - (n-2) < 2 min beta, here {:.3} >= {:.3})",
            kappa,
            2.0 * bmin
        ))
    }
}

fn stationary_resid(flow: &DiscreteConeFlow, u: &[f64]) -> f64 {
    let (_, g) = flow.geom(u);
    let a_s: f64 = g.av.iter().zip(&flow.smooth).filter(|(_, &s)| s).map(|(a, _)| a).sum();
    let rho = flow.c / a_s;
    let mut r = 0.0f64;
    for v in 0..u.len() {
        if flow.smooth[v] {
            r = r.max((rho * g.av[v] - g.k[v]).abs());
        }
    }
    for (ci, &c) in flow.mesh.cones.iter().enumerate() {
        r = r.max((g.k[c] - flow.cone_k[ci]).abs());
    }
    r
}

/// Uniformized log-scale u (area gauge 4 pi) from the far field `v_asymp`.
pub fn uniformize(flow: &mut DiscreteConeFlow) -> Result<(Vec<f64>, UnifStats), String> {
    let mut u = flow.smart_u().map_err(|e| format!("far-field start failed: {e}"))?;
    let mut st = UnifStats::default();
    let mut dt = 0.05f64;
    const DT_MAX: f64 = 1e10;
    const TOL: f64 = 1e-10;
    let newton0 = flow.stats.newton;
    for _ in 0..400 {
        let (_, g) = flow.geom(&u);
        // Mn (u - u_prev)/dt - (rho A - K) = 0: the residual is O(1) at every dt, so the absolute
        // Newton tolerance stays meaningful as dt -> infinity.
        let mn: Vec<f64> = g.av.iter().map(|a| a / dt).collect();
        let (un, ok, _) = flow.newton(&u, &u, &mn, 1.0, 1.0);
        if !ok {
            st.rejected += 1;
            dt /= 4.0;
            if dt < 1e-9 {
                return Err("uniformization failed: implicit continuation step did not converge".into());
            }
            continue;
        }
        u = flow.gauge(&un);
        st.steps += 1;
        let r = stationary_resid(flow, &u);
        st.resid = r;
        st.dt = dt;
        if r < TOL && dt >= 1e3 {
            st.newton = flow.stats.newton - newton0;
            return Ok((u, st));
        }
        dt = (dt * 8.0).min(DT_MAX);
    }
    st.newton = flow.stats.newton - newton0;
    if st.resid < 1e-7 {
        return Ok((u, st));
    }
    Err(format!("uniformization did not converge (stationary residual {:.1e})", st.resid))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn admissibility_rules() {
        assert!(admissibility(&[0.6, 0.6]).is_ok());
        assert!(admissibility(&[0.6, 0.7]).is_err());
        assert!(admissibility(&[0.9, 0.9, 0.2]).is_err()); // kappa 1.0 > 0.4
        assert!(admissibility(&[0.6, 0.7, 0.5]).is_ok());
        assert!(admissibility(&[0.5, 0.5, 0.5, 0.5]).is_err()); // flat
        assert!(admissibility(&[0.3, 0.3, 0.3, 0.3, 0.3]).is_ok()); // rho < 0
    }
}
