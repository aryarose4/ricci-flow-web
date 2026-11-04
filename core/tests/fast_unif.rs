//! Fast Uniformizer: direct implicit solve from the far field + one 3D realization.
//! Certificate: Gauss-Bonnet (sum of defects = 4 pi), cone defects, K = rho pointwise.
use ricci_core::flow::Seed;
use ricci_core::pipeline::{Session, SessionConfig};
use ricci_core::seedmesh::page_level;
use std::f64::consts::PI;
use std::time::Instant;

fn ngon(n: usize) -> Vec<[f64; 2]> {
    (0..n).map(|i| { let a = 2.0 * PI * i as f64 / n as f64; [0.7 * a.cos(), 0.7 * a.sin()] }).collect()
}

fn run(betas: &[f64], level: &str) -> (Session, f64) {
    let cfg = SessionConfig { positions: ngon(betas.len()), betas: betas.to_vec(), seed: Seed::Gluing, level: page_level(level).unwrap(), anti: Default::default() };
    let t0 = Instant::now();
    let mut s = Session::new(cfg).unwrap();
    s.set_fast();
    let rows = s.run_all().unwrap();
    assert_eq!(rows.len(), 1);
    (s, t0.elapsed().as_secs_f64())
}

fn certify(betas: &[f64], level: &str) {
    let (s, wall) = run(betas, level);
    let r = &s.rows()[0];
    let st = s.unif.clone().unwrap();
    let q = s.quality.clone().unwrap();
    println!(
        "FAST n={} betas={:?} {}: nv={} rho={:.5} gb={:.1e} cone={:.1e} |K/rho-1|max={:.1e} steps={} newton={} resid={:.1e} len_max={:.2e} rms={:.2e} deep={} label='{}' wall={:.2}s",
        betas.len(), betas, level, s.mesh.nv(), r.rho, r.gb_resid, r.cone_defect_err, r.k_dev_max, st.steps, st.newton, st.resid,
        r.len_rel_max, r.len_rel_rms, r.n_deep, q.label, wall
    );
    assert!(r.gb_resid.abs() < 1e-9, "Gauss-Bonnet");
    assert!(r.cone_defect_err < 1e-9, "cone defects");
    assert!(r.k_dev_max < 1e-8, "K = rho");
    let kappa = betas.iter().sum::<f64>() - (betas.len() as f64 - 2.0);
    assert!(r.rho.signum() == kappa.signum(), "sign(rho) = sign(kappa)");
    assert!(r.vertices.iter().all(|p| p.iter().all(|c| c.is_finite())));
}

#[test]
fn fast_gauss_bonnet_positive() {
    certify(&[0.6, 0.7, 0.5, 0.8], "coarse");
    certify(&[0.6, 0.7, 0.8], "coarse");
}

#[test]
fn fast_gauss_bonnet_n2_n6_n8() {
    certify(&[0.6, 0.6], "coarse");
    certify(&[0.9, 0.9, 0.9, 0.9, 0.9, 0.9], "coarse");
    certify(&[0.95; 8], "coarse");
}

#[test]
fn fast_gauss_bonnet_negative() {
    certify(&[0.3, 0.4, 0.5, 0.3, 0.4], "coarse");
    certify(&[0.6, 0.7, 0.5, 0.8, 0.7, 0.6, 0.5, 0.7], "coarse");
}

#[test]
fn inadmissible_refused() {
    let mk = |b: &[f64]| SessionConfig { positions: ngon(b.len()), betas: b.to_vec(), seed: Seed::Gluing, level: page_level("coarse").unwrap(), anti: Default::default() };
    assert!(Session::new(mk(&[0.6, 0.7])).is_err());
    assert!(Session::new(mk(&[0.9, 0.9, 0.2])).is_err());
}
