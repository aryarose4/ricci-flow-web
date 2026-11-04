//! rho < 0 True Ricci Flow (forward continuation + tracker) against Python `run_flow_true_negative`
//! goldens (`fixtures/*_neg_*_coarse`, exporter
//! `curvature-flow-and-geodesics/singular_ricci/dev/scratch/export_neg_flow_fixtures.py`).
//!
//! Tolerances:
//! * intrinsic (t, rho): rel 1e-8 (same BDF2 step control; Newton tol 1e-11; observed ~1e-15);
//!   K_defect: per-entry rel 1e-5 (Python streams `intr_defect` rounded to "%.6g");
//! * GB residual and cone defect error: magnitude < 1e-12 on both sides (exact by construction);
//! * embed_state sequence: exact;
//! * len_rel_rms / K_err_signed: loose.  The 3D shape comes from a nonconvex LM with fading bending
//!   (different linear algebra: faer vs scipy), so after a few warm-started frames the realized
//!   surfaces differ at the fit's own noise level, ~1e-4 relative edge error; on followed frames we
//!   require |rms_rust - rms_py| <= 0.5 * max(rms_py, 1e-3) + 2e-4 and |dK_signed| <= 0.03 (the
//!   tracker thresholds are 5e-3 and 0.10, so state decisions are not near the noise).  Lost frames
//!   keep different last-good shapes and are only reported.
#[path = "../src/fixtures.rs"]
mod fixtures;

use fixtures::Fixture;
use ricci_core::flow::Seed;
use ricci_core::pipeline::{Session, SessionConfig};
use ricci_core::seedmesh::page_level;
use std::time::Instant;

fn check(name: &str) -> Vec<String> {
    let fx = Fixture::open(name);
    let js: serde_json::Value = serde_json::from_str(&fx.json()).unwrap();
    let req = &js["meta"]["request"];
    let positions = req["positions"].as_array().unwrap().iter().map(|p| [p[0].as_f64().unwrap(), p[1].as_f64().unwrap()]).collect();
    let betas = req["betas"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
    let cfg = SessionConfig {
        positions,
        betas,
        seed: Seed::parse(req["seed_surface"].as_str().unwrap()).unwrap(),
        level: page_level(req["resolution"].as_str().unwrap()).unwrap(),
        anti: Default::default(),
    };
    let fr = js["scalars"]["frames"].as_array().unwrap();
    let kd = fx.f64("K_defect").data;
    let mut bad = vec![];
    let t0 = Instant::now();
    let mut s = Session::new(cfg).unwrap();
    assert!(s.neg);
    let mut rows = vec![];
    let mut tlast = Instant::now();
    let mut times = vec![];
    while !s.is_done() {
        if let Some(f) = s.next_frame().unwrap() {
            times.push(tlast.elapsed().as_secs_f64());
            tlast = Instant::now();
            rows.push(f);
        }
    }
    let wall = t0.elapsed().as_secs_f64();
    let nv = s.mesh.nv();
    if nv as u64 != js["meta"]["n_vertices"].as_u64().unwrap() {
        bad.push(format!("{name}: nv {nv}"));
    }
    if rows.len() != fr.len() {
        bad.push(format!("{name}: {} frames vs {}", rows.len(), fr.len()));
        return bad;
    }
    eprintln!("== {name}: nv={nv} frames={} wall={wall:.2}s (python {:.1}s)", rows.len(), js["meta"]["wall_s"].as_f64().unwrap());
    let (mut ek, mut et, mut er, mut drms, mut dks) = (0f64, 0f64, 0f64, 0f64, 0f64);
    for (i, (r, p)) in rows.iter().zip(fr).enumerate() {
        let pf = |k: &str| p[k].as_f64().unwrap_or(f64::NAN);
        et = et.max((r.t - pf("t")).abs());
        er = er.max(((r.rho - pf("rho")) / pf("rho")).abs());
        ek = ek.max(r.k_defect.iter().zip(&kd[i * nv..(i + 1) * nv]).fold(0f64, |m, (a, b)| m.max((a - b).abs() / b.abs().max(1e-12))));
        if r.gb_resid.abs() > 1e-12 || pf("gb_resid").abs() > 1e-12 {
            bad.push(format!("{name} frame {i}: GB {:.2e}", r.gb_resid));
        }
        if r.cone_defect_err > 1e-12 || pf("cone_defect_err") > 1e-12 {
            bad.push(format!("{name} frame {i}: cone defect err {:.2e}", r.cone_defect_err));
        }
        let ps = p["embed_state"].as_str().unwrap();
        if r.embed_state != ps {
            bad.push(format!("{name} frame {i}: state {} vs {ps}", r.embed_state));
        }
        if r.embed_lost != p["embed_lost"].as_bool().unwrap() {
            bad.push(format!("{name} frame {i}: embed_lost {}", r.embed_lost));
        }
        let (dr, dk) = ((r.len_rel_rms - pf("len_rel_rms")).abs(), (r.k_err_signed - pf("K_err_signed")).abs());
        let followed = ps == "following";
        if followed {
            drms = drms.max(dr);
            dks = dks.max(dk);
            if dr > 0.5 * pf("len_rel_rms").max(1e-3) + 2e-4 {
                bad.push(format!("{name} frame {i}: len_rel_rms {:.3e} vs {:.3e}", r.len_rel_rms, pf("len_rel_rms")));
            }
            if dk > 0.03 {
                bad.push(format!("{name} frame {i}: K_err_signed {:.4} vs {:.4}", r.k_err_signed, pf("K_err_signed")));
            }
        }
        // the page's overlays need these on every frame
        if r.lengths.len() != s.mesh.e.len() || r.a_vertex.len() != nv || r.cf_w.len() != 64 * 128 || r.vertices.len() != nv {
            bad.push(format!("{name} frame {i}: missing per-frame arrays"));
        }
        eprintln!(
            "  {i:2} t={:8.4} {:9} rms rust {:.3e} py {:.3e} | Kerr rust {:+.4} py {:+.4} | start {} | {:.3}s (py {:.3}s)",
            r.t, r.embed_state, r.len_rel_rms, pf("len_rel_rms"), r.k_err_signed, pf("K_err_signed"), r.start, times[i], pf("frame_s")
        );
    }
    eprintln!("   max |dt| {et:.1e} rel rho {er:.1e} rel Kdef {ek:.1e} | followed frames: max |d rms| {drms:.2e} |d Kerr| {dks:.3}");
    eprintln!("   notice rust: {:?}", s.embed_notice());
    eprintln!("   notice py  : {:?}", js["scalars"]["embed_notice"].as_str());
    eprintln!("   quality: {:?}", s.quality.as_ref().map(|q| &q.label));
    if et > 1e-12 || er > 1e-8 || ek > 1e-5 {
        bad.push(format!("{name}: intrinsic dt {et:.1e} rho {er:.1e} Kdef {ek:.1e}"));
    }
    let pn = js["scalars"]["embed_notice"].as_str();
    if s.embed_notice().is_some() != pn.is_some() {
        bad.push(format!("{name}: notice presence differs"));
    }
    if let (Some(a), Some(b)) = (s.embed_notice(), pn) {
        // same failure time (prefix up to "; it is")
        if a.split("; it is").next() != b.split("; it is").next() {
            bad.push(format!("{name}: notice {a:?}"));
        }
    }
    bad
}

#[test]
fn negative_fixtures() {
    let mut bad = vec![];
    for n in ["n4_neg_square_coarse", "n3_neg_small_coarse"] {
        bad.extend(check(n));
    }
    assert!(bad.is_empty(), "{bad:#?}");
}
