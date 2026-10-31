//! End-to-end: from fixture INPUTS only (positions, betas, seed, level) run the full Rust pipeline
//! (`pipeline::Session`) and compare to the rho > 0 golden fixtures.
#[path = "../src/fixtures.rs"]
mod fixtures;

use fixtures::Fixture;
use ricci_core::flow::Seed;
use ricci_core::pipeline::{Session, SessionConfig};
use ricci_core::realize::geom;
use ricci_core::seedmesh::page_level;
use std::time::Instant;

fn nrel(a: &[f64], b: &[f64]) -> f64 {
    let num = a.iter().zip(b).fold(0.0f64, |m, (x, y)| m.max((x - y).abs()));
    num / b.iter().fold(0.0f64, |m, y| m.max(y.abs())).max(1e-300)
}
fn rms(v: &[f64]) -> f64 {
    (v.iter().map(|x| x * x).sum::<f64>() / v.len() as f64).sqrt()
}

pub fn config_of(fx: &Fixture) -> (SessionConfig, serde_json::Value) {
    let js: serde_json::Value = serde_json::from_str(&fx.json()).unwrap();
    let req = js["meta"]["request"].clone();
    let positions = req["positions"].as_array().unwrap().iter().map(|p| [p[0].as_f64().unwrap(), p[1].as_f64().unwrap()]).collect();
    let betas = req["betas"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
    let seed = Seed::parse(req["seed_surface"].as_str().unwrap()).unwrap();
    let level = page_level(req["resolution"].as_str().unwrap()).unwrap();
    (SessionConfig { positions, betas, seed, level }, js)
}

fn check(name: &str) -> Vec<String> {
    let fx = Fixture::open(name);
    let (cfg, js) = config_of(&fx);
    let mut bad = vec![];
    let t0 = Instant::now();
    let mut s = Session::new(cfg).unwrap();
    let rows = s.run_all().unwrap();
    let wall = t0.elapsed().as_secs_f64();
    let (nv, ne) = (s.mesh.nv(), s.mesh.e.len());
    let nf = rows.len();
    let (u, l, k, rho) = (fx.f64("u").data, fx.f64("L").data, fx.f64("K_defect").data, fx.f64("rho_frame").data);
    let xf = fx.f64("vertices3d").data;
    let cf = fx.f64("cf_w").data;
    let nd_fix = fx.i64("n_deep").data;
    let ng = 64 * 128;
    let (mut eu, mut el, mut ek, mut er, mut gb, mut ecf, mut worst_len) = (0f64, 0f64, 0f64, 0f64, 0f64, 0f64, 0f64);
    eprintln!("== {name}: nv={nv} frames={nf} wall={wall:.2}s");
    for (i, r) in rows.iter().enumerate() {
        eu = eu.max(nrel(&r.u, &u[i * nv..(i + 1) * nv]));
        el = el.max(nrel(&r.lengths, &l[i * ne..(i + 1) * ne]));
        ek = ek.max(nrel(&r.k_defect, &k[i * nv..(i + 1) * nv]));
        er = er.max(((r.rho - rho[i]) / rho[i]).abs());
        gb = gb.max(r.gb_resid.abs());
        ecf = ecf.max(r.cf_w.iter().zip(&cf[i * ng..(i + 1) * ng]).fold(0f64, |m, (a, b)| m.max((a - b).abs())));
        let xfix: Vec<[f64; 3]> = xf[i * nv * 3..(i + 1) * nv * 3].chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let lk = &l[i * ne..(i + 1) * ne];
        let (rr, rf) = (rms(&geom::len_rel(&r.vertices, &s.mesh.e, lk)), rms(&geom::len_rel(&xfix, &s.mesh.e, lk)));
        let dev = (rr - rf).abs() / rf;
        worst_len = worst_len.max(dev);
        let fold_ok = r.n_deep as i64 == nd_fix[i];
        eprintln!("  {:2} t={:8.4} rms rust {:.4e} fix {:.4e} (dev {:5.2}%) deep {}/{}{} start={}{}", i, r.t, rr, rf, 100.0 * dev, r.n_deep, nd_fix[i],
            if fold_ok { "" } else { " !" }, r.start, if r.polytope { " polytope" } else { "" });
        if dev > 0.10 { bad.push(format!("{name} frame {i}: len rms dev {:.1}%", 100.0 * dev)); }
        if !fold_ok { bad.push(format!("{name} frame {i}: deep folds {} vs {}", r.n_deep, nd_fix[i])); }
    }
    let q = s.quality.clone().unwrap();
    let label_fix = js["scalars"]["quality"]["label"].as_str().unwrap();
    let st = &js["meta"]["solver"]["stats_after_run"];
    let (sf, stot) = (s.stats_flow.clone().unwrap(), s.stats_total().unwrap());
    eprintln!("   rel u {eu:.2e} L {el:.2e} K {ek:.2e} rho {er:.2e} | GB {gb:.2e} | cf_w max abs {ecf:.2e} | worst len-rms dev {:.2}%", 100.0 * worst_len);
    eprintln!("   label rust={:?} fixture={:?} | stats frames steps/newton {}/{} total {}/{} fixture {}/{} | quality detail: {}",
        q.label, label_fix, sf.steps, sf.newton, stot.steps, stot.newton, st["steps"], st["newton"], q.detail);
    for (lab, v) in [("u", eu), ("L", el), ("K", ek), ("rho", er)] {
        if v > 1e-8 { bad.push(format!("{name}: rel {lab} {v:.2e}")); }
    }
    if gb > 1e-12 { bad.push(format!("{name}: GB {gb:.2e}")); }
    if ecf > 1e-3 { bad.push(format!("{name}: cf_w {ecf:.2e}")); }
    if q.label != label_fix { bad.push(format!("{name}: label {}", q.label)); }
    if sf.steps as u64 != st["steps"].as_u64().unwrap() || sf.newton as u64 != st["newton"].as_u64().unwrap()
        || sf.rejected as u64 != st["rejected"].as_u64().unwrap() || sf.newton_fail as u64 != st["newton_fail"].as_u64().unwrap() {
        bad.push(format!("{name}: stats {sf:?}"));
    }
    bad
}

const POS: [&str; 5] = ["n3_pos_gluing", "n4_pos_gluing", "n4_pos_polytope", "n4_pos_smart", "n5_pos_gluing"];

#[test]
fn pipeline_rho_pos_fixtures() {
    let only = std::env::var("PIPELINE_FIXTURE").ok();
    let mut bad = vec![];
    for name in POS {
        if only.as_deref().map_or(true, |o| o == name) {
            bad.extend(check(name));
        }
    }
    assert!(bad.is_empty(), "{bad:#?}");
}

#[test]
fn rho_neg_is_refused() {
    let fx = Fixture::open("n4_neg_gluing");
    let (cfg, _) = config_of(&fx);
    let e = Session::new(cfg).err().expect("rho<0 must be refused");
    assert!(e.contains("rho > 0"), "{e}");
}

#[test]
fn overlapping_caps_and_bad_beta() {
    // two close cones: caps get clamped, the session still runs to a quality label
    let cfg = SessionConfig { positions: vec![[0.1, 0.0], [0.25, 0.0], [-0.8, 0.5]], betas: vec![0.5, 0.6, 0.7], seed: Seed::Smart, level: page_level("smoke").unwrap() };
    let mut s = Session::new(cfg).unwrap();
    assert!(!s.flow_mesh.seed.overlaps.is_empty());
    let rows = s.run_all().unwrap();
    assert_eq!(rows.len(), 21);
    assert!(rows.iter().all(|r| r.gb_resid.abs() < 1e-12));
    eprintln!("overlap config: quality {:?}", s.quality.as_ref().unwrap().label);
    let bad = SessionConfig { positions: vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], betas: vec![1.2, 0.6, 0.7], seed: Seed::Smart, level: page_level("smoke").unwrap() };
    assert!(Session::new(bad).is_err());
}
