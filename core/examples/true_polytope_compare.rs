//! True Polytope vs Polytope vs Cone Gluing on the lane's test configs.
//! `cargo run --release --example true_polytope_compare -- <out_dir>`: prints a table, writes
//! `<out_dir>/<config>_<seed>.json` (per-frame metrics + vertices + faces) for the renderer.
use ricci_core::flow::Seed;
use ricci_core::pipeline::{Session, SessionConfig};
use ricci_core::realize::geom::{self, V3};
use ricci_core::seedmesh::page_level;
use serde_json::json;
use std::time::Instant;

fn ring(n: usize, r: f64, ph: f64) -> Vec<[f64; 2]> {
    (0..n).map(|k| {
        let a = ph + 2.0 * std::f64::consts::PI * k as f64 / n as f64;
        [r * a.cos(), r * a.sin()]
    }).collect()
}

fn bbox_diag(x: &[V3]) -> f64 {
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for p in x {
        for c in 0..3 {
            lo[c] = lo[c].min(p[c]);
            hi[c] = hi[c].max(p[c]);
        }
    }
    geom::norm(geom::sub(hi, lo))
}

fn main() {
    let out = std::path::PathBuf::from(std::env::args().nth(1).unwrap_or("artifacts/true-polytope".into()));
    std::fs::create_dir_all(&out).unwrap();
    let d4 = vec![[0.7, 0.0], [-0.7, 0.0], [0.0, 0.7], [0.0, -0.7]];
    let configs: Vec<(&str, Vec<[f64; 2]>, Vec<f64>, &str)> = vec![
        ("n4_default_coarse", d4.clone(), vec![0.6, 0.7, 0.5, 0.8], "coarse"),
        ("n4_default_medium", d4.clone(), vec![0.6, 0.7, 0.5, 0.8], "medium"),
        ("n3_coarse", ring(3, 0.7, 0.5 * std::f64::consts::PI), vec![0.5, 0.6, 0.7], "coarse"),
        ("n5_pentagon_coarse", ring(5, 0.7, 0.5 * std::f64::consts::PI), vec![0.6, 0.7, 0.5, 0.8, 0.7], "coarse"),
        ("n3_sharp_coarse", ring(3, 0.7, 0.5 * std::f64::consts::PI), vec![0.3, 0.35, 0.4], "coarse"),
        ("n3_closepair_medium", vec![[0.06, 0.0], [-0.06, 0.0], [0.0, 0.7]], vec![0.8, 0.8, 0.8], "medium"),
    ];
    let only = std::env::var("ONLY").ok();
    println!("config | seed | nv | wall s | frame0 len_rms | frame0 reflex/deep | max len_rms (k>0) | max deep (k>0) | max jump | max gb | max cone err | crease exact | sharp0 | final label");
    for (name, pos, betas, lvl) in &configs {
        if let Some(o) = &only {
            if !name.contains(o.as_str()) {
                continue;
            }
        }
        for seed in ["true_polytope", "polytope", "gluing"] {
            let cfg = SessionConfig { positions: pos.clone(), betas: betas.clone(), seed: Seed::parse(seed).unwrap(), level: page_level(lvl).unwrap() };
            let t0 = Instant::now();
            let mut s = match Session::new(cfg) {
                Ok(s) => s,
                Err(e) => {
                    println!("{name} | {seed} | FAILED: {e}");
                    continue;
                }
            };
            let rows = match s.run_all() {
                Ok(r) => r,
                Err(e) => {
                    println!("{name} | {seed} | FAILED: {e}");
                    continue;
                }
            };
            let wall = t0.elapsed().as_secs_f64();
            let mut jumps = vec![0.0];
            for k in 1..rows.len() {
                let a = geom::kabsch_align(&rows[k].vertices, &rows[k - 1].vertices);
                let m = a.iter().zip(&rows[k - 1].vertices).map(|(p, q)| geom::norm(geom::sub(*p, *q))).fold(0.0, f64::max);
                jumps.push(m / bbox_diag(&rows[k].vertices));
            }
            let f0 = rows[0].faces.clone().unwrap_or(s.mesh.f.clone());
            let dh0 = geom::dihedrals(&rows[0].vertices, &f0);
            let sharp0 = dh0.iter().filter(|d| d.abs() > 0.1).count();
            let crease_exact = match &s.true_poly {
                Some(tp) => {
                    let te = geom::face_edge_table(&s.mesh.f);
                    let mut folded: Vec<[usize; 2]> = te.iter().zip(&dh0).filter(|(_, d)| d.abs() > 1e-7).map(|(e, _)| [e.0[0].min(e.0[1]), e.0[0].max(e.0[1])]).collect();
                    folded.sort();
                    let mut cr: Vec<[usize; 2]> = tp.crease_edges.iter().map(|&k| s.mesh.e[k]).collect();
                    cr.sort();
                    format!("{} ({} crease edges)", folded == cr, cr.len())
                }
                None => "-".into(),
            };
            let mx = |f: &dyn Fn(&ricci_core::pipeline::FrameOut) -> f64, from: usize| rows[from..].iter().map(|r| f(r)).fold(0.0, f64::max);
            let q = s.quality.clone().unwrap();
            println!(
                "{name} | {seed} | {} | {wall:.2} | {:.2e} | {}/{} | {:.2e} | {} | {:.3} | {:.1e} | {:.1e} | {crease_exact} | {sharp0} | {}",
                s.mesh.nv(),
                rows[0].len_rel_rms,
                rows[0].n_reflex,
                rows[0].n_deep,
                mx(&|r| r.len_rel_rms, 1),
                mx(&|r| r.n_deep as f64, 1) as usize,
                jumps.iter().cloned().fold(0.0, f64::max),
                mx(&|r| r.gb_resid.abs(), 0),
                mx(&|r| r.cone_defect_err.abs(), 0),
                q.label
            );
            let fr: Vec<_> = rows
                .iter()
                .zip(&jumps)
                .map(|(r, j)| {
                    json!({"k": r.k, "t": r.t, "len_rel_rms": r.len_rel_rms, "len_rel_max": r.len_rel_max, "n_reflex": r.n_reflex,
                        "n_deep": r.n_deep, "jump": j, "gb": r.gb_resid, "cone_err": r.cone_defect_err, "seed_folds": r.seed_folds,
                        "vertices": r.vertices, "faces": r.faces})
                })
                .collect();
            let crease: Vec<[usize; 2]> = s.true_poly.as_ref().map(|tp| tp.crease_edges.iter().map(|&k| s.mesh.e[k]).collect()).unwrap_or_default();
            let js = json!({"config": name, "seed": seed, "wall": wall, "label": q.label, "detail": q.detail, "faces": s.mesh.f,
                "cones": s.mesh.cones, "corners": s.seed_corners, "crease": crease, "frames": fr});
            std::fs::write(out.join(format!("{name}_{seed}.json")), js.to_string()).unwrap();
        }
    }
}
