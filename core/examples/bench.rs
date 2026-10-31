//! Native wall time of the full Session per rho>0 fixture: `cargo run --release --example bench`
use ricci_core::flow::Seed;
use ricci_core::pipeline::{Phase, Session, SessionConfig};
use ricci_core::seedmesh::page_level;
use std::time::Instant;

fn main() {
    let fx = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("fixtures");
    for name in ["n3_pos_gluing", "n4_pos_gluing", "n4_pos_polytope", "n4_pos_smart", "n5_pos_gluing"] {
        let js: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(fx.join(name).join("fixture.json")).unwrap()).unwrap();
        let r = &js["meta"]["request"];
        let cfg = SessionConfig {
            positions: r["positions"].as_array().unwrap().iter().map(|p| [p[0].as_f64().unwrap(), p[1].as_f64().unwrap()]).collect(),
            betas: r["betas"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect(),
            seed: Seed::parse(r["seed_surface"].as_str().unwrap()).unwrap(),
            level: page_level(&std::env::var("BENCH_LEVEL").unwrap_or(r["resolution"].as_str().unwrap().to_string())).unwrap(),
        };
        let t0 = Instant::now();
        let mut s = Session::new(cfg).unwrap();
        let mut per: Vec<(Phase, f64)> = vec![];
        let mut max_chunk = 0f64;
        while !s.is_done() {
            let ph = s.phase();
            let c = Instant::now();
            s.next_frame().unwrap();
            let dt = c.elapsed().as_secs_f64();
            max_chunk = max_chunk.max(dt);
            match per.last_mut() {
                Some((p, t)) if *p == ph => *t += dt,
                _ => per.push((ph, dt)),
            }
        }
        let ph: Vec<String> = per.iter().map(|(p, t)| format!("{p:?} {t:.2}")).collect();
        println!("{name}: wall {:.2} s (max chunk {:.0} ms) [{}] python fixture {:.2} s", t0.elapsed().as_secs_f64(), 1e3 * max_chunk, ph.join(", "),
            js["meta"]["wall_s_total"].as_f64().unwrap_or(f64::NAN));
    }
}
