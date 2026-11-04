//! Probe: strategies for realizing the final (uniformized) frame when its cold fit folds.
//! `cargo run --release --example final_probe -- SEED LEVEL "betas" "positions"`
use ricci_core::flow::{DiscreteConeFlow, FlowOptions, Seed};
use ricci_core::realize::{self, FastRealizer, RealizeMesh, Realized};
use ricci_core::seedmesh::{self, page_level};
use std::time::Instant;

fn rep(name: &str, r: &Realized, f: &[[usize; 3]], dt: f64) {
    println!("{name:28} deep {:3} rms {:.2e} max {:.2e} ok {} [{:.2}s] {}", realize::n_deep_folds(&r.x, f, 0.15), r.len_rel_rms, r.len_rel_max, r.ok, dt, r.start);
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let seed = Seed::parse(&a[1]).unwrap();
    let level = page_level(&a[2]).unwrap();
    let betas: Vec<f64> = a[3].split(',').map(|v| v.parse().unwrap()).collect();
    let pos: Vec<[f64; 2]> = a[4].split(';').map(|p| { let v: Vec<f64> = p.split(',').map(|z| z.parse().unwrap()).collect(); [v[0], v[1]] }).collect();
    let (mesh, _) = seedmesh::build_mesh(&pos, &betas, level).unwrap();
    let flow = DiscreteConeFlow::new(mesh.clone(), betas.clone(), FlowOptions { seed, ..Default::default() });
    let (frames, _fl) = flow.frames().unwrap();
    let rm = RealizeMesh::new(mesh.f.clone(), mesh.e.clone(), mesh.cones.clone(), mesh.nv());
    let (x0, kind) = realize::fast_seed(&mesh.y, &mesh.f, &mesh.cones, &betas);
    println!("seed kind {:?} nv {}", kind, mesh.nv());
    let rlz = FastRealizer::new(rm, x0.clone(), kind);
    let lt = &frames.last().unwrap().l;
    let c = Instant::now();
    let cold = rlz.realize(lt, None);
    rep("cold", &cold, &mesh.f, c.elapsed().as_secs_f64());
    for s in [0.9, 0.8, 1.1, 0.7] {
        let xs: Vec<[f64; 3]> = x0.iter().map(|p| [p[0] * s, p[1] * s, p[2] * s]).collect();
        let c = Instant::now();
        let r = rlz.realize_cold_from(lt, &xs);
        rep(&format!("cold seed x{s}"), &r, &mesh.f, c.elapsed().as_secs_f64());
    }
    // bare sphere seed
    let c = Instant::now();
    let r = rlz.realize_cold_from(lt, &mesh.y);
    rep("cold bare sphere", &r, &mesh.f, c.elapsed().as_secs_f64());
    // forward chain from frame 0 cold through every frame
    let c = Instant::now();
    let mut cur = rlz.realize(&frames[0].l, None);
    for fr in &frames[1..] {
        let r = rlz.realize(&fr.l, Some(&cur));
        if r.ok {
            cur = r;
        }
    }
    rep("forward chain all frames", &cur, &mesh.f, c.elapsed().as_secs_f64());
    // chain on the late frames only: from frame m cold
    for m in [14usize, 17] {
        let c = Instant::now();
        let mut cur = rlz.realize(&frames[m].l, None);
        let d0 = realize::n_deep_folds(&cur.x, &mesh.f, 0.15);
        for fr in &frames[m + 1..] {
            let r = rlz.realize(&fr.l, Some(&cur));
            if r.ok {
                cur = r;
            }
        }
        rep(&format!("chain from {m} (cold deep {d0})"), &cur, &mesh.f, c.elapsed().as_secs_f64());
    }
}
