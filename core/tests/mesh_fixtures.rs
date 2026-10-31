//! Mesh construction from fixture INPUTS (positions, betas, level) reproduces the fixture mesh.
#[path = "../src/fixtures.rs"]
mod fixtures;

use fixtures::Fixture;
use ricci_core::seedmesh::{build_mesh, page_level};

pub const ALL: [&str; 6] = ["n3_pos_gluing", "n4_pos_gluing", "n4_pos_polytope", "n4_pos_smart", "n5_pos_gluing", "n4_neg_gluing"];

fn check(name: &str) {
    let fx = Fixture::open(name);
    let js: serde_json::Value = serde_json::from_str(&fx.json()).unwrap();
    let req = &js["meta"]["request"];
    let pos: Vec<[f64; 2]> = req["positions"].as_array().unwrap().iter().map(|p| [p[0].as_f64().unwrap(), p[1].as_f64().unwrap()]).collect();
    let betas: Vec<f64> = req["betas"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
    let level = page_level(req["resolution"].as_str().unwrap()).unwrap();
    let (mesh, _) = build_mesh(&pos, &betas, level).unwrap();
    let y = fx.f64("sphere_vertices").data;
    let f = fx.i64("sphere_faces").data;
    let cones: Vec<usize> = fx.i64("cone_ids").data.iter().map(|&c| c as usize).collect();
    assert_eq!(mesh.y.len() * 3, y.len(), "{name}: nv");
    let dy = mesh.y.iter().flatten().zip(&y).fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
    let ff: Vec<i64> = mesh.f.iter().flatten().map(|&v| v as i64).collect();
    let same_order = ff == f;
    let mut a: Vec<[i64; 3]> = ff.chunks(3).map(|c| { let mut t = [c[0], c[1], c[2]]; let r = t.iter().position(|v| v == t.iter().min().unwrap()).unwrap(); t.rotate_left(r); t }).collect();
    let mut b: Vec<[i64; 3]> = f.chunks(3).map(|c| { let mut t = [c[0], c[1], c[2]]; let r = t.iter().position(|v| v == t.iter().min().unwrap()).unwrap(); t.rotate_left(r); t }).collect();
    a.sort(); b.sort();
    eprintln!("{name}: nv={} nf={} max|dY|={dy:.2e} faces identical order={same_order} same set={}", mesh.y.len(), mesh.f.len(), a == b);
    assert_eq!(mesh.cones, cones);
    assert!(dy <= 1e-12, "{name}: vertex error {dy}");
    assert!(same_order, "{name}: faces differ");
}

#[test]
fn mesh_from_inputs() {
    for n in ALL { check(n); }
}
