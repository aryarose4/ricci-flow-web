//! Golden-fixture tests for the intrinsic flow (mesh -> seed -> BDF2 frames).
#[path = "../src/fixtures.rs"]
mod fixtures;

use fixtures::Fixture;
use ricci_core::flow::{default_times, DiscreteConeFlow, FlowOptions, Seed};
use ricci_core::mesh::Mesh;

const REL_TOL: f64 = 1e-8;

fn rows3<T: Copy>(d: &[T]) -> Vec<[T; 3]> {
    d.chunks(3).map(|c| [c[0], c[1], c[2]]).collect()
}

/// normwise relative error max|a-b| / max|b|
fn nrel(a: &[f64], b: &[f64]) -> f64 {
    let num = a.iter().zip(b).fold(0.0f64, |m, (x, y)| m.max((x - y).abs()));
    let den = b.iter().fold(0.0f64, |m, y| m.max(y.abs()));
    num / den.max(1e-300)
}

fn seed_of(fx: &Fixture) -> Seed {
    let j = fx.json();
    for s in ["polytope", "smart", "gluing"] {
        if j.contains(&format!("\"seed_surface\": \"{s}\"")) {
            return Seed::parse(s).unwrap();
        }
    }
    panic!("seed_surface not found")
}

fn betas_of(fx: &Fixture) -> Vec<f64> {
    let j = fx.json();
    let a = j.find("\"betas\"").unwrap();
    let s = &j[a..];
    let (l, r) = (s.find('[').unwrap(), s.find(']').unwrap());
    s[l + 1..r].split(',').map(|x| x.trim().parse::<f64>().unwrap()).collect()
}

/// Returns (max rel u, max rel L, max rel K, max rel rho, max |GB|).
fn run_fixture(name: &str) -> (f64, f64, f64, f64, f64) {
    let fx = Fixture::open(name);
    let y = rows3(&fx.f64("sphere_vertices").data);
    let f: Vec<[usize; 3]> = rows3(&fx.i64("sphere_faces").data).into_iter().map(|t| [t[0] as usize, t[1] as usize, t[2] as usize]).collect();
    let cones: Vec<usize> = fx.i64("cone_ids").data.iter().map(|&c| c as usize).collect();
    let mesh = Mesh::from_parts(y, f, cones);
    // mesh combinatorics + base lengths
    let e_ref = fx.i64("edges").data;
    let e: Vec<i64> = mesh.e.iter().flat_map(|e| [e[0] as i64, e[1] as i64]).collect();
    assert_eq!(e, e_ref, "{name}: edges differ");
    let l0_ref = fx.f64("L0").data;
    assert!(nrel(&mesh.chord_lengths(), &l0_ref) < 1e-12, "{name}: L0 differs");

    let t_ref = fx.f64("t").data;
    let times = default_times(16.0, 12, 8, 0.05, 1e-3);
    assert!(nrel(&times, &t_ref) < 1e-14, "{name}: times differ {times:?}");
    let opts = FlowOptions { seed: seed_of(&fx), times: t_ref.clone(), ..Default::default() };
    let flow = DiscreteConeFlow::new(mesh, betas_of(&fx), opts);
    let (frames, flow) = flow.frames().unwrap_or_else(|e| panic!("{name}: {e}"));
    let nv = flow.nv();
    let ne = flow.mesh.e.len();
    let (u_ref, l_ref, k_ref, rho_ref) = (fx.f64("u").data, fx.f64("L").data, fx.f64("K_defect").data, fx.f64("rho_frame").data);
    let a_ref = fx.f64("A_vertex").data;
    assert_eq!(frames.len(), t_ref.len());
    let (mut eu, mut el, mut ek, mut er, mut gb, mut ea) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
    for (i, fr) in frames.iter().enumerate() {
        let du = nrel(&fr.u, &u_ref[i * nv..(i + 1) * nv]);
        let dl = nrel(&fr.l, &l_ref[i * ne..(i + 1) * ne]);
        let dk = nrel(&fr.k_defect, &k_ref[i * nv..(i + 1) * nv]);
        let da = nrel(&fr.a_vertex, &a_ref[i * nv..(i + 1) * nv]);
        let dr = ((fr.diag.rho - rho_ref[i]) / rho_ref[i]).abs();
        eu = eu.max(du);
        el = el.max(dl);
        ek = ek.max(dk);
        ea = ea.max(da);
        er = er.max(dr);
        gb = gb.max(fr.diag.gb_resid.abs());
    }
    eprintln!(
        "{name}: seed={:?} steps={} newton={} rejected={} fail={} corners={:?} | rel u {eu:.2e} L {el:.2e} K {ek:.2e} A {ea:.2e} rho {er:.2e} | GB max {gb:.2e}",
        flow.opts.seed, flow.stats.steps, flow.stats.newton, flow.stats.rejected, flow.stats.newton_fail, flow.seed_corners
    );
    (eu, el, ek, er, gb)
}

fn check(name: &str) {
    let (eu, el, ek, er, gb) = run_fixture(name);
    assert!(gb <= 1e-12, "{name}: GB {gb}");
    for (lab, v) in [("u", eu), ("L", el), ("K", ek), ("rho", er)] {
        assert!(v <= REL_TOL, "{name}: rel {lab} error {v:.3e} > {REL_TOL}");
    }
}

#[test]
fn n3_pos_gluing() {
    check("n3_pos_gluing");
}
#[test]
fn n4_pos_gluing() {
    check("n4_pos_gluing");
}
#[test]
fn n5_pos_gluing() {
    check("n5_pos_gluing");
}
#[test]
fn n4_pos_smart() {
    check("n4_pos_smart");
}
#[test]
fn n4_pos_polytope() {
    check("n4_pos_polytope");
}
/// rho < 0: same intrinsic DAE (C < 0); the 3D part (stage 2) is out of scope here.
#[test]
fn n4_neg_gluing_intrinsic() {
    check("n4_neg_gluing");
}
