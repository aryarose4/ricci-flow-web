//! rho < 0 True Polytope probe: `cargo run --release --example tp_neg -- LEVEL "b1,.." "x,y;.."`
use ricci_core::mesh::defects_and_areas;
use ricci_core::seedmesh::page_level;
use ricci_core::truepoly;
use std::f64::consts::PI;
use std::time::Instant;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let betas: Vec<f64> = a[2].split(',').map(|v| v.parse().unwrap()).collect();
    let pos: Vec<[f64; 2]> = a[3].split(';').map(|p| { let v: Vec<f64> = p.split(',').map(|z| z.parse().unwrap()).collect(); [v[0], v[1]] }).collect();
    for conf in [false, true] {
        let t0 = Instant::now();
        let r = if conf { truepoly::build(&pos, &betas, page_level(&a[1]).unwrap()) } else { truepoly::build_preview(&pos, &betas, page_level(&a[1]).unwrap()) };
        let tp = match r { Ok(t) => t, Err(e) => { println!("conf={conf} FAIL {e}"); continue; } };
        let g = defects_and_areas(&tp.mesh.f, &tp.mesh.opp, &tp.l0, tp.mesh.nv());
        let n = betas.len();
        let cerr = (0..n).map(|i| (g.k[i] - 2.0 * PI * (1.0 - betas[i])).abs()).fold(0.0, f64::max);
        let gb = (g.k.iter().sum::<f64>() - 4.0 * PI).abs();
        let flat = (n + tp.corners.len()..tp.mesh.nv()).map(|v| g.k[v].abs()).fold(0.0, f64::max);
        let ck: Vec<String> = tp.corners.iter().map(|&c| format!("{:.2}", g.k[c])).collect();
        let si = truepoly::self_intersections(&tp.x, &tp.mesh.f);
        println!("conf={conf} {:.2}s nv={} corners={} radii={:?} cone_err={:.1e} gb={:.1e} flat={:.1e} selfx={} match={:.1e} fit={:.1e} warn={:?} cornerK=[{}]",
            t0.elapsed().as_secs_f64(), tp.mesh.nv(), tp.corners.len(), tp.radii.iter().map(|r| (r * 100.0).round() / 100.0).collect::<Vec<_>>(), cerr, gb, flat, si, tp.match_err, tp.fit_resid, tp.warning, ck.join(","));
    }
}
