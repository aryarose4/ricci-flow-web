//! Diagnostic: roughness of the plotted conformal factor for True Polytope vs other seeds.
//! `cargo run --release --example cf_rough -- <coarse|medium>`
use ricci_core::flow::Seed;
use ricci_core::pipeline::{Session, SessionConfig};
use ricci_core::seedmesh::page_level;

fn rough(w: &[f64], e: &[[usize; 2]], skip: &[bool]) -> (f64, f64, usize) {
    let n = w.len();
    let (mut s, mut c) = (vec![0.0; n], vec![0usize; n]);
    for &[i, j] in e {
        s[i] += w[j];
        c[i] += 1;
        s[j] += w[i];
        c[j] += 1;
    }
    let (mut ss, mut m, mut arg, mut k) = (0.0, 0.0, 0, 0);
    for v in 0..n {
        if skip[v] || c[v] == 0 {
            continue;
        }
        let d = (w[v] - s[v] / c[v] as f64).abs();
        ss += d * d;
        k += 1;
        if d > m {
            m = d;
            arg = v;
        }
    }
    ((ss / k as f64).sqrt(), m, arg)
}

fn main() {
    let lvl = std::env::args().nth(1).unwrap_or("coarse".into());
    let pos: Vec<[f64; 2]> = (0..5)
        .map(|k| {
            let a = 0.5 * std::f64::consts::PI + 2.0 * std::f64::consts::PI * k as f64 / 5.0;
            [0.7 * a.cos(), 0.7 * a.sin()]
        })
        .collect();
    let betas = vec![0.6, 0.7, 0.5, 0.8, 0.7];
    for seed in ["true_polytope", "polytope", "gluing"] {
        let cfg = SessionConfig {
            positions: pos.clone(),
            betas: betas.clone(),
            seed: Seed::parse(seed).unwrap(),
            level: page_level(&lvl).unwrap(),
            anti: Default::default(),
        };
        let mut s = Session::new(cfg).unwrap();
        let rows = s.run_all().unwrap();
        let m = &s.mesh;
        // skip cones and their 1-ring (genuine log singularity)
        let mut skip = vec![false; m.nv()];
        for &c in &m.cones {
            skip[c] = true;
        }
        for &[i, j] in &m.e {
            if m.cones.contains(&i) {
                skip[j] = true;
            }
            if m.cones.contains(&j) {
                skip[i] = true;
            }
        }
        let ncorner = s.true_poly.as_ref().map(|t| t.corners.len()).unwrap_or(0);
        let last = rows.last().unwrap();
        let (r_u, mx_u, a_u) = rough(&last.u, &m.e, &skip);
        print!(
            "{seed:14} nv={} corners={ncorner} final u: rms {r_u:.2e} max {mx_u:.2e} @v{a_u}",
            m.nv()
        );
        if let Some(tp) = &s.true_poly {
            let (r_o, mx_o, a_o) = rough(&tp.u_off, &m.e, &skip);
            let tot: Vec<f64> = last.u.iter().zip(&tp.u_off).map(|(a, b)| a + b).collect();
            let (r_t, mx_t, a_t) = rough(&tot, &m.e, &skip);
            let is_corner = |v: usize| tp.corners.contains(&v);
            // residual of the u_off fit per edge
            let ch = m.chord_lengths();
            let res: Vec<f64> =
                m.e.iter()
                    .enumerate()
                    .map(|(k, &[i, j])| (tp.l0[k] / ch[k]).ln() - 0.5 * (tp.u_off[i] + tp.u_off[j]))
                    .collect();
            let rr = (res.iter().map(|x| x * x).sum::<f64>() / res.len() as f64).sqrt();
            print!(" | u_off: rms {r_o:.2e} max {mx_o:.2e} @v{a_o}(corner {}) | plotted: rms {r_t:.2e} max {mx_t:.2e} @v{a_t}(corner {}) | u_off fit residual rms {rr:.2e}", is_corner(a_o), is_corner(a_t));
        }
        let vx = &last.vertices;
        let c = &m.cones;
        let d = |a: usize, b: usize| {
            let (p, q) = (vx[c[a]], vx[c[b]]);
            ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
        };
        let ds: Vec<f64> = (0..c.len())
            .flat_map(|a| ((a + 1)..c.len()).map(move |b| (a, b)))
            .map(|(a, b)| d(a, b))
            .collect();
        let sc = ds.iter().cloned().fold(0.0, f64::max);
        println!(
            "
   final cone-pair 3D distances / max: {:?}",
            ds.iter()
                .map(|x| format!("{:.3}", x / sc))
                .collect::<Vec<_>>()
        );
    }
}
