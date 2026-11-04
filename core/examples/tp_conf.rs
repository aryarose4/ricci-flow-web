//! True Polytope conformal-coordinate diagnostic over several configs.
//! `cargo run --release --example tp_conf -- <coarse|medium> [iters]`
use ricci_core::flow::Seed;
use ricci_core::pipeline::{Session, SessionConfig};
use ricci_core::seedmesh::page_level;
use std::f64::consts::PI;
use std::time::Instant;

fn rough(w: &[f64], e: &[[usize; 2]], skip: &[bool]) -> f64 {
    let n = w.len();
    let (mut s, mut c) = (vec![0.0; n], vec![0usize; n]);
    for &[i, j] in e {
        s[i] += w[j];
        c[i] += 1;
        s[j] += w[i];
        c[j] += 1;
    }
    let (mut ss, mut k) = (0.0, 0);
    for v in 0..n {
        if skip[v] || c[v] == 0 {
            continue;
        }
        let d = w[v] - s[v] / c[v] as f64;
        ss += d * d;
        k += 1;
    }
    (ss / k as f64).sqrt()
}

fn ring(n: usize, r: f64) -> Vec<[f64; 2]> {
    (0..n)
        .map(|k| {
            let a = 0.5 * PI + 2.0 * PI * k as f64 / n as f64;
            [r * a.cos(), r * a.sin()]
        })
        .collect()
}

fn dists(s: &Session, last: &ricci_core::pipeline::FrameOut) -> Vec<f64> {
    let vx = &last.vertices;
    let c = &s.mesh.cones;
    let mut ds = vec![];
    for a in 0..c.len() {
        for b in a + 1..c.len() {
            let (p, q) = (vx[c[a]], vx[c[b]]);
            ds.push(((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt());
        }
    }
    let m = ds.iter().cloned().fold(0.0, f64::max);
    ds.iter().map(|x| x / m).collect()
}

fn main() {
    let lvl = std::env::args().nth(1).unwrap_or("coarse".into());
    let only = std::env::args().nth(2);
    // cube: 8 vertices of a cube on the sphere, stereographic
    let cube: Vec<[f64; 2]> = {
        let s = 1.0 / 3f64.sqrt();
        let mut v = vec![];
        for &x in &[-s, s] {
            for &y in &[-s, s] {
                for &z in &[-s, s] {
                    // tilt a little so no vertex sits at the pole
                    let (y2, z2) = (y * 0.96 + z * 0.28, -y * 0.28 + z * 0.96);
                    v.push([x / (1.0 - z2), y2 / (1.0 - z2)]);
                }
            }
        }
        v
    };
    let cfgs: Vec<(&str, Vec<[f64; 2]>, Vec<f64>)> = vec![
        (
            "n4 default",
            vec![[0.7, 0.0], [-0.7, 0.0], [0.0, 0.7], [0.0, -0.7]],
            vec![0.6, 0.7, 0.5, 0.8],
        ),
        ("n5 ring", ring(5, 0.7), vec![0.6, 0.7, 0.5, 0.8, 0.7]),
        (
            "n3 close pair",
            vec![[0.15, 0.0], [-0.15, 0.0], [0.0, 1.2]],
            vec![0.8, 0.8, 0.8],
        ),
        ("n8 cube", cube, vec![0.85; 8]),
    ];
    for (name, pos, betas) in cfgs {
        if let Some(o) = &only {
            if !name.starts_with(o.as_str()) {
                continue;
            }
        }
        let mut ref_d = None;
        for seed in ["gluing", "true_polytope"] {
            let cfg = SessionConfig {
                positions: pos.clone(),
                betas: betas.clone(),
                seed: Seed::parse(seed).unwrap(),
                level: page_level(&lvl).unwrap(),
                anti: Default::default(),
            };
            let t0 = Instant::now();
            let mut s = match Session::new(cfg) {
                Ok(s) => s,
                Err(e) => {
                    println!("{name} {seed}: ERR {e}");
                    continue;
                }
            };
            let tb = t0.elapsed().as_secs_f64();
            let rows = s.run_all().unwrap();
            let last = rows.last().unwrap();
            let d = dists(&s, last);
            let m = &s.mesh;
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
            match &s.true_poly {
                None => {
                    println!("{name:14} gluing nv={} build {tb:.2}s", m.nv());
                    ref_d = Some(d);
                }
                Some(tp) => {
                    let tot: Vec<f64> = last.u.iter().zip(&tp.u_off).map(|(a, b)| a + b).collect();
                    let dd = ref_d
                        .as_ref()
                        .map(|r| {
                            r.iter()
                                .zip(&d)
                                .map(|(a, b)| (a - b).abs())
                                .fold(0.0, f64::max)
                        })
                        .unwrap_or(f64::NAN);
                    println!("{name:14} TP nv={} build {tb:.2}s | fit resid {:.2e} | plotted rough {:.2e} | match err {:.1e} iters {} | max |dist ratio - gluing| {dd:.4}",
                        m.nv(), tp.fit_resid, rough(&tot, &m.e, &skip), tp.match_err, tp.match_iters);
                }
            }
        }
    }
}
