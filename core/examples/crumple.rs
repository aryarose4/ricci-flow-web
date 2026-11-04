//! Anti-crumple measurement harness (lane anti-crumple).
//! `cargo run --release --example crumple -- OUT.json MODE SEED LEVEL "b1,b2,.." "x1,y1;x2,y2;.."`
//! MODE: base | A | B | AB.  Overrides via env: AC_STEP, AC_MAXSUB, AC_CHEAP (0/1), AC_W0, AC_TAU, AC_HINGE (0/1).
//! Writes per-frame metrics + vertices (JSON) and prints a summary line.
use ricci_core::flow::Seed;
use ricci_core::hull;
use ricci_core::pipeline::{AntiCrumple, Session, SessionConfig};
use ricci_core::realize::geom::{self, V3};
use ricci_core::realize::{Adaptive, FoldPrior};
use ricci_core::seedmesh::page_level;
use std::fmt::Write as _;
use std::time::Instant;

fn envf(k: &str) -> Option<f64> {
    std::env::var(k).ok().and_then(|v| v.parse().ok())
}

struct M {
    n_reflex: usize,
    n_deep: usize,
    min_dih: f64,
    dent_depth: f64,
    dent_frac: f64,
    h_neg_frac: f64,
    jump: f64,
}

fn metrics(x: &[V3], f: &[[usize; 3]], cones: &[usize], prev: Option<&[V3]>) -> M {
    let nv = x.len();
    let mut lo = x[0];
    let mut hi = x[0];
    for p in x {
        for c in 0..3 {
            lo[c] = lo[c].min(p[c]);
            hi[c] = hi[c].max(p[c]);
        }
    }
    let diag = geom::norm(geom::sub(hi, lo));
    let fet = geom::face_edge_table(f);
    let d = geom::dihedrals(x, f);
    let n_reflex = d.iter().filter(|&&v| v < -1e-9).count();
    let n_deep = d.iter().filter(|&&v| v < -0.15).count();
    let min_dih = d.iter().cloned().fold(f64::INFINITY, f64::min);
    let mut h = vec![0.0; nv];
    for ((e, _), th) in fet.iter().zip(&d) {
        let l = geom::norm(geom::sub(x[e[0]], x[e[1]]));
        h[e[0]] += l * th;
        h[e[1]] += l * th;
    }
    let mut isc = vec![false; nv];
    for &c in cones {
        isc[c] = true;
    }
    let used: Vec<bool> = {
        let mut u = vec![false; nv];
        for t in f {
            for &v in t {
                u[v] = true;
            }
        }
        u
    };
    let nonc: Vec<usize> = (0..nv).filter(|&v| !isc[v] && used[v]).collect();
    let h_neg_frac = nonc.iter().filter(|&&v| h[v] < 0.0).count() as f64 / nonc.len() as f64;
    let (mut dent_depth, mut ndent) = (0.0f64, 0usize);
    if let Some(hl) = hull::convex_hull(x) {
        for &v in &nonc {
            let mut dmin = f64::INFINITY;
            for (n, b) in hl.normals.iter().zip(&hl.offsets) {
                dmin = dmin.min(b - geom::dot(*n, x[v]));
            }
            let dep = dmin.max(0.0) / diag;
            dent_depth = dent_depth.max(dep);
            if dep > 0.01 {
                ndent += 1;
            }
        }
    }
    let jump = match prev {
        Some(p) if p.len() == nv => {
            let a = geom::kabsch_align(x, p);
            a.iter().zip(p).map(|(q, r)| geom::norm(geom::sub(*q, *r))).fold(0.0, f64::max) / diag
        }
        _ => f64::NAN,
    };
    M { n_reflex, n_deep, min_dih, dent_depth, dent_frac: ndent as f64 / nonc.len() as f64, h_neg_frac, jump }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (out, mode, seed, level) = (&a[1], a[2].as_str(), &a[3], &a[4]);
    let betas: Vec<f64> = a[5].split(',').map(|v| v.parse().unwrap()).collect();
    let positions: Vec<[f64; 2]> = a[6]
        .split(';')
        .map(|p| {
            let v: Vec<f64> = p.split(',').map(|z| z.parse().unwrap()).collect();
            [v[0], v[1]]
        })
        .collect();
    let mut anti = AntiCrumple::default();
    if mode.contains('A') {
        let mut ad = Adaptive::default();
        if let Some(v) = envf("AC_STEP") {
            ad.step_rel = v;
        }
        if let Some(v) = envf("AC_MAXSUB") {
            ad.max_sub = v as usize;
        }
        if let Some(v) = envf("AC_CHEAP") {
            ad.cheap = v != 0.0;
        }
        anti.adaptive = Some(ad);
    }
    if mode.contains('C') {
        anti.rescue_final = envf("AC_RESCUE").map_or(4, |v| v as usize);
    }
    if mode.contains('F') {
        let mut fo = ricci_core::realize::ForwardOpts::default();
        if let Some(v) = envf("FW_STEP") {
            fo.step_rel = v;
        }
        if let Some(v) = envf("FW_MAXSUB") {
            fo.max_sub = v as usize;
        }
        if let Some(v) = envf("FW_W") {
            fo.w_guard = v;
        }
        if let Some(v) = envf("FW_HINGE") {
            fo.hinge = v != 0.0;
        }
        if mode.contains('H') {
            fo.hybrid_rescue = envf("AC_RESCUE").map_or(4, |v| v as usize);
        }
        anti.forward = Some(fo);
    }
    if mode.contains('B') {
        let mut p = FoldPrior::default();
        if let Some(v) = envf("AC_W0") {
            p.w0 = v;
        }
        if let Some(v) = envf("AC_TAU") {
            p.tau = v;
        }
        if let Some(v) = envf("AC_HINGE") {
            p.hinge = v != 0.0;
        }
        anti.prior = Some(p);
    }
    let cfg = SessionConfig { positions, betas, seed: Seed::parse(seed).unwrap(), level: page_level(level).unwrap(), anti: anti.clone() };
    let t0 = Instant::now();
    let mut s = Session::new(cfg).unwrap();
    let frames = s.run_all().unwrap();
    let wall = t0.elapsed().as_secs_f64();
    let mut js = String::new();
    let _ = write!(
        js,
        "{{\"mode\":\"{mode}\",\"anti\":\"{:?}\",\"wall\":{wall},\"n_hidden\":{},\"gated\":{},\"rho\":{},\"quality\":\"{}\",\"notice\":{:?},\"cones\":{:?},\"faces\":{:?},\"frames\":[",
        anti,
        s.n_hidden,
        s.prior_gated,
        s.rho(),
        s.quality.as_ref().map_or(String::new(), |q| q.label.clone()),
        s.embed_notice().unwrap_or_default(),
        s.mesh.cones,
        s.mesh.f
    );
    let mut prev: Option<Vec<V3>> = None;
    let (mut mx_deep, mut mx_dent, mut mx_jump, mut mx_rms, mut sum_hneg, mut mx_reflex) = (0usize, 0f64, 0f64, 0f64, 0f64, 0usize);
    let mut sum_jump = 0.0;
    let mut nj = 0;
    let mut failed_t = f64::NAN;
    for (i, fr) in frames.iter().enumerate() {
        let f = fr.faces.clone().unwrap_or_else(|| s.mesh.f.clone());
        let m = metrics(&fr.vertices, &f, &s.mesh.cones, prev.as_deref());
        if fr.embed_state == "failed" {
            failed_t = fr.t;
        }
        // frame 0 (seed / polytope) excluded from the maxima, as the quality rule exempts it
        if i > 0 {
            mx_deep = mx_deep.max(m.n_deep);
            mx_reflex = mx_reflex.max(m.n_reflex);
            mx_dent = mx_dent.max(m.dent_depth);
            mx_rms = mx_rms.max(fr.len_rel_rms);
            sum_hneg += m.h_neg_frac;
            if m.jump.is_finite() && fr.faces.is_none() && frames[i - 1].faces.is_none() {
                mx_jump = mx_jump.max(m.jump);
                sum_jump += m.jump;
                nj += 1;
            }
        }
        let v: Vec<String> = fr.vertices.iter().flatten().map(|z| format!("{:.5}", z)).collect();
        let fs = match &fr.faces {
            Some(ff) => format!("{:?}", ff),
            None => "null".into(),
        };
        let _ = write!(
            js,
            "{}{{\"k\":{},\"t\":{},\"len_rel_rms\":{},\"len_rel_max\":{},\"n_reflex\":{},\"n_deep\":{},\"min_dih\":{},\"dent_depth\":{},\"dent_frac\":{},\"h_neg_frac\":{},\"jump\":{},\"state\":\"{}\",\"start\":\"{}\",\"faces\":{},\"x\":[{}]}}",
            if i > 0 { "," } else { "" },
            fr.k,
            fr.t,
            fr.len_rel_rms,
            fr.len_rel_max,
            m.n_reflex,
            m.n_deep,
            m.min_dih,
            m.dent_depth,
            m.dent_frac,
            m.h_neg_frac,
            if m.jump.is_finite() { m.jump } else { -1.0 },
            fr.embed_state,
            fr.start,
            fs,
            v.join(",")
        );
        prev = Some(fr.vertices.clone());
    }
    js.push_str("]}");
    std::fs::write(out, js).unwrap();
    let nf = (frames.len() - 1).max(1) as f64;
    println!(
        "{mode}: wall {wall:.2}s hidden {} gated {} rescued {} from_bw {} | len_max {:.2e} | max deep {mx_deep} reflex {mx_reflex} | max dent {mx_dent:.4} | mean Hneg {:.3} | jump max {mx_jump:.4} mean {:.4} | max rms {:.2e} | {} failed_t {failed_t}",
        s.n_hidden,
        s.prior_gated,
        s.rescued,
        s.n_from_backward,
        frames.iter().map(|f| f.len_rel_max).fold(0.0, f64::max),
        sum_hneg / nf,
        sum_jump / nj.max(1) as f64,
        mx_rms,
        s.quality.as_ref().map_or(String::new(), |q| q.label.clone())
    );
}
