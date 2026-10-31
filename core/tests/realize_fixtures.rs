//! Golden-fixture test of the rho > 0 realization pipeline (`realize_backward` + `clean_start` +
//! Polytope frame 0 + quality note) against `fixtures/<name>/` (exported from the Python /flow).
//!
//! 3D vertices are path dependent, so invariants are compared: per-frame edge-length relative RMS
//! (vs the fixture's own, computed from its `vertices3d`), deep-fold counts, the quality label.
//! Kabsch-aligned vertex RMS vs the fixture is reported, not gated.
//!
//! The intrinsic flow is not run here: per-frame lengths come from the fixture, and the hidden
//! in-between metrics (`solver.between`) are approximated by interpolating u linearly in t between
//! the two bracketing output frames (the realized output frames are still fit to the exact lengths).

use ricci_core::quality::{flow_quality, mark_seed_folds, FrameRow};
use ricci_core::realize::{self, geom, FlowSource, MetricFrame, RealizeMesh};
use std::collections::HashMap;
use std::io::Read;
use std::path::PathBuf;
use std::time::Instant;

// ------------------------------------------------------------------ minimal npz loader
#[derive(Debug)]
struct Arr {
    shape: Vec<usize>,
    f: Vec<f64>,
}

fn u16le(b: &[u8], o: usize) -> usize {
    u16::from_le_bytes([b[o], b[o + 1]]) as usize
}
fn u32le(b: &[u8], o: usize) -> usize {
    u32::from_le_bytes([b[o], b[o + 1], b[o + 2], b[o + 3]]) as usize
}

fn parse_npy(b: &[u8]) -> Arr {
    assert_eq!(&b[..6], b"\x93NUMPY");
    let (hl, h0) = if b[6] == 1 { (u16le(b, 8), 10) } else { (u32le(b, 8), 12) };
    let hdr = std::str::from_utf8(&b[h0..h0 + hl]).unwrap();
    let descr = hdr.split("'descr':").nth(1).unwrap().split('\'').nth(1).unwrap().to_string();
    assert!(hdr.contains("'fortran_order': False"));
    let sh = hdr.split("'shape':").nth(1).unwrap();
    let sh = &sh[sh.find('(').unwrap() + 1..sh.find(')').unwrap()];
    let shape: Vec<usize> = sh.split(',').filter(|s| !s.trim().is_empty()).map(|s| s.trim().parse().unwrap()).collect();
    let d = &b[h0 + hl..];
    let n: usize = shape.iter().product();
    let f = match descr.as_str() {
        "<f8" => (0..n).map(|i| f64::from_le_bytes(d[8 * i..8 * i + 8].try_into().unwrap())).collect(),
        "<i8" => (0..n).map(|i| i64::from_le_bytes(d[8 * i..8 * i + 8].try_into().unwrap()) as f64).collect(),
        "<i4" => (0..n).map(|i| i32::from_le_bytes(d[4 * i..4 * i + 4].try_into().unwrap()) as f64).collect(),
        x => panic!("dtype {x}"),
    };
    Arr { shape, f }
}

fn load_npz(path: &PathBuf) -> HashMap<String, Arr> {
    let mut b = Vec::new();
    std::fs::File::open(path).unwrap().read_to_end(&mut b).unwrap();
    // end of central directory
    let eocd = (0..b.len() - 21).rev().find(|&i| &b[i..i + 4] == b"PK\x05\x06").unwrap();
    let n = u16le(&b, eocd + 10);
    let mut o = u32le(&b, eocd + 16);
    let mut out = HashMap::new();
    for _ in 0..n {
        assert_eq!(&b[o..o + 4], b"PK\x01\x02");
        let method = u16le(&b, o + 10);
        let csize = u32le(&b, o + 20);
        let (fl, xl, cl) = (u16le(&b, o + 28), u16le(&b, o + 30), u16le(&b, o + 32));
        let lho = u32le(&b, o + 42);
        let name = String::from_utf8(b[o + 46..o + 46 + fl].to_vec()).unwrap();
        o += 46 + fl + xl + cl;
        let (lfl, lxl) = (u16le(&b, lho + 26), u16le(&b, lho + 28));
        let data = &b[lho + 30 + lfl + lxl..lho + 30 + lfl + lxl + csize];
        let raw = match method {
            0 => data.to_vec(),
            8 => miniz_oxide::inflate::decompress_to_vec(data).unwrap(),
            m => panic!("zip method {m}"),
        };
        out.insert(name.trim_end_matches(".npy").to_string(), parse_npy(&raw));
    }
    out
}

fn v3s(a: &Arr) -> Vec<[f64; 3]> {
    a.f.chunks(3).map(|c| [c[0], c[1], c[2]]).collect()
}

// ------------------------------------------------------------------ flow stand-in
struct FixtureFlow {
    t: Vec<f64>,
    u: Vec<Vec<f64>>,
    l0: Vec<f64>,
    edges: Vec<[usize; 2]>,
}

impl FixtureFlow {
    fn lengths(&self, u: &[f64]) -> Vec<f64> {
        self.edges.iter().zip(&self.l0).map(|(e, l)| (0.5 * (u[e[0]] + u[e[1]])).exp() * l).collect()
    }
}

impl FlowSource for FixtureFlow {
    fn between(&self, frame: &MetricFrame, t_mid: f64) -> Vec<f64> {
        let k = frame.k;
        let (t0, t1) = (self.t[k], self.t[k + 1]);
        let s = (t_mid - t0) / (t1 - t0);
        let u: Vec<f64> = self.u[k].iter().zip(&self.u[k + 1]).map(|(a, b)| (1.0 - s) * a + s * b).collect();
        self.lengths(&u)
    }
    fn blend(&self, s: f64) -> Vec<f64> {
        let last = self.u.len() - 1;
        let u: Vec<f64> = self.u[0].iter().zip(&self.u[last]).map(|(a, b)| s * a + (1.0 - s) * b).collect();
        self.lengths(&u)
    }
}

fn rms(v: &[f64]) -> f64 {
    (v.iter().map(|x| x * x).sum::<f64>() / v.len() as f64).sqrt()
}

fn run_fixture(name: &str) -> bool {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("fixtures").join(name);
    let z = load_npz(&dir.join("arrays.npz"));
    let js: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(dir.join("fixture.json")).unwrap()).unwrap();
    let meta = &js["meta"];
    let betas: Vec<f64> = meta["request"]["betas"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
    let seed = meta["solver"]["seed"].as_str().unwrap().to_string();
    let corners: Vec<usize> = meta["solver"]["seed_corners"].as_array().unwrap().iter().map(|v| v.as_u64().unwrap() as usize).collect();
    let label_fix = js["scalars"]["quality"]["label"].as_str().unwrap().to_string();

    let y = v3s(&z["sphere_vertices"]);
    let faces: Vec<[usize; 3]> = z["sphere_faces"].f.chunks(3).map(|c| [c[0] as usize, c[1] as usize, c[2] as usize]).collect();
    let edges: Vec<[usize; 2]> = z["edges"].f.chunks(2).map(|c| [c[0] as usize, c[1] as usize]).collect();
    let cones: Vec<usize> = z["cone_ids"].f.iter().map(|&c| c as usize).collect();
    let (nf, nv, ne) = (z["t"].shape[0], y.len(), edges.len());
    let ts = z["t"].f.clone();
    let lf: Vec<Vec<f64>> = z["L"].f.chunks(ne).map(|c| c.to_vec()).collect();
    let uf: Vec<Vec<f64>> = z["u"].f.chunks(nv).map(|c| c.to_vec()).collect();
    let xf: Vec<Vec<[f64; 3]>> = z["vertices3d"].f.chunks(3 * nv).map(|c| c.chunks(3).map(|p| [p[0], p[1], p[2]]).collect()).collect();
    let nd_fix: Vec<usize> = z["n_deep"].f.iter().map(|&v| v as usize).collect();
    assert_eq!(geom::edges_of(&faces), edges, "edge order");

    let t0 = Instant::now();
    let mesh = RealizeMesh::new(faces.clone(), edges.clone(), cones.clone(), nv);
    let (x0, kind) = realize::fast_seed(&y, &faces, &cones, &betas);
    let rlz = realize::FastRealizer::new(mesh.clone(), x0, kind);
    let src = FixtureFlow { t: ts.clone(), u: uf.clone(), l0: z["L0"].f.clone(), edges: edges.clone() };
    let start = realize::clean_start(&rlz, &src, &lf[0], 6);
    let frs: Vec<MetricFrame> = (0..nf).map(|k| MetricFrame { k, t: ts[k], lengths: lf[k].clone() }).collect();
    let res = realize::realize_backward(&rlz, &src, &frs, start, &mut |_, _| {});
    let mut rows: Vec<FrameRow> = res
        .iter()
        .map(|r| FrameRow {
            len_rel_max: Some(r.len_rel_max),
            len_rel_rms: None, // FastRealizer diag has no len_rel_rms (as in Python)
            n_deep: realize::n_deep_folds(&r.x, &faces, 0.15),
            ..Default::default()
        })
        .collect();
    let mut xs: Vec<Vec<[f64; 3]>> = res.iter().map(|r| r.x.clone()).collect();
    let mut poly_faces = None;
    if seed == "polytope" {
        rows[0].q_len_rel_max = rows[0].len_rel_max;
        rows[0].q_len_rel_rms = rows[0].len_rel_rms;
        if let Some(p) = realize::polytope_vertices(&mesh, &res[0].x, &[cones.clone(), corners.clone()].concat(), &lf[0]) {
            rows[0].len_rel_max = Some(p.len_rel_max);
            xs[0] = p.x;
            poly_faces = Some(p.faces);
        }
    }
    mark_seed_folds(&mut rows, seed == "polytope");
    let q = flow_quality(&rows).unwrap();
    let wall = t0.elapsed().as_secs_f64();

    let mut ok = q.label == label_fix;
    eprintln!("== {name}: seed={seed} {:?}  label rust={:?} fixture={:?}  wall={:.2}s  polytope_faces={}", kind, q.label, label_fix, wall,
        poly_faces.as_ref().map_or(0, |f: &Vec<[usize; 3]>| f.len()));
    eprintln!("   k  t        rms_rust   rms_fix    max_rust   max_fix    deep r/f  kabsch_rms  start");
    for k in 0..nf {
        let rr = geom::len_rel(&xs[k], &edges, &lf[k]);
        let rf = geom::len_rel(&xf[k], &edges, &lf[k]);
        let (r_rms, f_rms) = (rms(&rr), rms(&rf));
        let r_max = rr.iter().cloned().fold(0.0, f64::max);
        let f_max = rf.iter().cloned().fold(0.0, f64::max);
        let al = geom::kabsch_align(&xs[k], &xf[k]);
        let krms = (al.iter().zip(&xf[k]).map(|(a, b)| geom::dot(geom::sub(*a, *b), geom::sub(*a, *b))).sum::<f64>() / nv as f64).sqrt();
        let len_ok = (r_rms - f_rms).abs() <= 0.10 * f_rms; // relative 10 %, no absolute slack
        let fold_ok = rows[k].n_deep == nd_fix[k];
        ok &= len_ok && fold_ok;
        eprintln!(
            "  {:2} {:8.4} {:10.3e} {:10.3e} {:10.3e} {:10.3e}  {}/{}{}  {:10.3e}  {}{}",
            k, ts[k], r_rms, f_rms, r_max, f_max, rows[k].n_deep, nd_fix[k], if fold_ok { " " } else { "!" }, krms, res[k].start,
            if len_ok { "" } else { "  LEN-FAIL" }
        );
    }
    ok
}

const POS: [&str; 5] = ["n3_pos_gluing", "n4_pos_gluing", "n4_pos_polytope", "n4_pos_smart", "n5_pos_gluing"];

#[test]
fn realize_rho_pos_fixtures() {
    let only = std::env::var("REALIZE_FIXTURE").ok();
    let mut bad = vec![];
    for name in POS {
        if only.as_deref().map_or(false, |o| o != name) {
            continue;
        }
        if !run_fixture(name) {
            bad.push(name);
        }
    }
    assert!(bad.is_empty(), "fixtures failing invariants: {bad:?}");
}
