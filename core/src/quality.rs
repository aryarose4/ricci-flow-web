//! Status-bar quality note of the /flow page: port of `viewer/backend/flow_true.py`
//! `frame_quality` / `flow_quality`, plus the Polytope-seed opening-fold exemption of
//! `run_flow_true` (`seed_folds`).

pub const QUALITY_LEN_RMS: f64 = 0.01;
pub const QUALITY_LEN_MAX: f64 = 0.05;

/// The per-frame fields `frame_quality` reads.  `q_len_*` (when set) override `len_*`
/// (Polytope frame 0: quality reads the smooth fit, the display shows the polytope).
#[derive(Clone, Debug, Default)]
pub struct FrameRow {
    pub len_rel_max: Option<f64>,
    pub len_rel_rms: Option<f64>,
    pub q_len_rel_max: Option<f64>,
    pub q_len_rel_rms: Option<f64>,
    pub n_deep: usize,
    pub seed_folds: bool,
    pub embed_lost: bool,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Quality {
    pub label: String,
    pub ok: bool,
    pub detail: String,
}

/// `frame_quality(r)` -> (passes, reasons).
pub fn frame_quality(r: &FrameRow) -> (bool, Vec<&'static str>) {
    let lmax = r.q_len_rel_max.or(r.len_rel_max);
    let lrms = r.q_len_rel_rms.or(r.len_rel_rms);
    let mut reasons = Vec::new();
    if r.n_deep > 0 && !r.seed_folds {
        reasons.push("crumpled");
    }
    let len_ok = match lmax {
        Some(m) => m <= QUALITY_LEN_MAX && lrms.map_or(true, |s| s <= QUALITY_LEN_RMS),
        None => false,
    };
    if r.embed_lost || !len_ok {
        reasons.push("can't reach intrinsic metric");
    }
    (reasons.is_empty(), reasons)
}

fn pct(x: f64) -> String {
    // Python "%g" of 1.0 / 5.0
    let s = format!("{}", x);
    s
}

/// `flow_quality(rows)`.
pub fn flow_quality(rows: &[FrameRow]) -> Option<Quality> {
    if rows.is_empty() {
        return None;
    }
    let per: Vec<_> = rows.iter().map(frame_quality).collect();
    let (rms, mx) = (pct(100.0 * QUALITY_LEN_RMS), pct(100.0 * QUALITY_LEN_MAX));
    if per.iter().all(|p| p.0) {
        return Some(Quality {
            label: "Certified Ricci Flow".into(),
            ok: true,
            detail: format!("Every frame matches its metric (edge lengths within {} % rms / {} % max, no folds)", rms, mx),
        });
    }
    let last = per.last().unwrap();
    if last.0 {
        let bad = per.iter().filter(|p| !p.0).count();
        return Some(Quality {
            label: "Certified Target".into(),
            ok: true,
            detail: format!("The final shape matches the uniformized metric; {} of {} frames along the way do not", bad, per.len()),
        });
    }
    Some(Quality {
        label: format!("Best Attempt ({})", last.1.join(", ")),
        ok: false,
        detail: format!("The final shape does not meet the quality threshold ({} % rms / {} % max edge length, no folds)", rms, mx),
    })
}

/// `run_flow_true`'s Polytope-seed rule: the opening frames' folds are exempt (`seed_folds`)
/// while the fold count never grows; the first frame with 0 folds or more folds ends the run.
pub fn mark_seed_folds(rows: &mut [FrameRow], polytope_seed: bool) {
    let mut in_run = polytope_seed;
    let mut prev: Option<usize> = None;
    for r in rows.iter_mut() {
        if !in_run {
            break;
        }
        let nd = r.n_deep;
        if nd == 0 || prev.map_or(false, |p| nd > p) {
            in_run = false;
        } else {
            r.seed_folds = true;
            prev = Some(nd);
        }
    }
}
