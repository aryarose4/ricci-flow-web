//! Status-bar quality note of the /flow page: port of `viewer/backend/flow_true.py`
//! `frame_quality` / `flow_quality`, plus the Polytope-seed opening-fold exemption of
//! `run_flow_true` (`seed_folds`).

pub const QUALITY_LEN_RMS: f64 = 0.01;
pub const QUALITY_LEN_MAX: f64 = 0.05;

/// The per-frame fields `frame_quality` reads, with Python dict semantics:
/// `None` = key present with value None (or absent: `r.get` gives None either way for `len_*`).
/// `has_q` = the `q_len_rel_*` keys are PRESENT in the row (Polytope frame 0); then they are read
/// verbatim, even when None (`r.get("q_len_rel_rms", ...)` returns a present None, no fallback).
/// Without `has_q` the `len_rel_*` values are read.
#[derive(Clone, Debug, Default)]
pub struct FrameRow {
    pub len_rel_max: Option<f64>,
    pub len_rel_rms: Option<f64>,
    pub q_len_rel_max: Option<f64>,
    pub q_len_rel_rms: Option<f64>,
    pub has_q: bool,
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
    let (lmax, lrms) = if r.has_q { (r.q_len_rel_max, r.q_len_rel_rms) } else { (r.len_rel_max, r.len_rel_rms) };
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

#[cfg(test)]
mod tests {
    use super::*;
    fn row(max: f64, rms: Option<f64>) -> FrameRow {
        FrameRow { len_rel_max: Some(max), len_rel_rms: rms, ..Default::default() }
    }
    #[test]
    fn rms_limit_and_present_none() {
        // passes max (4 %) but fails rms (2 % > 1 %) on the final frame -> Best Attempt
        let q = flow_quality(&[row(0.01, Some(0.001)), row(0.04, Some(0.02))]).unwrap();
        assert_eq!(q.label, "Best Attempt (can't reach intrinsic metric)");
        // rms failing only on an earlier frame -> Certified Target
        let q = flow_quality(&[row(0.04, Some(0.02)), row(0.01, Some(0.001))]).unwrap();
        assert_eq!(q.label, "Certified Target");
        // rms None (FastRealizer diag has no len_rel_rms) -> only max is checked
        assert!(frame_quality(&row(0.04, None)).0);
        // present-None q keys: the row's own (failing) rms is NOT read (no fallback)
        let mut r = row(0.12, Some(0.02));
        r.has_q = true;
        r.q_len_rel_max = Some(0.003);
        r.q_len_rel_rms = None;
        assert!(frame_quality(&r).0);
        // present q rms that fails the limit -> fails
        r.q_len_rel_rms = Some(0.02);
        assert!(!frame_quality(&r).0);
    }
}
