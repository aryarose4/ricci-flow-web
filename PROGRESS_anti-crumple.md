# PROGRESS — lane anti-crumple (claude)

## 2026-11-03 — step 1: opt-in flags wired
- Worktree `ricci-flow-web-anti-crumple`, branch `lane/anti-crumple`, GNU toolchain override.
- `SessionConfig.anti: AntiCrumple { adaptive: Option<realize::Adaptive>, prior: Option<realize::FoldPrior> }`, default off.
- (A) `Backward` hidden-metric rule configurable (port rule 0.08 / 4); `cheap` = WARM-only LM for hidden metrics. rho<0: hidden metrics between consecutive frames via `flow.between`.
- (B) `FadeSetup.prior` (+w, hinge) adds w * sum (theta - theta_ref)^2 to fade_lm; reference = realized final frame (rho>0) / previous accepted frame (rho<0).
- `cargo test --release` all pass with defaults.

## 2026-11-03 — step 2: harness + first measurements (coarse unless noted; runs in parallel, so wall times here are NOT the timing numbers)
- Harness `core/examples/crumple.rs` (per-frame len rms/max, reflex, deep (< -0.15), dent depth below convex hull / bbox diag, frac non-cone vertices with H<0, rigid-aligned max jump / bbox diag; JSON with vertices). Renderer `artifacts/anti-crumple/tools/render.py`.
- Root cause found: on several configs the cold fit of the UNIFORMIZED END folds (port has no BI convex anchor). Case 1 medium baseline: final frame 60 deep folds, max len err 8.2 %, and the chain carries a dent (2 deep at t=.025, min dih -.18). p3 (polytope, case-1 positions, coarse): final 7 deep.
- A (step .02, cap 16, cheap): no measurable change on any case (folds, dents, jumps identical to ~1e-3). Metric jumps between frames are not the driver here.
- New opt-in C `rescue_final`: cold-fit frames n-2.. (<=4) until clean, chain warm to the end. Fixes p3 and case 1 coarse (0 deep on every frame, Certified). Failed on medium case 1 (all cold fits fold) -> added bounce fallback (warm backward from the folded end until clean, then forward). Rerunning.
- B (w0 1e-3, tau 1): x1 (n3 b .25/.4/.5) 4 -> 0 deep, early rms .03 -> .16 %, max .1 -> .9 % (w0 3e-4: 0 deep, max .6 %). "Release" variant (refit without prior) restores the folds -> dropped.
- B with a FOLDED final reference (x5, x6 sharp cones) is harmful: drags folds to early frames, max len err ~1-1.5 %. -> gate: B is off when the final frame has deep folds (`prior_gated`).
- rho<0 (c4): B with reference = previous accepted frame moves tracker failure t .330 -> .718 (one more frame followed). A: no change.
- Sharp cones (beta <= .2): late folds sit on cone-incident edges (outside the prior's quads); none of A/B/C fix them (known gap).
