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

## 2026-11-03 — step 3: C tail repair, B tuning, serial final matrix (coarse)
- C: cold rescue failed on medium case 1 (all cold fits of frames 16-19 fold; the folded end is not `ok`, so "warm from it" is cold too). Replaced the bounce with a tail repair: the first clean backward frame within the last 4 re-fits the later frames forward warm; a gated B is then switched on with the repaired end as reference. Medium case 1: end 60 deep / 8.2 % -> 0 deep / .07 %.
- B w0 1e-3 pushed max edge error over 5*len_tol (1.5 %) -> frame not ok -> next frame cold-started (broken chain). Default w0 now 3e-4 (max err <= ~1.2 %), plus a safeguard: a not-ok prior fit is redone without the prior.
- rho<0 at w0 3e-4: no gain (failed_t .330 unchanged); at w0 1e-3: .330 -> .718.
- Serial matrix (bin/crumple7.exe) in artifacts/anti-crumple/final/matrix.log; GIFs/PNGs in final/.

## 2026-11-03 — step 4: final matrix done (candidate — pending user verification)
- Tables: artifacts/anti-crumple/final/SUMMARY.md (serial wall times, bin/crumple7.exe). GIFs + worst-frame PNGs in final/.
- A: no fold/dent/jump gain on any case; +0-65 % wall (c1 medium +24 %, c3b +67 %). Recommend leaving A off.
- B (w0 3e-4): c2a 4 -> 0 deep, Certified; c3a dents .013 -> .009, H<0 .015 -> .006; max len err up to ~1.2-1.5 %. Gated off on folded ends. rho<0 needs w0 1e-3 to gain (failed_t .330 -> .718).
- C: the big win. c1 medium end 60 deep / 8.2 % -> 0, Best Attempt -> Certified Target (+35 % wall with B). c3b polytope 7 -> 0, Certified Ricci Flow.
- Not fixed: c1 medium dent at t .012-.025 (2 deep at .025, from the forward fallback chain); sharp cones beta <= .2 (c2b) fold on cone-incident edges.
