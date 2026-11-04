# PROGRESS — lane forward-flow (claude)

## 2026-11-03 — step 0: setup + web
- Worktree `ricci-flow-web-forward-flow`, branch `lane/forward-flow` from main + merge lane/anti-crumple (conflict in pipeline.rs Session fields: kept both sides; true_polytope_compare example got `anti: Default::default()`). cargo test --release all pass.
- web: "True Polytope" seed radio (greyed, tooltip "rho > 0 only", falls back to Polytope when sum(1-beta) > 2); WasmSolver sets rescue_final(4) + fold prior (3e-4, 1, false) for rho>0. tsc + vite build OK.

## 2026-11-03 — step 1: opt-in forward realization (F)
- `realize::Forward` + `ForwardOpts {step_rel .04, max_sub 8, w_guard, hinge, hybrid_rescue}`; `AntiCrumple.forward` (default None). March warm from clean_start with hidden flow metrics; a fit that adds deep folds is retried with a hinge fold prior toward the previous shape and with 2x substeps; cleaner kept.
- Coarse first look: forward has NO early dent on c1 (0 deep for t<=3.39; backward's dent at t .012-.025 came from the backward/forward-fallback mix); with w_guard 1e-3 the late frames (t 7.4, 16) fold (3/5 deep). w_guard 1e-2: c1 coarse 0 deep, Certified.
- Measured (coarse, parallel runs, wall times inflated): c1 F w1e-2 0 deep, dent .031, Certified RF; F step .01 (max_sub 32) 0 deep, dent .013 (= backward BC .012) at ~2x wall; c2a F 0 deep Certified RF; c3b polytope F 0 deep Certified RF; true_polytope c1 F 0 deep, dent .014 vs backward BC .025 (forward better); c2b sharp F w1e-2 3 deep (backward 6 deep on 9 frames) -> still crumpled.
- Hybrid FH (forward, then backward+rescue and per-frame cleaner merge) implemented; c1 coarse never triggered it (forward already clean). c2b FH run lost (killed by mistake).
- Default w_guard -> 1e-2. Everything opt-in (`AntiCrumple.forward = None`); tests pass. Web not switched to forward (no wasm setter yet).
- c1 MEDIUM forward (w1e-2): 0 deep on every frame; early dent max .018 at t .012 (backward BC: 2 deep at t .025, dent .047). See SUMMARY in report.
