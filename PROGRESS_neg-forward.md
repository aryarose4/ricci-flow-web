# PROGRESS — lane neg-forward (claude)

## 2026-11-04 — guarded forward march for rho < 0 (candidate — pending user verification)
- Worktree `ricci-flow-web-neg-forward`, branch `lane/neg-forward`, GNU toolchain override.
- `pipeline.rs::neg_forward_step` (on when `anti.forward` is Some and rho < 0): each frame is fitted WARM from the
  previously SHOWN shape through hidden flow metrics (step_rel .04, cheap warm fits that always advance), also past
  tracker failure (no cold restarts, no frozen "last good" shape). If the fit adds deep folds, it is re-fitted with a
  hinge fold prior toward the previous shape at w_guard x {1, 10, 100} (2x substeps after the first); the fewest deep
  folds win, ties by lower len rms. The prior is NOT dropped when the fit misses the length tolerance: fold-free beats
  smaller length error. Row len/K errors are those of the shown shape. Tracker states unchanged (honest); after failure
  the notice gets "Shown instead: the fold-free shape closest to the flow, now X% off in edge lengths (max Y%) and Z% off
  in curvature at t = ..." updated every frame.
- Web: wasmSolver calls `set_forward(1e-2, 0)` for rho < 0 (page default); `embedTag` says
  "(3D not faithful: best fold-free shape, X% length error)". curvature.ts unchanged: the Length Fit overlay already
  uses the frame's own vertices + Kdef/Avert, so it now compares the shown shape with its own frame's target.
- Baseline honesty bug found: the old tracker path shows the frozen last-good vertices but reports len/K errors of a
  hidden attempt (a different shape) — numbers and overlay disagree.
- Coarse results (tables: artifacts/neg-forward/TABLES.md; GIFs c4.gif n3.gif mx.gif):
  - c4 n4 b.3 square: base following to t .07, failed_t .33; new: following to t .33, failed_t 1.56; 0 deep both;
    final len rms 1.15 % / max 4.65 %, K rms 50 % (shown live). Wall 27 s -> 41 s (parallel runs).
  - n3 small (.2/.3/.4): both fail at t .003 by the K criterion (+10 %; likely coarse-ball noise). New shows the live
    fold-free shape: len rms <= .22 %, K err signed grows to +21..27 % late (base reports ~0-1 % for a hidden attempt
    while showing the t .002 shape). 0 deep both. Wall 11 -> 24 s.
  - mixed n4 (.25/.45/.35/.3): base accepts a 1-deep-fold frame at t .152 as "following" (failed_t .718); new keeps
    t .152 fold-free (0 deep) at .52 % rms, so it fails earlier (failed_t .33) — honest. New still has 1 deep fold on
    4 late frames (even w 100x guard). Wall 22 -> 96 s (guard escalation on every late frame; slow).
- Not done: Medium run, browser click-through on :4393 (only tsc + vite build). cargo test --release all pass.
