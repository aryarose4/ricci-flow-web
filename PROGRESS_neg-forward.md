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

## 2026-11-04 — verifier follow-up: exact commands, re-run, recomputed errors, smart/polytope
- Exact commands (repo root of the worktree; `cargo build --release --example crumple`; no env vars; mode `base` = current
  tracker, mode `F` = `anti.forward = ForwardOpts::default()` (w_guard 1e-2, step_rel .04, max_sub 8, hinge)):
  - c4:  `core/target/release/examples/crumple.exe OUT.json {base|F} gluing coarse ".3,.3,.3,.3" "0.7,0;0,0.7;-0.7,0;0,-0.7"`
  - n3:  `... {base|F} gluing coarse ".2,.3,.4" "0.7,0;-0.35,0.606218;-0.35,-0.606218"`
  - mx:  `... {base|F} gluing coarse ".25,.45,.35,.3" "0.6,0.1;-0.2,0.7;-0.6,-0.3;0.3,-0.5"`
  - smart / polytope variants: same as c4 with seed `smart` / `polytope`.
- Re-run reproduces c4 exactly: base failed_t .330, F failed_t 1.56, 0 deep on every frame both.
  Run JSONs committed: artifacts/neg-forward/runs/c4_base.json, c4_F.json.
- Harness now writes `len_rms_recomp` / `len_max_recomp` = edge error of the EMITTED vertices vs that frame's own
  target lengths. c4 F, k 17 (t 1.56, failed): row 1.151 % / 4.661 %, recomputed 1.151 % / 4.661 % (identical); K err
  row +4.54 % = mean of the emitted-vertex Length-Fit values +4.54 %. Baseline k 17: row says 1.15 % rms but the emitted
  (frozen) vertices are 30.5 % rms / 86.7 % max off this frame's lengths, K 25 % vs reported 4.4 % -> the old display
  under-reports by ~25x.
- smart / polytope square: the t = 0 cold start never passes (6 / 14 deep folds at frame 0, state "searching" on EVERY
  frame in both modes; failed_t NaN because tracking never started). First F version marched warm from that folded
  start, and the guard only forbids ADDING folds, so F ended with 2-3 deep vs base 0 (base cold-starts each searching
  frame). Fix: the march starts only after a frame has passed; while searching F uses the tracker's cold start.
  Now smart/polytope F == base per frame (identical states/deep). F is not better there; neither is faithful. The
  gluing seed is the only one of the three that starts clean for this config.
- After the fix: c4 / n3 / mx F unchanged (same states, deep, failed_t). cargo test --release passes.
