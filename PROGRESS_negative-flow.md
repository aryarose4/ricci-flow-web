# PROGRESS — lane negative-flow (rho < 0 True Ricci Flow in the Rust core)

## 2026-11-03 — port complete (candidate — pending user verification)
- Worktree `ricci-flow-web-negative-flow`, branch `lane/negative-flow`, toolchain override stable-x86_64-pc-windows-gnu.
- Fixtures: new exporter `curvature-flow-and-geodesics/singular_ricci/dev/scratch/export_neg_flow_fixtures.py`
  (runs `flow_stream.run_flow`, realizer true_flow, resolution coarse, seed gluing; records streamed frames).
  `fixtures/n4_neg_square_coarse` (betas .3 x4 at (+-0.7,0),(0,+-0.7)) and `fixtures/n3_neg_small_coarse`
  (betas .2/.3/.4 on a radius-0.7 triangle). Python wall 19.4 s / 16.0 s.
- `core/src/tracker.rs`: Tracker (3-ring ball K, cone+1-ring exclusion, LEN_RMS_FAIL 5e-3, K_SIGNED_FAIL 0.10,
  CONFIRM 2, SEED_SCALES), `mismatch`, `cold_start`, `embed_notice` (+ Python `%g` formatter). No
  self-intersection count (off the criterion by default in Python).
- `realize.rs`: `FastRealizer::realize_cold_from(l, x0)` (cold fit from a scaled seed); rho>0 path unchanged.
- `pipeline.rs`: `Phase::Track` — one intrinsic frame + forward realization per chunk, frame returned
  immediately; quality note from the forward rows. rho = 0 still refused (`NEG_REFUSAL` reworded).
- `wasm.rs`: FlowFrame.embed_state / embed_lost / k_err_signed / k_err_rms; FlowSession.embed_notice / negative.
- Results (release): every frame's state sequence, len_rel_rms and K_err_signed equal Python to ~1e-14
  (same LM path); t 9e-16, rho 2e-15; GB, cone defect < 1e-12; notices identical; the old smoke
  `n4_neg_gluing` golden: states + notice + quality label identical.
- Timing per frame (coarse): n4 square native 0.94 s (0.16-0.56 s following, ~2.4 s once warm_full+polish),
  node/wasm 0.93 s; n3 small native 0.38 s, wasm 0.35 s. Python 0.92 / 0.76 s.
- `cargo test --release`: all suites pass (lib 4, intrinsic 6, mesh 1, realize 1, pipeline 3, negative 1).
  `rho_neg_is_refused` replaced by `rho_neg_runs_and_flat_is_refused`.
