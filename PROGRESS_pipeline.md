# PROGRESS_pipeline (lane rust-pipeline, harness: claude, tier H)

## 2026-10-30 — candidate, pending user verification

- `core/src/seedmesh.rs` ports `build_seed`. That covers stereographic_to_sphere, the cap overlap clamp,
  Fibonacci points, rings, tip heights, `_tips_valid`, and the exclusion factors with the 0.9 shrink fallback.
- `seed_hull` is an exact port of `seed_surface._hull`. It emulates CPython's set iteration order (tuple hash
  plus set probing and resizing), so the faces come out in the same order as in Python.
- The same file also has `build_flow_mesh` (power-graded rings, `_orient_out`), `build_mesh` and `page_level`
  (PAGE_MESH).
- Mesh test: the meshes of all 6 fixtures are rebuilt from their inputs alone.
  - Vertices differ by at most 2.2e-16.
  - Faces are identical, including their order.
  - Cone ids match.
- `core/src/pipeline.rs` adds `Session`, a steppable port of `run_flow_true`. It contains:
  - `FlowSrc`, which implements FlowSource. `between` calls `flow.between`. `blend` returns the lengths of
    s u0 + (1-s) u(T). As in Python's `clean_start`, u(T) comes from one integration from 0 to 16, done lazily.
  - the clean start;
  - steppable `realize::Backward`, whose forward chain advances one fit per step;
  - Polytope frame 0, `mark_seed_folds` and `flow_quality`;
  - the cf rasters.

  rho <= 0 returns an Err instead of panicking.
- `core/src/cfgrid.rs` ports `smooth_w_grid` and `conformal_grid`. It uses a cubic RBF with a linear tail. The
  dense LU is factored once per point set. It handles the Polytope corner singularities at t = 0 and the
  keep-subset used when nv > 1500.
- `core/src/wasm.rs` adds `FlowSession` and `FlowFrame`. The API is documented in `core/API.md`.
- Changes to `realize.rs`:
  - `clean_start` now gets u(T) lazily through `src.blend(0.0)`. The result is exact.
  - `geom::edges_of` and `geom::signed_volume` delegate to `mesh`.
  - The realize test's length gate is now |rms - rms_fix| <= 10 % of rms_fix, with no absolute slack.
- Release test results: intrinsic 6/6, mesh 1/1, pipeline 3/3, realize 1/1. The 3 pipeline tests are the
  5-fixture end-to-end run, the rho<0 refusal, and the cap-overlap / bad-beta case.
- End-to-end results against the fixtures:
  - Relative error: u <= 2.8e-15, L <= 3.9e-15, K <= 1.6e-14, rho <= 1.8e-15.
  - GB <= 1.07e-13.
  - Length-RMS deviation is 0.00 % on every frame.
  - Deep-fold counts and quality labels are equal.
  - Step and Newton counts equal the fixture's `stats_after_run`, which counts the intrinsic frames only.
  - cf_w max abs difference is 5.0e-5, within the fixture's rounding.
- Wall time per smoke fixture, full session including cf rasters:

  | Fixture | Native (s) | wasm, node 24 (s) |
  |---|---|---|
  | n3_pos_gluing | 0.63 | 0.63 |
  | n4_pos_gluing | 1.17 | 1.09 |
  | n4_pos_polytope | 1.81 | 1.79 |
  | n4_pos_smart | 0.66 | 0.62 |
  | n5_pos_gluing | 3.23 | 3.18 |

  The wasm binary is 554 KB. On the coarse level, native takes 3.2 to 7.5 s, and 26.6 s for n5.
- Both sparse solvers are kept. The flow's Newton step uses faer's sparse LU: the matrix is nonsymmetric, with a
  Sherman-Morrison update, and this takes about 5 % of the wall time. The LM normal equations use the skyline
  Cholesky, which needs an SPD matrix. One cannot replace the other unless a nonsymmetric skyline LU is written.

## 2026-10-30 — quality present-None fix (verifier follow-up)
- quality.rs: `FrameRow.has_q` models Python key presence: with q_ keys present they are read verbatim (present None = no fallback); `.or()` fallback removed. Unit test `quality::tests::rms_limit_and_present_none` (rms-only failure on last frame -> Best Attempt; on earlier frame -> Certified Target; present-None q rms ignores the row rms).
- Fast-realizer rows keep len_rel_rms = None on purpose: Python FastRealizer diag (true_flow.py:569) has no len_rel_rms key (the :464 one is the GN realizer), so Python never applies the rms limit in fast mode; fixtures store NaN. All 5 labels still match.
