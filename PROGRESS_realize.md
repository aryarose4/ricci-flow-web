# PROGRESS_realize (lane rust-realize, harness: claude)

## 2026-10-30 — candidate — pending user verification
Ported the rho>0 3D realization of the /flow True Ricci Flow to Rust (`core/src/{lm,hull,realize,quality}.rs`).

- `lm.rs`: envelope (skyline) Cholesky + RCM (no deps, wasm-safe); `fade_lm` (exact port: stages, mu schedule
  x8 / /4 floor 1e-7, 10 retries, volume-sign rejection, tol 1e-6, central-FD dihedral Jacobian as `bend_gn`);
  `least_squares_jac_scaled` replacing `least_squares(method="trf", x_scale="jac", max_nfev=30, ftol=xtol=gtol=1e-8,
  linear loss)` in `polytope_fit` (trust-region LM with growing jac column scaling, exact normal-equation steps
  instead of scipy's LSMR subproblem; analytic Jacobian instead of 2-point FD).
- `hull.rs`: incremental 3D hull (simplices + outward unit planes).
- `realize.rs`: geom helpers (`_dihedrals` engine value reproduced exactly incl. its unnormalised-edge atan2,
  Kabsch via Horn quaternion), `build_seed_on_mesh`/`fast_seed` (brentq -> bisection to 1e-15), `FastRealizer`,
  `n_deep_folds`, `clean_start`, `realize_backward`, `polytope_snap`, `polytope_fit`, `hull_faces`, `polytope_vertices`.
  Flow access via `FlowSource` trait (`between`, `blend`, optional `convex_anchor` — BI not ported, returns None).
- `quality.rs`: `frame_quality`, `flow_quality`, `mark_seed_folds` (Polytope opening-fold exemption).
- trf_track not ported (rho<0 path only).

Fixture test `core/tests/realize_fixtures.rs` (own npz loader via miniz_oxide; `between` = linear-in-t u interp
between output frames): all 5 rho>0 fixtures pass — rms equal to fixture to 3-4 digits on every frame, deep-fold
counts identical (incl. polytope 1-fold opening run and n5 frame 15's 2 folds), quality labels identical.
Kabsch vertex RMS vs fixture 1e-13..1e-3 (largest n5 frame 16, path difference from interpolated substeps).
Native release wall: n3 0.40 s, n4_gluing 0.93 s, n4_polytope 1.42 s, n4_smart 0.51 s, n5 2.95 s.
Debug `cargo test` ~75 s. wasm32 release build OK.
