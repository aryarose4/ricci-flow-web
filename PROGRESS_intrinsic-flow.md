# PROGRESS — lane rust-intrinsic-flow (harness: claude, tier H)

## 2026-10-30 — intrinsic flow ported (candidate — pending user verification)
- core/src/sparse.rs: COO assembly + faer 0.24 sparse LU (default-features off: std + sparse-linalg, no rayon).
- core/src/mesh.rs: edges_of, edge_index (opp), signed_volume, Mesh::{from_sphere, from_parts, chord_lengths},
  face_geometry, defects_and_areas, cot_laplacian, area_jacobian (port of discrete_conformal.py).
- core/src/flow.rs: DiscreteConeFlow (smart/gluing/polytope seeds, Newton w/ exact sparse Jacobian + Sherman-Morrison
  rho term, gauge, step, adaptive BDF2 integrate, diagnostics, between, frames), FlowRun streaming iterator,
  WasmFlow wasm-bindgen wrapper.
- core/src/fixtures.rs: test-only npz loader (zip[deflate] + npyz; npyz's npz feature needs a C toolchain).
- core/tests/intrinsic_fixtures.rs: all 6 fixtures (incl. rho<0 n4_neg_gluing intrinsic) match per-frame u, L, K, A, rho
  to <= 1.7e-14 normwise rel (target 1e-8); GB max 9.2e-14; step/Newton counts identical to Python stats.
- wasm32-unknown-unknown release build OK (493 KB with faer linked via WasmFlow).
- Deferred: build_mesh mesh *construction* (seed_surface.build_seed needs a 3D convex hull of Fibonacci+ring points;
  hull.rs belongs to the realization lane). Mesh::from_sphere does the orientation/edges part; meshes come from fixtures.
- Not ported (display/3D only): intrinsic_vertex_data rounding (Frame.k_defect/a_vertex carry the raw values), cf_w raster.
