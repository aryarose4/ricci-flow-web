# PROGRESS — lane tp-negative (True Polytope for rho < 0)

Status: candidate — pending user verification.

## Approach
- A convex polytope can't carry more than 4 pi of curvature, so for rho < 0 the seed is a STAR-SHAPED (radial-graph) closed
  polyhedron, which is embedded by construction. `truepoly::solve_star` does this:
  - Directions = the cone directions + 4 farthest-point corners, plus the midpoint of every cone–cone edge of the spherical
    Delaunay (the hull of the unit vectors), repeated until no two cones are adjacent.
  - Corners stay at r = 1. Each cone's link is then corners only, so each cone radius is an independent 1D equation. The
    pyramid apex defect goes to 2 pi as r goes to infinity; it is solved by bracketing plus bisection, so the cone error is
    at machine precision.
  - The corners take the remaining 4 pi − Σ < 0 (saddles), so some crease edges are concave by construction.
- `build_raw` uses that polyhedron when Σ 2π(1−β) > 4π. There is no coplanar facet merge, and the facets are meshed by the
  existing crease-respecting lattice/Delaunay mesher. The initial sphere coordinates map each face barycentrically onto the
  flat triangle of its directions (instead of the radial projection), then project onto the sphere.
- pipeline.rs:
  - The rho < 0 refusal is removed.
  - In both rho < 0 paths (tracker and forward march), frame 0 is the polyhedron scaled to the gauge (`true_poly_start`,
    start = "true_polytope") instead of a cold start. The flow continues warm from it.
  - When True Polytope is selected, `n_deep` excludes crease edges (concave creases are geometry, not folds).
  - If the build fails for rho < 0, the session falls back to the Polytope seed and sets `warning` ("True Polytope
    unavailable for these cones (...); started from the Polytope seed instead"), which the UI shows as a notice.
- web/src/main.ts (seed selector only): True Polytope is no longer greyed out or replaced for rho < 0.

## Results
Built with `build_preview`, coarse. "Flat err" is the largest |defect| at the flat vertices.

| Config | Cone err | GB err | Flat err | Self-intersections |
|---|---|---|---|---|
| n3 (.2/.3/.4) | ≤ 8e-15 | ≤ 1.2e-13 | ≤ 8e-15 | 0 |
| n3 (.05 × 3, radii ~15.8) | ≤ 8e-15 | ≤ 1.2e-13 | ≤ 8e-15 | 0 |
| c4 (.3 × 4) | ≤ 8e-15 | ≤ 1.2e-13 | ≤ 8e-15 | 0 |
| mixed n4 | ≤ 8e-15 | ≤ 1.2e-13 | ≤ 8e-15 | 0 |
| n4 (.1/.1/.1/.95) | ≤ 8e-15 | ≤ 1.2e-13 | ≤ 8e-15 | 0 |
| n5 (.2 × 5) | ≤ 8e-15 | ≤ 1.2e-13 | ≤ 8e-15 | 0 |
| n6 (.1 … .35) | ≤ 8e-15 | ≤ 1.2e-13 | ≤ 8e-15 | 0 |
| n6 medium (2116 vertices) | ≤ 8e-15 | ≤ 1.2e-13 | ≤ 8e-15 | 0 |

Full flows (crumple harness, mode F, coarse). The gluing-seed numbers are the neg-forward baselines.

| Config | Frame 0 len rms | Frame 0 state | Where tracking ends | Gluing seed, same config | Final len rms | Wall |
|---|---|---|---|---|---|---|
| n3 | 1.5e-14 | following | follows to t = .152 | fails at t .003 | 8e-4 | 19 s |
| c4 | 3.8e-14 | following | fails at t .0029 | fails at t .33 | 1.6 % | 76 s |
| mx | 1.4e-14 | following | fails at t .0029 | fails at t .33 | 1.2 % | 153 s |

- The GB residual and the cone defect error are below 1e-10 on every frame (pipeline test).
- `cargo test --release`: all pass. New tests:
  - `truepoly::negative_star_polyhedra` (7 configs, 3–6 cones)
  - `truepoly::self_intersection_detects_overlap`
  - `pipeline_fixtures::true_polytope_rho_neg_frame0_is_the_polyhedron`
- `wasm-pack build core --target web --release` passes; `npm run build` passes.

## Limitations
- Conformal coordinates (`conformalize`) fail for most rho < 0 configs: the sphere layout flips 4–258 faces. In a probe on
  n5 the conformal factor u went down to −12; vertex-scaling uniformization of sharp cone tips degenerates.
  - The build then uses the direction-map positions and sets a warning. Cones sit exactly at the user's positions, but their
    true conformal positions are approximate, so for n ≥ 4 the end shape can differ from other seeds.
  - Only the c4 square (match 5e-7), the n3 small config and the mixed n4 configs got conformal coordinates. A real fix needs
    graded refinement at the tips or Delaunay-flip (Penner) discrete conformal maps.
- For c4 and mx the sharp concave creases make the forward tracker fail early (t .0029). The fold-free best shape is still
  shown, with an honest notice.
- `solve_star` fails if a cone is too flat for its link (β near 1 next to a large link curvature, r < 0.3 bracket), or if
  cones are antipodal. Both cases go to the Polytope fallback with a notice. No such admissible config was hit in the tests.
- No browser click-through was done; only tsc + vite.

## Images
artifacts/tp-negative/{n3,c4,mx}_k{0,5,10,20}.png and {n3,c4,mx}.gif. Run JSONs are in the same folder.
