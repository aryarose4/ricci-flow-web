# PROGRESS — lane/tp-conformal (True Polytope conformal coordinates)

Status: candidate — pending user verification.

## Problem
`truepoly::build` used the RADIAL projection of the polytope as `mesh.y`. These are not
conformal coordinates of the polytope metric, so (1) the display offset `u_off` (LS fit of
log(l0/chord)) had per-edge residual rms ~0.1–0.2 and the plotted conformal factor was rough
(visible ring along the crease edges), and (2) the cones' conformal positions differed from the
user's, so True Polytope's final shape differed from gluing/polytope seeds.

## Fix (core/src/truepoly.rs only; pipeline.rs untouched)
- A. `conformalize`: uniformize the polytope metric with NO cones (every vertex smooth, existing
  `fastunif::uniformize`, pseudo-transient implicit Newton), lay the near-round lengths out on the
  unit sphere (`sphere_layout`, sparse Gauss-Newton from the radial positions), Möbius-normalise
  (`fit_mobius`/`mobius_apply`, LM on 6 real params, homogeneous coords) so cones best match the
  user's positions; `mesh.y` := these coordinates, `u_off` refit on them.
- B. `build_opts`: fixed point on the polytope's cone directions p_i ← normalize(p_i + (P_i − q_i))
  (q = conformal cone positions), ≤12 iterations, stop at 1e-6 or when stalled (<10% gain,
  err<1e-2: mesh re-generation noise floor). Over-relaxation 2.0 diverged (measured) → kept 1.0.
  Returns the best iterate. New TruePoly fields: fit_resid, match_err, match_iters.
- Crease test (`page_default_polytope`) still passes: t=0 metric is still exactly the polytope.

## Measurements (`cargo run --release --example tp_conf -- coarse`, before = main)
| config | fit resid before→after | plotted rough before→after | match err (iters) | max|dist ratio − gluing| before→after | build time before→after |
|---|---|---|---|---|---|
| n4 default | 1.20e-1 → 1.04e-3 | 4.74e-2 → 1.85e-2 | 4.3e-4 (9) | 0.0128 → 0.0096 | 0.01 → 1.9 s |
| n5 ring .6/.7/.5/.8/.7 r.7 | 1.98e-1 → 9.1e-4 | 4.44e-2 → 1.60e-2 | 1.7e-3 (11) | 0.0146 → 0.0029 | 0.01 → 3.4 s |
| n3 close pair .8×3 | 1.27e-1 → 8.3e-4 | 6.16e-2 → 1.54e-2 | 1e-16 (0, Möbius exact) | 0.0508 → 0.0508 | 0.01 → 0.2 s |
| n8 cube .85×8 | 6.13e-2 → 7.3e-4 | 3.46e-2 → 6.7e-3 | 3.4e-7 (6) | 0.0108 → 0.0014 | 0.01 → 2.8 s |
| n5 medium | 1.97e-1 → 2.3e-4 | 3.46e-2 → 8.3e-3 | 4.5e-4 (12) | 0.0135 → 0.0017 | 0.02 → 11.7 s |
Reference: gluing seed final-u roughness on n5 coarse is 2.3e-2, so TP's plotted factor is now
as smooth as the other seeds.

Honest caveats:
- B converges linearly (~0.7/step) to a ~1e-3 floor (n5 coarse 1.7e-3, n4 4e-4); the floor is
  re-meshing noise (facet meshes change with p). Remaining final-shape differences vs gluing
  (0.3–1%) are at the discretization level (they shrink coarse→medium: 0.0029→0.0017).
- n3 close pair: 3 cones are always Möbius-equivalent, so positions match exactly; the 5%
  distance-ratio gap vs gluing is unchanged by this lane — it is coarse-mesh discretization of
  the close pair (TP mesh ≠ gluing mesh), not a conformal-position error.
- Extra build time 0.2–3.4 s coarse, ~12 s medium (blocking in Session::new). Could be cut
  by fewer iterations (`MATCH_ITERS`) or warm-starting u across iterations.

## Screenshots (n5 ring, coarse, final frame)
- before: C:\Users\aryay\Desktop\ricci-flow-web-tp-conformal\web\review\tp_cf_before_disc.png
- after:  C:\Users\aryay\Desktop\ricci-flow-web-tp-conformal\web\review\tp_cf_after_disc.png
- gluing reference: C:\Users\aryay\Desktop\ricci-flow-web-tp-conformal\web\review\tp_cf_gluing_disc.png
(full pages: tp_cf_before.png / tp_cf_after.png / tp_cf_gluing.png in the same folder)

## Checks
cargo test --release: all pass. web `npm run build` (tsc + vite): pass. Browser (port 4396):
n5 TP "Certified Ricci Flow", rho 0.15090.

## Follow-up (verifier pass_with_caveats → merge blocker fixed)
- Merged main 49ef3d3 (clean).
- Preview path: `truepoly::build_preview` (radial, no conformal map / match loop) via
  `Session::new_preview` and wasm `FlowSession.preview(...)`; `previewOf` in web/src/wasmSolver.ts
  uses it. Only Run pays for the conformal build. Browser (port 4396, coarse, setConfig + 4 drag
  steps, `scripts/tp-preview.mjs`): TP n4 7–14 ms, n5 10–23 ms; gluing n4 16–66 ms, n5 22–43 ms.
- Full build speed: warm-starting uniformize from the radial fit was tried (dt0 = 1) — 7 vs 8
  steps and the next iteration failed to converge → reverted. Instead MATCH_ITERS 12 → 8 and stall
  stop at <15% gain. Coarse build: n4 1.6 s, n5 2.1 s, n3 0.2 s, n8 1.5 s (was 1.9/3.4/0.2/2.8).
  Cost: n5 match err 1.7e-3 → 2.7e-3, dist gap vs gluing 0.0029 → 0.0044 (still 3× better than main).
- Fallback: conformalize rejects non-finite layout, |ad−bc| ≤ 1e-8, any flipped sphere face, or
  non-finite u_off; on failure at the first iteration the radial build is returned with
  `TruePoly.warning` → `Session.warning` (surfaced by wasm `warning()`); on a later iteration the
  best conformal iterate is kept. Session::new no longer errors on a conformal-map failure.
- New test `conformal_vs_preview`. cargo test --release all pass; tsc + vite pass; browser full
  run n5 TP Certified, rho 0.15090.
- Medium n5 full build 11.7 → 7.6 s; match err 2.0e-3, dist gap 0.0018, fit resid 2.35e-4.
