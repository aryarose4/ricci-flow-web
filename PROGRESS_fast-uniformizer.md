# PROGRESS — lane fast-uniformizer (Opus, 2026-11-04)

## What was ported
Python "Fast Uniformizer" (`flow_stream.py` fast path: metric solve -> GN/LM fading-bending realization; rho<0 `fast_neg`)
is ported in spirit, reusing Rust pieces:
- `core/src/fastunif.rs` (new): `uniformize()` = pseudo-transient continuation on the cone flow: backward-Euler
  steps with dt x8 per step (0.05 -> 1e10), each solved by the existing exact damped Newton (`flow.rs`, sparse
  Jacobian + Sherman-Morrison rho term), started from the far field `v_asymp` (`smart_u`). rho is an output; the
  mass term removes the dilation null direction. Stops at stationary residual max|rho A_v - K_v| < 1e-10.
  No explicit Euler anywhere. `admissibility()` ports `core/model.py::check_admissibility` (n=1 no; n=2 equal betas
  (Troyanov football); kappa=0 refused; kappa<0 ok; kappa>0 Luo-Tian kappa < 2 min beta).
- `pipeline.rs`: `Session::set_fast()` -> phases FastSolve/FastRealize, one output frame. rho>0: cold Fast LM fit
  from the cone-cap seed, plus sphere-scale restarts if it folds (cleanest kept, `cleaner_picks_b` now pub).
  rho<0: `tracker::cold_start` multi-start + curvature mismatch; label/notice honest (best attempt, not an embedding).
  Quality: "Certified Target" when the single shape passes the /flow thresholds, else "Best Attempt (...)".
  Admissibility guard now also applies to True Ricci Flow sessions (inadmissible configs are refused with a message).
- `wasm.rs`: `set_fast()`, `fast_stats()`.
- Web: Mode radio enabled (`mode` in FlowRequest), Initial Condition greyed in Fast mode, status line reports
  steps/Newton/residual/GB/cone err/time, tval "uniformized (Fast)". Number of Cones 2..8; n=2 sliders are linked
  (only equal angles admissible); preview failures (inadmissible) now show the core's message in #warn.
- Python's rho>0 "M9 + discrete-conformal correction" is NOT ported (not needed: the discrete metric is solved
  on the flow mesh directly). Python's `fast_neg` saddle realizer is NOT ported; rho<0 uses the TRF tracker cold fit.

## Tests
`cargo test --release` (GNU toolchain) all pass; new `core/tests/fast_unif.rs` (GB < 1e-9, cone err < 1e-9,
|K/rho-1| < 1e-8, n = 2,3,4,6,8, rho<0 n=5,8; inadmissible refused) + unit test in fastunif.rs.
tsc + vite build pass. Browser (`web/scripts/fast.mjs`, coarse, 0 page errors). Port 4392 was occupied by a
foreign server, so the check ran on `vite preview` port 4393.

| mode | n | betas | wall s | quality | len rms |
|---|---|---|---|---|---|
| fast | 2 | .6,.6 | 0.32 | Certified Target | 0.03 % |
| fast | 3 | .6,.7,.8 | 0.37 | Certified Target | 0.03 % |
| fast | 6 | .9,.85,.9,.8,.9,.85 | 1.11 | Certified Target | 0.02 % |
| fast | 8 | .95/.9 alternating | 1.89 | Certified Target | 0.02 % |
| fast | 5 (rho<0) | .3,.4,.5,.3,.4 | 10.5 | Best Attempt (crumpled) | 1.37 % |
| TRF | 2 | same | 5.03 | Certified Ricci Flow | 0.03 % |
| TRF | 3 | same | 7.14 | Certified Ricci Flow | 0.03 % |
| TRF | 6 | same | 13.8 | Certified Ricci Flow | 0.02 % |
| TRF | 8 | same | 21.0 | Certified Ricci Flow | 0.02 % |
Intrinsic certificate (all fast runs): GB ~1e-13, cone err ~1e-15, max|K/rho-1| <= 1.3e-10, 6-7 steps / 14-20 Newton.
Caveat: GB is exact by construction for any closed triangulation; the substantive checks are cone err and K = rho.
