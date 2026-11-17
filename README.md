# ricci-flow-web

**Live:** https://aryarose4.github.io/ricci-flow-web/

Ricci flow on spheres with conic singularities, running entirely in the browser.
Place cone points, set their angles, and watch the metric flow toward constant
curvature, realized as a 3D surface at every frame.

## What it does

- **True Ricci Flow** — the cone-preserving Ricci flow, solved intrinsically on a
  triangulated sphere, with each frame realized in 3D and checked for isometry and folds.
- **Fast Uniformizer** — builds the final constant-curvature surface directly.
- **Initial conditions** — Cone Gluing, Smart Seed (flat-cone far field), or Polytope
  (a flat-faced polyhedron with corners at the cone points: convex for ρ > 0,
  star-shaped for ρ < 0).
- Works for both signs of ρ = 2 − Σ(1 − βᵢ).

## Layout

- `core/` — Rust solver compiled to WebAssembly (`wasm-pack`).
- `web/` — Vite + TypeScript + Three.js front end; the solver runs in a Web Worker.
- `fixtures/` — golden outputs from the original Python pipeline; the Rust port is
  tested against them.

## Build

```sh
wasm-pack build core --target web --release
cd web && npm ci && npm run dev      # or: npm run build
cargo test --manifest-path core/Cargo.toml
```

Pushes to `main` deploy to GitHub Pages via `.github/workflows/pages.yml`.
