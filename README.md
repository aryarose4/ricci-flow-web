# ricci-flow-web

Static, in-browser port of the `/flow` curvature-flow page from
`curvature-flow-and-geodesics`: Ricci flow on spheres with conic singularities.

- `core/` — Rust solver compiled to WebAssembly (`wasm-pack`).
- `web/` — Vite + TypeScript + Three.js front end; solver runs in a Web Worker.
- `fixtures/` — golden reference outputs exported from the Python pipeline; the
  Rust port is tested against them.

Deployed to GitHub Pages by `.github/workflows/pages.yml`.
