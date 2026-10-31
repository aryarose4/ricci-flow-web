# PROGRESS web-ui (harness: claude, tier M)
- Vite 8.3.5 + TS 5.9.3 + three 0.187.0 (pinned), base './'. web/src: solver.ts (Solver interface), worker.ts, fixtureSolver.ts, wasmSolver.ts (stub), view3d.ts, cfview.ts, shading.ts, main.ts.
- web/scripts/fixtures-to-json.mjs (fflate + minimal npy reader) runs as prebuild -> web/public/fixtures/*.json (gitignored, regenerated).
- npm ci && npm run build OK. npm run screenshots (headless Chromium via playwright, swiftshader WebGL) plays all 6 fixtures, 21 frames each, writes web/review/*.png (light+dark).
- Preview port 4391 (4173 is taken by an unrelated local server on this machine).
- Known gaps: conformal-factor discs are vertex-coloured mesh (jagged equator on the coarse smoke mesh) rather than the Python page's raster; positions inputs and Mesh Resolution not wired (fixture mode); last-frame edge accuracy is "-" when fixture value is NaN.
