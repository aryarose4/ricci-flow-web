// Wall time of the full FlowSession per rho>0 fixture in node (wasm32, --target web build).
// Usage (repo root): wasm-pack build core --target web --release && node core/bench/wasm_bench.mjs
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";
import { pathToFileURL } from "node:url";
const here = dirname(fileURLToPath(import.meta.url));
const pkg = join(here, "..", "pkg");
const { default: init, FlowSession } = await import(pathToFileURL(join(pkg, "ricci_core.js")).href);
await init({ module_or_path: readFileSync(join(pkg, "ricci_core_bg.wasm")) });
const fx = join(here, "..", "..", "fixtures");
for (const name of ["n3_pos_gluing", "n4_pos_gluing", "n4_pos_polytope", "n4_pos_smart", "n5_pos_gluing", "n4_neg_gluing"]) {
  const j = JSON.parse(readFileSync(join(fx, name, "fixture.json"), "utf8"));
  const r = j.meta.request;
  const t0 = performance.now();
  let s;
  try {
    s = new FlowSession(new Float64Array(r.positions.flat()), new Float64Array(r.betas), r.seed_surface, r.resolution);
  } catch (e) {
    console.log(`${name}: refused: ${e}`);
    continue;
  }
  let frames = 0, chunks = 0, maxChunk = 0, gbMax = 0;
  while (!s.done()) {
    const c0 = performance.now();
    const f = s.next_frame();
    maxChunk = Math.max(maxChunk, performance.now() - c0);
    chunks++;
    if (f) { frames++; gbMax = Math.max(gbMax, Math.abs(f.gb_resid())); f.free(); }
  }
  const wall = (performance.now() - t0) / 1000;
  console.log(`${name}: wall ${wall.toFixed(2)} s, frames ${frames}, chunks ${chunks}, max chunk ${maxChunk.toFixed(0)} ms, rho ${s.rho().toFixed(6)} (fixture ${j.scalars.rho.toFixed ? j.scalars.rho.toFixed(6) : j.scalars.rho}), GB ${gbMax.toExponential(1)}, quality "${s.quality_label()}" (fixture "${j.scalars.quality.label}")`);
  s.free();
}
