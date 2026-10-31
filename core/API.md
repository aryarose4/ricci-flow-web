# ricci_core WASM API

Build: `wasm-pack build core --target web --release` (from the repo root).

A streaming /flow "True Ricci Flow" session for rho > 0. It runs these steps in order:

1. Build the mesh from the cone positions.
2. Run the intrinsic cone-preserving flow (BDF2).
3. Compute the clean start.
4. Realize the frames in 3D, backward from the uniformized end (Fast LM).
5. Apply the Polytope frame 0, the seed-fold exemption and the quality note.
6. Emit per-frame rows with the conformal-factor raster.

It ports `viewer/backend/flow_true.run_flow_true` with the default `trf_realizer="fast"`. The native Rust API is
`ricci_core::pipeline::{Session, SessionConfig, FrameOut}`.

```ts
import init, { FlowSession, FlowFrame } from "./pkg/ricci_core.js";
await init();                                   // or init({ module_or_path: bytes })

type Seed = "gluing" | "smart" | "polytope";
type Level = "smoke" | "coarse" | "medium" | "fine" | "extra";   // page default: "coarse"

declare class FlowSession {
  /** positions: flat [re0, im0, re1, im1, ...] cone points in the complex (stereographic) plane;
   *  betas: cone angle / 2pi, each in [0.05, 0.99].  Throws a string for rho <= 0 (stage 2, not
   *  ported), betas out of range, length mismatch, unknown seed/level, or mesh failure. */
  constructor(positions: Float64Array, betas: Float64Array, seed: Seed, level: Level);

  /** One work chunk. Returns the next frame once the 3D realization is complete, else undefined. */
  next_frame(): FlowFrame | undefined;
  /** One work chunk without taking frames; false once everything is computed. */
  step(): boolean;
  done(): boolean;            // all frames taken
  progress(): number;         // 0..1
  phase(): string;            // human-readable phase
  n_frames(): number;         // 21 (default output times, t in [0, 16])

  rho(): number;              // uniformized curvature (NaN until the intrinsic flow finished)
  curvature_total(): number;  // C = 4pi - sum of cone defects (> 0)
  n_vertices(): number;
  sphere_vertices(): Float64Array;   // nv*3 unit-sphere solver mesh
  preview_vertices(): Float64Array;  // nv*3 seed polyhedron (sphere + tangent cone caps)
  faces(): Uint32Array;              // nf*3, outward
  edges(): Uint32Array;              // ne*2, a < b; order of FlowFrame.lengths()
  cones(): Uint32Array;              // vertex id of cone i
  seed_corners(): Uint32Array;       // Polytope seed's 4 extra corners (else empty)
  seed_betas(): Float64Array;        // cap betas actually used (differ if caps were clamped)
  clamped(): boolean;                // overlapping caps were shrunk

  /** The t = 0 start shape (the Python "seed" event), once the clean-start phase is over. */
  seed_frame(): FlowFrame | undefined;
  quality_label(): string | undefined;   // "Certified Ricci Flow" | "Certified Target" | "Best Attempt (...)"
  quality_detail(): string | undefined;
  quality_ok(): boolean | undefined;
  warning(): string | undefined;
  flow_stats(): Uint32Array;             // [steps, newton, rejected, newton_fail] (intrinsic frames)
  free(): void;
}

declare class FlowFrame {
  k(): number; t(): number;
  u(): Float64Array;          // log conformal factor per vertex (metric e^{2u} x round chord metric)
  k_defect(): Float64Array;   // angle defect per vertex; K = k_defect / a_vertex
  a_vertex(): Float64Array;   // barycentric vertex area
  lengths(): Float64Array;    // intrinsic edge lengths
  vertices(): Float64Array;   // nv*3 realized 3D shape
  faces(): Uint32Array | undefined;  // only when the triangulation differs (Polytope frame 0 hull)
  len_rel_max(): number; len_rel_rms(): number;   // 3D edge length vs this frame's metric
  n_deep(): number;           // edges folded inward by more than 0.15 rad
  n_reflex(): number; min_dihedral(): number;
  gb_resid(): number;         // sum K - 4pi (|.| <= 1e-12)
  cone_defect_err(): number; rho(): number; k_dev_max(): number;
  seed_folds(): boolean;      // folds exempt from quality (Polytope opening run)
  polytope(): boolean;        // frame shows the honest polytope (Polytope seed, t = 0)
  warning(): string | undefined;
  cf_w(): Float64Array;       // conformal-factor raster cf_lat x cf_lon, row-major
  cf_lat(): number;           // 64: rows lat in [0, pi] from the north pole (cell centres)
  cf_lon(): number;           // 128: cols lon in [-pi, pi)
  free(): void;
}
```

## Worker loop

```ts
const s = new FlowSession(pos, betas, "gluing", "coarse");
while (!s.done() && !stopRequested) {
  const f = s.next_frame();
  if (f) { postMessage({ frame: { t: f.t(), vertices: f.vertices(), cf: f.cf_w() /* ... */ } }); f.free(); }
  else postMessage({ progress: s.progress(), phase: s.phase() });
}
postMessage({ done: true, rho: s.rho(), quality: s.quality_label(), detail: s.quality_detail() });
```

Each call to `next_frame()` or `step()` does one of these chunks:

- one intrinsic frame;
- the clean start (up to about 8 LM fits);
- one backward-realized frame, together with its hidden in-between metrics;
- one fit of the forward chain;
- the Polytope fit;
- one frame emission, which includes computing the cf raster.

No frame is emitted until the whole backward realization has finished. This matches Python, which realizes the
frames starting from the uniformized end. The longest chunks measured were 0.8 s on smoke and about 7.6 s on
coarse (n5, a single LM fit). To stop in the middle of a chunk, terminate the Worker.

## Notes

- `cf_w` is the raster from `flow_true.conformal_grid`: the cone model plus a cubic RBF of the smooth
  remainder. To draw the north and south discs, sample the rows with lat < pi/2 and lat > pi/2.
- Mesh levels as (N, m, K):

  | Level | (N, m, K) |
  |---|---|
  | smoke | (150, 12, 2) |
  | coarse | (400, 16, 4) |
  | medium | (1200, 24, 8) |
  | fine | (2500, 32, 12) |
  | extra | (4000, 36, 14) |

- `WasmFlow` (intrinsic flow only, mesh supplied from JS) is still exported for compatibility.
