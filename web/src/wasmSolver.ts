import init, { FlowSession, type FlowFrame as WFrame } from '../../core/pkg/ricci_core.js';
import type { FlowFrame, FlowMeta, FlowRequest, PreviewData, Quality, Solver } from './solver';

let ready: Promise<unknown> | null = null;

/** Cheap preview: build the session (mesh + cap seed) and read the preview polyhedron, without running the flow. */
export async function previewOf(req: FlowRequest): Promise<PreviewData> {
  await (ready ??= init());
  const t0 = performance.now();
  const s = new FlowSession(new Float64Array(req.positions.flat()), new Float64Array(req.betas), req.seed, req.level ?? 'coarse');
  try {
    return {
      vertices: Float32Array.from(s.preview_vertices()), faces: s.faces(), edges: s.edges(), coneIds: Array.from(s.cones()),
      rho: (req.betas.reduce((a, b) => a + b, 0) - (req.betas.length - 2)) / 2, clamped: s.clamped(), ms: performance.now() - t0,
    };
  } finally { s.free(); }
}

/** Live solver: drives a core FlowSession one chunk at a time inside the worker (terminate the worker to stop). */
export class WasmSolver implements Solver {
  async run(req: FlowRequest, onMeta: (m: FlowMeta) => void, onFrame: (f: FlowFrame) => void, signal: AbortSignal,
            onProgress?: (p: number, phase: string) => void): Promise<void> {
    onProgress?.(0, 'Loading solver');
    await (ready ??= init());
    const pos = new Float64Array(req.positions.flat());
    const s = new FlowSession(pos, new Float64Array(req.betas), req.seed, req.level ?? 'coarse');
    try {
      const faces = s.faces();
      let sent = false, faces0: Uint32Array | null = null;
      while (!s.done() && !signal.aborted) {
        const f = s.next_frame();
        if (!f) { onProgress?.(s.progress(), s.phase()); continue; }
        if (!sent) {
          sent = true;
          faces0 = faces;
          const sv = s.sphere_vertices();
          onMeta({
            request: req, rho: s.rho(), coneIds: Array.from(s.cones()),
            sphereVertices: Float32Array.from(sv), sphereFaces: faces, edges: s.edges(),
            nFrames: s.n_frames(), label: undefined, notice: s.warning() ?? null,
          });
        }
        const last = s.done();
        const lbl = s.quality_label();
        const quality: Quality | null = last && lbl ? { label: lbl, ok: s.quality_ok() ?? false, detail: s.quality_detail() ?? undefined } : null;
        onFrame(convert(f, faces0!, s.progress(), quality));
        f.free();
      }
    } finally { s.free(); }
  }
}

function convert(f: WFrame, sessionFaces: Uint32Array, progress: number, quality: Quality | null): FlowFrame {
  const kd = f.k_defect(), av = f.a_vertex(), K = new Float32Array(kd.length);
  for (let i = 0; i < kd.length; i++) K[i] = kd[i] / av[i];
  const ny = f.cf_lat(), nx = f.cf_lon();
  return {
    t: f.t(), faces: f.faces() ?? sessionFaces, vertices: Float32Array.from(f.vertices()), u: Float32Array.from(f.u()),
    K, Kdef: Float32Array.from(kd), Avert: Float32Array.from(av), rho: f.rho(), polytope: f.polytope(),
    wGrid: { nx, ny, values: Float32Array.from(f.cf_w()), extent: [-Math.PI, Math.PI, 0, Math.PI] },
    L: Float32Array.from(f.lengths()), gbError: f.gb_resid(), lenRelRms: f.len_rel_rms(), coneDefectErr: f.cone_defect_err(),
    quality, progress,
  };
}
