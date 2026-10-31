// Compute interface. Solvers run inside the Web Worker (worker.ts) and stream frames to the UI.

export type SeedSurface = 'gluing' | 'smart' | 'polytope';

export interface FlowRequest {
  n: number;
  betas: number[];            // cone angle / 2pi
  positions: [number, number][]; // stereographic-plane cone positions
  seed: SeedSurface;
  /** FixtureSolver only: force a specific fixture name. */
  fixture?: string;
}

export interface Quality {
  label: string;   // "Certified Ricci Flow" | "Certified Target" | "Best Attempt (...)"
  ok: boolean;
  detail?: string;
}

/** Static data sent once, before the first frame. */
export interface FlowMeta {
  request: FlowRequest;
  rho: number;
  coneIds: number[];
  sphereVertices: Float32Array; // nv*3, unit sphere (used for the conformal-factor discs)
  sphereFaces: Uint32Array;     // nf*3 (triangulation of the sphere mesh)
  edges: Uint32Array;           // ne*2
  nFrames: number | null;       // null if unknown
  label?: string;
  notice?: string | null;
}

export interface FlowFrame {
  t: number;
  faces: Uint32Array;     // nf*3 for this frame's 3D surface
  vertices: Float32Array; // nv*3, 3D realization
  u: Float32Array;        // nv, log conformal factor
  K: Float32Array;        // nv, Gaussian curvature estimate
  L?: Float32Array;       // ne, target edge lengths (for the Length Fit overlay)
  gbError: number;        // sum(K_defect) - 4pi
  lenRelRms?: number;     // 3D edge length vs target, relative RMS
  coneDefectErr?: number;
  quality: Quality | null; // non-null only on the last frame
  progress: number;       // 0..1
}

export interface Solver {
  /** Resolve once the whole flow has been streamed (or aborted). */
  run(req: FlowRequest, onMeta: (m: FlowMeta) => void, onFrame: (f: FlowFrame) => void, signal: AbortSignal): Promise<void>;
}

// worker protocol
export type ToWorker =
  | { cmd: 'start'; id: number; req: FlowRequest; solver: 'fixture' | 'wasm'; baseUrl: string }
  | { cmd: 'stop'; id: number };
export type FromWorker =
  | { type: 'meta'; id: number; meta: FlowMeta }
  | { type: 'frame'; id: number; frame: FlowFrame }
  | { type: 'done'; id: number; stopped: boolean }
  | { type: 'error'; id: number; message: string };
