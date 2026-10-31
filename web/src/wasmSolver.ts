import type { FlowFrame, FlowMeta, FlowRequest, Solver } from './solver';

/**
 * ===== WASM SOLVER STUB (to be filled in by the core lane) =====
 * Will wrap the wasm-pack `--target web` output in ../../core/pkg:
 *
 *   import init, { /* solver entry points *\/ } from '../../core/pkg/<crate>.js';
 *   await init();
 *   // build mesh, run BDF2 flow, call onFrame per output time with {t, faces, vertices, u, K, ...}
 *
 * The import is intentionally absent so the front end builds without core/pkg.
 * Once core/pkg exists, add it here and switch the default in worker.ts / main.ts.
 */
export class WasmSolver implements Solver {
  async run(_req: FlowRequest, _onMeta: (m: FlowMeta) => void, _onFrame: (f: FlowFrame) => void, _signal: AbortSignal): Promise<void> {
    throw new Error('WasmSolver is not wired up yet (core/pkg not built); use the fixture solver.');
  }
}
