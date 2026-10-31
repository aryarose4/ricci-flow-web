/// <reference lib="webworker" />
import type { FromWorker, ToWorker } from './solver';
import { FixtureSolver } from './fixtureSolver';
import { WasmSolver, previewOf } from './wasmSolver';

const ctl = new Map<number, AbortController>();
const post = (m: FromWorker, tr: Transferable[] = []) => (self as unknown as Worker).postMessage(m, tr);

self.onmessage = async (ev: MessageEvent<ToWorker>) => {
  const m = ev.data;
  if (m.cmd === 'preview') {
    try { post({ type: 'preview', id: m.id, data: await previewOf(m.req) }); }
    catch (e) { post({ type: 'preview-fail', id: m.id, message: String((e as Error).message ?? e) }); }
    return;
  }
  if (m.cmd === 'stop') { ctl.get(m.id)?.abort(); return; }
  const ac = new AbortController();
  ctl.set(m.id, ac);
  try {
    const solver = m.solver === 'wasm' ? new WasmSolver() : new FixtureSolver(m.baseUrl);
    await solver.run(m.req,
      meta => post({ type: 'meta', id: m.id, meta }),
      frame => post({ type: 'frame', id: m.id, frame }),
      ac.signal,
      (progress, phase) => post({ type: 'progress', id: m.id, progress, phase }));
    post({ type: 'done', id: m.id, stopped: ac.signal.aborted });
  } catch (e) {
    if (ac.signal.aborted) post({ type: 'done', id: m.id, stopped: true });
    else post({ type: 'error', id: m.id, message: String((e as Error).message ?? e) });
  } finally { ctl.delete(m.id); }
};
