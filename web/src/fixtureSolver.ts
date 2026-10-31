import type { FlowFrame, FlowMeta, FlowRequest, Solver } from './solver';

interface FixtureIndexEntry { name: string; n: number; betas: number[]; seed: string; rho: number; quality: string }

/** Streams precomputed frames (public/fixtures/*.json) with an artificial delay to mimic a live solve. */
export class FixtureSolver implements Solver {
  constructor(private baseUrl: string, private frameDelayMs = 140) {}

  async pick(req: FlowRequest): Promise<string> {
    if (req.fixture) return req.fixture;
    const idx: FixtureIndexEntry[] = await (await fetch(this.baseUrl + 'index.json')).json();
    let best = idx[0], bd = Infinity;
    for (const e of idx) {
      const sameRho = (e.rho > 0) === (req.betas.reduce((a, b) => a + b, 0) - (req.n - 2) > 0);
      let d = (e.n === req.n ? 0 : 100) + (e.seed === req.seed ? 0 : 10) + (sameRho ? 0 : 50);
      if (e.n === req.n) d += e.betas.reduce((a, b, i) => a + Math.abs(b - (req.betas[i] ?? 0)), 0);
      if (d < bd) { bd = d; best = e; }
    }
    return best.name;
  }

  async run(req: FlowRequest, onMeta: (m: FlowMeta) => void, onFrame: (f: FlowFrame) => void, signal: AbortSignal) {
    const name = await this.pick(req);
    const doc = await (await fetch(this.baseUrl + name + '.json', { signal })).json();
    const sv = new Float32Array(doc.sphereVertices);
    const faceSets: Record<string, Uint32Array> = {};
    for (const k of Object.keys(doc.faceSets)) faceSets[k] = new Uint32Array(doc.faceSets[k]);
    const nf = doc.frames.length;
    onMeta({
      request: { ...doc.request, betas: doc.request.betas, seed: doc.request.seed_surface, positions: doc.request.positions, fixture: name },
      rho: doc.rho, coneIds: doc.coneIds, sphereVertices: sv,
      sphereFaces: faceSets['0'] ?? Object.values(faceSets)[0], edges: new Uint32Array(doc.edges),
      nFrames: nf, label: doc.label, notice: doc.embedNotice,
    });
    for (let k = 0; k < nf; k++) {
      if (signal.aborted) return;
      await new Promise(r => setTimeout(r, this.frameDelayMs));
      if (signal.aborted) return;
      const f = doc.frames[k];
      onFrame({
        t: f.t, faces: faceSets[f.faceSet], vertices: new Float32Array(f.vertices), u: new Float32Array(f.u),
        K: new Float32Array(f.K), Kdef: new Float32Array(f.Kdef), Avert: new Float32Array(f.Avert), rho: f.rho ?? doc.rho,
        polytope: k === 0 && doc.request.seed_surface === 'polytope',
        wGrid: f.cfW && doc.cfGrid ? { nx: doc.cfGrid.nx, ny: doc.cfGrid.ny, values: new Float32Array(f.cfW), extent: doc.cfGrid.extent } : null,
        L: new Float32Array(f.L), gbError: f.gbError ?? 0,
        lenRelRms: f.lenRelRms ?? undefined, coneDefectErr: f.coneDefectErr ?? undefined,
        quality: k === nf - 1 ? doc.quality : null, progress: (k + 1) / nf,
      });
    }
  }
}
