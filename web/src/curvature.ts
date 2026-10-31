// Port of viewer/frontend/flow_curvature.py (CURV_JS): Length Fit / Gaussian Curvature overlay estimators.
import { pickLut } from './shading';
import type { FlowFrame } from './solver';

const KR = 3;            // ball radius in rings
const RQ = 0.01, BUF = 0.15;
export const SENSITIVITY: [string, string][] = [['0.65', 'Low'], ['0.4', 'Medium'], ['0.15', 'High']];
export const SENSITIVITY_DEFAULT = '0.4';

interface Topo {
  faces: Uint32Array; nv: number; tips: number[];
  ex: Uint8Array; ball: (Int32Array | null)[]; vf: number[][]; fok: Uint8Array; edges?: { a: number; b: number; fs: number[] }[];
}
let TOPO: Topo | null = null;

function topology(faces: Uint32Array, nv: number, tips: number[]): Topo {
  if (TOPO && TOPO.faces === faces && TOPO.nv === nv && TOPO.tips === tips) return TOPO;
  const nf = faces.length / 3, nb: number[][] = Array.from({ length: nv }, () => []);
  for (let q = 0; q < nf; q++) for (let k = 0; k < 3; k++) { const a = faces[3 * q + k], b = faces[3 * q + (k + 1) % 3]; nb[a].push(b); nb[b].push(a); }
  for (let i = 0; i < nv; i++) nb[i] = Array.from(new Set(nb[i]));
  const ex = new Uint8Array(nv);
  for (const t of tips) { ex[t] = 1; for (const u of nb[t]) ex[u] = 1; }
  const ball: (Int32Array | null)[] = new Array(nv), stamp = new Int32Array(nv).fill(-1);
  for (let v = 0; v < nv; v++) {
    if (ex[v]) { ball[v] = null; continue; }
    const res = [v]; stamp[v] = v; let cur = [v];
    for (let d = 0; d < KR; d++) {
      const nxt: number[] = [];
      for (const c of cur) for (const w of nb[c]) if (stamp[w] !== v) { stamp[w] = v; nxt.push(w); res.push(w); }
      cur = nxt;
    }
    ball[v] = Int32Array.from(res);
  }
  const vf: number[][] = Array.from({ length: nv }, () => []);
  for (let q = 0; q < nf; q++) for (let k = 0; k < 3; k++) vf[faces[3 * q + k]].push(q);
  const fok = new Uint8Array(nf);
  for (let q = 0; q < nf; q++) fok[q] = (ex[faces[3 * q]] || ex[faces[3 * q + 1]] || ex[faces[3 * q + 2]]) ? 0 : 1;
  TOPO = { faces, nv, tips, ex, ball, vf, fok };
  return TOPO;
}

/** Angle defect over barycentric area, 3-ring ball averaged (NaN on the cone vertices and their 1-ring). */
export function estimateIntr(V: Float32Array, faces: Uint32Array, tips: number[], s: number): Float64Array {
  const nv = V.length / 3, T = topology(faces, nv, tips), nf = faces.length / 3;
  const defect = new Float64Array(nv).fill(2 * Math.PI), area = new Float64Array(nv), K = new Float64Array(nv);
  for (let q = 0; q < nf; q++) {
    const f = [faces[3 * q], faces[3 * q + 1], faces[3 * q + 2]];
    let ar = 0;
    for (let k = 0; k < 3; k++) {
      const A = 3 * f[k], B = 3 * f[(k + 1) % 3], C = 3 * f[(k + 2) % 3];
      const ux = V[B] - V[A], uy = V[B + 1] - V[A + 1], uz = V[B + 2] - V[A + 2], vx = V[C] - V[A], vy = V[C + 1] - V[A + 1], vz = V[C + 2] - V[A + 2];
      const cx = uy * vz - uz * vy, cy = uz * vx - ux * vz, cz = ux * vy - uy * vx, cr = Math.sqrt(cx * cx + cy * cy + cz * cz), dt = ux * vx + uy * vy + uz * vz;
      defect[f[k]] -= Math.atan2(cr, dt);
      if (k === 0) ar = 0.5 * cr / 3;
      area[f[k]] += ar;
    }
  }
  for (let i = 0; i < nv; i++) {
    const bl = T.ball[i]; if (!bl) { K[i] = NaN; continue; }
    let sd = 0, sa = 0;
    for (const u of bl) { if (T.ex[u]) continue; sd += defect[u]; sa += area[u]; }
    K[i] = sa > 0 ? sd * s * s / sa : NaN;
  }
  return K;
}

/** Gauss-map area ratio over the 3-ring ball (signed; saddles negative). */
export function estimateExtr(V: Float32Array, faces: Uint32Array, tips: number[], s: number): Float64Array {
  const nv = V.length / 3, nf = faces.length / 3, T = topology(faces, nv, tips);
  const Nn = new Float64Array(3 * nv), om = new Float64Array(nf), ar = new Float64Array(nf), K = new Float64Array(nv);
  for (let q = 0; q < nf; q++) {
    const A = 3 * faces[3 * q], B = 3 * faces[3 * q + 1], C = 3 * faces[3 * q + 2];
    const fx = (V[B + 1] - V[A + 1]) * (V[C + 2] - V[A + 2]) - (V[B + 2] - V[A + 2]) * (V[C + 1] - V[A + 1]);
    const fy = (V[B + 2] - V[A + 2]) * (V[C] - V[A]) - (V[B] - V[A]) * (V[C + 2] - V[A + 2]);
    const fz = (V[B] - V[A]) * (V[C + 1] - V[A + 1]) - (V[B + 1] - V[A + 1]) * (V[C] - V[A]);
    ar[q] = 0.5 * Math.sqrt(fx * fx + fy * fy + fz * fz);
    for (let i = 0; i < 3; i++) { const v = faces[3 * q + i]; Nn[3 * v] += fx; Nn[3 * v + 1] += fy; Nn[3 * v + 2] += fz; }
  }
  for (let i = 0; i < nv; i++) {
    const l = Math.sqrt(Nn[3 * i] ** 2 + Nn[3 * i + 1] ** 2 + Nn[3 * i + 2] ** 2) || 1e-300;
    Nn[3 * i] /= l; Nn[3 * i + 1] /= l; Nn[3 * i + 2] /= l;
  }
  for (let q = 0; q < nf; q++) {
    if (!T.fok[q]) continue;
    const a = 3 * faces[3 * q], b = 3 * faces[3 * q + 1], c = 3 * faces[3 * q + 2];
    const ax = Nn[a], ay = Nn[a + 1], az = Nn[a + 2], bx = Nn[b], by = Nn[b + 1], bz = Nn[b + 2], cx = Nn[c], cy = Nn[c + 1], cz = Nn[c + 2];
    const num = ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) + az * (bx * cy - by * cx);
    const den = 1 + (ax * bx + ay * by + az * bz) + (bx * cx + by * cy + bz * cz) + (cx * ax + cy * ay + cz * az);
    om[q] = 2 * Math.atan2(num, den);
  }
  const mark = new Int32Array(nv).fill(-1), fst = new Int32Array(nf).fill(-1);
  for (let i = 0; i < nv; i++) {
    const bl = T.ball[i]; if (!bl) { K[i] = NaN; continue; }
    let sw = 0, sa = 0;
    for (const u of bl) mark[u] = i;
    for (const u of bl) for (const fq of T.vf[u]) {
      if (fst[fq] === i || !T.fok[fq]) continue; fst[fq] = i;
      if (mark[faces[3 * fq]] === i && mark[faces[3 * fq + 1]] === i && mark[faces[3 * fq + 2]] === i) { sw += om[fq]; sa += ar[fq]; }
    }
    K[i] = sa > 0 ? sw * s * s / sa : NaN;
  }
  return K;
}

function pct(a: number[], q: number) { const i = (a.length - 1) * q, lo = Math.floor(i), f = i - lo; return a[lo] + (a[Math.min(a.length - 1, lo + 1)] - a[lo]) * f; }

export interface OverlayResult {
  /** per-face value in [-1,1] (NaN where grey) */
  faceT: Float32Array;
  /** per-vertex value in [-1,1] (NaN where grey) */
  vertT: Float32Array;
  lo: number; hi: number; // legend half-ranges (value of -1 and +1)
}

export class CurvatureOverlay {
  mode: 'none' | 'length' | 'gauss' = 'none';
  srel = parseFloat(SENSITIVITY_DEFAULT);
  rho = 1;
  private range: [number, number] | null = null;
  private rangeKey = '';
  private key = '';
  private cache: OverlayResult | null = null;
  private cacheFrame: FlowFrame | null = null;

  invalidate() { this.cache = null; this.cacheFrame = null; }
  resetRange() { this.range = null; this.rangeKey = ''; this.invalidate(); }

  private estimate(V: Float32Array, faces: Uint32Array, tips: number[]) { return this.mode === 'length' ? estimateIntr(V, faces, tips, 1) : estimateExtr(V, faces, tips, 1); }

  /** Range from the seed frame: p1/p99 of measured K widened to include rho, plus BUF of the span on each side (fixed for the flow). */
  private fitRange(f0: FlowFrame, tips: number[]) {
    const key = this.mode + '|' + this.rho;
    if (this.range && this.rangeKey === key) return;
    const K = this.estimate(f0.vertices, f0.faces, tips), e: number[] = [];
    for (const k of K) if (k === k) e.push(k);
    if (!e.length) this.range = null;
    else {
      e.sort((a, b) => a - b);
      const lo = Math.min(pct(e, RQ), this.rho), hi = Math.max(pct(e, 1 - RQ), this.rho), sp = Math.max(hi - lo, 1e-3 * Math.abs(this.rho));
      this.range = [lo - BUF * sp, hi + BUF * sp];
    }
    this.rangeKey = key; this.cache = null;
  }

  compute(f: FlowFrame, f0: FlowFrame, tips: number[]): OverlayResult | null {
    if (this.mode === 'none') return null;
    this.fitRange(f0, tips);
    const key = `${this.mode}|${this.srel}|${this.rho}`;
    if (this.cache && this.cacheFrame === f && this.key === key) return this.cache;
    const nv = f.vertices.length / 3, nf = f.faces.length / 3, rho = this.rho;
    const fr = this.srel / 0.5, r = this.range ?? [rho - 0.3 * Math.abs(rho), rho + 0.3 * Math.abs(rho)];
    let dl = Math.max(fr * (rho - r[0]), 1e-3 * Math.abs(rho)), dh = Math.max(fr * (r[1] - rho), 1e-3 * Math.abs(rho));
    const K = this.estimate(f.vertices, f.faces, tips);
    let Kt: Float64Array | null = null;
    if (this.mode === 'length' && f.Kdef && f.Avert) {
      const T = topology(f.faces, nv, tips); Kt = new Float64Array(nv);
      for (let i = 0; i < nv; i++) {
        const bl = T.ball[i]; if (!bl) { Kt[i] = NaN; continue; }
        let sd = 0, sa = 0;
        for (const u of bl) { if (T.ex[u]) continue; sd += f.Kdef[u]; sa += f.Avert[u]; }
        Kt[i] = sa > 0 ? sd / sa : NaN;
      }
      const rr = f.rho ?? rho;
      dl = dh = Math.max(this.srel * Math.abs(rr), 1e-6);
    }
    const polyW = this.mode === 'length' && !!f.polytope; // honest polytope frame: no fair length target
    const vertT = new Float32Array(nv);
    for (let i = 0; i < nv; i++) {
      if (polyW) { vertT[i] = 0; continue; }
      if (K[i] !== K[i] || (Kt && Kt[i] !== Kt[i])) { vertT[i] = NaN; continue; }
      const e = K[i] - (Kt ? Kt[i] : rho);
      vertT[i] = Math.max(-1, Math.min(1, e < 0 ? e / dl : e / dh));
    }
    const faceT = new Float32Array(nf);
    for (let q = 0; q < nf; q++) faceT[q] = (vertT[f.faces[3 * q]] + vertT[f.faces[3 * q + 1]] + vertT[f.faces[3 * q + 2]]) / 3;
    this.cache = { faceT, vertT, lo: dl, hi: dh }; this.cacheFrame = f; this.key = key;
    return this.cache;
  }

  /** Mean of the adjacent faces' values (NaN if all grey): used for wireframe edge colours. */
  edgeColors(f: FlowFrame, res: OverlayResult): { a: number; b: number; rgb: number[] }[] {
    const T = topology(f.faces, f.vertices.length / 3, TOPO ? TOPO.tips : []);
    if (!T.edges) {
      const m = new Map<number, number>(), list: { a: number; b: number; fs: number[] }[] = [], nv = T.nv, nf = f.faces.length / 3;
      for (let q = 0; q < nf; q++) for (let k = 0; k < 3; k++) {
        const a = f.faces[3 * q + k], b = f.faces[3 * q + (k + 1) % 3], key = a < b ? a * nv + b : b * nv + a, e = m.get(key);
        if (e === undefined) { m.set(key, list.length); list.push({ a: Math.min(a, b), b: Math.max(a, b), fs: [q] }); } else list[e].fs.push(q);
      }
      T.edges = list;
    }
    return T.edges.map(E => {
      let s = 0, n = 0;
      for (const q of E.fs) { const t = res.faceT[q]; if (t === t) { s += t; n++; } }
      const c = pickLut(n ? s / n : NaN);
      return { a: E.a, b: E.b, rgb: [Math.round(c[0]), Math.round(c[1]), Math.round(c[2])] };
    });
  }
}
