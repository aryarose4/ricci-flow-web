// Port of viewer/frontend/flow_playback.py (FP_JS): display-side "Smooth Playback".
// Between consecutive TRUE frames A, B the page shows ticks of a C1 cubic Hermite interpolant through P, A, B, N
// (rigidly de-rotated; time reparametrised by jump size; Bessel tangents with a monotone limiter). The solver output
// is untouched. Same constants as Python.

export const TARGET = 0.004, MAXSUB = 40, TICK_MS = 40;

type V3 = number[];
type M3 = number[][];
type Cloud = V3[];
interface Prep { R: M3; ca: V3; cb: V3; Bp: Cloud; k: V3; a: number; J: number; n: number }
interface Plan { n: number; m: number; yA: Float64Array; yB: Float64Array; m0: Float64Array; m1: Float64Array }

/** A true frame as the player sees it (V = vertex triples); `fp` is set on synthesized sub-frames. */
export interface PFrame { V: Cloud; fp?: { from: PFrame; to: PFrame; i: number; n: number }; fpxA?: Cloud; fpxP?: Prep }
export interface Host {
  frames(): PFrame[];
  set(fr: PFrame): void;        // display this frame (true or synthesized)
  setSlider(k: number): void;   // move the slider to true frame k
  follow(): boolean;            // keep the slider following playback
  smooth(): boolean;
  setPlayLabel(text: string): void;
  sliderValue(): number;
}

const ctr = (V: Cloud): V3 => { const c = [0, 0, 0], n = V.length; for (const v of V) { c[0] += v[0]; c[1] += v[1]; c[2] += v[2]; } return [c[0] / n, c[1] / n, c[2] / n]; };
function diag(V: Cloud) {
  const lo = [1e300, 1e300, 1e300], hi = [-1e300, -1e300, -1e300];
  for (const v of V) for (let k = 0; k < 3; k++) { if (v[k] < lo[k]) lo[k] = v[k]; if (v[k] > hi[k]) hi[k] = v[k]; }
  return Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);
}
const det3 = (M: M3) => M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
const cof = (M: M3): M3 => [[M[1][1] * M[2][2] - M[1][2] * M[2][1], M[1][2] * M[2][0] - M[1][0] * M[2][2], M[1][0] * M[2][1] - M[1][1] * M[2][0]],
  [M[0][2] * M[2][1] - M[0][1] * M[2][2], M[0][0] * M[2][2] - M[0][2] * M[2][0], M[0][1] * M[2][0] - M[0][0] * M[2][1]],
  [M[0][1] * M[1][2] - M[0][2] * M[1][1], M[0][2] * M[1][0] - M[0][0] * M[1][2], M[0][0] * M[1][1] - M[0][1] * M[1][0]]];
const I3: M3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

function rigid(A: Cloud, B: Cloud, ca: V3, cb: V3): M3 {
  const M: M3 = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < A.length; i++) for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) M[r][c] += (B[i][r] - cb[r]) * (A[i][c] - ca[c]);
  if (det3(M) <= 1e-18) return I3;
  let nf = 0; for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) nf += M[r][c] * M[r][c];
  nf = Math.sqrt(nf);
  let Q: M3 = M.map(row => row.map(x => x / nf));
  for (let it = 0; it < 60; it++) {
    const d = det3(Q); if (Math.abs(d) < 1e-14) return I3;
    const C = cof(Q); let mx = 0; const Qn: M3 = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) { Qn[r][c] = 0.5 * (Q[r][c] + C[r][c] / d); mx = Math.max(mx, Math.abs(Qn[r][c] - Q[r][c])); }
    Q = Qn; if (mx < 1e-13) break;
  }
  if (det3(Q) < 0) return I3;
  const ang = Math.acos(Math.max(-1, Math.min(1, (Q[0][0] + Q[1][1] + Q[2][2] - 1) / 2)));
  return ang > 2.5 ? I3 : Q;
}
function axang(R: M3) {
  const ang = Math.acos(Math.max(-1, Math.min(1, (R[0][0] + R[1][1] + R[2][2] - 1) / 2)));
  if (ang < 1e-9) return { k: [0, 0, 1], a: 0 };
  const k = [R[2][1] - R[1][2], R[0][2] - R[2][0], R[1][0] - R[0][1]], n = Math.hypot(k[0], k[1], k[2]);
  return { k: [k[0] / n, k[1] / n, k[2] / n], a: ang };
}
function prep(A: Cloud, Bf: PFrame): Prep {
  const B = Bf.V;
  if (Bf.fpxA === A && Bf.fpxP) return Bf.fpxP;
  const ca = ctr(A), cb = ctr(B), R = rigid(A, B, ca, cb), aa = axang(R), n = A.length, Bp: Cloud = new Array(n);
  let resid = 0, rmax = 0;
  for (let i = 0; i < n; i++) {
    const x = B[i][0] - cb[0], y = B[i][1] - cb[1], z = B[i][2] - cb[2];
    Bp[i] = [R[0][0] * x + R[1][0] * y + R[2][0] * z, R[0][1] * x + R[1][1] * y + R[2][1] * z, R[0][2] * x + R[1][2] * y + R[2][2] * z];
    let r = Math.hypot(Bp[i][0] - (A[i][0] - ca[0]), Bp[i][1] - (A[i][1] - ca[1]), Bp[i][2] - (A[i][2] - ca[2])); if (r > resid) resid = r;
    r = Math.hypot(A[i][0] - ca[0], A[i][1] - ca[1], A[i][2] - ca[2]); if (r > rmax) rmax = r;
  }
  const J = (resid + aa.a * rmax + Math.hypot(cb[0] - ca[0], cb[1] - ca[1], cb[2] - ca[2])) / Math.max(diag(A), 1e-12);
  const p: Prep = { R, ca, cb, Bp, k: aa.k, a: aa.a, J, n: Math.min(Math.max(Math.ceil(J / TARGET), 1), MAXSUB) };
  Bf.fpxA = A; Bf.fpxP = p; return p;
}
function quat(R: M3): number[] {
  const t = R[0][0] + R[1][1] + R[2][2]; let s: number, q: number[];
  if (t > 0) { s = 2 * Math.sqrt(t + 1); q = [s / 4, (R[2][1] - R[1][2]) / s, (R[0][2] - R[2][0]) / s, (R[1][0] - R[0][1]) / s]; }
  else if (R[0][0] >= R[1][1] && R[0][0] >= R[2][2]) { s = 2 * Math.sqrt(1 + R[0][0] - R[1][1] - R[2][2]); q = [(R[2][1] - R[1][2]) / s, s / 4, (R[0][1] + R[1][0]) / s, (R[0][2] + R[2][0]) / s]; }
  else if (R[1][1] >= R[2][2]) { s = 2 * Math.sqrt(1 + R[1][1] - R[0][0] - R[2][2]); q = [(R[0][2] - R[2][0]) / s, (R[0][1] + R[1][0]) / s, s / 4, (R[1][2] + R[2][1]) / s]; }
  else { s = 2 * Math.sqrt(1 + R[2][2] - R[0][0] - R[1][1]); q = [(R[1][0] - R[0][1]) / s, (R[0][2] + R[2][0]) / s, (R[1][2] + R[2][1]) / s, s / 4]; }
  const nn = Math.hypot(q[0], q[1], q[2], q[3]); return [q[0] / nn, q[1] / nn, q[2] / nn, q[3] / nn];
}
function qmat(q: number[]): M3 {
  const nn = Math.hypot(q[0], q[1], q[2], q[3]), w = q[0] / nn, x = q[1] / nn, y = q[2] / nn, z = q[3] / nn;
  return [[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)], [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)], [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]];
}
const qsign = (q: number[], r: number[]) => (q[0] * r[0] + q[1] * r[1] + q[2] * r[2] + q[3] * r[3] < 0 ? q.map(x => -x) : q);
// knot vector layout: [shape 3m | quaternion 4 | centroid 3]
function knot(S: Cloud, q: number[], c: V3) {
  const m = S.length, y = new Float64Array(3 * m + 7);
  for (let i = 0; i < m; i++) { y[3 * i] = S[i][0]; y[3 * i + 1] = S[i][1]; y[3 * i + 2] = S[i][2]; }
  for (let i = 0; i < 4; i++) y[3 * m + i] = q[i];
  for (let i = 0; i < 3; i++) y[3 * m + 4 + i] = c[i];
  return y;
}
function tang(d1: Float64Array, d2: Float64Array, h1: number, h2: number, m: number) {
  const out = new Float64Array(d1.length), mm = [0, 0, 0, 0];
  const groups: [number, number, number][] = [[3, m, 0], [4, 1, 3 * m], [3, 1, 3 * m + 4]];
  for (const [dim, cnt, o] of groups) for (let a = 0; a < cnt; a++) {
    const b = o + a * dim; let dot = 0, n1 = 0, n2 = 0, nm = 0;
    for (let k = 0; k < dim; k++) { mm[k] = (h2 * d1[b + k] + h1 * d2[b + k]) / (h1 + h2); dot += d1[b + k] * d2[b + k]; n1 += d1[b + k] ** 2; n2 += d2[b + k] ** 2; nm += mm[k] ** 2; }
    let sc = 0; if (dot > 0) sc = Math.min(1, 3 * Math.sqrt(Math.min(n1, n2)) / Math.max(Math.sqrt(nm), 1e-300));
    for (let k = 0; k < dim; k++) out[b + k] = mm[k] * sc;
  }
  return out;
}
const slope = (y1: Float64Array, y0: Float64Array, h: number) => { const o = new Float64Array(y1.length); for (let i = 0; i < o.length; i++) o[i] = (y1[i] - y0[i]) / h; return o; };

function seg(P: PFrame | null, Af: PFrame, Bf: PFrame, Nf: PFrame | null): Plan {
  const A = Af.V, B = Bf.V, pa = prep(A, Bf), m = A.length, n = pa.n, ca = pa.ca, cb = pa.cb, R = pa.R, S: Cloud = new Array(m);
  for (let i = 0; i < m; i++) S[i] = [A[i][0] - ca[0], A[i][1] - ca[1], A[i][2] - ca[2]];
  const qI = [1, 0, 0, 0], qB = qsign(quat(R), qI), yA = knot(S, qI, ca), yB = knot(pa.Bp, qB, cb), dA = slope(yB, yA, n);
  let m0 = new Float64Array(yA.length), m1 = new Float64Array(yA.length);
  if (P) {
    const p1 = prep(P.V, Af), R1 = p1.R, cp = p1.ca, Pt: Cloud = new Array(m);
    for (let i = 0; i < m; i++) {
      const x = P.V[i][0] - cp[0], y = P.V[i][1] - cp[1], z = P.V[i][2] - cp[2];
      Pt[i] = [R1[0][0] * x + R1[0][1] * y + R1[0][2] * z, R1[1][0] * x + R1[1][1] * y + R1[1][2] * z, R1[2][0] * x + R1[2][1] * y + R1[2][2] * z];
    }
    const R1T: M3 = [[R1[0][0], R1[1][0], R1[2][0]], [R1[0][1], R1[1][1], R1[2][1]], [R1[0][2], R1[1][2], R1[2][2]]];
    const yP = knot(Pt, qsign(quat(R1T), qI), cp);
    m0 = tang(slope(yA, yP, p1.n), dA, p1.n, n, m);
  }
  if (Nf) {
    const p2 = prep(B, Nf), R2 = p2.R, Nt: Cloud = new Array(m);
    for (let i = 0; i < m; i++) {
      const x = p2.Bp[i][0], y = p2.Bp[i][1], z = p2.Bp[i][2];
      Nt[i] = [R[0][0] * x + R[1][0] * y + R[2][0] * z, R[0][1] * x + R[1][1] * y + R[2][1] * z, R[0][2] * x + R[1][2] * y + R[2][2] * z];
    }
    const RR: M3 = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) for (let k = 0; k < 3; k++) RR[r][c] += R2[r][k] * R[k][c];
    const yN = knot(Nt, qsign(quat(RR), qB), p2.cb);
    m1 = tang(dA, slope(yN, yB, p2.n), n, p2.n, m);
  }
  return { n, m, yA, yB, m0, m1 };
}
function segEval(s: Plan, t: number): Cloud {
  const n = s.n, m = s.m, t2 = t * t, t3 = t2 * t, h00 = 2 * t3 - 3 * t2 + 1, h10 = t3 - 2 * t2 + t, h01 = -2 * t3 + 3 * t2, h11 = t3 - t2, L = s.yA.length, y = new Float64Array(L);
  for (let i = 0; i < L; i++) y[i] = h00 * s.yA[i] + h10 * n * s.m0[i] + h01 * s.yB[i] + h11 * n * s.m1[i];
  const Q = qmat([y[3 * m], y[3 * m + 1], y[3 * m + 2], y[3 * m + 3]]), cx = y[3 * m + 4], cy = y[3 * m + 5], cz = y[3 * m + 6], out: Cloud = new Array(m);
  for (let i = 0; i < m; i++) {
    const x = y[3 * i], yy = y[3 * i + 1], z = y[3 * i + 2];
    out[i] = [Q[0][0] * x + Q[0][1] * yy + Q[0][2] * z + cx, Q[1][0] * x + Q[1][1] * yy + Q[1][2] * z + cy, Q[2][0] * x + Q[2][1] * yy + Q[2][2] * z + cz];
  }
  return out;
}

export class Player {
  private lastReal: PFrame | null = null; private prevReal: PFrame | null = null;
  private queue: PFrame[] = []; private cur: { A: PFrame; B: PFrame; n: number; i: number; pl: Plan | null } | null = null;
  private timer: number | null = null; playing = false;
  constructor(private H: Host) {}

  private idx = (fr: PFrame) => this.H.frames().indexOf(fr);
  private anchor(fr: PFrame) { const F = this.H.frames(), k = F.indexOf(fr); this.lastReal = fr; this.prevReal = k > 0 ? F[k - 1] : null; }
  private synth(B: PFrame, V: Cloud, A: PFrame, i: number, n: number): PFrame { return { ...B, V, fp: { from: A, to: B, i, n } }; }
  clear() { if (this.timer != null) { clearInterval(this.timer); this.timer = null; } this.cur = null; this.queue = []; this.playing = false; this.H.setPlayLabel('Play'); }
  busy() { return this.timer != null; }

  private tick() {
    const H = this.H;
    if (!this.cur) {
      if (!this.queue.length) { this.clear(); if (this.lastReal && H.follow()) H.setSlider(this.idx(this.lastReal)); return; }
      const B = this.queue.shift()!, A = this.lastReal!; let n = 1, pl: Plan | null = null;
      if (H.smooth() && A && A.V.length === B.V.length && (this.playing || this.queue.length <= 4)) { // backlog: one tick per frame
        const P = this.prevReal && this.prevReal.V.length === A.V.length ? this.prevReal : null;
        const N = this.queue.length && this.queue[0].V.length === B.V.length ? this.queue[0] : null;
        pl = seg(P, A, B, N); n = pl.n;
      }
      this.cur = { A, B, n, i: 0, pl };
    }
    const c = this.cur; c.i++;
    const ia = this.idx(c.A), ib = this.idx(c.B);
    if (c.i >= c.n || !c.pl) { this.anchor(c.B); H.set(c.B); if (H.follow()) H.setSlider(ib); this.cur = null; return; }
    H.set(this.synth(c.B, segEval(c.pl, c.i / c.n), c.A, c.i, c.n));
    if (H.follow()) H.setSlider(2 * c.i >= c.n ? ib : ia); // nearest TRUE frame, monotone
  }
  private pump() { if (this.timer == null) this.timer = window.setInterval(() => this.tick(), TICK_MS); }

  /** A newly received TRUE frame (live stream). */
  show(fr: PFrame) {
    if (!this.H.smooth() || !this.lastReal || (!this.cur && !this.queue.length && this.lastReal.V.length !== fr.V.length)) { this.anchor(fr); this.H.set(fr); return; }
    if (!this.cur && !this.queue.length && this.timer == null && prep(this.lastReal.V, fr).n <= 1) { this.anchor(fr); this.H.set(fr); return; }
    this.queue.push(fr); this.pump();
  }
  stop(snap: boolean) {
    const had = this.timer != null, q = this.queue.slice();
    if (this.cur) q.unshift(this.cur.B);
    this.clear();
    if (snap && had) { const f = q.length ? q[q.length - 1] : null; if (f) { this.anchor(f); this.H.set(f); } }
    if (!snap) { this.lastReal = null; this.prevReal = null; }
  }
  play() {
    const H = this.H, fr = H.frames(); if (fr.length < 2) return;
    if (this.playing) { // pause: rest on a TRUE frame (nearest), never mid-motion
      const c = this.cur; this.clear();
      if (c) { const f = 2 * c.i >= c.n ? c.B : c.A; this.anchor(f); H.set(f); H.setSlider(this.idx(f)); }
      return;
    }
    let k = H.sliderValue(); if (!(k >= 0) || k >= fr.length - 1) k = 0;
    this.clear(); this.anchor(fr[k]); H.set(fr[k]); H.setSlider(k); this.queue = fr.slice(k + 1); this.playing = true; H.setPlayLabel('Pause'); this.pump();
  }
  /** The user moved the slider: drop any animation and re-anchor on that true frame. */
  seek(f: PFrame | undefined) { this.clear(); if (f) this.anchor(f); else { this.lastReal = null; this.prevReal = null; } }
}
