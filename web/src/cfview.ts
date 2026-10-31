import * as THREE from 'three';
import { magma, RAW_COLOR_SPACE } from './shading';
import type { FlowFrame, FlowMeta, WGrid } from './solver';

// Port of viewer/frontend/conformal_view.py (W3D_JS): height graph of the conformal factor w over the north / south
// hemisphere seen as a stereographic disc, smooth (bilinear sampling of the w raster on a fine mesh clustered around
// the cones, cone model added back analytically so the tips are sharp cusps), magma colours, faint Lambert shading.
// The raster is {nx, ny, values, extent}: row 0 at colatitude extent[2], column 0 at longitude extent[0].

type V3 = [number, number, number];
interface Hem { verts: V3[][]; NI: number; NJ: number }
interface Mk { i: number; u: number; ph: number; v: number }

const unitOf = (p: number[]): V3 => { const d = 1 + p[0] * p[0] + p[1] * p[1]; return [2 * p[0] / d, 2 * p[1] / d, (p[0] * p[0] + p[1] * p[1] - 1) / d]; };
const dirOf = (th: number, ln: number): V3 => [Math.sin(th) * Math.cos(ln), Math.sin(th) * Math.sin(ln), Math.cos(th)];
function cmodel(X: V3, U: V3[], B: number[], fl: number) {
  let s = 0;
  for (let k = 0; k < U.length; k++) { const dx = X[0] - U[k][0], dy = X[1] - U[k][1], dz = X[2] - U[k][2]; s += (B[k] - 1) * Math.log(Math.max(Math.sqrt(dx * dx + dy * dy + dz * dz), fl)); }
  return s;
}
function lines(n: number, T: number, extra: number[]) {
  const a: number[] = [];
  for (let i = 0; i <= n; i++) a.push(T * i / n);
  for (const e of extra) if (e > 0 && e < T) a.push(e);
  a.sort((p, q) => p - q);
  const o = [a[0]];
  for (let i = 1; i < a.length; i++) if (a[i] - o[o.length - 1] > 1e-6) o.push(a[i]);
  return o;
}

const TWO = 2 * Math.PI, HP = Math.PI / 2;

export class CfDiscs {
  private renderers: THREE.WebGLRenderer[];
  private scenes = [new THREE.Scene(), new THREE.Scene()];
  private cam = new THREE.OrthographicCamera(-1.25, 1.25, 1.25, -1.25, -10, 10);
  private meshes: (THREE.Mesh | null)[] = [null, null];
  private betas: number[] = []; private positions: number[][] = [];
  private range: [number, number] | null = null;
  private grid: WGrid | null = null;
  private hems: Hem[] | null = null;
  private dark = false;
  // view state (shared by both plots, same defaults and gestures as the Python page)
  private yaw = 0.6; private elev = 1.1; private zoom = 1;
  // draggable cone markers (shown while no flow is displayed, like the Python page's live preview)
  private mk: Mk[][] = [[], []]; private mkEls: HTMLDivElement[][] = [[], []];
  markersOn = false;
  setMarkers(on: boolean) { this.markersOn = on; this.draw(); }
  /** Called while a marker is dragged (fin = false) and on release (fin = true), with the new stereographic position. */
  onDragCone: ((i: number, pos: [number, number], fin: boolean) => void) | null = null;

  constructor(private canvases: HTMLCanvasElement[]) {
    this.renderers = canvases.map(c => {
      const r = new THREE.WebGLRenderer({ canvas: c, antialias: true, preserveDrawingBuffer: true });
      r.outputColorSpace = RAW_COLOR_SPACE; return r;
    });
    canvases.forEach((c, k) => {
      new ResizeObserver(() => this.resize()).observe(c.parentElement!);
      let drag: { x: number; y: number } | null = null;
      c.style.cursor = 'grab'; c.style.touchAction = 'none';
      let cone = -1;
      const pix = (e: PointerEvent) => { const r = c.getBoundingClientRect(); return [(e.clientX - r.left), (e.clientY - r.top)]; };
      const hit = (e: PointerEvent) => {
        if (!this.markersOn) return -1;
        const [px, py] = pix(e); let best = -1, bd = 18;
        for (const m of this.mk[k]) { const q = this.screenOf(m, k); const d = Math.hypot(q[0] - px, q[1] - py); if (d < bd) { bd = d; best = m.i; } }
        return best;
      };
      c.addEventListener('pointerdown', e => { cone = e.button === 0 ? hit(e) : -1; drag = { x: e.clientX, y: e.clientY }; c.style.cursor = 'grabbing'; c.setPointerCapture(e.pointerId); });
      c.addEventListener('pointermove', e => {
        if (!drag) { c.style.cursor = hit(e) >= 0 ? 'pointer' : 'grab'; return; }
        if (cone >= 0) { // drag the marker: inverse of the plot projection at the marker's own height
          const m = this.mk[k].find(z => z.i === cone); if (!m || !this.range) return;
          const [px, py] = pix(e), q0 = this.screenOf(m, k), s = c.clientWidth, h = 1.25 / this.zoom;
          const se = Math.max(0.3, Math.sin(this.elev)), cy = Math.cos(this.yaw), sy = Math.sin(this.yaw);
          const dx = (px - q0[0]) / s * 2 * h, dy = -(py - q0[1]) / s * 2 * h;
          const x1 = m.u * Math.cos(m.ph) * cy - m.u * Math.sin(m.ph) * sy + dx, y1 = m.u * Math.cos(m.ph) * sy + m.u * Math.sin(m.ph) * cy + dy / se;
          let qx = x1 * cy + y1 * sy, qy = -x1 * sy + y1 * cy; const u = Math.hypot(qx, qy);
          if (u > 0.995) { qx *= 0.995 / u; qy *= 0.995 / u; }
          if (u < 1e-4) { qx = 1e-4; qy = 0; }
          this.onDragCone?.(cone, this.diskToPlane(qx, qy, k), false); return;
        }
        this.yaw += (e.clientX - drag.x) * 0.012;
        this.elev = Math.max(0.1, Math.min(1.5, this.elev + (e.clientY - drag.y) * 0.008));
        drag = { x: e.clientX, y: e.clientY }; this.draw();
      });
      const up = () => { const was = cone; drag = null; cone = -1; c.style.cursor = 'grab'; if (was >= 0) this.onDragCone?.(was, this.positions[was] as [number, number], true); };
      c.addEventListener('pointerup', up); c.addEventListener('pointercancel', up);
      c.addEventListener('wheel', e => { e.preventDefault(); this.zoom = Math.max(0.25, Math.min(12, this.zoom * Math.exp(-e.deltaY * 0.0015))); this.draw(); }, { passive: false });
      c.addEventListener('dblclick', () => { this.zoom = 1; this.draw(); });
    });
    canvases.forEach((c, k) => {
      for (let i = 0; i < 5; i++) { const d = document.createElement('div'); d.className = 'cone'; d.hidden = true; d.style.pointerEvents = 'none'; c.parentElement!.appendChild(d); this.mkEls[k].push(d); }
    });
    this.resize();
  }

  private diskToPlane(qx: number, qy: number, k: number): [number, number] {
    const a = 2 * Math.atan(Math.hypot(qx, qy)), ph = Math.atan2(qy, qx), th = k ? Math.PI - a : a, ln = k ? Math.PI - ph : ph - Math.PI;
    const d = dirOf(th, ln); return [d[0] / (1 - d[2]), d[1] / (1 - d[2])];
  }

  private screenOf(m: Mk, k: number): [number, number] {
    const range = this.range!, lo = range[0], span = (range[1] - range[0]) || 1, s = this.canvases[k].clientWidth, h = 1.25 / this.zoom;
    const cy = Math.cos(this.yaw), sy = Math.sin(this.yaw), ce = Math.cos(this.elev), se = Math.sin(this.elev);
    const px = m.u * Math.cos(m.ph), py = m.u * Math.sin(m.ph), v = (m.v - lo) / span - 0.5, x1 = px * cy - py * sy, y1 = px * sy + py * cy;
    return [(x1 / h + 1) / 2 * s, (1 - (v * 0.8 * ce + y1 * se) / h) / 2 * s];
  }

  /** Live preview before any flow: the analytic cone-model height graph for the given cones. */
  setPreview(betas: number[], positions: number[][]) {
    this.betas = betas.slice(); this.positions = positions.map(p => [...p]);
    const ny = 64, nx = 128;
    this.range = null; this.grid = null; this.hems = this.build({ nx, ny, values: new Float32Array(nx * ny), extent: [-Math.PI, Math.PI, 0, Math.PI] });
    this.draw();
  }

  resize() {
    this.canvases.forEach((c, k) => {
      const p = c.parentElement!, s = Math.max(120, Math.floor(p.clientWidth));
      this.renderers[k].setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
      this.renderers[k].setSize(s, s, false);
      c.style.width = s + 'px'; c.style.height = s + 'px';
    });
    this.draw();
  }

  setTheme(dark: boolean) { this.dark = dark; this.draw(); }

  reset(meta: FlowMeta) {
    this.betas = meta.request.betas.slice(); this.positions = meta.request.positions.map(p => [...p]);
    this.range = null; this.grid = null; this.hems = null; this.markersOn = false;
    this.draw();
  }

  /** The colour/height scale is shared by both plots and grows over the run, so frames stay comparable. */
  setFrame(f: FlowFrame) {
    if (!f.wGrid) return;
    this.grid = f.wGrid; this.hems = this.build(f.wGrid); this.draw();
  }

  private build(g: WGrid): Hem[] {
    const { nx: W, ny: H, values } = g, [lon0, lon1, th0, th1] = g.extent, dLon = lon1 - lon0, dTh = th1 - th0;
    const U = this.positions.map(unitOf), B = this.betas, flG = 0.5 * Math.PI / H, flR = 0.002;
    // smooth remainder at the raster cell centres: subtract the cone model, add it back at every render vertex
    const R = new Float64Array(W * H);
    for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) {
      let v0 = values[y * W + x]; if (!isFinite(v0)) v0 = 0;
      v0 -= cmodel(dirOf(th0 + (y + 0.5) * dTh / H, lon0 + (x + 0.5) * dLon / W), U, B, flG);
      R[y * W + x] = v0;
    }
    const samp = (th: number, ln: number) => {
      const py = Math.max(0, Math.min(H - 1, (th - th0) / dTh * H - 0.5)), px = (ln - lon0) / dLon * W - 0.5;
      const y0 = Math.floor(py), ty = py - y0, y1 = Math.min(H - 1, y0 + 1), x0 = Math.floor(px), tx = px - x0, xa = ((x0 % W) + W) % W, xb = (xa + 1) % W;
      return (R[y0 * W + xa] * (1 - tx) + R[y0 * W + xb] * tx) * (1 - ty) + (R[y1 * W + xa] * (1 - tx) + R[y1 * W + xb] * tx) * ty;
    };
    this.mk = [[], []];
    const offs = [0, 0.0025, 0.007, 0.016, 0.035, 0.07], hem: Hem[] = [];
    let mn = Infinity, mx = -Infinity;
    for (let k = 0; k < 2; k++) { // 0: north, 1: south (viewed from the south pole)
      const eT: number[] = [], eP: number[] = [];
      for (const u of U) {
        const ti = Math.acos(Math.max(-1, Math.min(1, u[2]))), li = Math.atan2(u[1], u[0]), ai = k ? Math.PI - ti : ti, pj = k ? Math.PI - li : li + Math.PI, sn = Math.max(0.05, Math.sin(ti));
        for (const o of offs) { eT.push(ai + o, ai - o); const q = o / sn; for (const z of [pj + q, pj - q]) eP.push(((z % TWO) + TWO) % TWO); }
      }
      const tl = lines(44, HP, eT), pl = lines(192, TWO, eP), NI = tl.length - 1, NJ = pl.length - 1, verts: V3[][] = [];
      for (let i = 0; i <= NI; i++) {
        const row: V3[] = [], a = tl[i], u = Math.tan(a / 2);
        for (let j = 0; j <= NJ; j++) {
          const ph = pl[j], th = k ? Math.PI - a : a, ln = k ? Math.PI - ph : ph - Math.PI;
          const v = samp(th, ln) + cmodel(dirOf(th, ln), U, B, flR);
          if (v < mn) mn = v; if (v > mx) mx = v;
          row.push([u * Math.cos(ph), u * Math.sin(ph), v]);
        }
        verts.push(row);
      }
      hem.push({ verts, NI, NJ });
      U.forEach((u, i) => {
        const ti = Math.acos(Math.max(-1, Math.min(1, u[2]))), li = Math.atan2(u[1], u[0]), a = k ? Math.PI - ti : ti;
        if (a > HP) return;
        const ph = k ? Math.PI - li : li + Math.PI;
        this.mk[k].push({ i, u: Math.tan(a / 2), ph, v: samp(ti, li) + cmodel(dirOf(ti, li), U, B, flR) });
      });
    }
    this.range = this.range ? [Math.min(this.range[0], mn), Math.max(this.range[1], mx)] : [mn, mx];
    return hem;
  }

  private draw() {
    const h = 1.25 / this.zoom;
    this.cam.left = -h; this.cam.right = h; this.cam.top = h; this.cam.bottom = -h; this.cam.updateProjectionMatrix();
    for (let k = 0; k < 2; k++) {
      const r = this.renderers[k];
      r.setClearColor(this.dark ? 0x121417 : 0xf4f4f2, 1);
      const old = this.meshes[k]; if (old) { this.scenes[k].remove(old); old.geometry.dispose(); this.meshes[k] = null; }
      if (this.hems && this.range) { const m = this.mesh(this.hems[k], this.range); this.meshes[k] = m; this.scenes[k].add(m); }
      r.render(this.scenes[k], this.cam);
      this.mkEls[k].forEach((el, i) => {
        const m = this.markersOn && this.range ? this.mk[k].find(z => z.i === i) : undefined;
        el.hidden = !m; if (!m) return;
        el.textContent = String(i + 1);
        const q = this.screenOf(m, k); el.style.transform = `translate(${q[0]}px, ${q[1]}px) translate(-50%,-50%)`;
      });
    }
  }

  private mesh(H: Hem, range: [number, number]) {
    const lo = range[0], span = (range[1] - range[0]) || 1;
    const cy = Math.cos(this.yaw), sy = Math.sin(this.yaw), ce = Math.cos(this.elev), se = Math.sin(this.elev);
    // P = [screen x, screen y, depth (larger = farther, as in the Python painter)]
    const P: V3[][] = H.verts.map(r => r.map(p => {
      const v = (p[2] - lo) / span - 0.5, x1 = p[0] * cy - p[1] * sy, y1 = p[0] * sy + p[1] * cy;
      return [x1, v * 0.8 * ce + y1 * se, y1 * ce - v * 0.8 * se] as V3;
    }));
    const nq = H.NI * H.NJ, pos = new Float32Array(nq * 18), col = new Float32Array(nq * 18);
    let o = 0;
    for (let i = 0; i < H.NI; i++) for (let j = 0; j < H.NJ; j++) {
      const a = P[i][j], b = P[i][j + 1], c = P[i + 1][j + 1], d = P[i + 1][j];
      const ux = c[0] - a[0], uy = c[1] - a[1], uz = c[2] - a[2], wx = d[0] - b[0], wy = d[1] - b[1], wz = d[2] - b[2];
      const nx = uy * wz - uz * wy, ny = uz * wx - ux * wz, nz = ux * wy - uy * wx, nl = Math.hypot(nx, ny, nz) || 1;
      const sh = Math.abs((-0.3 * nx + 0.4 * ny - 0.85 * nz) / nl);
      const hv = (H.verts[i][j][2] + H.verts[i][j + 1][2] + H.verts[i + 1][j + 1][2] + H.verts[i + 1][j][2]) / 4;
      const cc = magma((hv - lo) / span), kk = 0.78 + 0.3 * sh, r = Math.min(1, cc[0] * kk), g = Math.min(1, cc[1] * kk), bl = Math.min(1, cc[2] * kk);
      for (const p of [a, b, c, a, c, d]) { pos[o] = p[0]; pos[o + 1] = p[1]; pos[o + 2] = -p[2]; col[o] = r; col[o + 1] = g; col[o + 2] = bl; o += 3; }
    }
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    g.setAttribute('color', new THREE.BufferAttribute(col, 3));
    return new THREE.Mesh(g, new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide }));
  }
}
