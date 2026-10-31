import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { BLUES_LUT, shadeVal, diverge } from './shading';
import type { FlowFrame, FlowMeta } from './solver';

export type Shading = 'reflective' | 'plain' | 'wireframe';
export type Overlay = 'none' | 'length' | 'gauss';

export class View3D {
  private renderer: THREE.WebGLRenderer;
  private scene = new THREE.Scene();
  private camera = new THREE.PerspectiveCamera(32, 1, 0.01, 100);
  private controls: OrbitControls;
  private mesh = new THREE.Mesh(new THREE.BufferGeometry(), new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide }));
  private wire = new THREE.LineSegments(new THREE.BufferGeometry(), new THREE.LineBasicMaterial({ transparent: true, opacity: 0.15 }));
  private frame: FlowFrame | null = null;
  private meta: FlowMeta | null = null;
  private center = new THREE.Vector3();
  private fitted = false;
  private coneEls: HTMLElement[] = [];
  shading: Shading = 'reflective';
  overlay: Overlay = 'none';
  legendRange = 0;
  onLegend: (() => void) | null = null;

  constructor(private canvas: HTMLCanvasElement, private labels: HTMLElement) {
    this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, preserveDrawingBuffer: true });
    this.scene.add(this.mesh, this.wire);
    this.camera.position.set(0, 0, 4.2);
    this.controls = new OrbitControls(this.camera, canvas);
    this.controls.enableDamping = false;
    this.controls.addEventListener('change', () => this.render());
    new ResizeObserver(() => this.resize()).observe(canvas.parentElement!);
    this.resize();
  }

  setTheme(dark: boolean) {
    (this.wire.material as THREE.LineBasicMaterial).color.set(dark ? 0xffffff : 0x000000);
    this.renderer.setClearColor(dark ? 0x121417 : 0xf4f4f2, 1);
    this.render();
  }

  resize() {
    const p = this.canvas.parentElement!, s = Math.max(200, Math.floor(p.clientWidth));
    this.renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
    this.renderer.setSize(s, s, false);
    this.canvas.style.width = s + 'px'; this.canvas.style.height = s + 'px';
    this.camera.aspect = 1; this.camera.updateProjectionMatrix();
    this.render();
  }

  reset(meta: FlowMeta) {
    this.meta = meta; this.frame = null; this.fitted = false;
    this.labels.textContent = '';
    this.coneEls = meta.coneIds.map((_, i) => {
      const d = document.createElement('div'); d.className = 'cone'; d.textContent = String(i + 1);
      this.labels.appendChild(d); return d;
    });
    this.camera.position.set(0, 0, 4.2); this.controls.target.set(0, 0, 0); this.controls.update();
  }

  setFrame(f: FlowFrame) {
    this.frame = f;
    if (!this.fitted) {
      const v = f.vertices, nv = v.length / 3;
      let cx = 0, cy = 0, cz = 0;
      for (let i = 0; i < nv; i++) { cx += v[3 * i]; cy += v[3 * i + 1]; cz += v[3 * i + 2]; }
      this.center.set(cx / nv, cy / nv, cz / nv); this.fitted = true;
      let R = 0;
      for (let i = 0; i < nv; i++) R = Math.max(R, Math.hypot(v[3 * i] - this.center.x, v[3 * i + 1] - this.center.y, v[3 * i + 2] - this.center.z));
      const dist = (R * 1.12) / Math.sin((this.camera.fov * Math.PI) / 360);
      this.camera.position.set(0, 0, dist); this.camera.far = dist * 10; this.camera.updateProjectionMatrix();
      this.controls.target.set(0, 0, 0); this.controls.update();
    }
    this.render();
  }

  render() {
    const f = this.frame;
    if (!f) { this.renderer.render(this.scene, this.camera); return; }
    this.camera.updateMatrixWorld();
    const V = f.vertices, F = f.faces, nF = F.length / 3, c0 = this.center;
    const pos = new Float32Array(nF * 9), col = new Float32Array(nF * 9);
    const view = this.camera.matrixWorldInverse.elements; // column-major
    const faceVal = this.overlayValues(f);
    const tmp = [0, 0, 0];
    for (let k = 0; k < nF; k++) {
      const a = F[3 * k], b = F[3 * k + 1], c = F[3 * k + 2];
      const pa = [V[3 * a] - c0.x, V[3 * a + 1] - c0.y, V[3 * a + 2] - c0.z];
      const pb = [V[3 * b] - c0.x, V[3 * b + 1] - c0.y, V[3 * b + 2] - c0.z];
      const pc = [V[3 * c] - c0.x, V[3 * c + 1] - c0.y, V[3 * c + 2] - c0.z];
      pos.set(pa, 9 * k); pos.set(pb, 9 * k + 3); pos.set(pc, 9 * k + 6);
      const ux = pb[0] - pa[0], uy = pb[1] - pa[1], uz = pb[2] - pa[2], wx = pc[0] - pa[0], wy = pc[1] - pa[1], wz = pc[2] - pa[2];
      const nx = uy * wz - uz * wy, ny = uz * wx - ux * wz, nz = ux * wy - uy * wx;
      // normal into view space (rotation part of the view matrix); two-sided
      let vx = view[0] * nx + view[4] * ny + view[8] * nz, vy = view[1] * nx + view[5] * ny + view[9] * nz, vz = view[2] * nx + view[6] * ny + view[10] * nz;
      if (vz < 0) { vx = -vx; vy = -vy; vz = -vz; }
      let r: number, g: number, bl: number;
      if (this.overlay !== 'none') {
        diverge(faceVal[k], tmp);
        const s = this.shading === 'plain' ? 1 : 0.55 + 0.45 * shadeVal(vx, vy, vz);
        r = tmp[0] * s; g = tmp[1] * s; bl = tmp[2] * s;
      } else if (this.shading === 'plain') { r = 0.45; g = 0.62; bl = 0.82; }
      else { const i = Math.round(255 * shadeVal(vx, vy, vz)) * 3; r = BLUES_LUT[i]; g = BLUES_LUT[i + 1]; bl = BLUES_LUT[i + 2]; }
      for (let q = 0; q < 3; q++) { col[9 * k + 3 * q] = r; col[9 * k + 3 * q + 1] = g; col[9 * k + 3 * q + 2] = bl; }
    }
    const g = this.mesh.geometry;
    g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    g.setAttribute('color', new THREE.BufferAttribute(col, 3));
    this.mesh.visible = this.shading !== 'wireframe' || this.overlay !== 'none';
    this.wire.geometry.dispose();
    this.wire.geometry = new THREE.WireframeGeometry(g);
    (this.wire.material as THREE.LineBasicMaterial).opacity = this.shading === 'wireframe' ? 0.8 : 0.15;
    this.renderer.render(this.scene, this.camera);
    this.placeCones(f);
  }

  /** Per-face diverging value in [-1,1] for the overlay; sets legendRange (the value of +-1). */
  private overlayValues(f: FlowFrame): Float32Array {
    const nF = f.faces.length / 3, out = new Float32Array(nF);
    if (this.overlay === 'none' || !this.meta) { this.legendRange = 0; return out; }
    const nv = f.vertices.length / 3, acc = new Float32Array(nv), cnt = new Float32Array(nv);
    if (this.overlay === 'length') {
      const E = this.meta.edges, V = f.vertices, L = f.L;
      if (!L) return out;
      for (let e = 0; e < E.length / 2; e++) {
        const a = E[2 * e], b = E[2 * e + 1];
        const d = Math.hypot(V[3 * a] - V[3 * b], V[3 * a + 1] - V[3 * b + 1], V[3 * a + 2] - V[3 * b + 2]);
        const x = d / L[e];
        acc[a] += x; acc[b] += x; cnt[a]++; cnt[b]++;
      }
      // the 3D shape is only meaningful up to the global area scale: normalise by the median ratio
      const r = Array.from(acc, (s, i) => s / (cnt[i] || 1)).sort((p, q) => p - q), med = r[r.length >> 1] || 1;
      for (let i = 0; i < nv; i++) acc[i] = acc[i] / (cnt[i] || 1) / med - 1;
    } else {
      const cones = new Set(this.meta.coneIds);
      const sm = Array.from(f.K).filter((_, i) => !cones.has(i)).sort((p, q) => p - q);
      const med = sm[sm.length >> 1] || 0;
      for (let i = 0; i < nv; i++) acc[i] = f.K[i] - med;
    }
    const sorted = Array.from(acc, Math.abs).sort((p, q) => p - q);
    const m = Math.max(1e-6, sorted[Math.floor(0.98 * (nv - 1))]);
    this.legendRange = m;
    const F = f.faces;
    for (let k = 0; k < nF; k++) out[k] = (acc[F[3 * k]] + acc[F[3 * k + 1]] + acc[F[3 * k + 2]]) / 3 / m;
    this.onLegend?.();
    return out;
  }

  private placeCones(f: FlowFrame) {
    if (!this.meta) return;
    const w = this.canvas.clientWidth, h = this.canvas.clientHeight, v = new THREE.Vector3();
    this.meta.coneIds.forEach((id, i) => {
      v.set(f.vertices[3 * id] - this.center.x, f.vertices[3 * id + 1] - this.center.y, f.vertices[3 * id + 2] - this.center.z).project(this.camera);
      const el = this.coneEls[i]; if (!el) return;
      el.style.transform = `translate(${(v.x * 0.5 + 0.5) * w}px, ${(-v.y * 0.5 + 0.5) * h}px) translate(-50%,-50%)`;
    });
  }
}
