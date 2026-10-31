import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { BLUES_LUT, shadeVal, pickLut, RAW_COLOR_SPACE } from './shading';
import { CurvatureOverlay } from './curvature';
import type { FlowFrame, FlowMeta, PreviewData } from './solver';

export type Shading = 'reflective' | 'plain' | 'wireframe';
export type Overlay = 'none' | 'length' | 'gauss';

export class View3D {
  private renderer: THREE.WebGLRenderer;
  private scene = new THREE.Scene();
  // The Python page draws an orthographic projection (flat painter); match it.
  private camera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0.01, 100);
  private controls: OrbitControls;
  private mesh = new THREE.Mesh(new THREE.BufferGeometry(), new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.FrontSide }));
  private wire = new THREE.LineSegments(new THREE.BufferGeometry(), new THREE.LineBasicMaterial({ transparent: true, opacity: 0.15 }));
  private wfLines = new THREE.LineSegments(new THREE.BufferGeometry(), new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, depthTest: false }));
  private frame: FlowFrame | null = null;
  private frame0: FlowFrame | null = null;
  private meta: FlowMeta | null = null;
  private center = new THREE.Vector3(); // fixed at the origin, like the Python page
  private radius = 1;
  private fitted = false;
  /** Python page's initial orientation: rotX(-0.45) . rotY(pi) (it looks at the surface from the far side), applied to the data. */
  private Rm = (() => { const c = Math.cos(0.45), sn = Math.sin(-0.45); return [[-1, 0, 0], [0, c, sn], [0, sn, -c]]; })();
  private dark = false;
  private coneEls: HTMLElement[] = [];
  readonly curv = new CurvatureOverlay();
  shading: Shading = 'reflective';
  overlay: Overlay = 'none';
  legendRange: [number, number] = [0, 0];
  onLegend: (() => void) | null = null;

  constructor(private canvas: HTMLCanvasElement, private labels: HTMLElement) {
    this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, preserveDrawingBuffer: true });
    this.renderer.outputColorSpace = RAW_COLOR_SPACE;
    this.scene.add(this.mesh, this.wire, this.wfLines);
    this.camera.position.set(0, 0, 10);
    this.controls = new OrbitControls(this.camera, canvas);
    this.controls.enableDamping = false;
    this.controls.addEventListener('change', () => this.render());
    new ResizeObserver(() => this.resize()).observe(canvas.parentElement!);
    this.resize();
  }

  setTheme(dark: boolean) {
    this.dark = dark;
    this.renderer.setClearColor(dark ? 0x121417 : 0xf4f4f2, 1);
    this.render();
  }

  setOverlay(o: Overlay) { this.overlay = o; this.curv.mode = o; this.curv.invalidate(); this.render(); }
  setSensitivity(v: number) { this.curv.srel = v; this.curv.invalidate(); this.render(); }

  resize() {
    const p = this.canvas.parentElement!, s = Math.max(200, Math.floor(p.clientWidth));
    this.renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
    this.renderer.setSize(s, s, false);
    this.canvas.style.width = s + 'px'; this.canvas.style.height = s + 'px';
    this.render();
  }

  reset(meta: FlowMeta) {
    this.previewMode = false; this.keepOrbit = false; this.meta = meta; this.frame = null; this.frame0 = null; this.fitted = false;
    this.curv.rho = meta.rho; this.curv.resetRange();
    this.labels.textContent = '';
    this.coneEls = meta.coneIds.map((_, i) => {
      const d = document.createElement('div'); d.className = 'cone'; d.textContent = String(i + 1);
      this.labels.appendChild(d); return d;
    });
    this.camera.position.set(0, 0, 10); this.camera.zoom = 1; this.camera.updateProjectionMatrix();
    this.controls.target.set(0, 0, 0); this.controls.update();
  }

  /** Live seed preview (cap polyhedron): same framing, shading and cone markers as a flow frame. Keeps the user's orbit. */
  setPreview(p: PreviewData, request: FlowMeta['request']) {
    const first = !this.previewMode;
    const nv = p.vertices.length / 3, z = new Float32Array(nv);
    this.meta = { request, rho: p.rho, coneIds: p.coneIds, sphereVertices: p.vertices, sphereFaces: p.faces, edges: p.edges, nFrames: null };
    this.curv.rho = p.rho;
    if (first || this.coneEls.length !== p.coneIds.length) {
      this.labels.textContent = '';
      this.coneEls = p.coneIds.map((_, i) => { const d = document.createElement('div'); d.className = 'cone'; d.textContent = String(i + 1); this.labels.appendChild(d); return d; });
    }
    if (first) { this.curv.resetRange(); this.camera.position.set(0, 0, 10); this.camera.zoom = 1; this.controls.target.set(0, 0, 0); }
    this.previewMode = true;
    const f: FlowFrame = { t: 0, faces: p.faces, vertices: p.vertices, u: z, K: z, gbError: 0, quality: null, progress: 0 };
    this.fitted = false; this.keepOrbit = !first;
    this.setFrame(f, f);
  }
  private previewMode = false;
  private keepOrbit = false;

  /** Back to the empty state (no flow shown). */
  clear() {
    this.previewMode = false; this.frame = null; this.frame0 = null; this.meta = null; this.fitted = false; this.labels.textContent = ''; this.coneEls = [];
    this.mesh.visible = false; this.wire.visible = false; this.wfLines.visible = false; this.onLegend?.();
    this.renderer.render(this.scene, this.camera);
  }

  setFrame(f: FlowFrame, f0?: FlowFrame) {
    this.frame = f; this.frame0 = f0 ?? this.frame0 ?? f;
    if (!this.fitted) {
      // Python page: raw origin as centre, scale = 0.38 * canvas / min(mr, 2.6) per unit of the frame-0 median radius
      const v = f.vertices, nv = v.length / 3, rr: number[] = [];
      for (let i = 0; i < nv; i++) rr.push(Math.hypot(v[3 * i], v[3 * i + 1], v[3 * i + 2]));
      rr.sort((p, q) => p - q);
      const med = rr[nv >> 1] || 1;
      let mr = 1;
      for (const id of this.meta?.coneIds ?? []) mr = Math.max(mr, rr.length ? Math.hypot(v[3 * id], v[3 * id + 1], v[3 * id + 2]) / med : 1);
      this.radius = med; this.fitted = true;
      const h = 0.5 / 0.38 * Math.min(Math.max(1, mr * 0.8), 2.6) * med; // seed tips (not in the fixtures) are nearer than the final cones
      this.camera.left = -h; this.camera.right = h; this.camera.top = h; this.camera.bottom = -h;
      this.camera.far = med * 40;
      if (this.keepOrbit) { const d = this.camera.position.length() || 1; this.camera.position.multiplyScalar(med * 12 / d); this.camera.updateProjectionMatrix(); }
      else { this.camera.position.set(0, 0, med * 12); this.camera.zoom = 1; this.camera.updateProjectionMatrix(); this.controls.target.set(0, 0, 0); this.controls.update(); }
    }
    this.render();
  }

  render() {
    const f = this.frame;
    if (!f || !this.meta) { this.renderer.render(this.scene, this.camera); return; }
    this.camera.updateMatrixWorld();
    const V = f.vertices, F = f.faces, nF = F.length / 3;
    const view = this.camera.matrixWorldInverse.elements; // column-major
    const res = this.curv.compute(f, this.frame0 ?? f, this.meta.coneIds);
    if (res) this.legendRange = [res.lo, res.hi];
    const wf = this.shading === 'wireframe';
    const pos = new Float32Array(nF * 9), col = new Float32Array(nF * 9);
    for (let k = 0; k < nF; k++) {
      const a = F[3 * k], b = F[3 * k + 1], c = F[3 * k + 2];
      const pa = this.rot(V, a), pb = this.rot(V, b), pc = this.rot(V, c);
      pos.set(pa, 9 * k); pos.set(pb, 9 * k + 3); pos.set(pc, 9 * k + 6);
      const ux = pb[0] - pa[0], uy = pb[1] - pa[1], uz = pb[2] - pa[2], wx = pc[0] - pa[0], wy = pc[1] - pa[1], wz = pc[2] - pa[2];
      const nx = uy * wz - uz * wy, ny = uz * wx - ux * wz, nz = ux * wy - uy * wx;
      // normal into view space (rotation part of the view matrix)
      const vx = view[0] * nx + view[4] * ny + view[8] * nz, vy = view[1] * nx + view[5] * ny + view[9] * nz, vz = view[2] * nx + view[6] * ny + view[10] * nz;
      const si = Math.round(255 * shadeVal(vx, vy, vz));
      let r: number, g: number, bl: number;
      if (res) {
        // flow_curvature.py face(): Picnic colour x lighting (+ white specular highlight when Reflective)
        const cc = pickLut(res.faceT[k]);
        let lum = 0.62 + 0.38 * si / 255, hi = 0;
        if (this.shading === 'reflective') { const d = Math.max(0, Math.min(1, (si / 255 - 0.35) / 0.65)); lum = 0.28 + 0.72 * d; hi = 110 * Math.pow(d, 12); }
        r = Math.min(255, cc[0] * lum + hi) / 255; g = Math.min(255, cc[1] * lum + hi) / 255; bl = Math.min(255, cc[2] * lum + hi) / 255;
      } else {
        // Plain = the same lit faces compressed toward mid grey-blue (no specular range)
        const sj = this.shading === 'plain' ? Math.round(110 + 0.45 * (si - 110)) : si;
        r = BLUES_LUT[3 * sj]; g = BLUES_LUT[3 * sj + 1]; bl = BLUES_LUT[3 * sj + 2];
      }
      for (let q = 0; q < 3; q++) { col[9 * k + 3 * q] = r; col[9 * k + 3 * q + 1] = g; col[9 * k + 3 * q + 2] = bl; }
    }
    const g = this.mesh.geometry;
    g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    g.setAttribute('color', new THREE.BufferAttribute(col, 3));
    this.mesh.visible = !wf; this.wire.visible = !wf; this.wfLines.visible = wf;
    if (!wf) {
      this.wire.geometry.dispose();
      this.wire.geometry = new THREE.WireframeGeometry(g);
      (this.wire.material as THREE.LineBasicMaterial).color.set(this.dark ? 0xffffff : 0x000000);
    } else this.buildWireframe(f, res);
    this.renderer.render(this.scene, this.camera);
    this.placeCones(f);
    this.onLegend?.();
  }

  private rot(V: Float32Array, i: number): [number, number, number] {
    const x = V[3 * i], y = V[3 * i + 1], z = V[3 * i + 2], R = this.Rm;
    return [R[0][0] * x + R[0][1] * y + R[0][2] * z, R[1][0] * x + R[1][1] * y + R[1][2] * z, R[2][0] * x + R[2][1] * y + R[2][2] * z];
  }

  /** Wireframe: edges only, see-through, back edges dimmer; with an overlay each edge takes the mean of its faces. */
  private buildWireframe(f: FlowFrame, res: ReturnType<CurvatureOverlay['compute']>) {
    const V = f.vertices, view = this.camera.matrixWorldInverse.elements, M = this.meta!;
    const base = this.dark ? [110, 175, 225] : [8, 81, 156];
    const E = res ? this.curv.edgeColors(f, res)
      : Array.from({ length: M.edges.length / 2 }, (_, e) => ({ a: M.edges[2 * e], b: M.edges[2 * e + 1], rgb: base }));
    const zOf = (i: number) => { const p = this.rot(V, i); return view[2] * p[0] + view[6] * p[1] + view[10] * p[2]; };
    const zmx = Math.max(1e-9, this.radius);
    const pos = new Float32Array(E.length * 6), col = new Float32Array(E.length * 8);
    E.forEach((e, k) => {
      const zq = Math.max(0, Math.min(1, 0.5 + 0.5 * (zOf(e.a) + zOf(e.b)) / (2 * zmx))), al = 0.28 + 0.67 * zq;
      [e.a, e.b].forEach((v, j) => {
        const q = this.rot(V, v); pos[6 * k + 3 * j] = q[0]; pos[6 * k + 3 * j + 1] = q[1]; pos[6 * k + 3 * j + 2] = q[2];
        col[8 * k + 4 * j] = e.rgb[0] / 255; col[8 * k + 4 * j + 1] = e.rgb[1] / 255; col[8 * k + 4 * j + 2] = e.rgb[2] / 255; col[8 * k + 4 * j + 3] = al;
      });
    });
    const g = this.wfLines.geometry;
    g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    g.setAttribute('color', new THREE.BufferAttribute(col, 4));
  }

  /** Cone markers: red disc with a white ring and the cone number; cones on the far side are small, faded and unnumbered. */
  private placeCones(f: FlowFrame) {
    if (!this.meta) return;
    const w = this.canvas.clientWidth, h = this.canvas.clientHeight, v = new THREE.Vector3();
    const view = this.camera.matrixWorldInverse.elements;
    this.meta.coneIds.forEach((id, i) => {
      const [x, y, z] = this.rot(f.vertices, id);
      const depth = (view[2] * x + view[6] * y + view[10] * z) / this.radius; // positive toward the viewer
      v.set(x, y, z).project(this.camera);
      const el = this.coneEls[i]; if (!el) return;
      el.classList.toggle('back', depth <= -0.05);
      el.style.transform = `translate(${(v.x * 0.5 + 0.5) * w}px, ${(-v.y * 0.5 + 0.5) * h}px) translate(-50%,-50%)`;
    });
  }
}
