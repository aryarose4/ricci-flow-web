import * as THREE from 'three';
import { magma } from './shading';
import type { FlowFrame, FlowMeta } from './solver';

/** Conformal factor u on the north / south hemisphere, drawn as a stereographic disc (vertex-coloured, magma). */
export class CfView {
  private renderer: THREE.WebGLRenderer;
  private scene = new THREE.Scene();
  private cam = new THREE.OrthographicCamera(-1.05, 1.05, 1.05, -1.05, -1, 1);
  private mesh: THREE.Mesh | null = null;
  private nv = 0;

  constructor(private canvas: HTMLCanvasElement, private north: boolean) {
    this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, preserveDrawingBuffer: true });
    new ResizeObserver(() => this.resize()).observe(canvas.parentElement!);
    this.resize();
  }
  resize() {
    const p = this.canvas.parentElement!, s = Math.max(120, p.clientWidth);
    this.renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
    this.renderer.setSize(s, s, false);
    this.canvas.style.width = s + 'px'; this.canvas.style.height = s + 'px';
    this.renderer.render(this.scene, this.cam);
  }
  setTheme(dark: boolean) { this.renderer.setClearColor(dark ? 0x121417 : 0xf4f4f2, 1); this.renderer.render(this.scene, this.cam); }

  reset(meta: FlowMeta) {
    if (this.mesh) { this.scene.remove(this.mesh); this.mesh.geometry.dispose(); this.mesh = null; }
    const S = meta.sphereVertices, F = meta.sphereFaces, nv = S.length / 3, sgn = this.north ? 1 : -1;
    this.nv = nv;
    const pos = new Float32Array(nv * 3);
    for (let i = 0; i < nv; i++) {
      const z = S[3 * i + 2] * sgn, d = Math.max(1 + z, 1e-6); // project from the opposite pole; equator -> r = 1
      pos[3 * i] = S[3 * i] / d; pos[3 * i + 1] = S[3 * i + 1] * sgn / d; pos[3 * i + 2] = 0;
    }
    const keep: number[] = [];
    for (let k = 0; k < F.length / 3; k++) {
      const a = F[3 * k], b = F[3 * k + 1], c = F[3 * k + 2];
      if (((S[3 * a + 2] + S[3 * b + 2] + S[3 * c + 2]) / 3) * sgn >= -0.02) keep.push(a, b, c);
    }
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    g.setAttribute('color', new THREE.BufferAttribute(new Float32Array(nv * 3), 3));
    g.setIndex(new THREE.BufferAttribute(new Uint32Array(keep), 1));
    this.mesh = new THREE.Mesh(g, new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide }));
    this.scene.add(this.mesh);
    this.renderer.render(this.scene, this.cam);
  }

  /** [lo,hi]: shared range over all frames seen so far, so colours are comparable during playback. */
  setFrame(f: FlowFrame, lo: number, hi: number) {
    if (!this.mesh) return;
    const col = this.mesh.geometry.getAttribute('color') as THREE.BufferAttribute, tmp = [0, 0, 0], span = hi - lo || 1;
    for (let i = 0; i < this.nv; i++) { magma((f.u[i] - lo) / span, tmp); col.setXYZ(i, tmp[0], tmp[1], tmp[2]); }
    col.needsUpdate = true;
    this.renderer.render(this.scene, this.cam);
  }
}
