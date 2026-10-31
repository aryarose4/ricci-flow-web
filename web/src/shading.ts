// Port of viewer/frontend/shading.py: matplotlib Blues ramp + view-space light.
import * as THREE from 'three';
// Colours below are the page's final sRGB values; keep three.js from re-encoding them.
THREE.ColorManagement.enabled = false;
export const RAW_COLOR_SPACE = THREE.LinearSRGBColorSpace;
const BLUES: number[][] = [
  [247,251,255],[241,247,253],[234,243,251],[228,239,249],[222,235,247],[216,231,245],[210,227,243],[204,223,241],
  [198,219,239],[188,215,235],[178,210,232],[168,206,228],[157,202,225],[145,195,222],[132,188,219],[119,181,217],
  [106,174,214],[96,167,210],[86,160,206],[75,152,202],[65,145,198],[57,137,193],[49,129,189],[41,121,185],
  [32,112,180],[26,104,174],[20,96,168],[14,88,162],[8,80,155],[8,72,142],[8,64,130],[8,55,118],[8,48,107]];
export const LIGHT_VIEW = [0.4, -0.5, 0.77];

/** LUT of 256 rgb triples in 0..1. */
export const BLUES_LUT: Float32Array = (() => {
  const out = new Float32Array(256 * 3), n = BLUES.length - 1;
  for (let i = 0; i < 256; i++) {
    const x = (i / 255) * n, j = Math.min(n - 1, Math.floor(x)), t = x - j, a = BLUES[j], b = BLUES[j + 1];
    for (let c = 0; c < 3; c++) out[i * 3 + c] = (a[c] + t * (b[c] - a[c])) / 255;
  }
  return out;
})();

export function shadeVal(nx: number, ny: number, nz: number): number {
  const nl = Math.hypot(nx, ny, nz) || 1;
  const d = (nx * LIGHT_VIEW[0] + ny * LIGHT_VIEW[1] + nz * LIGHT_VIEW[2]) / nl;
  return 0.35 + 0.65 * Math.min(1, Math.max(0, d));
}

// matplotlib magma (7 stops, as in conformal_view.py)
const MAGMA = [[0,0,4],[40,17,89],[110,31,129],[182,54,121],[241,96,93],[254,159,109],[252,253,191]];
export function magma(t: number, out: number[] = [0, 0, 0]): number[] {
  t = Math.max(0, Math.min(1, t)); t = 0.03 + 0.97 * t; t *= MAGMA.length - 1;
  const j = Math.min(MAGMA.length - 2, Math.floor(t)), f = t - j, a = MAGMA[j], b = MAGMA[j + 1];
  for (let c = 0; c < 3; c++) out[c] = (a[c] + f * (b[c] - a[c])) / 255;
  return out;
}

// flow_curvature.py: overlay colormap (DIVERGING_CMAP = "Picnic"), stops listed LOW -> HIGH, neutral white at the centre.
export const PICNIC: number[][] = [[0,0,255],[51,153,255],[102,204,255],[153,204,255],[204,204,255],[255,255,255],[255,204,255],[255,153,255],[255,102,204],[255,51,153],[255,0,0]];
export const PICNIC_CSS = 'linear-gradient(90deg,' + PICNIC.map(c => `rgb(${c.join(',')})`).join(',') + ')';
export const OVERLAY_GREY = [120, 120, 120];
/** 257-entry LUT (0..255 rgb) indexed by round(128 + 128 t), t in [-1,1]. */
export const PICNIC_LUT: number[][] = (() => {
  const out: number[][] = [], n = PICNIC.length - 1;
  for (let i = 0; i <= 256; i++) {
    const x = (i / 256) * n, j = Math.min(n - 1, Math.floor(x)), f = x - j, a = PICNIC[j], b = PICNIC[j + 1];
    out.push([a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f]);
  }
  return out;
})();
export function pickLut(t: number): number[] { return Number.isNaN(t) ? OVERLAY_GREY : PICNIC_LUT[Math.round(128 + 128 * Math.max(-1, Math.min(1, t)))]; }
