// Port of viewer/frontend/shading.py: matplotlib Blues ramp + view-space light.
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

/** Diverging blue - white - red, s in [-1, 1]. */
export function diverge(s: number, out: number[] = [0, 0, 0]): number[] {
  s = Math.max(-1, Math.min(1, s));
  if (s < 0) { const k = 1 + s; out[0] = 0.1 + 0.85 * k; out[1] = 0.25 + 0.7 * k; out[2] = 1; }
  else { out[0] = 1; out[1] = 0.95 - 0.8 * s; out[2] = 0.95 - 0.7 * s; }
  return out;
}
