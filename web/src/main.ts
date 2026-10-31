import './style.css';
import { View3D, type Overlay, type Shading } from './view3d';
import { CfDiscs } from './cfview';
import { PICNIC_CSS } from './shading';
import { SENSITIVITY, SENSITIVITY_DEFAULT } from './curvature';
import { Player, type PFrame } from './playback';
import type { FlowFrame, FlowMeta, FlowRequest, FromWorker, SeedSurface, ToWorker } from './solver';

const $ = <T extends HTMLElement>(id: string) => document.getElementById(id) as T;
const params = new URLSearchParams(location.search);
const fixturesMode = params.get('fixtures') === '1'; // testing only: stream golden fixtures instead of solving live
const fixturesBase = new URL('fixtures/', document.baseURI).href;

type Level = 'coarse' | 'medium' | 'fine' | 'extra';
/** A true frame as displayed: the solver data plus its vertex triples for the smooth-playback player. */
type Fr = FlowFrame & PFrame;

// ---------- theme ----------
let dark = true; // the Python page defaults to dark; the choice is remembered under the same key
try { const s = localStorage.getItem('sr-theme'); if (s) dark = s === 'dark'; } catch { /* ignore */ }

const view = new View3D($('cv') as HTMLCanvasElement, $('lbls'));
const cf = new CfDiscs([$('cf-n') as HTMLCanvasElement, $('cf-s') as HTMLCanvasElement]);
$('lg-bar').style.background = PICNIC_CSS;

function applyTheme() {
  document.documentElement.dataset.theme = dark ? 'dark' : 'light';
  $('theme').textContent = dark ? 'Light mode' : 'Dark mode';
  view.setTheme(dark); cf.setTheme(dark);
  try { localStorage.setItem('sr-theme', dark ? 'dark' : 'light'); } catch { /* ignore */ }
  redraw();
}
$('theme').onclick = () => { dark = !dark; applyTheme(); };

// ---------- radio lists ----------
type Opt = [value: string, text: string, disabledTip?: string];
function radios(host: HTMLElement, name: string, opts: Opt[], init: string, onChange: (v: string) => void) {
  host.textContent = '';
  for (const [v, text, off] of opts) {
    const l = document.createElement('label'), r = document.createElement('input'), s = document.createElement('span');
    r.type = 'radio'; r.name = name; r.value = v; r.checked = v === init; s.textContent = text;
    if (off) { r.disabled = true; l.title = off; l.style.cursor = 'not-allowed'; l.style.opacity = '.55'; }
    r.onchange = () => { if (r.checked) onChange(v); };
    l.append(r, s); host.appendChild(l);
  }
  return (v: string) => { host.querySelectorAll('input').forEach(i => { i.checked = i.value === v; }); };
}
let seed: SeedSurface = 'gluing', level: Level = 'coarse', smoothOn = true;
const setSeedUI = radios($('seed'), 'seed', [['gluing', 'Cone Gluing'], ['smart', 'Smart Seed'], ['polytope', 'Polytope']], seed, v => { seed = v as SeedSurface; configChanged(); });
radios($('realizer'), 'realizer', [['true_flow', 'True Ricci Flow'], ['fast', 'Fast Uniformizer', 'not yet ported']], 'true_flow', () => { /* only True Ricci Flow is available */ });
const setLevelUI = radios($('ressel'), 'ressel', [['coarse', 'Coarse'], ['medium', 'Medium'], ['fine', 'Fine'], ['extra', 'Extra Fine']], level, v => { level = v as Level; configChanged(); });
const setSmoothUI = radios($('fpsmooth'), 'fpsmooth', [['off', 'Off'], ['on', 'On']], 'on', v => { smoothOn = v === 'on'; });
radios($('shade'), 'shade', [['reflective', 'Reflective'], ['plain', 'Plain'], ['wireframe', 'Wireframe']], 'reflective', v => { view.shading = v as Shading; view.render(); });
radios($('overlay'), 'overlay', [['none', 'None'], ['length', 'Length Fit'], ['gauss', 'Gaussian Curvature']], 'none', v => { $('sens').hidden = v === 'none'; view.setOverlay(v as Overlay); updateLegend(); });
radios($('curvscale'), 'curvscale', SENSITIVITY.map(([v, n]) => [v, n] as Opt), SENSITIVITY_DEFAULT, v => view.setSensitivity(parseFloat(v)));

function updateLegend() {
  const lg = $('legend'), on = view.overlay !== 'none';
  lg.hidden = !on;
  if (!on) return;
  const [mlo, mhi] = view.legendRange;
  $('lg-title').textContent = view.overlay === 'length' ? 'Length Fit' : 'Gaussian Curvature';
  $('lg-mid').textContent = 'Target';
  $('lg-lo').textContent = '−' + mlo.toPrecision(2); $('lg-hi').textContent = '+' + mhi.toPrecision(2);
}
view.onLegend = updateLegend;

// ---------- cone controls ----------
let betas: number[] = [], positions: [number, number][] = [];
const defaultBetas = [0.6, 0.7, 0.5, 0.8, 0.7];
function ngon(n: number): [number, number][] { return Array.from({ length: n }, (_, i) => [0.7 * Math.cos(2 * Math.PI * i / n), 0.7 * Math.sin(2 * Math.PI * i / n)] as [number, number]); }
function rhoOf() { return (betas.reduce((a, b) => a + b, 0) - (betas.length - 2)) / 2; }
function updateRho() {
  const r = rhoOf();
  $('rhoval').textContent = (r >= 0 ? '+' : '−') + Math.abs(r).toFixed(3);
  $('rhonote').textContent = r > 0 ? '(sphere-like)' : r < 0 ? '(hyperbolic-like)' : '(flat)';
}
/** Any edit of the cones cancels a running flow and returns to the live preview (as on the Python page). */
function configChanged() {
  if (running || frames.length || meta) toPreview();
  updateRho(); cf.setPreview(betas, positions); cf.setMarkers(true);
  requestPreview(false);
}
let sliderEls: { sl: HTMLInputElement; lab: HTMLLabelElement }[] = [];
function buildSliders(n: number) {
  const host = $('sliders'); host.textContent = ''; sliderEls = [];
  betas = Array.from({ length: n }, (_, i) => betas[i] ?? defaultBetas[i % 5]);
  positions = ngon(n);
  for (let i = 0; i < n; i++) {
    const row = document.createElement('div'); row.className = 'row';
    const lab = document.createElement('label'); lab.textContent = `β = ${betas[i].toFixed(3)}`;
    const sl = document.createElement('input'); sl.type = 'range'; sl.min = '0.05'; sl.max = '0.95'; sl.step = '0.01'; sl.value = String(betas[i]);
    sl.oninput = () => { betas[i] = +sl.value; lab.textContent = `β = ${betas[i].toFixed(3)}`; configChanged(); };
    row.append(lab, sl); host.appendChild(row); sliderEls.push({ sl, lab });
  }
  configChanged();
}
$('nsel').onchange = () => { betas = []; buildSliders(+($('nsel') as HTMLSelectElement).value); };
cf.onDragCone = (i, pos, fin) => {
  positions[i] = [pos[0], pos[1]];
  if (running || frames.length) toPreview();
  cf.setPreview(betas, positions); cf.setMarkers(true);
  requestPreview(!fin);   // live: throttled while dragging, latest wins
};

// ---------- solver worker ----------
const spawn = () => new Worker(new URL('./worker.ts', import.meta.url), { type: 'module' });
let worker = spawn();
let runId = 0, running = false, meta: FlowMeta | null = null, t0 = 0, wallMs = 0, result: { rho: number; label: string; detail: string; ok: boolean } | null = null;
let frames: Fr[] = [];
let cur = 0;
let cfShown: Fr | null = null;
worker.onmessage = onWorker;

// ---------- live seed preview (own worker, so a running flow is never disturbed) ----------
const pworker = spawn();
let pId = 0, pBusy = false, pDirty = false, pTimer: number | null = null, pLast = 0, pReq: FlowRequest | null = null;
let previewMs = NaN;
function currentReq(): FlowRequest { return { n: betas.length, betas: betas.slice(), positions: positions.map(p => [...p] as [number, number]), seed, level }; }
/** Debounced (150 ms) for sliders/radios, throttled (~60 ms) during cone drags; one request in flight, latest wins. */
function requestPreview(drag: boolean) {
  if (running) return;
  if (pTimer != null) clearTimeout(pTimer);
  const wait = drag ? Math.max(0, 60 - (performance.now() - pLast)) : 150;
  pTimer = window.setTimeout(() => { pTimer = null; sendPreview(); }, wait);
}
function sendPreview() {
  if (running) return;
  if (pBusy) { pDirty = true; return; }
  pBusy = true; pDirty = false; pReq = currentReq(); pLast = performance.now();
  pworker.postMessage({ cmd: 'preview', id: ++pId, req: pReq } satisfies ToWorker);
}
pworker.onmessage = (ev: MessageEvent<FromWorker>) => {
  const m = ev.data;
  if (m.type !== 'preview' && m.type !== 'preview-fail') return;
  pBusy = false;
  if (pDirty) { sendPreview(); }
  if (running || frames.length || meta) return;
  if (m.type === 'preview') {
    previewMs = m.data.ms; view.setPreview(m.data, pReq!);
    $('st-v').textContent = String(m.data.vertices.length / 3); $('st-e').textContent = String(m.data.edges.length / 2); $('st-f').textContent = String(m.data.faces.length / 3);
    $('warn').textContent = m.data.clamped ? 'Overlapping cone caps were shrunk' : '';
  } else { view.clear(); $('warn').textContent = ''; }
};

/** A chunk of the solver cannot be interrupted, so Stop terminates the worker and starts a fresh one. */
function killWorker() { worker.terminate(); worker = spawn(); worker.onmessage = onWorker; }

function triples(v: Float32Array): number[][] { const o: number[][] = new Array(v.length / 3); for (let i = 0; i < o.length; i++) o[i] = [v[3 * i], v[3 * i + 1], v[3 * i + 2]]; return o; }
function flat(V: number[][]): Float32Array { const o = new Float32Array(V.length * 3); V.forEach((p, i) => { o[3 * i] = p[0]; o[3 * i + 1] = p[1]; o[3 * i + 2] = p[2]; }); return o; }

function onWorker(ev: MessageEvent<FromWorker>) {
  const m = ev.data;
  if (m.id !== runId) return;
  if (m.type === 'meta') {
    meta = m.meta; frames = []; cfShown = null;
    view.reset(m.meta); cf.reset(m.meta);
    $('st-v').textContent = String(m.meta.sphereVertices.length / 3);
    $('st-e').textContent = String(m.meta.edges.length / 2);
    $('st-f').textContent = String(m.meta.sphereFaces.length / 3);
    $('notice').hidden = true;
    if (m.meta.notice) $('warn').textContent = m.meta.notice;
    setStatus(`Flowing${m.meta.request.fixture ? ` (${m.meta.request.fixture})` : ''}: rho = ${m.meta.rho.toPrecision(4)}`);
  } else if (m.type === 'progress') {
    setStatus(m.phase); $('loadpct').textContent = `${Math.round(100 * m.progress)}%`;
  } else if (m.type === 'frame') {
    const f = m.frame as Fr; f.V = triples(f.vertices); frames.push(f);
    const sl = $('play') as HTMLInputElement; sl.max = String(frames.length - 1);
    $('loadpct').textContent = `${Math.round(100 * f.progress)}%`;
    ($('playbtn') as HTMLButtonElement).disabled = false; sl.disabled = false;
    player.show(f);
    if (f.quality) {
      setQuality(f.quality.label, f.quality.ok, f.quality.detail);
      result = { rho: meta?.rho ?? NaN, label: f.quality.label, detail: f.quality.detail ?? '', ok: f.quality.ok };
    }
  } else if (m.type === 'done') {
    finish(m.stopped);
  } else if (m.type === 'error') {
    running = false; wallMs = performance.now() - t0; lock(frames.length > 0);
    setStatus('FAILED - ' + m.message);
  }
}

function setStatus(s: string) { $('status').textContent = s; }
function setQuality(label: string, ok: boolean, detail?: string) {
  const q = $('qnote'); q.textContent = label; q.className = ok ? 'good' : 'bad'; q.title = detail ?? '';
}
function fmt(x: number | undefined | null) { return x == null || !isFinite(x) ? '' : x === 0 ? '0' : Math.abs(x) < 1e-3 || Math.abs(x) >= 1e4 ? x.toExponential(1) : x.toPrecision(3); }

// ---------- display + smooth playback ----------
/** Port of frameData() in embed_flow.py: a smooth sub-frame blends rho, t, Kdef and Avert with the same fraction as its geometry. Cached per sub-frame. */
type Sub = Fr & { _sh?: FlowFrame };
function subFrame(fr: Fr, to: Fr): FlowFrame {
  const sf = fr as Sub; if (sf._sh) return sf._sh;
  const A = fr.fp!.from as Fr, s = fr.fp!.i / fr.fp!.n;
  const mix = (a: number | undefined, b: number | undefined) => (a == null || b == null ? b : (1 - s) * a + s * b);
  const mix1 = (a: Float32Array | undefined, b: Float32Array | undefined) => {
    if (!a || !b || a.length !== b.length) return b; const o = new Float32Array(b.length); for (let i = 0; i < o.length; i++) o[i] = (1 - s) * a[i] + s * b[i]; return o;
  };
  return sf._sh = { ...to, vertices: flat(fr.V), rho: mix(A.rho, to.rho), t: mix(A.t, to.t)!, Kdef: mix1(A.Kdef, to.Kdef), Avert: mix1(A.Avert, to.Avert) };
}
function display(p: PFrame) {
  const fr = p as Fr, to = (fr.fp ? fr.fp.to : fr) as Fr;
  if (!meta) return;
  const shown: FlowFrame = fr.fp ? subFrame(fr, to) : to;
  view.setFrame(shown, frames[0]);
  if (cfShown !== to) { cf.setFrame(to); cfShown = to; }
  $('st-l').textContent = to.lenRelRms == null ? '' : (100 * to.lenRelRms).toFixed(2) + ' %';
  $('st-c').textContent = fmt(to.coneDefectErr);
  if (!fr.fp) {
    cur = Math.max(0, frames.indexOf(to));
    const sl = $('play') as HTMLInputElement; sl.value = String(cur);
    $('tval').textContent = `t = ${to.t.toPrecision(3)}`;
    $('playlbl').textContent = `frame ${cur + 1} / ${meta.nFrames ?? frames.length}`;
  }
  updateLegend();
}
const player = new Player({
  frames: () => frames, set: display,
  setSlider(k) { if (k < 0) return; cur = k; ($('play') as HTMLInputElement).value = String(k); const f = frames[k]; if (f) { $('tval').textContent = `t = ${f.t.toPrecision(3)}`; $('playlbl').textContent = `frame ${k + 1} / ${meta?.nFrames ?? frames.length}`; } },
  follow: () => true, smooth: () => smoothOn,
  setPlayLabel: t => { $('playbtn').textContent = t; },
  sliderValue: () => parseInt(($('play') as HTMLInputElement).value, 10),
});
function showFrame() { const f = frames[cur]; if (f) display(f); }
function redraw() { if (frames.length) showFrame(); }

function lock(haveFrames: boolean) {
  ($('play') as HTMLInputElement).disabled = !haveFrames;
  ($('playbtn') as HTMLButtonElement).disabled = !haveFrames;
  ($('stop') as HTMLButtonElement).disabled = !running;
  ($('gen') as HTMLButtonElement).disabled = running;
}
function finish(stopped: boolean) {
  running = false; wallMs = performance.now() - t0; lock(frames.length > 0);
  if (!frames.length) { setStatus(stopped ? 'Stopped.' : 'Done: no frames.'); return; }
  setStatus(stopped ? `Stopped at frame ${cur + 1}.` : `Done: ${frames.length} frames, t = 0 .. ${frames[frames.length - 1].t.toPrecision(3)}.`);
  if (!stopped) $('loadpct').textContent = '100%';
  if (stopped) player.stop(true);
  else if (player.busy()) { // let the smooth tail finish, but never rest mid-run
    const last = frames.length; setTimeout(() => { if (frames.length === last && player.busy()) player.stop(true); showEnd(); }, 2500);
  } else showEnd();
}
/** The Python page ends a run on the final frame. */
function showEnd() {
  const k = frames.length - 1; if (k < 0 || running || cur === k) return;
  player.seek(frames[k]); cur = k; ($('play') as HTMLInputElement).value = String(k); showFrame();
}

($('play') as HTMLInputElement).oninput = e => { const k = +(e.target as HTMLInputElement).value; player.seek(frames[k]); cur = k; showFrame(); };
$('playbtn').onclick = () => player.play();

function clearStats() { for (const id of ['st-v', 'st-e', 'st-f', 'st-l', 'st-c']) $(id).textContent = ''; }
/** Drop the shown flow and return to the cone preview. */
function toPreview(keepView = false) {
  player.stop(false);
  if (running) { runId++; killWorker(); running = false; }
  frames = []; meta = null; cur = 0; cfShown = null; result = null;
  if (!keepView) { view.clear(); clearStats(); requestPreview(false); }
  $('loadpct').textContent = ''; $('qnote').textContent = ''; $('qnote').title = ''; $('tval').textContent = ''; $('playlbl').textContent = ''; $('notice').hidden = true; $('warn').textContent = '';
  ($('play') as HTMLInputElement).max = '0'; ($('play') as HTMLInputElement).value = '0';
  setStatus(''); lock(false);
}

function start() {
  toPreview(true);
  const req: FlowRequest = { n: betas.length, betas: betas.slice(), positions: positions.map(p => [...p] as [number, number]), seed, level };
  cf.setMarkers(false);
  runId++; running = true; t0 = performance.now(); wallMs = 0; lock(false);
  ($('stop') as HTMLButtonElement).disabled = false;
  $('loadpct').textContent = '0%';
  setStatus('Starting…');
  worker.postMessage({ cmd: 'start', id: runId, req, solver: fixturesMode ? 'fixture' : 'wasm', baseUrl: fixturesBase } satisfies ToWorker);
}
$('gen').onclick = start;
$('stop').onclick = () => {
  if (!running) return;
  runId++; killWorker(); // a chunk cannot be interrupted: terminate and respawn
  finish(true);
};
$('reset').onclick = () => {
  toPreview();
  seed = 'gluing'; setSeedUI(seed); level = 'coarse'; setLevelUI(level); smoothOn = true; setSmoothUI('on');
  ($('nsel') as HTMLSelectElement).value = '4'; betas = []; buildSliders(4);
};

// test hook for screenshots / automation
function setConfig(c: { positions?: number[][]; betas?: number[]; seed?: SeedSurface; level?: Level; smooth?: boolean }) {
  if (c.betas) {
    ($('nsel') as HTMLSelectElement).value = String(c.betas.length); betas = c.betas.slice(); buildSliders(c.betas.length);
    betas.forEach((b, i) => { sliderEls[i].sl.value = String(b); sliderEls[i].lab.textContent = `β = ${b.toFixed(3)}`; });
  }
  if (c.positions) positions = c.positions.map(p => [p[0], p[1]] as [number, number]);
  if (c.seed) { seed = c.seed; setSeedUI(seed); }
  if (c.level) { level = c.level; setLevelUI(level); }
  if (c.smooth != null) { smoothOn = c.smooth; setSmoothUI(smoothOn ? 'on' : 'off'); }
  configChanged();
}
(window as unknown as Record<string, unknown>).__flow = {
  frames: () => frames.length, running: () => running, setConfig, generate: start,
  seek: (k: number) => { player.seek(frames[Math.min(k, frames.length - 1)]); cur = Math.min(k, frames.length - 1); showFrame(); },
  result: () => result, previewMs: () => previewMs, previewReady: () => !pBusy && pTimer == null && !pDirty, dragCone: (i: number, p: [number, number], fin: boolean) => cf.onDragCone?.(i, p, fin),
  /** test: max |delta faceT| between consecutive linear sub-frames of frames k..k+1, old (B's Kdef) vs blended */
  faceTProbe: (k: number, n: number) => {
    const A = frames[k], B = frames[k + 1]; view.setOverlay('length');
    const run = (blend: boolean) => {
      let prev: Float32Array | null = null, mx = 0;
      for (let i = 1; i <= n; i++) {
        const V = A.V.map((a, j) => a.map((x, c) => (1 - i / n) * x + (i / n) * B.V[j][c]));
        const fr = { ...B, V, fp: { from: A, to: B, i, n } } as Fr;
        const sh = blend ? subFrame(fr, B) : { ...B, vertices: flat(V) };
        const r = view.curv.compute(sh, frames[0], meta!.coneIds); if (!r) continue;
        if (prev) for (let q = 0; q < prev.length; q++) { const d = Math.abs(r.faceT[q] - prev[q]); if (d === d) mx = Math.max(mx, d); }
        prev = Float32Array.from(r.faceT);
      }
      return mx;
    };
    return { before: run(false), after: run(true) };
  }, wallMs: () => wallMs, status: () => $('status').textContent, positions: () => positions, betas: () => betas,
};

buildSliders(4); lock(false);
applyTheme();
