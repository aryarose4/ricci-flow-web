import './style.css';
import { View3D, type Overlay, type Shading } from './view3d';
import { CfView } from './cfview';
import type { FlowFrame, FlowMeta, FlowRequest, FromWorker, SeedSurface, ToWorker } from './solver';

const $ = <T extends HTMLElement>(id: string) => document.getElementById(id) as T;
const fixturesBase = new URL('fixtures/', document.baseURI).href;

// ---------- theme ----------
let dark = matchMedia('(prefers-color-scheme: dark)').matches;
try { const s = localStorage.getItem('theme'); if (s) dark = s === 'dark'; } catch { /* ignore */ }

const view = new View3D($('cv') as HTMLCanvasElement, $('lbls'));
const cfN = new CfView($('cf-n') as HTMLCanvasElement, true);
const cfS = new CfView($('cf-s') as HTMLCanvasElement, false);

function applyTheme() {
  document.documentElement.dataset.theme = dark ? 'dark' : 'light';
  $('theme').textContent = dark ? 'Light mode' : 'Dark mode';
  view.setTheme(dark); cfN.setTheme(dark); cfS.setTheme(dark);
  try { localStorage.setItem('theme', dark ? 'dark' : 'light'); } catch { /* ignore */ }
  redraw();
}
$('theme').onclick = () => { dark = !dark; applyTheme(); };

// ---------- radio lists ----------
function radios(host: HTMLElement, name: string, opts: [string, string][], init: string, onChange: (v: string) => void) {
  host.textContent = '';
  for (const [v, text] of opts) {
    const l = document.createElement('label'), r = document.createElement('input'), s = document.createElement('span');
    r.type = 'radio'; r.name = name; r.value = v; r.checked = v === init; s.textContent = text;
    r.onchange = () => { if (r.checked) onChange(v); };
    l.append(r, s); host.appendChild(l);
  }
  return (v: string) => { host.querySelectorAll('input').forEach(i => { i.checked = i.value === v; }); };
}
let seed: SeedSurface = 'gluing';
const setSeedUI = radios($('seed'), 'seed', [['gluing', 'Cone Gluing'], ['smart', 'Smart Seed'], ['polytope', 'Polytope']], seed, v => { seed = v as SeedSurface; });
radios($('shade'), 'shade', [['reflective', 'Reflective'], ['plain', 'Plain'], ['wireframe', 'Wireframe']], 'reflective', v => { view.shading = v as Shading; view.render(); });
radios($('overlay'), 'overlay', [['none', 'None'], ['length', 'Length Fit'], ['gauss', 'Gaussian Curvature']], 'none', v => { view.overlay = v as Overlay; view.render(); updateLegend(); });

function updateLegend() {
  const lg = $('legend'), on = view.overlay !== 'none';
  lg.hidden = !on;
  if (!on) return;
  const m = view.legendRange, f = (x: number) => (x >= 0 ? '+' : '−') + Math.abs(x).toPrecision(2);
  $('lg-title').textContent = view.overlay === 'length' ? 'Length Fit' : 'Gaussian Curvature';
  $('lg-mid').textContent = view.overlay === 'length' ? 'Target' : 'ρ';
  $('lg-lo').textContent = f(-m); $('lg-hi').textContent = f(m);
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
function buildSliders(n: number) {
  const host = $('sliders'); host.textContent = '';
  betas = Array.from({ length: n }, (_, i) => betas[i] ?? defaultBetas[i % 5]);
  positions = ngon(n);
  for (let i = 0; i < n; i++) {
    const row = document.createElement('div'); row.className = 'row';
    const lab = document.createElement('label'); lab.textContent = `β${i + 1} = ${betas[i].toFixed(2)}`;
    const sl = document.createElement('input'); sl.type = 'range'; sl.min = '0.05'; sl.max = '0.95'; sl.step = '0.01'; sl.value = String(betas[i]);
    sl.oninput = () => { betas[i] = +sl.value; lab.textContent = `β${i + 1} = ${betas[i].toFixed(2)}`; updateRho(); };
    const px = document.createElement('input'), py = document.createElement('input');
    for (const [el, k] of [[px, 0], [py, 1]] as const) {
      el.type = 'number'; el.step = '0.05'; el.value = positions[i][k].toFixed(2); el.title = k ? 'y position (stereographic plane)' : 'x position (stereographic plane)';
      el.onchange = () => { positions[i][k] = +el.value; };
    }
    row.append(lab, sl, px, py); host.appendChild(row);
  }
  updateRho();
}
$('nsel').onchange = () => { betas = []; buildSliders(+($('nsel') as HTMLSelectElement).value); };

// ---------- fixture list ----------
interface Idx { name: string; n: number; betas: number[]; seed: SeedSurface; rho: number; quality: string }
let fixtures: Idx[] = [];
const fixSel = $('fixture') as HTMLSelectElement;
async function loadIndex() {
  fixtures = await (await fetch(fixturesBase + 'index.json')).json();
  fixSel.textContent = '';
  const auto = new Option('Nearest match', ''); fixSel.add(auto);
  for (const f of fixtures) fixSel.add(new Option(`${f.name}  (${f.quality.split(' (')[0]})`, f.name));
}
fixSel.onchange = () => {
  const f = fixtures.find(x => x.name === fixSel.value); if (!f) return;
  ($('nsel') as HTMLSelectElement).value = String(f.n);
  betas = f.betas.slice(); seed = f.seed; setSeedUI(seed);
  buildSliders(f.n);
};

// ---------- solver worker ----------
const worker = new Worker(new URL('./worker.ts', import.meta.url), { type: 'module' });
let runId = 0, running = false, meta: FlowMeta | null = null;
let frames: FlowFrame[] = [];
let lo = Infinity, hi = -Infinity, cur = 0;
let timer: number | null = null;

worker.onmessage = (ev: MessageEvent<FromWorker>) => {
  const m = ev.data;
  if (m.id !== runId) return;
  if (m.type === 'meta') {
    meta = m.meta; frames = []; lo = Infinity; hi = -Infinity;
    view.reset(m.meta); cfN.reset(m.meta); cfS.reset(m.meta);
    $('st-v').textContent = String(m.meta.sphereVertices.length / 3);
    $('st-e').textContent = String(m.meta.edges.length / 2);
    $('st-f').textContent = String(m.meta.sphereFaces.length / 3);
    const n = $('notice'); n.hidden = true;
    setStatus(`Flowing (${m.meta.request.fixture ?? 'live'}): ${m.meta.label ?? ''}`);
  } else if (m.type === 'frame') {
    const f = m.frame; frames.push(f);
    for (const x of f.u) { if (x < lo) lo = x; if (x > hi) hi = x; }
    cur = frames.length - 1; showFrame();
    $('loadpct').textContent = `${Math.round(100 * f.progress)}%`;
    if (f.quality) setQuality(f.quality.label, f.quality.ok, f.quality.detail);
  } else if (m.type === 'done') {
    finish(m.stopped);
  } else if (m.type === 'error') {
    $('warn').textContent = m.message; finish(true);
  }
};

function setStatus(s: string) { $('status').textContent = s; }
function setQuality(label: string, ok: boolean, detail?: string) {
  const q = $('qnote'); q.textContent = label; q.className = ok ? 'good' : 'bad'; q.title = detail ?? '';
}
function fmt(x: number | undefined | null) { return x == null || !isFinite(x) ? '-' : x === 0 ? '0' : Math.abs(x) < 1e-3 || Math.abs(x) >= 1e4 ? x.toExponential(1) : x.toPrecision(3); }

function showFrame() {
  const f = frames[cur]; if (!f || !meta) return;
  view.setFrame(f); cfN.setFrame(f, lo, hi); cfS.setFrame(f, lo, hi);
  const sl = $('play') as HTMLInputElement; sl.max = String(Math.max(0, frames.length - 1)); sl.value = String(cur);
  $('tval').textContent = `t = ${f.t.toPrecision(3)}`;
  $('playlbl').textContent = `frame ${cur + 1} / ${meta.nFrames ?? frames.length}`;
  $('st-l').textContent = f.lenRelRms == null ? '-' : (100 * f.lenRelRms).toFixed(2) + ' %';
  $('st-c').textContent = fmt(f.coneDefectErr);
  $('st-g').textContent = fmt(f.gbError);
  const fl = $('notice');
  if (meta.notice && cur === frames.length - 1 && !running) { fl.textContent = meta.notice; fl.hidden = false; } else fl.hidden = true;
  updateLegend();
}
function redraw() { if (frames.length) showFrame(); }

function lock(locked: boolean) {
  ($('play') as HTMLInputElement).disabled = locked;
  ($('playbtn') as HTMLButtonElement).disabled = locked;
  ($('stop') as HTMLButtonElement).disabled = !running;
  ($('gen') as HTMLButtonElement).disabled = running;
}
function finish(stopped: boolean) {
  running = false; lock(false);
  if (!frames.length) { lock(true); return; }
  setStatus(stopped ? `Stopped at frame ${cur + 1}.` : `Done: ${frames.length} frames, t = 0 .. ${frames[frames.length - 1].t.toPrecision(3)}.`);
  if (!stopped) $('loadpct').textContent = '100%';
  showFrame();
}

function stopPlay() { if (timer != null) { clearInterval(timer); timer = null; } $('playbtn').textContent = 'Play'; }
($('play') as HTMLInputElement).oninput = e => { stopPlay(); cur = +(e.target as HTMLInputElement).value; showFrame(); };
$('playbtn').onclick = () => {
  if (timer != null) { stopPlay(); return; }
  if (cur >= frames.length - 1) cur = 0;
  $('playbtn').textContent = 'Pause';
  timer = window.setInterval(() => {
    if (cur >= frames.length - 1) { stopPlay(); return; }
    cur++; showFrame();
  }, 220);
};

function start() {
  stopPlay();
  if (running) worker.postMessage({ cmd: 'stop', id: runId } satisfies ToWorker);
  $('warn').textContent = ''; $('qnote').textContent = ''; $('loadpct').textContent = '0%';
  const req: FlowRequest = { n: betas.length, betas: betas.slice(), positions: positions.map(p => [...p] as [number, number]), seed, fixture: fixSel.value || undefined };
  runId++; running = true; frames = []; cur = 0; lock(true);
  ($('stop') as HTMLButtonElement).disabled = false;
  setStatus('Starting…');
  worker.postMessage({ cmd: 'start', id: runId, req, solver: 'fixture', baseUrl: fixturesBase } satisfies ToWorker);
}
$('gen').onclick = start;
$('stop').onclick = () => { worker.postMessage({ cmd: 'stop', id: runId } satisfies ToWorker); };
$('reset').onclick = () => {
  stopPlay(); if (running) worker.postMessage({ cmd: 'stop', id: runId } satisfies ToWorker);
  runId++; running = false; frames = []; meta = null; cur = 0;
  fixSel.value = ''; seed = 'gluing'; setSeedUI(seed); ($('nsel') as HTMLSelectElement).value = '4'; betas = []; buildSliders(4);
  for (const id of ['st-v', 'st-e', 'st-f', 'st-l', 'st-c', 'st-g']) $(id).textContent = '-';
  $('loadpct').textContent = ''; $('qnote').textContent = ''; $('tval').textContent = ''; $('playlbl').textContent = ''; $('notice').hidden = true;
  setStatus('Ready.'); lock(true); ($('gen') as HTMLButtonElement).disabled = false;
};

// test hook for screenshots / automation
(window as unknown as Record<string, unknown>).__flow = {
  frames: () => frames.length, running: () => running, seek: (k: number) => { stopPlay(); cur = Math.min(k, frames.length - 1); showFrame(); },
  setFixture: (name: string) => { fixSel.value = name; fixSel.onchange?.(new Event('change')); },
};

buildSliders(4); lock(true); ($('gen') as HTMLButtonElement).disabled = false;
applyTheme();
loadIndex().catch(e => { $('warn').textContent = 'Could not load fixtures: ' + e; });
