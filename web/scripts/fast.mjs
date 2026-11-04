// Fast Uniformizer + 2..8 cones browser check: node scripts/fast.mjs [port]
import { chromium } from 'playwright';
const port = process.argv[2] ?? '4393', out = new URL('../../artifacts/fast-uniformizer/', import.meta.url);
const b = await chromium.launch(); const p = await b.newPage({ viewport: { width: 1400, height: 900 } });
const errs = []; p.on('pageerror', e => errs.push(String(e))); p.on('console', m => { if (m.type() === 'error') errs.push(m.text()); });
await p.goto(`http://127.0.0.1:${port}/`); await p.waitForFunction(() => window.__flow);
const cfgs = { 2: [0.6, 0.6], 3: [0.6, 0.7, 0.8], 6: [0.9, 0.85, 0.9, 0.8, 0.9, 0.85], 8: [0.95, 0.9, 0.95, 0.9, 0.95, 0.9, 0.95, 0.9], '5neg': [0.3, 0.4, 0.5, 0.3, 0.4] };
const rows = [];
for (const mode of ['fast', 'true_flow']) for (const [n, betas] of Object.entries(cfgs)) {
  if (mode === 'true_flow' && n === '5neg') continue;
  await p.evaluate(c => window.__flow.setConfig(c), { betas, mode, level: 'coarse', seed: 'gluing' });
  await p.waitForTimeout(400);
  await p.evaluate(() => window.__flow.generate());
  await p.waitForFunction(() => !window.__flow.running(), null, { timeout: 600000 });
  await p.waitForTimeout(800);
  const r = await p.evaluate(() => ({ frames: window.__flow.frames(), wall: window.__flow.wallMs(), status: window.__flow.status(), q: window.__flow.result(), warn: document.getElementById('warn').textContent, lrr: document.getElementById('st-l').textContent }));
  rows.push({ mode, n, betas, ...r });
  await p.screenshot({ path: new URL(`${mode}_n${n}.png`, out).pathname.slice(1) });
}
// inadmissible: Luo-Tian violation
await p.evaluate(c => window.__flow.setConfig(c), { betas: [0.9, 0.9, 0.2], mode: 'fast' }); await p.waitForTimeout(1500);
rows.push({ inadmissible: await p.evaluate(() => document.getElementById('warn').textContent) });
await p.screenshot({ path: new URL('inadmissible.png', out).pathname.slice(1) });
console.log(JSON.stringify({ rows, errs }, null, 1)); await b.close();
