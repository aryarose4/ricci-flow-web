// Measures real per-true-frame display intervals across a full replay.
// Usage: node scripts/playback-timing.mjs <port> [label]   (dev server already running)  -> prints JSON + table
import { createRequire } from 'node:module';
const require = createRequire(new URL('../package.json', import.meta.url));
const { chromium } = require('playwright');
const port = process.argv[2] || '4391', label = process.argv[3] || 'run';
const levels = (process.env.LEVELS || 'coarse,medium').split(','), ovs = (process.env.OVS || 'gauss,length').split(','), smooths = [true, false];
const browser = await chromium.launch({ args: ['--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--ignore-gpu-blocklist'] });
const errs = [], out = [];
for (const level of levels) {
  const page = await (await browser.newContext({ viewport: { width: 1240, height: 1000 } })).newPage();
  page.on('pageerror', e => errs.push(String(e))); page.on('console', m => { if (m.type() === 'error') errs.push(m.text()); });
  await page.goto(`http://127.0.0.1:${port}/`);
  await page.waitForFunction(() => window.__flow && window.__flow.previewReady());
  await page.evaluate(l => window.__flow.setConfig({ level: l, smooth: true }), level);
  await page.waitForFunction(() => window.__flow.previewReady());
  await page.evaluate(() => window.__flow.generate());
  await page.waitForFunction(() => !window.__flow.running() && window.__flow.frames() > 0, null, { timeout: 600000 });
  const nF = await page.evaluate(() => window.__flow.frames());
  for (const ov of ovs) for (const smooth of smooths) {
    await page.evaluate(([o, s]) => { document.querySelector(`#fpsmooth input[value="${s ? 'on' : 'off'}"]`).click(); const r = document.querySelector(`#overlay input[value="${o}"]`); r.click(); window.__flow.seek(0); }, [ov, smooth]);
    await page.evaluate(() => { window.__ev = []; window.__onDisplay = (t, k, sub) => window.__ev.push([t, k, sub ? 1 : 0]); });
    await page.waitForTimeout(500); await page.evaluate(() => { window.__ev = []; }); await page.click('#playbtn');
    await page.waitForFunction(() => document.querySelector('#playbtn').textContent === 'Play', null, { timeout: 900000, polling: 200 });
    const st = await page.evaluate(() => window.__flow.playStats ? window.__flow.playStats() : null);
    const ev = await page.evaluate(() => { window.__onDisplay = null; return window.__ev; });
    const arr = {}; let subs = 0; for (const [t, k, s] of ev) { if (s) subs++; else if (!(k in arr)) arr[k] = t; }
    const ks = Object.keys(arr).map(Number).sort((a, b) => a - b), d = [];
    for (let i = 1; i < ks.length; i++) if (ks[i] === ks[i - 1] + 1 && arr[ks[i]] > arr[ks[i - 1]]) d.push(arr[ks[i]] - arr[ks[i - 1]]);
    const half = Math.floor(d.length / 2), m = a => a.length ? a.reduce((x, y) => x + y, 0) / a.length : NaN, mx = a => Math.max(...a);
    // render-time proxy: gaps between consecutive displays
    const gaps = []; for (let i = 1; i < ev.length; i++) gaps.push(ev[i][0] - ev[i - 1][0]);
    const av = a => a && a.length ? +(a.reduce((x, y) => x + y, 0) / a.length).toFixed(1) : null, mxx = a => a && a.length ? +Math.max(...a).toFixed(1) : null;
    const r = { label, subCostMean: st && av(st.dispSub), subCostMax: st && mxx(st.dispSub), planMax: st && mxx(st.planMs), nEvalMax: st && mxx(st.nEval), level, ov, smooth, nF, trueShown: ks.length, subs, meanAll: m(d), maxAll: mx(d), meanEarly: m(d.slice(0, half)), meanLate: m(d.slice(half)), maxEarly: mx(d.slice(0, half)), maxLate: mx(d.slice(half)), maxGap: mx(gaps), total: ev.length ? ev[ev.length - 1][0] - ev[0][0] : 0 };
    out.push(r); console.log(JSON.stringify(Object.fromEntries(Object.entries(r).map(([k, v]) => [k, typeof v === 'number' ? +v.toFixed(1) : v]))));
  }
  await page.context().close();
}
await browser.close();
console.log('ERRORS', JSON.stringify(errs));
