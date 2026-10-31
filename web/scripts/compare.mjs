// Side-by-side screenshots: the real Python /flow viewer (http://127.0.0.1:8765, started separately with
// `python -m singular_ricci.dev.viewer.app` in the curvature repo) vs this port, at a fixture's configuration.
// Usage: npm run build && node scripts/compare.mjs [fixture ...]   -> review/compare_<fixture>_<theme>.png
import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import { mkdirSync, readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join, resolve } from 'node:path';

const root = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const outDir = join(root, 'review'); mkdirSync(outDir, { recursive: true });
const names = process.argv.slice(2).length ? process.argv.slice(2) : ['n4_pos_gluing', 'n4_pos_polytope'];
const extra = process.env.OVERLAY ? [process.env.OVERLAY] : [null];
const srv = spawn(process.execPath, [join(root, 'node_modules', 'vite', 'bin', 'vite.js'), 'preview', '--host', '127.0.0.1', '--port', '4392', '--strictPort'], { cwd: root, stdio: 'ignore' });
await new Promise(r => setTimeout(r, 2500));
const browser = await chromium.launch({ args: ['--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--ignore-gpu-blocklist'] });
const opts = { viewport: { width: 1240, height: 1500 } };
const init = t => { try { localStorage.setItem('sr-theme', t); } catch (e) {} };
let code = 0;
try {
  for (const name of names) for (const theme of ['light', 'dark']) for (const ov of extra) {
    const req = JSON.parse(readFileSync(join(root, '..', 'fixtures', name, 'fixture.json'), 'utf8')).meta.request;
    const sfx = ov ? '_' + ov : '';
    // --- port
    let ctx = await browser.newContext(opts); await ctx.addInitScript(init, theme);
    let page = await ctx.newPage();
    await page.goto('http://127.0.0.1:4392/');
    await page.evaluate(r => window.__flow.setConfig({ positions: r.positions, betas: r.betas, seed: r.seed_surface, level: 'coarse' }), req);
    const t0 = Date.now();
    await page.evaluate(() => window.__flow.generate());
    await page.waitForFunction(() => window.__flow.frames() > 0 && !window.__flow.running(), null, { timeout: 280000 });
    console.log(name, theme, 'port wall ms', Date.now() - t0, JSON.stringify(await page.evaluate(() => window.__flow.result())));
    await page.waitForTimeout(4000);
    if (ov) await page.click(`#overlay input[value=${ov}]`);
    await page.waitForTimeout(800);
    const port = await page.screenshot({ fullPage: true }); await ctx.close();
    // --- python
    ctx = await browser.newContext(opts); await ctx.addInitScript(init, theme);
    page = await ctx.newPage();
    await page.goto('http://127.0.0.1:8765/');
    await page.waitForSelector('body[data-ready="1"]');
    await page.evaluate(c => window.__flow.setConfig({ positions: c.positions, betas: c.betas }), req);
    await page.evaluate(s => { const e = document.getElementById('seedsurf'); e.value = s; e.dispatchEvent(new Event('change')); }, req.seed_surface);
    await page.evaluate(() => window.__flow.setResolution('coarse'));
    await page.waitForTimeout(1500);
    await page.evaluate(() => window.__flow.generate());
    await page.waitForFunction(() => window.__flow.playCount() > 5 && /done|FAILED/i.test(window.__flow.phaseText()), null, { timeout: 280000 });
    await page.waitForTimeout(1500);
    await page.evaluate(() => window.__flow.showFrame(window.__flow.playCount() - 1));
    if (ov) await page.evaluate(o => window.__flow.overlay(o), ov === 'length' ? 'length' : 'gauss');
    await page.waitForTimeout(1500);
    const py = await page.screenshot({ fullPage: true }); await ctx.close();
    // --- compose
    ctx = await browser.newContext({ viewport: { width: 2500, height: 1600 } });
    page = await ctx.newPage();
    const b64 = b => 'data:image/png;base64,' + b.toString('base64');
    await page.setContent(`<body style="margin:0;background:#888;font:bold 18px sans-serif;color:#fff"><div style="display:flex;gap:20px;padding:10px"><div>Python (current page)<br><img src="${b64(py)}"></div><div>Port<br><img src="${b64(port)}"></div></div>`);
    await page.waitForTimeout(500);
    await page.screenshot({ path: join(outDir, `live_compare_${name}_${theme}${sfx}.png`), fullPage: true });
    await ctx.close();
    console.log('wrote', `live_compare_${name}_${theme}${sfx}.png`);
  }
} catch (e) { console.error(e); code = 1; } finally { await browser.close(); srv.kill(); }
process.exit(code);
