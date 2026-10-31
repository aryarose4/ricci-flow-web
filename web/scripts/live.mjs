// Live acceptance: runs the wasm solver in headless Chromium for fixture configurations at Coarse and compares the
// final rho / quality label with the golden fixture. Writes review/live_<fixture>_<theme>.png.
// Usage: npm run build && node scripts/live.mjs [fixture ...]
import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import { mkdirSync, readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join, resolve } from 'node:path';

const root = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const outDir = join(root, 'review'); mkdirSync(outDir, { recursive: true });
const names = process.argv.slice(2).length ? process.argv.slice(2) : ['n4_pos_gluing', 'n4_pos_polytope'];
const srv = spawn(process.execPath, [join(root, 'node_modules', 'vite', 'bin', 'vite.js'), 'preview', '--host', '127.0.0.1', '--port', '4393', '--strictPort'], { cwd: root, stdio: 'ignore' });
await new Promise(r => setTimeout(r, 2500));
const browser = await chromium.launch({ args: ['--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--ignore-gpu-blocklist'] });
let code = 0;
try {
  for (const name of names) for (const theme of ['light', 'dark']) {
    const fx = JSON.parse(readFileSync(join(root, '..', 'fixtures', name, 'fixture.json'), 'utf8'));
    const req = fx.meta.request;
    const ctx = await browser.newContext({ viewport: { width: 1240, height: 1500 } });
    await ctx.addInitScript(t => { try { localStorage.setItem('sr-theme', t); } catch (e) {} }, theme);
    const page = await ctx.newPage();
    page.on('pageerror', e => { console.error('pageerror', e.message); code = 1; });
    await page.goto('http://127.0.0.1:4393/');
    await page.evaluate(r => window.__flow.setConfig({ positions: r.positions, betas: r.betas, seed: r.seed_surface, level: 'coarse' }), req);
    const t0 = Date.now();
    await page.evaluate(() => window.__flow.generate());
    await page.waitForFunction(() => window.__flow.frames() > 0 && !window.__flow.running(), null, { timeout: 280000 });
    const wall = Date.now() - t0;
    await page.waitForTimeout(4000); // smooth playback settles
    const res = await page.evaluate(() => window.__flow.result());
    const ok = res && res.label.split(' (')[0] === fx.scalars.quality.label.split(' (')[0] && Math.abs(res.rho - fx.scalars.rho) < 0.03 * Math.abs(fx.scalars.rho);
    console.log(JSON.stringify({ name, theme, wall_ms: wall, frames: await page.evaluate(() => window.__flow.frames()), live: res && { rho: res.rho, label: res.label }, fixture: { rho: fx.scalars.rho, label: fx.scalars.quality.label }, match: ok }));
    if (!ok) code = 1;
    await page.screenshot({ path: join(outDir, `live_${name}_${theme}.png`), fullPage: true });
    await ctx.close();
  }
} catch (e) { console.error(e); code = 1; } finally { await browser.close(); srv.kill(); }
process.exit(code);
