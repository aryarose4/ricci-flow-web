// True Polytope preview latency (slider / drag path) for n4 / n5 at Coarse.
// Usage: npm run build && node scripts/tp-preview.mjs
import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { dirname, join, resolve } from 'node:path';

const root = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const srv = spawn(process.execPath, [join(root, 'node_modules', 'vite', 'bin', 'vite.js'), 'preview', '--host', '127.0.0.1', '--port', '4396', '--strictPort'], { cwd: root, stdio: 'ignore' });
await new Promise(r => setTimeout(r, 2500));
const browser = await chromium.launch({ args: ['--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--ignore-gpu-blocklist'] });
const ring = (n, r) => [...Array(n)].map((_, k) => { const a = Math.PI / 2 + 2 * Math.PI * k / n; return [r * Math.cos(a), r * Math.sin(a)]; });
const cfgs = [
  ['n4', [[0.7, 0], [-0.7, 0], [0, 0.7], [0, -0.7]], [0.6, 0.7, 0.5, 0.8]],
  ['n5', ring(5, 0.7), [0.6, 0.7, 0.5, 0.8, 0.7]],
];
let code = 0;
try {
  const page = await browser.newPage({ viewport: { width: 1240, height: 1500 } });
  page.on('pageerror', e => { console.error('pageerror', e.message); code = 1; });
  await page.goto('http://127.0.0.1:4396/');
  for (const seed of ['gluing', 'true_polytope']) for (const [name, positions, betas] of cfgs) {
    await page.evaluate(c => window.__flow.setConfig(c), { positions, betas, seed, level: 'coarse' });
    await page.waitForTimeout(400);
    await page.waitForFunction(() => window.__flow.previewReady(), null, { timeout: 60000 });
    const ms = [await page.evaluate(() => window.__flow.previewMs())];
    for (let k = 1; k <= 4; k++) {
      const p = [positions[0][0] + 0.03 * k, positions[0][1] - 0.02 * k];
      await page.evaluate(([p, fin]) => window.__flow.dragCone(0, p, fin), [p, k === 4]);
      await page.waitForTimeout(300);
      await page.waitForFunction(() => window.__flow.previewReady(), null, { timeout: 60000 });
      ms.push(await page.evaluate(() => window.__flow.previewMs()));
    }
    console.log(JSON.stringify({ seed, name, preview_ms: ms.map(x => Math.round(x)) }));
  }
} catch (e) { console.error(e); code = 1; } finally { await browser.close(); srv.kill(); }
process.exit(code);
