// True Polytope conformal-factor screenshot (n5 ring .6/.7/.5/.8/.7, r 0.7, coarse).
// Usage: npm run build && node scripts/tp-shot.mjs <tag>   -> review/tp_cf_<tag>.png (+ _disc crop)
import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import { mkdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join, resolve } from 'node:path';

const root = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const outDir = join(root, 'review'); mkdirSync(outDir, { recursive: true });
const tag = process.argv[2] || 'after';
const seed = process.argv[3] || 'true_polytope';
const srv = spawn(process.execPath, [join(root, 'node_modules', 'vite', 'bin', 'vite.js'), 'preview', '--host', '127.0.0.1', '--port', '4396', '--strictPort'], { cwd: root, stdio: 'ignore' });
await new Promise(r => setTimeout(r, 2500));
const browser = await chromium.launch({ args: ['--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--ignore-gpu-blocklist'] });
let code = 0;
try {
  const positions = [...Array(5)].map((_, k) => { const a = Math.PI / 2 + 2 * Math.PI * k / 5; return [0.7 * Math.cos(a), 0.7 * Math.sin(a)]; });
  const ctx = await browser.newContext({ viewport: { width: 1240, height: 1500 } });
  const page = await ctx.newPage();
  page.on('pageerror', e => { console.error('pageerror', e.message); code = 1; });
  await page.goto('http://127.0.0.1:4396/');
  await page.evaluate(c => window.__flow.setConfig(c), { positions, betas: [0.6, 0.7, 0.5, 0.8, 0.7], seed, level: 'coarse' });
  await page.evaluate(() => window.__flow.generate());
  await page.waitForFunction(() => window.__flow.frames() > 0 && !window.__flow.running(), null, { timeout: 280000 });
  await page.waitForTimeout(5000);
  console.log(JSON.stringify(await page.evaluate(() => window.__flow.result())));
  await page.screenshot({ path: join(outDir, `tp_cf_${tag}.png`), fullPage: true });
  const el = await page.$('#cf-n');
  if (el) await (await el.evaluateHandle(e => e.parentElement.parentElement)).asElement().screenshot({ path: join(outDir, `tp_cf_${tag}_disc.png`) });
  await ctx.close();
} catch (e) { console.error(e); code = 1; } finally { await browser.close(); srv.kill(); }
process.exit(code);
