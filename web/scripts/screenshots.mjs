// Plays back every fixture in headless Chromium and writes review/<fixture>_<theme>.png.
// Usage: npm run build && npm run screenshots   (starts its own `vite preview`).
import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import { mkdirSync, readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join, resolve } from 'node:path';

const root = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const outDir = join(root, 'review');
mkdirSync(outDir, { recursive: true });
const names = JSON.parse(readFileSync(join(root, 'public', 'fixtures', 'index.json'), 'utf8')).map(e => e.name);

const srv = spawn(process.execPath, [join(root, 'node_modules', 'vite', 'bin', 'vite.js'), 'preview', '--host', '127.0.0.1', '--port', '4391', '--strictPort'], { cwd: root, stdio: 'ignore' });
await new Promise(r => setTimeout(r, 2500));
let code = 0;
try {
  const browser = await chromium.launch({ args: ['--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--ignore-gpu-blocklist'] });
  for (const theme of ['light', 'dark']) {
    const ctx = await browser.newContext({ viewport: { width: 1240, height: 1500 } });
    await ctx.addInitScript(t => { try { localStorage.setItem('sr-theme', t); } catch (e) {} }, theme);
    const page = await ctx.newPage();
    page.on('pageerror', e => { console.error('pageerror', e.message); code = 1; });
    await page.goto('http://127.0.0.1:4391/');
    await page.waitForFunction(() => document.querySelectorAll('#fixture option').length > 1);
    for (const name of names) {
      await page.evaluate(n => window.__flow.setFixture(n), name);
      await page.click('#gen');
      await page.waitForFunction(() => window.__flow.frames() > 0 && !window.__flow.running(), null, { timeout: 60000 });
      const nfr = await page.evaluate(() => window.__flow.frames());
      const q = await page.textContent('#qnote');
      console.log(theme, name, 'frames', nfr, 'quality', q);
      if (nfr !== 21) code = 1;
      await page.screenshot({ path: join(outDir, `${name}_${theme}.png`), fullPage: true });
    }
    await ctx.close();
  }
  await browser.close();
} finally { srv.kill(); }
process.exit(code);
