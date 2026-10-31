import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
const srv = spawn(process.execPath, ['node_modules/vite/bin/vite.js','preview','--host','127.0.0.1','--port','4394','--strictPort'], { stdio: 'ignore' });
await new Promise(r => setTimeout(r, 2500));
const b = await chromium.launch({ args: ['--use-angle=swiftshader','--enable-unsafe-swiftshader','--ignore-gpu-blocklist'] });
try {
  const p = await b.newPage({ viewport: { width: 1240, height: 1500 } });
  p.on('console', m => console.log('console:', m.type(), m.text()));
  p.on('pageerror', e => console.log('PAGEERROR', e.message));
  await p.goto('http://127.0.0.1:4394/');
  const run = async tag => {
    await p.click('#gen');
    await p.waitForFunction(() => !window.__flow.running(), null, { timeout: 120000 });
    console.log(tag, 'frames', await p.evaluate(() => window.__flow.frames()), 'status', await p.evaluate(() => window.__flow.status()));
    await p.screenshot({ path: `review/repro_${tag}.png` });
  };
  await run('first');
  await p.click('#reset'); await p.waitForTimeout(500); await p.screenshot({ path: 'review/repro_after_reset.png' });
  await run('second');
  await p.reload(); await p.waitForTimeout(1500); await p.screenshot({ path: 'review/repro_after_reload.png' });
  await run('third');
} finally { await b.close(); srv.kill(); }
