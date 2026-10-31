import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
const srv = spawn(process.execPath, ['node_modules/vite/bin/vite.js','preview','--host','127.0.0.1','--port','4395','--strictPort'], { stdio: 'ignore' });
await new Promise(r => setTimeout(r, 2500));
const b = await chromium.launch({ args: ['--use-angle=swiftshader','--enable-unsafe-swiftshader','--ignore-gpu-blocklist'] });
try {
  const p = await b.newPage({ viewport: { width: 1240, height: 1500 } });
  p.on('pageerror', e => console.log('PAGEERROR', e.message));
  await p.goto('http://127.0.0.1:4395/');
  const ready = () => p.waitForFunction(() => window.__flow.previewReady(), null, { timeout: 60000 }).then(() => p.waitForTimeout(400));
  await ready(); await p.screenshot({ path: 'review/preview_load.png' });
  console.log('preview ms (coarse, load)', await p.evaluate(() => window.__flow.previewMs()));
  for (const sd of ['smart', 'polytope', 'gluing']) {
    await p.evaluate(s => window.__flow.setConfig({ seed: s }), sd); await ready();
    await p.screenshot({ path: `review/preview_seed_${sd}.png` });
  }
  // mid-drag: move cone 0 in steps, shoot before release
  const pos = await p.evaluate(() => window.__flow.positions());
  for (let k = 1; k <= 6; k++) {
    await p.evaluate(([x, y]) => window.__flow.dragCone(0, [x, y], false), [pos[0][0] + 0.08 * k, pos[0][1] + 0.05 * k]);
    await p.waitForTimeout(120);
  }
  await ready(); await p.screenshot({ path: 'review/preview_middrag.png' });
  console.log('mid-drag preview ms', await p.evaluate(() => window.__flow.previewMs()));
  await p.evaluate(() => window.__flow.dragCone(0, window.__flow.positions()[0], true)); await ready();
  await p.evaluate(() => window.__flow.setConfig({ betas: [0.8, 0.8, 0.8, 0.8, 0.6] })); await ready();
  await p.screenshot({ path: 'review/preview_n5.png' });
  const run = async tag => {
    await p.click('#gen');
    await p.waitForFunction(() => !window.__flow.running(), null, { timeout: 180000 });
    await p.waitForTimeout(4000);
    console.log(tag, await p.evaluate(() => ({ f: window.__flow.frames(), s: window.__flow.status(), pl: document.getElementById('playlbl').textContent, t: document.getElementById('tval').textContent })));
    await p.screenshot({ path: `review/preview_flow_${tag}.png` });
  };
  await run('first');
  await p.click('#reset'); await ready(); await p.screenshot({ path: 'review/preview_after_reset.png' });
  await run('second');
  await p.reload(); await ready(); await p.screenshot({ path: 'review/preview_after_reload.png' });
  await run('third');
} finally { await b.close(); srv.kill(); }
