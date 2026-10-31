import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
const srv = spawn(process.execPath, ['node_modules/vite/bin/vite.js','preview','--host','127.0.0.1','--port','4396','--strictPort'], { stdio: 'ignore' });
await new Promise(r => setTimeout(r, 2500));
const b = await chromium.launch({ args: ['--use-angle=swiftshader','--enable-unsafe-swiftshader','--ignore-gpu-blocklist'] });
try {
  const p = await b.newPage();
  await p.goto('http://127.0.0.1:4396/');
  await p.click('#gen'); await p.waitForFunction(() => !window.__flow.running(), null, { timeout: 180000 });
  for (const k of [2, 8, 14]) console.log('k', k, await p.evaluate(k => window.__flow.faceTProbe(k, 10), k));
} finally { await b.close(); srv.kill(); }
