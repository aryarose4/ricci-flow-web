import { defineConfig } from 'vite';
// base './' so the build works under any GitHub Pages subpath.
export default defineConfig({ base: './', worker: { format: 'es' }, build: { target: 'es2022' } });
