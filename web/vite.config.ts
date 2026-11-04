import { defineConfig } from 'vite';
// base './' so the build works under any GitHub Pages subpath.
// fs.allow '..' lets the dev server serve ../core/pkg (the wasm build).
export default defineConfig({
  base: './',
  worker: { format: 'es' },
  build: { target: 'es2022' },
  server: { fs: { allow: ['..'] } },
});
