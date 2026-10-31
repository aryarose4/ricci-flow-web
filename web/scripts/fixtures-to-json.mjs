// Converts ../fixtures/*/{arrays.npz,fixture.json} into public/fixtures/*.json (+ index.json).
import { readFileSync, writeFileSync, mkdirSync, readdirSync, existsSync } from 'node:fs';
import { unzipSync } from 'fflate';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';

const here = dirname(fileURLToPath(import.meta.url));
const src = join(here, '..', '..', 'fixtures');
const out = join(here, '..', 'public', 'fixtures');
mkdirSync(out, { recursive: true });

function parseNpy(buf) {
  const dv = new DataView(buf.buffer, buf.byteOffset, buf.byteLength);
  const major = buf[6];
  const hl = major === 1 ? dv.getUint16(8, true) : dv.getUint32(8, true);
  const off = major === 1 ? 10 : 12;
  const hdr = new TextDecoder().decode(buf.subarray(off, off + hl));
  const descr = /'descr':\s*'([^']+)'/.exec(hdr)[1];
  const shape = (/'shape':\s*\(([^)]*)\)/.exec(hdr)[1]).split(',').map(s => s.trim()).filter(Boolean).map(Number);
  if (/'fortran_order':\s*True/.test(hdr)) throw new Error('fortran order unsupported');
  const n = shape.reduce((a, b) => a * b, 1);
  const body = buf.slice(off + hl); // copy -> aligned
  let data;
  if (descr === '<f8') data = new Float64Array(body.buffer, 0, n);
  else if (descr === '<i8') data = Array.from(new BigInt64Array(body.buffer, 0, n), Number);
  else if (descr === '<i4') data = new Int32Array(body.buffer, 0, n);
  else if (descr === '<f4') data = new Float32Array(body.buffer, 0, n);
  else throw new Error('unsupported dtype ' + descr);
  return { shape, data };
}
const r = (x, d = 6) => (Number.isFinite(x) ? Math.round(x * 10 ** d) / 10 ** d : null);
const rows = (a, k) => { const w = a.data.length / a.shape[0]; return Array.from({ length: a.shape[0] }, (_, i) => Array.from(a.data.subarray ? a.data.subarray(i * w, (i + 1) * w) : a.data.slice(i * w, (i + 1) * w))); };
const flat = (x, d) => Array.from(x, v => r(v, d));

const names = readdirSync(src, { withFileTypes: true }).filter(e => e.isDirectory() && existsSync(join(src, e.name, 'arrays.npz'))).map(e => e.name).sort();
const index = [];
for (const name of names) {
  const files = unzipSync(new Uint8Array(readFileSync(join(src, name, 'arrays.npz'))));
  const A = {};
  for (const [k, v] of Object.entries(files)) A[k.replace(/\.npy$/, '')] = parseNpy(v);
  const fx = JSON.parse(readFileSync(join(src, name, 'fixture.json'), 'utf8'));
  const nf = A.t.data.length, nv = A.sphere_vertices.shape[0];
  const faceSets = {};
  for (const k of Object.keys(A)) if (/^faces3d_\d+$/.test(k)) faceSets[k.slice(8)] = Array.from(A[k].data);
  const frames = [];
  for (let k = 0; k < nf; k++) {
    const sl = (key, w) => A[key].data.slice(k * w, (k + 1) * w);
    const sc = key => A[key].data[k];
    frames.push({
      t: r(sc('t'), 9), faceSet: String(A.face_idx.data[k]),
      vertices: flat(sl('vertices3d', nv * 3), 5), u: flat(sl('u', nv), 6),
      K: flat(Array.from(sl('K_defect', nv), (d, i) => d / sl('A_vertex', nv)[i]), 5),
      Kdef: flat(sl('K_defect', nv), 6), Avert: flat(sl('A_vertex', nv), 8),
      cfW: A.cf_w ? flat(sl('cf_w', A.cf_w.shape[1] * A.cf_w.shape[2]), 4) : null,
      L: flat(sl('L', A.L.shape[1]), 6),
      gbError: sc('gb_resid'), lenRelRms: r(sc('len_rel_rms'), 8), lenRelMax: r(sc('len_rel_max'), 8),
      coneDefectErr: sc('cone_defect_err'), nDeep: sc('n_deep'), rho: r(sc('rho_frame'), 8),
    });
  }
  const m = fx.meta;
  const doc = {
    name, request: m.request, rho: fx.scalars.rho, quality: fx.scalars.quality, label: fx.scalars.label,
    embedNotice: fx.scalars.embed_notice ?? null, embedFailure: fx.scalars.embed_failure ?? null,
    coneIds: Array.from(A.cone_ids.data), sphereVertices: flat(A.sphere_vertices.data, 6),
    edges: Array.from(A.edges.data), faceSets, frames,
    // smooth_w raster per frame (cf_w): ny rows over colatitude [0, pi], nx columns over longitude [-pi, pi]
    cfGrid: A.cf_w ? { ny: A.cf_w.shape[1], nx: A.cf_w.shape[2], extent: [-Math.PI, Math.PI, 0, Math.PI] } : null,
  };
  writeFileSync(join(out, name + '.json'), JSON.stringify(doc));
  index.push({ name, n: m.request.n, betas: m.request.betas, seed: m.request.seed_surface, rho: fx.scalars.rho, quality: fx.scalars.quality.label });
}
writeFileSync(join(out, 'index.json'), JSON.stringify(index));
console.log('wrote', names.length, 'fixtures to', out);
