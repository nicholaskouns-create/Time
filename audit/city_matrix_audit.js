#!/usr/bin/env node
// Precision audit of THE MATRIX two-site MPS worker (E47-Kartekeya website/interfaces/matrix/mps-worker.js).
// Runs the worker unmodified and compares each truncation with the Schmidt-optimal one at the same bond.
//   node audit/city_matrix_audit.js [path/to/mps-worker.js]   (default: fetched from GitHub)
import { readFileSync } from 'node:fs';
import vm from 'node:vm';
import { cmat, cmul, cadj, hermEig } from '../engine/linalg.js';

const SRC_URL = 'https://raw.githubusercontent.com/nicholaskouns-create/E47-Kartekeya/main/website/interfaces/matrix/mps-worker.js';
async function fetchText(url) { const r = await fetch(url); if (!r.ok) throw Error(`${url}: HTTP ${r.status}`); return r.text(); }
const src = process.argv[2] ? readFileSync(process.argv[2], 'utf8') : await fetchText(SRC_URL);
const ctx = { postMessage() {}, Math, Float64Array, Array, Number, Object, JSON };
vm.createContext(ctx);
vm.runInContext(src, ctx);
const n = 8, layers = 24, phi = 0.7, N = 1 << n;

// exact dense reference, written independently of the worker
function step(re, im, k, chi, log) {
  const s = Math.SQRT1_2, c = Math.cos(phi), sn = Math.sin(phi);
  for (let q = 0; q < n; q++) {
    const b = 1 << (n - 1 - q);
    for (let i = 0; i < N; i++) if (!(i & b)) { const j = i | b, ar = re[i], ai = im[i], br = re[j], bi = im[j]; re[i] = s * (ar + br); im[i] = s * (ai + bi); re[j] = s * (ar - br); im[j] = s * (ai - bi); }
    for (let i = 0; i < N; i++) if (i & b) { const r = re[i], m = im[i]; re[i] = c * r - sn * m; im[i] = sn * r + c * m; }
  }
  for (let q = k % 2; q < n - 1; q += 2) {
    const a = 1 << (n - 1 - q), b = 1 << (n - 2 - q);
    for (let i = 0; i < N; i++) if ((i & a) && !(i & b)) { const j = i | b; let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
    if (chi) schmidtTruncate(re, im, q, chi, log);
  }
}
// optimal rank-chi truncation across bond (q | q+1): project onto the top chi singular vectors of the smaller side
function reshape(re, im, q) { const r = 1 << (q + 1), c = N / r, M = cmat(r, c); M.re.set(re); M.im.set(im); return M; }
function schmidtTruncate(re, im, q, chi, log) {
  const M = reshape(re, im, q), left = M.n <= M.m, G = left ? cmul(M, cadj(M)) : cmul(cadj(M), M), { values, vectors } = hermEig(G), d = G.n;
  const order = values.map((v, k) => [v, k]).sort((x, y) => y[0] - x[0]), keep = order.slice(0, chi).map(x => x[1]);
  const tail = order.slice(chi).reduce((s, x) => s + Math.max(0, x[0]), 0), total = values.reduce((s, v) => s + Math.max(0, v), 0);
  const Pk = cmat(d); for (const k of keep) for (let i = 0; i < d; i++) for (let j = 0; j < d; j++) { const ar = vectors.re[i * d + k], ai = vectors.im[i * d + k], br = vectors.re[j * d + k], bi = -vectors.im[j * d + k]; Pk.re[i * d + j] += ar * br - ai * bi; Pk.im[i * d + j] += ar * bi + ai * br; }
  const T = left ? cmul(Pk, M) : cmul(M, Pk);   // Pk = Σ v v† over the kept singular vectors (Hermitian)
  re.set(T.re); im.set(T.im);
  if (log) log.push({ bond: q, discarded_rel: total ? tail / total : 0 });
}
const run = chi => { const re = new Float64Array(N), im = new Float64Array(N), log = []; re[0] = 1; for (let k = 0; k < layers; k++) step(re, im, k, chi, log); return { re, im, log }; };
const exact = run(0);
const fidelity = (re, im) => { let nr = 0, or = 0, oi = 0; for (let i = 0; i < N; i++) { nr += re[i] ** 2 + im[i] ** 2; or += exact.re[i] * re[i] + exact.im[i] * im[i]; oi += exact.re[i] * im[i] - exact.im[i] * re[i]; } return { fidelity: (or * or + oi * oi) / nr, norm2: nr }; };

// the worker, instrumented through its own functions
const probe = vm.runInContext(`(function(n,layers,phi,chi){
  const nrm2=v=>v.reduce((s,z)=>s+z[0]*z[0]+z[1]*z[1],0); let m=init(n),steps=[],disc=0;
  for(let k=0;k<layers;k++){ for(let q=0;q<n;q++){m[q]=one(m[q],H());m[q]=one(m[q],R(phi))}
    for(let q=k%2;q<n-1;q+=2){ const before=dense(m), mm=m.slice(), [ea,eb]=two(mm[q],mm[q+1],CNOT,1<<n); mm[q]=ea; mm[q+1]=eb; const ex=dense(mm);
      const [a,b,t]=two(m[q],m[q+1],CNOT,chi); m[q]=a; m[q+1]=b; disc+=t.discarded; const after=dense(m);
      let d=0; for(let i=0;i<after.length;i++) d+=(after[i][0]-ex[i][0])**2+(after[i][1]-ex[i][1])**2;
      steps.push({bond:q,reported_rel:t.discarded,norm2_before:nrm2(before),norm2_after:nrm2(after),err2:d,exact_after:ex}); } }
  return {steps,final:dense(m),disc};
})`, ctx);

const out = { audit_code: 'MC-MATRIX-PRECISION-AUDIT-20261019-001', object: 'THE MATRIX two-site MPS worker (E47-Kartekeya)', circuit: { n, layers, phi, entangler: 'brickwork' }, runs: [] };
for (const chi of [16, 8, 4]) {
  const w = probe(n, layers, phi, chi), re = Float64Array.from(w.final.map(z => z[0])), im = Float64Array.from(w.final.map(z => z[1]));
  const fw = fidelity(re, im), trunc = w.steps.filter(s => s.reported_rel > 1e-14);
  let worse = 0, maxRatio = 0, lossRatio = [Infinity, 0];
  for (const s of trunc) {
    const eRe = Float64Array.from(s.exact_after.map(z => z[0])), eIm = Float64Array.from(s.exact_after.map(z => z[1])), lg = [];
    const absNorm = eRe.reduce((t, x, i) => t + x * x + eIm[i] * eIm[i], 0); schmidtTruncate(eRe, eIm, s.bond, chi, lg);
    const optimal = lg[0].discarded_rel * absNorm;
    if (s.err2 > optimal * (1 + 1e-9) + 1e-15) worse++;
    if (optimal > 1e-15) maxRatio = Math.max(maxRatio, s.err2 / optimal);
    const r = (s.norm2_before - s.norm2_after) / (s.reported_rel * s.norm2_before); lossRatio = [Math.min(lossRatio[0], r), Math.max(lossRatio[1], r)];
  }
  const canon = chi < 16 ? run(chi) : null, fc = canon ? fidelity(canon.re, canon.im) : null, dc = canon ? canon.log.reduce((s, x) => s + x.discarded_rel, 0) : null;
  out.runs.push({ chi, worker: { fidelity: fw.fidelity, norm2: fw.norm2, discarded_sum: w.disc, window_if_canonical: [1 - w.disc, Math.exp(-w.disc)], inside_window: fw.norm2 >= 1 - w.disc - 1e-12 && fw.norm2 <= Math.exp(-w.disc) + 1e-12 },
    truncating_steps: trunc.length, steps_worse_than_schmidt_optimum: worse, max_error_ratio_to_optimum: maxRatio, actual_over_reported_loss: trunc.length ? lossRatio : null,
    canonical_tebd: canon ? { fidelity: fc.fidelity, norm2: fc.norm2, discarded_sum: dc, inside_window: fc.norm2 >= 1 - dc - 1e-12 && fc.norm2 <= Math.exp(-dc) + 1e-12 } : null });
}
out.findings = [
  { id: 'M1', severity: 'high', title: 'Truncations are not made in canonical gauge, so the discarded weight is not the state’s truncation error and the truncation is not optimal', fix: 'Keep the MPS in Vidal Γλ form, or move the orthogonality centre to the active bond (QR sweeps) before each SVD. The untruncated χ = 16 parity certificate is unaffected.' },
  { id: 'M2', severity: 'medium', title: 'The circuit diagram draws a T1/Tφ noise box after every layer, but the worker applies the noise channel once, after the last layer', fix: 'Apply the trajectory channel per layer, or draw it once at the end.' },
  { id: 'M3', severity: 'low', title: 'In controlled-truncation mode the validation passes whenever the error is finite', fix: 'After M1, compare the fidelity with Π(1 − ε_k) instead.' },
  { id: 'M4', severity: 'low', title: 'Validation runs only for the brickwork entangler; the chain entangler has no dense reference', fix: 'Add the chain ordering to exact().' },
  { id: 'M5', severity: 'note', title: 'hermitianEig is defined but unused', fix: 'Remove it or use it.' },
];
console.log(JSON.stringify(out, null, 1));
