#!/usr/bin/env node
// Quantum Clock command line: the same engine the page runs, headless.
//   node tools/qclock-cli.js certificate [--out data/certificate.json]
//   node tools/qclock-cli.js sweep --param gamma --from 5 --to 80 --steps 16 [--c 1 --omega 2 --clock linear] [--csv]
//   node tools/qclock-cli.js run --tau 50 [--omega 0.785 --c 1 --gamma 20 --clock matched --eta 1 --v 0.3,0,0]
//   node tools/qclock-cli.js verify-receipt path/to/receipt.json
import { readFileSync, writeFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { runCertificate } from '../engine/certificate.js';
import { sweep } from '../engine/automate.js';
import { NodeClock } from '../engine/sim.js';
import { verifyMatrixReceipt, verifyClockReceipt, MATRIX_SCHEMA } from '../engine/receipt.js';

const here = dirname(fileURLToPath(import.meta.url)), root = resolve(here, '..');
const [cmd, ...rest] = process.argv.slice(2);
const opt = {}; for (let k = 0; k < rest.length; k++) if (rest[k].startsWith('--')) { const key = rest[k].slice(2), val = rest[k + 1] && !rest[k + 1].startsWith('--') ? rest[++k] : true; opt[key] = val; } else (opt._ ||= []).push(rest[k]);
const num = (x, d) => x === undefined ? d : Number(x);

if (cmd === 'certificate') {
  const golden = JSON.parse(readFileSync(resolve(root, 'data/golden.json'), 'utf8'));
  const cert = await runCertificate({ golden, environment: { runtime: `node ${process.version}`, platform: process.platform }, onProgress: (k, t) => process.stderr.write(`\r${String(k).padStart(2)} ${t.slice(0, 90).padEnd(90)}`) });
  process.stderr.write('\n');
  for (const c of cert.checks) console.log(`${c.pass ? 'PASS' : 'FAIL'}  ${c.id.padEnd(3)} ${c.title}${c.residual !== null && typeof c.residual === 'number' ? `  [residual ${c.residual.toExponential(2)}]` : ''}`);
  console.log(`${cert.certificate_code}: ${cert.status} in ${cert.elapsed_ms} ms · digest ${cert.digest}`);
  const out = opt.out ? resolve(process.cwd(), opt.out) : resolve(root, 'data/certificate.json');
  writeFileSync(out, JSON.stringify(cert, null, 1) + '\n');
  process.exit(cert.passed === cert.total ? 0 : 1);
} else if (cmd === 'sweep') {
  const res = sweep({ param: opt.param || 'gamma', from: num(opt.from, 5), to: num(opt.to, 80), steps: num(opt.steps, 16), base: { omega: num(opt.omega, 2), c: num(opt.c, 1), gamma: num(opt.gamma, 20), eta: num(opt.eta, 1), clock: opt.clock || 'linear' } });
  if (opt.csv) { console.log(res.columns.join(',')); for (const r of res.rows) console.log(res.columns.map(k => r[k]).join(',')); } else console.log(JSON.stringify(res, null, 1));
} else if (cmd === 'run') {
  const v = (opt.v || '0,0,0').split(',').map(Number);
  const clk = new NodeClock({ omega: num(opt.omega, 2 * Math.PI / 8), c: num(opt.c, 1), gamma: num(opt.gamma, 20), clock: opt.clock || 'matched', eta: num(opt.eta, 1), velocity: v, seed: num(opt.seed, 47125) });
  const tau = num(opt.tau, 50); clk.advance(tau / clk.N);
  console.log(JSON.stringify(clk.state(), null, 1));
} else if (cmd === 'verify-receipt') {
  const x = JSON.parse(readFileSync(opt._?.[0], 'utf8'));
  const r = x.schema === MATRIX_SCHEMA ? await verifyMatrixReceipt(x) : await verifyClockReceipt(x);
  if (r.M) delete r.M;
  console.log(JSON.stringify(r, null, 1));
  process.exit(r.status ? (['REPLAY-VERIFIED', 'MEASURED-CONSISTENT'].includes(r.status) ? 0 : 2) : (r.ok && r.digestOk ? 0 : 2));
} else {
  console.log(readFileSync(fileURLToPath(import.meta.url), 'utf8').split('\n').slice(1, 7).map(l => l.replace(/^\/\/ ?/, '')).join('\n'));
  process.exit(cmd ? 1 : 0);
}
