#!/usr/bin/env node
// Precision audit of The Matrix Quantum Simulator (matrix-quantum-simulator/ on nicholaskouns-create.github.io).
// Runs the app's own engine and receipt code, unmodified, in a Node vm with a minimal DOM stub.
//   node audit/site_matrix_audit.js [path/to/matrix-quantum-simulator]   (default: its current source, fetched from GitHub)
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import vm from 'node:vm';
import { webcrypto, createHash } from 'node:crypto';

const SRC_URL = 'https://raw.githubusercontent.com/nicholaskouns-create/nicholaskouns-create.github.io/main/matrix-quantum-simulator/';
const appDir = process.argv[2] ? resolve(process.argv[2]) : null;
async function source(name) {
  if (appDir) return readFileSync(resolve(appDir, name), 'utf8');
  const r = await fetch(SRC_URL + name);
  if (!r.ok) throw Error(`${SRC_URL + name}: HTTP ${r.status}`);
  return r.text();
}
const html = await source('app.html'), receiptSrc = await source('receipt.js');
const sha256 = s => createHash('sha256').update(s, 'utf8').digest('hex');
const script = html.slice(html.indexOf('const C=(r=0,i=0)'), html.indexOf("draw();render();log('> MATRIX QUANTUM ENGINE ONLINE')"));

function boot() {
  const els = new Map(), el = id => { if (!els.has(id)) els.set(id, { id, textContent: '', innerHTML: '', value: '', className: '', classList: { add() {}, remove() {} }, append() {}, appendChild() {}, onclick: null, onchange: null, click() { this.onclick?.(); }, files: [], setAttribute() {}, removeAttribute() {}, scrollTop: 0, scrollHeight: 0, dataset: {} }); return els.get(id); };
  const document = { getElementById: el, createElement: () => ({ ...el('tmp' + Math.random()), click() {} }), querySelectorAll: () => [] };
  const ctx = { document, console, Math, Array, Object, JSON, Number, Set, Map, Error, TextEncoder, crypto: webcrypto, URL: { createObjectURL: () => '', revokeObjectURL() {} }, Blob: class {}, navigator: {}, setTimeout, window: {} };
  vm.createContext(ctx);
  vm.runInContext(script + ';globalThis.engine={get ops(){return ops},set ops(v){ops=v},get psi(){return psi},set psi(v){psi=v},setPick(g){pick=g},place,run,P,C};', ctx);
  vm.runInContext(receiptSrc, ctx);
  ctx.MatrixReceipt = ctx.window.MatrixReceipt;
  return ctx;
}
const amp = psi => psi.map(z => [z.r, z.i]);
const dist = (a, b) => Math.sqrt(a.reduce((s, z, k) => s + (z[0] - b[k][0]) ** 2 + (z[1] - b[k][1]) ** 2, 0));
const findings = [];
const finding = (id, title, severity, evidence, fix, open = true) => findings.push({ id, title, severity, status: open ? 'open' : 'resolved', evidence, fix });

// F1: the importer rejects the exporter's receipts
{
  const ctx = boot(), E = ctx.engine, out = {};
  for (const [name, ops] of Object.entries({ empty: [], bell: [{ g: 'H', a: 0, b: null, col: 0 }, { g: 'CX', a: 0, b: 1, col: 1 }] })) {
    E.ops = ops; E.run();
    const r = await ctx.MatrixReceipt.export();
    try { await ctx.MatrixReceipt.import(r); out[name] = 'imported'; } catch (e) { out[name] = 'rejected: ' + e.message; }
  }
  finding('F1', 'Receipt import rejects every receipt the app exports with at least one gate', 'high',
    { import_results: out, cause: 'valid() checks o.g, but receipt() writes each gate as o.gate; read() then uses o.gate' },
    'In receipt.js valid(), read the gate as o.gate (the two occurrences of o.g).', out.bell !== 'imported');
}
// F2: execution order is placement order, not diagram order
{
  const ctx = boot(), E = ctx.engine;
  E.ops = [{ g: 'H', a: 0, b: null, col: 0 }, { g: 'CX', a: 0, b: 1, col: 1 }]; E.run(); const bell = amp(E.psi);
  E.setPick('H'); E.place(0, 0); E.run(); const replaced = amp(E.psi);
  E.ops = []; E.setPick('S'); E.place(0, 1); E.setPick('H'); E.place(0, 0); E.run(); const hs = amp(E.psi);
  const expectHS = [[Math.SQRT1_2, 0], [0, 0], [0, 0], [0, 0], [0, Math.SQRT1_2], [0, 0], [0, 0], [0, 0]];
  finding('F2', 'run() applies gates in the order they were placed, not the column order the circuit shows', 'high',
    { bell_then_replace_H: { before: bell, after: replaced, distance: dist(bell, replaced) }, H_col0_S_col1: { computed: hs, diagram_says: expectHS, distance: dist(hs, expectHS) } },
    'In run(), iterate over [...ops].sort((x, y) => x.col - y.col); gates sharing a column act on different wires and commute.', dist(bell, replaced) > 1e-12 || dist(hs, expectHS) > 1e-12);
}
// F3: a two-qubit gate clears its whole column
{
  const ctx = boot(), E = ctx.engine;
  E.ops = []; E.setPick('X'); E.place(2, 3); const before = JSON.parse(JSON.stringify(E.ops));
  E.setPick('CX'); E.place(0, 3); E.place(1, 3); const after = JSON.parse(JSON.stringify(E.ops));
  finding('F3', 'Placing a two-qubit gate deletes every gate in that column, including on the untouched wire', 'medium',
    { before, after }, 'In place(), remove only gates in that column that touch the two chosen qubits.', !after.some(o => o.g === 'X' && o.a === 2));
}
// F4: export can pair a circuit with a stale statevector
{
  const ctx = boot(), E = ctx.engine;
  E.ops = []; E.run(); E.setPick('H'); E.place(0, 0);
  const r = await ctx.MatrixReceipt.export();
  finding('F4', 'Placing a gate does not recompute the state, so EXPORT can pair a circuit with an older statevector', 'medium',
    { circuit: r.circuit, statevector: r.statevector, note: 'H is in the circuit; the statevector is still |000>' },
    'Call run() at the start of receipt(), or recompute psi in place().', Math.abs(r.statevector[0][0] - 1) < 1e-12);
}
// F5: measurement sampling round-off edge
{
  const ctx = boot(), E = ctx.engine;
  const shot = html.slice(html.indexOf("document.getElementById('shot').onclick=()=>{"), html.indexOf("draw();render();log('> MATRIX QUANTUM ENGINE ONLINE')"));
  finding('F5', 'MEASURE falls through to |111> when round-off leaves the cumulative probability below the random draw', 'low',
    { code: shot.slice(0, 140), note: 'with Σ|ψ|² = 1 − 1e−16 and u above it, k reaches 7 whatever its probability' },
    'Stop at the last basis state with nonzero probability.', shot.includes('while(k<7&&(u-=P(psi[k]))>=0)k++'));
}
// F6: digest is integrity, not authenticity
finding('F6', 'The SHA-256 digest detects accidental edits but not deliberate ones: editing amplitudes and recomputing the digest passes import', 'note',
  { note: 'The Quantum Clock bridge replays the circuit to catch this (status REPLAY-MISMATCH).' }, 'Replay the circuit on import and compare to the statevector.');

const report = { audit_code: 'MC-MATRIX-PRECISION-AUDIT-20261019-001', object: 'The Matrix Quantum Simulator (site app)', path: 'matrix-quantum-simulator/', date: '2026-10-19', method: 'the app’s engine and receipt code executed unmodified in a Node vm',
  source: appDir || SRC_URL, sha256: { 'app.html': sha256(html), 'receipt.js': sha256(receiptSrc) }, findings };
console.log(JSON.stringify(report, null, 1));
