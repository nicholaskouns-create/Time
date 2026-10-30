// node --test tests/engine.test.js   (Node 22+, no dependencies)
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { runCertificate, CERTIFICATE_ID } from '../engine/certificate.js';
import { NodeClock } from '../engine/sim.js';
import { Network, xorshift } from '../engine/network.js';
import { sweep, encodeHash, decodeHash } from '../engine/automate.js';
import * as Rc from '../engine/receipt.js';
import { radixString } from '../engine/clocks.js';

const root = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const golden = JSON.parse(readFileSync(resolve(root, 'data/golden.json'), 'utf8'));

test('precision certificate: every check passes', { timeout: 120000 }, async () => {
  const cert = await runCertificate({ golden, environment: { runtime: `node ${process.version}` } });
  const failed = cert.checks.filter(c => !c.pass).map(c => `${c.id} ${c.title}: ${JSON.stringify(c.value)}`);
  assert.deepEqual(failed, []);
  assert.equal(cert.certificate_code, CERTIFICATE_ID);
  assert.ok(cert.total >= 40, `only ${cert.total} checks`);
});

test('committed certificate matches a fresh run (same checks, same verdicts)', { timeout: 120000 }, async () => {
  const committed = JSON.parse(readFileSync(resolve(root, 'data/certificate.json'), 'utf8'));
  const fresh = await runCertificate({ golden });
  assert.deepEqual(fresh.checks.map(c => [c.id, c.pass]), committed.checks.map(c => [c.id, c.pass]));
  assert.equal(fresh.status, committed.status);
});

test('running clock: matched hand reads proper time on a moving worldline, through clicks', () => {
  const clk = new NodeClock({ omega: 2 * Math.PI / 8, c: 1, gamma: 20, clock: 'matched', eta: 1, velocity: [0.3, 0, 0], seed: 7 });
  clk.advance(60);
  const s = clk.state();
  assert.ok(Math.abs(s.lapse - Math.sqrt(0.82)) < 1e-15);
  assert.ok(Math.abs(s.tau - s.lapse * 60) < clk.delta + 1e-9);
  assert.ok(Math.abs(s.phaseTime - s.tau) < 1e-9, `deterministic hand ${s.phaseTime} vs tau ${s.tau}`);
  assert.ok(Math.abs(s.phaseTimeConditional - s.tau) < 1e-9, `conditional hand ${s.phaseTimeConditional} vs tau ${s.tau}`);
  assert.ok(s.clicks > 0, 'expected recovery clicks');
});

test('running clock: linear clock trajectories pick up phase kicks at clicks', () => {
  const clk = new NodeClock({ omega: 2, c: 1, gamma: 4, clock: 'linear', eta: 1, seed: 11 });
  clk.advance(200);
  const s = clk.state();
  assert.ok(s.clicks > 10);
  assert.ok(Math.abs(s.phaseTimeConditional - s.tau) > 1e-4, 'linear-clock trajectory should not read τ exactly');
});

test('eta = 0: the conditional state is the deterministic state', () => {
  const clk = new NodeClock({ omega: 2, c: 1, gamma: 20, clock: 'linear', eta: 0 });
  clk.advance(25);
  const s = clk.state();
  assert.equal(s.clicks, 0);
  assert.ok(Math.abs(s.conditional.retention - s.deterministic.retention) < 1e-13);
  assert.ok(Math.hypot(s.conditional.coherence[0] - s.deterministic.coherence[0], s.conditional.coherence[1] - s.deterministic.coherence[1]) < 1e-13);
});

test('network trajectories: three nodes, clicks are recorded per node and states stay normalised', () => {
  const net = new Network({ n: 3, edges: [[0, 1], [1, 2]], omega: 2, c: 1, gamma: 10, J: 0.2, clock: 'matched', lapses: [1, 0.9, 0.8] });
  const psi0 = Network.productState([Network.clockVector(), Network.clockVector(), Network.clockVector()]);
  const { psi, clicks } = net.trajectory(psi0, 20, 0.01, xorshift(3));
  let nn = 0; for (let i = 0; i < psi.re.length; i++) nn += psi.re[i] ** 2 + psi.im[i] ** 2;
  assert.ok(Math.abs(nn - 1) < 1e-12);
  assert.ok(clicks.length > 0 && clicks.every(c => c.node >= 0 && c.node < 3));
});

test('sweep: numeric and closed-form columns agree', () => {
  const res = sweep({ param: 'gamma', from: 5, to: 60, steps: 6, base: { omega: 2, c: 1, clock: 'linear' } });
  for (const r of res.rows) {
    assert.ok(Math.abs(r.retention_a - r.retention_a_closed) < 1e-13);
    assert.ok(Math.abs(r.retention_b - r.retention_b_closed) < 1e-13);
    assert.ok(r.dephasing > 0);
  }
  const m = sweep({ param: 'c', from: 0.1, to: 3, steps: 5, base: { omega: 2, gamma: 20, clock: 'matched' } });
  for (const r of m.rows) assert.ok(Math.abs(r.dephasing) < 1e-12);
});

test('permalink round trip', () => {
  const p = { omega: 0.785, c: 1, gamma: 20, eta: 0.5, clock: 'matched', velocity: [0.3, 0, 0.1], metric: 'e47', theta: 1.2, phase: 0.4, seed: 9 };
  assert.deepEqual(decodeHash(encodeHash(p)), p);
});

test('radix dialects: exact repetends', () => {
  assert.equal(radixString(1, 4, 12), '0.3');
  assert.equal(radixString(1, 3, 60), '0;20');
  assert.equal(radixString(47, 125, 5), '0.142');
  assert.equal(radixString(1, 7, 10), '0.(142857)');
});

test('Matrix bridge: tick receipt for 1..8 ticks uses the shortest T word', async () => {
  assert.deepEqual([1, 2, 3, 4, 5, 6, 7, 8].map(k => Rc.tickWord(k).join('')), ['T', 'S', 'ST', 'Z', 'ZT', 'ZS', 'ZST', '']);
  const circuit = [{ gate: 'H', targets: [0], column: 0 }];
  const src = await Rc.seal({ schema: Rc.MATRIX_SCHEMA, source: 'test', exported_at: 'x', qubits: 3, basis: 'q0q1q2', initial_state: '000', circuit, statevector: Rc.matrixReplay(circuit), measurement: { collapsed: false, outcome: null }, norm: 1 });
  const { receipt } = await Rc.matrixReceiptWithTicks(src, 0, 2);
  const v = await Rc.verifyMatrixReceipt(receipt);
  assert.equal(v.status, 'REPLAY-VERIFIED');
  await assert.rejects(Rc.matrixReceiptWithTicks({ ...src, circuit: [...circuit, { gate: 'X', targets: [1], column: 6 }] }, 0, 7), /8 columns/);
});


test('MODSYNC link preserves private validator and provides a labeled public record', () => {
  const html = readFileSync(resolve(root, 'index.html'), 'utf8');
  const line = html.split('\n').find((s) => s.includes('MC-E47-MODSYNC-20261015-001'));
  assert.ok(line, 'MODSYNC certificate reference is present');
  assert.ok(line.includes('https://github.com/nicholaskouns-create/Time/blob/main/README.md'), 'public README record remains navigable');
  assert.ok(line.includes('https://drive.google.com/file/d/1-bCyEuWcN1rtRb1IJF02O8tTe7RbmaEC/view'), 'original private validator is preserved');
  assert.ok(line.includes('validator on Drive (restricted)'), 'restricted permission is disclosed');
  assert.ok(line.includes('not a substitute for the private validator'), 'public record is not presented as validator');
});
