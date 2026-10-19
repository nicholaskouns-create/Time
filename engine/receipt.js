// Receipts: the Matrix Quantum Simulator's state receipt (city.quantum.state-receipt.v1),
// verified by digest AND by replaying its circuit with the Matrix's own gate arithmetic,
// and the Quantum Clock's own receipt (city.quantum.clock-receipt.v1).
import { cmat, hermEig } from './linalg.js';

export const MATRIX_SCHEMA = 'city.quantum.state-receipt.v1';
export const CLOCK_SCHEMA = 'city.quantum.clock-receipt.v1';
export const GATES = ['H', 'X', 'Y', 'Z', 'S', 'T', 'CX', 'SW'];

// sorted-json-v1, identical to matrix-quantum-simulator/receipt.js (number formatting is ECMAScript's, as in RFC 8785)
export function canonical(x) {
  if (Array.isArray(x)) return '[' + x.map(canonical).join(',') + ']';
  if (x && typeof x === 'object') return '{' + Object.keys(x).sort().map(k => JSON.stringify(k) + ':' + canonical(x[k])).join(',') + '}';
  return JSON.stringify(x);
}
export async function sha256hex(text) {
  const d = await globalThis.crypto.subtle.digest('SHA-256', new TextEncoder().encode(text));
  return [...new Uint8Array(d)].map(v => v.toString(16).padStart(2, '0')).join('');
}
export const digestOf = payload => sha256hex(canonical(payload));
export async function seal(payload) { return { ...payload, integrity: { algorithm: 'SHA-256', canonicalization: 'sorted-json-v1', digest: await digestOf(payload) } }; }

// ---------------------------------------------------------------------------
// The Matrix engine, transcribed operation for operation from matrix-quantum-simulator/app.html,
// so a consistent receipt replays bit for bit.
const C = (r = 0, i = 0) => ({ r, i }), A = (a, b) => C(a.r + b.r, a.i + b.i), M = (a, b) => C(a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r), P = a => a.r * a.r + a.i * a.i, S = Math.SQRT1_2;
const U = { H: [[C(S), C(S)], [C(S), C(-S)]], X: [[C(), C(1)], [C(1), C()]], Y: [[C(), C(0, -1)], [C(0, 1), C()]], Z: [[C(1), C()], [C(), C(-1)]], S: [[C(1), C()], [C(), C(0, 1)]], T: [[C(1), C()], [C(), C(S, S)]] };
const init = () => Array.from({ length: 8 }, (_, i) => C(i ? 0 : 1));
const bit = q => 1 << (2 - q);
function one(s, g, q) { const o = s.map(z => C(z.r, z.i)), m = bit(q); for (let i = 0; i < 8; i++) if (!(i & m)) { const j = i | m, a = s[i], b = s[j]; o[i] = A(M(U[g][0][0], a), M(U[g][0][1], b)); o[j] = A(M(U[g][1][0], a), M(U[g][1][1], b)); } return o; }
function two(s, g, a, b) { const o = s.map(z => C(z.r, z.i)), ma = bit(a), mb = bit(b); if (g === 'CX') { for (let i = 0; i < 8; i++) if ((i & ma) && !(i & mb)) { const j = i | mb; [o[i], o[j]] = [o[j], o[i]]; } } else { for (let i = 0; i < 8; i++) if (!(i & ma) && (i & mb)) { const j = (i | ma) & ~mb; [o[i], o[j]] = [o[j], o[i]]; } } return o; }
export function matrixApply(psi, gate, targets) { return targets.length === 2 ? two(psi, gate, targets[0], targets[1]) : one(psi, gate, targets[0]); }
export function matrixReplay(circuit, order = 'column', start = null) {
  const ops = order === 'column' ? circuit.map((o, k) => [o, k]).sort((x, y) => x[0].column - y[0].column || x[1] - y[1]).map(x => x[0]) : circuit;
  let psi = start ? start.map(([r, i]) => C(r, i)) : init();
  for (const o of ops) psi = matrixApply(psi, o.gate, o.targets);
  return psi.map(z => [z.r, z.i]);
}
const l2 = (x, y) => Math.sqrt(x.reduce((s, z, k) => s + (z[0] - y[k][0]) ** 2 + (z[1] - y[k][1]) ** 2, 0));

// The importer of matrix-quantum-simulator/receipt.js as deployed, transcribed verbatim. Note that it
// validates each gate as o.g while the exporter writes o.gate, so it rejects every receipt the Matrix
// exports with a non-empty circuit. Receipts written here carry both keys so the deployed importer accepts them.
export function matrixAppValid(x) {
  if (x.schema !== MATRIX_SCHEMA || x.qubits !== 3 || x.basis !== 'q0q1q2' || !Array.isArray(x.statevector) || x.statevector.length !== 8 || !Array.isArray(x.circuit)) throw Error('Unsupported receipt');
  if (!x.statevector.every(z => Array.isArray(z) && z.length === 2 && z.every(Number.isFinite))) throw Error('Invalid amplitudes');
  const norm = x.statevector.reduce((s, z) => s + z[0] ** 2 + z[1] ** 2, 0); if (Math.abs(norm - 1) > 1e-8) throw Error('State is not normalized');
  if (x.circuit.length > 1000 || !x.circuit.every(o => GATES.includes(o.g) && Array.isArray(o.targets) && o.targets.length === (['CX', 'SW'].includes(o.g) ? 2 : 1) && o.targets.every(q => Number.isInteger(q) && q >= 0 && q < 3) && new Set(o.targets).size === o.targets.length && Number.isInteger(o.column) && o.column >= 0 && o.column < 8)) throw Error('Invalid circuit');
  return true;
}

export function checkMatrixShape(x) { // the Matrix's validation, reading the gate from o.gate as its exporter writes it
  if (!x || x.schema !== MATRIX_SCHEMA || x.qubits !== 3 || x.basis !== 'q0q1q2' || !Array.isArray(x.statevector) || x.statevector.length !== 8 || !Array.isArray(x.circuit)) throw Error('Unsupported receipt (expects city.quantum.state-receipt.v1, 3 qubits, basis q0q1q2)');
  if (!x.statevector.every(z => Array.isArray(z) && z.length === 2 && z.every(Number.isFinite))) throw Error('Invalid amplitudes');
  const norm = x.statevector.reduce((s, z) => s + z[0] ** 2 + z[1] ** 2, 0);
  if (Math.abs(norm - 1) > 1e-8) throw Error('State is not normalized');
  if (x.circuit.length > 1000 || !x.circuit.every(o => GATES.includes(o.gate) && Array.isArray(o.targets) && o.targets.length === (['CX', 'SW'].includes(o.gate) ? 2 : 1) && o.targets.every(q => Number.isInteger(q) && q >= 0 && q < 3) && new Set(o.targets).size === o.targets.length && Number.isInteger(o.column) && o.column >= 0 && o.column < 8)) throw Error('Invalid circuit');
  return norm;
}

// Full verification: digest, then replay in diagram (column) order and in listed order.
export async function verifyMatrixReceipt(x, tol = 1e-12) {
  const norm = checkMatrixShape(x);
  if (!x.integrity || x.integrity.algorithm !== 'SHA-256') throw Error('Missing SHA-256 digest');
  const { integrity, ...payload } = x, digest = await digestOf(payload), digestOk = digest === integrity.digest;
  const occupancy = new Map(); let collisions = 0;
  for (const o of x.circuit) for (const q of o.targets) { const k = o.column + ':' + q; occupancy.set(k, (occupancy.get(k) || 0) + 1); if (occupancy.get(k) === 2) collisions++; }
  const col = matrixReplay(x.circuit, 'column'), ins = matrixReplay(x.circuit, 'insertion');
  const rCol = l2(x.statevector, col), rIns = l2(x.statevector, ins), orderGap = l2(col, ins);
  let status, detail;
  if (!digestOk) { status = 'DIGEST-MISMATCH'; detail = 'The SHA-256 digest does not match the canonical payload.'; }
  else if (x.measurement?.collapsed) {
    const k = x.measurement.outcome, basis = Number.isInteger(k) && Math.abs(x.statevector[k]?.[0] ** 2 + x.statevector[k]?.[1] ** 2 - 1) < 1e-9;
    const pCol = basis ? col[k][0] ** 2 + col[k][1] ** 2 : 0, pIns = basis ? ins[k][0] ** 2 + ins[k][1] ** 2 : 0;
    status = basis && pCol > 1e-12 ? 'MEASURED-CONSISTENT' : basis && pIns > 1e-12 ? 'MEASURED-ORDER-DEPENDENT' : 'MEASURED-INCONSISTENT';
    detail = `Collapsed to ${Number.isInteger(k) && k >= 0 && k < 8 ? '|' + k.toString(2).padStart(3, '0') + '⟩' : 'an invalid outcome'}; probability of that outcome under the circuit: ${pCol.toExponential(3)} (diagram order), ${pIns.toExponential(3)} (listed order).`;
  } else if (rCol <= tol) { status = 'REPLAY-VERIFIED'; detail = 'The circuit, run in diagram order from |000⟩, reproduces the statevector.'; }
  else if (rIns <= tol) { status = 'ORDER-DEPENDENT'; detail = 'The statevector matches the gates in the order they were placed, not the order the diagram shows.'; }
  else { status = 'REPLAY-MISMATCH'; detail = 'Neither order of the circuit reproduces the statevector (state not re-run after editing, or amplitudes changed).'; }
  return { status, detail, digest, digestOk, norm, residualColumnOrder: rCol, residualListedOrder: rIns, orderGap, columnCollisions: collisions, gates: x.circuit.length, maxColumn: Math.max(-1, ...x.circuit.map(o => o.column)) };
}

// ---------------------------------------------------------------------------
// Bridge: one Matrix qubit becomes the clock's logical qubit, the other two ride along as spectators.
// Encoding: |0> -> |b> (C = 30), |1> -> |a> (C = 6). With it, one forward tick ωτ = π/4 is the Matrix T gate.
export const ENCODING = { 0: 'b (C = 30)', 1: 'a (C = 6)', note: 'one forward tick of ωτ = π/4 equals the Matrix T gate up to a global phase' };
const spectators = k => [0, 1, 2].filter(q => q !== k);
// 3-qubit amplitude index i (q0 most significant) <-> joint index (2l + f)*4 + s
export function encodeQubit(statevector, k, flag = 0) {
  const [s0, s1] = spectators(k), out = new Array(16).fill(0).map(() => [0, 0]);
  for (let i = 0; i < 8; i++) {
    const x = q => (i >> (2 - q)) & 1, l = 1 - x(k), s = x(s0) * 2 + x(s1);
    out[(2 * l + flag) * 4 + s] = statevector[i].slice();
  }
  return out;
}
export function densityOf(vec) { const d = vec.length, R = cmat(d); for (let i = 0; i < d; i++) for (let j = 0; j < d; j++) { R.re[i * d + j] = vec[i][0] * vec[j][0] + vec[i][1] * vec[j][1]; R.im[i * d + j] = vec[i][1] * vec[j][0] - vec[i][0] * vec[j][1]; } return R; }
// logical 3-qubit density (8x8, Matrix basis) = partial trace of the flag
export function logicalDensity(rho16, k) {
  const [s0, s1] = spectators(k), L = cmat(8);
  const jointIndex = (i, f) => { const x = q => (i >> (2 - q)) & 1; return (2 * (1 - x(k)) + f) * 4 + x(s0) * 2 + x(s1); };
  for (let i = 0; i < 8; i++) for (let j = 0; j < 8; j++) for (const f of [0, 1]) {
    const a = jointIndex(i, f), b = jointIndex(j, f);
    L.re[i * 8 + j] += rho16.re[a * 16 + b]; L.im[i * 8 + j] += rho16.im[a * 16 + b];
  }
  return L;
}
export function expectationIn(Lrho, vec) { // <ψ|ρ|ψ>
  let r = 0;
  for (let i = 0; i < 8; i++) for (let j = 0; j < 8; j++) {
    const ar = vec[i][0], ai = -vec[i][1], br = vec[j][0], bi = vec[j][1], mr = Lrho.re[i * 8 + j], mi = Lrho.im[i * 8 + j];
    r += (ar * (mr * br - mi * bi) - ai * (mr * bi + mi * br));
  }
  return r;
}
export function purity(R) { let s = 0; for (let i = 0; i < R.n; i++) for (let j = 0; j < R.n; j++) s += R.re[i * R.n + j] * R.re[j * R.n + i] - R.im[i * R.n + j] * R.im[j * R.n + i]; return s; }
export const tickWord = k => { const m = ((k % 8) + 8) % 8, w = []; if (m & 4) w.push('Z'); if (m & 2) w.push('S'); if (m & 1) w.push('T'); return w; };

// A Matrix receipt carrying `ticks` forward clock ticks on qubit k, appended after the existing gates.
export async function matrixReceiptWithTicks(source, k, ticks) {
  const word = tickWord(ticks), start = Math.max(-1, ...source.circuit.map(o => o.column)) + 1;
  if (start + word.length > 8) throw Error(`The Matrix grid has 8 columns; ${word.length} tick gate(s) after column ${start - 1} do not fit.`);
  let psi = source.statevector.map(z => z.slice());
  // each gate is written under both `gate` (what the Matrix exporter and reader use) and `g` (what its validator checks)
  const circuit = source.circuit.map(o => ({ gate: o.gate, g: o.gate, targets: o.targets.slice(), column: o.column }));
  word.forEach((g, n) => { circuit.push({ gate: g, g, targets: [k], column: start + n }); psi = matrixApply(psi.map(([r, i]) => ({ r, i })), g, [k]).map(z => [z.r, z.i]); });
  const payload = { schema: MATRIX_SCHEMA, source: 'quantum-clock', exported_at: new Date().toISOString(), qubits: 3, basis: 'q0q1q2', initial_state: '000', circuit, statevector: psi, measurement: { collapsed: false, outcome: null }, norm: psi.reduce((s, z) => s + z[0] ** 2 + z[1] ** 2, 0) };
  return { receipt: await seal(payload), word };
}

export function checkClockReceiptShape(x) {
  if (!x || x.schema !== CLOCK_SCHEMA) throw Error('Not a city.quantum.clock-receipt.v1');
  const R = x.state?.density, d = x.state?.dim;
  if (!Array.isArray(R) || !Number.isInteger(d) || R.length !== d || !R.every(row => Array.isArray(row) && row.length === d && row.every(z => Array.isArray(z) && z.length === 2 && z.every(Number.isFinite)))) throw Error('Invalid density matrix');
  const M = cmat(d); for (let i = 0; i < d; i++) for (let j = 0; j < d; j++) { M.re[i * d + j] = R[i][j][0]; M.im[i * d + j] = R[i][j][1]; }
  let herm = 0, tr = 0; for (let i = 0; i < d; i++) { tr += M.re[i * d + i]; for (let j = 0; j < d; j++) herm = Math.max(herm, Math.hypot(M.re[i * d + j] - M.re[j * d + i], M.im[i * d + j] + M.im[j * d + i])); }
  const minEig = Math.min(...hermEig(M).values);
  return { hermiticity: herm, trace: tr, minEigenvalue: minEig, ok: herm < 1e-10 && Math.abs(tr - 1) < 1e-10 && minEig > -1e-10, M };
}
export async function verifyClockReceipt(x) {
  const shape = checkClockReceiptShape(x);
  if (!x.integrity || x.integrity.algorithm !== 'SHA-256') throw Error('Missing SHA-256 digest');
  const { integrity, ...payload } = x, digest = await digestOf(payload);
  return { ...shape, digest, digestOk: digest === integrity.digest };
}
