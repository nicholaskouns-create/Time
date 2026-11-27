// A Dial of Forty-Seven Hours: the proof's finite claims, recomputed on the Quantum Clock's own engine.
// Evidence classes as in engine/certificate.js: E0 exact proof, E1 executable reconstruction.
// Runs in the browser (dial/index.html) and in Node (tests/dial.test.js); no dependencies.
import * as L from '../engine/linalg.js';
import { Carrier, SPECTRUM, siteM } from '../engine/carrier.js';
import { primeClock, primeSpectrumResidual } from '../engine/clocks.js';

export const DIAL_CERTIFICATE_ID = 'MC-E47-DIAL-20261126-001';
const N = 125;
const comm = (A, B) => L.radd(L.rmul(A, B), L.rmul(B, A), -1);
const acomm = (A, B) => L.radd(L.rmul(A, B), L.rmul(B, A));
const totalM = n => siteM(n).reduce((s, x) => s + x, 0);
const outer = (u, v) => { const O = L.rmat(N); for (let i = 0; i < N; i++) if (u[i]) for (let j = 0; j < N; j++) O.a[i * N + j] = u[i] * v[j]; return O; };
const col = (B, c) => Float64Array.from({ length: B.n }, (_, r) => B.a[r * B.m + c]);

// Single-site spin-2 operators on V2⊗V2⊗V2, in the site convention of engine/carrier.js.
function siteOperators() {
  const m = [2, 1, 0, -1, -2], Jp = L.rmat(5), Jz = L.rmat(5);
  for (let i = 0; i < 4; i++) Jp.a[i * 5 + i + 1] = Math.sqrt(6 - m[i + 1] * (m[i + 1] + 1));
  for (let i = 0; i < 5; i++) Jz.a[i * 5 + i] = m[i];
  const embed = (A, slot) => {
    const out = L.rmat(N);
    for (let r = 0; r < N; r++) {
      const d = [Math.floor(r / 25), Math.floor(r / 5) % 5, r % 5];
      for (let k = 0; k < 5; k++) { const v = A.a[d[slot] * 5 + k]; if (!v) continue; const e = d.slice(); e[slot] = k; out.a[r * N + e[0] * 25 + e[1] * 5 + e[2]] += v; }
    }
    return out;
  };
  return [0, 1, 2].map(s => ({ z: embed(Jz, s), p: embed(Jp, s), m: embed(L.rtrans(Jp), s) }));
}

// Total pair exchange X = J1·J2 + J2·J3 + J3·J1, with Jx Jx' + Jy Jy' = (J+ J-' + J- J+')/2 across sites.
export function exchangeOperator() {
  const s = siteOperators(); let X = L.rmat(N);
  for (const [k, l] of [[0, 1], [1, 2], [2, 0]]) {
    X = L.radd(X, L.rmul(s[k].z, s[l].z));
    X = L.radd(X, L.radd(L.rmul(s[k].p, s[l].m), L.rmul(s[k].m, s[l].p)), 0.5);
  }
  return X;
}

// Multiplets of the sector C = lam with spin J: one orthonormal ladder (2J+1 vectors) per copy.
export function multiplets(car, lam, J) {
  const B = car.basis(lam), tops = [];
  for (let c = 0; c < B.m; c++) {
    const v = col(B, c), r0 = v.findIndex(x => Math.abs(x) > 1e-9);
    if (totalM(r0) === J) tops.push(v);
  }
  return tops.map(v0 => {
    const ladder = [v0]; let v = v0;
    for (let k = 0; k < 2 * J; k++) { v = L.rmatvec(car.Jm, v); const nv = Math.hypot(...v); v = v.map(x => x / nv); ladder.push(v); }
    return ladder;
  });
}

// A pairing D across the grading with rank r: u_i <-> w_i for i < r (u in H+, w in H-).
export function pairing(car, r) {
  const Bp = car.basis(6), Bm = car.basis(30); let D = L.rmat(N);
  for (let i = 0; i < r; i++) { const u = col(Bp, i), w = col(Bm, i); D = L.radd(L.radd(D, outer(w, u)), outer(u, w)); }
  return D;
}

export function runDialCertificate(car = new Carrier()) {
  const checks = [];
  const add = (id, title, evidence, { value, target, residual = null, tolerance = null, pass }) =>
    checks.push({ id, title, evidence, value, target, residual, tolerance, pass: !!pass });
  const I = L.reye(N), C = car.C, P = car.P, X = exchangeOperator();
  const J = [car.Jz, car.Jp, car.Jm];
  const commJ = A => Math.max(...J.map(G => L.rmaxabs(comm(A, G))));

  const m = car.multiplicities();
  add('D1', 'Sectors of V₂⊗³: 1, 9, 25, 28, 27, 22, 13 states; the kernel takes 25 + 22 = 47', 'E0+E1',
    { value: m, target: SPECTRUM, pass: Object.keys(SPECTRUM).every(k => m[k] === SPECTRUM[k]) && m[6] + m[30] === 47 });

  const rX = L.rmaxabs(L.radd(L.radd(C, I, -18), X, -2));
  add('D2', 'C = 18I + 2X, with X = J₁·J₂ + J₂·J₃ + J₃·J₁ built site by site', 'E0+E1',
    { value: rX, target: 0, residual: rX, tolerance: 1e-11, pass: rX < 1e-11 });

  const rK = L.rmaxabs(L.radd(car.K, L.radd(L.rmul(X, X), I, -36), -4));
  add('D3', 'K = (C − 6I)(C − 30I) = 4(X² − 36I)', 'E0+E1', { value: rK, target: 0, residual: rK, tolerance: 1e-9, pass: rK < 1e-9 });

  const S = L.rscale(L.rmul(L.radd(L.rscale(I, 18), C, -1), P), 1 / 12);
  const rS = Math.max(L.rmaxabs(L.radd(S, L.rmul(X, P), 1 / 6)), L.rmaxabs(L.radd(L.rmul(S, S), P, -1)), L.rmaxabs(L.radd(S, L.rtrans(S), -1)));
  const trS = L.rtrace(S);
  add('D4', 'Grading S = (18I − C)P/12 = −XP/6: S† = S, S² = P, Tr S = 3', 'E0+E1',
    { value: { trace: trS }, target: { trace: 3 }, residual: Math.max(rS, Math.abs(trS - 3)), tolerance: 1e-10, pass: rS < 1e-10 && Math.abs(trS - 3) < 1e-10 });

  const Pp = L.rscale(L.radd(P, S), 0.5), Pm = L.rscale(L.radd(P, S, -1), 0.5);
  const rP = Math.max(L.rmaxabs(L.radd(Pp, L.rscale(L.rmul(L.radd(L.rscale(I, 30), C, -1), P), 1 / 24), -1)),
    L.rmaxabs(L.radd(Pm, L.rscale(L.rmul(L.radd(C, I, -6), P), 1 / 24), -1)), L.rmaxabs(L.rmul(Pp, Pm)));
  add('D5', 'P₊ = (P + S)/2 = (30I − C)P/24, P₋ = (P − S)/2 = (C − 6I)P/24, P₊P₋ = 0; ranks 25 and 22', 'E0+E1',
    { value: [L.rtrace(Pp), L.rtrace(Pm)], target: [25, 22], residual: rP, tolerance: 1e-10, pass: rP < 1e-10 && Math.abs(L.rtrace(Pp) - 25) < 1e-9 && Math.abs(L.rtrace(Pm) - 22) < 1e-9 });

  // restriction to E47 in the basis [H+ | H-]
  const Bp = car.basis(6), Bm = car.basis(30), B = L.rmat(N, 47);
  for (let r = 0; r < N; r++) { for (let c = 0; c < 25; c++) B.a[r * 47 + c] = Bp.a[r * 25 + c]; for (let c = 0; c < 22; c++) B.a[r * 47 + 25 + c] = Bm.a[r * 22 + c]; }
  const restrict = A => L.rmul(L.rtrans(B), L.rmul(A, B));
  const S47 = restrict(S);
  const pairStats = r => {
    const D = pairing(car, r), { values, vectors } = L.rsymEig(restrict(D));
    const cnt = [-1, 0, 1].map(t => values.filter(v => Math.abs(v - t) < 1e-9).length);
    let chir = 0; values.forEach((v, c) => { if (Math.abs(v) < 1e-9) { const x = col(vectors, c), y = L.rmatvec(S47, x); chir += x.reduce((s, xi, i) => s + xi * y[i], 0); } });
    return { D, cnt, chir, anti: L.rmaxabs(acomm(S, D)) };
  };
  const full = pairStats(22);
  add('D6', 'A full pairing D (AA† = I₂₂) anticommutes with S; spectrum on E₄₇ is −1, 0, +1 with multiplicities 22, 3, 22', 'E0+E1',
    { value: full.cnt, target: [22, 3, 22], residual: full.anti, tolerance: 1e-10, pass: full.anti < 1e-10 && full.cnt.join() === '22,3,22' });

  const part = pairStats(15);
  add('D7', 'The three that remain: zero modes carry Tr S = 3 for every rank (rank 22: 3 zero modes; rank 15: 17 = 3 + 2·7)', 'E0+E1',
    { value: { rank22: [full.cnt[1], full.chir], rank15: [part.cnt[1], part.chir] }, target: { rank22: [3, 3], rank15: [17, 3] },
      residual: Math.max(Math.abs(full.chir - 3), Math.abs(part.chir - 3)), tolerance: 1e-9,
      pass: full.cnt[1] === 3 && part.cnt[1] === 17 && Math.abs(full.chir - 3) < 1e-9 && Math.abs(part.chir - 3) < 1e-9 });

  const copies = [...multiplets(car, 6, 2), ...multiplets(car, 30, 5)];
  const Pi = copies.map(ladder => ladder.reduce((A, v) => L.radd(A, outer(v, v)), L.rmat(N)));
  const sumPi = Pi.reduce((A, Q) => L.radd(A, Q), L.rmat(N));
  const invRes = Math.max(...Pi.map(commJ)), evenRes = Math.max(...Pi.map(Q => L.rmaxabs(comm(Q, S))));
  const brk = commJ(full.D);
  add('D8', 'Symmetry forbids pairing: SU(2)-invariant operators on E₄₇ are even ([·, S] = 0), and the pairing D fails to commute with SU(2)', 'E0+E1',
    { value: { invariantEven: evenRes, pairingCommutator: brk }, target: { invariantEven: 0, pairingCommutator: '> 0' }, residual: evenRes, tolerance: 1e-10,
      pass: evenRes < 1e-10 && brk > 1e-3 });

  const ranks = Pi.map(Q => Math.round(L.rtrace(Q)));
  const rSum = L.rmaxabs(L.radd(sumPi, P, -1));
  add('D9', 'The symmetric dial attains 7 hours: 5 copy projectors of rank 5 and 2 of rank 11 commute with SU(2) and sum to P', 'E0+E1',
    { value: ranks, target: [5, 5, 5, 5, 5, 11, 11], residual: Math.max(invRes, rSum), tolerance: 1e-10,
      pass: ranks.join() === '5,5,5,5,5,11,11' && invRes < 1e-10 && rSum < 1e-10 });

  const f = car.functionOfC(lam => lam / 24);              // the linear clock f = ωC/24, ω = 1
  const rf = Math.max(L.rmaxabs(L.radd(L.rmul(f, P), L.radd(L.rscale(P, 0.75), S, -0.5), -1)), commJ(f));
  add('D10', 'The Quantum Clock\'s hand is the grading: f(C)P = ω(¾P − ½S), [f(C), SU(2)] = 0; one 1/47 turn is U₀ up to a global phase', 'E0+E1',
    { value: { hours: 2, relativePhasePerTurnOver47: (30 - 6) / 24 }, target: { hours: 2, relativePhasePerTurnOver47: 1 }, residual: rf, tolerance: 1e-10, pass: rf < 1e-10 });

  const pc = primeClock(car), res47 = primeSpectrumResidual(pc);
  const keep = L.rmaxabs(comm(pc.U, P)), brk47 = commJ(pc.U), crossS = L.rmaxabs(comm(pc.U, S));
  add('D11', 'The prime 47-cycle U₄₇ keeps E₄₇ and shows all 47 hours, so it cannot commute with SU(2)', 'E0+E1',
    { value: { keepsKernel: keep, spectrumResidual: res47, su2Commutator: brk47, gradingCommutator: crossS }, target: { keepsKernel: 0, spectrumResidual: 0, su2Commutator: '> 0' },
      residual: Math.max(keep, res47), tolerance: 1e-9, pass: keep < 1e-9 && res47 < 1e-9 && brk47 > 1e-3 });

  const F = [0, 1], Lu = [2, 1]; for (let i = 2; i < 20; i++) { F.push(F[i - 1] + F[i - 2]); Lu.push(Lu[i - 1] + Lu[i - 2]); }
  const tower = [3, 5, 7, 9, 11, 13].map(k => ({ k, spinTimesTwo: F[k] - 1, plus: F[k] ** 2, minus: F[k - 2] * Lu[k], index: F[k] ** 2 - F[k - 2] * Lu[k] }));
  add('D12', 'Index along the Fibonacci–Lucas tower: F_k² − F_{k−2}L_k = F_{k−3}F_{k−1} = 0, 3, 24, 168, 1155, 7920', 'E0',
    { value: tower.map(t => t.index), target: [0, 3, 24, 168, 1155, 7920], pass: tower.every(t => t.index === F[t.k - 3] * F[t.k - 1]) && tower.map(t => t.index).join() === '0,3,24,168,1155,7920' });

  return { id: DIAL_CERTIFICATE_ID, checks, passed: checks.filter(c => c.pass).length, total: checks.length, tower };
}
