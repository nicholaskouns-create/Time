// The E47 clock-memory node (the Quantum Clock's unit cell).
//
// Two clock states inside E47: |a> = pair singlet ⊗ |m3 = 0> (C = 6) and |b> = P30|38>/‖·‖ (C = 30).
// Two syndrome states outside E47: |a'> = P_s1|38>/‖·‖, |b'> = P_s2|38>/‖·‖ (default s = (12, 42)).
// V = [a, b], W = [a', b'].  Generator on the 125-dim carrier, in proper time τ:
//     H = f(C) + c (W V† + V W†),     R = √γ V W†       (leak and logical-blind recovery)
// f(C) is the clock: linear f = ωC/24 ticks at f(30) − f(6) = ω inside E47.
// Every operator maps span{a, a', b, b'} into itself, so the node is exactly 4-dimensional:
// index 2l + f with l = 0 (a), 1 (b) and f = 0 (protected), 1 (leaked).
import { cmat, ceye, cmul, cadj, ckron, caxpy, cscale, csub, cexpm, csolve, ctrace, cfromRows, smallEigvals, rmat, rmul, rtrans, radd, rmaxabs } from './linalg.js';
import { Carrier, SPECTRUM, REF_SITE, pairSingletCopy } from './carrier.js';

export const NONSELECTED = [0, 2, 12, 20, 42];

export function clockFunction(omega, clock = 'linear', syndrome = [12, 42]) {
  const lin = {}; for (const lam of Object.keys(SPECTRUM).map(Number)) lin[lam] = omega * lam / 24;
  if (clock === 'linear') return lin;
  if (clock === 'matched') { const f = { ...lin }; f[syndrome[1]] = lin[syndrome[0]] + (lin[30] - lin[6]); return f; }
  throw Error('clock must be linear or matched');
}

export function detunings(omega, clock = 'linear', syndrome = [12, 42]) {
  const f = clockFunction(omega, clock, syndrome);
  return [f[syndrome[0]] - f[6], f[syndrome[1]] - f[30]];
}

// 4x4 reduced (H, R) on logical ⊗ flag, index 2l + f
export function reduced({ omega, c, gamma, clock = 'linear', syndrome = [12, 42] }) {
  const f = clockFunction(omega, clock, syndrome), sec = [6, syndrome[0], 30, syndrome[1]];
  const H = cmat(4), R = cmat(4);
  for (let k = 0; k < 4; k++) H.re[k * 4 + k] = f[sec[k]];
  for (const l of [0, 1]) { H.re[(2 * l) * 4 + 2 * l + 1] += c; H.re[(2 * l + 1) * 4 + 2 * l] += c; R.re[(2 * l) * 4 + 2 * l + 1] = Math.sqrt(gamma); }
  return { H, R };
}

// 2x2 branch (H, R) on the flag of logical branch l
export function branch(l, { omega, c, gamma, clock = 'linear', syndrome = [12, 42] }) {
  const f = clockFunction(omega, clock, syndrome), sec = l === 0 ? [6, syndrome[0]] : [30, syndrome[1]];
  return { H: cfromRows([[f[sec[0]], c], [c, f[sec[1]]]]), R: cfromRows([[0, Math.sqrt(gamma)], [0, 0]]) };
}

// ---------------------------------------------------------------------------
// Liouville space (row-major vectorisation: vec(A X B) = (A ⊗ Bᵀ) vec(X))

const transpose = A => { const B = cmat(A.m, A.n); for (let i = 0; i < A.n; i++) for (let j = 0; j < A.m; j++) { B.re[j * A.n + i] = A.re[i * A.m + j]; B.im[j * A.n + i] = A.im[i * A.m + j]; } return B; };
const conj = A => ({ n: A.n, m: A.m, re: A.re.slice(), im: A.im.map(x => -x) });

export function liouvillian(H, Ls) {
  const d = H.n, I = ceye(d);
  const Lv = cscale(csub(ckron(H, I), ckron(I, transpose(H))), 0, -1);
  for (const L of Ls) {
    const LdL = cmul(cadj(L), L);
    caxpy(1, 0, ckron(L, conj(L)), Lv);
    caxpy(-0.5, 0, ckron(LdL, I), Lv);
    caxpy(-0.5, 0, ckron(I, transpose(LdL)), Lv);
  }
  return Lv;
}
export const jumpSuper = L => ckron(L, conj(L));

export function vecOf(rho) { const v = cmat(rho.n * rho.m, 1); v.re.set(rho.re); v.im.set(rho.im); return v; }
export function unvec(v, d) { const r = cmat(d); r.re.set(v.re); r.im.set(v.im); return r; }
export function applySuper(S, rho) { return unvec(cmul(S, vecOf(rho)), rho.n); }

export function steadyState(Lv, d) {
  // (L + u tᵀ) ρ = u with u = vec(I/d), t = vec(I): unique stationary state when it exists
  const M = { n: Lv.n, m: Lv.m, re: Lv.re.slice(), im: Lv.im.slice() }, u = cmat(d * d, 1);
  for (let i = 0; i < d; i++) u.re[i * d + i] = 1 / d;
  for (let r = 0; r < d * d; r++) for (let i = 0; i < d; i++) M.re[r * d * d + i * d + i] += u.re[r];
  const x = csolve(M, u), rho = unvec(x, d), [tr] = ctrace(rho);
  for (let k = 0; k < rho.re.length; k++) { rho.re[k] /= tr; rho.im[k] /= tr; }
  return rho;
}

// Long-time counting statistics of the jumps of channel L (detected with efficiency η):
// K1 = mean rate, K2 = variance rate (full counting statistics to second order), F = K2/K1.
export function countingStats(H, Ls, counted, eta = 1) {
  const d = H.n, Lv = liouvillian(H, Ls), rho = steadyState(Lv, d), J = jumpSuper(Ls[counted]);
  const t = v => { let s = 0; for (let i = 0; i < d; i++) s += v.re[i * d + i]; return s; };
  const rv = vecOf(rho), Jr = cmul(J, rv), K1 = eta * t(Jr);
  const M = { n: Lv.n, m: Lv.m, re: Lv.re.slice(), im: Lv.im.slice() };
  for (let r = 0; r < d * d; r++) for (let i = 0; i < d; i++) { M.re[r * d * d + i * d + i] += rv.re[r]; M.im[r * d * d + i * d + i] += rv.im[r]; }
  const b = cmat(d * d, 1); for (let k = 0; k < d * d; k++) { b.re[k] = -(eta * Jr.re[k] - K1 * rv.re[k]); b.im[k] = -(eta * Jr.im[k] - K1 * rv.im[k]); }
  const rho1 = csolve(M, b), K2 = K1 + 2 * eta * t(cmul(J, rho1));
  return { rate: K1, variance: K2, fano: K2 / K1, steady: rho };
}

// ---------------------------------------------------------------------------
// closed forms (proper time)

export const stationaryLeak = (c, gamma, delta = 0) => 4 * c * c / (gamma * gamma + 8 * c * c + 4 * delta * delta);
export const resonantFano = (c, gamma) => 1 - 24 * c * c * gamma * gamma / (gamma * gamma + 8 * c * c) ** 2;
export const fanoAtEfficiency = (F, eta) => 1 + eta * (F - 1);
export const clickTimePrecision = (F, rate, tau) => Math.sqrt(F / (rate * tau));

// rate λ of ρ_ab = <a|ρ|b>:  λ = iω + i c²(z_a − z̄_b) + γ c² z_a z̄_b,  z_l = 1/(Δ_l − iγ/2)   (second order in c/γ)
export function adiabaticCoherenceRate(c, gamma, da, db, omega) {
  const z = d => { const den = d * d + gamma * gamma / 4; return [d / den, (gamma / 2) / den]; }; // 1/(d − iγ/2) = (d + iγ/2)/(d² + γ²/4)
  const [ar, ai] = z(da), [br, bi] = z(db);
  // i c²(z_a − z̄_b) = −c²(ai + bi) + i c²(ar − br);   γc² z_a z̄_b = γc²(ar br + ai bi) + i γc²(ai br − ar bi)
  return [-c * c * (ai + bi) + gamma * c * c * (ar * br + ai * bi), omega + c * c * (ar - br) + gamma * c * c * (ai * br - ar * bi)];
}
export function adiabaticDephasing(c, gamma, da, db) {
  const z = d => { const den = d * d + gamma * gamma / 4; return [d / den, (gamma / 2) / den]; };
  const [ar, ai] = z(da), [br, bi] = z(db);
  return 0.5 * gamma * c * c * ((ar - br) ** 2 + (ai - bi) ** 2);
}

// exact rate of the slowest coherence mode: eigenvalues of the ρ_ba block generator
export function coherenceModeRates(params) {
  const A = branch(0, params), B = branch(1, params), I = ceye(2);
  const R = A.R, RdR = cmul(cadj(R), R);
  // d X/dt = −i(H_b X − X H_a) + R X R† − ½{R†R, X}   for X = ρ_ba (2x2 on the flag)
  const Lba = cscale(csub(ckron(B.H, I), ckron(I, transpose(A.H))), 0, -1);
  caxpy(1, 0, ckron(R, conj(R)), Lba); caxpy(-0.5, 0, ckron(RdR, I), Lba); caxpy(-0.5, 0, ckron(I, transpose(RdR)), Lba);
  const ev = smallEigvals(Lba).sort((x, y) => y[0] - x[0]);
  return { eigenvalues: ev, slowest: ev[0], generator: Lba };
}

// ---------------------------------------------------------------------------
// the full 125-dim node, for the exact-reduction certificate

export class ClockMemoryNode {
  constructor(carrier = null, syndrome = [12, 42]) {
    if (syndrome.some(s => s === 6 || s === 30 || !(s in SPECTRUM)) || syndrome[0] === syndrome[1]) throw Error('syndrome sectors must be two distinct Casimir eigenvalues outside E47');
    const car = carrier || new Carrier();
    this.car = car; this.syndrome = syndrome;
    const e = new Float64Array(125); e[REF_SITE] = 1;
    const proj = lam => { const P = car.projector(lam), v = new Float64Array(125); for (let i = 0; i < 125; i++) v[i] = P.a[i * 125 + REF_SITE]; const n = Math.hypot(...v); if (n < 1e-3) throw Error(`site ${REF_SITE} has no C = ${lam} component`); return v.map(x => x / n); };
    this.a = pairSingletCopy(2);
    this.b = proj(30); this.as = proj(syndrome[0]); this.bs = proj(syndrome[1]);
    this.sector = [6, syndrome[0], 30, syndrome[1]];
    this.U = rmat(125, 4);
    [this.a, this.as, this.b, this.bs].forEach((v, k) => { for (let r = 0; r < 125; r++) this.U.a[r * 4 + k] = v[r]; });
  }
  full({ omega, c, gamma, clock = 'linear' }) {
    const f = clockFunction(omega, clock, this.syndrome), H = this.car.functionOfC(lam => f[lam]), R = rmat(125);
    const add = (X, x, y, s) => { for (let i = 0; i < 125; i++) if (x[i] !== 0) for (let j = 0; j < 125; j++) X.a[i * 125 + j] += s * x[i] * y[j]; };
    add(H, this.as, this.a, c); add(H, this.a, this.as, c); add(H, this.bs, this.b, c); add(H, this.b, this.bs, c);
    add(R, this.a, this.as, Math.sqrt(gamma)); add(R, this.b, this.bs, Math.sqrt(gamma));
    return { H, R };
  }
  // residuals of the exact reduction: span invariance and U†XU = X_red
  reductionResiduals(params) {
    const { H, R } = this.full(params), Ut = rtrans(this.U), Pi = rmul(this.U, Ut);
    const leak = X => { const XU = rmul(X, this.U); return rmaxabs(radd(XU, rmul(Pi, XU), -1)); };
    const red = reduced({ ...params, syndrome: this.syndrome });
    const diff = (X, Xr) => { const Y = rmul(Ut, rmul(X, this.U)); let m = 0; for (let k = 0; k < 16; k++) m = Math.max(m, Math.abs(Y.a[k] - Xr.re[k]), Math.abs(Xr.im[k])); return m; };
    return {
      orthonormality: rmaxabs(radd(rmul(Ut, this.U), { n: 4, m: 4, a: Float64Array.from([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]) }, -1)),
      H_invariance: leak(H), R_invariance: leak(R), Rdag_invariance: leak(rtrans(R)),
      H_reduced: diff(H, red.H), R_reduced: diff(R, red.R),
    };
  }
}

// ---------------------------------------------------------------------------
// states and readouts (index 2l + f)

export function clockState(theta = Math.PI / 2, phase = 0, flag = 0) {
  // cos(θ/2)|a> + e^{iφ} sin(θ/2)|b>, as a density matrix
  const psi = [[0, 0], [0, 0], [0, 0], [0, 0]];
  psi[flag] = [Math.cos(theta / 2), 0];
  psi[2 + flag] = [Math.cos(phase) * Math.sin(theta / 2), Math.sin(phase) * Math.sin(theta / 2)];
  const rho = cmat(4);
  for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) {
    rho.re[i * 4 + j] = psi[i][0] * psi[j][0] + psi[i][1] * psi[j][1];
    rho.im[i * 4 + j] = psi[i][1] * psi[j][0] - psi[i][0] * psi[j][1];
  }
  return rho;
}
export function readout(rho) {
  const g = (i, j) => [rho.re[i * 4 + j], rho.im[i * 4 + j]];
  const retention = rho.re[0] + rho.re[10];
  const [c1r, c1i] = g(2, 0), [c2r, c2i] = g(3, 1);       // Tr(|a><b| ⊗ I ρ) = ρ_{b0,a0} + ρ_{b1,a1}
  return { retention, leaked: rho.re[5] + rho.re[15], coherence: [c1r + c2r, c1i + c2i], pa: rho.re[0] + rho.re[5], pb: rho.re[10] + rho.re[15] };
}

export { Carrier };
