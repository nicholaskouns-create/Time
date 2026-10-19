// Discrete clocks of the corpus: the prime 47-cycle on E47, the 5x5x5 torus clock,
// the 235 joint period, and the E47-modular synchronisation clock z -> z^2 mod 5875 (period 220).
import { rmat, rmul, rtrans, radd, reye, rmaxabs, rsymEig } from './linalg.js';
import { FILTER_DENOM } from './carrier.js';

// ---------------------------------------------------------------- prime clock U47
export function primeClock(car) {
  const B6 = car.basis(6), B30 = car.basis(30), B = rmat(125, 47);
  for (let r = 0; r < 125; r++) {
    for (let c = 0; c < 25; c++) B.a[r * 47 + c] = B6.a[r * 25 + c];
    for (let c = 0; c < 22; c++) B.a[r * 47 + 25 + c] = B30.a[r * 22 + c];
  }
  const S = rmat(47); for (let k = 0; k < 47; k++) S.a[((k + 1) % 47) * 47 + k] = 1;   // S e_k = e_{k+1}
  const Bt = rtrans(B), U = radd(radd(reye(125), rmul(B, Bt), -1), rmul(B, rmul(S, Bt)));
  return { B, S, U };
}

export function orderOf(U, kmax, tol = 1e-9) {
  const n = U.n, I = reye(n); let X = reye(n);
  const dist = [];
  for (let k = 1; k <= kmax; k++) { X = rmul(X, U); const d = rmaxabs(radd(X, I, -1)); dist.push(d); if (d < tol) return { order: k, dist }; }
  return { order: null, dist };
}

// eigenvectors of the cycle: f_r = Σ_k ω^{-rk} e_k has S f_r = ω^r f_r;  residual of U (B f_r) = ω^r (B f_r)
export function primeSpectrumResidual({ B, U }) {
  let worst = 0;
  for (let r = 0; r < 47; r++) {
    const vr = new Float64Array(125), vi = new Float64Array(125);
    for (let k = 0; k < 47; k++) { const a = -2 * Math.PI * r * k / 47, cr = Math.cos(a) / Math.sqrt(47), ci = Math.sin(a) / Math.sqrt(47); for (let x = 0; x < 125; x++) { vr[x] += B.a[x * 47 + k] * cr; vi[x] += B.a[x * 47 + k] * ci; } }
    const wr = Math.cos(2 * Math.PI * r / 47), wi = Math.sin(2 * Math.PI * r / 47);
    for (let x = 0; x < 125; x++) {
      let ur = 0, ui = 0; for (let y = 0; y < 125; y++) { ur += U.a[x * 125 + y] * vr[y]; ui += U.a[x * 125 + y] * vi[y]; }
      worst = Math.max(worst, Math.hypot(ur - (wr * vr[x] - wi * vi[x]), ui - (wr * vi[x] + wi * vr[x])));
    }
  }
  return worst;
}

// ---------------------------------------------------------------- torus clock on Z5^3 sites
export function torusClock() {
  const U = rmat(125);
  for (let n = 0; n < 125; n++) { const a = Math.floor(n / 25), b = Math.floor(n / 5) % 5, c = n % 5; const m = ((a + 1) % 5) * 25 + ((b + 1) % 5) * 5 + (c + 1) % 5; U.a[m * 125 + n] = 1; }
  return U;
}
export const permutationOrder = perm => { // order of a permutation given as an index map
  let k = 1, cur = perm.slice(); const id = cur.every((v, i) => v === i);
  if (id) return 1;
  while (true) { cur = cur.map(i => perm[i]); k++; if (cur.every((v, i) => v === i)) return k; if (k > 100000) return null; }
};
export const lcm = (a, b) => { const g = (x, y) => y ? g(y, x % y) : x; return a / g(a, b) * b; };

// ---------------------------------------------------------------- modular synchronisation clock
export const MODSYNC_M = 47 * 125; // 5875
export function multOrder(a, n) { let k = 1, x = a % n; while (x !== 1) { x = x * a % n; k++; if (k > n) return null; } return k; }
export function modsyncOrbit(len = 1200) { const z = [2]; for (let i = 0; i < len; i++) z.push(z[z.length - 1] * z[z.length - 1] % MODSYNC_M); return z; }
export function modsyncTransientPeriod(z) { const first = new Map(); for (let n = 0; n < z.length; n++) { if (first.has(z[n])) return { mu: first.get(z[n]), period: n - first.get(z[n]) }; first.set(z[n], n); } return null; }
export function periodMod(z, q) { for (let p = 1; p < 400; p++) { let ok = true; for (let n = 2; n < 600 && ok; n++) if (z[n + p] % q !== z[n] % q) ok = false; if (ok) return p; } return null; }
// z_n = 2^(2^n) mod 5875 with exact BigInt modular exponentiation
export function powmodBig(base, e, m) { let r = 1n, b = BigInt(base) % m; while (e > 0n) { if (e & 1n) r = r * b % m; b = b * b % m; e >>= 1n; } return r; }

// Γ in M-blocks (block-diagonal, exactly), and Γ^n − P spectral norms
export function gammaBlocks(car) {
  return car.blocks.map(blk => {
    const b = blk.sites.length, G = rmat(b), P = rmat(b);
    for (let i = 0; i < b; i++) for (let j = 0; j < b; j++) { G.a[i * b + j] = car.Gamma.a[blk.sites[i] * 125 + blk.sites[j]]; P.a[i * b + j] = car.P.a[blk.sites[i] * 125 + blk.sites[j]]; }
    return { sites: blk.sites, G, P };
  });
}
export function contractionNorms(car, ns) {
  const blocks = gammaBlocks(car), want = new Set(ns), out = {}, pw = blocks.map(b => reye(b.G.n));
  const nmax = Math.max(...ns);
  for (let n = 1; n <= nmax; n++) {
    for (let k = 0; k < blocks.length; k++) pw[k] = rmul(pw[k], blocks[k].G);
    if (want.has(n)) {
      let s = 0;
      for (let k = 0; k < blocks.length; k++) { const D = radd(pw[k], blocks[k].P, -1); const { values } = rsymEig(D); for (const v of values) s = Math.max(s, Math.abs(v)); }
      out[n] = s;
    }
  }
  return out;
}

// Joint operator U = Γ ⊗ T stepped one application at a time with the clock (T|z_n> = |z_{n+1}>),
// on a real 125 x 220 state: returns residuals of the synchronisation theorem.
export function modsyncJoint(car, { seed = 47, steps = 220 } = {}) {
  const blocks = gammaBlocks(car), W = 220;
  let s = seed >>> 0; const rnd = () => { s ^= s << 13; s >>>= 0; s ^= s >>> 17; s ^= s << 5; s >>>= 0; return s / 4294967296 - 0.5; };
  let X = new Float64Array(125 * W); for (let k = 0; k < X.length; k++) X[k] = rnd();
  let nn = 0; for (const x of X) nn += x * x; nn = Math.sqrt(nn); X = X.map(x => x / nn);
  const applyBlocks = (Y, mats) => { const out = new Float64Array(Y.length); mats.forEach((M, k) => { const st = blocks[k].sites, b = st.length; for (let i = 0; i < b; i++) for (let j = 0; j < b; j++) { const g = M.a[i * b + j]; if (g === 0) continue; const ro = st[i] * W, co = st[j] * W; for (let w = 0; w < W; w++) out[ro + w] += g * Y[co + w]; } }); return out; };
  const roll = (Y, r) => { const out = new Float64Array(Y.length); for (let x = 0; x < 125; x++) for (let w = 0; w < W; w++) out[x * W + ((w + r) % W + W) % W] = Y[x * W + w]; return out; };
  const PX = applyBlocks(X, blocks.map(b => b.P));
  let state = X, snaps = {};
  for (let n = 1; n <= steps; n++) { state = roll(applyBlocks(state, blocks.map(b => b.G)), 1); if (n === 110 || n === 220 || n === 20) snaps[n] = state; }
  // Γ^220 X by block powers (no clock)
  let pw = blocks.map(b => reye(b.G.n)); for (let n = 0; n < 220; n++) pw = pw.map((M, k) => rmul(M, blocks[k].G));
  const G220X = applyBlocks(X, pw);
  const dist = (A, B) => { let d = 0; for (let k = 0; k < A.length; k++) d += (A[k] - B[k]) ** 2; return Math.sqrt(d); };
  const dot = (A, B) => { let d = 0; for (let k = 0; k < A.length; k++) d += A[k] * B[k]; return d; };
  const w0 = dot(PX, PX), s220 = dot(snaps[220], snaps[220]);
  return {
    U220_minus_Gamma220: dist(snaps[220], G220X),
    U220_minus_P: dist(snaps[220], PX),
    U110_minus_PT110: dist(snaps[110], roll(PX, 110)),
    U110_minus_P: dist(snaps[110], PX),
    protectedWeight0: w0, success220: s220,
    returnFidelity220: dot(PX, snaps[220]) ** 2 / (w0 * s220),
  };
}

// ---------------------------------------------------------------- radix dialects
const DIGITS = { 2: '01', 5: '01234', 10: '0123456789', 12: '0123456789XE', 20: '0123456789ABCDEFGHIJ' };
// exact expansion of the rational p/q in base b: integer part, preperiod digits, repetend digits
export function radixExact(p, q, b) {
  let ip = Math.floor(p / q), r = p - ip * q; const seen = new Map(), digits = [];
  while (r !== 0 && !seen.has(r)) { seen.set(r, digits.length); r *= b; digits.push(Math.floor(r / q)); r %= q; }
  if (r === 0) return { int: ip, pre: digits, rep: [] };
  const k = seen.get(r); return { int: ip, pre: digits.slice(0, k), rep: digits.slice(k) };
}
export function radixString(p, q, b, maxRep = 24) {
  const { int, pre, rep } = radixExact(p, q, b);
  if (b === 60) { const f = a => a.join(','); return `${int};${f(pre)}${rep.length ? (pre.length ? ',' : '') + '(' + f(rep.slice(0, maxRep)) + (rep.length > maxRep ? ',…' : '') + ')' : ''}`; }
  const D = DIGITS[b], f = a => a.map(d => D[d]).join('');
  return `${int}.${f(pre)}${rep.length ? '(' + f(rep.slice(0, maxRep)) + (rep.length > maxRep ? '…' : '') + ')' : ''}`;
}
export function radixFloat(x, b, n = 6) { // first n fractional digits of x ∈ [0,1) in base b
  const out = []; let y = x; for (let k = 0; k < n; k++) { y *= b; const d = Math.floor(y + 1e-12); out.push(Math.min(b - 1, d)); y -= d; }
  return b === 60 ? '0;' + out.join(',') : '0.' + out.map(d => DIGITS[b][d]).join('');
}

// ---------------------------------------------------------------- readout helpers
export function unwrap(angles) { const out = []; let off = 0; for (let k = 0; k < angles.length; k++) { if (k) { const d = angles[k] + off - out[k - 1]; if (d > Math.PI) off -= 2 * Math.PI; else if (d < -Math.PI) off += 2 * Math.PI; } out.push(angles[k] + off); } return out; }
export const phaseTime = (coherences, omega) => unwrap(coherences.map(([r, i]) => -Math.atan2(i, r))).map(a => a / omega);
export const clickTime = (n, rate) => n / rate;
export { FILTER_DENOM };
