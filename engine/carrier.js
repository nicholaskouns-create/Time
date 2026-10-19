// The E47 carrier V2⊗V2⊗V2 (dim 125), built from scratch.
// Sites n = 25a + 5b + c carry m = (2-a, 2-b, 2-c) (descending Condon–Shortley order),
// the convention of the canonical spacetime runner and of the TQGNN reference.
import { rmat, reye, rmul, radd, rscale, rsymEig, rnormF, rmaxabs, rtrans, cmat, cmul, cadj, cnormF, hermEig, caxpy } from './linalg.js';

export const SPECTRUM = { 0: 1, 2: 9, 6: 25, 12: 28, 20: 27, 30: 22, 42: 13 };
export const SELECTED = [6, 30];
export const FILTER_DENOM = 99144;            // Γ = I − K²/99144,  ε* = 1/99144
export const REF_SITE = 38;                   // m = (1, 0, −1), total M = 0
export const M_SINGLE = [2, 1, 0, -1, -2];

export const siteM = n => [2 - Math.floor(n / 25), 2 - Math.floor(n / 5) % 5, 2 - n % 5];

// spin-2 ladder operator J+ in descending-m order: (J+)_{i,i+1} = sqrt(j(j+1) − m_{i+1}(m_{i+1}+1))
function jplus() {
  const J = rmat(5);
  for (let i = 0; i < 4; i++) { const mm = M_SINGLE[i + 1]; J.a[i * 5 + i + 1] = Math.sqrt(6 - mm * (mm + 1)); }
  return J;
}

function embed3(A, slot) { // A acting on factor `slot` of V2⊗V2⊗V2 (real 5x5 -> real 125x125)
  const out = rmat(125);
  for (let r = 0; r < 125; r++) {
    const d = [Math.floor(r / 25), Math.floor(r / 5) % 5, r % 5];
    for (let k = 0; k < 5; k++) {
      const v = A.a[d[slot] * 5 + k]; if (v === 0) continue;
      const e = d.slice(); e[slot] = k;
      out.a[r * 125 + (e[0] * 25 + e[1] * 5 + e[2])] += v;
    }
  }
  return out;
}

export class Carrier {
  constructor() {
    const Jp = jplus(), Jm = rtrans(Jp), Jz = rmat(5);
    for (let i = 0; i < 5; i++) Jz.a[i * 5 + i] = M_SINGLE[i];
    const tot = A => radd(radd(embed3(A, 0), embed3(A, 1)), embed3(A, 2));
    this.Jp = tot(Jp); this.Jm = tot(Jm); this.Jz = tot(Jz);
    // C = Jz² + (J+J− + J−J+)/2 : real symmetric and exactly M-conserving
    const C = radd(rmul(this.Jz, this.Jz), rscale(radd(rmul(this.Jp, this.Jm), rmul(this.Jm, this.Jp)), 0.5));
    for (let i = 0; i < 125; i++) for (let j = i + 1; j < 125; j++) { const s = 0.5 * (C.a[i * 125 + j] + C.a[j * 125 + i]); C.a[i * 125 + j] = s; C.a[j * 125 + i] = s; }
    this.C = C;
    this.dim = 125;
    this._spectral();
    const I = reye(125);
    this.K = rmul(radd(C, I, -6), radd(C, I, -30));
    this.K2 = rmul(this.K, this.K);
    this.P = this.projector(...SELECTED);
    this.Gamma = radd(I, this.K2, -1 / FILTER_DENOM);
  }

  // Block-diagonalise C by total M (13 blocks = 13 graph components), Jacobi in each block.
  _spectral() {
    const blocks = new Map();
    for (let n = 0; n < 125; n++) { const M = siteM(n).reduce((s, x) => s + x, 0); if (!blocks.has(M)) blocks.set(M, []); blocks.get(M).push(n); }
    this.blocks = [...blocks.entries()].sort((x, y) => y[0] - x[0]).map(([M, sites]) => ({ M, sites }));
    const vecs = rmat(125), vals = [];
    let col = 0, maxOffBlock = 0;
    for (let i = 0; i < 125; i++) for (let j = 0; j < 125; j++) {
      const Mi = siteM(i).reduce((s, x) => s + x, 0), Mj = siteM(j).reduce((s, x) => s + x, 0);
      if (Mi !== Mj) maxOffBlock = Math.max(maxOffBlock, Math.abs(this.C.a[i * 125 + j]));
    }
    this.maxOffBlock = maxOffBlock;
    for (const blk of this.blocks) {
      const b = blk.sites.length, B = rmat(b);
      for (let i = 0; i < b; i++) for (let j = 0; j < b; j++) B.a[i * b + j] = this.C.a[blk.sites[i] * 125 + blk.sites[j]];
      const { values, vectors } = rsymEig(B);
      blk.values = values;
      for (let c = 0; c < b; c++) {
        for (let i = 0; i < b; i++) vecs.a[blk.sites[i] * 125 + col] = vectors.a[i * b + c];
        vals.push(values[c]); col++;
      }
    }
    this.evals = vals;
    this.labels = vals.map(v => Math.round(v));
    this.integralityResidual = Math.max(...vals.map(v => Math.abs(v - Math.round(v))));
    this.evecs = vecs;
  }

  multiplicities() { const m = {}; for (const l of this.labels) m[l] = (m[l] || 0) + 1; return m; }

  basis(lam) { // 125 x d real basis of the eigenspace C = lam
    const cols = this.labels.map((l, k) => l === lam ? k : -1).filter(k => k >= 0), B = rmat(125, cols.length);
    cols.forEach((k, c) => { for (let r = 0; r < 125; r++) B.a[r * cols.length + c] = this.evecs.a[r * 125 + k]; });
    return B;
  }

  projector(...lams) {
    const P = rmat(125);
    for (const lam of lams) { const B = this.basis(lam); const BBt = rmul(B, rtrans(B)); for (let k = 0; k < P.a.length; k++) P.a[k] += BBt.a[k]; }
    return P;
  }

  functionOfC(f) { // f: eigenvalue -> number
    const out = rmat(125);
    for (const lam of Object.keys(SPECTRUM).map(Number)) { const P = this.projector(lam), v = f(lam); for (let k = 0; k < out.a.length; k++) out.a[k] += v * P.a[k]; }
    return out;
  }

  // graph view: C as a weighted graph on the 125 sites
  components() {
    const seen = new Array(125).fill(false), comps = [];
    for (let s = 0; s < 125; s++) {
      if (seen[s]) continue;
      const comp = [], queue = [s]; seen[s] = true;
      while (queue.length) {
        const u = queue.shift(); comp.push(u);
        for (let v = 0; v < 125; v++) if (v !== u && !seen[v] && Math.abs(this.C.a[u * 125 + v]) > 1e-12) { seen[v] = true; queue.push(v); }
      }
      comps.push(comp.sort((x, y) => x - y));
    }
    return comps;
  }
}

export function pairSingletCopy(k = 2) { // |pair singlet of spins 1,2> ⊗ |m3 = M_SINGLE[k]>  (inside C = 6)
  const v = new Float64Array(125);
  for (let ia = 0; ia < 5; ia++) {
    const ma = M_SINGLE[ia], ib = M_SINGLE.indexOf(-ma);
    v[(ia * 5 + ib) * 5 + k] = (-1) ** (2 - ma) / Math.sqrt(5);
  }
  return v;
}

// ---------------------------------------------------------------------------
// Background metric g = B4†(P6 − P0)B4 (port of the canonical runner's construction)

export function pullbackMetric(car) {
  // W: 125 x 5, columns pair-singlet ⊗ e_k
  const W = cmat(125, 5);
  for (let k = 0; k < 5; k++) { const v = pairSingletCopy(k); for (let r = 0; r < 125; r++) W.re[r * 5 + k] = v[r]; }
  // u: the Casimir singlet
  const u0 = car.basis(0); const u = new Float64Array(125); for (let r = 0; r < 125; r++) u[r] = u0.a[r];
  // so(3) generators on R^3 and the real spin-2 model Sym_0(3,R)
  const Qm = [[1, 0, 0], [0, 0, 0], [0, 0, -1]];
  const A = [
    [[0, 0, 0], [0, 0, -1], [0, 1, 0]],
    [[0, 0, 1], [0, 0, 0], [-1, 0, 0]],
    [[0, -1, 0], [1, 0, 0], [0, 0, 0]],
  ];
  const mm = (X, Y) => X.map((row, i) => Y[0].map((_, j) => row.reduce((s, x, k) => s + x * Y[k][j], 0)));
  const sub = (X, Y) => X.map((row, i) => row.map((x, j) => x - Y[i][j]));
  const tr = X => X[0][0] + X[1][1] + X[2][2];
  const D = A.map(Ai => sub(mm(Ai, Qm), mm(Qm, Ai)));
  const s2 = Math.SQRT2, s6 = Math.sqrt(6);
  const E = [
    [[1 / s2, 0, 0], [0, -1 / s2, 0], [0, 0, 0]],
    [[1 / s6, 0, 0], [0, 1 / s6, 0], [0, 0, -2 / s6]],
    [[0, 1 / s2, 0], [1 / s2, 0, 0], [0, 0, 0]],
    [[0, 0, 1 / s2], [0, 0, 0], [1 / s2, 0, 0]],
    [[0, 0, 0], [0, 0, 1 / s2], [0, 1 / s2, 0]],
  ];
  // rho[i][a][b] = Tr(E_a (A_i E_b − E_b A_i))
  const rho = A.map(Ai => E.map(Ea => E.map(Eb => tr(mm(Ea, sub(mm(Ai, Eb), mm(Eb, Ai)))))));
  // JW_i = W† J_i W with Jx = (J+ + J−)/2, Jy = (J+ − J−)/(2i), Jz
  const Jx = cmat(125), Jy = cmat(125), Jz = cmat(125);
  for (let k = 0; k < 125 * 125; k++) {
    Jx.re[k] = 0.5 * (car.Jp.a[k] + car.Jm.a[k]);
    Jy.im[k] = -0.5 * (car.Jp.a[k] - car.Jm.a[k]);
    Jz.re[k] = car.Jz.a[k];
  }
  const Wd = cadj(W);
  const JW = [Jx, Jy, Jz].map(J => cmul(Wd, cmul(J, W)));
  // Solve Phi rho_i + i JW_i Phi = 0 for Phi (5x5), vec column-major: vec(Phi) index = a*5 + r  (Phi[r][a])
  // Row block i: (rho_i^T ⊗ I5) + (I5 ⊗ i JW_i)
  const Ms = cmat(75, 25);
  for (let i = 0; i < 3; i++) for (let a = 0; a < 5; a++) for (let r = 0; r < 5; r++) {
    const row = i * 25 + a * 5 + r;
    for (let b = 0; b < 5; b++) Ms.re[row * 25 + b * 5 + r] += rho[i][b][a];          // (rho^T)_{a b}
    for (let s = 0; s < 5; s++) {                                                      // i JW_{r s}
      const jr = JW[i].re[r * 5 + s], ji = JW[i].im[r * 5 + s];
      Ms.re[row * 25 + a * 5 + s] += -ji; Ms.im[row * 25 + a * 5 + s] += jr;
    }
  }
  const G = cmul(cadj(Ms), Ms);
  const { values, vectors } = hermEig(G);
  const nullity = values.filter(v => Math.abs(v) < 1e-10).length;
  const phi = cmat(5);
  for (let a = 0; a < 5; a++) for (let r = 0; r < 5; r++) { phi.re[r * 5 + a] = vectors.re[(a * 5 + r) * 25 + 0]; phi.im[r * 5 + a] = vectors.im[(a * 5 + r) * 25 + 0]; }
  // intertwiner residual
  let resid = 0;
  for (let i = 0; i < 3; i++) {
    const R = cmat(5);
    const rhoC = cmat(5); for (let x = 0; x < 5; x++) for (let y = 0; y < 5; y++) rhoC.re[x * 5 + y] = rho[i][x][y];
    caxpy(1, 0, cmul(phi, rhoC), R);
    caxpy(0, 1, cmul(JW[i], phi), R);
    resid = Math.max(resid, cnormF(R));
  }
  let tn = 0; for (let k = 0; k < 25; k++) tn += phi.re[k] ** 2 + phi.im[k] ** 2;
  const sc = Math.sqrt(5 / tn); for (let k = 0; k < 25; k++) { phi.re[k] *= sc; phi.im[k] *= sc; }
  const dvecs = D.map(Di => E.map(Ea => tr(mm(Ea, Di))));
  const vs = dvecs.map(dv => { // W (Phi dv)
    const pr = new Float64Array(5), pi = new Float64Array(5);
    for (let r = 0; r < 5; r++) for (let a = 0; a < 5; a++) { pr[r] += phi.re[r * 5 + a] * dv[a]; pi[r] += phi.im[r * 5 + a] * dv[a]; }
    const vr = new Float64Array(125), vi = new Float64Array(125);
    for (let x = 0; x < 125; x++) for (let r = 0; r < 5; r++) { vr[x] += W.re[x * 5 + r] * pr[r]; vi[x] += W.re[x * 5 + r] * pi[r]; }
    return [vr, vi];
  });
  const B = [[u, new Float64Array(125)], ...vs];
  const P0 = car.projector(0), P6 = car.projector(6);
  const Dm = radd(P6, P0, -1);
  const g = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];
  let gim = 0;
  for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) {
    const [xr, xi] = B[i], [yr, yi] = B[j];
    let sr = 0, si = 0;
    for (let p = 0; p < 125; p++) {
      let tr_ = 0, ti_ = 0;
      for (let q = 0; q < 125; q++) { const d = Dm.a[p * 125 + q]; if (d !== 0) { tr_ += d * yr[q]; ti_ += d * yi[q]; } }
      sr += xr[p] * tr_ + xi[p] * ti_; si += xr[p] * ti_ - xi[p] * tr_;
    }
    g[i][j] = sr; gim = Math.max(gim, Math.abs(si));
  }
  return { g, intertwinerResidual: resid, nullity, imaginaryResidual: gim };
}

export const METRIC = [[-1, 0, 0, 0], [0, 2, 0, 0], [0, 0, 8, 0], [0, 0, 0, 2]];

// N = dτ/dt = sqrt(−g(ẋ, ẋ)) for ẋ = (1, v); null when the worldline is not timelike
export function lapse(v, g = METRIC) {
  const x = [1, v[0], v[1], v[2]];
  let s = 0; for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) s += x[i] * g[i][j] * x[j];
  return -s > 0 ? Math.sqrt(-s) : null;
}
export const MINKOWSKI = [[-1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];

export { rnormF, rmaxabs };
