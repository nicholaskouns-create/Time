// Temporal quantum graph neural network (TQGNN) of E47 clock-memory nodes.
// Each node carries the 4-dim unit (logical ⊗ flag) and its own worldline lapse N_v = dτ_v/dt.
// Coordinate-time generator:  H_t = Σ_v N_v H_v + Σ_e N_e H_e,   L_{v,t} = √N_v R_v,   N_e = (N_j + N_k)/2,
// H_e = J (σ+_j σ−_k + h.c.) + Jzz Z_j Z_k  (QGNN message passing between clocks).
// With equal lapses this is the time-chain identity dρ/dt = N L_τ ρ.
import { cmat, ceye, cmul, cadj, ckron, caxpy, cscale, cexpm, cfromRows } from './linalg.js';
import { reduced, liouvillian } from './node.js';

const A_TO_B = cfromRows([[0, 0, 0, 0], [0, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0]]); // |b><a| ⊗ I_flag  (index 2l+f)
const ZL = cfromRows([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]);
export function embed(op, v, n) { let out = ceye(1); for (let w = 0; w < n; w++) out = ckron(out, w === v ? op : ceye(4)); return out; }

function toCSR(A) {
  const rowPtr = [0], col = [], re = [], im = [];
  for (let i = 0; i < A.n; i++) { for (let j = 0; j < A.m; j++) { const k = i * A.m + j; if (A.re[k] !== 0 || A.im[k] !== 0) { col.push(j); re.push(A.re[k]); im.push(A.im[k]); } } rowPtr.push(col.length); }
  return { n: A.n, rowPtr: Int32Array.from(rowPtr), col: Int32Array.from(col), re: Float64Array.from(re), im: Float64Array.from(im) };
}

export class Network {
  constructor({ n = 2, edges = [[0, 1]], omega = 2, c = 1, gamma = 20, J = 0, Jzz = 0, clock = 'matched', syndrome = [12, 42], lapses = null }) {
    this.n = n; this.edges = edges; this.omega = omega; this.gamma = gamma; this.D = 4 ** n;
    this.lapses = lapses || new Array(n).fill(1);
    const { H, R } = reduced({ omega, c, gamma, clock, syndrome });
    const Nv = this.lapses, edgeLapse = (j, k) => 0.5 * (Nv[j] + Nv[k]);
    let Ht = cmat(this.D);
    for (let v = 0; v < n; v++) caxpy(Nv[v], 0, embed(H, v, n), Ht);
    this.Hedge = edges.map(([j, k]) => {
      const sp_j = embed(A_TO_B, j, n), sp_k = embed(A_TO_B, k, n), hop = cmul(sp_j, cadj(sp_k));
      const He = cscale(hop, J); caxpy(J, 0, cadj(hop), He); caxpy(Jzz, 0, cmul(embed(ZL, j, n), embed(ZL, k, n)), He);
      return He;
    });
    edges.forEach(([j, k], e) => caxpy(edgeLapse(j, k), 0, this.Hedge[e], Ht));
    this.H = Ht;
    this.L = gamma > 0 ? [...Array(n).keys()].map(v => cscale(embed(R, v, n), Math.sqrt(Nv[v]))) : [];
    // H_eff = H − (i/2) Σ L†L ; each L has at most one nonzero per row: store (k(i), value)
    this.Heff = cmat(this.D); this.Heff.re.set(Ht.re); this.Heff.im.set(Ht.im);
    for (const L of this.L) { const LdL = cmul(cadj(L), L); caxpy(0, -0.5, LdL, this.Heff); }
    this.HeffCSR = toCSR(this.Heff);
    this.jumps = this.L.map(L => { const k = new Int32Array(this.D).fill(-1), vr = new Float64Array(this.D), vi = new Float64Array(this.D); for (let i = 0; i < this.D; i++) for (let j = 0; j < this.D; j++) { const q = i * this.D + j; if (L.re[q] !== 0 || L.im[q] !== 0) { k[i] = j; vr[i] = L.re[q]; vi[i] = L.im[q]; } } return { k, vr, vi }; });
  }

  // GKSL right-hand side for ρ given as {re, im} Float64Array(D*D)
  rhs(rho) {
    const D = this.D, { rowPtr, col, re, im } = this.HeffCSR, Xr = new Float64Array(D * D), Xi = new Float64Array(D * D);
    for (let i = 0; i < D; i++) for (let p = rowPtr[i]; p < rowPtr[i + 1]; p++) {
      const k = col[p], hr = re[p], hi = im[p], io = i * D, ko = k * D;
      for (let j = 0; j < D; j++) { const xr = rho.re[ko + j], xi = rho.im[ko + j]; Xr[io + j] += hr * xr - hi * xi; Xi[io + j] += hr * xi + hi * xr; }
    }
    const outR = new Float64Array(D * D), outI = new Float64Array(D * D);
    for (let i = 0; i < D; i++) for (let j = 0; j < D; j++) { // −i(X − X†)
      const a = i * D + j, b = j * D + i, dr = Xr[a] - Xr[b], di = Xi[a] + Xi[b];
      outR[a] = di; outI[a] = -dr;
    }
    for (const { k, vr, vi } of this.jumps) for (let i = 0; i < D; i++) { if (k[i] < 0) continue; for (let j = 0; j < D; j++) { if (k[j] < 0) continue; // r_i ρ[k_i,k_j] conj(r_j)
      const q = k[i] * D + k[j], xr = rho.re[q], xi = rho.im[q], ar = vr[i] * xr - vi[i] * xi, ai = vr[i] * xi + vi[i] * xr;
      outR[i * D + j] += ar * vr[j] + ai * vi[j]; outI[i * D + j] += ai * vr[j] - ar * vi[j]; } }
    return { re: outR, im: outI };
  }
  rk4(rho, h) {
    const ax = (y, s, k) => ({ re: y.re.map((x, q) => x + s * k.re[q]), im: y.im.map((x, q) => x + s * k.im[q]) });
    const k1 = this.rhs(rho), k2 = this.rhs(ax(rho, h / 2, k1)), k3 = this.rhs(ax(rho, h / 2, k2)), k4 = this.rhs(ax(rho, h, k3));
    return { re: rho.re.map((x, q) => x + h / 6 * (k1.re[q] + 2 * k2.re[q] + 2 * k3.re[q] + k4.re[q])), im: rho.im.map((x, q) => x + h / 6 * (k1.im[q] + 2 * k2.im[q] + 2 * k3.im[q] + k4.im[q])) };
  }
  stableStep() { const scale = Math.max(Math.abs(this.omega) * 2 * Math.max(...this.lapses), this.gamma * Math.max(...this.lapses), 1); return 0.05 / scale; }
  // an upper bound on ‖L‖ in the max-row-sum sense, for choosing Taylor sub-steps
  normBound() {
    let hrow = 0; const { rowPtr, re, im } = this.HeffCSR;
    for (let i = 0; i < this.D; i++) { let s = 0; for (let p = rowPtr[i]; p < rowPtr[i + 1]; p++) s += Math.hypot(re[p], im[p]); hrow = Math.max(hrow, s); }
    let jump = 0; for (const { vr, vi } of this.jumps) { let m = 0; for (let i = 0; i < this.D; i++) m = Math.max(m, vr[i] ** 2 + vi[i] ** 2); jump += m; }
    return 2 * hrow + jump;
  }
  // ρ(t) = e^{L t} ρ by Taylor series on sub-steps of norm ≤ 1/2 (exact to round-off; L applied through rhs)
  propagate(rho, t) {
    if (t === 0) return { re: rho.re.slice(), im: rho.im.slice() };
    const s = Math.max(1, Math.ceil(2 * this.normBound() * Math.abs(t))), h = t / s;
    let cur = { re: rho.re.slice(), im: rho.im.slice() };
    for (let step = 0; step < s; step++) {
      let term = cur, sum = { re: cur.re.slice(), im: cur.im.slice() };
      for (let k = 1; k <= 80; k++) {
        const r = this.rhs(term); let tmax = 0, smax = 0;
        for (let q = 0; q < r.re.length; q++) { r.re[q] *= h / k; r.im[q] *= h / k; sum.re[q] += r.re[q]; sum.im[q] += r.im[q]; tmax = Math.max(tmax, Math.abs(r.re[q]), Math.abs(r.im[q])); smax = Math.max(smax, Math.abs(sum.re[q]), Math.abs(sum.im[q])); }
        term = r;
        if (tmax <= 1e-18 * smax) break;
      }
      cur = sum;
    }
    return cur;
  }
  // exact coordinate-time Liouvillian (for small networks)
  liouvillianDense() { return liouvillian(this.H, this.L); }

  // readouts: per node retention, coherence Tr(|a><b|_v ρ), coordinate-time click intensity
  readout(rho) {
    const D = this.D, n = this.n, out = [];
    for (let v = 0; v < n; v++) {
      const w = 4 ** (n - 1 - v); let ret = 0, cr = 0, ci = 0, leaked = 0;
      for (let i = 0; i < D; i++) {
        const digit = Math.floor(i / w) % 4, l = digit >> 1, f = digit & 1;
        if (f === 0) ret += rho.re[i * D + i]; else leaked += rho.re[i * D + i];
        if (l === 0) { const ib = i + 2 * w; cr += rho.re[ib * D + i]; ci += rho.im[ib * D + i]; }
      }
      out.push({ retention: ret, leaked, coherence: [cr, ci], clicks: this.lapses[v] * this.gamma * leaked });
    }
    return out;
  }

  static productState(nodeVectors) { // nodeVectors: arrays of 4 [re,im]
    let psi = [[1, 0]];
    for (const v of nodeVectors) { const out = []; for (const x of psi) for (const y of v) out.push([x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0]]); psi = out; }
    const nn = Math.sqrt(psi.reduce((s, z) => s + z[0] ** 2 + z[1] ** 2, 0));
    return psi.map(z => [z[0] / nn, z[1] / nn]);
  }
  static clockVector(theta = Math.PI / 2, phase = 0) { return [[Math.cos(theta / 2), 0], [0, 0], [Math.cos(phase) * Math.sin(theta / 2), Math.sin(phase) * Math.sin(theta / 2)], [0, 0]]; }
  static density(psi) { const D = psi.length, re = new Float64Array(D * D), im = new Float64Array(D * D); for (let i = 0; i < D; i++) for (let j = 0; j < D; j++) { re[i * D + j] = psi[i][0] * psi[j][0] + psi[i][1] * psi[j][1]; im[i * D + j] = psi[i][1] * psi[j][0] - psi[i][0] * psi[j][1]; } return { re, im }; }

  // pure-state jump trajectories (η = 1): exact no-jump propagator over δ, jumps located by bisection
  trajectory(psi0, tEnd, delta, rng) {
    const D = this.D, E = cexpm(cscale(this.Heff, 0, -delta));
    let psi = { re: Float64Array.from(psi0.map(z => z[0])), im: Float64Array.from(psi0.map(z => z[1])) }, t = 0, r = rng(), clicks = [];
    const mv = (M, x) => { const yr = new Float64Array(D), yi = new Float64Array(D); for (let i = 0; i < D; i++) { let sr = 0, si = 0; for (let j = 0; j < D; j++) { const mr = M.re[i * D + j], mi = M.im[i * D + j]; sr += mr * x.re[j] - mi * x.im[j]; si += mr * x.im[j] + mi * x.re[j]; } yr[i] = sr; yi[i] = si; } return { re: yr, im: yi }; };
    const nrm = x => { let s = 0; for (let i = 0; i < D; i++) s += x.re[i] ** 2 + x.im[i] ** 2; return s; };
    while (t < tEnd - 1e-12) {
      const nxt = mv(E, psi);
      if (nrm(nxt) >= r) { psi = nxt; t += delta; continue; }
      let lo = 0, hi = delta, best = null;
      for (let it = 0; it < 40; it++) { const mid = 0.5 * (lo + hi), y = mv(cexpm(cscale(this.Heff, 0, -mid)), psi); if (nrm(y) >= r) lo = mid; else { hi = mid; best = y; } }
      const y = best || nxt; t += hi;
      const w = this.jumps.map(({ k, vr, vi }) => { let s = 0; for (let i = 0; i < D; i++) if (k[i] >= 0) s += (vr[i] ** 2 + vi[i] ** 2) * (y.re[k[i]] ** 2 + y.im[k[i]] ** 2); return s; });
      let u = rng() * w.reduce((s, x) => s + x, 0), ch = 0; while (ch < w.length - 1 && (u -= w[ch]) > 0) ch++;
      const { k, vr, vi } = this.jumps[ch], z = { re: new Float64Array(D), im: new Float64Array(D) };
      for (let i = 0; i < D; i++) if (k[i] >= 0) { const xr = y.re[k[i]], xi = y.im[k[i]]; z.re[i] = vr[i] * xr - vi[i] * xi; z.im[i] = vr[i] * xi + vi[i] * xr; }
      const zn = Math.sqrt(nrm(z)); for (let i = 0; i < D; i++) { z.re[i] /= zn; z.im[i] /= zn; }
      psi = z; clicks.push({ t, node: ch }); r = rng();
    }
    const nn = Math.sqrt(nrm(psi)); for (let i = 0; i < D; i++) { psi.re[i] /= nn; psi.im[i] /= nn; }
    return { psi, clicks };
  }
}

// An incrementally advanced pure-state trajectory of a network (η = 1), for animation.
export class NetworkTrajectory {
  constructor(net, psi0, delta, rng) {
    this.net = net; this.delta = delta; this.rng = rng; this.D = net.D;
    this.E = cexpm(cscale(net.Heff, 0, -delta));
    this.psi = { re: Float64Array.from(psi0.map(z => z[0])), im: Float64Array.from(psi0.map(z => z[1])) };
    this.r = rng(); this.t = 0; this.acc = 0; this.clicks = [];
  }
  _mv(M, x) { const D = this.D, yr = new Float64Array(D), yi = new Float64Array(D); for (let i = 0; i < D; i++) { let sr = 0, si = 0; const o = i * D; for (let j = 0; j < D; j++) { const mr = M.re[o + j], mi = M.im[o + j]; sr += mr * x.re[j] - mi * x.im[j]; si += mr * x.im[j] + mi * x.re[j]; } yr[i] = sr; yi[i] = si; } return { re: yr, im: yi }; }
  _nrm(x) { let s = 0; for (let i = 0; i < this.D; i++) s += x.re[i] ** 2 + x.im[i] ** 2; return s; }
  advance(dt) {
    let left = dt + this.acc; const fresh = [];
    while (left >= this.delta - 1e-15) {
      const nxt = this._mv(this.E, this.psi);
      if (this._nrm(nxt) >= this.r) { this.psi = nxt; }
      else {
        let lo = 0, hi = this.delta, at = nxt;
        for (let it = 0; it < 30; it++) { const mid = 0.5 * (lo + hi), y = this._mv(cexpm(cscale(this.net.Heff, 0, -mid)), this.psi); if (this._nrm(y) >= this.r) lo = mid; else { hi = mid; at = y; } }
        const w = this.net.jumps.map(({ k, vr, vi }) => { let s = 0; for (let i = 0; i < this.D; i++) if (k[i] >= 0) s += (vr[i] ** 2 + vi[i] ** 2) * (at.re[k[i]] ** 2 + at.im[k[i]] ** 2); return s; });
        let u = this.rng() * w.reduce((s, x) => s + x, 0), ch = 0; while (ch < w.length - 1 && (u -= w[ch]) > 0) ch++;
        const { k, vr, vi } = this.net.jumps[ch], z = { re: new Float64Array(this.D), im: new Float64Array(this.D) };
        for (let i = 0; i < this.D; i++) if (k[i] >= 0) { const xr = at.re[k[i]], xi = at.im[k[i]]; z.re[i] = vr[i] * xr - vi[i] * xi; z.im[i] = vr[i] * xi + vi[i] * xr; }
        const zn = Math.sqrt(this._nrm(z)); for (let i = 0; i < this.D; i++) { z.re[i] /= zn; z.im[i] /= zn; }
        const rest = this.delta - hi; this.psi = rest > 1e-15 ? this._mv(cexpm(cscale(this.net.Heff, 0, -rest)), z) : z;
        const click = { t: this.t + hi, node: ch }; this.clicks.push(click); fresh.push(click); this.r = this.rng();
      }
      this.t += this.delta; left -= this.delta;
    }
    this.acc = left;
    return fresh;
  }
  density() { // normalised |ψ><ψ| in the {re, im} layout used by Network.readout
    const D = this.D, n = Math.sqrt(this._nrm(this.psi)), re = new Float64Array(D * D), im = new Float64Array(D * D), a = this.psi.re.map(x => x / n), b = this.psi.im.map(x => x / n);
    for (let i = 0; i < D; i++) for (let j = 0; j < D; j++) { re[i * D + j] = a[i] * a[j] + b[i] * b[j]; im[i * D + j] = b[i] * a[j] - a[i] * b[j]; }
    return { re, im };
  }
}

export function xorshift(seed = 47125) { let x = (seed >>> 0) || 1; return () => { x ^= x << 13; x >>>= 0; x ^= x >>> 17; x ^= x << 5; x >>>= 0; return x / 4294967296; }; }
