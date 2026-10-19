// The running Quantum Clock: one E47 clock-memory node on a worldline, propagated exactly.
// Deterministic state: ρ(τ + δ) = e^{Lδ} ρ(τ).  Stochastic state (detection efficiency η):
// ρ̃ evolves under e^{(L − ηJ)δ} (its trace is the no-click probability); a detected click applies J.
// At η = 0 nothing is recorded and the conditional state equals ρ: stochastic becomes dynamic.
import { cmat, cmul, cexpm, cscale, csub, caxpy } from './linalg.js';
import { reduced, liouvillian, jumpSuper, applySuper, clockState, readout, branch, countingStats, detunings, stationaryLeak } from './node.js';
import { lapse, METRIC, MINKOWSKI } from './carrier.js';
import { xorshift } from './network.js';

export const DEFAULTS = { omega: 2 * Math.PI / 8, c: 1, gamma: 20, clock: 'matched', syndrome: [12, 42], eta: 1, velocity: [0, 0, 0], metric: 'e47', theta: Math.PI / 2, phase: 0, seed: 47125 };

const traceRe = rho => rho.re[0] + rho.re[5] + rho.re[10] + rho.re[15];
const scaled = (rho, s) => { const r = cmat(4); for (let k = 0; k < 16; k++) { r.re[k] = rho.re[k] * s; r.im[k] = rho.im[k] * s; } return r; };

export class NodeClock {
  constructor(params = {}) { this.configure(params); this.reset(); }

  configure(params) {
    const p = { ...DEFAULTS, ...(this.p || {}), ...params };
    const g = p.metric === 'minkowski' ? MINKOWSKI : METRIC, N = lapse(p.velocity, g);
    if (N === null) throw Error('This velocity is not timelike for the chosen metric.');
    if (!(p.omega > 0) || !(p.gamma >= 0) || !(p.c >= 0) || !(p.eta >= 0 && p.eta <= 1)) throw Error('Parameters out of range: ω > 0, γ ≥ 0, c ≥ 0, 0 ≤ η ≤ 1.');
    this.p = p; this.N = N;
    const { H, R } = reduced(p);
    this.H = H; this.R = R; this.Lv = liouvillian(H, [R]); this.J = jumpSuper(R);
    const period = 2 * Math.PI / Math.max(1e-9, Math.abs(p.omega));
    this.delta = Math.min(0.02, Math.max(1e-4, period / 720));
    this.Sd = cexpm(cscale(this.Lv, this.delta));
    this.Lne = csub(this.Lv, cscale(this.J, p.eta));
    this.Sne = cexpm(cscale(this.Lne, this.delta));
    const [da, db] = detunings(p.omega, p.clock, p.syndrome);
    this.detuning = [da, db];
    this.branchStats = [0, 1].map(l => { const b = branch(l, p); return countingStats(b.H, [b.R], 0, p.eta); });
    this.leak = [stationaryLeak(p.c, p.gamma, da), stationaryLeak(p.c, p.gamma, db)];
  }

  reset() {
    const p = this.p;
    this.rho = clockState(p.theta, p.phase, 0);
    this.rhoc = clockState(p.theta, p.phase, 0);
    this.tau = 0; this.t = 0; this.acc = 0; this.clicks = []; this.rng = xorshift(p.seed); this.r = this.rng();
    const c0 = readout(this.rho).coherence; this.phi0 = -Math.atan2(c0[1], c0[0]); this.phiDet = 0; this.phiSto = 0;
    this.lastDet = c0; this.lastSto = c0;
  }

  // advance by coordinate time dt (proper time N dt); returns clicks recorded in this call
  advance(dt) {
    let dtau = this.N * dt + this.acc; const n = Math.floor(dtau / this.delta + 1e-12); this.acc = dtau - n * this.delta;
    const fresh = [];
    for (let s = 0; s < n; s++) {
      this.rho = applySuper(this.Sd, this.rho);
      const next = applySuper(this.Sne, this.rhoc);
      if (this.p.eta > 0 && traceRe(next) < this.r) {
        // locate the click inside the step: trace(e^{(L−ηJ)x} ρ̃) = r
        let lo = 0, hi = this.delta, at = next;
        for (let it = 0; it < 32; it++) { const mid = 0.5 * (lo + hi), y = applySuper(cexpm(cscale(this.Lne, mid)), this.rhoc); if (traceRe(y) >= this.r) lo = mid; else { hi = mid; at = y; } }
        const jumped = applySuper(this.J, at), tr = traceRe(jumped);
        this.rhoc = scaled(jumped, 1 / tr);
        const rest = this.delta - hi; if (rest > 0) this.rhoc = applySuper(cexpm(cscale(this.Lne, rest)), this.rhoc);
        const click = { tau: this.tau + hi, t: (this.tau + hi) / this.N }; this.clicks.push(click); fresh.push(click);
        this.r = this.rng();
      } else this.rhoc = next;
      this.tau += this.delta;
      this._track();
    }
    this.t += dt;
    return fresh;
  }

  _track() { // continuous unwrapping of both clock hands
    const unwrapStep = (prev, cur, acc) => { const d = Math.atan2(prev[0] * cur[1] - prev[1] * cur[0], prev[0] * cur[0] + prev[1] * cur[1]); return acc - d; };
    const cd = readout(this.rho).coherence, cs = this.conditional().coherence;
    this.phiDet = unwrapStep(this.lastDet, cd, this.phiDet); this.lastDet = cd;
    this.phiSto = unwrapStep(this.lastSto, cs, this.phiSto); this.lastSto = cs;
  }

  conditional() { const tr = traceRe(this.rhoc); return readout(scaled(this.rhoc, 1 / tr)); }

  state() {
    const det = readout(this.rho), sto = this.conditional(), p = this.p;
    const rate = det.pa * this.branchStats[0].rate + det.pb * this.branchStats[1].rate;   // detected clicks per unit proper time
    return {
      t: this.t, tau: this.tau, lapse: this.N,
      deterministic: det, conditional: sto,
      phaseTime: this.phiDet / p.omega, phaseTimeConditional: this.phiSto / p.omega,
      clicks: this.clicks.length, clickRate: rate, clickTime: rate > 0 ? this.clicks.length / rate : null,
      hand: this.phi0 + this.phiDet, handConditional: this.phi0 + this.phiSto,
    };
  }
}
