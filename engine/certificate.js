// The Quantum Clock precision certificate: every claim the page makes, recomputed.
// Evidence classes follow the E47 convention: E0 exact proof, E1 executable reconstruction,
// E2 simulation. Each check records value, target, residual, tolerance.
import * as L from './linalg.js';
import { Carrier, SPECTRUM, pullbackMetric, lapse, METRIC } from './carrier.js';
import * as Nd from './node.js';
import { Network, xorshift } from './network.js';
import * as K from './clocks.js';
import * as Rc from './receipt.js';

export const ENGINE_VERSION = 'qclock-engine 1.0.0';
export const CERTIFICATE_ID = 'MC-QCLOCK-PRECISION-20261019-001';

const fr = (p, q) => { const g = (a, b) => b ? g(b, a % b) : Math.abs(a); const d = g(p, q); return [p / d, q / d]; };

export async function runCertificate({ golden = null, onProgress = null, environment = {} } = {}) {
  const checks = [];
  const t0 = (globalThis.performance || Date).now();
  const add = (group, id, title, evidence, { value, target, residual = null, tolerance = null, pass }) => {
    checks.push({ group, id, title, evidence, value, target, residual, tolerance, pass: !!pass });
    if (onProgress) onProgress(checks.length, title);
  };
  const tick = () => new Promise(r => setTimeout(r, 0));

  // ---------------------------------------------------------------- A. carrier
  const car = new Carrier();
  const mult = car.multiplicities();
  add('carrier', 'A1', 'Casimir spectrum of V₂⊗³: multiplicities 1, 9, 25, 28, 27, 22, 13', 'E0+E1',
    { value: mult, target: SPECTRUM, pass: Object.keys(SPECTRUM).every(k => mult[k] === SPECTRUM[k]) });
  add('carrier', 'A2', 'Casimir eigenvalues are the integers j(j+1)', 'E1', { value: car.integralityResidual, target: 0, residual: car.integralityResidual, tolerance: 1e-12, pass: car.integralityResidual < 1e-12 });
  const comps = car.components().map(c => c.length);
  add('carrier', 'A3', 'C is a graph operator with 13 components (total-M blocks)', 'E0+E1',
    { value: comps.join(','), target: '1,3,6,10,15,18,19,18,15,10,6,3,1', residual: car.maxOffBlock, tolerance: 0, pass: comps.join(',') === '1,3,6,10,15,18,19,18,15,10,6,3,1' && car.maxOffBlock === 0 });
  const P = car.P, idem = L.rmaxabs(L.radd(L.rmul(P, P), P, -1)), sym = L.rmaxabs(L.radd(P, L.rtrans(P), -1)), trP = L.rtrace(P);
  add('carrier', 'A4', 'P47 is the orthogonal projector of rank 47', 'E0+E1', { value: { trace: trP, idempotence: idem, symmetry: sym }, target: { trace: 47 }, residual: Math.max(Math.abs(trP - 47), idem, sym), tolerance: 1e-12, pass: Math.abs(trP - 47) < 1e-12 && idem < 1e-13 && sym < 1e-14 });
  const KP = L.rmaxabs(L.rmul(car.K, P));
  add('carrier', 'A5', 'K P47 = 0 with K = (C − 6I)(C − 30I)', 'E0+E1', { value: KP, target: 0, residual: KP, tolerance: 1e-11, pass: KP < 1e-11 });
  const qvals = car.evals.filter((v, k) => ![6, 30].includes(car.labels[k])).map(v => ((v - 6) * (v - 30)) ** 2);
  const qmin = Math.min(...qvals), qmax = Math.max(...qvals);
  add('carrier', 'A6', 'Q = K² has nonzero spectrum from 11664 to 186624', 'E0+E1', { value: [qmin, qmax], target: [11664, 186624], residual: Math.max(Math.abs(qmin - 11664), Math.abs(qmax - 186624)), tolerance: 1e-8, pass: Math.abs(qmin - 11664) < 1e-8 && Math.abs(qmax - 186624) < 1e-8 });
  // exact rational: ε* = 2/(qmin + qmax), ρ* = (qmax − qmin)/(qmax + qmin)
  const eps = fr(2, 11664 + 186624), rhoStar = fr(186624 - 11664, 186624 + 11664);
  const rhoNum = Math.max(...qvals.map(q => Math.abs(1 - q / 99144)));
  add('carrier', 'A7', 'Optimal step ε* = 1/99144 and contraction ρ* = 15/17 (exact rationals)', 'E0',
    { value: { epsilon: `${eps[0]}/${eps[1]}`, rho: `${rhoStar[0]}/${rhoStar[1]}`, rho_numeric: rhoNum }, target: { epsilon: '1/99144', rho: '15/17' }, residual: Math.abs(rhoNum - 15 / 17), tolerance: 1e-12, pass: eps[0] === 1 && eps[1] === 99144 && rhoStar[0] === 15 && rhoStar[1] === 17 && Math.abs(rhoNum - 15 / 17) < 1e-12 });
  const ns = [1, 5, 20, 100, 220], norms = K.contractionNorms(car, ns);
  const contr = ns.map(n => ({ n, matrix_power: norms[n], exact: (15 / 17) ** n, diff: Math.abs(norms[n] - (15 / 17) ** n) }));
  add('carrier', 'A8', '‖Γⁿ − P‖₂ = (15/17)ⁿ by matrix powers, n = 1 … 220', 'E0+E1', { value: contr, target: '(15/17)^n', residual: Math.max(...contr.map(r => r.diff)), tolerance: '1e-13 + 1e-10·(15/17)^n', pass: contr.every(r => r.diff < 1e-13 + 1e-10 * r.exact) });
  add('carrier', 'A9', 'Coherence fraction Ωc = 47/125 = 0.376', 'E0', { value: trP / 125, target: 0.376, residual: Math.abs(trP / 125 - 0.376), tolerance: 1e-14, pass: Math.abs(trP / 125 - 0.376) < 1e-14 });
  await tick();

  // ---------------------------------------------------------------- B. spacetime
  const pm = pullbackMetric(car), gTarget = METRIC;
  const gRes = Math.max(...pm.g.flatMap((row, i) => row.map((x, j) => Math.abs(x - gTarget[i][j]))));
  add('spacetime', 'B1', 'Pullback metric g = B4†(P6 − P0)B4 = diag(−1, 2, 8, 2), built from the carrier', 'E0+E1', { value: pm.g, target: gTarget, residual: Math.max(gRes, pm.intertwinerResidual, pm.imaginaryResidual), tolerance: 1e-12, pass: gRes < 1e-12 && pm.nullity === 1 && pm.intertwinerResidual < 1e-12 });
  const vs = [[0.3, 0, 0], [0, 0.2, 0], [0.1, 0.1, 0.2], [0.8, 0, 0]];
  const lap = vs.map(v => lapse(v)), lapExact = vs.map(v => { const s = 1 - 2 * v[0] ** 2 - 8 * v[1] ** 2 - 2 * v[2] ** 2; return s > 0 ? Math.sqrt(s) : null; });
  add('spacetime', 'B2', 'Lapse N = √(1 − 2v₁² − 8v₂² − 2v₃²); spacelike velocities are refused', 'E0', { value: lap, target: lapExact, residual: Math.max(...lap.slice(0, 3).map((x, k) => Math.abs(x - lapExact[k]))), tolerance: 1e-15, pass: lap[3] === null && lap.slice(0, 3).every((x, k) => Math.abs(x - lapExact[k]) < 1e-15) });
  if (golden) {
    const gp = Math.max(...pm.g.flatMap((row, i) => row.map((x, j) => Math.abs(x - golden.metric.g[i][j]))));
    add('spacetime', 'B3', 'Parity with the Python reference: metric', 'E1', { value: gp, target: 0, residual: gp, tolerance: 1e-12, pass: gp < 1e-12 });
  }
  await tick();

  // ---------------------------------------------------------------- C. node
  const node = new Nd.ClockMemoryNode(car);
  const base = { omega: 2, c: 1, gamma: 20 };
  const redL = node.reductionResiduals({ ...base, clock: 'linear' }), redM = node.reductionResiduals({ ...base, clock: 'matched' });
  add('node', 'C1', 'Clock-memory basis [a, a′, b, b′] is orthonormal', 'E1', { value: redL.orthonormality, target: 0, residual: redL.orthonormality, tolerance: 1e-14, pass: redL.orthonormality < 1e-14 });
  const Pv = v => L.rmatvec(P, v), dist = (x, y) => Math.sqrt(x.reduce((s, z, k) => s + (z - y[k]) ** 2, 0)), nrm = x => Math.sqrt(x.reduce((s, z) => s + z * z, 0));
  const cval = v => { const Cv = L.rmatvec(car.C, v); return v.reduce((s, z, k) => s + z * Cv[k], 0); };
  const inE = Math.max(dist(Pv(node.a), node.a), dist(Pv(node.b), node.b)), outE = Math.max(nrm(Pv(node.as)), nrm(Pv(node.bs)));
  const cv = [node.a, node.as, node.b, node.bs].map(cval);
  add('node', 'C2', 'a, b lie in E47 (C = 6, 30); a′, b′ lie outside (C = 12, 42)', 'E0+E1', { value: { casimir: cv, in_E47: inE, outside: outE }, target: { casimir: [6, 12, 30, 42] }, residual: Math.max(inE, outE, ...cv.map((x, k) => Math.abs(x - [6, 12, 30, 42][k]))), tolerance: 1e-12, pass: inE < 1e-13 && outE < 1e-13 && cv.every((x, k) => Math.abs(x - [6, 12, 30, 42][k]) < 1e-12) });
  const inv = Math.max(redL.H_invariance, redL.R_invariance, redL.Rdag_invariance, redM.H_invariance, redM.R_invariance, redM.Rdag_invariance);
  add('node', 'C3', 'H, R and R† map span{a, a′, b, b′} into itself (linear and matched clocks)', 'E0+E1', { value: inv, target: 0, residual: inv, tolerance: 1e-13, pass: inv < 1e-13 });
  const red = Math.max(redL.H_reduced, redL.R_reduced, redM.H_reduced, redM.R_reduced);
  add('node', 'C4', 'The 125-dim generator restricts exactly to the 4-dim node: U†HU = H_red, U†RU = R_red', 'E0+E1', { value: red, target: 0, residual: red, tolerance: 1e-13, pass: red < 1e-13 });
  { // Clock Math: kernel-aligned decoherence-free theorem on the carrier
    const rng = xorshift(4747), B = car.basis(6), B2 = car.basis(30), Bm = L.rmat(125, 47);
    for (let r = 0; r < 125; r++) { for (let c = 0; c < 25; c++) Bm.a[r * 47 + c] = B.a[r * 25 + c]; for (let c = 0; c < 22; c++) Bm.a[r * 47 + 25 + c] = B2.a[r * 22 + c]; }
    const X = L.rmat(47); for (let k = 0; k < X.a.length; k++) X.a[k] = rng() - 0.5;
    let rho = L.rmul(Bm, L.rmul(L.rmul(X, L.rtrans(X)), L.rtrans(Bm))); const tr = L.rtrace(rho); rho = L.rscale(rho, 1 / tr);
    const D = L.radd(L.rmul(car.K, L.rmul(rho, car.K)), L.rscale(L.radd(L.rmul(car.K2, rho), L.rmul(rho, car.K2)), 0.5), -1);
    const Hc = L.rscale(car.C, 2 / 24), comm = L.rmaxabs(L.radd(L.rmul(Hc, P), L.rmul(P, Hc), -1));
    const dmax = L.rmaxabs(D);
    add('node', 'C5', 'Clock Math theorem: with L = K the dissipator vanishes on ρ = P47ρP47, and [ωC/24, P47] = 0', 'E0+E1', { value: { dissipator: dmax, commutator: comm }, target: 0, residual: Math.max(dmax, comm), tolerance: 1e-8, pass: dmax < 1e-8 && comm < 1e-13 });
  }
  await tick();

  // ---------------------------------------------------------------- D. clock law
  const evolve = (p, rho0, tau) => Nd.applySuper(L.cexpm(L.cscale(Nd.liouvillian(Nd.reduced(p).H, [Nd.reduced(p).R]), tau)), rho0);
  { const p = { omega: 2, c: 0, gamma: 0, clock: 'linear' }, rho0 = Nd.clockState(Math.PI / 2, 0); let worst = 0;
    for (const tau of [0.1, 1.3, 7.7, 31.4, 100]) { const r = Nd.readout(evolve(p, rho0, tau)); const O = 2 * r.coherence[0]; worst = Math.max(worst, Math.abs(O - Math.cos(2 * tau))); }
    add('clock', 'D1', 'Linear clock H = ωC/24 on E47: ⟨O⟩ = cos ωτ', 'E0+E1', { value: worst, target: 0, residual: worst, tolerance: 1e-12, pass: worst < 1e-12 }); }
  { const p = { omega: 2, c: 1, gamma: 20, clock: 'matched' }, rho0 = Nd.clockState(Math.PI / 2, 0.3), k0 = Nd.readout(rho0).coherence; let worst = 0, retMin = 1;
    for (const tau of [0.37, 3.1, 25, 200]) { const r = Nd.readout(evolve(p, rho0, tau)), ph = -2 * tau, er = k0[0] * Math.cos(ph) - k0[1] * Math.sin(ph), ei = k0[0] * Math.sin(ph) + k0[1] * Math.cos(ph); worst = Math.max(worst, Math.hypot(r.coherence[0] - er, r.coherence[1] - ei)); retMin = Math.min(retMin, r.retention); }
    add('clock', 'D2', 'Matched clock under continuous leak and recovery keeps κ(τ) = e^(−iωτ) κ(0) exactly', 'E0+E1', { value: { coherence_error: worst, min_retention: retMin }, target: 0, residual: worst, tolerance: 1e-12, pass: worst < 1e-12 && retMin < 0.999 }); }
  { const errs = [0.1, 0.01].map(c => { const p = { omega: 2, c, gamma: 20, clock: 'linear' }, s = Nd.coherenceModeRates(p).slowest, [da, db] = Nd.detunings(2, 'linear'), ad = Nd.adiabaticCoherenceRate(c, 20, da, db, 2); return { c, numeric: -s[0], adiabatic: -ad[0], rel: Math.abs(s[0] / ad[0] - 1) }; });
    const ratio = errs[0].rel / errs[1].rel;
    add('clock', 'D3', 'Linear clock is dephased by recovery: γφ = ½γc²|z_a − z_b|² (second order, error ∝ c²)', 'E0+E1', { value: errs, target: 'relative error → 0 as c², ratio 100', residual: errs[1].rel, tolerance: 1e-5, pass: errs.every(e => e.numeric > 0) && errs[1].rel < 1e-5 && ratio > 90 && ratio < 110 }); }
  { const s = Nd.coherenceModeRates({ omega: 2, c: 1, gamma: 20, clock: 'matched' }).slowest;
    add('clock', 'D4', 'Matched clock: the coherence mode is exactly −iω (zero dephasing)', 'E0+E1', { value: s, target: [0, -2], residual: Math.max(Math.abs(s[0]), Math.abs(s[1] + 2)), tolerance: 1e-12, pass: Math.abs(s[0]) < 1e-12 && Math.abs(s[1] + 2) < 1e-12 }); }
  { const lam = Object.keys(SPECTRUM).map(Number), gaps = []; for (const x of lam) for (const y of lam) if (y - x === 24) gaps.push([x, y]);
    const syn = [0, 2, 12, 20, 42], blind = []; for (const s1 of syn) for (const s2 of syn) if (s1 !== s2 && s2 - s1 === 24) blind.push([s1, s2]);
    add('clock', 'D5', 'Only (6, 30) has Casimir gap 24, so no syndrome pair makes the linear clock logical-blind', 'E0', { value: { gap24: gaps, blind_syndromes: blind }, target: { gap24: [[6, 30]], blind_syndromes: [] }, pass: gaps.length === 1 && gaps[0][0] === 6 && gaps[0][1] === 30 && blind.length === 0 }); }
  if (golden) { let worst = 0; for (const cs of golden.node) { const r = evolve({ omega: cs.omega, c: cs.c, gamma: cs.gamma, clock: cs.clock }, Nd.clockState(cs.theta, cs.phase), cs.tau); for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) worst = Math.max(worst, Math.hypot(r.re[i * 4 + j] - cs.rho[i][j][0], r.im[i * 4 + j] - cs.rho[i][j][1])); }
    add('clock', 'D6', `Parity with the Python reference: node density matrices (${golden.node.length} cases, τ ≤ 10)`, 'E1', { value: worst, target: 0, residual: worst, tolerance: 1e-12, pass: worst < 1e-12 }); }
  await tick();

  // ---------------------------------------------------------------- E. recovery statistics
  { let worst = 0; for (const c of [0.5, 1, 2]) for (const g of [5, 20, 40]) for (const d of [0, 0.5, 1]) { const H = L.cfromRows([[0, c], [c, d]]), R = L.cfromRows([[0, Math.sqrt(g)], [0, 0]]), st = Nd.countingStats(H, [R], 0); const q = Nd.stationaryLeak(c, g, d); worst = Math.max(worst, Math.abs(st.steady.re[3] - q) / q, Math.abs(st.rate - g * q) / (g * q)); }
    add('recovery', 'E1', 'Stationary leak 4c²/(γ² + 8c² + 4Δ²) and click rate γq (27 parameter points)', 'E0+E1', { value: worst, target: 0, residual: worst, tolerance: 1e-12, pass: worst < 1e-12 }); }
  { const row = g => { const H = L.cfromRows([[0, 1], [1, 0]]), R = L.cfromRows([[0, Math.sqrt(g)], [0, 0]]), st = Nd.countingStats(H, [R], 0); return { retention: 1 - st.steady.re[3], clicks: st.rate }; };
    const r20 = row(20), r40 = row(40), ex = [[101 / 102, 10 / 51], [401 / 402, 20 / 201]];
    const res = Math.max(Math.abs(r20.retention - ex[0][0]), Math.abs(r20.clicks - ex[0][1]), Math.abs(r40.retention - ex[1][0]), Math.abs(r40.clicks - ex[1][1]));
    add('recovery', 'E2', 'Review table: γ = 20 → 99.019608 %, 0.19607843; γ = 40 → 99.751244 %, 0.09950249 (exactly 101/102, 10/51, 401/402, 20/201)', 'E0+E1', { value: { gamma20: r20, gamma40: r40 }, target: { gamma20: '101/102, 10/51', gamma40: '401/402, 20/201' }, residual: res, tolerance: 1e-14, pass: res < 1e-14 && (100 * r20.retention).toFixed(6) === '99.019608' && r20.clicks.toFixed(8) === '0.19607843' && (100 * r40.retention).toFixed(6) === '99.751244' && r40.clicks.toFixed(8) === '0.09950249' }); }
  { let worst = 0, subP = true; for (const c of [0.5, 1, 2]) for (const g of [5, 20, 40]) { const H = L.cfromRows([[0, c], [c, 0]]), R = L.cfromRows([[0, Math.sqrt(g)], [0, 0]]), st = Nd.countingStats(H, [R], 0), F = Nd.resonantFano(c, g); worst = Math.max(worst, Math.abs(st.fano - F)); subP = subP && st.fano < 1; }
    add('recovery', 'E3', 'Recovery clicks are sub-Poissonian: F = 1 − 24c²γ²/(γ² + 8c²)² (full counting statistics)', 'E0+E1', { value: worst, target: 0, residual: worst, tolerance: 1e-12, pass: worst < 1e-12 && subP }); }
  { let worst = 0; const H = L.cfromRows([[0, 1], [1, 0.5]]), R = L.cfromRows([[0, Math.sqrt(20)], [0, 0]]), F1 = Nd.countingStats(H, [R], 0, 1).fano; for (const eta of [0.1, 0.5, 0.9]) worst = Math.max(worst, Math.abs(Nd.countingStats(H, [R], 0, eta).fano - Nd.fanoAtEfficiency(F1, eta)));
    add('recovery', 'E4', 'Detection efficiency thins the clicks: F_η = 1 + η(F − 1)', 'E0+E1', { value: worst, target: 0, residual: worst, tolerance: 1e-12, pass: worst < 1e-12 }); }
  if (golden) { let worst = 0; for (const f of golden.fcs) { const H = L.cfromRows([[0, f.c], [f.c, f.delta]]), R = L.cfromRows([[0, Math.sqrt(f.gamma)], [0, 0]]), st = Nd.countingStats(H, [R], 0); worst = Math.max(worst, Math.abs(st.rate - f.rate_fd), Math.abs(st.fano - f.fano_fd)); }
    add('recovery', 'E5', `Parity with the Python reference: counting statistics (${golden.fcs.length} cases, finite-difference reference)`, 'E1', { value: worst, target: 0, residual: worst, tolerance: 1e-8, pass: worst < 1e-8 }); }
  { // Monte Carlo of click counts, branch at Δ = 0: mean and variance vs K1 τ, K2 τ
    const c = 1, g = 5, tau = 80, Ntr = 800, H = L.cfromRows([[0, c], [c, 0]]), R = L.cfromRows([[0, Math.sqrt(g)], [0, 0]]), st = Nd.countingStats(H, [R], 0);
    // exact finite-time mean from |0>: E N(τ) = γ ∫0^τ <1|e^{Lt}ρ0|1> dt (Simpson on an exact propagator grid)
    const Lq = Nd.liouvillian(H, [R]), Mg = 8000, hq = tau / Mg, Pq = L.cexpm(L.cscale(Lq, hq)); let rq = Nd.vecOf(L.cfromRows([[1, 0], [0, 0]])), integ = 0;
    for (let k = 0; k <= Mg; k++) { const w = k === 0 || k === Mg ? 1 : k % 2 ? 4 : 2; integ += w * rq.re[3]; rq = L.cmul(Pq, rq); }
    const meanExact = g * integ * hq / 3;
    const Heff = L.cfromRows([[0, c], [c, [0, -g / 2]]]), dlt = 0.002, E = L.cexpm(L.cscale(Heff, 0, -dlt)), rng = xorshift(2026);
    const counts = [];
    for (let k = 0; k < Ntr; k++) {
      let pr = [1, 0], pi = [0, 0], t = 0, r = rng(), n = 0;
      while (t < tau) {
        const nr = [E.re[0] * pr[0] - E.im[0] * pi[0] + E.re[1] * pr[1] - E.im[1] * pi[1], E.re[2] * pr[0] - E.im[2] * pi[0] + E.re[3] * pr[1] - E.im[3] * pi[1]];
        const ni = [E.re[0] * pi[0] + E.im[0] * pr[0] + E.re[1] * pi[1] + E.im[1] * pr[1], E.re[2] * pi[0] + E.im[2] * pr[0] + E.re[3] * pi[1] + E.im[3] * pr[1]];
        const nn = nr[0] ** 2 + ni[0] ** 2 + nr[1] ** 2 + ni[1] ** 2;
        if (nn < r) { n++; pr = [1, 0]; pi = [0, 0]; r = rng(); } else { pr = nr; pi = ni; }
        t += dlt;
      }
      counts.push(n);
    }
    const mean = counts.reduce((s, x) => s + x, 0) / Ntr, vr = counts.reduce((s, x) => s + (x - mean) ** 2, 0) / (Ntr - 1);
    const zMean = (mean - meanExact) / Math.sqrt(st.variance * tau / Ntr), zVar = (vr - st.variance * tau) / (st.variance * tau * Math.sqrt(2 / (Ntr - 1)));
    add('recovery', 'E6', 'Monte Carlo jump trajectories reproduce the click mean K1τ and variance K2τ (800 runs to τ = 80)', 'E2', { value: { mean, variance: vr, fano_sample: vr / mean, z_mean: zMean, z_variance: zVar }, target: { mean: meanExact, mean_asymptotic: st.rate * tau, variance: st.variance * tau, fano: st.fano }, residual: Math.max(Math.abs(zMean), Math.abs(zVar)), tolerance: 'z < 4', pass: Math.abs(zMean) < 4 && Math.abs(zVar) < 4 });
  }
  await tick();

  // ---------------------------------------------------------------- F. time-chain identity
  { const p = { omega: 2, c: 1, gamma: 20, clock: 'linear' }, { H, R } = Nd.reduced(p), Lv = Nd.liouvillian(H, [R]), rho0 = Nd.clockState(1.1, 0.4), Nl = Math.sqrt(0.74), t = 3.7;
    const a = Nd.applySuper(L.cexpm(L.cscale(Lv, Nl * t)), rho0), b = Nd.applySuper(L.cexpm(L.cscale(Nd.liouvillian(L.cscale(H, Nl), [L.cscale(R, Math.sqrt(Nl))]), t)), rho0);
    const res = L.cmaxabs(L.csub(a, b));
    add('time', 'F1', 'Time-chain identity: generator (N H, √N R) over t equals (H, R) over τ = N t', 'E0+E1', { value: res, target: 0, residual: res, tolerance: 1e-13, pass: res < 1e-13 }); }
  { const p = { omega: 2, c: 1, gamma: 20, clock: 'linear' }, { H, R } = Nd.reduced(p), Lv = Nd.liouvillian(H, [R]), rho0 = Nd.clockState(1.1, 0.4);
    const Nf = t => Math.sqrt(1 - 2 * (0.3 * Math.sin(t)) ** 2), T = 4, h = 2e-4;
    let x = Nd.vecOf(rho0); const f = (t, y) => L.cscale(L.cmul(Lv, y), Nf(t));
    for (let t = 0; t < T - 1e-12; t += h) { const k1 = f(t, x), k2 = f(t + h / 2, L.caxpy(h / 2, 0, k1, L.cclone(x))), k3 = f(t + h / 2, L.caxpy(h / 2, 0, k2, L.cclone(x))), k4 = f(t + h, L.caxpy(h, 0, k3, L.cclone(x))); x = L.caxpy(h / 6, 0, k1, x); L.caxpy(h / 3, 0, k2, x); L.caxpy(h / 3, 0, k3, x); L.caxpy(h / 6, 0, k4, x); }
    let tau = 0; const M = 200000; for (let k = 0; k <= M; k++) { const w = k === 0 || k === M ? 1 : k % 2 ? 4 : 2; tau += w * Nf(T * k / M); } tau *= T / M / 3;
    const ex = Nd.vecOf(Nd.applySuper(L.cexpm(L.cscale(Lv, tau)), rho0)), res = L.cmaxabs(L.csub(x, ex));
    add('time', 'F2', 'Time-dependent lapse: dρ/dt = N(t)L_τρ integrates to e^(τ(t) L_τ) ρ(0)', 'E0+E1', { value: { tau, residual: res }, target: 0, residual: res, tolerance: 1e-9, pass: res < 1e-9 });
    const tt = 2.2, hh = 1e-4, rhoAt = s => { let tq = 0; const Mq = 20000; for (let k = 0; k <= Mq; k++) { const w = k === 0 || k === Mq ? 1 : k % 2 ? 4 : 2; tq += w * Nf(s * k / Mq); } tq *= s / Mq / 3; return { tq, v: Nd.vecOf(Nd.applySuper(L.cexpm(L.cscale(Lv, tq)), rho0)) }; };
    const plus = rhoAt(tt + hh), minus = rhoAt(tt - hh), mid = rhoAt(tt), dt = L.cscale(L.csub(plus.v, minus.v), 1 / (2 * hh));
    const dtau = L.cmul(Lv, mid.v), chain = L.cmaxabs(L.csub(L.cscale(dt, 1 / Nf(tt)), dtau));
    add('time', 'F3', 'C0-TIME-CHAIN: dτ = N dt implies d/dτ = N⁻¹ d/dt (not N d/dt)', 'E0+E1', { value: chain, target: 0, residual: chain, tolerance: 1e-6, pass: chain < 1e-6 && L.cmaxabs(L.csub(L.cscale(dt, Nf(tt)), dtau)) > 1e-3 }); }
  { const Nl = Math.sqrt(0.82), t = 1.9, mk = lapses => new Network({ n: 2, edges: [[0, 1]], omega: 2, c: 0.5, gamma: 10, J: 0.3, Jzz: 0.2, clock: 'linear', lapses });
    const psi = Network.productState([Network.clockVector(Math.PI / 2, 0), Network.clockVector(Math.PI / 3, 0.7)]), r0 = Network.density(psi);
    const a = mk([Nl, Nl]).propagate(r0, t), b = mk([1, 1]).propagate(r0, Nl * t);
    const res = Math.max(...a.re.map((x, q) => Math.hypot(x - b.re[q], a.im[q] - b.im[q])));
    add('time', 'F4', 'Network time-chain identity: equal lapses N give dρ/dt = N L_τ ρ', 'E0+E1', { value: res, target: 0, residual: res, tolerance: 1e-12, pass: res < 1e-12 }); }
  { const N1 = 1, N2 = Math.sqrt(0.74), t = 5.3, net = new Network({ n: 2, edges: [[0, 1]], omega: 2, c: 1, gamma: 20, J: 0, Jzz: 0, clock: 'matched', lapses: [N1, N2] });
    const r = net.readout(net.propagate(Network.density(Network.productState([Network.clockVector(), Network.clockVector()])), t));
    const ph = r.map(x => -Math.atan2(x.coherence[1], x.coherence[0])), diff = Math.atan2(Math.sin(ph[0] - ph[1]), Math.cos(ph[0] - ph[1])), ex = Math.atan2(Math.sin(2 * (N1 - N2) * t), Math.cos(2 * (N1 - N2) * t));
    add('time', 'F5', 'Two matched clocks on different worldlines drift apart at exactly (N₁ − N₂)ω', 'E0+E1', { value: diff, target: ex, residual: Math.abs(diff - ex), tolerance: 1e-12, pass: Math.abs(diff - ex) < 1e-12 }); }
  if (golden) { const cfg = golden.network; let worst = 0;
    for (const cs of cfg.cases) { const net = new Network({ n: 2, edges: cfg.edges, omega: cfg.omega, c: cfg.c, gamma: cfg.gamma, J: cfg.J, Jzz: cfg.Jzz, clock: cs.clock, lapses: cfg.velocities.map(v => lapse(v)) });
      const ro = net.readout(net.propagate(Network.density(Network.productState(cfg.initial.map(([th, ph]) => Network.clockVector(th, ph)))), cs.t));
      ro.forEach((r, v) => { worst = Math.max(worst, Math.abs(r.retention - cs.retention[v]), Math.hypot(r.coherence[0] - cs.coherence[v][0], r.coherence[1] - cs.coherence[v][1]), Math.abs(r.clicks - cs.clicks[v])); }); }
    add('time', 'F6', `Parity with the Python reference: two-node TQGNN readouts (${cfg.cases.length} cases)`, 'E1', { value: worst, target: 0, residual: worst, tolerance: 1e-12, pass: worst < 1e-12 }); }
  await tick();

  // ---------------------------------------------------------------- G. prime and synchronisation clocks
  const pc = K.primeClock(car);
  { const orth = L.rmaxabs(L.radd(L.rmul(L.rtrans(pc.U), pc.U), L.reye(125), -1)), borth = L.rmaxabs(L.radd(L.rmul(L.rtrans(pc.B), pc.B), L.reye(47), -1));
    add('prime', 'G1', 'The 47-cycle U47 on an orthonormal basis of E47 is orthogonal', 'E1', { value: { orthogonality: orth, basis: borth }, target: 0, residual: Math.max(orth, borth), tolerance: 1e-13, pass: orth < 1e-13 && borth < 1e-13 }); }
  { const o = K.orderOf(pc.U, 60), minBelow = Math.min(...o.dist.slice(0, 46));
    add('prime', 'G2', 'ord(U47) = 47 (no smaller power is the identity; 47 is prime)', 'E0+E1', { value: o.order, target: 47, residual: o.dist[46], tolerance: 1e-12, pass: o.order === 47 && minBelow > 0.5 }); }
  { const sr = K.primeSpectrumResidual(pc);
    add('prime', 'G3', 'spec(U47) on E47 is every 47th root of unity', 'E0+E1', { value: sr, target: 0, residual: sr, tolerance: 1e-13, pass: sr < 1e-13 }); }
  { const perm = [...Array(125).keys()].map(n => { const a = Math.floor(n / 25), b = Math.floor(n / 5) % 5, c = n % 5; return ((a + 1) % 5) * 25 + ((b + 1) % 5) * 5 + (c + 1) % 5; });
    const o = K.permutationOrder(perm); let p47 = [...Array(125).keys()], p2 = [...Array(125).keys()]; for (let k = 0; k < 47; k++) p47 = p47.map(i => perm[i]); for (let k = 0; k < 2; k++) p2 = p2.map(i => perm[i]);
    const same = p47.every((v, i) => v === p2[i]), notId = p47.some((v, i) => v !== i);
    add('prime', 'G4', 'Torus clock U_Δ = U_xU_yU_z on 5×5×5 has order 5, so U_Δ⁴⁷ = U_Δ² ≠ I', 'E0+E1', { value: { order: o, U47_equals_U2: same }, target: { order: 5 }, pass: o === 5 && same && notId });
    add('prime', 'G5', 'Joint clock U_Δ ⊗ U47 has period lcm(5, 47) = 235', 'E0', { value: K.lcm(o, 47), target: 235, pass: K.lcm(o, 47) === 235 }); }
  { const z = K.modsyncOrbit(1200), ords = [K.multOrder(2, 47), K.multOrder(2, 125), K.multOrder(2, 23), K.multOrder(2, 25)];
    let big = true; for (let n = 0; n <= 300; n++) if (Number(K.powmodBig(2, 2n ** BigInt(n), 5875n)) !== z[n]) { big = false; break; }
    const tp = K.modsyncTransientPeriod(z), p47 = K.periodMod(z, 47), p125 = K.periodMod(z, 125), distinct = new Set(z.slice(2, 222)).size;
    add('prime', 'G6', 'Synchronisation clock z ↦ z² mod 5875: transient 2, period 220 = lcm(11, 20)', 'E0', { value: { orders: ords, transient: tp.mu, period: tp.period, mod47: p47, mod125: p125, distinct, bigint: big }, target: { orders: [23, 100, 11, 20], transient: 2, period: 220, mod47: 11, mod125: 20, distinct: 220 }, pass: ords.join() === '23,100,11,20' && tp.mu === 2 && tp.period === 220 && p47 === 11 && p125 === 20 && distinct === 220 && big }); }
  await tick();
  { const j = K.modsyncJoint(car), r220 = (15 / 17) ** 220, r110 = (15 / 17) ** 110;
    add('prime', 'G7', 'Joint operator U = Γ ⊗ T stepped with the clock: U²²⁰ = Γ²²⁰ ⊗ I, within (15/17)²²⁰ of P ⊗ I', 'E0+E1', { value: j, target: { U220_minus_P: `≤ ${r220.toExponential(4)}` }, residual: j.U220_minus_Gamma220, tolerance: 1e-13,
      pass: j.U220_minus_Gamma220 < 1e-13 && j.U220_minus_P <= r220 + 1e-14 && j.U110_minus_PT110 <= r110 + 1e-14 && j.U110_minus_P > 0.1 && Math.abs(j.success220 - j.protectedWeight0) < 1e-12 && Math.abs(j.returnFidelity220 - 1) < 1e-12 }); }
  await tick();

  // ---------------------------------------------------------------- H. the Matrix bridge
  { const d = await Rc.sha256hex('abc');
    add('matrix', 'H1', 'SHA-256 test vector (FIPS 180-2, "abc")', 'E1', { value: d, target: 'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad', pass: d === 'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad' }); }
  { const a = { b: 1, a: [2, { y: 1, x: 0 }] }, b = { a: [2, { x: 0, y: 1 }], b: 1 };
    add('matrix', 'H2', 'sorted-json-v1 canonicalisation is independent of key order', 'E1', { value: Rc.canonical(a), target: Rc.canonical(b), pass: Rc.canonical(a) === Rc.canonical(b) && (await Rc.digestOf(a)) === (await Rc.digestOf(b)) }); }
  { const s = [[0.5, 0.1], [0.3, -0.2], [0.1, 0.4], [-0.2, 0.2], [0.3, 0.1], [0.1, -0.3], [0.2, 0.2], [-0.1, 0.25]], nn = Math.sqrt(s.reduce((a, z) => a + z[0] ** 2 + z[1] ** 2, 0)), v = s.map(z => [z[0] / nn, z[1] / nn]);
    const run = gs => gs.reduce((acc, [g, t]) => Rc.matrixApply(acc.map(([r, i]) => ({ r, i })), g, t).map(z => [z.r, z.i]), v.map(z => z.slice()));
    const d = (x, y) => Math.sqrt(x.reduce((a, z, k) => a + (z[0] - y[k][0]) ** 2 + (z[1] - y[k][1]) ** 2, 0));
    const T8 = d(run(Array(8).fill(['T', [1]])), v), T2 = d(run([['T', [1]], ['T', [1]]]), run([['S', [1]]])), T4 = d(run(Array(4).fill(['T', [1]])), run([['Z', [1]]])), H2 = d(run([['H', [0]], ['H', [0]]]), v), CX2 = d(run([['CX', [0, 2]], ['CX', [0, 2]]]), v), SW2 = d(run([['SW', [1, 2]], ['SW', [1, 2]]]), v);
    const worst = Math.max(T8, T2, T4, H2, CX2, SW2);
    add('matrix', 'H3', 'Matrix gate arithmetic: T⁸ = I, T² = S, T⁴ = Z, H² = CX² = SW² = I', 'E0+E1', { value: { T8, T2, T4, H2, CX2, SW2 }, target: 0, residual: worst, tolerance: 1e-15, pass: worst < 1e-15 }); }
  { // one forward clock tick (ωτ = π/4) on an encoded Matrix qubit equals the T gate, under leak and recovery
    const omega = 2 * Math.PI / 8, p = { omega, c: 1, gamma: 20, clock: 'matched' }, { H, R } = Nd.reduced(p), Lv = Nd.liouvillian(H, [R]);
    const ghz = [[Math.SQRT1_2, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [Math.SQRT1_2, 0]];
    let worst = 0, purityMin = 1;
    for (const k of [0, 1, 2]) for (const ticks of [1, 3, 8]) {
      const rho16 = Rc.densityOf(Rc.encodeQubit(ghz, k)), S1 = L.cexpm(L.cscale(Lv, (Math.PI / 4) * ticks / omega));
      const out = L.cmat(16); // (S ⊗ id) on 4x4 blocks
      for (let s = 0; s < 4; s++) for (let s2 = 0; s2 < 4; s2++) { const X = L.cmat(4); for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) { X.re[i * 4 + j] = rho16.re[(i * 4 + s) * 16 + j * 4 + s2]; X.im[i * 4 + j] = rho16.im[(i * 4 + s) * 16 + j * 4 + s2]; } const Y = Nd.applySuper(S1, X); for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) { out.re[(i * 4 + s) * 16 + j * 4 + s2] = Y.re[i * 4 + j]; out.im[(i * 4 + s) * 16 + j * 4 + s2] = Y.im[i * 4 + j]; } }
      const Lr = Rc.logicalDensity(out, k);
      let psiT = ghz.map(z => z.slice()); for (let n = 0; n < ticks; n++) psiT = Rc.matrixApply(psiT.map(([r, i]) => ({ r, i })), 'T', [k]).map(z => [z.r, z.i]);
      worst = Math.max(worst, Math.abs(1 - Rc.expectationIn(Lr, psiT))); purityMin = Math.min(purityMin, Rc.purity(Lr));
    }
    add('matrix', 'H4', 'One clock tick ωτ = π/4 on an encoded Matrix qubit is exactly the T gate, under leak and recovery', 'E0+E1', { value: { infidelity: worst, purity_min: purityMin }, target: 0, residual: worst, tolerance: 1e-12, pass: worst < 1e-12 && Math.abs(1 - purityMin) < 1e-12 }); }
  { // receipt verification statuses on synthetic receipts
    const mk = async (circuit, statevector, measurement = { collapsed: false, outcome: null }) => Rc.seal({ schema: Rc.MATRIX_SCHEMA, source: 'matrix-quantum-simulator', exported_at: '2026-10-19T00:00:00.000Z', qubits: 3, basis: 'q0q1q2', initial_state: '000', circuit, statevector, measurement, norm: statevector.reduce((s, z) => s + z[0] ** 2 + z[1] ** 2, 0) });
    const bell = [{ gate: 'H', targets: [0], column: 0 }, { gate: 'CX', targets: [0, 1], column: 1 }];
    const ok = await Rc.verifyMatrixReceipt(await mk(bell, Rc.matrixReplay(bell, 'column')));
    const swapped = [{ gate: 'S', targets: [0], column: 1 }, { gate: 'H', targets: [0], column: 0 }];
    const ord = await Rc.verifyMatrixReceipt(await mk(swapped, Rc.matrixReplay(swapped, 'insertion')));
    const stale = await Rc.verifyMatrixReceipt(await mk(bell, Rc.matrixReplay([bell[0]], 'column')));
    const good = await mk(bell, Rc.matrixReplay(bell, 'column')), flipped = good.statevector.map((z, k) => k === 6 ? [-z[0], -z[1]] : z.slice());
    const tampered = await Rc.verifyMatrixReceipt({ ...good, statevector: flipped });            // digest kept
    const resealed = await Rc.verifyMatrixReceipt(await mk(bell, flipped));                     // digest recomputed
    const meas = await Rc.verifyMatrixReceipt(await mk(bell, [[0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [1, 0], [0, 0]], { collapsed: true, outcome: 6 }));
    const st = { ok: ok.status, order: ord.status, stale: stale.status, tampered: tampered.status, resealed: resealed.status, measured: meas.status };
    const want = { ok: 'REPLAY-VERIFIED', order: 'ORDER-DEPENDENT', stale: 'REPLAY-MISMATCH', tampered: 'DIGEST-MISMATCH', resealed: 'REPLAY-MISMATCH', measured: 'MEASURED-CONSISTENT' };
    add('matrix', 'H5', 'Receipt replay tells apart verified, order-dependent, stale, tampered, re-sealed and measured receipts', 'E1', { value: st, target: want, pass: Object.keys(want).every(k => st[k] === want[k]) });
    const ghzC = [{ gate: 'H', targets: [0], column: 0 }, { gate: 'CX', targets: [0, 1], column: 1 }, { gate: 'CX', targets: [1, 2], column: 2 }];
    const src = await mk(ghzC, Rc.matrixReplay(ghzC, 'column')), { receipt, word } = await Rc.matrixReceiptWithTicks(src, 2, 7);
    Rc.checkMatrixShape(receipt); const back = await Rc.verifyMatrixReceipt(receipt); let appOk = false; try { appOk = Rc.matrixAppValid(receipt); } catch { appOk = false; }
    add('matrix', 'H6', 'A clock-tick receipt exported back to the Matrix passes the deployed importer and replays exactly', 'E1', { value: { status: back.status, word: word.join(' '), residual: back.residualColumnOrder, deployed_importer: appOk }, target: { status: 'REPLAY-VERIFIED', word: 'Z S T', deployed_importer: true }, residual: back.residualColumnOrder, tolerance: 0, pass: back.status === 'REPLAY-VERIFIED' && back.digestOk && word.join(' ') === 'Z S T' && back.residualColumnOrder === 0 && appOk === true }); }

  const passed = checks.filter(c => c.pass).length;
  const payload = {
    certificate_code: CERTIFICATE_ID, title: 'Quantum Clock — E47 clock-memory unit, TQGNN network and Matrix bridge (executed)', engine: ENGINE_VERSION,
    status: `${passed}/${checks.length} PASS`, passed, total: checks.length,
    evidence_classes: { E0: 'exact proof', E1: 'executable reconstruction', E2: 'simulation' },
    claim_boundary: 'The E47 Quantum Clock: a finite-dimensional quantum clock with phase-based elapsed-time readout, a 4-dim clock-memory node, matched GKSL dynamics, a two- or three-node TQGNN, and the 3-qubit Matrix receipt protocol. Evidence comprises mathematical derivation and numerical quantum simulation. Laboratory hardware realization and calibration to physical seconds remain open; the metric-to-phase bridge (ω·Δτ per tick) specifies the readout convention.',
    golden: golden ? { schema: golden.schema, generator_sha256: golden.generator_sha256, environment: golden.environment } : null,
    environment, elapsed_ms: Math.round((globalThis.performance || Date).now() - t0), checks,
  };
  return { ...payload, digest: await Rc.digestOf({ ...payload, environment: null, elapsed_ms: null }) };
}
