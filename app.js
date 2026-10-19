// The Quantum Clock page.
import * as Q from './engine/index.js';
import { palette, fitCanvas, drawFace, drawStrip, plot, bars, ringFigure } from './ui.js';

const { NodeClock, Network, NetworkTrajectory, xorshift, receipts: R, clocks: K } = Q;
const $ = id => document.getElementById(id);
const reduced = matchMedia('(prefers-reduced-motion: reduce)').matches;
let P = palette();
matchMedia('(prefers-color-scheme: light)').addEventListener?.('change', () => { P = palette(); });
const fmt = (x, d = 4) => Number.isFinite(x) ? x.toFixed(d) : '—';
const sci = (x, d = 2) => Number.isFinite(x) ? (x === 0 ? '0' : x.toExponential(d)) : '—';
const download = (name, text, type = 'application/json') => { const a = document.createElement('a'), u = URL.createObjectURL(new Blob([text], { type })); a.href = u; a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(u), 2000); };
const say = (id, text, cls = '') => { const e = $(id); e.innerHTML = ''; const s = document.createElement('span'); if (cls) s.className = cls; s.textContent = text; e.append(s); };

// ======================================================================= the clock
const init = { ...Q.DEFAULTS, speed: 1, ...Q.decodeHash(location.hash) };
let speed = init.speed; delete init.speed;
let clock;
try { clock = new NodeClock(init); } catch { clock = new NodeClock(Q.DEFAULTS); }
let running = !reduced, flash = 0, hist = [], lastHistTau = -1;

const sliders = { omega: ['in-omega', 'o-omega', 3], c: ['in-c', 'o-c', 2], gamma: ['in-gamma', 'o-gamma', 1], eta: ['in-eta', 'o-eta', 2] };
const vsl = [['in-v1', 'o-v1'], ['in-v2', 'o-v2'], ['in-v3', 'o-v3']];
function syncControls() {
  const p = clock.p;
  for (const [k, [inp, out, d]] of Object.entries(sliders)) { $(inp).value = p[k]; $(out).textContent = (+p[k]).toFixed(d); }
  vsl.forEach(([inp, out], k) => { $(inp).value = p.velocity[k]; $(out).textContent = (+p.velocity[k]).toFixed(3); });
  $('in-speed').value = speed; $('o-speed').textContent = speed.toFixed(1) + ' t/s';
  document.querySelectorAll('[data-clock]').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.clock === p.clock)));
  document.querySelectorAll('[data-metric]').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.metric === p.metric)));
  explain();
}
function reconfigure(change, { reset = false } = {}) {
  try { clock.configure(change); if (reset) { clock.reset(); hist = []; lastHistTau = -1; } say('clock-msg', ''); if ('omega' in change && netObj) buildNetwork(); }
  catch (e) { say('clock-msg', e.message, 'bad'); }
  syncControls(); updateBridgeClock();
}
for (const [k, [inp]] of Object.entries(sliders)) $(inp).addEventListener('input', e => reconfigure({ [k]: +e.target.value }));
const weights = () => clock.p.metric === 'minkowski' ? [1, 1, 1] : [2, 8, 2];
vsl.forEach(([inp], k) => $(inp).addEventListener('input', e => {
  const v = clock.p.velocity.slice(), w = weights(); v[k] = +e.target.value;
  const others = w.reduce((s, wi, i) => i === k ? s : s + wi * v[i] * v[i], 0), vmax = Math.sqrt(Math.max(0, (0.98 - others) / w[k]));
  if (v[k] > vmax) { v[k] = Math.floor(vmax * 1000) / 1000; say('clock-msg', clock.p.metric === 'minkowski' ? 'Held inside the light cone: v₁² + v₂² + v₃² < 1.' : 'Held inside the light cone: 2v₁² + 8v₂² + 2v₃² < 1 for the E47 metric.'); }
  reconfigure({ velocity: v });
}));
$('in-speed').addEventListener('input', e => { speed = +e.target.value; $('o-speed').textContent = speed.toFixed(1) + ' t/s'; });
document.querySelectorAll('[data-clock]').forEach(b => b.addEventListener('click', () => reconfigure({ clock: b.dataset.clock }, { reset: true })));
document.querySelectorAll('[data-metric]').forEach(b => b.addEventListener('click', () => {
  const v = clock.p.velocity.slice(), w = b.dataset.metric === 'minkowski' ? [1, 1, 1] : [2, 8, 2];
  const s = w.reduce((a, wi, i) => a + wi * v[i] * v[i], 0); if (s >= 0.98) { const f = Math.sqrt(0.97 / s); for (let i = 0; i < 3; i++) v[i] *= f; }
  reconfigure({ metric: b.dataset.metric, velocity: v });
}));
const PRESETS = {
  standard: { clock: 'matched', omega: 2 * Math.PI / 8, c: 1, gamma: 20, eta: 1, velocity: [0, 0, 0] },
  kicks: { clock: 'linear', omega: 2, c: 1, gamma: 4, eta: 1, velocity: [0, 0, 0] },
  zeno: { clock: 'linear', omega: 2, c: 1, gamma: 60, eta: 1, velocity: [0, 0, 0] },
  deterministic: { eta: 0 },
  moving: { velocity: [0.3, 0, 0], metric: 'e47' },
};
document.querySelectorAll('[data-preset]').forEach(b => b.addEventListener('click', () => reconfigure(PRESETS[b.dataset.preset], { reset: true })));
$('play').addEventListener('click', () => { running = !running; $('play').textContent = running ? 'Pause' : 'Play'; $('play').setAttribute('aria-pressed', String(running)); });
$('reset').addEventListener('click', () => { clock.reset(); hist = []; lastHistTau = -1; });
$('permalink').addEventListener('click', async () => {
  const h = Q.encodeHash({ ...clock.p, speed }); history.replaceState(null, '', h);
  try { await navigator.clipboard.writeText(location.href); say('clock-msg', 'Link to this clock copied.', 'ok'); } catch { say('clock-msg', 'Link is in the address bar.'); }
});
$('save-clock').addEventListener('click', async () => {
  const s = clock.state(), rho = Q.linalg.ctoRows(clock.rho);
  const payload = { schema: R.CLOCK_SCHEMA, source: 'quantum-clock', engine: Q.ENGINE_VERSION, exported_at: new Date().toISOString(),
    params: { omega: clock.p.omega, c: clock.p.c, gamma: clock.p.gamma, clock: clock.p.clock, syndrome: clock.p.syndrome, eta: clock.p.eta, velocity: clock.p.velocity, metric: clock.p.metric },
    time: { t: s.t, tau: s.tau, lapse: s.lapse }, state: { basis: 'node index 2l+f: a, a′, b, b′', dim: 4, density: rho },
    readouts: { retention: s.deterministic.retention, coherence: s.deterministic.coherence, phase_time: s.phaseTime, clicks: s.clicks } };
  download(`quantum-clock-tau-${s.tau.toFixed(2)}.json`, JSON.stringify(await R.seal(payload), null, 2));
});

function explain() {
  const p = clock.p, [da, db] = clock.detuning;
  const txt = p.clock === 'matched'
    ? `Matched clock: the syndrome sector C = 42 is retuned so both logical branches leak with the same detuning, Δ = ${fmt(da, 3)}. Recovery cannot tell |a⟩ from |b⟩, so the phase is exact on every run, through every click. The leak still costs population: ${fmt(100 * (1 - clock.leak[0]), 4)} % stays in E47.`
    : `Linear clock f = ωC/24: the branches leak with different detunings, Δa = ${fmt(da, 3)} and Δb = ${fmt(db, 3)}, because only (6, 30) has Casimir gap 24. Each click carries which-branch information, so a monitored run's phase jumps at clicks and the average dephases. The leak also pulls the clock frequency to ω + δω with δω = ${sci(-Q.coherenceModeRates(p).slowest[1] - p.omega, 3)}, so even the gold hand drifts unless that shift is calibrated out.`;
  $('clock-explain').textContent = txt + (p.eta === 0 ? ' With η = 0 nothing is recorded: the green hand is the gold one.' : '');
}

function clockReadouts() {
  const s = clock.state(), p = clock.p, d = s.deterministic;
  $('r-t').textContent = fmt(s.t, 2);
  $('r-tau').textContent = fmt(s.tau, 3);
  $('r-lapse').textContent = `N = ${fmt(s.lapse, 6)} · τ = N t`;
  $('r-phase').textContent = fmt(s.phaseTime, 4);
  $('r-phase-err').textContent = `error ${sci(s.phaseTime - s.tau)}`;
  $('r-phase-c').textContent = fmt(s.phaseTimeConditional, 4);
  $('r-phase-c-err').textContent = `error ${sci(s.phaseTimeConditional - s.tau)}`;
  $('r-clicks').textContent = String(s.clicks);
  const prec = s.clickRate > 0 && s.tau > 0 ? Math.sqrt(clock.branchStats[0].fano / (s.clickRate * s.tau)) : NaN;
  $('r-clicktime').textContent = s.clickTime !== null && s.clicks ? `click clock τ̂ = ${fmt(s.clickTime, 2)} (±${fmt(100 * prec, 1)} %)` : 'click clock τ̂ = n/λ';
  $('r-ret').textContent = fmt(d.retention, 6);
  const stat = 1 - (d.pa * clock.leak[0] + d.pb * clock.leak[1]);
  $('r-ret-closed').textContent = `stationary ${fmt(stat, 6)}`;
  $('r-coh').textContent = fmt(Math.hypot(...d.coherence), 6);
  const mode = Q.coherenceModeRates(p).slowest;
  $('r-deph').textContent = p.clock === 'matched' ? `dephasing ${sci(Math.abs(mode[0]))} (exactly 0)` : `dephasing γφ = ${sci(-mode[0], 3)}`;
  $('r-rate').textContent = fmt(s.clickRate, 5);
  $('r-fano').textContent = `Fano F_η = ${fmt(clock.branchStats[0].fano, 5)}`;
  const frac = ((s.hand / (2 * Math.PI)) % 1 + 1) % 1, prime = Math.floor(frac * 47 + 1e-9) % 47;
  $('r-radix').textContent = `${K.radixFloat(frac, 60, 3)}₆₀ · ${K.radixFloat(frac, 12, 4)}₁₂ · ${K.radixFloat(frac, 5, 5)}₅`;
  $('r-prime').textContent = `prime tick ${prime}/47 = ${K.radixString(prime, 47, 60, 4)}₆₀`;
}

// ======================================================================= the network
const net = { n: 3, topo: 'path', mode: 'flow', clock: 'matched', J: 0, Jzz: 0, c: 0.5, gamma: 10, v: [0, 0.3, Math.sqrt(0.13)], running: !reduced };
let netObj, netRho, netTraj, netHist = [], netFlash = [], netUnwrap = [];
const edgesOf = () => net.topo === 'ring' && net.n === 3 ? [[0, 1], [1, 2], [0, 2]] : [...Array(net.n - 1).keys()].map(k => [k, k + 1]);
function buildNetwork() {
  const lapses = net.v.slice(0, net.n).map(v => Q.lapse([v, 0, 0]));
  netObj = new Network({ n: net.n, edges: edgesOf(), omega: clock.p.omega, c: net.c, gamma: net.gamma, J: net.J, Jzz: net.Jzz, clock: net.clock, lapses });
  const psi = Network.productState([...Array(net.n)].map(() => Network.clockVector(Math.PI / 2, 0)));
  netRho = Network.density(psi);
  netTraj = new NetworkTrajectory(netObj, psi, Math.min(0.01, 0.2 / Math.max(1, net.gamma)), xorshift(4747));
  netHist = []; netFlash = new Array(net.n).fill(0); netUnwrap = [...Array(net.n)].map(() => ({ last: [0.5, 0], acc: 0 })); net.t = 0;
  renderNodeControls(); say('net-msg', '');
}
function renderNodeControls() {
  const host = $('node-ctl'); host.innerHTML = '';
  for (let v = 0; v < net.n; v++) {
    const box = document.createElement('div'); box.className = 'box';
    box.innerHTML = `<div class="ctl"><label for="nv-${v}">clock ${v + 1} · velocity v₁ <output id="onv-${v}">${net.v[v].toFixed(3)}</output></label><input id="nv-${v}" type="range" min="0" max="0.7" step="0.001" value="${net.v[v]}"></div><span class="mono muted" id="nl-${v}"></span>`;
    host.append(box);
    box.querySelector('input').addEventListener('input', e => { net.v[v] = Math.min(+e.target.value, 0.7); buildNetwork(); });
    $(`nl-${v}`).textContent = `N = ${fmt(Q.lapse([net.v[v], 0, 0]), 6)}`;
  }
  const ro = $('net-readout'); ro.innerHTML = '';
  for (let v = 0; v < net.n; v++) ro.insertAdjacentHTML('beforeend', `<div class="cell"><span class="k">clock ${v + 1} reads</span><span class="val" id="nr-${v}">0</span><span class="sub" id="ns-${v}">τ = 0</span></div>`);
}
const bindNet = (id, out, key, d) => $(id).addEventListener('input', e => { net[key] = +e.target.value; $(out).textContent = net[key].toFixed(d); buildNetwork(); });
bindNet('n-J', 'o-nJ', 'J', 2); bindNet('n-Jzz', 'o-nJzz', 'Jzz', 2); bindNet('n-c', 'o-nc', 'c', 2); bindNet('n-gamma', 'o-ngamma', 'gamma', 1);
['o-nJ', 'o-nJzz', 'o-nc', 'o-ngamma'].forEach((o, k) => { $(o).textContent = [net.J, net.Jzz, net.c, net.gamma][k].toFixed(k === 3 ? 1 : 2); });
const toggle = (attr, key, cast = x => x) => document.querySelectorAll(`[data-${attr}]`).forEach(b => b.addEventListener('click', () => { net[key] = cast(b.dataset[attr]); document.querySelectorAll(`[data-${attr}]`).forEach(x => x.setAttribute('aria-pressed', String(x === b))); buildNetwork(); }));
toggle('nodes', 'n', Number); toggle('topo', 'topo'); toggle('netmode', 'mode'); toggle('netclock', 'clock');
$('net-play').addEventListener('click', () => { net.running = !net.running; $('net-play').textContent = net.running ? 'Pause' : 'Play'; $('net-play').setAttribute('aria-pressed', String(net.running)); });
$('net-reset').addEventListener('click', buildNetwork);
$('net-tqgnn').addEventListener('click', () => { net.n = 3; net.v = [0, 0.3, Math.sqrt(0.13)]; document.querySelectorAll('[data-nodes]').forEach(x => x.setAttribute('aria-pressed', String(x.dataset.nodes === '3'))); buildNetwork(); say('net-msg', 'Velocities 0, 0.3 and √0.13 along v₁ give the TQGNN lapses 1, √0.82, √0.74 in the E47 metric.'); });
$('net-learn').addEventListener('click', () => {
  if (netHist.length < 20) { say('net-msg', 'Let the network run a little first.'); return; }
  const ts = netHist.map(h => h.t), est = [];
  for (let v = 0; v < net.n; v++) { const ys = netHist.map(h => h.phase[v]), n = ts.length, mt = ts.reduce((s, x) => s + x, 0) / n, my = ys.reduce((s, x) => s + x, 0) / n; let sxy = 0, sxx = 0; for (let k = 0; k < n; k++) { sxy += (ts[k] - mt) * (ys[k] - my); sxx += (ts[k] - mt) ** 2; } est.push(sxy / sxx / clock.p.omega); }
  const truth = netObj.lapses, err = Math.max(...est.map((e, k) => Math.abs(e - truth[k])));
  say('net-msg', `Recovered from the clock phases: ${est.map(e => e.toFixed(9)).join(', ')} · metric lapses ${truth.map(e => e.toFixed(9)).join(', ')} · max error ${sci(err)}${net.J || net.Jzz ? ' (with coupling on, the phases also carry the exchange, so the slopes are only approximately the lapses)' : ''}${net.clock === 'linear' && net.c > 0 ? ' (the linear clock’s frequency is shifted by the leak, so its slope is N(ω + δω))' : ''}.`, err < 1e-6 ? 'ok' : '');
});

function netAdvance(dt) {
  if (net.mode === 'flow') netRho = netObj.propagate(netRho, dt);
  else { const fresh = netTraj.advance(dt); for (const c of fresh) netFlash[c.node] = 1; }
  net.t += dt;
  const ro = netObj.readout(net.mode === 'flow' ? netRho : netTraj.density());
  ro.forEach((r, v) => { const u = netUnwrap[v], cur = r.coherence, d = Math.atan2(u.last[0] * cur[1] - u.last[1] * cur[0], u.last[0] * cur[0] + u.last[1] * cur[1]); u.acc -= d; u.last = cur; });
  if (!netHist.length || net.t - netHist[netHist.length - 1].t > 0.05) netHist.push({ t: net.t, phase: netUnwrap.map(u => u.acc) });
  if (netHist.length > 1500) netHist.splice(0, netHist.length - 1500);
  return ro;
}
let lastNetRo = null;
function drawNetwork(ro) {
  const cv = $('net'), { ctx, w, h } = fitCanvas(cv);
  ctx.clearRect(0, 0, w, h); ctx.fillStyle = P.slab; ctx.fillRect(0, 0, w, h);
  const n = net.n, R = Math.min(h * 0.31, w / (n * 2.6)), cxs = [...Array(n)].map((_, k) => w * (k + 0.5) / n), cy = Math.max(R * 1.45, h * 0.42);
  for (const [j, k] of edgesOf()) {
    ctx.beginPath(); ctx.strokeStyle = net.J || net.Jzz ? P.gold : P.rule; ctx.lineWidth = 1 + 5 * Math.max(net.J, net.Jzz); ctx.setLineDash(net.J || net.Jzz ? [] : [4, 5]);
    if (Math.abs(j - k) === 1) { ctx.moveTo(cxs[j] + R, cy); ctx.lineTo(cxs[k] - R, cy); }
    else { const top = Math.max(3, cy - R * 1.42); ctx.moveTo(cxs[j], cy - R); ctx.bezierCurveTo(cxs[j], top, cxs[k], top, cxs[k], cy - R); }
    ctx.stroke(); ctx.setLineDash([]);
  }
  ro.forEach((r, v) => { drawFace(ctx, cxs[v], cy, R, { hand: netUnwrap[v].acc, retention: r.retention, flash: netFlash[v], small: true, caption: `τ = ${(netObj.lapses[v] * net.t).toFixed(2)}`, sub: `N = ${netObj.lapses[v].toFixed(4)}` }, P); netFlash[v] = Math.max(0, netFlash[v] - 0.05); });
}
function netReadouts(ro) {
  ro.forEach((r, v) => { const el = $(`nr-${v}`); if (!el) return; el.textContent = fmt(netUnwrap[v].acc / clock.p.omega, 4); $(`ns-${v}`).textContent = `τ = ${fmt(netObj.lapses[v] * net.t, 4)} · in E47 ${fmt(r.retention, 5)}`; });
  if (netHist.length > 1) {
    const series = [];
    for (let v = 1; v < net.n; v++) {
      series.push({ pts: netHist.map(hh => [hh.t, hh.phase[v] - hh.phase[0]]), kind: 'line', cls: v === 1 ? 'ln-gold' : 'ln-ink' });
      series.push({ pts: netHist.map(hh => [hh.t, (netObj.lapses[v] - netObj.lapses[0]) * clock.p.omega * hh.t]), kind: 'line', cls: 'ln-verd' });
    }
    plot($('net-plot'), { series, xlabel: 'coordinate time t', ylabel: 'phase − phase of clock 1', width: 720, height: 220 });
  }
}

// ======================================================================= construct (built here)
function construct() {
  const t0 = performance.now(), car = new Q.Carrier(), mult = car.multiplicities(), lam = [0, 2, 6, 12, 20, 30, 42];
  bars($('fig-spec'), { labels: lam, values: lam.map(l => mult[l]), selected: [6, 30], xlabel: 'Casimir eigenvalue j(j+1)', ylabel: 'multiplicity' });
  const ns = [1, 2, 5, 10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220], norms = K.contractionNorms(car, ns);
  plot($('fig-contr'), { logY: true, xlabel: 'n', ylabel: '‖Γⁿ − P‖₂', series: [{ pts: [...Array(221).keys()].map(n => [n, (15 / 17) ** n]), kind: 'line', cls: 'ln-ink' }, { pts: ns.map(n => [n, norms[n]]), kind: 'dots', cls: 'dot-gold' }] });
  const Pm = car.P, idem = Q.linalg.rmaxabs(Q.linalg.radd(Q.linalg.rmul(Pm, Pm), Pm, -1)), KP = Q.linalg.rmaxabs(Q.linalg.rmul(car.K, Pm));
  const node = new Q.ClockMemoryNode(car), red = node.reductionResiduals({ omega: 2, c: 1, gamma: 20, clock: 'matched' }), pm = Q.pullbackMetric(car);
  const cval = v => { const Cv = Q.linalg.rmatvec(car.C, v); return v.reduce((s, z, k) => s + z * Cv[k], 0); };
  const cells = [
    ['carrier', '125 = 5 × 5 × 5', `13 components · ${car.components().map(c => c.length).join(' ')}`],
    ['projector P47', `rank ${Q.linalg.rtrace(Pm).toFixed(12)}`, `‖P² − P‖ = ${sci(idem)} · ‖KP‖ = ${sci(KP)}`],
    ['Q = K² gap', '11664 … 186624', 'ε* = 1/99144 · ρ* = 15/17'],
    ['clock states', [node.a, node.as, node.b, node.bs].map(v => cval(v).toFixed(10)).map(x => (+x).toString()).join(' · '), 'C of a, a′, b, b′'],
    ['exact reduction', sci(Math.max(red.H_invariance, red.R_invariance, red.Rdag_invariance, red.H_reduced)), '125-dim generator = 4-dim node'],
    ['metric g', pm.g.map((r, i) => r[i].toFixed(12)).map(x => (+x).toString()).join(', '), `diagonal of B4†(P6 − P0)B4 · off-diagonal ≤ ${sci(Math.max(...pm.g.flatMap((r, i) => r.map((x, j) => i === j ? 0 : Math.abs(x)))))}`],
  ];
  $('con-readout').innerHTML = cells.map(([k, v, s]) => `<div class="cell"><span class="k">${k}</span><span class="val mono">${v}</span><span class="sub">${s}</span></div>`).join('');
  $('cap-spec').textContent = `Casimir spectrum of V₂⊗³, built in ${Math.round(performance.now() - t0)} ms. The gold bands, C = 6 and C = 30, are E47: 25 + 22 = 47 of 125. Eigenvalues are integers to ${sci(car.integralityResidual)}.`;
}
const orbit = K.modsyncOrbit(240);
function primeFig() {
  const n = +$('prime-n').value; $('o-prime-n').textContent = n;
  ringFigure($('fig-prime'), { rings: [{ n: 47, r: 105, active: n % 47 }, { n: 5, r: 52, active: n % 5, cls: 'dot-verd', dot: 3 }],
    center: [`n = ${n}`, `outer · 47-cycle: ${n % 47}`, `inner · torus: ${n % 5}`, n > 0 && n % 235 === 0 ? 'both return together' : n % 47 === 0 && n ? '47-cycle back, torus not' : n % 5 === 0 && n ? 'torus back, 47-cycle not' : '', 'lcm(5, 47) = 235'] });
}
function syncFig() {
  const n = +$('sync-n').value; $('o-sync-n').textContent = n; const z = orbit[n];
  ringFigure($('fig-sync'), { rings: [{ n: 220, r: 105, active: n >= 2 ? (n - 2) % 220 : -1, dot: 1.4 }, { n: 11, r: 62, active: n >= 2 ? (n - 2) % 11 : -1, cls: 'dot-verd', dot: 2.6 }, { n: 20, r: 34, active: n >= 2 ? (n - 2) % 20 : -1, cls: 'dot-rose', dot: 2 }],
    center: [`z${n < 10 ? '₀₁₂₃₄₅₆₇₈₉'[n] : '_' + n} = ${z}`, n < 2 ? 'transient step' : `outer · 220-cycle: ${(n - 2) % 220}`, `middle · mod 47 (11): ${z % 47}`, `inner · mod 125 (20): ${z % 125}`, `(15/17)ⁿ = ${((15 / 17) ** n).toExponential(4)}`, n === 222 ? 'z₂₂₂ = z₂: returned' : ''] });
}
$('prime-n').addEventListener('input', primeFig); $('sync-n').addEventListener('input', syncFig);

// ======================================================================= certificate
let lastCert = null, committed = null;
function certRow(c) {
  const res = typeof c.residual === 'number' ? sci(c.residual) : (c.residual ?? '—'), tol = typeof c.tolerance === 'number' ? sci(c.tolerance, 0) : (c.tolerance ?? 'exact');
  const ev = c.evidence.includes('E0') ? 'e0' : c.evidence === 'E2' ? 'e2' : 'e1';
  return `<tr><td class="mono">${c.id}</td><td style="font-family:var(--serif);font-size:.98rem">${c.title}</td><td><span class="ev ${ev}">${c.evidence}</span></td><td class="num">${res}</td><td class="num">${tol}</td><td>${c.pass ? '<span class="pass">PASS</span>' : '<span class="fail">FAIL</span>'}</td></tr>`;
}
function showCert(cert, where) {
  lastCert = cert;
  const ok = cert.passed === cert.total;
  $('cert-status').textContent = `${cert.status} in ${(cert.elapsed_ms / 1000).toFixed(1)} s${where}`; $('cert-status').className = ok ? '' : 'bad';
  $('cert-badge').textContent = cert.status; $('cert-badge').className = 'badge ' + (ok ? 'ok' : 'bad');
  $('cert-meta').textContent = `${cert.engine} · digest ${cert.digest.slice(0, 16)}… · ${cert.environment?.runtime || ''}`;
  $('cert-table').tBodies[0].innerHTML = cert.checks.map(certRow).join('');
  $('cert-dl').disabled = false; $('cert-bar').style.width = '100%';
  compareCommitted();
  if (window.parent !== window) window.parent.postMessage({ type: 'QCLOCK_CERTIFICATE', status: cert.status, digest: cert.digest, code: cert.certificate_code }, '*');
}
async function compareCommitted() {
  try {
    committed = committed || await fetch('data/certificate.json').then(r => r.json());
    const a = lastCert.checks.map(c => c.id + ':' + c.pass).join(), b = committed.checks.map(c => c.id + ':' + c.pass).join();
    say('cert-compare', a === b ? `Same ${lastCert.total} checks and verdicts as the committed certificate (digest ${committed.digest.slice(0, 12)}…, ${committed.environment?.runtime || 'node'}).` : 'Differs from the committed certificate: see the verdict column.', a === b ? 'ok' : 'bad');
  } catch { say('cert-compare', 'Committed certificate not reachable from this copy of the page.'); }
}
function runCert() {
  $('cert-status').textContent = 'running…'; $('cert-status').className = ''; $('cert-badge').textContent = 'running'; $('cert-badge').className = 'badge'; $('cert-bar').style.width = '0';
  const environment = { runtime: (navigator.userAgentData?.brands?.map(b => b.brand + ' ' + b.version).join(', ') || navigator.userAgent).slice(0, 160), threads: navigator.hardwareConcurrency || null };
  const total = 48;
  try {
    const w = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
    w.onmessage = e => { const m = e.data; if (m.type === 'progress') { $('cert-bar').style.width = `${Math.min(100, 100 * m.k / total)}%`; } else if (m.type === 'done') { showCert(m.cert, ' in this browser'); w.terminate(); } else { $('cert-status').textContent = 'error: ' + m.message; $('cert-status').className = 'bad'; } };
    w.onerror = () => { w.terminate(); runCertMain(environment); };
    w.postMessage({ environment });
  } catch { runCertMain(environment); }
}
async function runCertMain(environment) {
  const golden = await fetch('data/golden.json').then(r => r.json()).catch(() => null);
  const cert = await Q.runCertificate({ golden, environment, onProgress: k => { $('cert-bar').style.width = `${100 * k / 48}%`; } });
  showCert(cert, ' in this browser');
}
$('cert-run').addEventListener('click', runCert);
$('cert-dl').addEventListener('click', () => lastCert && download(`${lastCert.certificate_code}.json`, JSON.stringify(lastCert, null, 1)));

// ======================================================================= the Matrix bridge
let mx = { receipt: null, report: null, result: null };
function updateBridgeClock() { $('mx-clock').textContent = `${clock.p.clock} · c = ${clock.p.c} · γ = ${clock.p.gamma} · ω = ${fmt(clock.p.omega, 3)}`; }
const badgeCls = s => ['REPLAY-VERIFIED', 'MEASURED-CONSISTENT'].includes(s) ? 'ok' : ['ORDER-DEPENDENT', 'MEASURED-ORDER-DEPENDENT'].includes(s) ? 'warn' : 'bad';
async function loadReceipt(text) {
  mx = { receipt: null, report: null, result: null }; $('mx-result').textContent = ''; $('mx-amps').innerHTML = '';
  ['mx-run', 'mx-export', 'mx-export-clock'].forEach(id => $(id).disabled = true);
  try {
    if (text.length > 100000) throw Error('Receipt too large');
    const x = JSON.parse(text), rep = await R.verifyMatrixReceipt(x);
    mx.receipt = x; mx.report = rep;
    reportLine(rep.status, badgeCls(rep.status), `${rep.detail}\nSHA-256 ${rep.digestOk ? 'verified' : 'does not match'} · ${rep.digest.slice(0, 16)}… · ${rep.gates} gate(s), last column ${rep.maxColumn} · replay residual ${sci(rep.residualColumnOrder)} in diagram order, ${sci(rep.residualListedOrder)} in placement order${rep.orderGap > 1e-12 ? ` · the two orders differ by ${sci(rep.orderGap)}` : ''}${rep.columnCollisions ? ` · ${rep.columnCollisions} column collision(s)` : ''}`);
    $('mx-run').disabled = false;
  } catch (e) { reportLine('REJECTED', 'bad', String(e.message).slice(0, 300)); }
}
function reportLine(status, cls, text) { // receipt-derived text only ever goes in as text, never as markup
  const host = $('mx-report'); host.innerHTML = '';
  const b = document.createElement('span'); b.className = 'badge ' + cls; b.textContent = status;
  const t = document.createElement('span'); t.style.whiteSpace = 'pre-line'; t.textContent = ' ' + text;
  host.append(b, t);
}
$('mx-file').addEventListener('change', async e => { const f = e.target.files[0]; if (f) await loadReceipt(await f.text()); e.target.value = ''; });
$('mx-paste').addEventListener('click', async () => { try { await loadReceipt(await navigator.clipboard.readText()); } catch (e) { $('mx-report').textContent = 'Clipboard unavailable: ' + e.message; } });
const sample = async (circuit, order) => R.seal({ schema: R.MATRIX_SCHEMA, source: 'matrix-quantum-simulator', exported_at: new Date().toISOString(), qubits: 3, basis: 'q0q1q2', initial_state: '000', circuit, statevector: R.matrixReplay(circuit, order), measurement: { collapsed: false, outcome: null }, norm: 1 });
$('mx-sample').addEventListener('click', async () => loadReceipt(JSON.stringify(await sample([{ gate: 'H', targets: [0], column: 0 }, { gate: 'CX', targets: [0, 1], column: 1 }, { gate: 'CX', targets: [1, 2], column: 2 }], 'column'))));
$('mx-sample-bad').addEventListener('click', async () => loadReceipt(JSON.stringify(await sample([{ gate: 'S', targets: [0], column: 1 }, { gate: 'H', targets: [0], column: 0 }, { gate: 'H', targets: [2], column: 0 }], 'insertion'))));
$('mx-ticks').addEventListener('input', e => { $('o-mx-ticks').textContent = e.target.value; });
const ampTable = (title, vec) => `<div><div class="eyebrow" style="margin-bottom:.4rem">${title}</div><div class="table-box"><table><tbody>${vec.map((z, i) => `<tr><td class="mono">|${i.toString(2).padStart(3, '0')}⟩</td><td class="num">${z[0].toFixed(6)} ${z[1] < 0 ? '−' : '+'} ${Math.abs(z[1]).toFixed(6)}i</td><td class="num">${(100 * (z[0] ** 2 + z[1] ** 2)).toFixed(3)} %</td></tr>`).join('')}</tbody></table></div></div>`;
$('mx-run').addEventListener('click', async () => {
  const k = +$('mx-q').value, ticks = +$('mx-ticks').value, p = { ...clock.p }, tau = ticks * (Math.PI / 4) / p.omega;
  const { H, R: Rn } = Q.reduced(p), S = Q.linalg.cexpm(Q.linalg.cscale(Q.liouvillian(H, [Rn]), tau));
  const rho16 = R.densityOf(R.encodeQubit(mx.receipt.statevector, k)), out = Q.linalg.cmat(16);
  for (let s = 0; s < 4; s++) for (let s2 = 0; s2 < 4; s2++) { const X = Q.linalg.cmat(4); for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) { X.re[i * 4 + j] = rho16.re[(i * 4 + s) * 16 + j * 4 + s2]; X.im[i * 4 + j] = rho16.im[(i * 4 + s) * 16 + j * 4 + s2]; } const Y = Q.applySuper(S, X); for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) { out.re[(i * 4 + s) * 16 + j * 4 + s2] = Y.re[i * 4 + j]; out.im[(i * 4 + s) * 16 + j * 4 + s2] = Y.im[i * 4 + j]; } }
  const L = R.logicalDensity(out, k), pur = R.purity(L);
  let psiT = mx.receipt.statevector.map(z => z.slice()); for (let n = 0; n < ticks; n++) psiT = R.matrixApply(psiT.map(([r, i]) => ({ r, i })), 'T', [k]).map(z => [z.r, z.i]);
  const fid = R.expectationIn(L, psiT); let inE = 0; for (let i = 0; i < 16; i++) if ((Math.floor(i / 4) & 1) === 0) inE += out.re[i * 16 + i];
  mx.result = { k, ticks, tau, purity: pur, fidelity: fid, rho: out, psiT, inE47: inE };
  const exact = 1 - fid < 1e-10 && 1 - pur < 1e-10, verified = ['REPLAY-VERIFIED', 'ORDER-DEPENDENT', 'MEASURED-CONSISTENT', 'MEASURED-ORDER-DEPENDENT'].includes(mx.report.status);
  const word = R.tickWord(ticks), fits = mx.report.maxColumn + 1 + word.length <= 8;
  $('mx-result').innerHTML = `Clock ran for proper time τ = ${fmt(tau, 6)} (${ticks} tick${ticks > 1 ? 's' : ''} of π/4 at ω = ${fmt(p.omega, 4)}) on q${k}, ${p.clock} clock, c = ${p.c}, γ = ${p.gamma}. In E47 at the end: ${fmt(inE, 6)}. Logical state purity <b>${fmt(pur, 15)}</b>, overlap with the Matrix’s T<sup>${ticks}</sup> on q${k}: <b>${fmt(fid, 15)}</b> (infidelity ${sci(1 - fid)}). ` +
    (exact ? `Exact: these ticks are the gate word <b>${word.join(' ') || 'I'}</b>.` : p.clock === 'linear' ? 'The linear clock dephases under recovery, so the result is mixed and cannot be written as Matrix gates.' : 'Not exact.') +
    (exact && !fits ? ' The Matrix grid has no free columns left for that word.' : '') + (exact && !verified ? ' The receipt did not replay, so it is not sent back.' : '');
  $('mx-amps').innerHTML = ampTable('receipt statevector', mx.receipt.statevector) + ampTable(`after ${ticks} tick(s) = T${ticks > 1 ? '^' + ticks : ''} on q${k}`, psiT);
  $('mx-export').disabled = !(exact && verified && fits && ticks % 8 !== 0); $('mx-export-clock').disabled = false;
});
$('mx-export').addEventListener('click', async () => {
  const { receipt, word } = await R.matrixReceiptWithTicks(mx.receipt, mx.result.k, mx.result.ticks);
  download('matrix-quantum-receipt-clock-ticks.json', JSON.stringify(receipt, null, 2));
  say('mx-result', `Exported: ${word.join(' ')} appended on q${mx.result.k}. In The Matrix choose IMPORT JSON ↑. Each gate is written as both gate and g, which is what the Matrix importer’s validator currently reads.`, 'ok');
});
$('mx-export-clock').addEventListener('click', async () => {
  const r = mx.result, payload = { schema: R.CLOCK_SCHEMA, source: 'quantum-clock', engine: Q.ENGINE_VERSION, exported_at: new Date().toISOString(),
    params: { omega: clock.p.omega, c: clock.p.c, gamma: clock.p.gamma, clock: clock.p.clock, syndrome: clock.p.syndrome },
    time: { tau: r.tau, ticks: r.ticks }, link: { matrix_digest: mx.report.digest, matrix_status: mx.report.status, qubit: r.k, encoding: R.ENCODING },
    state: { basis: 'joint index (2l+f)·4 + s; l: a=0, b=1; s: the other two qubits in order', dim: 16, density: Q.linalg.ctoRows(r.rho) },
    readouts: { logical_purity: r.purity, overlap_with_T_word: r.fidelity, in_E47: r.inE47 } };
  download('quantum-clock-receipt.json', JSON.stringify(await R.seal(payload), null, 2));
});

// ======================================================================= automate
let lastSweep = null;
$('sw-param').addEventListener('change', e => { const d = { gamma: [2, 60], c: [0.05, 3], omega: [0.2, 4], eta: [0, 1] }[e.target.value]; $('sw-from').value = d[0]; $('sw-to').value = d[1]; });
$('sw-run').addEventListener('click', () => {
  try {
    const param = $('sw-param').value, from = +$('sw-from').value, to = +$('sw-to').value, steps = Math.max(2, Math.min(400, Math.round(+$('sw-steps').value))), log = $('sw-log').value === '1';
    if (!(Number.isFinite(from) && Number.isFinite(to))) throw Error('from and to must be numbers');
    if (log && (from <= 0 || to <= 0)) throw Error('logarithmic spacing needs positive limits');
    if (param === 'eta' && (from < 0 || to > 1)) throw Error('η lies in [0, 1]');
    const t0 = performance.now(), base = { omega: clock.p.omega, c: clock.p.c, gamma: clock.p.gamma, eta: clock.p.eta, clock: $('sw-clock').value };
    const res = Q.sweep({ param, from, to, steps, log, base }); lastSweep = res;
    const xs = res.rows.map(r => r[param]);
    plot($('sw-fig1'), { xlabel: Q.SWEEP_PARAMS[param], ylabel: 'population in E47', series: [
      { pts: res.rows.map((r, k) => [xs[k], r.retention_a_closed]), kind: 'line', cls: 'ln-gold' }, { pts: res.rows.map((r, k) => [xs[k], r.retention_a]), kind: 'dots', cls: 'dot-gold' },
      { pts: res.rows.map((r, k) => [xs[k], r.retention_b_closed]), kind: 'line', cls: 'ln-ink' }, { pts: res.rows.map((r, k) => [xs[k], r.retention_b]), kind: 'dots', cls: 'dot-ink' }] });
    const matched = base.clock === 'matched';
    plot($('sw-fig2'), { xlabel: Q.SWEEP_PARAMS[param], ylabel: matched ? 'dephasing (absolute)' : 'dephasing γφ', logY: !matched, series: matched
      ? [{ pts: res.rows.map((r, k) => [xs[k], Math.abs(r.dephasing)]), kind: 'dots', cls: 'dot-verd' }]
      : [{ pts: res.rows.map((r, k) => [xs[k], r.dephasing_adiabatic]), kind: 'line', cls: 'ln-gold' }, { pts: res.rows.map((r, k) => [xs[k], r.dephasing]), kind: 'dots', cls: 'dot-verd' }] });
    const worst = Math.max(...res.rows.map(r => Math.max(Math.abs(r.retention_a - r.retention_a_closed), Math.abs(r.retention_b - r.retention_b_closed))));
    say('sw-msg', `${steps} points in ${Math.round(performance.now() - t0)} ms · steady state vs closed form: worst ${sci(worst)}${matched ? ` · matched dephasing ≤ ${sci(Math.max(...res.rows.map(r => Math.abs(r.dephasing))))}` : ''}`, 'ok');
    $('sw-csv').disabled = false; $('sw-json').disabled = false;
  } catch (e) { say('sw-msg', e.message, 'bad'); }
});
$('sw-csv').addEventListener('click', () => lastSweep && download(`qclock-sweep-${lastSweep.param}.csv`, Q.toCSV(lastSweep), 'text/csv'));
$('sw-json').addEventListener('click', () => lastSweep && download(`qclock-sweep-${lastSweep.param}.json`, JSON.stringify(lastSweep, null, 1)));

// ======================================================================= scripting API
window.QuantumClock = { ...Q, clock: () => clock, network: () => netObj };
window.addEventListener('message', e => {
  const m = e.data; if (!m || typeof m !== 'object' || typeof m.type !== 'string' || !m.type.startsWith('QCLOCK_')) return;
  const reply = msg => { try { e.source?.postMessage(msg, e.origin && e.origin !== 'null' ? e.origin : '*'); } catch { /* window gone */ } };
  if (m.type === 'QCLOCK_SET' && m.params && typeof m.params === 'object') {
    const allowed = {}; for (const k of ['omega', 'c', 'gamma', 'eta']) if (Number.isFinite(m.params[k])) allowed[k] = m.params[k];
    if (['linear', 'matched'].includes(m.params.clock)) allowed.clock = m.params.clock;
    if (Array.isArray(m.params.velocity) && m.params.velocity.length === 3 && m.params.velocity.every(Number.isFinite)) allowed.velocity = m.params.velocity;
    reconfigure(allowed, { reset: !!m.reset }); reply({ type: 'QCLOCK_STATE', state: clock.state(), params: clock.p });
  } else if (m.type === 'QCLOCK_GET') reply({ type: 'QCLOCK_STATE', state: clock.state(), params: clock.p });
});

// ======================================================================= loop
const face = $('face');
let lastFrame = performance.now(), lastReadout = 0;
function frame(now) {
  const dt = Math.min(0.05, (now - lastFrame) / 1000); lastFrame = now;
  if (running) {
    const fresh = clock.advance(dt * speed); if (fresh.length) flash = 1;
    const s = clock.state();
    if (s.tau - lastHistTau > 0.04 || fresh.length) { hist.push({ tau: s.tau, err: s.phaseTime - s.tau, errC: s.phaseTimeConditional - s.tau, ret: s.deterministic.retention, click: fresh.length > 0 }); lastHistTau = s.tau; if (hist.length > 3000) hist.splice(0, hist.length - 3000); }
  }
  const s = clock.state(), { ctx, w } = fitCanvas(face);
  ctx.clearRect(0, 0, w, w);
  drawFace(ctx, w / 2, w / 2, w / 2 - 4, { hand: s.hand, handC: clock.p.eta > 0 ? s.handConditional : null, retention: s.deterministic.retention, flash, caption: `τ = ${s.tau.toFixed(2)}` }, P);
  flash = Math.max(0, flash - dt * 2.5);
  if (net.running && netObj) lastNetRo = netAdvance(dt * speed);
  if (netObj && lastNetRo) drawNetwork(lastNetRo);
  if (now - lastReadout > 120) { clockReadouts(); drawStrip($('strip'), hist, P); if (lastNetRo) netReadouts(lastNetRo); lastReadout = now; }
  requestAnimationFrame(frame);
}

// ======================================================================= start
syncControls(); updateBridgeClock(); buildNetwork(); lastNetRo = netObj.readout(netRho);
if (!running) { $('play').textContent = 'Play'; $('play').setAttribute('aria-pressed', 'false'); }
if (!net.running) { $('net-play').textContent = 'Play'; $('net-play').setAttribute('aria-pressed', 'false'); }
requestAnimationFrame(frame);
setTimeout(() => { try { construct(); } catch (e) { $('cap-spec').textContent = 'Construction failed: ' + e.message; } primeFig(); syncFig(); runCert(); $('sw-run').click(); }, 60);
