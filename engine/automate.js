// Automation: parameter sweeps (closed form beside numeric), permalinks, exports.
import { branch, countingStats, detunings, stationaryLeak, resonantFano, coherenceModeRates, adiabaticCoherenceRate, adiabaticDephasing, fanoAtEfficiency } from './node.js';

export const SWEEP_PARAMS = { gamma: 'recovery rate γ', c: 'leak coupling c', omega: 'clock frequency ω', eta: 'detection efficiency η' };

export function sweepPoint(p) {
  const [da, db] = detunings(p.omega, p.clock, p.syndrome);
  const ba = branch(0, p), bb = branch(1, p);
  const sa = countingStats(ba.H, [ba.R], 0, p.eta), sb = countingStats(bb.H, [bb.R], 0, p.eta);
  const mode = coherenceModeRates(p).slowest, ad = adiabaticCoherenceRate(p.c, p.gamma, da, db, p.omega);
  return {
    delta_a: da, delta_b: db,
    retention_a: 1 - sa.steady.re[3], retention_a_closed: 1 - stationaryLeak(p.c, p.gamma, da),
    retention_b: 1 - sb.steady.re[3], retention_b_closed: 1 - stationaryLeak(p.c, p.gamma, db),
    click_rate_a: sa.rate, click_rate_b: sb.rate,
    fano_a: sa.fano, fano_b: sb.fano, fano_resonant_closed: fanoAtEfficiency(resonantFano(p.c, p.gamma), p.eta),
    dephasing: -mode[0], dephasing_adiabatic: adiabaticDephasing(p.c, p.gamma, da, db),
    clock_frequency: -mode[1], clock_frequency_adiabatic: ad[1],
    click_clock_rel_error_per_sqrt_tau: Math.sqrt(sa.fano / sa.rate),
  };
}

export function sweep({ param = 'gamma', from = 5, to = 80, steps = 16, log = false, base = {} }) {
  const p0 = { omega: 2, c: 1, gamma: 20, eta: 1, clock: 'linear', syndrome: [12, 42], ...base };
  const rows = [];
  for (let k = 0; k < steps; k++) {
    const s = steps === 1 ? 0 : k / (steps - 1), x = log ? from * (to / from) ** s : from + (to - from) * s;
    rows.push({ [param]: x, ...sweepPoint({ ...p0, [param]: x }) });
  }
  return { schema: 'QCLOCK-SWEEP-1.0', param, base: p0, columns: Object.keys(rows[0]), rows };
}

export function toCSV({ columns, rows }) { return [columns.join(','), ...rows.map(r => columns.map(k => r[k]).join(','))].join('\n') + '\n'; }

// permalinks: the page state lives in the URL hash
const KEYS = { w: 'omega', c: 'c', g: 'gamma', e: 'eta', k: 'clock', v: 'velocity', m: 'metric', th: 'theta', ph: 'phase', s: 'seed', sp: 'speed', view: 'view' };
export function encodeHash(p) {
  const out = [];
  for (const [short, long] of Object.entries(KEYS)) if (p[long] !== undefined) out.push(`${short}=${encodeURIComponent(Array.isArray(p[long]) ? p[long].join(',') : p[long])}`);
  return '#' + out.join('&');
}
export function decodeHash(hash) {
  const p = {};
  for (const part of hash.replace(/^#/, '').split('&')) {
    const [short, raw] = part.split('='); const long = KEYS[short]; if (!long || raw === undefined) continue;
    const val = decodeURIComponent(raw);
    if (long === 'velocity') { const v = val.split(',').map(Number); if (v.length === 3 && v.every(Number.isFinite)) p.velocity = v; }
    else if (long === 'clock') { if (['linear', 'matched'].includes(val)) p.clock = val; }
    else if (long === 'metric') { if (['e47', 'minkowski'].includes(val)) p.metric = val; }
    else if (long === 'view') { if (/^[a-z-]{1,24}$/.test(val)) p.view = val; }
    else if (Number.isFinite(Number(val))) p[long] = Number(val);
  }
  return p;
}
