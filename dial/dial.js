// A Dial of Forty-Seven Hours: the instrument and its live certificate.
// The dial is drawn from the proof; the certificate runs on the Quantum Clock's own engine (../engine).
import { Carrier } from '../engine/carrier.js';
import { runDialCertificate } from './dial-certificate.js';

const NS = 'http://www.w3.org/2000/svg';
const H = 47, CX = 280, CY = 280, R0 = 228;            // hours, centre, radius of the state points
const PLUS = new Set([...Array(13).keys(), ...Array.from({ length: 12 }, (_, i) => 35 + i)]);   // 25 gold hours, centred on 12 o'clock
const RIGHT = Array.from({ length: 11 }, (_, i) => [12 - i, 13 + i]);      // nested around 3 o'clock
const LEFT = Array.from({ length: 11 }, (_, i) => [35 + i, 34 - i]);       // nested around 9 o'clock
const PAIRS = RIGHT.flatMap((p, i) => [p, LEFT[i]]);                       // 22 pairs, outermost first
const hourOfBasis = k => (k < 25 ? (35 + k) % H : 13 + (k - 25));          // prime-clock basis e_k on the rim, clockwise

const ang = h => Math.PI / 2 - h * 2 * Math.PI / H;
const pt = (r, h) => [CX + r * Math.cos(ang(h)), CY - r * Math.sin(ang(h))];
const el = (tag, attrs = {}, parent) => { const e = document.createElementNS(NS, tag); for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v); if (parent) parent.appendChild(e); return e; };

function geodesic(h1, h2, n = 56) {                      // arc orthogonal to the rim between two hours
  const a = ang(h1), b = ang(h2);
  const p1 = [R0 * Math.cos(a), R0 * Math.sin(a)], p2 = [R0 * Math.cos(b), R0 * Math.sin(b)];
  let d = Math.abs(a - b) % (2 * Math.PI); if (d > Math.PI) d = 2 * Math.PI - d;
  const ux = p1[0] + p2[0], uy = p1[1] + p2[1], un = Math.hypot(ux, uy);
  const c = [ux / un * R0 / Math.cos(d / 2), uy / un * R0 / Math.cos(d / 2)], rho = R0 * Math.tan(d / 2);
  const f1 = Math.atan2(p1[1] - c[1], p1[0] - c[0]), f2 = Math.atan2(p2[1] - c[1], p2[0] - c[0]);
  let best = null;
  for (const sgn of [1, -1]) {
    const span = sgn === 1 ? ((f2 - f1) % (2 * Math.PI) + 2 * Math.PI) % (2 * Math.PI) : -(((f1 - f2) % (2 * Math.PI) + 2 * Math.PI) % (2 * Math.PI));
    const mx = c[0] + rho * Math.cos(f1 + span / 2), my = c[1] + rho * Math.sin(f1 + span / 2);
    if (!best || Math.hypot(mx, my) < best[0]) best = [Math.hypot(mx, my), span];
  }
  let s = '';
  for (let t = 0; t <= n; t++) { const f = f1 + best[1] * t / n; s += (t ? ' L' : 'M') + (CX + c[0] + rho * Math.cos(f)).toFixed(1) + ' ' + (CY - c[1] - rho * Math.sin(f)).toFixed(1); }
  return s;
}

// ---------------------------------------------------------------- the face
const svg = document.getElementById('dial');
el('circle', { cx: CX, cy: CY, r: 272, class: 'rim' }, svg);
el('circle', { cx: CX, cy: CY, r: 238, class: 'rim', 'stroke-width': 1, opacity: 0.6 }, svg);
const ticks = Array.from({ length: H }, (_, h) => { const [x1, y1] = pt(246, h), [x2, y2] = pt(264, h); return el('line', { x1, y1, x2, y2, class: 'tick' }, svg); });
el('circle', { cx: CX, cy: CY, r: 3, class: 'hub' }, svg);
const layer = el('g', {}, svg);

const state = { mode: 'grading', rank: 22, beat: 1, hours: [0, 9, 19, 28, 38, 14, 33], timer: null };
const reroll = () => { const pool = [...Array(H).keys()].sort(() => Math.random() - 0.5); state.hours = pool.slice(0, 7); };

function spoke(h, weight, cls, width) {
  const [x1, y1] = pt(R0, h), [x2, y2] = pt(R0 - 6 - 3.2 * weight, h);
  el('line', { x1, y1, x2, y2, class: 'spoke ' + cls, 'stroke-width': width }, layer);
  el('circle', { cx: x1, cy: y1, r: 4, class: 'spoke ' + cls, fill: 'currentColor', style: `fill: var(--${cls === 'ink' ? 'ink-2' : cls})` }, layer);
}
const label = (y, text, cls) => { const t = el('text', { x: CX, y, 'text-anchor': 'middle', class: cls }, layer); t.textContent = text; };

const MODES = {
  grading: {
    note: 'S = (18I − C)P/12 = −XP/6 marks the two sides: +1 on the 25 states of total spin 2, −1 on the 22 of total spin 5. It is a function of the Casimir, so it commutes with every SU(2) rotation. Its trace, 25 − 22 = 3, is the imbalance the rest of the proof turns on.',
    draw() { label(262, '25 | 22', 'big'); label(304, 'TR S = 3', 'small'); },
    cells: () => [['Operator', 'S, the grading'], ['Values', '+1 on 25 · −1 on 22'], ['Commutes with SU(2)', 'yes', 'ok'], ['Tr S', '3', 'hl'], ['S²', 'P'], ['Trace of P', '47']],
    ticks: 'sector',
  },
  pairing: {
    note: 'A pairing D swaps the two sides: D(u, v) = (A†v, Au) with A from the 25 to the 22. Rank r joins r states on each side; 47 − 2r = 3 + 2(22 − r) stay unpaired. Whatever r is, the unpaired states carry Tr S = 3: at full rank exactly three remain, all on the gold side. Any nonzero D breaks SU(2).',
    draw() {
      const used = PAIRS.slice(0, state.rank);
      for (const [g, b] of used) el('path', { d: geodesic(g, b), class: 'arc' }, layer);
      const paired = new Set(used.flat());
      for (let h = 0; h < H; h++) if (!paired.has(h)) { const [x, y] = pt(R0, h); el('circle', { cx: x, cy: y, r: 7.5, class: 'free ' + (PLUS.has(h) ? 'plus' : 'minus') }, layer); }
      label(132, String(H - 2 * state.rank), 'big'); label(156, 'UNPAIRED', 'small');
    },
    cells: () => [['Rank r', String(state.rank)], ['Paired states', String(2 * state.rank)], ['Unpaired', `${H - 2 * state.rank} = 3 + 2(22 − ${state.rank})`, 'hl'],
      ['Tr S on the unpaired', '3', 'hl'], ['Spectrum of D', state.rank ? `−1 ×${state.rank} · 0 ×${H - 2 * state.rank} · +1 ×${state.rank}` : '0 ×47'], ['Commutes with SU(2)', state.rank ? 'no' : 'yes (D = 0)', state.rank ? 'no' : 'ok']],
    ticks: 'sector', rank: true,
  },
  symmetric: {
    note: 'A dial that commutes with every SU(2) rotation is U = (B ⊗ I₅) ⊕ (F ⊗ I₁₁): one phase per copy of spin 2 (five copies, five states each) and per copy of spin 5 (two copies, eleven each). It can show at most seven hours. Play it: the hours hop as the beats go by, and never more than seven are lit.',
    draw() {
      const lit = state.hours.map((h, a) => [(h * state.beat) % H, a < 5 ? 5 : 11, a < 5 ? 'plus' : 'minus']);
      for (const [h, w, c] of lit) spoke(h, w, c, 4.5);
      label(300, String(new Set(lit.map(x => x[0])).size), 'big'); label(326, 'HOURS LIT', 'small');
    },
    cells: () => { const lit = new Set(state.hours.map(h => (h * state.beat) % H)); return [['Hours lit', String(lit.size), 'hl'], ['States per hour', '5 · 5 · 5 · 5 · 5 · 11 · 11'], ['Commutes with SU(2)', 'yes', 'ok'], ['Beat', `${state.beat % H} / 47`], ['Most hours possible', '5 + 2 = 7'], ['States', '25 + 22 = 47']]; },
    ticks: 'neutral', clock: true, reroll: true,
  },
  u0: {
    note: 'U₀ = ζ on the 25 and 1 on the 22, with ζ = e^(2πi/47). Its gold hand turns one hour per beat and comes home after 47 beats; the blue stays at the top. It commutes with SU(2), and it is exactly one 1/47 turn of the Quantum Clock’s own hand, up to a global phase.',
    draw() { spoke(0, 22, 'minus', 5); spoke(state.beat % H, 25, 'plus', 5); label(300, String(state.beat % H ? 2 : 1), 'big'); label(326, 'HOURS LIT', 'small'); },
    cells: () => [['Hours lit', String(state.beat % H ? 2 : 1), 'hl'], ['States per hour', '25 · 22'], ['Commutes with SU(2)', 'yes', 'ok'], ['Beat', `${state.beat % H} / 47`], ['Order', '47'], ['Hand', 'the grading S']],
    ticks: 'neutral', clock: true,
  },
  prime: {
    note: 'The prime 47-cycle moves e₀ → e₁ → … → e₄₆ → e₀ through a basis of E₄₇. Its spectrum is every 47th root of unity, all 47 hours lit at once, so by step 7 it cannot commute with SU(2). Watch its hand: it walks through the 25, crosses into the 22, and comes back.',
    draw() {
      for (let h = 0; h < H; h++) spoke(h, 1, 'ink', 2.2);
      const k = state.beat % H, hh = hourOfBasis(k), [x, y] = pt(R0 - 12, hh);
      el('line', { x1: CX, y1: CY, x2: x, y2: y, class: 'hand ' + (k < 25 ? 'plus' : 'minus') }, layer);
      label(300, '47', 'big'); label(326, 'HOURS LIT', 'small');
    },
    cells: () => { const k = state.beat % H; return [['Hours lit', '47', 'hl'], ['States per hour', '1 each'], ['Commutes with SU(2)', 'no', 'no'], ['Beat', `${k} / 47`], ['Hand at', `e${k} · ${k < 25 ? 'gold side, S = +1' : 'blue side, S = −1'}`], ['Keeps E₄₇', 'yes']]; },
    ticks: 'neutral', clock: true,
  },
};

const readout = document.getElementById('readout'), note = document.getElementById('mode-note');
function render() {
  const m = MODES[state.mode];
  layer.replaceChildren();
  ticks.forEach((t, h) => t.setAttribute('class', 'tick' + (m.ticks === 'sector' ? (PLUS.has(h) ? ' plus' : ' minus') : '')));
  m.draw();
  readout.replaceChildren(...m.cells().map(([k, v, cls]) => { const d = document.createElement('div'); d.className = 'cell' + (cls ? ' ' + cls : ''); d.innerHTML = '<span class="k"></span><span class="val"></span>'; d.children[0].textContent = k; d.children[1].textContent = v; return d; }));
  note.textContent = m.note;
  document.getElementById('rank-ctl').hidden = !m.rank;
  document.getElementById('clock-ctl').hidden = !m.clock;
  document.getElementById('reroll').hidden = !m.reroll;
  svg.setAttribute('aria-label', `The dial of 47 hours, ${state.mode} view. ` + m.cells().map(c => c[0] + ': ' + c[1]).join('; '));
}

function stop() { clearInterval(state.timer); state.timer = null; document.getElementById('play').textContent = 'Play'; }
document.querySelectorAll('.chip').forEach(b => b.addEventListener('click', () => {
  stop(); state.mode = b.dataset.mode; state.beat = 1;
  document.querySelectorAll('.chip').forEach(c => c.setAttribute('aria-pressed', String(c === b)));
  render();
}));
const rank = document.getElementById('rank'), rankOut = document.getElementById('rank-out');
rank.addEventListener('input', () => { state.rank = +rank.value; rankOut.value = rank.value; render(); });
document.getElementById('step').addEventListener('click', () => { state.beat = (state.beat + 1) % H; render(); });
document.getElementById('reroll').addEventListener('click', () => { reroll(); render(); });
document.getElementById('play').addEventListener('click', e => {
  if (state.timer) return stop();
  e.currentTarget.textContent = 'Pause';
  state.timer = setInterval(() => { state.beat = (state.beat + 1) % H; render(); }, 420);
});
render();

// ---------------------------------------------------------------- the live certificate
const fmt = c => {
  const v = Array.isArray(c.value) ? c.value.map(x => (Number.isInteger(x) ? x : +x.toFixed(6))).join(', ') : null;
  const r = typeof c.residual === 'number' ? (c.residual === 0 ? '0' : c.residual.toExponential(1)) : null;
  if (v && r) return `${v} · ${r}`;
  return r || v || Object.values(c.value).join(' · ');
};
setTimeout(() => {
  const bar = document.getElementById('cert-bar'), status = document.getElementById('cert-status');
  try {
    const t0 = performance.now(), cert = runDialCertificate(new Carrier()), ms = Math.round(performance.now() - t0);
    bar.style.width = '100%';
    status.textContent = `${cert.passed}/${cert.total} PASS · ${ms} ms`;
    if (cert.passed !== cert.total) status.classList.add('bad');
    document.getElementById('cert-body').replaceChildren(...cert.checks.map(c => {
      const tr = document.createElement('tr');
      tr.innerHTML = '<td></td><td></td><td></td><td class="num"></td><td></td>';
      const td = tr.children;
      td[0].textContent = c.id; td[1].textContent = c.title; td[2].textContent = c.evidence; td[3].textContent = fmt(c);
      td[4].innerHTML = c.pass ? '<span class="pass">PASS</span>' : '<span class="fail">FAIL</span>';
      return tr;
    }));
    document.getElementById('tower-body').replaceChildren(...cert.tower.map(t => {
      const tr = document.createElement('tr'); if (t.k === 5) tr.className = 'hl';
      const half = n => (n % 2 ? `${n}/2` : String(n / 2));
      tr.innerHTML = '<td></td><td></td><td></td><td class="num"></td><td class="num"></td><td class="num"></td><td class="num"></td><td class="num"></td>';
      [t.k, half(t.spinTimesTwo), half(Math.round((Math.sqrt(1 + 4 * 5 * ((t.spinTimesTwo / 2) * (t.spinTimesTwo / 2 + 1))) - 1))), t.plus, t.minus, t.plus + t.minus, (t.spinTimesTwo + 1) ** 3, t.index]
        .forEach((v, i) => { tr.children[i].textContent = String(v); });
      return tr;
    }));
  } catch (err) {
    status.textContent = 'did not run: ' + err.message; status.classList.add('bad');
  }
}, 30);
