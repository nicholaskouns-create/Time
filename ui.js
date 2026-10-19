// Drawing for the Quantum Clock page: clock faces (canvas), the history strip (canvas), small SVG plots.

export function palette() {
  const cs = getComputedStyle(document.documentElement), v = n => cs.getPropertyValue(n).trim();
  return { ground: v('--ground'), slab: v('--slab'), slab2: v('--slab-2'), rule: v('--rule'), ink: v('--ink'), ink2: v('--ink-2'), ink3: v('--ink-3'), gold: v('--gold'), goldSoft: v('--gold-soft'), verd: v('--verd'), verdSoft: v('--verd-soft'), rose: v('--rose'), mono: v('--mono') || 'monospace', serif: v('--serif') || 'serif' };
}

export function fitCanvas(cv) {
  const dpr = Math.min(3, window.devicePixelRatio || 1), r = cv.getBoundingClientRect();
  const w = Math.max(1, Math.round(r.width * dpr)), h = Math.max(1, Math.round(r.height * dpr));
  if (cv.width !== w || cv.height !== h) { cv.width = w; cv.height = h; }
  const ctx = cv.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx, w: r.width, h: r.height };
}

export const GATE_WORDS = ['I', 'T', 'S', 'ST', 'Z', 'ZT', 'ZS', 'ZST'];
const polar = (cx, cy, r, a) => [cx + r * Math.sin(a), cy - r * Math.cos(a)];
const mod2pi = a => ((a % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);

// st: { hand, handC, retention, flash, tau, caption, small }
export function drawFace(ctx, cx, cy, R, st, P) {
  ctx.save();
  ctx.beginPath(); ctx.arc(cx, cy, R, 0, 2 * Math.PI); ctx.fillStyle = P.slab; ctx.fill();
  ctx.lineWidth = 1; ctx.strokeStyle = P.rule; ctx.stroke();
  // prime 47-cycle rim (declared bridge: one prime tick per 2π/47 of phase)
  const prime = Math.floor(mod2pi(st.hand) / (2 * Math.PI / 47) + 1e-9) % 47;
  for (let k = 0; k < 47; k++) {
    const a = 2 * Math.PI * k / 47, [x1, y1] = polar(cx, cy, R * 0.955, a), [x2, y2] = polar(cx, cy, R * (k === prime ? 0.885 : 0.925), a);
    ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x2, y2); ctx.strokeStyle = k === prime ? P.gold : P.ink3; ctx.lineWidth = k === prime ? 2 : 0.8; ctx.stroke();
  }
  // eight T ticks with their Matrix gate words
  for (let k = 0; k < 8; k++) {
    const a = k * Math.PI / 4, [x1, y1] = polar(cx, cy, R * 0.86, a), [x2, y2] = polar(cx, cy, R * 0.76, a);
    ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x2, y2); ctx.strokeStyle = P.ink2; ctx.lineWidth = 1.6; ctx.stroke();
    if (!st.small) {
      const [tx, ty] = polar(cx, cy, R * 0.66, a);
      ctx.fillStyle = P.ink2; ctx.font = `500 ${Math.max(9, R * 0.062)}px ${P.mono}`; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.fillText(GATE_WORDS[k], tx, ty);
    }
  }
  // population in E47: arc of length retention, the leaked remainder in rose
  const rr = R * (st.small ? 0.5 : 0.47), ret = Math.max(0, Math.min(1, st.retention));
  ctx.lineWidth = Math.max(2, R * 0.022);
  ctx.beginPath(); ctx.arc(cx, cy, rr, -Math.PI / 2, -Math.PI / 2 + 2 * Math.PI * ret); ctx.strokeStyle = P.ink3; ctx.stroke();
  if (ret < 1) { ctx.beginPath(); ctx.arc(cx, cy, rr, -Math.PI / 2 + 2 * Math.PI * ret, 3 * Math.PI / 2); ctx.strokeStyle = P.rose; ctx.stroke(); }
  // click flash
  if (st.flash > 0) { ctx.beginPath(); ctx.arc(cx, cy, R * (0.12 + 0.4 * (1 - st.flash)), 0, 2 * Math.PI); ctx.strokeStyle = P.verd; ctx.globalAlpha = st.flash; ctx.lineWidth = 2; ctx.stroke(); ctx.globalAlpha = 1; }
  // hands
  const hand = (a, len, w, col) => { const [x, y] = polar(cx, cy, len, a), [bx, by] = polar(cx, cy, -len * 0.12, a); ctx.beginPath(); ctx.moveTo(bx, by); ctx.lineTo(x, y); ctx.strokeStyle = col; ctx.lineWidth = w; ctx.lineCap = 'round'; ctx.stroke(); };
  hand(st.hand, R * 0.84, Math.max(2.5, R * 0.03), P.gold);
  if (st.handC !== undefined && st.handC !== null) hand(st.handC, R * 0.74, Math.max(1.4, R * 0.012), P.verd);
  ctx.beginPath(); ctx.arc(cx, cy, Math.max(3, R * 0.03), 0, 2 * Math.PI); ctx.fillStyle = P.gold; ctx.fill();
  if (st.caption && !st.small) { ctx.fillStyle = P.ink; ctx.font = `500 ${Math.max(12, R * 0.085)}px ${P.serif}`; ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.fillText(st.caption, cx, cy + R * 0.24); }
  if (st.small) { // captions under small faces
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    if (st.caption) { ctx.fillStyle = P.ink; ctx.font = `500 ${Math.max(12, R * 0.16)}px ${P.serif}`; ctx.fillText(st.caption, cx, cy + R + 6); }
    if (st.sub) { ctx.fillStyle = P.ink3; ctx.font = `400 ${Math.max(10, R * 0.1)}px ${P.mono}`; ctx.fillText(st.sub, cx, cy + R + 8 + Math.max(14, R * 0.2)); }
  }
  ctx.restore();
  return prime;
}

// history strip: hand errors (top), population in E47 (middle), clicks (bottom)
export function drawStrip(cv, hist, P, window = 40) {
  const { ctx, w, h } = fitCanvas(cv);
  ctx.clearRect(0, 0, w, h); ctx.fillStyle = P.slab; ctx.fillRect(0, 0, w, h);
  if (!hist.length) return;
  const tEnd = hist[hist.length - 1].tau, t0 = Math.max(0, tEnd - window), pad = 8, lane1 = [pad, h * 0.52], lane2 = [h * 0.58, h * 0.8], lane3 = [h * 0.86, h - pad];
  const X = t => 46 + (w - 54) * (t - t0) / Math.max(1e-9, window);
  const pts = hist.filter(p => p.tau >= t0);
  let emax = 1e-12; for (const p of pts) emax = Math.max(emax, Math.abs(p.err), Math.abs(p.errC));
  const scale = 10 ** Math.ceil(Math.log10(emax * 1.05));
  const Y1 = e => (lane1[0] + lane1[1]) / 2 - (lane1[1] - lane1[0]) / 2 * e / scale;
  ctx.strokeStyle = P.rule; ctx.lineWidth = 1; ctx.setLineDash([2, 4]);
  ctx.beginPath(); ctx.moveTo(46, Y1(0)); ctx.lineTo(w - 8, Y1(0)); ctx.stroke(); ctx.setLineDash([]);
  ctx.fillStyle = P.ink3; ctx.font = `400 10px ${P.mono}`; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
  ctx.fillText('+' + scale.toExponential(0), 42, lane1[0] + 5); ctx.fillText('−' + scale.toExponential(0), 42, lane1[1] - 5); ctx.fillText('error', 42, Y1(0));
  const line = (get, col, Y, dash = []) => { ctx.beginPath(); pts.forEach((p, k) => { const x = X(p.tau), y = Y(get(p)); k ? ctx.lineTo(x, y) : ctx.moveTo(x, y); }); ctx.strokeStyle = col; ctx.lineWidth = 1.4; ctx.setLineDash(dash); ctx.stroke(); ctx.setLineDash([]); };
  line(p => p.err, P.gold, Y1); line(p => p.errC, P.verd, Y1, [4, 3]);
  let rmin = 1; for (const p of pts) rmin = Math.min(rmin, p.ret); rmin = Math.min(rmin, 0.999);
  const Y2 = r => lane2[1] - (lane2[1] - lane2[0]) * (r - rmin) / (1 - rmin);
  ctx.fillText('in E47', 42, (lane2[0] + lane2[1]) / 2);
  line(p => p.ret, P.ink2, Y2);
  ctx.fillText('clicks', 42, (lane3[0] + lane3[1]) / 2);
  ctx.strokeStyle = P.verd; ctx.lineWidth = 1.2;
  for (const p of pts) if (p.click) { const x = X(p.tau); ctx.beginPath(); ctx.moveTo(x, lane3[0]); ctx.lineTo(x, lane3[1]); ctx.stroke(); }
}

// ---------------------------------------------------------------- SVG plots
const NS = 'http://www.w3.org/2000/svg';
const el = (tag, attrs = {}, text) => { const e = document.createElementNS(NS, tag); for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v); if (text !== undefined) e.textContent = text; return e; };
const niceTicks = (lo, hi, n = 5) => { const span = hi - lo || 1, step0 = span / n, mag = 10 ** Math.floor(Math.log10(step0)), step = [1, 2, 5, 10].map(s => s * mag).find(s => span / s <= n) || 10 * mag; const out = []; for (let k = Math.ceil(lo / step - 1e-9); k * step <= hi + 1e-9 * span; k++) { const v = k * step; out.push(Math.abs(v) < step * 1e-9 ? 0 : +v.toPrecision(12)); } return out; };
const fmt = v => { const a = Math.abs(v); return a !== 0 && (a < 1e-3 || a >= 1e5) ? v.toExponential(0) : String(+v.toPrecision(4)); };

// series: [{ pts: [[x, y]], kind: 'line' | 'dots', cls }]
export function plot(svg, { series, xlabel = '', ylabel = '', logY = false, yDomain = null, xDomain = null, width = 520, height = 240 }) {
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
  const m = { l: 58, r: 14, t: 12, b: 38 }, W = width - m.l - m.r, H = height - m.t - m.b;
  const all = series.flatMap(s => s.pts).filter(([x, y]) => Number.isFinite(x) && Number.isFinite(y) && (!logY || y > 0));
  if (!all.length) return;
  let [x0, x1] = xDomain || [Math.min(...all.map(p => p[0])), Math.max(...all.map(p => p[0]))];
  let ys = all.map(p => logY ? Math.log10(p[1]) : p[1]), [y0, y1] = yDomain ? (logY ? yDomain.map(Math.log10) : yDomain) : [Math.min(...ys), Math.max(...ys)];
  if (y1 - y0 < 1e-12) { y0 -= 0.5; y1 += 0.5; } else if (!yDomain) { const pd = (y1 - y0) * 0.06; y0 -= pd; y1 += pd; }
  if (x1 - x0 < 1e-12) { x0 -= 0.5; x1 += 0.5; }
  const X = x => m.l + W * (x - x0) / (x1 - x0), Y = y => m.t + H - H * ((logY ? Math.log10(y) : y) - y0) / (y1 - y0);
  svg.append(el('line', { class: 'ax', x1: m.l, y1: m.t + H, x2: m.l + W, y2: m.t + H }), el('line', { class: 'ax', x1: m.l, y1: m.t, x2: m.l, y2: m.t + H }));
  for (const t of niceTicks(x0, x1, 6)) { svg.append(el('line', { class: 'gridl', x1: X(t), y1: m.t, x2: X(t), y2: m.t + H }), el('text', { class: 'tk', x: X(t), y: m.t + H + 14, 'text-anchor': 'middle' }, fmt(t))); }
  if (logY) { for (let e = Math.ceil(y0); e <= Math.floor(y1); e += Math.max(1, Math.ceil((y1 - y0) / 6))) { const y = m.t + H - H * (e - y0) / (y1 - y0); svg.append(el('line', { class: 'gridl', x1: m.l, y1: y, x2: m.l + W, y2: y }), el('text', { class: 'tk', x: m.l - 6, y: y + 3, 'text-anchor': 'end' }, `1e${e}`)); } }
  else for (const t of niceTicks(y0, y1, 5)) { svg.append(el('line', { class: 'gridl', x1: m.l, y1: Y(t), x2: m.l + W, y2: Y(t) }), el('text', { class: 'tk', x: m.l - 6, y: Y(t) + 3, 'text-anchor': 'end' }, fmt(t))); }
  svg.append(el('text', { class: 'lbl', x: m.l + W / 2, y: height - 6, 'text-anchor': 'middle' }, xlabel), el('text', { class: 'lbl', x: 12, y: m.t + H / 2, 'text-anchor': 'middle', transform: `rotate(-90 12 ${m.t + H / 2})` }, ylabel));
  for (const s of series) {
    const pts = s.pts.filter(([x, y]) => Number.isFinite(x) && Number.isFinite(y) && (!logY || y > 0));
    if (s.kind === 'dots') for (const [x, y] of pts) svg.append(el('circle', { class: s.cls || 'dot-gold', cx: X(x), cy: Y(y), r: s.r || 3 }));
    else if (pts.length) svg.append(el('path', { class: s.cls || 'ln-gold', d: pts.map(([x, y], k) => `${k ? 'L' : 'M'}${X(x).toFixed(2)},${Y(y).toFixed(2)}`).join('') }));
  }
}

export function bars(svg, { labels, values, selected, xlabel, ylabel, width = 520, height = 240 }) {
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const m = { l: 48, r: 10, t: 16, b: 40 }, W = width - m.l - m.r, H = height - m.t - m.b, vmax = Math.max(...values) * 1.15, bw = W / labels.length;
  svg.append(el('line', { class: 'ax', x1: m.l, y1: m.t + H, x2: m.l + W, y2: m.t + H }));
  labels.forEach((lab, k) => {
    const h = H * values[k] / vmax, x = m.l + k * bw + bw * 0.18, sel = selected.includes(lab);
    svg.append(el('rect', { class: 'bar' + (sel ? ' sel' : ''), x, y: m.t + H - h, width: bw * 0.64, height: h }));
    svg.append(el('text', { class: 'tk', x: x + bw * 0.32, y: m.t + H - h - 5, 'text-anchor': 'middle' }, String(values[k])));
    svg.append(el('text', { class: 'tk', x: x + bw * 0.32, y: m.t + H + 14, 'text-anchor': 'middle' }, String(lab)));
  });
  svg.append(el('text', { class: 'lbl', x: m.l + W / 2, y: height - 6, 'text-anchor': 'middle' }, xlabel), el('text', { class: 'lbl', x: 12, y: m.t + H / 2, 'text-anchor': 'middle', transform: `rotate(-90 12 ${m.t + H / 2})` }, ylabel));
}

export function ringFigure(svg, { width = 520, height = 260, rings, center }) {
  // rings: [{ n, r, active, cls, label }]
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const cx = width / 2 - 80, cy = height / 2;
  for (const ring of rings) {
    svg.append(el('circle', { cx, cy, r: ring.r, fill: 'none', class: 'ax' }));
    for (let k = 0; k < ring.n; k++) {
      const a = 2 * Math.PI * k / ring.n, x = cx + ring.r * Math.sin(a), y = cy - ring.r * Math.cos(a), on = ring.active === k;
      svg.append(el('circle', { cx: x, cy: y, r: on ? 5 : ring.dot || 2.2, class: on ? (ring.cls || 'dot-gold') : 'dot-ink', opacity: on ? 1 : 0.55 }));
    }
  }
  (center || []).forEach((line, k) => svg.append(el('text', { class: k ? 'tk' : 'lbl', x: width - 150, y: 40 + 22 * k, 'text-anchor': 'start' }, line)));
}
