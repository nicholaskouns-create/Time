// Dense real and complex linear algebra on row-major Float64Array storage.
// Complex matrices are {n, m, re, im}; real matrices are {n, m, a}.
// Everything here is small (≤ 256) and written for exactness first, speed second.

export const cmat = (n, m = n) => ({ n, m, re: new Float64Array(n * m), im: new Float64Array(n * m) });
export const rmat = (n, m = n) => ({ n, m, a: new Float64Array(n * m) });

export function ceye(n) { const A = cmat(n); for (let i = 0; i < n; i++) A.re[i * n + i] = 1; return A; }
export function reye(n) { const A = rmat(n); for (let i = 0; i < n; i++) A.a[i * n + i] = 1; return A; }

export function cclone(A) { return { n: A.n, m: A.m, re: A.re.slice(), im: A.im.slice() }; }
export function rclone(A) { return { n: A.n, m: A.m, a: A.a.slice() }; }

export function cfromReal(R) { const A = cmat(R.n, R.m); A.re.set(R.a); return A; }

export function cfromRows(rows) {
  // rows: array of arrays of numbers or [re, im] pairs
  const n = rows.length, m = rows[0].length, A = cmat(n, m);
  for (let i = 0; i < n; i++) for (let j = 0; j < m; j++) {
    const z = rows[i][j];
    if (Array.isArray(z)) { A.re[i * m + j] = z[0]; A.im[i * m + j] = z[1]; } else A.re[i * m + j] = z;
  }
  return A;
}

export function ctoRows(A) {
  const out = [];
  for (let i = 0; i < A.n; i++) { const r = []; for (let j = 0; j < A.m; j++) r.push([A.re[i * A.m + j], A.im[i * A.m + j]]); out.push(r); }
  return out;
}

export function cmul(A, B) {
  if (A.m !== B.n) throw Error(`cmul shape ${A.n}x${A.m} · ${B.n}x${B.m}`);
  const n = A.n, k = A.m, m = B.m, C = cmat(n, m);
  const ar = A.re, ai = A.im, br = B.re, bi = B.im, cr = C.re, ci = C.im;
  for (let i = 0; i < n; i++) {
    const io = i * k, co = i * m;
    for (let p = 0; p < k; p++) {
      const xr = ar[io + p], xi = ai[io + p];
      if (xr === 0 && xi === 0) continue;
      const po = p * m;
      for (let j = 0; j < m; j++) {
        const yr = br[po + j], yi = bi[po + j];
        cr[co + j] += xr * yr - xi * yi;
        ci[co + j] += xr * yi + xi * yr;
      }
    }
  }
  return C;
}

export function rmul(A, B) {
  if (A.m !== B.n) throw Error('rmul shape');
  const n = A.n, k = A.m, m = B.m, C = rmat(n, m), a = A.a, b = B.a, c = C.a;
  for (let i = 0; i < n; i++) {
    const io = i * k, co = i * m;
    for (let p = 0; p < k; p++) {
      const x = a[io + p];
      if (x === 0) continue;
      const po = p * m;
      for (let j = 0; j < m; j++) c[co + j] += x * b[po + j];
    }
  }
  return C;
}

// Y <- Y + (sr + i si) X   (returns Y)
export function caxpy(sr, si, X, Y) {
  const xr = X.re, xi = X.im, yr = Y.re, yi = Y.im;
  for (let k = 0; k < xr.length; k++) { yr[k] += sr * xr[k] - si * xi[k]; yi[k] += sr * xi[k] + si * xr[k]; }
  return Y;
}
export function cadd(A, B) { return caxpy(1, 0, B, cclone(A)); }
export function csub(A, B) { return caxpy(-1, 0, B, cclone(A)); }
export function cscale(A, sr, si = 0) { const B = cmat(A.n, A.m); return caxpy(sr, si, A, B); }
export function radd(A, B, s = 1) { const C = rclone(A); for (let k = 0; k < C.a.length; k++) C.a[k] += s * B.a[k]; return C; }
export function rscale(A, s) { const C = rclone(A); for (let k = 0; k < C.a.length; k++) C.a[k] *= s; return C; }

export function cadj(A) {
  const B = cmat(A.m, A.n);
  for (let i = 0; i < A.n; i++) for (let j = 0; j < A.m; j++) { B.re[j * A.n + i] = A.re[i * A.m + j]; B.im[j * A.n + i] = -A.im[i * A.m + j]; }
  return B;
}
export function rtrans(A) { const B = rmat(A.m, A.n); for (let i = 0; i < A.n; i++) for (let j = 0; j < A.m; j++) B.a[j * A.n + i] = A.a[i * A.m + j]; return B; }

export function ckron(A, B) {
  const n = A.n * B.n, m = A.m * B.m, C = cmat(n, m);
  for (let i = 0; i < A.n; i++) for (let j = 0; j < A.m; j++) {
    const ar = A.re[i * A.m + j], ai = A.im[i * A.m + j];
    if (ar === 0 && ai === 0) continue;
    for (let k = 0; k < B.n; k++) for (let l = 0; l < B.m; l++) {
      const br = B.re[k * B.m + l], bi = B.im[k * B.m + l], o = (i * B.n + k) * m + (j * B.m + l);
      C.re[o] = ar * br - ai * bi; C.im[o] = ar * bi + ai * br;
    }
  }
  return C;
}
export function rkron(A, B) {
  const n = A.n * B.n, m = A.m * B.m, C = rmat(n, m);
  for (let i = 0; i < A.n; i++) for (let j = 0; j < A.m; j++) {
    const x = A.a[i * A.m + j]; if (x === 0) continue;
    for (let k = 0; k < B.n; k++) for (let l = 0; l < B.m; l++) C.a[(i * B.n + k) * m + (j * B.m + l)] = x * B.a[k * B.m + l];
  }
  return C;
}

export function ctrace(A) { let r = 0, i = 0; for (let k = 0; k < A.n; k++) { r += A.re[k * A.m + k]; i += A.im[k * A.m + k]; } return [r, i]; }
export function rtrace(A) { let r = 0; for (let k = 0; k < A.n; k++) r += A.a[k * A.m + k]; return r; }
export function cnormF(A) { let s = 0; for (let k = 0; k < A.re.length; k++) s += A.re[k] * A.re[k] + A.im[k] * A.im[k]; return Math.sqrt(s); }
export function rnormF(A) { let s = 0; for (let k = 0; k < A.a.length; k++) s += A.a[k] * A.a[k]; return Math.sqrt(s); }
export function cnorm1(A) {
  let best = 0;
  for (let j = 0; j < A.m; j++) { let s = 0; for (let i = 0; i < A.n; i++) s += Math.hypot(A.re[i * A.m + j], A.im[i * A.m + j]); if (s > best) best = s; }
  return best;
}
export function cmaxabs(A) { let b = 0; for (let k = 0; k < A.re.length; k++) b = Math.max(b, Math.hypot(A.re[k], A.im[k])); return b; }
export function rmaxabs(A) { let b = 0; for (let k = 0; k < A.a.length; k++) b = Math.max(b, Math.abs(A.a[k])); return b; }

export function cmatvec(A, xr, xi) {
  const yr = new Float64Array(A.n), yi = new Float64Array(A.n);
  for (let i = 0; i < A.n; i++) {
    let sr = 0, si = 0; const o = i * A.m;
    for (let j = 0; j < A.m; j++) { const ar = A.re[o + j], ai = A.im[o + j]; sr += ar * xr[j] - ai * xi[j]; si += ar * xi[j] + ai * xr[j]; }
    yr[i] = sr; yi[i] = si;
  }
  return [yr, yi];
}
export function rmatvec(A, x) { const y = new Float64Array(A.n); for (let i = 0; i < A.n; i++) { let s = 0; const o = i * A.m; for (let j = 0; j < A.m; j++) s += A.a[o + j] * x[j]; y[i] = s; } return y; }

// Solve A X = B (complex) by LU with partial pivoting. A is n x n, B is n x m.
export function csolve(A, B) {
  const n = A.n, m = B.m;
  const ar = A.re.slice(), ai = A.im.slice(), xr = B.re.slice(), xi = B.im.slice();
  for (let c = 0; c < n; c++) {
    let p = c, best = -1;
    for (let r = c; r < n; r++) { const v = Math.hypot(ar[r * n + c], ai[r * n + c]); if (v > best) { best = v; p = r; } }
    if (best === 0) throw Error('csolve: singular matrix');
    if (p !== c) {
      for (let j = 0; j < n; j++) { let t = ar[c * n + j]; ar[c * n + j] = ar[p * n + j]; ar[p * n + j] = t; t = ai[c * n + j]; ai[c * n + j] = ai[p * n + j]; ai[p * n + j] = t; }
      for (let j = 0; j < m; j++) { let t = xr[c * m + j]; xr[c * m + j] = xr[p * m + j]; xr[p * m + j] = t; t = xi[c * m + j]; xi[c * m + j] = xi[p * m + j]; xi[p * m + j] = t; }
    }
    const pr = ar[c * n + c], pi = ai[c * n + c], den = pr * pr + pi * pi;
    for (let r = c + 1; r < n; r++) {
      const vr = ar[r * n + c], vi = ai[r * n + c];
      if (vr === 0 && vi === 0) continue;
      const fr = (vr * pr + vi * pi) / den, fi = (vi * pr - vr * pi) / den;
      for (let j = c; j < n; j++) { const yr = ar[c * n + j], yi = ai[c * n + j]; ar[r * n + j] -= fr * yr - fi * yi; ai[r * n + j] -= fr * yi + fi * yr; }
      for (let j = 0; j < m; j++) { const yr = xr[c * m + j], yi = xi[c * m + j]; xr[r * m + j] -= fr * yr - fi * yi; xi[r * m + j] -= fr * yi + fi * yr; }
    }
  }
  for (let c = n - 1; c >= 0; c--) {
    const pr = ar[c * n + c], pi = ai[c * n + c], den = pr * pr + pi * pi;
    for (let j = 0; j < m; j++) {
      let sr = xr[c * m + j], si = xi[c * m + j];
      for (let k = c + 1; k < n; k++) { const ur = ar[c * n + k], ui = ai[c * n + k], vr = xr[k * m + j], vi = xi[k * m + j]; sr -= ur * vr - ui * vi; si -= ur * vi + ui * vr; }
      xr[c * m + j] = (sr * pr + si * pi) / den; xi[c * m + j] = (si * pr - sr * pi) / den;
    }
  }
  return { n, m, re: xr, im: xi };
}

// Matrix exponential: scaling and squaring with the degree-13 Padé approximant (Higham 2005).
const PADE13 = [64764752532480000, 32382376266240000, 7771770303897600, 1187353796428800, 129060195264000,
  10559470521600, 670442572800, 33522128640, 1323241920, 40840800, 960960, 16380, 182, 1];
const THETA13 = 5.371920351148152;
export function cexpm(A0) {
  const n = A0.n, nrm = cnorm1(A0);
  const s = nrm > THETA13 ? Math.max(0, Math.ceil(Math.log2(nrm / THETA13))) : 0;
  const A = s ? cscale(A0, 2 ** -s) : A0;
  const b = PADE13, I = ceye(n);
  const A2 = cmul(A, A), A4 = cmul(A2, A2), A6 = cmul(A4, A2);
  const lin = (terms) => { const R = cmat(n); for (const [c, X] of terms) caxpy(c, 0, X, R); return R; };
  const U1 = cmul(A6, lin([[b[13], A6], [b[11], A4], [b[9], A2]]));
  const U = cmul(A, caxpy(1, 0, U1, lin([[b[7], A6], [b[5], A4], [b[3], A2], [b[1], I]])));
  const V1 = cmul(A6, lin([[b[12], A6], [b[10], A4], [b[8], A2]]));
  const V = caxpy(1, 0, V1, lin([[b[6], A6], [b[4], A4], [b[2], A2], [b[0], I]]));
  let R = csolve(csub(V, U), cadd(V, U));
  for (let k = 0; k < s; k++) R = cmul(R, R);
  return R;
}

// Cyclic Jacobi eigen-decomposition of a real symmetric matrix. Returns {values, vectors (columns)}.
export function rsymEig(A0, { tol = 1e-15, maxSweeps = 60 } = {}) {
  const n = A0.n, A = A0.a.slice(), V = new Float64Array(n * n);
  for (let i = 0; i < n; i++) V[i * n + i] = 1;
  let off0 = 0; for (let k = 0; k < A.length; k++) off0 += A[k] * A[k];
  const scale = Math.sqrt(off0) || 1;
  for (let sweep = 0; sweep < maxSweeps; sweep++) {
    let off = 0;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += A[p * n + q] * A[p * n + q];
    if (Math.sqrt(off) <= tol * scale) break;
    for (let p = 0; p < n - 1; p++) for (let q = p + 1; q < n; q++) {
      const apq = A[p * n + q];
      if (Math.abs(apq) <= 1e-300) continue;
      const app = A[p * n + p], aqq = A[q * n + q];
      const theta = (aqq - app) / (2 * apq);
      const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
      const c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < n; k++) {
        const akp = A[k * n + p], akq = A[k * n + q];
        A[k * n + p] = c * akp - s * akq; A[k * n + q] = s * akp + c * akq;
      }
      for (let k = 0; k < n; k++) {
        const apk = A[p * n + k], aqk = A[q * n + k];
        A[p * n + k] = c * apk - s * aqk; A[q * n + k] = s * apk + c * aqk;
      }
      A[p * n + q] = 0; A[q * n + p] = 0;
      for (let k = 0; k < n; k++) {
        const vkp = V[k * n + p], vkq = V[k * n + q];
        V[k * n + p] = c * vkp - s * vkq; V[k * n + q] = s * vkp + c * vkq;
      }
    }
  }
  const idx = [...Array(n).keys()].sort((i, j) => A[i * n + i] - A[j * n + j]);
  const values = idx.map(i => A[i * n + i]);
  const vectors = rmat(n);
  for (let c = 0; c < n; c++) for (let r = 0; r < n; r++) vectors.a[r * n + c] = V[r * n + idx[c]];
  return { values, vectors };
}

// Hermitian eigen-decomposition through the real symmetric embedding [[X,-Y],[Y,X]].
// Each eigenvalue appears twice there; one vector per pair is kept by Gram–Schmidt in C^n.
export function hermEig(H) {
  const n = H.n, R = rmat(2 * n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    const x = H.re[i * n + j], y = H.im[i * n + j];
    R.a[i * 2 * n + j] = x; R.a[i * 2 * n + n + j] = -y;
    R.a[(n + i) * 2 * n + j] = y; R.a[(n + i) * 2 * n + n + j] = x;
  }
  const { values, vectors } = rsymEig(R);
  const outVals = [], cols = [];
  for (let c = 0; c < 2 * n && cols.length < n; c++) {
    let vr = new Float64Array(n), vi = new Float64Array(n);
    for (let r = 0; r < n; r++) { vr[r] = vectors.a[r * 2 * n + c]; vi[r] = vectors.a[(n + r) * 2 * n + c]; }
    for (const [ur, ui] of cols) { // remove components along kept vectors (complex inner product)
      let pr = 0, pi = 0;
      for (let r = 0; r < n; r++) { pr += ur[r] * vr[r] + ui[r] * vi[r]; pi += ur[r] * vi[r] - ui[r] * vr[r]; }
      for (let r = 0; r < n; r++) { vr[r] -= pr * ur[r] - pi * ui[r]; vi[r] -= pr * ui[r] + pi * ur[r]; }
    }
    let nn = 0; for (let r = 0; r < n; r++) nn += vr[r] * vr[r] + vi[r] * vi[r];
    if (nn < 1e-12) continue;
    nn = Math.sqrt(nn); for (let r = 0; r < n; r++) { vr[r] /= nn; vi[r] /= nn; }
    cols.push([vr, vi]); outVals.push(values[c]);
  }
  const V = cmat(n);
  cols.forEach(([vr, vi], c) => { for (let r = 0; r < n; r++) { V.re[r * n + c] = vr[r]; V.im[r * n + c] = vi[r]; } });
  return { values: outVals, vectors: V };
}

// Roots of a monic complex polynomial z^d + c[d-1] z^(d-1) + ... + c[0] (Durand–Kerner, then Newton polish).
export function polyRoots(cr, ci) {
  const d = cr.length;
  const evalP = (zr, zi) => { let pr = 1, pi = 0; for (let k = d - 1; k >= 0; k--) { const tr = pr * zr - pi * zi + cr[k], ti = pr * zi + pi * zr + ci[k]; pr = tr; pi = ti; } return [pr, pi]; };
  const evalDP = (zr, zi) => { // derivative
    let pr = d, pi = 0;
    for (let k = d - 1; k >= 1; k--) { const tr = pr * zr - pi * zi + k * cr[k], ti = pr * zi + pi * zr + k * ci[k]; pr = tr; pi = ti; }
    return [pr, pi];
  };
  let R = 1; for (let k = 0; k < d; k++) R = Math.max(R, 1 + Math.hypot(cr[k], ci[k]));
  let zr = [], zi = [];
  for (let k = 0; k < d; k++) { const a = 2 * Math.PI * k / d + 0.4; zr.push(R * 0.9 * Math.cos(a)); zi.push(R * 0.9 * Math.sin(a)); }
  for (let it = 0; it < 2000; it++) {
    let delta = 0;
    for (let k = 0; k < d; k++) {
      const [pr, pi] = evalP(zr[k], zi[k]);
      let qr = 1, qi = 0;
      for (let j = 0; j < d; j++) if (j !== k) { const ar = zr[k] - zr[j], ai = zi[k] - zi[j]; const tr = qr * ar - qi * ai, ti = qr * ai + qi * ar; qr = tr; qi = ti; }
      const den = qr * qr + qi * qi; if (den === 0) continue;
      const wr = (pr * qr + pi * qi) / den, wi = (pi * qr - pr * qi) / den;
      zr[k] -= wr; zi[k] -= wi; delta = Math.max(delta, Math.hypot(wr, wi));
    }
    if (delta < 1e-17 * R) break;
  }
  for (let k = 0; k < d; k++) for (let it = 0; it < 4; it++) { // Newton polish
    const [pr, pi] = evalP(zr[k], zi[k]), [dr, di] = evalDP(zr[k], zi[k]), den = dr * dr + di * di;
    if (den === 0) break;
    zr[k] -= (pr * dr + pi * di) / den; zi[k] -= (pi * dr - pr * di) / den;
  }
  return zr.map((r, k) => [r, zi[k]]);
}

// Eigenvalues of a small complex matrix via its characteristic polynomial (Faddeev–LeVerrier).
export function smallEigvals(A) {
  const n = A.n; let M = cmat(n); const cr = new Float64Array(n), ci = new Float64Array(n);
  let ckr = 1, cki = 0; // c_n = 1
  for (let k = 1; k <= n; k++) {
    // M_k = A M_{k-1} + c_{n-k+1} I ;  c_{n-k} = -tr(A M_k)/k
    const AM = cmul(A, M);
    for (let i = 0; i < n; i++) { AM.re[i * n + i] += ckr; AM.im[i * n + i] += cki; }
    M = AM;
    const [tr, ti] = ctrace(cmul(A, M));
    ckr = -tr / k; cki = -ti / k;
    cr[n - k] = ckr; ci[n - k] = cki;
  }
  return polyRoots(Array.from(cr), Array.from(ci));
}
