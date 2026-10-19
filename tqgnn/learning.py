"""QGRNN-style learning for the TQGNN: recover lapses and couplings from clock data.

Logical sector only (no leakage): each node is a clock qubit {|a>, |b>} with
energies (6, 30) omega/24; edges carry exchange J_e.  Coordinate-time data are
the node coherences Tr(|a><b|_v rho(t)).  The model generator is

    H_t = sum_v N_v omega_v h_v + sum_e N_e J_e x_e,     N_e = (N_j + N_k)/2.

Only the products N_v omega_v and N_e J_e enter, so coordinate-time data cannot
separate a lapse from a clock frequency; identical calibrated clocks
(omega_v = omega known) make the lapses identifiable.
"""
import numpy as np
from scipy.optimize import least_squares


def _ops(n):
    a2b = np.array([[0, 0], [1, 0]], dtype=complex)
    h = np.diag([6.0, 30.0]) / 24.0

    def emb(op, v):
        out = np.eye(1, dtype=complex)
        for w in range(n):
            out = np.kron(out, op if w == v else np.eye(2))
        return out

    hs = [emb(h, v) for v in range(n)]
    sp = [emb(a2b, v) for v in range(n)]
    coh = [emb(a2b.T, v) for v in range(n)]   # |a><b|
    return hs, sp, coh


def simulate(n, edges, N, omegas, Js, psi0, times):
    hs, sp, coh = _ops(n)
    H = sum(N[v] * omegas[v] * hs[v] for v in range(n))
    for e, (j, k) in enumerate(edges):
        Ne = 0.5 * (N[j] + N[k])
        hop = sp[j] @ sp[k].conj().T
        H = H + Ne * Js[e] * (hop + hop.conj().T)
    w, U = np.linalg.eigh(H)
    c0 = U.conj().T @ psi0
    out = []
    for t in times:
        psi = U @ (np.exp(-1j * w * t) * c0)
        out.append([np.vdot(psi, C @ psi) for C in coh])
    return np.array(out)   # (times, n) complex


def residuals(data, n, edges, N, omegas, Js, psi0, times):
    sim = simulate(n, edges, N, omegas, Js, psi0, times)
    d = sim - data
    return np.concatenate([d.real.ravel(), d.imag.ravel()])


def fit_calibrated(data, n, edges, omega, psi0, times, x0):
    """Unknowns: lapses N_v and couplings J_e; clocks calibrated (omega known)."""
    ne = len(edges)

    def r(x):
        return residuals(data, n, edges, x[:n], [omega] * n, x[n:n + ne], psi0, times)

    return least_squares(r, x0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=2000)


def phase_slope_lapses(data, times, omega, t_max=2.0):
    """Initial lapse estimate from each clock's early phase slope (phi_v ~ omega N_v t)."""
    sel = times <= t_max
    phi = np.unwrap(-np.angle(data[sel]), axis=0)
    return np.array([np.polyfit(times[sel], phi[:, v], 1)[0] / omega for v in range(data.shape[1])])


def train(data, n, edges, omega, psi0, times, J0=0.05, windows=(2.0, 5.0, 10.0, 20.0, None)):
    """Calibrated training: phase-slope initialisation, then least squares on
    growing time windows (warm start), as for QGRNN Hamiltonian learning."""
    x = np.concatenate([phase_slope_lapses(data, times, omega), np.full(len(edges), J0)])
    res = None
    for T in windows:
        sel = times <= (times[-1] if T is None else T)
        res = fit_calibrated(data[sel], n, edges, omega, psi0, times[sel], x)
        x = res.x
    return res


def jacobian_rank(n, edges, N, omegas, Js, psi0, times, h=1e-6, tol=1e-6):
    """Numerical rank of d(data)/d(N, omega, J) at the given point."""
    x0 = np.concatenate([N, omegas, Js])
    ne = len(edges)

    def r(x):
        return residuals(np.zeros((len(times), n)), n, edges, x[:n], x[n:2 * n], x[2 * n:2 * n + ne], psi0, times)

    cols = []
    for i in range(len(x0)):
        e = np.zeros_like(x0); e[i] = h
        cols.append((r(x0 + e) - r(x0 - e)) / (2 * h))
    Jm = np.column_stack(cols)
    s = np.linalg.svd(Jm, compute_uv=False)
    return int(np.sum(s > tol * s[0])), s
