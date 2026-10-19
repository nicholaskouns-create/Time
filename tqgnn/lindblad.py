"""GKSL utilities: right-hand side, integration, Liouvillian, counting statistics.

Vectorisation is row-major: vec(A X B) = (A kron B^T) vec(X).
"""
import numpy as np
from scipy.integrate import solve_ivp


def gksl_rhs(H, Ls, rho):
    """-i[H, rho] + sum_k (L rho L^dag - 1/2 {L^dag L, rho}) for Hermitian rho."""
    G = sum((L.conj().T @ L for L in Ls), np.zeros_like(H))
    X = (H - 0.5j * G) @ rho
    out = -1j * (X - X.conj().T)
    for L in Ls:
        out = out + L @ rho @ L.conj().T
    return out


def evolve(generator, rho0, t_eval, rtol=1e-11, atol=1e-13):
    """Integrate d rho/dt = L_t rho with DOP853; generator(t) -> (H, [L_k])."""
    d = rho0.shape[0]

    def f(t, y):
        rho = y.view(complex).reshape(d, d)
        H, Ls = generator(t)
        return gksl_rhs(H, Ls, rho).reshape(-1).view(float)

    y0 = np.ascontiguousarray(rho0, dtype=complex).reshape(-1).view(float).copy()
    sol = solve_ivp(f, (t_eval[0], t_eval[-1]), y0, method="DOP853",
                    t_eval=t_eval, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    return np.array([np.ascontiguousarray(sol.y[:, k]).view(complex).reshape(d, d)
                     for k in range(len(t_eval))])


def liouvillian(H, Ls):
    d = H.shape[0]
    I = np.eye(d)
    Lv = -1j * (np.kron(H, I) - np.kron(I, H.T))
    for L in Ls:
        LdL = L.conj().T @ L
        Lv = Lv + np.kron(L, L.conj()) - 0.5 * np.kron(LdL, I) - 0.5 * np.kron(I, LdL.T)
    return Lv


def jump_super(L):
    return np.kron(L, L.conj())


def steady_state(H, Ls):
    """Unique stationary state (null vector of the Liouvillian), normalised."""
    d = H.shape[0]
    ev, vec = np.linalg.eig(liouvillian(H, Ls))
    k = np.argmin(np.abs(ev))
    rho = vec[:, k].reshape(d, d)
    rho = rho / np.trace(rho)
    return 0.5 * (rho + rho.conj().T)


def counting_rates(H, Ls, counted, eta=1.0, h=1e-2):
    """Long-time mean rate and Fano factor of detected jumps.

    Tilted generator L_s = L + (e^s - 1) eta sum_{k in counted} J_k; theta(s) is its
    eigenvalue of largest real part; mean = theta'(0), variance rate = theta''(0).
    """
    Lv = liouvillian(H, Ls)
    Jsum = sum(jump_super(Ls[k]) for k in counted)

    def theta(s):
        ev = np.linalg.eigvals(Lv + (np.exp(s) - 1.0) * eta * Jsum)
        return ev[np.argmax(ev.real)].real

    t = {k: theta(k * h) for k in (-2, -1, 0, 1, 2)}
    d1 = (-t[2] + 8 * t[1] - 8 * t[-1] + t[-2]) / (12 * h)
    d2 = (-t[2] + 16 * t[1] - 30 * t[0] + 16 * t[-1] - t[-2]) / (12 * h * h)
    return d1, d2 / d1
