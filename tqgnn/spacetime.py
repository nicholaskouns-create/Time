"""Background metric and lapse.

pullback_metric() reproduces the canonical runner's construction
g = B4^dagger (P6 - P0) B4 for the four-frame B4 = [u, v1, v2, v3]
(u: the Casimir singlet, v_i: SO(3)-orbit tangents in the C = 6 copy W).
The metric is background structure for the network: it supplies each
node's lapse N = dtau/dt along a timelike worldline.
"""
import numpy as np
from scipy.integrate import quad

from .carrier import Carrier

METRIC = np.diag([-1.0, 2.0, 8.0, 2.0])


def pullback_metric(carrier=None):
    """Return (g, intertwiner_residual, u, W) from the Casimir-spectral pullback."""
    car = carrier or Carrier()
    j = 2
    m = car.m_single
    tol = 1e-8

    # W: the explicit V2 copy (pair singlet of spins 1,2) (x) spin 3, inside C = 6
    pair_singlet = np.zeros(25, dtype=complex)
    for ia, ma in enumerate(m):
        ib = int(np.where(m == -ma)[0][0])
        pair_singlet[ia * 5 + ib] = (-1) ** int(j - ma) / np.sqrt(2 * j + 1)
    W = np.column_stack([np.kron(pair_singlet, np.eye(5, dtype=complex)[:, k]) for k in range(5)])

    # u: the singlet (C = 0), the time direction of the frame
    u = car.basis(0)[:, 0].astype(complex)

    # real spin-2 model Sym_0(3,R) and its tangent vectors under so(3)
    Qm = np.diag([1.0, 0.0, -1.0])
    A1 = np.array([[0, 0, 0], [0, 0, -1], [0, 1, 0]], float)
    A2 = np.array([[0, 0, 1], [0, 0, 0], [-1, 0, 0]], float)
    A3 = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]], float)
    A = [A1, A2, A3]
    D = [Ai @ Qm - Qm @ Ai for Ai in A]

    E1 = np.zeros((3, 3)); E1[0, 0] = 1; E1[1, 1] = -1; E1 /= np.sqrt(2)
    E2 = np.zeros((3, 3)); E2[0, 0] = 1; E2[1, 1] = 1; E2[2, 2] = -2; E2 /= np.sqrt(6)
    E3 = np.zeros((3, 3)); E3[0, 1] = E3[1, 0] = 1; E3 /= np.sqrt(2)
    E4 = np.zeros((3, 3)); E4[0, 2] = E4[2, 0] = 1; E4 /= np.sqrt(2)
    E5 = np.zeros((3, 3)); E5[1, 2] = E5[2, 1] = 1; E5 /= np.sqrt(2)
    Eb = [E1, E2, E3, E4, E5]
    rho = np.array([[[np.trace(Ea @ (Ai @ Eb_ - Eb_ @ Ai)) for Eb_ in Eb] for Ea in Eb] for Ai in A])
    JW = [W.conj().T @ Jt @ W for Jt in (car.Jx, car.Jy, car.Jz)]

    # intertwiner Phi: Sym_0(3,R) -> W, solving Phi rho_i + i JW_i Phi = 0
    rows = [np.kron(rho[i].T, np.eye(5)) + np.kron(np.eye(5), 1j * JW[i]) for i in range(3)]
    Ms = np.vstack(rows)
    gew, gev = np.linalg.eigh(Ms.conj().T @ Ms)
    if not (gew[0] < 1e-10 and gew[1] > 1e-3):
        raise RuntimeError("intertwiner nullity failure")
    Phi = gev[:, 0].reshape(5, 5, order="F")
    resid = max(np.linalg.norm(Phi @ rho[i] + 1j * JW[i] @ Phi) for i in range(3))
    Phi = Phi * np.sqrt(5 / np.trace(Phi.conj().T @ Phi))

    dvecs = [np.array([np.trace(Ea @ Di) for Ea in Eb]) for Di in D]
    vs = [W @ (Phi @ dv) for dv in dvecs]
    P0 = car.projector(0)
    P6 = car.projector(6)
    B4 = np.column_stack([u] + vs)
    g = np.real_if_close(B4.conj().T @ (P6 - P0) @ B4, tol=1e6)
    return np.real(g), resid, u, W


def lapse(velocity, g=METRIC):
    """N = dtau/dt = sqrt(-g(xdot, xdot)) for xdot = (1, v); raises if not timelike."""
    xdot = np.concatenate([[1.0], np.asarray(velocity, dtype=float)])
    s = -xdot @ g @ xdot
    if s <= 0:
        raise ValueError("worldline is not timelike for this metric")
    return float(np.sqrt(s))


class Worldline:
    """A worldline given by its coordinate velocity v(t); supplies N(t) and tau(t)."""

    def __init__(self, velocity, g=METRIC):
        self._v = velocity if callable(velocity) else (lambda t, v=np.asarray(velocity, float): v)
        self.g = g

    def N(self, t):
        return lapse(self._v(t), self.g)

    def tau(self, t):
        val, _ = quad(self.N, 0.0, t, epsabs=1e-13, epsrel=1e-13, limit=200)
        return val
