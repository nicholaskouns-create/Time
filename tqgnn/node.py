"""E47 clock-memory node.

Two clock states inside E47 (|a> in C = 6, |b> in C = 30) and two syndrome
states outside it (|a'>, |b'>).  V = [a, b] and W = [a', b'] are isometries,
so span{V, W} = logical (x) flag with flag 0 = protected, flag 1 = leaked.

Generator on the full 125-dim carrier (proper time):
    H = f(C) + c (W V^dag + V W^dag),     R = sqrt(gamma) V W^dag
f(C) is the clock.  The review's linear clock is f(C) = omega C / 24, which
makes the clock tick at omega inside E47: f(30) - f(6) = omega.
Because every operator above maps span{V, W} into itself, the dynamics of
any state in that span is exactly the 4-dim reduced dynamics.
"""
import numpy as np

from .carrier import Carrier, SPECTRUM

REF_SITE = 38  # site with m = (1, 0, -1), total M = 0


def _pair_singlet_copy(car, k=2):
    """|a>: pair singlet of spins 1,2 (x) |m3>, inside C = 6 (canonical runner's W copy)."""
    m = car.m_single
    pair = np.zeros(25)
    for ia, ma in enumerate(m):
        ib = int(np.where(m == -ma)[0][0])
        pair[ia * 5 + ib] = (-1) ** int(2 - ma) / np.sqrt(5)
    return np.kron(pair, np.eye(5)[:, k])


class ClockMemoryNode:
    def __init__(self, carrier=None, syndrome=(12, 42)):
        car = carrier or Carrier()
        self.car = car
        self.syndrome = tuple(syndrome)
        if any(s in (6, 30) for s in syndrome) or any(s not in SPECTRUM for s in syndrome):
            raise ValueError("syndrome sectors must be Casimir eigenvalues outside E47")
        e = np.zeros(125)
        e[REF_SITE] = 1.0

        def proj(lam):
            v = car.projector(lam) @ e
            n = np.linalg.norm(v)
            if n < 1e-3:
                raise ValueError(f"reference site has no C = {lam} component")
            return v / n

        self.a = _pair_singlet_copy(car)
        self.b = proj(30)
        self.a_s = proj(syndrome[0])
        self.b_s = proj(syndrome[1])
        self.V = np.column_stack([self.a, self.b])
        self.W = np.column_stack([self.a_s, self.b_s])
        # reduced index 2*l + f : (a,0), (a,1), (b,0), (b,1)
        self.U = np.column_stack([self.a, self.a_s, self.b, self.b_s])
        self.sector_of = {0: 6, 1: syndrome[0], 2: 30, 3: syndrome[1]}

    # -- clocks ------------------------------------------------------------
    def clock_function(self, omega, clock="linear"):
        """Return f: eigenvalue -> energy.  'linear': omega*C/24.
        'matched': linear on every sector except the second syndrome sector, which is
        moved so that f(b') - f(a') = f(30) - f(6) (equal detunings for a and b)."""
        lin = {lam: omega * lam / 24.0 for lam in SPECTRUM}
        if clock == "linear":
            return lin
        if clock == "matched":
            f = dict(lin)
            f[self.syndrome[1]] = lin[self.syndrome[0]] + (lin[30] - lin[6])
            return f
        raise ValueError(clock)

    def detunings(self, omega, clock="linear"):
        f = self.clock_function(omega, clock)
        return f[self.syndrome[0]] - f[6], f[self.syndrome[1]] - f[30]

    # -- operators ---------------------------------------------------------
    def reduced(self, omega, c, gamma, clock="linear"):
        """4x4 (H, R) on logical (x) flag."""
        f = self.clock_function(omega, clock)
        H = np.diag([f[self.sector_of[k]] for k in range(4)]).astype(complex)
        X = np.array([[0, 1], [1, 0]], dtype=complex)
        down = np.array([[0, 1], [0, 0]], dtype=complex)   # |0><1| on the flag
        H = H + c * np.kron(np.eye(2), X)
        R = np.sqrt(gamma) * np.kron(np.eye(2), down)
        return H, R

    def branch(self, l, omega, c, gamma, clock="linear"):
        """2x2 (H, R) of logical branch l (0 = a, 1 = b) on the flag; the reduced
        operators are block diagonal in l, so each branch evolves on its own."""
        f = self.clock_function(omega, clock)
        H = np.diag([f[self.sector_of[2 * l]], f[self.sector_of[2 * l + 1]]]).astype(complex)
        H = H + c * np.array([[0, 1], [1, 0]], dtype=complex)
        R = np.sqrt(gamma) * np.array([[0, 1], [0, 0]], dtype=complex)
        return H, R

    def full(self, omega, c, gamma, clock="linear"):
        """125x125 (H, R) on the carrier."""
        f = self.clock_function(omega, clock)
        H = self.car.function_of_C(lambda lam: f[lam]).astype(complex)
        H = H + c * (self.W @ self.V.T + self.V @ self.W.T)
        R = np.sqrt(gamma) * (self.V @ self.W.T)
        return H, R.astype(complex)


# closed forms ---------------------------------------------------------------
def stationary_leak(c, gamma, delta=0.0):
    """Stationary leaked population of one branch: 4c^2 / (gamma^2 + 8c^2 + 4 delta^2)."""
    return 4 * c**2 / (gamma**2 + 8 * c**2 + 4 * delta**2)


def resonant_fano(c, gamma):
    """Asymptotic Fano factor of recovery clicks at zero detuning."""
    return 1 - 24 * c**2 * gamma**2 / (gamma**2 + 8 * c**2) ** 2


def adiabatic_coherence_rate(c, gamma, delta_a, delta_b, omega):
    """Second-order (c << gamma) rate of the clock coherence rho_ab = <a|rho|b>:

        d rho_ab/dt = lambda rho_ab,
        lambda = i omega + i c^2 (z_a - conj(z_b)) + gamma c^2 z_a conj(z_b),
        z_l = 1 / (Delta_l - i gamma/2).

    -Re(lambda) is the clock dephasing rate, Im(lambda) the clock frequency."""
    za = 1 / (delta_a - 0.5j * gamma)
    zb = 1 / (delta_b - 0.5j * gamma)
    return 1j * omega + 1j * c**2 * (za - np.conj(zb)) + gamma * c**2 * za * np.conj(zb)


def adiabatic_dephasing(c, gamma, delta_a, delta_b):
    """gamma_phi = gamma c^2 |z_a - z_b|^2 / 2; zero exactly when Delta_a = Delta_b."""
    za = 1 / (delta_a - 0.5j * gamma)
    zb = 1 / (delta_b - 0.5j * gamma)
    return 0.5 * gamma * c**2 * abs(za - zb) ** 2
