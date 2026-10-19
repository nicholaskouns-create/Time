"""The E47 carrier V2^{(x)3} read as a 125-node graph.

Sites n = 25a + 5b + c with a, b, c in {0,...,4} carry the magnetic numbers
m = (2-a, 2-b, 2-c) (descending Condon-Shortley order, as in the canonical
spacetime runner).  The total Casimir C = J_tot^2 is a real symmetric graph
operator on these 125 nodes; E47 is its spectral band {6, 30}.
"""
from collections import deque

import numpy as np

# Casimir eigenvalue -> eigenspace dimension
SPECTRUM = {0: 1, 2: 9, 6: 25, 12: 28, 20: 27, 30: 22, 42: 13}
SELECTED = (6, 30)
FILTER_DENOM = 99144  # Gamma = I - K^2 / 99144


def spin_matrices(j=2):
    """Spin-j matrices in descending-m order."""
    m = np.arange(j, -j - 1, -1, dtype=float)
    d = len(m)
    Jz = np.diag(m).astype(complex)
    Jp = np.zeros((d, d), dtype=complex)
    for i in range(d - 1):
        mm = m[i + 1]
        Jp[i, i + 1] = np.sqrt(j * (j + 1) - mm * (mm + 1))
    Jm = Jp.conj().T
    return (Jp + Jm) / 2, (Jp - Jm) / 2j, Jz, m


class Carrier:
    """Three spin-2 systems: Casimir, selector, projector and spectral filter."""

    def __init__(self):
        Jx, Jy, Jz, m = spin_matrices(2)
        I5 = np.eye(5)

        def tot(A):
            return (np.kron(np.kron(A, I5), I5) + np.kron(np.kron(I5, A), I5)
                    + np.kron(np.kron(I5, I5), A))

        self.Jx, self.Jy, self.Jz = tot(Jx), tot(Jy), tot(Jz)
        C = self.Jx @ self.Jx + self.Jy @ self.Jy + self.Jz @ self.Jz
        assert np.allclose(C.imag, 0.0)
        self.C = C.real
        self.dim = 125
        self.m_single = m
        w, U = np.linalg.eigh(self.C)
        self.evals, self.evecs = w, U
        self.labels = np.rint(w).astype(int)
        I = np.eye(125)
        self.K = (self.C - 6 * I) @ (self.C - 30 * I)
        self.P = self.projector(*SELECTED)
        self.Q = I - self.P
        self.Gamma = I - self.K @ self.K / FILTER_DENOM

    # -- spectral data ---------------------------------------------------
    def basis(self, lam):
        """Orthonormal basis (125 x dim) of the Casimir eigenspace C = lam."""
        return self.evecs[:, self.labels == lam]

    def projector(self, *lams):
        B = np.hstack([self.basis(lam) for lam in lams])
        return B @ B.conj().T

    def function_of_C(self, f):
        """f(C) via spectral projectors; f maps eigenvalue -> number."""
        return sum(f(lam) * self.projector(lam) for lam in SPECTRUM)

    # -- graph view ------------------------------------------------------
    def site_m(self):
        """(125, 3) array of magnetic numbers (m1, m2, m3) per site."""
        idx = np.arange(125)
        a, b, c = idx // 25, (idx // 5) % 5, idx % 5
        return np.stack([2 - a, 2 - b, 2 - c], axis=1)

    def adjacency(self, tol=1e-12):
        A = np.abs(self.C) > tol
        np.fill_diagonal(A, False)
        return A


def components(adj):
    """Connected components of a boolean adjacency matrix (BFS)."""
    n = adj.shape[0]
    seen = np.zeros(n, dtype=bool)
    comps = []
    for s in range(n):
        if seen[s]:
            continue
        comp, queue = [], deque([s])
        seen[s] = True
        while queue:
            u = queue.popleft()
            comp.append(u)
            for v in np.nonzero(adj[u])[0]:
                if not seen[v]:
                    seen[v] = True
                    queue.append(v)
        comps.append(sorted(comp))
    return comps


def distances(adj):
    """All-pairs graph distances (np.inf where disconnected)."""
    n = adj.shape[0]
    D = np.full((n, n), np.inf)
    nbrs = [np.nonzero(adj[u])[0] for u in range(n)]
    for s in range(n):
        D[s, s] = 0
        queue = deque([s])
        while queue:
            u = queue.popleft()
            for v in nbrs[u]:
                if D[s, v] == np.inf:
                    D[s, v] = D[s, u] + 1
                    queue.append(v)
    return D
