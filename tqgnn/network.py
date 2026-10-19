"""Temporal quantum graph neural network (TQGNN) of E47 clock-memory nodes.

Each graph node v carries the reduced E47 clock-memory unit (logical (x) flag,
dimension 4) and follows its own worldline with lapse N_v(t).  The proper-time
generator is graph-structured in the QGNN sense,

    H_tau = sum_v H_v + sum_{(j,k) in E} [ J_jk (s+_j s-_k + s-_j s+_k) + Z_jk z_j z_k ],
    jump operators R_v (one recovery channel per node),

and the coordinate-time generator rescales every term by the lapse of the
worldline it lives on:

    H_t = sum_v N_v(t) H_v + sum_e N_e(t) H_e,      R_{v,t} = sqrt(N_v(t)) R_v,

with N_e = (N_j + N_k)/2 by default (a modelling choice; any positive edge
lapse keeps the generator GKSL).  With equal lapses this is exactly the
time-chain identity d rho/dt = N(t) L_tau rho.  One layer of the network is
the CPTP map over one coordinate-time step; stacked layers form the
recurrent QGNN unrolled in time.
"""
import numpy as np

from .lindblad import evolve
from .node import ClockMemoryNode
from .spacetime import Worldline

# logical operators inside one node (index 2*l + f)
_A_TO_B = np.kron(np.array([[0, 0], [1, 0]], dtype=complex), np.eye(2))   # |b><a| (x) I
_ZL = np.kron(np.diag([1.0, -1.0]).astype(complex), np.eye(2))           # |a><a| - |b><b|
_PROT = np.kron(np.eye(2), np.diag([1.0, 0.0])).astype(complex)          # flag 0
_COH = np.kron(np.array([[0, 1], [0, 0]], dtype=complex), np.eye(2))     # |a><b| (x) I


class TQGNN:
    def __init__(self, n_nodes, edges, omega=2.0, c=0.0, gamma=0.0, J=0.0, Jzz=0.0,
                 worldlines=None, clock="linear", edge_lapse="mean", node=None):
        self.n = n_nodes
        self.edges = [tuple(e) for e in edges]
        self.omega = omega
        self.gamma = gamma
        self.node = node or ClockMemoryNode()
        Hn, Rn = self.node.reduced(omega, c, gamma, clock)
        self.dim = 4 ** n_nodes
        self.worldlines = worldlines or [Worldline([0.0, 0.0, 0.0]) for _ in range(n_nodes)]
        Jd = J if isinstance(J, dict) else {e: J for e in self.edges}
        Zd = Jzz if isinstance(Jzz, dict) else {e: Jzz for e in self.edges}
        self.edge_lapse_rule = edge_lapse

        self.H_node = [self.embed(Hn, v) for v in range(n_nodes)]
        self.R_node = [self.embed(Rn, v) for v in range(n_nodes)] if gamma > 0 else []
        self.H_edge = {}
        for (j, k) in self.edges:
            sp_j, sp_k = self.embed(_A_TO_B, j), self.embed(_A_TO_B, k)
            hop = sp_j @ sp_k.conj().T
            Hjk = Jd.get((j, k), 0.0) * (hop + hop.conj().T)
            Hjk = Hjk + Zd.get((j, k), 0.0) * self.embed(_ZL, j) @ self.embed(_ZL, k)
            self.H_edge[(j, k)] = Hjk
        self.P_node = [self.embed(_PROT, v) for v in range(n_nodes)]
        self.coh_node = [self.embed(_COH, v) for v in range(n_nodes)]

    # -- structure -----------------------------------------------------------
    def embed(self, op, v):
        out = np.eye(1, dtype=complex)
        for w in range(self.n):
            out = np.kron(out, op if w == v else np.eye(4))
        return out

    def lapses(self, t):
        return np.array([w.N(t) for w in self.worldlines])

    def edge_lapse(self, Nj, Nk):
        if self.edge_lapse_rule == "mean":
            return 0.5 * (Nj + Nk)
        if self.edge_lapse_rule == "geometric":
            return np.sqrt(Nj * Nk)
        raise ValueError(self.edge_lapse_rule)

    def generator(self, t, lapses=None):
        N = self.lapses(t) if lapses is None else lapses
        H = sum(N[v] * self.H_node[v] for v in range(self.n))
        for (j, k), He in self.H_edge.items():
            H = H + self.edge_lapse(N[j], N[k]) * He
        Ls = [np.sqrt(N[v]) * self.R_node[v] for v in range(len(self.R_node))]
        return H, Ls

    def proper_generator(self):
        return self.generator(0.0, lapses=np.ones(self.n))

    # -- execution -------------------------------------------------------------
    def run_dynamic(self, rho0, t_eval, rtol=1e-11, atol=1e-13):
        """Deterministic (dynamic) execution: the GKSL flow in coordinate time."""
        return evolve(self.generator, rho0, t_eval, rtol=rtol, atol=atol)

    def readout(self, rho, t=0.0):
        """Per-node retention Tr(P_v rho), clock coherence Tr(|a><b|_v rho) and
        coordinate-time click intensity N_v gamma Tr((1 - P_v) rho)."""
        N = self.lapses(t)
        ret = np.array([np.trace(P @ rho).real for P in self.P_node])
        coh = np.array([np.trace(Cm @ rho) for Cm in self.coh_node])
        clicks = N * self.gamma * (1 - ret)
        return ret, coh, clicks

    # -- states ------------------------------------------------------------------
    @staticmethod
    def product_state(node_vectors):
        psi = np.ones(1, dtype=complex)
        for v in node_vectors:
            psi = np.kron(psi, np.asarray(v, dtype=complex))
        return psi / np.linalg.norm(psi)

    @staticmethod
    def clock_state(theta=np.pi / 2, phase=0.0, flag=0):
        """cos(theta/2)|a> + e^{i phase} sin(theta/2)|b>, flag fixed."""
        logical = np.array([np.cos(theta / 2), np.exp(1j * phase) * np.sin(theta / 2)])
        f = np.zeros(2); f[flag] = 1
        return np.kron(logical, f)
