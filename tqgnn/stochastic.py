"""Stochastic execution: unravelling the recovery channels with detection efficiency eta.

Conditional state (Wiseman-Milburn jump unravelling, efficiency eta):

    between detected clicks   d rho~/dt = L_t rho~ - eta sum_k J_k rho~        (linear, trace decays)
    at a detected click in k  rho -> J_k rho / Tr(J_k rho),   J_k rho = R_k rho R_k^dag

Tr rho~ is the probability of no detected click, so clicks are drawn by the
waiting-time method.  For every eta, E[rho_c(t)] = rho(t) (the GKSL solution).
At eta = 0 nothing is recorded and rho_c(t) = rho(t) exactly: the stochastic
system has become the deterministic (dynamic) one.
"""
import numpy as np


def _rk4(f, y, t, dt):
    k1 = f(t, y)
    k2 = f(t + dt / 2, y + dt / 2 * k1)
    k3 = f(t + dt / 2, y + dt / 2 * k2)
    k4 = f(t + dt, y + dt * k3)
    return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def conditional_density(generator, rho0, t_grid, n_traj, eta, dt, rng):
    """Batched conditional density matrices.  Returns (states[t, traj], clicks[traj])."""
    d = rho0.shape[0]
    M = n_traj
    rho = np.repeat(np.asarray(rho0, dtype=complex)[None], M, axis=0)
    r = rng.random(M)
    clicks = [[] for _ in range(M)]
    out = np.empty((len(t_grid), M, d, d), dtype=complex)
    out[0] = rho

    def f(t, R):
        H, Ls = generator(t)
        G = sum((L.conj().T @ L for L in Ls), np.zeros((d, d), complex))
        X = (H - 0.5j * G) @ R
        res = -1j * (X - np.conj(np.swapaxes(X, 1, 2)))
        for L in Ls:
            res = res + (1 - eta) * (L @ R @ L.conj().T)
        return res

    t = t_grid[0]
    for i in range(1, len(t_grid)):
        n_sub = max(1, int(round((t_grid[i] - t_grid[i - 1]) / dt)))
        h = (t_grid[i] - t_grid[i - 1]) / n_sub
        for _ in range(n_sub):
            rho = _rk4(f, rho, t, h)
            t += h
            tr = np.einsum("mii->m", rho).real
            hit = np.nonzero(tr < r)[0]
            if hit.size:
                _, Ls = generator(t)
                for m in hit:
                    w = np.array([np.trace(L @ rho[m] @ L.conj().T).real for L in Ls]) * eta
                    k = rng.choice(len(Ls), p=w / w.sum())
                    new = Ls[k] @ rho[m] @ Ls[k].conj().T
                    rho[m] = new / np.trace(new).real
                    clicks[m].append((t, k))
                    r[m] = rng.random()
        tr = np.einsum("mii->m", rho).real
        out[i] = rho / tr[:, None, None]
    return out, clicks


def jump_trajectories(generator, psi0, t_grid, n_traj, dt, rng):
    """eta = 1 pure-state trajectories (batched).  Returns (states[t, traj], clicks[traj])."""
    d = psi0.shape[0]
    M = n_traj
    psi = np.repeat(np.asarray(psi0, dtype=complex)[None], M, axis=0)
    r = rng.random(M)
    clicks = [[] for _ in range(M)]
    out = np.empty((len(t_grid), M, d), dtype=complex)
    out[0] = psi

    def f(t, Y):
        H, Ls = generator(t)
        G = sum((L.conj().T @ L for L in Ls), np.zeros((d, d), complex))
        return Y @ (-1j * (H - 0.5j * G)).T

    t = t_grid[0]
    for i in range(1, len(t_grid)):
        n_sub = max(1, int(round((t_grid[i] - t_grid[i - 1]) / dt)))
        h = (t_grid[i] - t_grid[i - 1]) / n_sub
        for _ in range(n_sub):
            psi = _rk4(f, psi, t, h)
            t += h
            nrm = np.einsum("mi,mi->m", psi.conj(), psi).real
            hit = np.nonzero(nrm < r)[0]
            if hit.size:
                _, Ls = generator(t)
                for m in hit:
                    amps = [L @ psi[m] for L in Ls]
                    w = np.array([np.vdot(a, a).real for a in amps])
                    k = rng.choice(len(Ls), p=w / w.sum())
                    psi[m] = amps[k] / np.sqrt(w[k])
                    clicks[m].append((t, k))
                    r[m] = rng.random()
        nrm = np.einsum("mi,mi->m", psi.conj(), psi).real
        out[i] = psi / np.sqrt(nrm)[:, None]
    return out, clicks


def thin(clicks, eta, rng):
    """Keep each recorded click independently with probability eta."""
    return [[c for c in rec if rng.random() < eta] for rec in clicks]
