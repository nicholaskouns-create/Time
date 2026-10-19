"""Finite-dimensional quantum time and E47-preserving graph dynamics."""
import numpy as np

def commutator(a, b):
    return a @ b - b @ a

def lindblad_rhs(rho, hamiltonian, jumps=()):
    """Deterministic density-matrix derivative, with (rate, operator) jumps."""
    out = -1j * commutator(hamiltonian, rho)
    for rate, op in jumps:
        adj = op.conj().T
        gram = adj @ op
        out += rate * (op @ rho @ adj - 0.5 * (gram @ rho + rho @ gram))
    return out

def projected_hamiltonian(h, projector):
    """Remove off-code couplings while preserving within-code dynamics."""
    return projector @ h @ projector

def rk4_step(rho, t, dt, h_of_t, jumps_of_t=lambda t: ()):
    """Fourth-order deterministic integration; positivity is not guaranteed at finite dt."""
    f = lambda s, x: lindblad_rhs(x, h_of_t(s), jumps_of_t(s))
    k1 = f(t, rho)
    k2 = f(t + dt/2, rho + dt*k1/2)
    k3 = f(t + dt/2, rho + dt*k2/2)
    k4 = f(t + dt, rho + dt*k3)
    return rho + dt*(k1 + 2*k2 + 2*k3 + k4)/6
