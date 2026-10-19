"""Clock readouts: sector phase, recovery clicks, prime cycles."""
import numpy as np


def phase_time(coherence, omega):
    """Proper time read from the clock coherence Tr(|a><b| rho) ~ e^{-i omega tau}/2."""
    return np.unwrap(-np.angle(np.asarray(coherence))) / omega


def click_time(n_clicks, rate):
    """Proper time estimated by counting clicks of known mean rate."""
    return n_clicks / rate


def click_time_precision(fano, rate, tau):
    """Asymptotic relative standard deviation of the click-time estimate."""
    return np.sqrt(fano / (rate * tau))


def cycle_unitary(basis):
    """Cyclic shift |e_k> -> |e_{k+1 mod n}> on the columns of an orthonormal basis,
    identity on the orthogonal complement."""
    n = basis.shape[1]
    S = np.roll(np.eye(n), 1, axis=0)
    d = basis.shape[0]
    return basis @ S @ basis.conj().T + (np.eye(d) - basis @ basis.conj().T)


def matrix_order(U, kmax, tol=1e-9):
    """Smallest k <= kmax with U^k = I, else None."""
    X = np.eye(U.shape[0], dtype=complex)
    for k in range(1, kmax + 1):
        X = X @ U
        if np.linalg.norm(X - np.eye(U.shape[0])) < tol:
            return k
    return None
