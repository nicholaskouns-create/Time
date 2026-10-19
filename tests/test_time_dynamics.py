import unittest
import numpy as np
from src.time_dynamics import lindblad_rhs, projected_hamiltonian, rk4_step

class QuantumTimeTests(unittest.TestCase):
    def test_unitary_derivative(self):
        h = np.array([[0,1],[1,0]], dtype=complex)
        rho = np.diag([1.,0.]).astype(complex)
        d = lindblad_rhs(rho,h)
        self.assertAlmostEqual(np.trace(d).real,0.)
        self.assertLess(np.linalg.norm(d-d.conj().T),1e-12)

    def test_projected_invariance(self):
        p = np.diag([1.,0.]).astype(complex)
        h = np.array([[1.,2.],[2.,3.]], dtype=complex)
        hp = projected_hamiltonian(h,p)
        rho = p.copy()
        result = rk4_step(rho,0,0.01,lambda t:hp)
        self.assertLess(np.linalg.norm(result-p@result@p),1e-12)

if __name__ == "__main__":
    unittest.main()
