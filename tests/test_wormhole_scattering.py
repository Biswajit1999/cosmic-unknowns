"""Independent analytical and numerical benchmarks. Author: Biswajit Jana."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
from wormhole_scattering import scatter, potential, ground_energy, halfline_reflection


class ScatteringTests(unittest.TestCase):
    def test_transparent_critical_coupling(self):
        for k in (0.02, 0.3, 2.0):
            s = scatter(k, xi=0.5, domain=50)
            self.assertLess(abs(s.transmission-1), 2e-8)
            self.assertLess(s.R, 1e-16)

    def test_flux_balance(self):
        for xi in (0, 1/6, 0.45):
            s = scatter(0.3, xi=xi, domain=80)
            self.assertLess(s.flux_error, 2e-8)
            self.assertGreaterEqual(s.T, 0)
            self.assertLessEqual(s.T, 1+2e-8)

    def test_published_minimal_low_energy_limit(self):
        s = scatter(0.003, domain=320)
        coefficient = 4/np.pi**2
        self.assertLess(abs(s.T/s.k**2/coefficient-1), 0.002)

    def test_analytic_zero_energy_solution(self):
        # u=sqrt(1+X^2)*(pi/2-atan X), analytically u''=u/(1+X^2)^2.
        x = np.linspace(-5, 5, 101)
        h = 2e-4
        def u(z):
            return np.sqrt(1+z*z)*(np.pi/2-np.arctan(z))
        second = (u(x+h)-2*u(x)+u(x-h))/h**2
        np.testing.assert_allclose(second, potential(x)*u(x), rtol=4e-4, atol=2e-7)

    def test_bound_state_free_box(self):
        L = 20
        energy, _, _ = ground_energy(0.5, domain=L, spacing=0.025)
        self.assertLess(abs(energy/(np.pi/(2*L))**2-1), 2e-6)

    def test_massless_instability_above_threshold(self):
        energy, _, _ = ground_energy(0.6, domain=100, spacing=0.04)
        self.assertLess(energy, 0)

    def test_one_sided_response_equivalence(self):
        for k in (0.07, 0.4, 1.5):
            s = scatter(k, domain=60)
            r, _ = halfline_reflection(k, domain=60)
            self.assertLess(abs(r-s.reflection), 2e-8)

    def test_scattering_domain_convergence(self):
        a = scatter(0.05, domain=80)
        b = scatter(0.05, domain=160)
        self.assertLess(abs(a.T-b.T)/b.T, 5e-4)


if __name__ == "__main__":
    unittest.main()
