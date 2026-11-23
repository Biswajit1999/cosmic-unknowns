"""Independent limiting cases, flux and reciprocity. Author: Biswajit Jana."""
import sys
import unittest
from pathlib import Path
import numpy as np
from mpmath import mp
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from conformal_scattering import scatter, stress
from massive_conformal import scatter_massive, tail_data, tail, scalar_source_obstruction


class MassiveTests(unittest.TestCase):
    def test_zero_mass_reproduces_previous_solver(self):
        for ell in [0, 1]:
            a = scatter_massive(.8, .2, 0., ell)
            b = scatter(.8, ell, .2, form='reduced')
            self.assertLess(abs(a['t']-b['t']), 3e-8)
            self.assertLess(abs(a['r']-b['r']), 3e-8)

    def test_ultrastatic_mass_is_frequency_shift(self):
        a = scatter_massive(.7, 0., .3, 1)
        b = scatter(np.sqrt(.7**2-.3**2), 1, 0., form='reduced')
        self.assertLess(abs(a['t']-b['t']), 3e-8)
        self.assertLess(abs(a['r']-b['r']), 3e-8)

    def test_reciprocity_and_velocity_weighted_flux(self):
        a = scatter_massive(.6, -.2, .3)
        b = scatter_massive(.6, -.2, .3, incidence='left')
        self.assertLess(abs(a['t_flux']-b['t_flux']), 3e-8)
        self.assertLess(a['flux_error'], 3e-8)
        self.assertLess(b['flux_error'], 3e-8)
        self.assertGreater(abs(abs(a['t'])**2-a['T']), .005)

    def test_original_equation_and_convergence(self):
        a = scatter_massive(.54, -.2, .3, 1)
        b = scatter_massive(.54, -.2, .3, 1, form='radial', domain=2*a['domain'], order=16)
        self.assertLess(abs(a['r']-b['r']), 3e-8)
        self.assertLess(abs(a['t_flux']-b['t_flux']), 3e-8)

    def test_coulomb_leading_coefficient_and_current(self):
        for side in [-1, 1]:
            q, eta, v, a = tail_data(.7, .2, .3, 1, side)
            f0 = 1 if side == 1 else 1-np.pi*.2
            self.assertAlmostEqual(v[1], -side*2*.3**2*f0*.2)
            self.assertAlmostEqual(eta, v[1]/(2*q))
            u, du, _, _ = tail(500., .7, .2, .3, 1, side)
            self.assertLess(abs(np.imag(u.conjugate()*du)/q-1), 1e-12)

    def test_closed_channel_rejected(self):
        with self.assertRaises(ValueError):
            scatter_massive(.4, -.2, .3)

    def test_single_scalar_source_obstruction(self):
        x = np.linspace(-8, 8, 101)
        for beta in [-.2, 0., .2]:
            rho, pr, pt = stress(x, beta)
            np.testing.assert_allclose(rho+pt, scalar_source_obstruction(x, beta), atol=2e-14)
        self.assertLess(scalar_source_obstruction(1., .2), 0)
        self.assertGreater(scalar_source_obstruction(-1., .2), 0)


if __name__ == '__main__':
    unittest.main()
