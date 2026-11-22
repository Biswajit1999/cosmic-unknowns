"""Geometry and independent scalar-equation checks. Author: Biswajit Jana."""
import sys
import unittest
from pathlib import Path
import numpy as np
from mpmath import mp
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from conformal_scattering import geometry, metric_potential, common_potential, stress, throat, scatter


class ConformalTests(unittest.TestCase):
    def test_potential_identity_all_sampled_modes_and_couplings(self):
        x = np.linspace(-30, 30, 301)
        for beta in [-.2, 0., .2, .3]:
            for ell in [0, 1, 2, 5]:
                for xi in [0., 1/6, .4]:
                    np.testing.assert_allclose(metric_potential(x, beta, ell, xi), common_potential(x, ell, xi), atol=3e-14)

    def test_independent_radius_derivatives(self):
        with mp.workdps(50):
            for beta in [-.2, .2]:
                radius = lambda x: (1+mp.mpf(beta)*(mp.atan(x)-mp.pi/2))*mp.sqrt(1+x*x)
                for x in [-4., -.3, 0., 2.]:
                    values = geometry(x, beta)
                    self.assertAlmostEqual(values[4], float(mp.diff(radius, x)), places=13)
                    self.assertAlmostEqual(values[5], float(mp.diff(radius, x, 2)), places=13)

    def test_regular_throat_and_nec(self):
        for beta in [-.2, 0., .2, .3]:
            t = throat(beta)
            self.assertGreater(t['radius'], 0)
            self.assertGreater(t['flare_out'], 0)
            self.assertLess(t['radial_nec'], 0)
            self.assertAlmostEqual(float(geometry(t['X'], beta)[4]), 0, places=12)

    def test_einstein_trace_matches_curvature(self):
        x = np.linspace(-10, 10, 101)
        for beta in [-.2, 0, .2]:
            f = geometry(x, beta)[0]
            rho, pr, pt = stress(x, beta)
            np.testing.assert_allclose(rho-pr-2*pt, -2/(f*f*(1+x*x)**2), atol=2e-14)

    def test_original_field_and_reduced_scattering_agree(self):
        for ell in [0, 1, 2]:
            a = scatter(.8, ell, .2, form='radial')
            b = scatter(.8, ell, 0, form='reduced')
            self.assertLess(abs(a['r']-b['r']), 2e-8)
            self.assertLess(abs(a['t']-b['t']), 2e-8)
            self.assertLess(a['flux_error'], 2e-8)

    def test_transparent_exact_limit(self):
        a = scatter(.7, 0, .2, xi=.5)
        self.assertLess(abs(a['r']), 2e-8)
        self.assertLess(abs(a['t']-1), 2e-8)

    def test_excluded_degenerate_endpoint(self):
        with self.assertRaises(ValueError):
            geometry(0., 1/np.pi)


if __name__ == '__main__':
    unittest.main()
