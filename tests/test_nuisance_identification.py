"""Symmetry and nuisance-domain tests. Author: Biswajit Jana."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from massive_conformal import potential, scatter_massive
from nuisance_identification import centered_parameters, reversed_parameters, NullSpectrum


class NuisanceTests(unittest.TestCase):
    def test_reversed_potential_identity(self):
        x = np.linspace(-15, 15, 301)
        for b in [-.2, .2]:
            br, mr = reversed_parameters(b, .3)
            for ell in [0, 1]:
                np.testing.assert_allclose(potential(x, b, .3, ell), potential(-x, br, mr, ell), atol=3e-15)

    def test_centered_parity_and_flux(self):
        bp, mp = centered_parameters(.08, .3)
        bm, mm = centered_parameters(-.08, .3)
        x = np.linspace(-8, 8, 101)
        np.testing.assert_allclose(potential(x, bp, mp), potential(-x, bm, mm), atol=2e-15)
        a = scatter_massive(.7, bp, mp)
        b = scatter_massive(.7, bm, mm)
        self.assertLess(abs(a['T']-b['T']), 2e-9)

    def test_derivative_confounding(self):
        h = .0005; nu = .3; k = .8
        db = (scatter_massive(k, h, nu)['T']-scatter_massive(k, -h, nu)['T'])/(2*h)
        dm = (scatter_massive(k, 0., nu+h)['T']-scatter_massive(k, 0., nu-h)['T'])/(2*h)
        self.assertLess(abs(db+np.pi*nu*dm/2), 2e-6)

    def test_no_surrogate_extrapolation(self):
        rows = [dict(q=x, T=x/3, ell=l) for l in [0, 1] for x in [.2, 1., 2.2]]
        model = NullSpectrum(rows)
        with self.assertRaises(ValueError):
            model.predict([3.], [0], [1., .3, 0., 0.])
        value = model.predict([1.], [0], [1., .3, .01, .02])[0]
        self.assertAlmostEqual(value, .01+.97*np.sqrt(1-.09)/3)

    def test_end_exchange_is_involution(self):
        b, m = reversed_parameters(*reversed_parameters(.2, .3))
        self.assertAlmostEqual(b, .2)
        self.assertAlmostEqual(m, .3)


if __name__ == '__main__':
    unittest.main()
