"""Independent geometric and reduced-energy checks. Author: Biswajit Jana."""
import sys
import unittest
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from magnetic_energy import MagneticEnergy, field_coordinates, field_length, scale_diagnostics


class MagneticEnergyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.model = MagneticEnergy(128, 256)

    def test_field_line_flux_label_and_endpoints(self):
        for n in [.01, .3, 1., 1.9]:
            eta = np.linspace(-.99, .99, 103)
            z, r = field_coordinates(eta, n)
            flux = (z+.5)/np.sqrt((z+.5)**2+r*r)-(z-.5)/np.sqrt((z-.5)**2+r*r)
            np.testing.assert_allclose(flux, n, atol=2e-13, rtol=2e-13)
            z, r = field_coordinates(np.array([-1., 1.]), n)
            np.testing.assert_allclose(z, [-.5, .5], atol=1e-15)
            np.testing.assert_allclose(r, 0, atol=1e-15)

    def test_length_from_independent_polygonal_curve(self):
        eta = np.linspace(-1., 1., 20001)
        for n in [.1, .7, 1.9]:
            z, r = field_coordinates(eta, n)
            polygon = np.hypot(np.diff(z), np.diff(r)).sum()
            self.assertLess(abs(polygon/float(field_length(n))-1), 3e-8)

    def test_shortest_path_bound_and_straight_limit(self):
        self.assertEqual(float(field_length(2)), 1.)
        self.assertTrue(np.all(field_length(np.geomspace(1e-5, 1.999, 100)) > 1))
        constant = 2+np.arcsinh(np.sqrt(3))/np.sqrt(3)
        self.assertLess(abs(float(1e-7*field_length(1e-7))-constant), 2e-7)

    def test_independent_energy_minimizer(self):
        for s in [.01, 1., 100.]:
            row = self.model.minimum(s)
            opt = minimize_scalar(lambda logy: self.model.energy(np.exp(logy), s),
                                  bounds=(np.log(.5), np.log(1e4)), method='bounded',
                                  options={'xatol': 1e-12})
            self.assertTrue(opt.success)
            self.assertLess(abs(np.exp(opt.x)/row['y']-1), 1e-7)
            self.assertGreater(row['y'], max(1., 1/s))
            self.assertLess(row['energy'], 0)
            self.assertGreater(row['curvature'], 0)

    def test_energy_derivative_and_h_monotonicity(self):
        for s in [.1, 10.]:
            for y in [1., 3., 20.]:
                step = y*1e-4
                derivative = (self.model.energy(y+step, s)-self.model.energy(y-step, s))/(2*step)
                self.assertLess(abs(derivative-self.model.stationarity(y,s)/y**2), 5e-7)
                hp = (self.model.stationarity(y+step,s)-self.model.stationarity(y-step,s))/(2*step)
                self.assertLess(abs(hp/self.model.hprime(y,s)-1), 3e-8)

    def test_tail_omission_increases_length(self):
        cut = MagneticEnergy(128, 256, cutoff=.1)
        for s in [.01, 1., 100.]:
            base = self.model.minimum(s)['y']; shifted = cut.minimum(s)['y']
            self.assertGreater(shifted, base)
            self.assertLess(self.model.stationarity(shifted, s), .05+1e-10)

    def test_short_separation_scale_identity(self):
        s = 1e-9; q = 1e6; g = .3
        # Leading analytic y=1/s, e=-3s/8 gives |E|ell=q/16.
        scales = scale_diagnostics(q, g, s, 1/s, -3*s/8)
        self.assertAlmostEqual(scales['binding_times_ell'], q/16)
        self.assertAlmostEqual(scales['schwarzian_coefficient'], q/16)

    def test_invalid_domains(self):
        for n in [0, -1, 2.1]:
            with self.assertRaises(ValueError): field_length(n)
        with self.assertRaises(ValueError): self.model.minimum(0)
        with self.assertRaises(ValueError): MagneticEnergy(cutoff=1.5)


if __name__ == '__main__': unittest.main()
