"""Independent analytical and statistical checks. Author: Biswajit Jana."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
from transmission_inference import (wilson, pt_transmission, pt_parameters,
    numerical_pt, physical_parameters, deviance, fit_pt)


class InferenceTests(unittest.TestCase):
    def test_wilson_boundary_samples(self):
        low, high = wilson(0, 10)
        self.assertAlmostEqual(low, 0)
        self.assertAlmostEqual(high, 0.2775327998628892)
        a, b = wilson(10, 10)
        self.assertAlmostEqual(a, 1-high)
        self.assertAlmostEqual(b, 1)

    def test_pt_independent_ode_both_branches(self):
        for height, alpha in [(0.1, 1.), (1., 1.), (0.25, 1.)]:
            theta = pt_parameters(height, alpha)
            np.testing.assert_allclose(physical_parameters(theta), [height, alpha], atol=1e-13)
            for k in [0.1, 0.7, 1.5]:
                self.assertLess(abs(pt_transmission(k, theta)-numerical_pt(k, height, alpha)), 2e-9)

    def test_noise_free_fit(self):
        k = np.geomspace(0.05, 1.5, 24)
        expected = pt_parameters(0.8, 1.2)
        theta, _, bound, _ = fit_pt(k, pt_transmission(k, expected), multistart=True)
        np.testing.assert_allclose(theta, expected, atol=2e-5)
        self.assertFalse(bound)

    def test_deviance_edge_counts(self):
        self.assertAlmostEqual(deviance([0, 10], 10, [0.5, 0.5]), 40*np.log(2))
        self.assertAlmostEqual(deviance([2, 5], 10, [0.2, 0.5]), 0)

    def test_low_frequency_limit(self):
        theta = pt_parameters(1., 1.)
        c, b = np.exp(theta)
        k = 1e-5
        self.assertAlmostEqual(pt_transmission(k, theta)/k**2, b*b/c, places=8)


if __name__ == "__main__":
    unittest.main()
