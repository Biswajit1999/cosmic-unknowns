"""Validation response identities. Author: Biswajit Jana."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from validation_identification import drift_bounds, validation_regions
from calibrated_identification import physical_interval


class ValidationTests(unittest.TestCase):
    def test_interval_difference_bounds(self):
        rng = np.random.default_rng(205)
        for _ in range(100):
            r, v = rng.uniform(.1, .9, (2, 48))
            rl, ru = r-rng.uniform(0, .1, 48), r+rng.uniform(0, .1, 48)
            vl, vu = v-rng.uniform(0, .1, 48), v+rng.uniform(0, .1, 48)
            low, high = drift_bounds(rl, ru, vl, vu)
            self.assertLessEqual(low, np.max(abs(v-r)))
            self.assertGreaterEqual(high, np.max(abs(v-r)))

    def test_direct_region_nested_in_scalar_region(self):
        rng = np.random.default_rng(23)
        for _ in range(30):
            regions, _ = validation_regions(rng.binomial(1000, .4, 24), 1000,
                rng.binomial(500, .02, 24), rng.binomial(500, .02, 24), 500,
                rng.binomial(300, .03, 24), rng.binomial(300, .01, 24), 300)
            a, b = regions['direct_validation'], regions['scalar_transfer']
            self.assertTrue(np.all(a[0] >= b[0]-1e-10))
            self.assertTrue(np.all(a[1] <= b[1]+1e-10))

    def test_absent_validation_is_noninformative(self):
        regions, _ = validation_regions([50, 90], 100, [2, 3], [1, 3], 100, [0, 0], [0, 0], 0)
        for name in ['scalar_transfer', 'direct_validation', 'guarded_validation']:
            np.testing.assert_array_equal(regions[name][0], [0, 0])
            np.testing.assert_array_equal(regions[name][1], [1, 1])

    def test_guarded_region_contains_drifted_truth(self):
        rng = np.random.default_rng(4)
        t = rng.random(100); f = rng.uniform(.01, .04, 100); g = rng.uniform(.01, .04, 100)
        fs = f+rng.uniform(-.003, .003, 100); gs = g+rng.uniform(-.003, .003, 100)
        r = (1-fs)*t+gs*(1-t)
        lo, hi = physical_interval(r, r, f-.003, f+.003, g-.003, g+.003)
        self.assertTrue(np.all((lo <= t) & (t <= hi)))

    def test_known_nonblind_response_removes_exchange(self):
        t = np.linspace(.01, .99, 24); u = np.linspace(-1, 1, 24)
        f, g = .02-.0067*u, .02+.0067*u
        r = (1-f)*t+g*(1-t)
        np.testing.assert_allclose((r-g)/(1-f-g), t, atol=2e-16)

    def test_invalid_residual(self):
        with self.assertRaises(ValueError):
            validation_regions([5], 10, [1], [1], 10, [1], [1], 10, eta=-.01)


if __name__ == '__main__':
    unittest.main()
