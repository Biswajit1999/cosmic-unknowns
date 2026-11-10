"""Algebraic equivalence and interval checks. Author: Biswajit Jana."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
from calibration_limits import required_error, matching_errors, recorded_probability, clipped_probability, bernoulli_kl


class CalibrationTests(unittest.TestCase):
    def test_constructed_equivalence_random_probabilities(self):
        rng = np.random.default_rng(394)
        p, q = rng.uniform(.0001, .9999, (2, 1000))
        fp, fn = matching_errors(p, q)
        np.testing.assert_allclose(recorded_probability(q, fp, fn), p, atol=2e-16)
        np.testing.assert_allclose(np.maximum(fp, fn), required_error(p, q))

    def test_bound_is_necessary_and_sufficient(self):
        p = np.array([.2, .7])
        q = np.array([.3, .6])
        delta = float(max(required_error(p, q)))
        np.testing.assert_allclose(clipped_probability(p, q, delta), p)
        self.assertGreater(float(sum(bernoulli_kl(p, clipped_probability(p, q, delta*.99)))), 0)

    def test_no_noise_and_equal_models(self):
        p, q = np.array([.1, .8]), np.array([.2, .6])
        np.testing.assert_array_equal(clipped_probability(p, q, 0), q)
        np.testing.assert_array_equal(required_error(p, p), np.zeros(2))

    def test_extreme_targets(self):
        np.testing.assert_array_equal(required_error(np.array([0., 1.]), np.array([.3, .8])), np.ones(2))


if __name__ == "__main__":
    unittest.main()
