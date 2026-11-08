"""Physical consistency checks on the educational baseline."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from wormhole_baseline import radius, stresses, STRESS_FACTOR


class WormholeBaselineTests(unittest.TestCase):
    def test_proper_distance_metric_and_throat(self):
        r0 = 3.0
        l = np.array([-4.0, -1.0, 0.0, 1.0, 4.0])
        r = radius(l, r0)
        self.assertEqual(float(radius(0, r0)), r0)
        # (dr/dl)^2 = 1-b/r, with b=r0^2/r; avoids dividing by zero at the throat.
        np.testing.assert_allclose((l/r)**2, 1-r0**2/r**2, atol=1e-14)

    def test_local_stress_conservation(self):
        r0 = 2.0
        r = np.linspace(2.1, 20, 100)
        step = 1e-4
        derivative = (stresses(r+step, r0)[1]-stresses(r-step, r0)[1])/(2*step)
        _, pr, pt = stresses(r, r0)
        np.testing.assert_allclose(derivative, 2*(pt-pr)/r, rtol=1e-7)

    def test_throat_nec_from_independent_flare_out_expression(self):
        r0 = 4.0
        h = 1e-5
        b_prime = ((r0*r0/(r0+h))-(r0*r0/(r0-h)))/(2*h)
        self.assertLess(b_prime, 1)
        energy, radial, _ = stresses(r0, r0)
        expected = STRESS_FACTOR*(b_prime-1)/r0**2
        self.assertAlmostEqual(float((energy+radial)/expected), 1, places=9)

    def test_throat_stress_scaling(self):
        self.assertAlmostEqual(float(stresses(20, 20)[0]/stresses(2, 2)[0]), 0.01)

    def test_invalid_domain(self):
        with self.assertRaises(ValueError):
            radius(1, 0)
        with self.assertRaises(ValueError):
            stresses(0.5, 1)


if __name__ == "__main__":
    unittest.main()
