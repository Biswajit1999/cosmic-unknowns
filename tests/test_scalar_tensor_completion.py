"""Full background equations and kinetic sign. Author: Biswajit Jana."""
import sys
import unittest
from pathlib import Path
from mpmath import mp
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from scalar_tensor_completion import couplings, check_equations, instability_mode


class ScalarTensorTests(unittest.TestCase):
    def test_full_equations(self):
        with mp.workdps(60):
            for b, x in [('.2', '-.3'), ('-.2', '2'), ('0', '0')]:
                metric, scalar, kinetic = check_equations(mp.mpf(x), mp.mpf(b))
                self.assertLess(max(abs(r) for r in metric), mp.mpf('1e-50'))
                self.assertLess(abs(scalar), mp.mpf('1e-50'))
                self.assertLess(abs(kinetic-2), mp.mpf('1e-50'))

    def test_ellis_action_limit(self):
        f, h, z = couplings(mp.mpf('.4'), mp.mpf(0))
        self.assertEqual((f, h, z), (1, 1, 2))

    def test_singular_conformal_map_excluded(self):
        with self.assertRaises(ValueError): couplings(mp.mpf(-3), mp.mpf(1))

    def test_negative_mode_and_variational_bound(self):
        x, mode, row = instability_mode(20., .05)
        self.assertLess(row['eigenvalue'], -1.)
        self.assertLess(row['eigen_residual'], 1e-9)
        with mp.workdps(35):
            energy = mp.quad(lambda y: (y*y-3)/(1+y*y)**3, [-mp.inf, mp.inf])
            norm = mp.quad(lambda y: 1/(1+y*y), [-mp.inf, mp.inf])
            self.assertLess(abs(energy/norm+1), mp.mpf('1e-30'))


if __name__ == '__main__': unittest.main()
