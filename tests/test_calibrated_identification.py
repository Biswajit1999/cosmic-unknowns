"""Tests of coverage construction and continuum exclusion. Author: Biswajit Jana."""
import sys
from pathlib import Path
import unittest
import numpy as np
from scipy.stats import binom
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from calibrated_identification import cp_interval,physical_interval,continuum_test,joint_intervals
from transmission_inference import pt_transmission

class CalibrationTests(unittest.TestCase):
    def test_exact_small_n_coverage(self):
        n=20; y=np.arange(n+1); lo,hi=cp_interval(y,n,.025)
        for p in np.linspace(0,1,101):
            coverage=np.sum(binom.pmf(y,n,p)*((lo<=p)&(p<=hi)))
            self.assertGreaterEqual(coverage, .95-1e-12)

    def test_zero_errors_bound(self):
        lo,hi=cp_interval(np.array([0]),100,.001)
        self.assertEqual(lo[0],0)
        self.assertAlmostEqual(hi[0],1-.001**.01,places=10)

    def test_response_projection_contains_truth(self):
        rng=np.random.default_rng(12)
        q,f,g=rng.random((3,100))
        p=(1-f)*q+g*(1-q)
        lo,hi=physical_interval(np.maximum(0,p-.01),np.minimum(1,p+.01),
            np.maximum(0,f-.03),np.minimum(1,f+.03),np.maximum(0,g-.03),np.minimum(1,g+.03))
        self.assertTrue(np.all((lo<=q)&(q<=hi)))

    def test_known_comparator_never_excluded(self):
        k=np.geomspace(.05,1.5,24)
        for theta in [[2.5,.94],[-15,np.log(.11)],[24,np.log(29)]]:
            p=pt_transmission(k,theta)
            result=continuum_test(k,np.maximum(0,p-1e-4),np.minimum(1,p+1e-4))
            self.assertNotEqual(result['status'],'excluded')

    def test_decreasing_spectrum_excluded(self):
        k=np.geomspace(.05,1.5,24);p=np.linspace(.8,.2,24)
        r=continuum_test(k,p-.001,p+.001,return_certificate=True)
        self.assertEqual(r['status'],'excluded')
        intervals=sorted((v[0],v[1]) for v in r['leaves'])
        self.assertEqual(intervals[0][0],.1);self.assertEqual(intervals[-1][1],30.)
        for a,b in zip(intervals,intervals[1:]):self.assertEqual(a[1],b[0])

    def test_computational_limit_is_not_rejection(self):
        r=continuum_test(np.array([.1,.5]),np.array([.3,.2]),np.array([.31,.21]),max_nodes=0)
        self.assertEqual(r['status'],'undetermined')

    def test_no_calibration_leaves_full_physical_interval(self):
        lo,hi=joint_intervals(np.array([0,50,100]),100,np.zeros(3,int),np.zeros(3,int),0)
        np.testing.assert_array_equal(lo,np.zeros(3))
        np.testing.assert_array_equal(hi,np.ones(3))

    def test_known_blind_response_does_not_identify_q(self):
        lo,hi=physical_interval(np.array([.6]),np.array([.6]),.4,.4,.6,.6)
        self.assertEqual(lo[0],0);self.assertEqual(hi[0],1)

if __name__=='__main__':unittest.main()
