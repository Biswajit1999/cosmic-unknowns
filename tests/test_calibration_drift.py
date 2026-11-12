"""Calibration drift identities and confidence-region containment."""
import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from calibration_drift import drift_intervals,signed_drift
from calibrated_identification import joint_intervals,physical_interval

class DriftTests(unittest.TestCase):
    def test_zero_drift_recovers_previous_intervals(self):
        y=np.array([1,40,99]);fn=np.array([0,2,3]);fp=np.array([1,0,2])
        for a,b in zip(drift_intervals(y,100,fn,fp,200,0),joint_intervals(y,100,fn,fp,200)):
            np.testing.assert_array_equal(a,b)

    def test_nested_allowances(self):
        y=np.array([1,40,99]);fn=np.array([0,2,3]);fp=np.array([1,0,2])
        lo,hi=drift_intervals(y,100,fn,fp,200,.01)
        lo2,hi2=drift_intervals(y,100,fn,fp,200,.03)
        self.assertTrue(np.all(lo2<=lo));self.assertTrue(np.all(hi2>=hi))

    def test_signed_drift_identity(self):
        rng=np.random.default_rng(1);q=rng.random(50);u=rng.uniform(-1,1,50)
        for sign in [-1,1]:
            fn,fp=signed_drift(.02,u,.015,sign)
            np.testing.assert_allclose((1-fn)*q+fp*(1-q),.02+.96*q+sign*.015*u,atol=2e-16)

    def test_exact_exchange_construction(self):
        q=np.array([.1,.4,.9]);p=q+np.array([.002,-.004,.001]);gap=max(abs(p-q));u=(p-q)/gap
        fn,fp=signed_drift(.02,u,.96*gap)
        np.testing.assert_allclose((1-fn)*q+fp*(1-q),.02+.96*p,atol=1e-16)

    def test_bounded_response_rectangle_contains_true_model(self):
        rng=np.random.default_rng(9);q=rng.random(80);f=rng.uniform(.02,.05,80);g=rng.uniform(.02,.05,80)
        fs=f+rng.uniform(-.01,.01,80);gs=g+rng.uniform(-.01,.01,80);r=(1-fs)*q+gs*(1-q)
        lo,hi=physical_interval(r,r,f-.01,f+.01,g-.01,g+.01)
        self.assertTrue(np.all((lo<=q)&(q<=hi)))

    def test_two_models_share_midpoint_at_half_exchange(self):
        q=np.array([.1,.4,.9]);p=q+np.array([.002,-.004,.001])
        gap=max(abs(p-q));u=(p-q)/gap;eps=.96*gap/2
        fq,gq=signed_drift(.02,u,eps,1)
        fp,gp=signed_drift(.02,u,eps,-1)
        rq=(1-fq)*q+gq*(1-q);rp=(1-fp)*p+gp*(1-p)
        np.testing.assert_allclose(rq,rp,atol=2e-16)
        np.testing.assert_allclose(rq,.02+.96*(p+q)/2,atol=2e-16)
        # Below this allowance at least one channel's response intervals are disjoint.
        smaller=.999*eps
        self.assertTrue(np.any(.96*abs(p-q)>2*smaller))

    def test_invalid_allowance(self):
        with self.assertRaises(ValueError):drift_intervals([1],10,[0],[0],10,-.1)

if __name__=='__main__':unittest.main()
