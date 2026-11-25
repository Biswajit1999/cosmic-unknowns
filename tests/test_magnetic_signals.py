"""Independent finite-signal checks. Author: Biswajit Jana."""
import sys
import unittest
from pathlib import Path
import numpy as np
from scipy.integrate import quad
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from magnetic_signals import (amplitudes,matching_solve,local_transmission,mouth_probability,
    halfwidth,eventual_fraction,collected_fraction,collection_pulses,gaussian_average,
    gaussian_fourier,pulse_amplitude,compact_envelope,compact_autocorrelation)


class MagneticSignalTests(unittest.TestCase):
    def test_matching_and_flux(self):
        for x in [.005,.02,.05]:
            for z in [.2,.5,1.,1.001,1.1]:
                r,t=amplitudes(x,z); rr,tt=matching_solve(x,z)
                self.assertLess(max(abs(r-rr),abs(t-tt)),2e-12)
                self.assertAlmostEqual(abs(r)**2+abs(t)**2,1.,places=12)

    def test_narrow_resonance_and_halfwidth(self):
        for x in [1e-6,.005,.05]:
            self.assertAlmostEqual(local_transmission(x,0),-1,places=14)
            self.assertAlmostEqual(abs(local_transmission(x,halfwidth(x)))**2,.5,places=12)

    def test_geometric_series_reconstructs_amplitude(self):
        x=.1;p=mouth_probability(x);r=1-p
        for z in [.01,.2,.5]:
            m=np.arange(2000)
            series=-p*np.sum(r**m*np.exp(1j*(2*m+1)*np.pi*z))
            self.assertLess(abs(series-local_transmission(x,z)),1e-13)

    def test_gaussian_against_independent_fourier(self):
        for b in [.005,.1,1.]:
            for det in [0.,.2,.5]:
                v,error=gaussian_average(.02,b,det)
                self.assertLess(abs(v-gaussian_fourier(.02,b,det)),1e-10)

    def test_phase_average_direct_integral(self):
        for x in [.02,.1]:
            result=quad(lambda z:abs(local_transmission(x,z))**2,-.5,.5,points=[0],epsabs=1e-12)[0]
            self.assertAlmostEqual(result,float(eventual_fraction(x)),places=11)

    def test_pulse_energy_and_complete_sum(self):
        x=.05;p=mouth_probability(x);b=2.
        direct=quad(lambda u:float(abs(pulse_amplitude(np.array(u),x,b,1))**2),0,2*np.pi,epsabs=1e-13)[0]
        self.assertAlmostEqual(direct,p*p,places=12)
        N=1000; m=np.arange(N)
        self.assertAlmostEqual(float(collected_fraction(x,N)),float(np.sum(p*p*(1-p)**(2*m))),places=12)

    def test_collection_quantile(self):
        for x in [1e-6,.02,.1]:
            n=int(np.ceil(collection_pulses(x)))
            target=.9*eventual_fraction(x)
            self.assertGreaterEqual(collected_fraction(x,n),target*(1-1e-14))
            self.assertLess(collected_fraction(x,n-1),target)

    def test_domains(self):
        for x in [0,-1,.2]:
            with self.assertRaises(ValueError):mouth_probability(x)
        with self.assertRaises(ValueError):gaussian_average(.01,0)
        with self.assertRaises(ValueError):collection_pulses(.01,1)

    def test_compact_pulse_orthogonality(self):
        for duration in [.4,1.4]:
            for chirp in [0.,7.]:
                self.assertAlmostEqual(compact_autocorrelation(0,duration,chirp),1.,places=12)
                self.assertEqual(compact_autocorrelation(2*np.pi,duration,chirp),0.)

    def test_long_pulse_has_overlap(self):
        self.assertGreater(compact_autocorrelation(2*np.pi,8,0).real,0)
        self.assertGreater(abs(compact_autocorrelation(2*np.pi,8,7)),0)


if __name__=='__main__':unittest.main()
