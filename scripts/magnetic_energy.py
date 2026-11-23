"""MMP finite-separation reduced energy. Author: Biswajit Jana.

Reproduces arXiv:1807.04726v3 Eq. (5.47), with independently derived
prolate coordinates for the flat-exterior magnetic field-line lengths.
This is not a full semiclassical spacetime or perturbation solver.
"""
from functools import lru_cache
import numpy as np
from scipy.special import roots_legendre
from scipy.optimize import brentq


@lru_cache(maxsize=16)
def quadrature(order):
    return roots_legendre(order)


def field_coordinates(eta, nu):
    """Mouth centers z/d=+-1/2; eta in [-1,1], 0<nu<=2."""
    if not 0 < nu <= 2:
        raise ValueError('Field label must obey 0 < nu <= 2')
    eta = np.asarray(eta)
    if np.any(abs(eta) > 1):
        raise ValueError('Prolate eta must lie in [-1,1]')
    a = 1-eta*eta
    disc = np.sqrt(a*a+nu*nu*eta*eta)
    xm = a*(2-nu)/(nu+nu*nu*eta*eta/(disc+a))
    xi = 1+xm
    return xi*eta/2, np.sqrt(xm*(xi+1)*a)/2


def field_length(nu, order=512):
    """Full exterior length/d. Integrate the prolate flat-space line element."""
    values = np.asarray(nu, dtype=float)
    if np.any((values <= 0) | (values > 2)):
        raise ValueError('Field labels must obey 0 < nu <= 2')
    flat = values.ravel(); result = np.ones_like(flat)
    eta, weights = quadrature(order)
    mask = flat < 2
    n = flat[mask, None]; e = eta[None, :]; a = 1-e*e
    disc = np.sqrt(a*a+n*n*e*e)
    xm = a*(2-n)/(n+n*n*e*e/(disc+a)); xi = 1+xm
    derivative = e*(n-2*xi)/disc
    speed = .5*np.sqrt((xi*xi-e*e)*(derivative**2/(xm*(xi+1))+1/a))
    result[mask] = speed @ weights
    return result.reshape(values.shape)


class MagneticEnergy:
    """Dimensionless e=6dE/(pi q); s=d/d0; y=pi ell/d.

    cutoff removes 0<nu<cutoff without renormalizing the remaining measure.
    It is a sensitivity experiment, not a derived finite-charge correction.
    """
    def __init__(self, outer=512, inner=512, cutoff=0., shortest=False):
        if not 0 <= cutoff < 1.5:
            raise ValueError('Require cutoff < 1.5 for a negative-energy minimum')
        nodes, weights = quadrature(outer)
        self.nu = cutoff+(nodes+1)*(2-cutoff)/2
        self.weights = weights*(2-cutoff)/4
        self.f = np.ones_like(self.nu) if shortest else field_length(self.nu, inner)
        self.cutoff = cutoff

    def moments(self, y, power):
        if y <= 0: raise ValueError('y must be positive')
        return float(np.sum(self.weights/(y+self.f)**power))

    def energy(self, y, s):
        if s <= 0: raise ValueError('s must be positive')
        return 3/(8*s*y*y)+1/(4*y)-self.moments(y, 1)

    def stationarity(self, y, s=np.inf):
        if s <= 0 or y <= 0: raise ValueError('Require s,y>0')
        return float(np.sum(self.weights*(y/(y+self.f))**2))-.25-.75/(s*y)

    def hprime(self, y, s):
        return float(np.sum(self.weights*2*y*self.f/(y+self.f)**3))+.75/(s*y*y)

    def minimum(self, s=np.inf):
        if s <= 0: raise ValueError('s must be positive')
        upper = max(4., 4/s)
        while self.stationarity(upper, s) <= 0: upper *= 2
        y = brentq(lambda v: self.stationarity(v, s), .5, upper, xtol=1e-12, rtol=1e-13)
        return dict(s=float(s), y=float(y), throat_over_d0=float(s*y),
                    energy=self.energy(y, s), curvature=self.hprime(y, s)/y**2,
                    stationarity=self.stationarity(y, s))


def scale_diagnostics(q, g, s, y, energy):
    """Natural units hbar=c=1. These are necessary scale checks, not sufficiency."""
    if min(q, g, s, y) <= 0: raise ValueError('Scales must be positive')
    ell0_over_re = 16*np.pi*q/g**2
    ell_over_re = ell0_over_re*s*y
    d_over_re = np.pi*ell0_over_re*s
    return dict(ell_over_re=ell_over_re, d_over_re=d_over_re,
                schwarzian_coefficient=q/(16*s*y),
                orbital_omega_ell=2*ell_over_re/d_over_re**1.5,
                binding_times_ell=-q*y*energy/6)
