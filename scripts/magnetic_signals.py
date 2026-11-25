"""Study 014: low-frequency magnetic scalar signals. Author: Biswajit Jana.

Input: Freivogel/Fumagalli/Tomasevic arXiv:2606.12528v2 Eqs.20-22.
The frozen-mouth response is an approximation, not a full spacetime solver.
"""
import numpy as np
from scipy.integrate import quad
from scipy.special import erfcx


def mouth_probability(x):
    x = np.asarray(x, dtype=float)
    if np.any((x <= 0) | (x >= .2)):
        raise ValueError('Require 0 < omega*r_e < 0.2; low-frequency model only')
    return 4*x*x/(1+x*x)**2


def amplitudes(x, z):
    """Published matching amplitudes. z=omega*ell. No frozen-x assumption here."""
    x = np.asarray(x); z = np.asarray(z)
    mouth_probability(x)
    a = x*x; phase = np.pi*z
    denominator = 2j*a*np.cos(phase)+(1+a*a)*np.sin(phase)
    return (a*a-1)*np.sin(phase)/denominator, -2j*a/denominator


def local_transmission(x, detuning):
    """Frozen x; detuning z-n evaluated directly to resolve very narrow peaks.

    Integer carrier parity changes only an overall sign, omitted here.
    expm1 avoids cancellation; exact-resonance limit is -1.
    """
    p = mouth_probability(x); r = 1-p
    phase = np.pi*np.asarray(detuning)
    return -p*np.exp(1j*phase)/(p-r*np.expm1(2j*phase))


def matching_solve(x, z):
    """Independent four-equation throat/far-region matching, unit incident wave."""
    a = x*x; ep = np.exp(.5j*np.pi*z); em = ep.conjugate()
    matrix = np.array([[em, ep, -1j, 0], [1j*a*em, -1j*a*ep, -1, 0],
                       [ep, em, 0, -1j], [-1j*a*ep, 1j*a*em, 0, -1]], dtype=complex)
    cin, cout, reflected, transmitted = np.linalg.solve(matrix, [-1j, 1, 0, 0])
    return reflected, transmitted


def halfwidth(x):
    p = mouth_probability(x)
    return np.arcsin(p/(2*np.sqrt(1-p)))/np.pi


def eventual_fraction(x):
    p = mouth_probability(x)
    return p/(2-p)


def collected_fraction(x, number_of_pulses):
    """Separated pulses: fraction of incident energy after N complete pulses."""
    n = np.asarray(number_of_pulses)
    if np.any(n < 0): raise ValueError('Pulse count must be nonnegative')
    p = mouth_probability(x)
    return eventual_fraction(x)*(-np.expm1(2*n*np.log1p(-p)))


def collection_pulses(x, fraction=.9):
    """Continuous count threshold; ceil for number of complete pulses."""
    if not 0 < fraction < 1: raise ValueError('Require 0 < fraction < 1')
    p = mouth_probability(x)
    return np.log1p(-fraction)/(2*np.log1p(-p))


def gaussian_average(x, bandwidth, detuning=0., tolerance=1e-10):
    """Exact frozen-x Gaussian energy average via resonance-flattening variable.

    bandwidth is the standard deviation of incident energy in z=omega*ell.
    Periodize the Gaussian; map tan(pi delta)=k tan(theta), k=p/(2-p).
    Under this map |T|^2 d(delta)=k d(theta)/pi. Avoid missing narrow peaks.
    """
    if bandwidth <= 0: raise ValueError('Bandwidth must be positive')
    p = float(mouth_probability(x)); k = p/(2-p)
    center = (detuning+.5) % 1-.5
    count = max(2, int(np.ceil(9*bandwidth+abs(center)+1)))
    integers = np.arange(-count, count+1)
    norm = np.sqrt(2*np.pi)*bandwidth
    def integrand(theta):
        delta = np.arctan(k*np.tan(theta))/np.pi
        w = np.sum(np.exp(-.5*((delta+integers-center)/bandwidth)**2))/norm
        return k*w/np.pi
    points = [-np.pi/2, np.pi/2, 0.]
    # A packet centered on the phase-cell edge has two wrapped peaks.
    # Split around both; otherwise adaptive quadrature can miss one of them.
    for image in [-1, 0, 1]:
        for multiple in [-8, -3, -1, 0, 1, 3, 8]:
            delta = center+image+multiple*bandwidth
            if -.5 < delta < .5:
                points.append(float(np.arctan(np.tan(np.pi*delta)/k)))
    points = sorted(set(points)); answer = 0.; estimate = 0.
    for left, right in zip(points[:-1], points[1:]):
        v, e = quad(integrand, left, right, epsabs=tolerance*k,
                    epsrel=tolerance, limit=200)
        answer += v; estimate += e
    return answer, estimate


def gaussian_fourier(x, bandwidth, detuning=0.):
    """Independent Poisson-kernel Fourier sum; efficient for bandwidth >= .005."""
    if bandwidth < .005: raise ValueError('Use peak-resolving integration for smaller bandwidth')
    p = float(mouth_probability(x))
    count = int(np.ceil(np.sqrt(45/(2*np.pi**2))/bandwidth))+1
    k = np.arange(1, count+1)
    correction = np.exp(k*np.log1p(-p)-2*np.pi**2*bandwidth**2*k*k)*np.cos(2*np.pi*k*detuning)
    return float(eventual_fraction(x)*(1+2*np.sum(correction)))


def gaussian_lorentzian(x, bandwidth):
    """Single-resonance approximation, centered packet; not valid over many peaks."""
    p = mouth_probability(x); h = p/(2*np.pi*np.sqrt(1-p))
    ratio = h/bandwidth
    return np.sqrt(np.pi/2)*ratio*erfcx(ratio/np.sqrt(2))


def pulse_amplitude(u, x, bandwidth, pulses=12, detuning=0.):
    """Normalized energy envelope of the frozen-response delayed pulse train.

    u=(t-r0-r_observer)/ell, each pulse centered at (2m+1)pi.
    Carrier removed; independent round-trip phases retained if detuning != 0.
    Gaussian has noncompact tails. A hard causal turn-on is not claimed.
    """
    u = np.asarray(u); p = float(mouth_probability(x))
    norm = (2*bandwidth**2/np.pi)**.25
    result = np.zeros_like(u, dtype=complex)
    for m in range(pulses):
        delay = (2*m+1)*np.pi
        result -= p*np.exp(m*np.log1p(-p)+2j*np.pi*m*detuning)*norm*np.exp(-bandwidth**2*(u-delay)**2)
    return result


def compact_envelope(u, duration, chirp=0.):
    """Unit-energy compact cos-squared envelope; effective response test only.

    Strict time support implies spectral tails, so this is not an exactly
    band-limited physical solution of the low-frequency matching problem.
    """
    if duration <= 0: raise ValueError('Duration must be positive')
    u=np.asarray(u); inside=abs(u)<duration/2
    return np.where(inside, np.sqrt(8/(3*duration))*np.cos(np.pi*u/duration)**2
                    *np.exp(1j*chirp*(u/duration)**2), 0j)


def compact_autocorrelation(delay, duration, chirp=0.):
    lo=max(-duration/2,delay-duration/2); hi=min(duration/2,delay+duration/2)
    if lo >= hi:return 0j
    def product(u):return np.conj(compact_envelope(u,duration,chirp))*compact_envelope(u-delay,duration,chirp)
    return (quad(lambda u:float(product(u).real),lo,hi,epsabs=1e-12)[0]
            +1j*quad(lambda u:float(product(u).imag),lo,hi,epsabs=1e-12)[0])
