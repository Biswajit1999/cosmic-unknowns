"""Regular conformal metrics with a common scalar operator. Author: Biswajit Jana.

Units a=c=1. This constructs metrics and required Einstein stress, not a
supporting matter action or a claim of coupled dynamical stability.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


def geometry(x, beta):
    if not math.isfinite(beta) or beta >= 1/math.pi:
        raise ValueError('Require finite beta < 1/pi for two regular ends')
    x = np.asarray(x)
    root = np.sqrt(1+x*x)
    f = 1+beta*(np.arctan(x)-math.pi/2)
    fp = beta/(1+x*x); fpp = -2*beta*x/(1+x*x)**2
    r = f*root
    rp = fp*root+f*x/root
    rpp = fpp*root+2*fp*x/root+f/root**3
    return f, fp, fpp, r, rp, rpp


def metric_potential(x, beta, ell=0, xi=0.):
    f, fp, fpp, r, rp, rpp = geometry(x, beta)
    ricci = -2/(f*f*(1+np.asarray(x)**2)**2)
    return ell*(ell+1)*f*f/(r*r)+rpp/r+xi*f*f*ricci


def common_potential(x, ell=0, xi=0.):
    return ell*(ell+1)/(1+x*x)+(1-2*xi)/(1+x*x)**2


def stress(x, beta):
    """8 pi G a^2/c^4 times energy density and principal pressures."""
    f, fp, fpp, r, rp, rpp = geometry(x, beta)
    rl = rp/f; rll = rpp/f**2-rp*fp/f**3
    phil = fp/f**2; phill = fpp/f**3-2*fp**2/f**4
    rho = (1-rl*rl-2*r*rll)/r**2
    radial = (-1+rl*rl+2*r*rl*phil)/r**2
    transverse = phill+phil*phil+(phil*rl+rll)/r
    return rho, radial, transverse


def throat(beta):
    root = brentq(lambda x: float(geometry(x, beta)[4]), -100, 100, xtol=1e-14)
    f, fp, fpp, r, rp, rpp = geometry(root, beta)
    rho, pr, pt = stress(root, beta)
    return dict(beta=beta, X=float(root), radius=float(r), lapse_sqrt=float(f),
                left_lapse_sqrt=1-math.pi*beta, ricci=float(-2/(f*f*(1+root*root)**2)),
                radial_nec=float(rho+pr), flare_out=float(rpp/f**2))


def tail_coefficients(k, ell, order=12, xi=0.):
    """u_out=exp(ikx) sum a_n/x^n, from the exact common potential.

    U_{2j+2}=(-1)^j [ell(ell+1)-(1-2xi)j].
    The asymptotic series is truncated, not asserted to converge.
    """
    if k <= 0 or ell < 0 or int(ell) != ell or order < 1:
        raise ValueError('Invalid wave or tail parameters')
    a = np.zeros(order+1, dtype=complex); a[0] = 1
    angular = ell*(ell+1)
    for m in range(order):
        convolution = 0j
        for power in range(2, m+3, 2):
            j = (power-2)//2
            convolution += (-1)**j*(angular-(1-2*xi)*j)*a[m+2-power]
        a[m+1] = (m*(m+1)*a[m]-convolution)/(2j*k*(m+1))
    return a


def jost_tail(x, k, ell, order=12, xi=0.):
    a = tail_coefficients(k, ell, order, xi)
    n = np.arange(len(a)); powers = x**(-n.astype(float))
    series = np.dot(a, powers)
    derivative = np.dot(-n*a, powers)/x
    phase = np.exp(1j*k*x)
    return phase*series, phase*(1j*k*series+derivative)


def scatter(k, ell=0, beta=0., form='radial', domain=None, order=12,
            rtol=2e-10, atol=2e-12, xi=0.):
    """Right incidence; common tortoise-coordinate wave normalization.

    radial integrates the original field phi=u/r in the full metric.
    reduced independently integrates u with its metric-derived potential.
    Boundary asymptotics include the angular 1/x^2 tail.
    """
    if form not in ['radial', 'reduced']:
        raise ValueError('Unknown equation form')
    domain = max(80., 30/k) if domain is None else float(domain)
    j, jp = jost_tail(domain, k, ell, order, xi)
    left_u = j; left_du = -jp
    if form == 'radial':
        _, _, _, r, rp, _ = geometry(-domain, beta)
        initial = [left_u/r, (left_du-rp*left_u/r)/r]
        def rhs(x, y):
            f, fp, fpp, r, rp, rpp = geometry(x, beta)
            curvature_term = -2*xi/(1+x*x)**2
            return [y[1], -2*rp/r*y[1]-(k*k-ell*(ell+1)*f*f/r**2-curvature_term)*y[0]]
    else:
        initial = [left_u, left_du]
        def rhs(x, y):
            return [y[1], (metric_potential(x, beta, ell, xi)-k*k)*y[0]]
    sol = solve_ivp(rhs, (-domain, domain), initial, method='DOP853',
                    rtol=rtol, atol=atol, max_step=min(2., .8/k))
    if not sol.success:
        raise RuntimeError(sol.message)
    u, du = sol.y[:, -1]
    if form == 'radial':
        _, _, _, r, rp, _ = geometry(domain, beta)
        u, du = r*u, rp*u+r*du
    incoming, outgoing = np.linalg.solve(np.array([[j.conjugate(), j], [jp.conjugate(), jp]]), [u, du])
    transmission, reflection = 1/incoming, outgoing/incoming
    return dict(k=k, ell=ell, beta=beta, form=form, domain=domain, order=order, xi=xi,
                t=complex(transmission), r=complex(reflection), T=float(abs(transmission)**2),
                R=float(abs(reflection)**2), flux_error=float(abs(abs(transmission)**2+abs(reflection)**2-1)),
                evaluations=sol.nfev)
