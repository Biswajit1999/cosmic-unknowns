"""Massive scalar scattering with logarithmic tail phases. Author: Biswajit Jana.

Fixed prescribed geometries; no claim of physical supporting matter or stability.
Only two open asymptotic channels, away from their exact thresholds, are supported.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
from conformal_scattering import geometry, common_potential


def potential(x, beta, mass, ell=0):
    return common_potential(x, ell)+mass**2*geometry(x, beta)[0]**2


def tail_data(k, beta, mass, ell, side, order=12):
    """Potential U=V0+V1/y+..., y=abs(X), and outgoing series coefficients."""
    if side not in [-1, 1] or order < 1 or int(order) != order:
        raise ValueError('Require side +/-1 and positive integer tail order')
    if not all(math.isfinite(v) for v in [k, beta, mass]) or k <= 0 or mass < 0:
        raise ValueError('Require finite positive K and nonnegative mass')
    if ell < 0 or int(ell) != ell:
        raise ValueError('Require nonnegative integer angular mode')
    geometry(0., beta)
    f = np.zeros(order+3)
    f[0] = 1 if side == 1 else 1-math.pi*beta
    for n in range(1, len(f), 2):
        f[n] = -side*beta*(-1)**((n-1)//2)/n
    v = mass**2*np.convolve(f, f)[:len(f)]
    for n in range(2, len(v), 2):
        j = (n-2)//2
        v[n] += (-1)**j*(ell*(ell+1)-j)
    if k*k <= v[0]:
        raise ValueError('Both asymptotic channels must be open; threshold/closed channel unsupported')
    q = math.sqrt(k*k-v[0]); eta = v[1]/(2*q); alpha = -1j*eta
    a = np.zeros(order+1, dtype=complex); a[0] = 1
    for n in range(order):
        convolution = sum(v[s]*a[n+2-s] for s in range(2, n+3))
        a[n+1] = ((n*(n+1)-2*alpha*n+alpha*alpha-alpha)*a[n]-convolution)/(2j*q*(n+1))
    return q, eta, v, a


def tail(y, k, beta, mass, ell=0, side=1, order=12):
    q, eta, v, a = tail_data(k, beta, mass, ell, side, order)
    n = np.arange(len(a)); powers = y**(-n.astype(float))
    series = np.dot(a, powers); derivative = np.dot(-n*a, powers)/y
    phase = np.exp(1j*(q*y-eta*np.log(2*q*y)))
    return phase*series, phase*((1j*q-1j*eta/y)*series+derivative), q, eta


def scatter_massive(k, beta=0., mass=.3, ell=0, incidence='right', form='reduced',
                    domain=None, order=12, rtol=2e-10, atol=2e-12):
    if incidence not in ['right', 'left'] or form not in ['reduced', 'radial']:
        raise ValueError('Invalid incidence or equation form')
    qr = tail_data(k, beta, mass, ell, 1, order)[0]
    ql = tail_data(k, beta, mass, ell, -1, order)[0]
    length = max(100., 40/min(qr, ql)) if domain is None else float(domain)
    if not math.isfinite(length) or length <= 1:
        raise ValueError('Boundary radius must exceed one')
    incoming_side = 1 if incidence == 'right' else -1
    outgoing_side = -incoming_side
    out, dout, qout, _ = tail(length, k, beta, mass, ell, outgoing_side, order)
    incoming_basis, din, qin, _ = tail(length, k, beta, mass, ell, incoming_side, order)
    start, end = outgoing_side*length, incoming_side*length
    u, du = out, outgoing_side*dout
    if form == 'radial':
        f, fp, fpp, radius, rp, rpp = geometry(start, beta)
        initial = [u/radius, (du-rp*u/radius)/radius]
        def rhs(x, state):
            f, fp, fpp, radius, rp, rpp = geometry(x, beta)
            coefficient = k*k-f*f*(ell*(ell+1)/radius**2+mass**2)
            return [state[1], -2*rp/radius*state[1]-coefficient*state[0]]
    else:
        initial = [u, du]
        def rhs(x, state):
            return [state[1], (potential(x, beta, mass, ell)-k*k)*state[0]]
    sol = solve_ivp(rhs, (start, end), initial, method='DOP853', rtol=rtol, atol=atol,
                    max_step=min(2., .8/max(qr, ql)))
    if not sol.success:
        raise RuntimeError(sol.message)
    u, du = sol.y[:, -1]
    if form == 'radial':
        f, fp, fpp, radius, rp, rpp = geometry(end, beta)
        u, du = radius*u, rp*u+radius*du
    derivative = incoming_side*din
    ai, ao = np.linalg.solve([[incoming_basis.conjugate(), incoming_basis],
                             [derivative.conjugate(), derivative]], [u, du])
    t = 1/ai; r = ao/ai
    tf = math.sqrt(qout/qin)*t
    transmission = abs(tf)**2; reflection = abs(r)**2
    jr = np.imag(np.conj(incoming_basis)*din)
    jl = np.imag(np.conj(out)*dout)
    return dict(k=k, beta=beta, mass=mass, ell=ell, incidence=incidence, form=form,
        domain=length, order=order, q_right=qr, q_left=ql, r=complex(r), t=complex(t),
        t_flux=complex(tf), T=float(transmission), R=float(reflection),
        flux_error=float(abs(transmission+reflection-1)),
        tail_current_error=float(max(abs(jr/qin-1), abs(jl/qout-1))), evaluations=sol.nfev)


def scalar_source_obstruction(x, beta):
    """Required scaled rho+p_t; zero for a single static radial scalar P(chi,Z)."""
    f = geometry(x, beta)[0]
    return -2*beta*np.asarray(x)/(f**3*(1+np.asarray(x)**2)**2)
