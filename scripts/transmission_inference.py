"""Restricted scalar-barrier inference. Author: Biswajit Jana.

These utilities model ideal synthetic counting experiments, not telescope data.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from scipy.special import expit, xlogy
from scipy.stats import norm


def wilson(successes, trials, confidence=0.95):
    """Pointwise Wilson score interval; not a simultaneous confidence band."""
    s = np.asarray(successes, dtype=float)
    if trials <= 0 or np.any(s < 0) or np.any(s > trials) or not 0 < confidence < 1:
        raise ValueError("Invalid binomial sample or confidence level")
    z = norm.ppf((1+confidence)/2)
    p = s/trials
    d = 1+z*z/trials
    center = (p+z*z/(2*trials))/d
    half = z*np.sqrt(p*(1-p)/trials+z*z/(4*trials**2))/d
    return np.maximum(0, center-half), np.minimum(1, center+half)


def pt_logit(k, theta):
    """theta=(log C, log b); T=sinh(b K)^2/[sinh(b K)^2+C]."""
    k = np.asarray(k)
    if np.any(k <= 0):
        raise ValueError("Require positive wave numbers")
    z = np.exp(theta[1])*k
    # logsinh without overflow for large positive arguments.
    return 2*(z+np.log(-np.expm1(-2*z))-np.log(2))-theta[0]


def pt_transmission(k, theta):
    return expit(pt_logit(k, theta))


def physical_parameters(theta):
    """Return height H and inverse width alpha for H sech(alpha X)^2."""
    c, b = np.exp(theta)
    alpha = np.pi/b
    if c >= 1:
        factor = 1+(2*np.arccosh(np.sqrt(c))/np.pi)**2
    else:
        factor = 1-(2*np.arccos(np.sqrt(c))/np.pi)**2
    return float(alpha*alpha*factor/4), float(alpha)


def pt_parameters(height, alpha):
    if height <= 0 or alpha <= 0:
        raise ValueError("Require positive barrier height and inverse width")
    d = 4*height/alpha**2-1
    c = np.cosh(np.pi*np.sqrt(d)/2)**2 if d >= 0 else np.cos(np.pi*np.sqrt(-d)/2)**2
    return np.array([np.log(c), np.log(np.pi/alpha)])


def deviance(counts, trials, probability):
    """Binomial deviance against saturated frequencies, including zero counts."""
    y = np.asarray(counts)/trials
    p = np.clip(probability, 1e-15, 1-1e-15)
    terms = xlogy(y, y)-xlogy(y, p)+xlogy(1-y, 1-y)-xlogy(1-y, 1-p)
    return float(max(0., 2*trials*np.sum(terms)))


BOUNDS = [(-16., 25.), (np.log(0.1), np.log(30.))]


def fit_pt(k, fractions, start=None, multistart=False):
    """Binomial MLE for equal trial counts; return fit and optimization audit.

    Scaling by trial count does not change the optimizer. The constant saturated
    likelihood is omitted. Bounds are part of the restricted comparator family.
    """
    k, y = np.asarray(k), np.asarray(fractions)
    def objective(theta):
        z = np.exp(theta[1])*k
        eta = pt_logit(k, theta)
        p = expit(eta)
        jac = np.array([-np.ones_like(k), 2*z/np.tanh(z)])
        return float(np.sum(np.logaddexp(0, eta)-y*eta)), jac@(p-y)
    starts = [np.array(start) if start is not None else pt_parameters(1., 1.)]
    if multistart:
        starts += [pt_parameters(0.4, 0.5), pt_parameters(2., 2.)]
    fits = [minimize(objective, s, method="L-BFGS-B", jac=True, bounds=BOUNDS,
                     options={"ftol": 1e-14, "gtol": 1e-9, "maxiter": 400}) for s in starts]
    best = min(fits, key=lambda f: f.fun)
    active = any(abs(best.x[i]-edge) < 1e-5 for i, bounds in enumerate(BOUNDS) for edge in bounds)
    return best.x, bool(best.success), bool(active), float(np.max(np.abs(best.jac)))


def numerical_pt(k, height, alpha):
    """Independent ODE check of the analytic comparator."""
    length = 24/alpha
    def rhs(x, y):
        return [y[1], (height/np.cosh(alpha*x)**2-k*k)*y[0]]
    phase = np.exp(1j*k*length)
    sol = solve_ivp(rhs, (length, -length), [phase, 1j*k*phase], method="DOP853",
                    rtol=2e-11, atol=2e-13, max_step=min(0.3/alpha, 0.6/k))
    if not sol.success:
        raise RuntimeError(sol.message)
    u, du = sol.y[:, -1]
    incoming = 0.5*np.exp(1j*k*length)*(u+du/(1j*k))
    return float(abs(1/incoming)**2)
