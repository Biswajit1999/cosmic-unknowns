"""Bounded channel-dependent classification errors. Author: Biswajit Jana.

Identities concern ideal Bernoulli probabilities, not wormhole existence.
"""
import numpy as np
from scipy.optimize import minimize, differential_evolution
from scipy.special import xlogy
from transmission_inference import pt_transmission, pt_logit, BOUNDS


def required_error(p, q):
    """Smallest per-channel bound permitting q to be recorded as p."""
    p, q = np.asarray(p), np.asarray(q)
    if np.any((p < 0) | (p > 1)) or np.any((q <= 0) | (q >= 1)):
        raise ValueError("Require 0<=p<=1 and 0<q<1")
    return np.where(p >= q, (p-q)/(1-q), (q-p)/q)


def matching_errors(p, q):
    """Construct false-positive and false-negative probabilities explicitly."""
    p, q = np.asarray(p), np.asarray(q)
    required_error(p, q)
    false_positive = np.maximum(p-q, 0)/(1-q)
    false_negative = np.maximum(q-p, 0)/q
    return false_positive, false_negative


def recorded_probability(q, fp, fn):
    return (1-np.asarray(fn))*np.asarray(q)+np.asarray(fp)*(1-np.asarray(q))


def clipped_probability(p, q, delta):
    if not 0 <= delta < 1:
        raise ValueError("Require 0<=delta<1")
    return np.clip(p, (1-delta)*q, q+delta*(1-q))


def bernoulli_kl(p, q):
    p = np.asarray(p)
    q = np.clip(q, 1e-15, 1-1e-15)
    return xlogy(p, p)-xlogy(p, q)+xlogy(1-p, 1-p)-xlogy(1-p, 1-q)


def fit_error_budget(k, p, starts):
    """Best located feasible witness. No certified global optimality claim."""
    k, p = np.asarray(k), np.asarray(p)
    bounds = [*BOUNDS, (0, .5)]
    def constraints(v):
        q = pt_transmission(k, v[:2])
        delta = v[2]
        return np.r_[p-(1-delta)*q, q+delta*(1-q)-p]
    def jac(v):
        q = pt_transmission(k, v[:2])
        z = np.exp(v[1])*k
        dq = q[:, None]*(1-q)[:, None]*np.column_stack([-np.ones_like(k), 2*z/np.tanh(z)])
        return np.vstack([np.column_stack([-(1-v[2])*dq, q]),
                          np.column_stack([(1-v[2])*dq, 1-q])])
    attempts = []
    for start in starts:
        delta = float(np.max(required_error(p, pt_transmission(k, start))))
        result = minimize(lambda v: v[2], np.r_[start, min(delta+.001, .49)],
            jac=lambda v: np.array([0., 0., 1.]), method="SLSQP", bounds=bounds,
            constraints=[{"type": "ineq", "fun": constraints, "jac": jac}],
            options={"ftol": 2e-13, "maxiter": 1200})
        actual = float(np.max(required_error(p, pt_transmission(k, result.x[:2]))))
        attempts.append({"theta": result.x[:2].tolist(), "delta": actual,
                         "success": bool(result.success), "constraint_min": float(min(constraints(result.x)))})
    best = min(attempts, key=lambda r: r["delta"])
    return np.array(best["theta"]), best["delta"], attempts


def global_search_check(k, p, seed):
    """Independent stochastic search supports, but cannot certify, a minimum."""
    def objective(theta):
        q = np.clip(pt_transmission(k, theta), 1e-14, 1-1e-14)
        return float(np.max(required_error(p, q)))
    result = differential_evolution(objective, BOUNDS, seed=seed, popsize=20,
        maxiter=700, tol=1e-9, atol=1e-12, polish=False)
    return {"theta": result.x.tolist(), "delta": float(result.fun), "success": bool(result.success)}


def profile_kl(k, p, delta, starts):
    """Locally profiled KL with freely varying, bounded channel errors."""
    def objective(theta):
        q = pt_transmission(k, theta)
        return float(np.sum(bernoulli_kl(p, clipped_probability(p, q, delta))))
    results = [minimize(objective, start, method="Nelder-Mead", bounds=BOUNDS,
               options={"xatol": 1e-9, "fatol": 1e-14, "maxiter": 1500}) for start in starts]
    best = min(results, key=lambda r: r.fun)
    return best.x, max(0., float(best.fun)), bool(best.success)
