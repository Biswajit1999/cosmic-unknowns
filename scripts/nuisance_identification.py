"""Bounded deterministic nuisance fits. Author: Biswajit Jana.

Fit residuals are achievable witnesses, not certified global minima or p-values.
"""
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares

BOUNDS = dict(scale=(.95, 1.05), mass=(.10, .44), false_positive=(0., .03), false_negative=(0., .03))
STAGES = ['fixed', 'mass', 'mass_scale', 'mass_scale_response']


def centered_parameters(lam, centered_mass):
    denominator = 1+np.pi*lam/2
    if abs(lam) >= 2/np.pi or centered_mass < 0:
        raise ValueError('Require positive endpoint lapse and nonnegative centered mass')
    return lam/denominator, centered_mass*denominator


def reversed_parameters(beta, mass):
    factor = 1-np.pi*beta
    if factor <= 0 or mass < 0:
        raise ValueError('Require regular geometry and nonnegative mass')
    return -beta/factor, mass*factor


class NullSpectrum:
    def __init__(self, rows):
        self.functions = {}
        for ell in [0, 1]:
            rr = sorted([r for r in rows if int(r['ell']) == ell], key=lambda r: float(r['q']))
            self.functions[ell] = PchipInterpolator([float(r['q']) for r in rr], [float(r['T']) for r in rr], extrapolate=False)

    def predict(self, k, ell, parameters):
        k = np.asarray(k); ell = np.asarray(ell)
        scale, mass, b, c = parameters
        if mass < 0 or b < 0 or c < 0 or b+c > 1:
            raise ValueError('Unphysical mass or response probabilities')
        q2 = (scale*k)**2-mass**2
        if np.any(q2 <= 0):
            raise ValueError('Null channels must remain open')
        q = np.sqrt(q2); values = np.empty_like(k, dtype=float)
        for mode in np.unique(ell):
            selected = ell == mode
            values[selected] = self.functions[int(mode)](q[selected])
        if not np.all(np.isfinite(values)):
            raise ValueError('No template extrapolation permitted')
        return b+(1-b-c)*values

    def fit(self, k, ell, target, stage):
        if stage not in STAGES:
            raise ValueError('Unknown fit stage')
        base = np.array([1., .3, 0., 0.])
        active = {'fixed': [], 'mass': [1], 'mass_scale': [0, 1], 'mass_scale_response': [0, 1, 2, 3]}[stage]
        trials = []
        if not active:
            residual = self.predict(k, ell, base)-target
            return base, residual, [dict(start=0, cost=float(np.dot(residual, residual)), success=True, evaluations=1)]
        limits = np.array([BOUNDS[key] for key in ['scale', 'mass', 'false_positive', 'false_negative']])
        starts = [(s, m) for s in ([1.] if stage == 'mass' else [.96, 1., 1.04]) for m in [.18, .3, .42]]
        best = None
        for i, (s, m) in enumerate(starts):
            initial = np.array([s, m, .004, .004])
            def residual(x):
                parameters = base.copy(); parameters[active] = x
                return self.predict(k, ell, parameters)-target
            fit = least_squares(residual, initial[active], bounds=(limits[active, 0], limits[active, 1]),
                ftol=2e-12, xtol=2e-12, gtol=2e-12, max_nfev=3000)
            cost = float(np.dot(fit.fun, fit.fun))
            trials.append(dict(start=i, cost=cost, success=bool(fit.success), evaluations=fit.nfev))
            if best is None or cost < best[0]:
                parameters = base.copy(); parameters[active] = fit.x
                best = (cost, parameters, fit.fun.copy())
        return best[1], best[2], trials
