"""Independent validation with simultaneous confidence regions. Author: Biswajit Jana."""
import numpy as np
from calibrated_identification import cp_interval, physical_interval, continuum_test


def drift_bounds(rl, ru, vl, vu):
    """Simultaneous bounds on the maximum absolute change of supplied rates."""
    lower = np.maximum(0, np.maximum(vl-ru, rl-vu))
    upper = np.maximum(vu-rl, ru-vl)
    return float(np.max(lower)), float(np.max(upper))


def validation_regions(y, n, rf, rg, m, vf, vg, v, eta=.003, alpha=.05):
    """Five binomial intervals/channel share one fixed confidence budget.

    V probes each known class at the validation response. Only a prespecified
    residual bound eta connects validation to science. No shared-response
    assumption is inferred merely from contemporaneous timestamps.
    """
    if not 0 <= eta <= 1 or v < 0:
        raise ValueError('Invalid residual allowance or validation budget')
    channels = len(y)
    tail = alpha/(10*channels)
    pl, pu = cp_interval(y, n, tail)
    fl, fu = cp_interval(rf, m, tail)
    gl, gu = cp_interval(rg, m, tail)
    if v:
        vfl, vfu = cp_interval(vf, v, tail)
        vgl, vgu = cp_interval(vg, v, tail)
    else:
        vfl = vgl = np.zeros(channels)
        vfu = vgu = np.ones(channels)
    dl, du = drift_bounds(np.r_[fl, gl], np.r_[fu, gu],
                          np.r_[vfl, vgl], np.r_[vfu, vgu])
    rectangles = {
        'stale_reference': (fl, fu, gl, gu),
        'scalar_transfer': (np.maximum(0, fl-du), np.minimum(1, fu+du),
                            np.maximum(0, gl-du), np.minimum(1, gu+du)),
        'direct_validation': (vfl, vfu, vgl, vgu),
        'guarded_validation': (np.maximum(0, vfl-eta), np.minimum(1, vfu+eta),
                               np.maximum(0, vgl-eta), np.minimum(1, vgu+eta)),
    }
    regions = {name: physical_interval(pl, pu, *box) for name, box in rectangles.items()}
    return regions, (dl, du)


def evaluate(y, n, rf, rg, m, vf, vg, v, k, eta=.003):
    regions, bounds = validation_regions(y, n, rf, rg, m, vf, vg, v, eta)
    return {name: continuum_test(k, *region) for name, region in regions.items()}, bounds
