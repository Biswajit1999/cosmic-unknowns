"""Joint calibration confidence regions and continuum exclusion.

Author: Biswajit Jana. Statistical statements are model-conditional.
Interval certificates apply to the supplied numerical probability bounds.
"""
import numpy as np
from scipy.special import logit
from scipy.stats import beta
from mpmath import iv

iv.dps = 35


def cp_interval(y, n, tail):
    """Two-sided Clopper-Pearson interval, probability tail per side.

    Numerical beta quantiles are expanded by 1e-12; this is not a formal
    interval implementation of the inverse incomplete beta function.
    """
    y = np.asarray(y, dtype=int)
    if n <= 0 or np.any((y < 0) | (y > n)) or not 0 < tail < .5:
        raise ValueError('Invalid binomial experiment')
    lo = np.zeros(y.shape); hi = np.ones(y.shape)
    use = y > 0
    lo[use] = beta.ppf(tail, y[use], n-y[use]+1)
    use = y < n
    hi[use] = beta.ppf(1-tail, y[use]+1, n-y[use])
    return np.maximum(0, lo-1e-12), np.minimum(1, hi+1e-12)


def physical_interval(pl, pu, fl, fu, gl, gu):
    """Invert response-rectangle overlap, including slopes of either sign.

    recorded_min = gl+(1-fu-gl)q <= pu
    recorded_max = gu+(1-fl-gu)q >= pl.
    """
    args = [np.asarray(v, dtype=float) for v in [pl, pu, fl, fu, gl, gu]]
    pl, pu, fl, fu, gl, gu = np.broadcast_arrays(*args)
    low = np.zeros(pl.shape); high = np.ones(pl.shape)
    for slope, rhs in [(1-fu-gl, pu-gl), (-(1-fl-gu), gu-pl)]:
        positive, negative, zero = slope > 0, slope < 0, slope == 0
        high[positive] = np.minimum(high[positive], rhs[positive]/slope[positive])
        low[negative] = np.maximum(low[negative], rhs[negative]/slope[negative])
        high[zero & (rhs < 0)] = -1
    return np.maximum(0, low-1e-11), np.minimum(1, high+1e-11)


def joint_intervals(y, n, false_negatives, false_positives, m, alpha=.05):
    """72 intervals for 24 channels; union-bound failure <= alpha.

    M known-class injections PER class PER channel. Independent calibration
    probes must have the same response as the science experiment.
    """
    channels = len(y)
    tail = alpha/(6*channels)
    pl, pu = cp_interval(y, n, tail)
    if m == 0:
        return np.zeros(channels),np.ones(channels)
    fl, fu = cp_interval(false_negatives, m, tail)
    gl, gu = cp_interval(false_positives, m, tail)
    return physical_interval(pl, pu, fl, fu, gl, gu)


def _logsinh2(b, k):
    z = b*k
    return 2*(z + np.log(-np.expm1(-2*z))-np.log(2))


def _iv_logsinh2(b, k):
    z = iv.mpf(float(b))*iv.mpf(float(k))
    return 2*iv.log((iv.exp(z)-iv.exp(-z))/2)


def _iv_logit(q):
    v = iv.mpf(float(q))
    return iv.log(v/(1-v))


def continuum_test(k, qlo, qhi, max_nodes=12000, min_width=1e-11,
                   return_certificate=False):
    """Exclude all logC in [-16,25] and b in [0.1,30], or do not reject.

    S_i(b)=2 log sinh(b K_i) is increasing. A branch [l,r] requires
    max_i(S_i(l)-logit(qhi_i), -16) <=
    min_j(S_j(r)-logit(qlo_j), 25).
    Every excluded branch is rechecked with 35-digit interval arithmetic.
    Inconclusive computation NEVER counts as rejection.
    """
    k, qlo, qhi = [np.asarray(v, dtype=float) for v in [k, qlo, qhi]]
    if np.any(k <= 0) or np.any(qlo < 0) or np.any(qhi > 1):
        raise ValueError('Invalid probabilities or frequencies')
    if np.any(qlo > qhi):
        return {'status':'excluded', 'nodes':0, 'leaves':[], 'empty_channel':int(np.argmax(qlo-qhi))}
    lower_offset = -logit(qhi)
    upper_offset = -logit(qlo)
    stack = [(.1, 30.)]; leaves=[]; nodes=0
    while stack and nodes < max_nodes:
        bl, br = stack.pop(); nodes += 1
        lows = np.r_[-16., _logsinh2(bl, k)+lower_offset]
        highs = np.r_[25., _logsinh2(br, k)+upper_offset]
        i, j = int(np.argmax(lows)), int(np.argmin(highs))
        if lows[i] > highs[j]+1e-10:
            lv = iv.mpf(-16) if i == 0 else _iv_logsinh2(bl,k[i-1])-_iv_logit(qhi[i-1])
            uv = iv.mpf(25) if j == 0 else _iv_logsinh2(br,k[j-1])-_iv_logit(qlo[j-1])
            if bool((lv-uv).a > 0):
                leaves.append([bl,br,i,j])
                continue
        mid = (bl+br)/2
        s = _logsinh2(mid,k)
        clo=max(-16.,float(np.max(s+lower_offset)))
        chi=min(25.,float(np.min(s+upper_offset)))
        if clo <= chi:
            return {'status':'compatible','nodes':nodes,'theta':[.5*(clo+chi),float(np.log(mid))]}
        if br-bl < min_width:
            return {'status':'undetermined','nodes':nodes}
        stack.extend([(mid,br),(bl,mid)])
    status = 'excluded' if not stack else 'undetermined'
    result={'status':status,'nodes':nodes,'leaf_count':len(leaves)}
    if return_certificate: result['leaves']=leaves
    return result


def zero_error_interval(p, delta, allowance=0.):
    """Population target with symmetric classification allowance and padding."""
    p=np.asarray(p)
    return physical_interval(np.maximum(0,p-allowance),np.minimum(1,p+allowance),
                             0,delta,0,delta)
