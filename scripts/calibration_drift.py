"""Calibration-transfer uncertainty. Author: Biswajit Jana.

Drift is an absolute change in each class-error probability, not a relative
amplitude/phase calibration error of a gravitational-wave instrument.
"""
import numpy as np
from scipy.optimize import minimize,differential_evolution
from calibrated_identification import cp_interval,physical_interval
from transmission_inference import pt_transmission,BOUNDS


def drift_intervals(y,n,fn,fp,m,epsilon,alpha=.05):
    """Project simultaneous calibration regions enlarged by a fixed drift bound."""
    if not 0<=epsilon<=1:raise ValueError('Drift allowance must lie in [0,1]')
    count=len(y);tail=alpha/(6*count)
    pl,pu=cp_interval(y,n,tail)
    if m==0:return np.zeros(count),np.ones(count)
    fl,fu=cp_interval(fn,m,tail);gl,gu=cp_interval(fp,m,tail)
    return physical_interval(pl,pu,np.maximum(0,fl-epsilon),np.minimum(1,fu+epsilon),
        np.maximum(0,gl-epsilon),np.minimum(1,gu+epsilon))


def minimax_comparator(k,p,starts):
    """Locate a small sup-norm mismatch; interval exclusion supplies lower bounds."""
    def constraints(v):
        diff=p-pt_transmission(k,v[:2]);return np.r_[v[2]-diff,v[2]+diff]
    attempts=[]
    for start in starts:
        result=minimize(lambda v:v[2],np.r_[start,.03],method='SLSQP',
            jac=lambda v:np.array([0.,0.,1.]),bounds=[*BOUNDS,(0,1)],
            constraints={'type':'ineq','fun':constraints},options={'ftol':1e-13,'maxiter':2000})
        gap=float(np.max(np.abs(p-pt_transmission(k,result.x[:2]))))
        attempts.append({'theta':result.x[:2].tolist(),'gap':gap,'success':bool(result.success)})
    de=differential_evolution(lambda theta:np.max(np.abs(p-pt_transmission(k,theta))),BOUNDS,
        seed=2026111207,popsize=20,maxiter=700,tol=1e-10,atol=1e-13,polish=False)
    attempts.append({'theta':de.x.tolist(),'gap':float(de.fun),'success':bool(de.success),'method':'differential_evolution'})
    best=min(attempts,key=lambda r:r['gap'])
    return np.array(best['theta']),best['gap'],attempts


def signed_drift(reference_error,direction,amount,sign=1):
    """fn=e-sign*d*u, fp=e+sign*d*u: recorded probability shifts by sign*d*u."""
    u=np.asarray(direction)
    if amount<0 or np.any(abs(u)>1+1e-12):raise ValueError('Invalid drift pattern')
    fn=reference_error-sign*amount*u;fp=reference_error+sign*amount*u
    if np.any((fn<0)|(fn>1)|(fp<0)|(fp>1)):raise ValueError('Drift leaves probability range')
    return fn,fp
