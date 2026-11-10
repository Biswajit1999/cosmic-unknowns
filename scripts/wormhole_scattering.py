"""Scalar scattering on a prescribed Ellis geometry.

Author: Biswajit Jana
Coordinates X=x/a, K=k*a. This is a probe-field calculation, not an existence proof.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh_tridiagonal


def potential(x, xi=0.0, ell=0):
    """Dimensionless a^2 V; scalar mass is kept in the spectral parameter."""
    x = np.asarray(x)
    return ell*(ell+1)/(1+x*x) + (1-2*xi)/(1+x*x)**2


@dataclass
class Scattering:
    k: float
    xi: float
    ell: int
    domain: float
    transmission: complex
    reflection: complex
    evaluations: int

    @property
    def T(self):
        return float(abs(self.transmission)**2)

    @property
    def R(self):
        return float(abs(self.reflection)**2)

    @property
    def flux_error(self):
        return abs(self.T+self.R-1)


def scatter(k, xi=0.0, ell=0, domain=160.0, rtol=2e-10, atol=2e-12):
    """Right-outgoing Jost solution; incident amplitude normalized on the left.

    Plane-wave boundaries at +/-domain approximate infinity. Flux conservation
    alone does not measure that truncation error; compare domains independently.
    """
    if k <= 0 or domain <= 0 or ell < 0 or int(ell) != ell:
        raise ValueError("Require k>0, domain>0, and integer ell>=0.")
    def rhs(x, y):
        return [y[1], (potential(x, xi, ell)-k*k)*y[0]]
    phase = np.exp(1j*k*domain)
    sol = solve_ivp(rhs, (domain, -domain), [phase, 1j*k*phase],
                    method="DOP853", rtol=rtol, atol=atol,
                    max_step=min(2.0, 0.8/k))
    if not sol.success:
        raise RuntimeError(sol.message)
    u, du = sol.y[:, -1]
    incident = 0.5*np.exp(1j*k*domain)*(u+du/(1j*k))
    reflected = 0.5*np.exp(-1j*k*domain)*(u-du/(1j*k))
    return Scattering(k, xi, ell, domain, 1/incident, reflected/incident, sol.nfev)


def halfline_reflection(k, xi=0.0, domain=160.0, rtol=2e-10):
    """Replace the unobserved right side by its exact frequency-dependent load.

    Z=u'(0)/u(0) inherits the same outgoing boundary at +domain. This demonstrates
    exterior-response equivalence, not a realizable material substitute.
    """
    def rhs(x, y):
        return [y[1], (potential(x, xi)-k*k)*y[0]]
    phase = np.exp(1j*k*domain)
    right = solve_ivp(rhs, (domain, 0), [phase, 1j*k*phase],
                      method="DOP853", rtol=rtol, atol=rtol/100,
                      max_step=min(2.0, 0.8/k))
    if not right.success:
        raise RuntimeError(right.message)
    u0, du0 = right.y[:, -1]
    if abs(u0) < 1e-12:
        raise RuntimeError("Impedance pole: use an inverse load representation.")
    impedance = du0/u0
    left = solve_ivp(rhs, (0, -domain), [1+0j, impedance],
                     method="DOP853", rtol=rtol, atol=rtol/100,
                     max_step=min(2.0, 0.8/k))
    if not left.success:
        raise RuntimeError(left.message)
    u, du = left.y[:, -1]
    inc = 0.5*np.exp(1j*k*domain)*(u+du/(1j*k))
    ref = 0.5*np.exp(-1j*k*domain)*(u-du/(1j*k))
    return ref/inc, impedance


def ground_energy(xi, domain=120.0, spacing=0.025, ell=0):
    """Lowest Dirichlet finite-box eigenvalue of -d_X^2+V, second-order FD.

    Convergence requires BOTH enlarging the box and reducing the spacing.
    Positive finite-box energies are not proof of absence of a full-line bound state.
    """
    if domain <= 0 or spacing <= 0:
        raise ValueError("Require positive domain and spacing.")
    intervals = int(np.ceil(2*domain/spacing))
    dx = 2*domain/intervals
    x = np.linspace(-domain, domain, intervals+1)[1:-1]
    diag = 2/dx**2+potential(x, xi, ell)
    offdiag = np.full(len(x)-1, -1/dx**2)
    energy = eigh_tridiagonal(diag, offdiag, eigvals_only=True,
                              select="i", select_range=(0, 0),
                              lapack_driver="stebz", tol=1e-12)[0]
    return float(energy), float(dx), int(len(x))


def low_energy_coefficient(xi):
    """B in T~4K^2/B^2 for the s wave, away from zero-energy resonances."""
    if xi == 0:
        return float(np.pi)
    if xi > 0:
        s = np.sqrt(2*xi)
        return float(np.sin(np.pi*s)/s)
    s = np.sqrt(-2*xi)
    return float(np.sinh(np.pi*s)/s)
