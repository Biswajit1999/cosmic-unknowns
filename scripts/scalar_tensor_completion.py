"""Study 012: explicit but ghost-supported action. Author: Biswajit Jana.

This is an invertible field redefinition of the known Ellis phantom action.
It does not remove its ghost sign or its mapped linear instability.
"""
import csv
import hashlib
import json
from pathlib import Path
from mpmath import mp
import numpy as np
from scipy.linalg import eigh_tridiagonal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from verify_conformal_study import connection, tensor_check

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT/'results/study-012'; FIG = ROOT/'figures/study-012'


def couplings(phi, beta):
    f = 1+beta*(phi-mp.pi/2)
    if f <= 0: raise ValueError('Conformal transformation requires positive F')
    h = f**-2; z = 2*f**-2+6*beta**2*f**-4
    return f, h, z


def check_equations(x, beta):
    phi = mp.atan(x); f, h, z = couplings(phi, beta)
    metric, gamma = connection(x, mp.mpf('1.1'), beta)
    curvature, rho, pr, pt, pp = tensor_check(x, mp.mpf('1.1'), beta)
    einstein = [-rho, pr, pt, pp]
    h_of_x = lambda y: couplings(mp.atan(y), beta)[1]
    hp = mp.diff(h_of_x, x); hpp = mp.diff(h_of_x, x, 2)
    hessian = [((hpp if a == 1 else 0)-gamma[1][a][a]*hp)/metric[a] for a in range(4)]
    boxh = sum(hessian)
    phip = mp.diff(mp.atan, x); phipp = mp.diff(mp.atan, x, 2)
    kinetic = phip**2/metric[1]
    boxphi = sum(((phipp if a == 1 else 0)-gamma[1][a][a]*phip)/metric[a] for a in range(4))
    metric_residuals = [h*einstein[a]+boxh-hessian[a]+z*kinetic*((1 if a == 1 else 0)-mp.mpf('.5')) for a in range(4)]
    hphi = mp.diff(lambda v: couplings(v, beta)[1], phi)
    zphi = mp.diff(lambda v: couplings(v, beta)[2], phi)
    scalar_residual = hphi*curvature-zphi*kinetic-2*z*boxphi
    einstein_kinetic = z/h-mp.mpf('1.5')*(hphi/h)**2
    return metric_residuals, scalar_residual, einstein_kinetic


def instability_mode(length, spacing):
    """Published regular coupled Ellis master operator, gamma1=0.

    W=-3/(1+X^2)^2, Gonzalez/Guzman/Sarbach arXiv:0806.0608v2 Eq.32.
    Dirichlet finite box, centered second-order differences.
    """
    count = int(round(2*length/spacing))-1
    x = np.linspace(-length+spacing, length-spacing, count)
    diagonal = 2/spacing**2-3/(1+x*x)**2
    off = np.full(count-1, -1/spacing**2)
    values, vectors = eigh_tridiagonal(diagonal, off, select='i', select_range=(0, 0))
    eigenvalue = float(values[0]); mode = vectors[:, 0]
    mode *= np.sign(mode[count//2]); mode /= np.sqrt(np.sum(mode*mode)*spacing)
    applied = diagonal*mode
    applied[:-1] += off*mode[1:]; applied[1:] += off*mode[:-1]
    residual = float(max(abs(applied-eigenvalue*mode)))
    return x, mode, dict(domain=length, spacing=spacing, interior_points=count,
        eigenvalue=eigenvalue, growth_rate=float(np.sqrt(-eigenvalue)), eigen_residual=residual)


def main():
    OUT.mkdir(exist_ok=True, parents=True); FIG.mkdir(exist_ok=True, parents=True)
    checks = []; profiles = []
    with mp.workdps(70):
        for b in ['-0.2', '0', '0.1', '0.2', '0.3']:
            for location in ['-8', '-2', '-0.3', '0', '0.7', '3', '8']:
                residuals, scalar, coefficient = check_equations(mp.mpf(location), mp.mpf(b))
                assert max(abs(r) for r in residuals) < mp.mpf('1e-60')
                assert abs(scalar) < mp.mpf('1e-60')
                assert abs(coefficient-2) < mp.mpf('1e-60')
                checks.append(dict(beta=float(b), X=float(location),
                    metric_residuals=[str(r) for r in residuals], scalar_residual=str(scalar),
                    einstein_kinetic_coefficient=str(coefficient)))
        for beta in ['-0.2', '0', '0.2']:
            b = mp.mpf(beta)
            for x in np.linspace(-6, 6, 601):
                phi = mp.atan(float(x)); f, h, z = couplings(phi, b)
                hphi = mp.diff(lambda p: couplings(p, b)[1], phi)
                profiles.append(dict(beta=float(b), X=float(x), phi=float(phi), F=float(f),
                    curvature_coupling=float(h), kinetic_coupling=float(z),
                    einstein_kinetic=float(z/h-mp.mpf('1.5')*(hphi/h)**2)))
    with (OUT/'couplings.csv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(profiles[0]), lineterminator='\n'); writer.writeheader(); writer.writerows(profiles)
    stability = []
    for length in [10., 20., 30.]:
        for spacing in [.1, .05, .025, .0125]:
            x, mode, row = instability_mode(length, spacing); stability.append(row)
    mode_rows = [dict(X=float(xx), mode=float(yy), potential=float(-3/(1+xx*xx)**2)) for xx, yy in zip(x, mode)]
    for name, rows in [('instability', stability), ('unstable_mode', mode_rows)]:
        with (OUT/(name+'.csv')).open('w', encoding='utf-8', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n'); writer.writeheader(); writer.writerows(rows)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'svg.hashsalt': 'study012'})
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.6))
    for beta, color in [(-.2, '#26728e'), (0., '#313e50'), (.2, '#b45b38')]:
        rr = [r for r in profiles if r['beta'] == beta]
        for ax, column in zip(axes, ['curvature_coupling', 'kinetic_coupling', 'einstein_kinetic']):
            ax.plot([r['X'] for r in rr], [r[column] for r in rr], color=color, label=f'beta={beta:g}')
    axes[0].set(title='Nonminimal curvature coupling', ylabel='h(phi) = F^(-2)')
    axes[1].set(title='Jordan-frame kinetic coefficient', ylabel='z(phi) = 2 F^(-2) + 6 beta^2 F^(-4)')
    axes[2].set(title='Einstein-frame sign is unchanged', ylabel='C_E = z/h − (3/2)(h_phi/h)^2', ylim=(-2.6, 2.6))
    axes[2].axhline(-2, color='gray', ls='--', label='Healthy canonical sign (comparison)')
    axes[2].text(.5, .79, 'C_E = +2: phantom sign', transform=axes[2].transAxes, ha='center', fontsize=9)
    axes[0].legend(fontsize=8); axes[2].legend(fontsize=7, loc='lower center')
    for ax in axes: ax.set_xlabel('X'); ax.grid(alpha=.2)
    fig.suptitle('An explicit supporting action does not cure the ghost\nAction convention: h R + z (grad phi)^2; metric signature (−,+,+,+)')
    fig.tight_layout(rect=(0, 0, 1, .89))
    fig.savefig(FIG/'01-action-and-ghost.png', dpi=180, bbox_inches='tight')
    p = FIG/'01-action-and-ghost.svg'; fig.savefig(p, bbox_inches='tight', metadata={'Date': None, 'Creator': 'Biswajit Jana'})
    p.write_text('\n'.join(s.rstrip() for s in p.read_text(encoding='utf-8').splitlines())+'\n', encoding='utf-8', newline='\n'); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7))
    mask = abs(x) <= 6
    axes[0].plot(x[mask], -3/(1+x[mask]**2)**2, label='Published regular master potential W', color='#26728e')
    axes[0].plot(x[mask], mode[mask], label='Computed normalized unstable eigenfunction', color='#b45b38')
    axes[0].set(xlabel='X', ylabel='Dimensionless W and mode amplitude', title='Coupled metric–scalar instability, Ellis frame')
    axes[0].legend(fontsize=8)
    for length in [10., 20., 30.]:
        rr = [r for r in stability if r['domain'] == length]
        axes[1].plot([r['spacing']**2 for r in rr], [r['growth_rate'] for r in rr], 'o-', label=f'Half-domain L={length:g}')
    axes[1].set(xlabel='Grid spacing squared', ylabel='Growth rate gamma, units c/a', title='Grid refinement and finite-box check')
    axes[1].legend(fontsize=8)
    for ax in axes: ax.grid(alpha=.2)
    fig.suptitle('The explicit action inherits a growing mode\nPublished master equation reproduced; no nonlinear spacetime evolution')
    fig.tight_layout(rect=(0, 0, 1, .9))
    fig.savefig(FIG/'02-coupled-instability.png', dpi=180, bbox_inches='tight')
    p = FIG/'02-coupled-instability.svg'; fig.savefig(p, bbox_inches='tight', metadata={'Date': None, 'Creator': 'Biswajit Jana'})
    p.write_text('\n'.join(s.rstrip() for s in p.read_text(encoding='utf-8').splitlines())+'\n', encoding='utf-8', newline='\n'); plt.close(fig)
    source = ['scripts/scalar_tensor_completion.py', 'scripts/verify_conformal_study.py']
    summary = dict(author='Biswajit Jana', date='2026-11-23', metric_points=len(checks),
        metric_equations_checked=4*len(checks), scalar_equations_checked=len(checks), precision_digits=70,
        maximum_metric_residual=max(float(abs(mp.mpf(v))) for r in checks for v in r['metric_residuals']),
        maximum_scalar_residual=max(float(abs(mp.mpf(r['scalar_residual']))) for r in checks),
        coupling_profile_rows=len(profiles), checks=checks,
        stability_checks=stability, finest_growth_rate=stability[-1]['growth_rate'],
        richardson_growth_rate=(4*stability[-1]['growth_rate']-stability[-2]['growth_rate'])/3,
        maximum_eigen_residual=max(r['eigen_residual'] for r in stability),
        source_hashes={s: hashlib.sha256((ROOT/s).read_bytes().replace(b'\r\n', b'\n')).hexdigest() for s in source},
        conclusion='Explicit transformed phantom action; positive Einstein-frame kinetic coefficient +2. No stability cure.',
        provenance='Reproduces the published regular Ellis coupled-master eigenproblem; transformed-frame mapping is analytic. No nonlinear evolution.')
    (OUT/'verification.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8', newline='\n')
    paths = [*OUT.glob('*.csv'), OUT/'verification.json', *FIG.glob('*.png'), *FIG.glob('*.svg')]
    hashes = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    (OUT/'sha256.json').write_text(json.dumps(hashes, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['checks', 'source_hashes', 'stability_checks']}))


if __name__ == '__main__': main()
