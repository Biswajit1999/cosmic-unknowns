"""Study 009: identical scalar operators, distinct metrics. Author: Biswajit Jana."""
import csv
import hashlib
import json
import platform
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from conformal_scattering import geometry, stress, throat, scatter, metric_potential
from wormhole_scattering import scatter as plane_scatter

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/study-009'; FIG = ROOT/'figures/study-009'
BETAS = [-.2, 0., .2]


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n', encoding='utf-8', newline='\n')


def write_csv(path, rows):
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def flatten(result):
    result = result.copy()
    for name in ['r', 't']:
        number = result.pop(name)
        result[name+'_real'] = number.real; result[name+'_imag'] = number.imag
    return result


def task(spec):
    k, ell, beta = spec
    return flatten(scatter(k, ell, beta))


def audit_task(spec):
    k, ell, beta = spec
    base = scatter(k, ell, beta); length = base['domain']
    settings = [('double_domain', dict(domain=2*length)),
                ('order8', dict(order=8)), ('order16', dict(order=16)),
                ('tight_tolerance', dict(rtol=2e-12, atol=2e-14)),
                ('reduced_equation', dict(form='reduced'))]
    output = []
    for name, kwargs in settings:
        alternative = scatter(k, ell, beta, **kwargs)
        output.append(dict(k=k, ell=ell, beta=beta, check=name,
            delta_r=abs(base['r']-alternative['r']), delta_t=abs(base['t']-alternative['t']),
            delta_transmission_probability=abs(base['T']-alternative['T']), flux_error=alternative['flux_error'],
            **{key: alternative[key] for key in ['domain', 'order', 'form']},
            r_real=alternative['r'].real, r_imag=alternative['r'].imag,
            t_real=alternative['t'].real, t_imag=alternative['t'].imag))
    if beta == 0:
        plane = plane_scatter(k, ell=ell, domain=length)
        output.append(dict(k=k, ell=ell, beta=beta, check='plane_boundary',
            delta_r=abs(base['r']-plane.reflection), delta_t=abs(base['t']-plane.transmission),
            delta_transmission_probability=abs(base['T']-plane.T), flux_error=plane.flux_error, domain=length, order=0,
            form='prior_plane_solver', r_real=plane.reflection.real, r_imag=plane.reflection.imag,
            t_real=plane.transmission.real, t_imag=plane.transmission.imag))
    return output


def save(fig, name):
    fig.savefig(FIG/(name+'.png'), dpi=180, bbox_inches='tight')
    svg = FIG/(name+'.svg')
    fig.savefig(svg, bbox_inches='tight', metadata={'Date': None, 'Creator': 'Biswajit Jana'})
    svg.write_bytes(('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n').encode())
    plt.close(fig)


def plots(spectra, profiles, throats, audits):
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'svg.hashsalt': 'study009', 'font.family': 'DejaVu Sans'})
    colors = ['#3178a0', '#343b45', '#b85a3c']
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.8))
    for beta, color in zip(BETAS, colors):
        rr = [r for r in profiles if r['beta'] == beta]
        x = [r['X'] for r in rr]
        axes[0, 0].plot(x, [r['radius'] for r in rr], color=color, label=f'beta = {beta:g}')
        axes[0, 1].plot(x, [r['ricci'] for r in rr], color=color)
        axes[1, 0].plot(x, [r['rho_plus_pr'] for r in rr], color=color)
    axes[0, 0].set(xlabel='Tortoise coordinate X', ylabel='Areal radius / a', title='Different spatial geometry')
    axes[0, 0].legend(); axes[0, 1].set(xlabel='X', ylabel='a^2 Ricci scalar', title='Different curvature')
    axes[1, 0].set(xlabel='X', ylabel='8 pi G a^2 (rho + p_r) / c^4', title='Required radial null-energy combination')
    axes[1, 0].axhline(0, color='gray', ls=':')
    axes[1, 1].plot([r['beta'] for r in throats], [r['radius'] for r in throats], color='#147e83')
    axes[1, 1].set(xlabel='Conformal parameter beta', ylabel='Minimum areal radius / a', title='A family of distinct throats')
    for ax in axes.flat: ax.grid(alpha=.2)
    fig.suptitle('Regular metric constructions with identical massless-scalar radial potentials')
    fig.tight_layout(rect=(0, 0, 1, .95)); save(fig, '01-distinct-geometries')

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8))
    mode_colors = ['#116d89', '#b65b36', '#796092']
    styles = ['--', '-', ':']
    for ell, color in enumerate(mode_colors):
        base = sorted([r for r in spectra if r['ell'] == ell and r['beta'] == 0], key=lambda r: r['k'])
        for beta, style in zip(BETAS, styles):
            rr = sorted([r for r in spectra if r['ell'] == ell and r['beta'] == beta], key=lambda r: r['k'])
            k = np.array([r['k'] for r in rr]); t = np.array([complex(r['t_real'], r['t_imag']) for r in rr])
            reflection = np.array([complex(r['r_real'], r['r_imag']) for r in rr])
            axes[0, 0].plot(k, [r['T'] for r in rr], style, color=color, label=f'ell = {ell}' if beta == 0 else None)
            axes[0, 1].plot(k, np.unwrap(np.angle(t)), style, color=color)
            valid = np.array([r['R'] for r in rr]) > 1e-10
            axes[1, 0].plot(k[valid], np.unwrap(np.angle(reflection[valid])), style, color=color)
            if beta != 0:
                residual = [max(abs(complex(r['r_real'], r['r_imag'])-complex(b['r_real'], b['r_imag'])),
                                abs(complex(r['t_real'], r['t_imag'])-complex(b['t_real'], b['t_imag']))) for r, b in zip(rr, base)]
                axes[1, 1].semilogy(k, np.maximum(residual, 1e-16), style, color=color)
    axes[0, 0].set(yscale='log', ylim=(1e-16, 1.2), ylabel='Transmission probability T', title='All three geometries overlap')
    axes[0, 0].legend()
    axes[0, 1].set(ylabel='Unwrapped arg(t), radians', title='Transmission phase in shared convention')
    axes[1, 0].set(ylabel='Unwrapped arg(r), radians', title='Reflection phase; R > 1e-10')
    axes[1, 1].set(ylabel='max(|delta r|, |delta t|)', title='Numerical disagreement with beta = 0')
    for ax in axes.flat: ax.set_xlabel('K = omega a / c'); ax.grid(alpha=.2)
    fig.suptitle('Independent integrations of the original field equation\nLine styles: beta=-0.2 dashed, 0 solid, +0.2 dotted; no counting noise')
    fig.tight_layout(rect=(0, 0, 1, .94)); save(fig, '02-multimode-complex-scattering')

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    axes[0].plot([r['beta'] for r in throats], [r['lapse_sqrt'] for r in throats], label='Throat / right-infinity clock rate')
    axes[0].plot([r['beta'] for r in throats], [r['left_lapse_sqrt'] for r in throats], label='Left / right-infinity clock rate')
    axes[0].set(xlabel='beta', ylabel='Static proper-time rate ratio', title='A known clock supplies different information')
    axes[0].legend(fontsize=8)
    for beta, color in zip(BETAS, colors):
        rr = [r for r in profiles if r['beta'] == beta]
        axes[1].plot([r['X'] for r in rr], [r['massive_delta_U'] for r in rr], color=color, label=f'beta={beta:g}')
    axes[1].set(xlabel='X', ylabel='Delta U = (mu a)^2 (F^2 - 1)', title='A mass term breaks the radial-operator identity')
    axes[1].legend(fontsize=8)
    for ax in axes: ax.grid(alpha=.2)
    fig.suptitle('Conditional discriminants: calibrated clocks or a probe with known mass (mu a = 0.3)')
    fig.tight_layout(rect=(0, 0, 1, .94)); save(fig, '03-noninvariant-observables')

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    for name, label in [('double_domain', 'Double boundary radius'), ('order8', 'Tail order 8 instead of 12'),
                        ('order16', 'Tail order 16'), ('tight_tolerance', '100x tighter tolerance'),
                        ('reduced_equation', 'Independent reduced equation'), ('plane_boundary', 'Prior plane-wave boundaries')]:
        rr = [r for r in audits if r['beta'] == 0 and r['ell'] == 2 and r['check'] == name]
        axes[0].semilogy([r['k'] for r in rr], [max(r['delta_r'], r['delta_t'], 1e-16) for r in rr], 'o-', label=label)
    axes[0].set(xlabel='K', ylabel='Complex-amplitude difference from production', title='ell = 2 boundary and equation checks')
    axes[0].legend(fontsize=7)
    for ell, color in enumerate(mode_colors):
        rr = [r for r in audits if r['beta'] == 0 and r['ell'] == ell and r['check'] == 'plane_boundary']
        axes[1].loglog([r['flux_error'] for r in rr], [max(r['delta_r'], r['delta_t']) for r in rr], 'o', color=color, label=f'ell = {ell}')
    axes[1].set(xlabel='Plane-boundary flux error', ylabel='Difference from tail-corrected amplitudes', title='Flux balance alone does not certify phase')
    axes[1].legend()
    for ax in axes: ax.grid(alpha=.2)
    fig.suptitle('Deterministic convergence checks; differences are not statistical confidence intervals')
    fig.tight_layout(rect=(0, 0, 1, .94)); save(fig, '04-numerical-audit')


def main():
    OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    grid = np.unique(np.r_[np.geomspace(.08, .5, 21), np.linspace(.5, 3.5, 61)])
    specs = [(float(k), ell, beta) for ell in range(3) for beta in BETAS for k in grid]
    spectra = []; audits = []; start = time.time()
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(task, spec) for spec in specs]
        for i, future in enumerate(as_completed(futures), 1):
            spectra.append(future.result())
            if i % 100 == 0: print(f'spectra {i}/{len(specs)} {time.time()-start:.1f}s', flush=True)
        checks = [(k, ell, beta) for ell in range(3) for beta in [0., .2] for k in [.1, .3, .8, 1.7, 3.]]
        for rows in pool.map(audit_task, checks): audits.extend(rows)
    spectra.sort(key=lambda r: (r['ell'], r['beta'], r['k']))
    audits.sort(key=lambda r: (r['ell'], r['beta'], r['check'], r['k']))
    profiles = []
    for beta in BETAS:
        for x in np.linspace(-6, 6, 1201):
            f, fp, fpp, r, rp, rpp = geometry(x, beta)
            rho, pr, pt = stress(x, beta)
            profiles.append(dict(beta=beta, X=float(x), F=float(f), radius=float(r),
                ricci=float(-2/(f*f*(1+x*x)**2)), rho=float(rho), pr=float(pr), pt=float(pt),
                rho_plus_pr=float(rho+pr), massive_delta_U=float(.09*(f*f-1)),
                U0=float(metric_potential(x, beta, 0)), U1=float(metric_potential(x, beta, 1)), U2=float(metric_potential(x, beta, 2))))
    throats = [throat(float(b)) for b in np.linspace(-.2, .2, 101)]
    for name, rows in [('spectra', spectra), ('profiles', profiles), ('throats', throats), ('convergence', audits)]:
        write_csv(OUT/(name+'.csv'), rows)
    differences = []
    for ell in range(3):
        baseline = [r for r in spectra if r['ell'] == ell and r['beta'] == 0]
        for beta in [-.2, .2]:
            rr = [r for r in spectra if r['ell'] == ell and r['beta'] == beta]
            for a, b in zip(rr, baseline):
                differences.append(max(abs(complex(a['r_real'], a['r_imag'])-complex(b['r_real'], b['r_imag'])),
                                       abs(complex(a['t_real'], a['t_imag'])-complex(b['t_real'], b['t_imag']))))
    summary = dict(author='Biswajit Jana', spectral_integrations=len(spectra), convergence_comparisons=len(audits),
        maximum_geometry_amplitude_difference=max(differences), max_flux_error=max(r['flux_error'] for r in spectra),
        maximum_nonplane_convergence_difference=max(max(r['delta_r'], r['delta_t']) for r in audits if r['check'] != 'plane_boundary'),
        maximum_plane_boundary_difference=max(max(r['delta_r'], r['delta_t']) for r in audits if r['check'] == 'plane_boundary'))
    sources = ['scripts/conformal_scattering.py', 'scripts/run_conformal_study.py', 'scripts/wormhole_scattering.py']
    design = dict(author='Biswajit Jana', date='2026-11-22', betas=BETAS, modes=[0,1,2], K=grid.tolist(),
        field='massless minimally coupled scalar; theorem also covers fixed curvature coupling',
        phase_convention='Shared tortoise coordinate X; radiative amplitude u=r phi; not a fixed-areal-radius experiment',
        geometry='ds^2/a^2 = F^2[-d_tau^2+dX^2+(1+X^2)dOmega^2]; F=1+beta(atan X-pi/2)',
        baseline=dict(domain='max(80,30/K)', tail_order=12, rtol=2e-10, atol=2e-12),
        environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, matplotlib=matplotlib.__version__),
        source_hash_normalization='Python source CRLF normalized to LF',
        source_hashes={name: hashlib.sha256((ROOT/name).read_bytes().replace(b'\r\n', b'\n')).hexdigest() for name in sources})
    write_json(OUT/'design.json', design); write_json(OUT/'summary.json', summary)
    plots(spectra, profiles, throats, audits)
    paths = sorted([*OUT.glob('*.csv'), OUT/'design.json', OUT/'summary.json', *FIG.glob('*.png'), *FIG.glob('*.svg')])
    write_json(OUT/'sha256.json', {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
