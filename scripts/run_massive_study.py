"""Study 010 data, scientific illustrations and audit. Author: Biswajit Jana."""
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
from massive_conformal import scatter_massive, scalar_source_obstruction
from conformal_scattering import scatter, stress

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/study-010'; FIG = ROOT/'figures/study-010'
BETAS = [-.2, 0., .2]
COLORS = ['#26728e', '#313e50', '#b45b38']


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n', encoding='utf-8', newline='\n')


def write_csv(path, rows):
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def flatten(row):
    row = row.copy()
    for name in ['r', 't', 't_flux']:
        value = row.pop(name)
        row[name+'_real'] = value.real; row[name+'_imag'] = value.imag
    return row


def task(spec):
    k, beta, mass, ell = spec
    return flatten(scatter_massive(k, beta, mass, ell))


def audit(spec):
    k, beta, ell = spec
    base = scatter_massive(k, beta, .3, ell)
    settings = [('double_domain', dict(domain=2*base['domain'])),
                ('order8', dict(order=8)), ('order16', dict(order=16)),
                ('tight_tolerance', dict(rtol=2e-12, atol=2e-14)),
                ('original_equation', dict(form='radial')),
                ('reverse_incidence', dict(incidence='left'))]
    rows = []
    for check, kwargs in settings:
        alt = scatter_massive(k, beta, .3, ell, **kwargs)
        rows.append(dict(check=check, delta_t_flux=abs(base['t_flux']-alt['t_flux']),
            delta_transmission_probability=abs(base['T']-alt['T']),
            # Reverse reflection phase need not match for an asymmetric potential.
            delta_r=abs(base['r']-alt['r']) if check != 'reverse_incidence' else '',
            **flatten(alt)))
    return rows


def save(fig, name):
    fig.savefig(FIG/(name+'.png'), dpi=180, bbox_inches='tight', facecolor='white')
    p = FIG/(name+'.svg')
    fig.savefig(p, bbox_inches='tight', metadata={'Date': None, 'Creator': 'Biswajit Jana'})
    p.write_text('\n'.join(line.rstrip() for line in p.read_text(encoding='utf-8').splitlines())+'\n', encoding='utf-8', newline='\n')
    plt.close(fig)


def select(rows, beta, mass, ell=0):
    return sorted([r for r in rows if r['beta'] == beta and r['mass'] == mass and r['ell'] == ell], key=lambda r: r['k'])


def plots(spectra, audits, profiles, thresholds, meridian):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
        'axes.spines.top': False, 'axes.spines.right': False, 'svg.hashsalt': 'study010'})
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.8), sharex=True)
    for ell in [0, 1]:
        base = select(spectra, 0., .3, ell)
        for beta, color in zip(BETAS, COLORS):
            for mass, style in [(0., '--'), (.3, '-')]:
                rr = select(spectra, beta, mass, ell)
                axes[0, ell].plot([r['k'] for r in rr], [r['T'] for r in rr], style, color=color,
                    label=fr'$\beta={beta:g}$' if mass else None)
            rr = select(spectra, beta, .3, ell)
            axes[1, ell].plot([r['k'] for r in rr], [r['T']-b['T'] for r, b in zip(rr, base)], color=color)
        axes[0, ell].set(title=fr'Angular mode $\ell={ell}$', ylabel='Transmitted flux fraction T', ylim=(0, 1.02))
        axes[0, ell].legend()
        axes[1, ell].axhline(0, color='gray', lw=.7)
        axes[1, ell].set(xlabel=r'$K=\omega a/c$', ylabel=r'$T(\beta)-T(0)$, massive probe')
    for ax in axes.flat: ax.grid(alpha=.18)
    fig.suptitle('A known mass breaks the shared massless response\nDashed: massless (overlapping); solid: massive, mu a = 0.3', y=.99)
    fig.tight_layout(rect=(0, 0, 1, .93)); save(fig, '01-massive-discrimination')

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7))
    beta = np.array([r['beta'] for r in thresholds]); left = np.array([r['left_threshold'] for r in thresholds])
    axes[0].plot(beta, left, label='Left-end threshold', color=COLORS[0])
    axes[0].axhline(.3, label='Right-end threshold', color=COLORS[2])
    axes[0].fill_between(beta, np.maximum(.3, left), .8, color='#deedf3', label='Both channels open')
    axes[0].axhline(.52, color='#66727e', ls=':', label='Production grid minimum')
    axes[0].set(xlabel=r'$\beta$', ylabel='K', ylim=(0, .8), title='Propagation thresholds for mu a = 0.3')
    axes[0].legend(fontsize=8)
    for beta, color in zip(BETAS, COLORS):
        rr = [r for r in profiles if r['beta'] == beta]
        axes[1].plot([r['X'] for r in rr], [r['rho_plus_pt'] for r in rr], color=color, label=fr'$\beta={beta:g}$')
    axes[1].axhline(0, color='gray', ls=':')
    axes[1].set(xlabel='Tortoise coordinate X', ylabel=r'$8\pi G a^2(\rho+p_t)/c^4$', title='A single static radial scalar requires zero')
    axes[1].legend(fontsize=8)
    for ax in axes: ax.grid(alpha=.18)
    fig.suptitle('A distinguishable mathematical geometry still needs a viable source')
    fig.tight_layout(rect=(0, 0, 1, .93)); save(fig, '02-thresholds-and-matter')

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    checks = ['double_domain', 'order8', 'order16', 'tight_tolerance', 'original_equation', 'reverse_incidence']
    labels = ['2x boundary', 'Tail 8', 'Tail 16', '100x tighter', 'Original field', 'Reciprocity']
    for i, (check, label) in enumerate(zip(checks, labels)):
        rr = [r for r in audits if r['check'] == check]
        values = [max(r['delta_t_flux'], float(r['delta_r'] or 0), 1e-16) for r in rr]
        axes[0].scatter(np.full(len(values), i), values, s=22, alpha=.6)
    axes[0].set(yscale='log', ylabel='Absolute complex-amplitude difference', title='18 audit problems per check')
    axes[0].set_xticks(range(len(labels)), labels, rotation=25, ha='right')
    for beta, color in zip(BETAS, COLORS):
        rr = select(spectra, beta, .3)
        axes[1].plot([r['k'] for r in rr], [r['T'] for r in rr], color=color, label=fr'$\beta={beta:g}$')
        axes[1].plot([r['k'] for r in rr], [r['t_real']**2+r['t_imag']**2 for r in rr], '--', color=color)
    axes[1].set(xlabel='K', ylabel='Transmission', title='Solid: flux T; dashed: unweighted |t|^2')
    axes[1].legend(fontsize=8)
    for ax in axes: ax.grid(alpha=.18)
    fig.suptitle('Numerical validation and the unequal-speed flux correction\nConvergence differences are not statistical confidence intervals')
    fig.tight_layout(rect=(0, 0, 1, .91)); save(fig, '03-numerical-audit')

    # Geometric illustration: only the beta=0 equatorial spatial slice embeds this way.
    fig = plt.figure(figsize=(13.5, 5.9), facecolor='white')
    ax = fig.add_axes([.015, .17, .49, .69], projection='3d')
    radius = np.array([r['radius'] for r in meridian]); z = np.array([r['embedding_z'] for r in meridian])
    angle = np.linspace(0, 2*np.pi, 91)
    xx = radius[:, None]*np.cos(angle); yy = radius[:, None]*np.sin(angle)
    zz = np.broadcast_to(z[:, None], xx.shape)
    ax.plot_surface(xx, yy, zz, color='#73a9c2', alpha=.62, linewidth=0, rcount=81, ccount=61)
    ax.plot_wireframe(xx, yy, zz, color='#204c68', alpha=.35, rstride=10, cstride=9, linewidth=.5)
    ax.plot(np.cos(angle), np.sin(angle), np.zeros_like(angle), color='#ba6038', lw=2.4)
    ax.view_init(elev=20, azim=-52); ax.set_axis_off()
    ax.set_box_aspect((2*max(radius), 2*max(radius), max(z)-min(z)))
    ax.set_title('Calculated Ellis spatial embedding', fontsize=13, pad=0)
    fig.text(.055, .145, r'$r/a=\sqrt{1+X^2}$     $z/a=\mathrm{asinh}(X)$', fontsize=12, color='#263d50')
    fig.text(.055, .09, 'Orange ring: throat at r = a\nSpatial slice, not a telescope view or a physical extra dimension.', fontsize=9, color='#596977')
    for j, mass in enumerate([0., .3]):
        bx = fig.add_axes([.59, .56-j*.34, .365, .235])
        for beta, color in zip(BETAS, COLORS):
            rr = select(spectra, beta, mass)
            bx.plot([r['k'] for r in rr], [r['T'] for r in rr], color=color,
                ls=['--', '-', ':'][BETAS.index(beta)] if mass == 0 else '-', label=fr'$\beta={beta:g}$')
        bx.set(ylabel='Transmitted flux', ylim=(0, 1.03), title='Massless: three geometries overlap' if mass == 0 else 'Known massive probe: responses separate')
        bx.grid(alpha=.15)
        if mass == 0: bx.legend(fontsize=8, loc='lower right')
        else: bx.set_xlabel(r'$K=\omega a/c$;  massive panel: $\mu a=0.3$')
    fig.suptitle('COSMIC UNKNOWNS  /  Geometry, signals and what can be inferred', x=.04, ha='left', fontsize=17, color='#183d58', y=.99)
    fig.text(.59, .065, 'Fixed-background model calculations • Biswajit Jana\nNo astronomical detection or viable supporting matter established.', fontsize=9, color='#596977')
    save(fig, '00-research-overview')


def main():
    OUT.mkdir(exist_ok=True, parents=True); FIG.mkdir(exist_ok=True, parents=True)
    kgrid = np.unique(np.r_[np.linspace(.52, .9, 39), np.linspace(.9, 2., 43)])
    specs = [(float(k), b, m, ell) for ell in [0, 1] for m in [0., .3] for b in BETAS for k in kgrid]
    spectra = []; audits = []; start = time.time()
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(task, spec) for spec in specs]
        for i, future in enumerate(as_completed(futures), 1):
            spectra.append(future.result())
            if i % 100 == 0: print(f'spectra {i}/{len(specs)}, {time.time()-start:.1f}s', flush=True)
        for rows in pool.map(audit, [(k, b, ell) for ell in [0, 1] for b in BETAS for k in [.54, .85, 1.6]]):
            audits.extend(rows)
    spectra.sort(key=lambda r: (r['ell'], r['mass'], r['beta'], r['k']))
    profiles = []
    for b in BETAS:
        for x in np.linspace(-6, 6, 601):
            rho, pr, pt = stress(x, b)
            profiles.append(dict(beta=b, X=float(x), rho_plus_pt=float(rho+pt),
                analytic_obstruction=float(scalar_source_obstruction(x, b))))
    thresholds = [dict(beta=float(b), right_threshold=.3, left_threshold=float(.3*(1-np.pi*b))) for b in np.linspace(-.2, .2, 121)]
    meridian = [dict(X=float(x), radius=float(np.sqrt(1+x*x)), embedding_z=float(np.arcsinh(x))) for x in np.linspace(-3, 3, 161)]
    for name, rows in [('spectra', spectra), ('convergence', audits), ('source_profiles', profiles), ('thresholds', thresholds), ('embedding_meridian', meridian)]:
        write_csv(OUT/(name+'.csv'), rows)
    gaps = []
    for ell in [0, 1]:
        ref = select(spectra, 0., .3, ell)
        for b in [-.2, .2]:
            rr = select(spectra, b, .3, ell)
            gaps.extend(dict(ell=ell, beta=b, k=r['k'], delta_T=r['T']-s['T']) for r, s in zip(rr, ref))
    max_gap = max(gaps, key=lambda r: abs(r['delta_T']))
    summary = dict(author='Biswajit Jana', spectral_integrations=len(spectra), audit_comparisons=len(audits),
        maximum_flux_error=max(r['flux_error'] for r in spectra),
        maximum_tail_current_error=max(r['tail_current_error'] for r in spectra),
        maximum_audit_amplitude_difference=max(max(r['delta_t_flux'], float(r['delta_r'] or 0)) for r in audits),
        maximum_transmission_probability_gap=max_gap,
        maximum_reciprocity_difference=max(r['delta_t_flux'] for r in audits if r['check'] == 'reverse_incidence'),
        source_identity_error=max(abs(r['rho_plus_pt']-r['analytic_obstruction']) for r in profiles))
    sources = ['scripts/massive_conformal.py', 'scripts/run_massive_study.py', 'scripts/conformal_scattering.py']
    design = dict(author='Biswajit Jana', date='2026-11-22', betas=BETAS, masses=[0., .3], modes=[0, 1], K=kgrid.tolist(),
        domain='max(100,40/min(q_right,q_left))', order=12, rtol=2e-10, atol=2e-12,
        data_kind='deterministic model data; no observational measurements or statistical intervals',
        boundary_phase='exp(i(q*y-eta*log(2*q*y))) sum a_n/y^n; y=abs(X)',
        embedding=dict(beta=0, azimuthal_points=91, meaning='Equatorial spatial embedding, not a sky image'),
        environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, matplotlib=matplotlib.__version__),
        source_hash_normalization='CRLF normalized to LF',
        source_hashes={s: hashlib.sha256((ROOT/s).read_bytes().replace(b'\r\n', b'\n')).hexdigest() for s in sources})
    write_json(OUT/'summary.json', summary); write_json(OUT/'design.json', design)
    plots(spectra, audits, profiles, thresholds, meridian)
    paths = [*OUT.glob('*.csv'), OUT/'summary.json', OUT/'design.json', *FIG.glob('*.png'), *FIG.glob('*.svg')]
    write_json(OUT/'sha256.json', {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)})
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
