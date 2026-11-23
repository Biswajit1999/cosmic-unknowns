"""Study 011: fitting nuisance parameters and checking parity. Author: Biswajit Jana."""
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
from matplotlib.ticker import NullFormatter
from conformal_scattering import scatter
from massive_conformal import scatter_massive
from nuisance_identification import NullSpectrum, STAGES, BOUNDS, centered_parameters

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT/'results/study-011'; FIG = ROOT/'figures/study-011'


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8', newline='\n')


def write_csv(path, rows):
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def template_task(spec):
    q, ell = spec
    a = scatter(q, ell=ell, form='reduced')
    return dict(q=q, ell=ell, T=a['T'], flux_error=a['flux_error'])


def direct_task(row):
    q = np.sqrt((row['scale']*row['k'])**2-row['mass']**2)
    t = scatter(q, ell=row['ell'], form='reduced')['T']
    p = row['false_positive']+(1-row['false_positive']-row['false_negative'])*t
    return dict(**row, direct_T=t, direct_probability=p, direct_residual=p-row['truth'],
                interpolation_difference=p-row['interpolated_probability'])


def parity_task(spec):
    k, ell, lam = spec
    beta, mass = centered_parameters(lam, .3)
    a = scatter_massive(k, beta, mass, ell)
    return dict(k=k, ell=ell, lam=lam, beta=beta, mass=mass, T=a['T'], flux_error=a['flux_error'])


def derivative_task(spec):
    k, ell, h = spec
    bp = scatter_massive(k, h, .3, ell)['T']; bm = scatter_massive(k, -h, .3, ell)['T']
    mp = scatter_massive(k, 0., .3+h, ell)['T']; mm = scatter_massive(k, 0., .3-h, ell)['T']
    db = (bp-bm)/(2*h); dm = (mp-mm)/(2*h)
    return dict(k=k, ell=ell, step=h, T_beta_plus=bp, T_beta_minus=bm,
        T_mass_plus=mp, T_mass_minus=mm, dT_dbeta=db, dT_dmass=dm,
        predicted_dT_dbeta=-np.pi*.3*dm/2, identity_residual=db+np.pi*.3*dm/2)


def save(fig, name):
    fig.savefig(FIG/(name+'.png'), dpi=180, bbox_inches='tight')
    p = FIG/(name+'.svg')
    fig.savefig(p, bbox_inches='tight', metadata={'Date': None, 'Creator': 'Biswajit Jana'})
    p.write_text('\n'.join(s.rstrip() for s in p.read_text(encoding='utf-8').splitlines())+'\n', encoding='utf-8', newline='\n')
    plt.close(fig)


def plots(fits, predictions, parity, derivatives):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'svg.hashsalt': 'study011'})
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    labels = ['Fixed', 'Mass free', '+ scale', '+ response']
    for beta, color in [(-.2, '#26728e'), (.2, '#b45b38')]:
        for obs, style in [('ell0', '--'), ('joint', '-')]:
            rr = [next(r for r in fits if r['beta_truth'] == beta and r['observables'] == obs and r['stage'] == s) for s in STAGES]
            for ax, key in zip(axes, ['direct_rms', 'direct_maximum_residual']):
                ax.semilogy(range(4), [r[key] for r in rr], style+'o', color=color, label=f'beta={beta:g}, {obs}')
    axes[0].set(ylabel='RMS probability residual', title='Objective used for the fit')
    axes[1].set(ylabel='Largest absolute channel residual', title='Achieved residual; not a minimax optimum')
    for ax in axes: ax.set_xticks(range(4), labels); ax.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    fig.suptitle('The fixed-parameter spectral gap shrinks when nuisance parameters vary\nBounds: scale 0.95–1.05; mass 0.10–0.44; each response error 0–0.03')
    fig.tight_layout(rect=(0, 0, 1, .89)); save(fig, '01-nuisance-residuals')

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.7), sharex='col')
    for col, beta in enumerate([-.2, .2]):
        for ell, color in [(0, '#26728e'), (1, '#805c9a')]:
            rr = [r for r in predictions if r['beta_truth'] == beta and r['observables'] == 'joint'
                  and r['stage'] == 'mass_scale_response' and r['ell'] == ell]
            axes[0, col].plot([r['k'] for r in rr], [r['truth'] for r in rr], '-', color=color, label=f'truth ell={ell}')
            axes[0, col].plot([r['k'] for r in rr], [r['direct_probability'] for r in rr], 'x', color=color, label=f'fitted null ell={ell}')
            axes[1, col].plot([r['k'] for r in rr], [r['direct_residual'] for r in rr], 'o-', ms=3, color=color, label=f'ell={ell}')
        axes[0, col].set(title=f'Truth beta={beta:g}; fitted geometry beta=0', ylabel='Recorded probability')
        axes[0, col].legend(fontsize=8)
        axes[1, col].axhline(0, color='gray', lw=.7)
        axes[1, col].set(xlabel='Nominal K', ylabel='Direct ODE fit − truth')
    for ax in axes.flat: ax.grid(alpha=.2)
    fig.suptitle('A fitted alternative can imitate most of the massive-probe difference\n21 frequencies per mode; common mass, scale and response parameters')
    fig.tight_layout(rect=(0, 0, 1, .92)); save(fig, '02-fitted-spectra')

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.9))
    for ell, color in [(0, '#26728e'), (1, '#b45b38')]:
        rr = [r for r in derivatives if r['ell'] == ell and r['step'] == .0005]
        axes[0].plot([r['k'] for r in rr], [r['dT_dbeta'] for r in rr], 'o-', color=color, label=f'dT/dbeta, ell={ell}')
        axes[0].plot([r['k'] for r in rr], [r['predicted_dT_dbeta'] for r in rr], 'x--', color=color, label=f'−pi nu/2 × dT/dnu, ell={ell}')
        baseline = next(r['T'] for r in parity if r['ell'] == ell and r['k'] == .65 and r['lam'] == 0)
        rr = [r for r in parity if r['ell'] == ell and r['k'] == .65 and r['lam'] > 0]
        lam = np.array([r['lam'] for r in rr]); change = np.array([abs(r['T']-baseline) for r in rr])
        axes[1].loglog(lam, change, 'o-', color=color, label=f'ell={ell}, K=0.65')
        axes[1].loglog(lam, change[0]*(lam/lam[0])**2, ':', color=color)
    axes[0].set(xlabel='K', ylabel='Transmission derivative', title='Geometry and mass have collinear local effects')
    axes[1].set(xlabel='Centered asymmetry |lambda|', ylabel='|T(lambda) − T(0)| at fixed centered mass', title='Leading separation is second order near symmetry')
    axes[1].set_xticks([.005, .01, .02, .04, .08], ['0.005', '0.01', '0.02', '0.04', '0.08'])
    axes[1].xaxis.set_minor_formatter(NullFormatter())
    for ax in axes: ax.legend(fontsize=7); ax.grid(alpha=.2)
    fig.suptitle('An analytic reason the nuisance fit succeeds\nFixed centered mass = 0.3; dotted guides are proportional to lambda squared')
    fig.tight_layout(rect=(0, 0, 1, .89)); save(fig, '03-local-confounding')


def main():
    OUT.mkdir(exist_ok=True, parents=True); FIG.mkdir(exist_ok=True, parents=True)
    templates = []; start = time.time(); grid = np.linspace(.20, 2.12, 1001)
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(template_task, (float(q), ell)) for ell in [0, 1] for q in grid]
        for i, future in enumerate(as_completed(futures), 1):
            templates.append(future.result())
            if i % 400 == 0: print(f'templates {i}/2002, {time.time()-start:.1f}s', flush=True)
        templates.sort(key=lambda r: (r['ell'], r['q']))
        model = NullSpectrum(templates)
        indices = np.linspace(0, 999, 121).astype(int)
        checks = list(pool.map(template_task, [(float((grid[i]+grid[i+1])/2), ell) for ell in [0, 1] for i in indices]))
        for r in checks:
            r['interpolated_T'] = float(model.functions[r['ell']](r['q']))
            r['error'] = r['interpolated_T']-r['T']
        with (ROOT/'results/study-010/spectra.csv').open() as stream: parent = list(csv.DictReader(stream))
        all_k = sorted(set(float(r['k']) for r in parent)); chosen = np.array(all_k)[::4].tolist()
        fits = []; starts = []; requested = []
        for beta in [-.2, .2]:
            for obs in ['ell0', 'joint']:
                truth_rows = [r for r in parent if float(r['beta']) == beta and float(r['mass']) == .3 and float(r['k']) in chosen and (obs == 'joint' or int(r['ell']) == 0)]
                k = np.array([float(r['k']) for r in truth_rows]); ell = np.array([int(r['ell']) for r in truth_rows]); target = np.array([float(r['T']) for r in truth_rows])
                for stage in STAGES:
                    parameters, residual, runs = model.fit(k, ell, target, stage)
                    p = model.predict(k, ell, parameters)
                    names = ['scale', 'mass', 'false_positive', 'false_negative']; pars = dict(zip(names, parameters.tolist()))
                    fits.append(dict(beta_truth=beta, observables=obs, stage=stage, channels=len(k), **pars,
                        interpolated_rms=float(np.sqrt(np.mean(residual**2))), interpolated_maximum_residual=float(max(abs(residual)))))
                    starts.extend(dict(beta_truth=beta, observables=obs, stage=stage, **r) for r in runs)
                    requested.extend(dict(beta_truth=beta, observables=obs, stage=stage, k=float(ki), ell=int(li), truth=float(yi),
                        interpolated_probability=float(pi), **pars) for ki, li, yi, pi in zip(k, ell, target, p))
        predictions = list(pool.map(direct_task, requested))
        for row in fits:
            rr = [r for r in predictions if all(r[key] == row[key] for key in ['beta_truth', 'observables', 'stage'])]
            errors = np.array([r['direct_residual'] for r in rr])
            row['direct_rms'] = float(np.sqrt(np.mean(errors**2))); row['direct_maximum_residual'] = float(max(abs(errors)))
        lambdas = [-.08, -.04, -.02, -.01, -.005, 0., .005, .01, .02, .04, .08]
        parity = list(pool.map(parity_task, [(k, ell, lam) for ell in [0, 1] for k in [.65, .9, 1.4] for lam in lambdas]))
        derivatives = list(pool.map(derivative_task, [(k, ell, h) for ell in [0, 1] for k in [.55, .75, 1., 1.6] for h in [.001, .0005]]))
    parity_pairs = []
    for r in parity:
        if r['lam'] > 0:
            opposite = next(s for s in parity if s['k'] == r['k'] and s['ell'] == r['ell'] and s['lam'] == -r['lam'])
            parity_pairs.append(abs(r['T']-opposite['T']))
    summary = dict(author='Biswajit Jana', template_integrations=len(templates), heldout_integrations=len(checks),
        fitted_cases=len(fits), direct_fit_integrations=len(predictions), parity_integrations=len(parity),
        derivative_integrations=4*len(derivatives),
        maximum_heldout_interpolation_error=max(abs(r['error']) for r in checks),
        maximum_fit_interpolation_error=max(abs(r['interpolation_difference']) for r in predictions),
        maximum_centered_parity_difference=max(parity_pairs),
        maximum_derivative_identity_residual=max(abs(r['identity_residual']) for r in derivatives if r['step'] == .0005),
        fits=fits)
    for name, rows in [('templates', templates), ('heldout', checks), ('fits', fits), ('optimizer_starts', starts),
                       ('predictions', predictions), ('parity', parity), ('derivatives', derivatives)]:
        write_csv(OUT/(name+'.csv'), rows)
    sources = ['scripts/nuisance_identification.py', 'scripts/run_nuisance_study.py', 'scripts/massive_conformal.py', 'scripts/conformal_scattering.py']
    design = dict(author='Biswajit Jana', date='2026-11-23', bounds=BOUNDS, channels_per_mode=21, K=chosen,
        truth_betas=[-.2, .2], truth_mass=.3, response='b+(1-b-c)*T; common across all frequencies and retained modes',
        objective='Unweighted mean squared probability residual; multistart local optimization; no likelihood or global certificate',
        data_kind='Noise-free model probabilities; no confidence intervals or astronomical observations',
        environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
        parent_hashes={'results/study-010/spectra.csv': hashlib.sha256((ROOT/'results/study-010/spectra.csv').read_bytes()).hexdigest()},
        source_hash_normalization='CRLF normalized to LF',
        source_hashes={s: hashlib.sha256((ROOT/s).read_bytes().replace(b'\r\n', b'\n')).hexdigest() for s in sources})
    write_json(OUT/'design.json', design); write_json(OUT/'summary.json', summary)
    plots(fits, predictions, parity, derivatives)
    paths = [*OUT.glob('*.csv'), OUT/'design.json', OUT/'summary.json', *FIG.glob('*.png'), *FIG.glob('*.svg')]
    write_json(OUT/'sha256.json', {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)})
    print(json.dumps(summary), flush=True)


if __name__ == '__main__': main()
