"""Study 013: quantum-supported magnetic reduced energy. Author: Biswajit Jana."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from mpmath import mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from magnetic_energy import MagneticEnergy, field_coordinates, field_length, scale_diagnostics

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/study-013'; FIG = ROOT/'figures/study-013'
C_TAIL = 2+np.arcsinh(np.sqrt(3))/np.sqrt(3)


def csv_file(name, rows):
    with (OUT/name).open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def save(fig, name):
    fig.savefig(FIG/(name+'.png'), dpi=180, bbox_inches='tight', metadata={'Author': 'Biswajit Jana'})
    path = FIG/(name+'.svg')
    fig.savefig(path, bbox_inches='tight', metadata={'Date': None, 'Creator': 'Biswajit Jana'})
    path.write_text('\n'.join(line.rstrip() for line in path.read_text(encoding='utf-8').splitlines())+'\n', encoding='utf-8', newline='\n')
    plt.close(fig)


def angle_length(nu):
    """Independent 45-digit angle-coordinate arc integral, MMP Appendix D.

    Integrate only the first symmetric half to avoid the far-end angle singularity.
    Positions are rederived from rho(cot theta1-cot theta2)=1.
    """
    with mp.workdps(45):
        n = mp.mpf(str(nu))
        def point(t):
            t2 = mp.acos(mp.cos(t)-n)
            rho = mp.sin(t)*mp.sin(t2)/mp.sin(t2-t)
            z = rho*mp.cot(t)-mp.mpf('.5')
            return z, rho
        def speed(t):
            return mp.sqrt(sum(mp.diff(lambda u: point(u)[i], t)**2 for i in [0, 1]))
        return float(2*mp.quad(speed, [0, mp.acos(n/2)/2, mp.acos(n/2)]))


def adaptive_minimum(s, inner=256, tolerance=2e-11):
    """Endpoint-adapted loss integral for s->0, avoiding subtractive cancellation.

    This audits the nonuniform small-separation limit. QUADPACK estimates are
    numerical diagnostics, not validated bounds or physical uncertainties.
    """
    def loss(y):
        def integrand(n):
            f = float(field_length(n, inner))
            return .5*(f/(y+f))*(1+y/(y+f))
        split = sorted(set([0., 2., *[p/y for p in [.01, .1, 1, 10, 100] if p/y < 2]]))
        vals = [quad(integrand, a, b, epsabs=tolerance/y, epsrel=tolerance,
                     limit=150) for a, b in zip(split[:-1], split[1:])]
        return sum(v[0] for v in vals), sum(v[1] for v in vals)
    # Solve for delta=y-1/s so the first two terms cancel analytically.
    def equation(delta):
        y = 1/s+delta
        return .75*delta-y*loss(y)[0]
    hi = max(100., 20*np.log1p(1/s))
    while equation(hi) < 0: hi *= 2
    delta = brentq(equation, 0., hi, xtol=1e-10, rtol=1e-12)
    y = 1/s+delta; integral, err = loss(y)
    return dict(s=s, y=y, excess_over_one=s*delta,
                logarithmic_coefficient=delta/np.log(1/s),
                integral_error_estimate=err, stationarity_scaled=equation(delta))


def main():
    OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'svg.hashsalt': 'Biswajit-Jana-study013'})
    model = MagneticEnergy(512, 512); shortest = MagneticEnergy(64, 64, shortest=True)
    sweep = []
    for s in np.geomspace(.001, 1000, 241):
        row = model.minimum(float(s)); toy = shortest.minimum(float(s))
        row.update(shortest_y=toy['y'], shortest_energy=toy['energy'])
        sweep.append(row)
    csv_file('equilibrium.csv', sweep)
    asymptotic = model.minimum()['y']
    convergence = []
    for outer, inner in [(64, 64), (128, 128), (256, 256), (512, 512), (1024, 512), (512, 1024)]:
        candidate = MagneticEnergy(outer, inner)
        for s in [.001, .01, .1, 1., 10., 100., 1000.]:
            row = candidate.minimum(s); reference = model.minimum(s)
            row.update(outer=outer, inner=inner, relative_y_difference=abs(row['y']/reference['y']-1))
            convergence.append(row)
    csv_file('quadrature_audit.csv', convergence)
    field_audit = []
    for n in [.001, .01, .1, .3, .7, 1., 1.5, 1.9, 1.99]:
        reference = angle_length(n); production = float(field_length(n))
        field_audit.append(dict(nu=n, prolate_length=production, angle_length=reference,
                                relative_difference=abs(production/reference-1)))
    csv_file('field_length_audit.csv', field_audit)
    lengths = [dict(nu=float(n), length=float(f), nu_times_length=float(n*f))
               for n in np.geomspace(1e-7, 2, 401) for f in [field_length(n)]]
    csv_file('field_lengths.csv', lengths)
    coordinates = []
    for n in [.3, .5, .8, 1.2, 1.6, 1.9]:
        for eta in np.linspace(-1, 1, 601):
            z, rho = field_coordinates(eta, n)
            coordinates.append(dict(nu=n, eta=float(eta), z_over_d=float(z), rho_over_d=float(rho)))
    csv_file('field_lines.csv', coordinates)
    landscapes = []
    for s in [.03, .3, 3., 30.]:
        minimum = model.minimum(s)
        for y in np.geomspace(.5, 300, 501):
            landscapes.append(dict(s=s, y=float(y), energy=model.energy(y, s),
                                   energy_over_abs_min=model.energy(y, s)/abs(minimum['energy'])))
    csv_file('energy_landscapes.csv', landscapes)
    tail = [adaptive_minimum(float(s)) for s in np.geomspace(1e-2, 1e-7, 21)]
    csv_file('small_separation.csv', tail)
    tail_audit = []
    for s in [1e-2, 1e-4, 1e-7]:
        low = adaptive_minimum(s, 256); high = adaptive_minimum(s, 512, 2e-12)
        tail_audit.append(dict(s=s, coefficient_low=low['logarithmic_coefficient'],
            coefficient_high=high['logarithmic_coefficient'],
            coefficient_difference=abs(high['logarithmic_coefficient']-low['logarithmic_coefficient'])))
    csv_file('small_separation_audit.csv', tail_audit)
    cutoff_rows = []
    for cutoff in [.001, .01, .1]:
        reduced = MagneticEnergy(512, 512, cutoff=cutoff)
        for s in np.geomspace(.001, 1000, 61):
            baseline = model.minimum(s); candidate = reduced.minimum(s)
            # Omitted H weight <= cutoff/2. At the baseline root H_cut<=0.
            upper = brentq(lambda y: model.stationarity(y, s)-cutoff/2,
                           baseline['y'], 100*baseline['y'])
            cutoff_rows.append(dict(cutoff=cutoff, s=float(s), y_full=baseline['y'],
                y_cut=candidate['y'], bound_upper=float(upper),
                fractional_shift=candidate['y']/baseline['y']-1))
    csv_file('tail_omission.csv', cutoff_rows)
    scales = []
    for q in [1e3, 1e6, 1e9]:
        for row in sweep:
            scales.append(dict(q=q, g=.3, s=row['s'], **scale_diagnostics(q, .3, row['s'], row['y'], row['energy'])))
    csv_file('scale_diagnostics.csv', scales)
    # Four scientific figures; every plotted series is saved above.
    fig, ax = plt.subplots(1, 2, figsize=(12, 5.1), layout='constrained')
    colors = plt.cm.viridis(np.linspace(.1, .85, 6))
    for color, n in zip(colors, [.3, .5, .8, 1.2, 1.6, 1.9]):
        rr = [r for r in coordinates if r['nu'] == n]
        z = np.array([r['z_over_d'] for r in rr]); rho = np.array([r['rho_over_d'] for r in rr])
        ax[0].plot(z, rho, color=color, label=f'nu={n:g}')
        ax[0].plot(z, -rho, color=color)
    ax[0].scatter([-.5, .5], [0, 0], color='#ac4434', zorder=5)
    ax[0].set(xlabel='Axial position z/d', ylabel='Transverse position / d',
              title='Calculated exterior magnetic field lines', aspect='equal')
    ax[0].legend(fontsize=8, loc='upper left', bbox_to_anchor=(1.02, 1))
    for s in [.03, .3, 3., 30.]:
        rr = [r for r in landscapes if r['s'] == s]
        line, = ax[1].semilogx([r['y'] for r in rr], [r['energy_over_abs_min'] for r in rr], label=f'd/d0={s:g}')
        ax[1].plot(model.minimum(s)['y'], -1, 'o', color=line.get_color())
    ax[1].axhline(0, color='.6', lw=.7)
    ax[1].set(xlabel='Throat travel time / exterior travel time, y', ylabel='Reduced energy / |minimum energy|',
              ylim=(-1.12, 1.8), title='One minimum for each separation')
    ax[1].legend(fontsize=8); ax[1].grid(alpha=.2)
    fig.suptitle('Quantum-supported magnetic model: geometry and energy\nFlat-exterior approximation; not a sky image or a full stability calculation', fontsize=12)
    save(fig, '01-magnetic-support')
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.9), layout='constrained')
    s = np.array([r['s'] for r in sweep]); y = np.array([r['y'] for r in sweep])
    ax[0].loglog(s, s*y, color='#126a88', label='Integrated magnetic field lengths')
    ax[0].loglog(s, s*np.array([r['shortest_y'] for r in sweep]), '--', color='#b26937', label='All exterior paths = d (toy)')
    ax[0].loglog(s, s, ':', color='.4', label='Direct exterior travel time')
    ax[0].set(xlabel='Mouth separation d/d0', ylabel='Throat travel time / d0 (c = 1)', title='The throat is always the slower route')
    ax[0].legend(fontsize=8); ax[0].grid(alpha=.2)
    ax[1].semilogx(s, y, color='#126a88', label='Full field-line integral')
    ax[1].semilogx(s, [r['shortest_y'] for r in sweep], '--', color='#b26937', label='Shortest-path toy')
    ax[1].axhline(asymptotic, color='#126a88', ls=':', label=f'Large-separation limit = {asymptotic:.6f}')
    ax[1].axhline(1, color='.5', ls=':')
    ax[1].set(xlabel='Mouth separation d/d0', ylabel='Travel-time ratio y = pi ell / d',
              ylim=(.8, 8), title='Geometry changes the limiting ratio')
    ax[1].legend(fontsize=8); ax[1].grid(alpha=.2)
    fig.suptitle('Reproducing the MMP reduced-energy equilibrium\n241 separations; unique minimum proved within the stated approximation', fontsize=12)
    save(fig, '02-equilibrium')
    fig, ax = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    for n in [64, 128, 256, 1024]:
        rr = [r for r in convergence if r['outer'] == n]
        ax[0,0].loglog([r['s'] for r in rr], np.maximum([r['relative_y_difference'] for r in rr], 1e-16), 'o-', label=f'Outer order {n}')
    ax[0,0].set(xlabel='d/d0', ylabel='Relative change in y from 512 x 512', title='Quadrature resolution check')
    ax[0,0].legend(fontsize=8)
    ax[0,1].semilogx([r['nu'] for r in lengths], [r['nu_times_length'] for r in lengths], color='#126a88')
    ax[0,1].axhline(C_TAIL, ls='--', color='#b26937', label=f'C = {C_TAIL:.6f}')
    ax[0,1].set(xlabel='Magnetic field-line label nu', ylabel='nu f(nu)', title='Long exterior paths: f(nu) ~ C/nu'); ax[0,1].legend(fontsize=8)
    ax[1,0].plot(np.log10([1/r['s'] for r in tail]), [r['logarithmic_coefficient'] for r in tail], 'o-', color='#126a88', ms=3)
    ax[1,0].axhline(4*C_TAIL/3, ls='--', color='#b26937', label=f'Derived leading coefficient = {4*C_TAIL/3:.6f}')
    ax[1,0].set(xlabel='log10(d0/d)', ylabel='(pi ell/d0 - 1) / [s ln(1/s)]', title='Nonuniform short-separation correction'); ax[1,0].legend(fontsize=8)
    for cutoff in [.001, .01, .1]:
        rr = [r for r in cutoff_rows if r['cutoff'] == cutoff]
        ax[1,1].loglog([r['s'] for r in rr], [r['fractional_shift'] for r in rr], label=f'omit nu < {cutoff:g}')
    ax[1,1].set(xlabel='d/d0', ylabel='Fractional increase of minimizing y', title='Omitted-tail sensitivity (not an error bar)'); ax[1,1].legend(fontsize=8)
    for a in ax.flat: a.grid(alpha=.2)
    fig.suptitle('Numerical checks and the exterior long-path tail\nDeterministic calculations; no observational or Monte Carlo confidence intervals', fontsize=12)
    save(fig, '03-tail-and-numerics')
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.6), layout='constrained')
    for q in [1e3, 1e6, 1e9]:
        rr = [r for r in scales if r['q'] == q]
        for a, key in zip(ax, ['schwarzian_coefficient', 'orbital_omega_ell', 'binding_times_ell']):
            a.loglog([r['s'] for r in rr], [r[key] for r in rr], label=f'q={q:.0e}')
    ax[0].axhline(1, color='.5', ls=':'); ax[1].axhline(1, color='.5', ls=':')
    ax[0].set(ylabel='r_e^3 / (G ell)', title='Weak quantum fluctuation check\nRequires coefficient >> 1')
    ax[1].set(ylabel='Orbital angular frequency x ell', title='Slow exterior motion check\nRequires Omega ell << 1')
    ax[2].set(ylabel='|E_min| ell', title='Binding energy / (1/ell)\nEnergy scale, not packet capacity')
    for a in ax: a.set_xlabel('Mouth separation d/d0'); a.grid(alpha=.2); a.legend(fontsize=8)
    fig.suptitle('An energy minimum alone does not establish a controlled solution\nIllustrative charges and fixed g = 0.3; further quantum and matter constraints remain', fontsize=12)
    save(fig, '04-scale-checks')
    summary = dict(author='Biswajit Jana', date='2026-11-23', source='arXiv:1807.04726v3 Eq.5.47',
        model='Leading reduced energy with point-monopole flat exterior; no full field equations solved',
        outer_order=512, inner_order=512, separation_points=len(sweep), large_separation_ratio=asymptotic,
        long_path_coefficient=C_TAIL, short_separation_log_coefficient=4*C_TAIL/3,
        max_field_length_relative_difference=max(r['relative_difference'] for r in field_audit),
        max_refined_y_relative_difference=max(r['relative_y_difference'] for r in convergence if r['outer']==1024 or r['inner']==1024),
        max_stationarity_residual=max(abs(r['stationarity']) for r in sweep),
        min_curvature=min(r['curvature'] for r in sweep),
        max_tail_coefficient_refinement=max(r['coefficient_difference'] for r in tail_audit),
        source_hashes={p: hashlib.sha256((ROOT/p).read_bytes().replace(b'\r\n', b'\n')).hexdigest()
            for p in ['scripts/magnetic_energy.py', 'scripts/run_magnetic_study.py']})
    (OUT/'summary.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8', newline='\n')
    paths = [*OUT.glob('*.csv'), OUT/'summary.json', *FIG.glob('*.png'), *FIG.glob('*.svg')]
    (OUT/'sha256.json').write_text(json.dumps({p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__': main()
