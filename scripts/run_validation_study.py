"""Study 008: validation at fixed incident cost. Author: Biswajit Jana."""
import argparse
import csv
import gzip
import hashlib
import io
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
import mpmath
from calibration_drift import signed_drift
from validation_identification import evaluate
from transmission_inference import wilson

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/study-008'
FIG = ROOT/'figures/study-008'
METHODS = ['stale_reference', 'scalar_transfer', 'direct_validation', 'guarded_validation']
LABELS = ['Stale reference (misspecified)', 'Inferred scalar drift envelope', 'Channel-wise validation', 'Validation + 0.3 pp residual bound']
COLORS = ['#b94e37', '#85649d', '#087d8a', '#c29024']
FRACTIONS = [0, .002, .005, .01, .02, .05, .1, .2, .35, .5, .7, .9, .98]


def json_write(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n', encoding='utf-8', newline='\n')


def csv_write(path, rows):
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def tasks(design):
    output = []
    for budget in design['budgets_per_channel']:
        for fraction in design['validation_fractions']:
            v = int(round(fraction*(budget-2*design['M'])/2))
            output.append(dict(case=len(output), panel='allocation', B=budget, M=design['M'],
                               N=budget-2*design['M']-2*v, V=v, fraction=fraction, residual=0.))
    for residual in [.001, .003, .005]:
        budget = 1600000; v = int(.35*(budget-2*design['M'])/2)
        output.append(dict(case=len(output), panel='stress', B=budget, M=design['M'],
                           N=budget-2*design['M']-2*v, V=v, fraction=.35, residual=residual))
    return output


def case(task, design):
    index, n, m, v = [task[key] for key in ['case', 'N', 'M', 'V']]
    rng = np.random.default_rng(2026112200+index)
    k, p, q, u = [np.array(design[key]) for key in ['K', 'p_ellis', 'q_comparator', 'direction']]
    records = []; repetitions = design['repetitions']; d = design['reference_to_validation_drift']
    for group, truth, sign in [('comparator_null', q, 1), ('ellis_power', p, -1)]:
        fval, gval = signed_drift(.02, u, d, sign)
        fs, gs = signed_drift(.02, u, d+task['residual'], sign)
        for rep in range(repetitions):
            y = rng.binomial(n, (1-fs)*truth+gs*(1-truth))
            rf = rng.binomial(m, .02, 24); rg = rng.binomial(m, .02, 24)
            vf = rng.binomial(v, fval); vg = rng.binomial(v, gval)
            decisions, bounds = evaluate(y, n, rf, rg, m, vf, vg, v, k, design['guard_eta'])
            row = dict(group=group, replicate=rep, drift_low=bounds[0], drift_high=bounds[1])
            for prefix, array in [('y', y), ('rf', rf), ('rg', rg), ('vf', vf), ('vg', vg)]:
                row.update({f'{prefix}{i}': int(value) for i, value in enumerate(array)})
            for method in METHODS:
                row[method+'_status'] = decisions[method]['status']
                row[method+'_nodes'] = decisions[method]['nodes']
            records.append(row)
    name = f'case-{index:02d}.csv.gz'
    with (OUT/name).open('wb') as raw:
        with gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding='utf-8', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(records[0]), lineterminator='\n')
                writer.writeheader(); writer.writerows(records)
    rates = []; coverage = []
    for group in ['comparator_null', 'ellis_power']:
        subset = [r for r in records if r['group'] == group]
        cover = sum(r['drift_low'] <= d <= r['drift_high'] for r in subset)
        cl, ch = wilson(cover, repetitions)
        coverage.append(dict(**task, group=group, repetitions=repetitions, covered=cover,
            coverage=cover/repetitions, wilson_low=float(cl), wilson_high=float(ch),
            upper_median=float(np.median([r['drift_high'] for r in subset])),
            lower_median=float(np.median([r['drift_low'] for r in subset]))))
        for method in METHODS:
            rejected = sum(r[method+'_status'] == 'excluded' for r in subset)
            lo, hi = wilson(rejected, repetitions)
            valid = (method != 'stale_reference' and
                     (task['residual'] <= design['guard_eta'] if method == 'guarded_validation' else task['residual'] == 0))
            rates.append(dict(**task, group=group, method=method, assumed_transfer_valid=int(valid),
                seed=2026112200+index, repetitions=repetitions, rejected=rejected, rate=rejected/repetitions,
                wilson_low=float(lo), wilson_high=float(hi),
                undetermined=sum(r[method+'_status'] == 'undetermined' for r in subset),
                nodes=sum(r[method+'_nodes'] for r in subset), file=name))
    return rates, coverage


def save(fig, name):
    fig.savefig(FIG/(name+'.png'), dpi=180, bbox_inches='tight')
    svg = FIG/(name+'.svg')
    fig.savefig(svg, bbox_inches='tight', metadata={'Date': None, 'Creator': 'Biswajit Jana'})
    svg.write_bytes(('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n').encode())
    plt.close(fig)


def errorcurve(ax, rows, method, xkey, label, color, scale=1):
    subset = sorted([r for r in rows if r['method'] == method], key=lambda r: r[xkey])
    y = np.array([r['rate'] for r in subset])
    lo = np.array([r['wilson_low'] for r in subset]); hi = np.array([r['wilson_high'] for r in subset])
    ax.errorbar([scale*r[xkey] for r in subset], y, yerr=[y-lo, hi-y], marker='o', ms=4,
                capsize=2, color=color, label=label)


def plots(rates, coverage, design):
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'svg.hashsalt': 'study008', 'font.family': 'DejaVu Sans'})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    for col, budget in enumerate(design['budgets_per_channel']):
        for row, group in enumerate(['comparator_null', 'ellis_power']):
            ax = axes[row, col]
            rr = [r for r in rates if r['panel'] == 'allocation' and r['B'] == budget and r['group'] == group]
            for method, label, color in zip(METHODS, LABELS, COLORS):
                errorcurve(ax, rr, method, 'fraction', label, color, 100)
            ax.set_ylim(-.025, 1.025); ax.grid(alpha=.2)
            if row == 0: ax.axhline(.05, ls=':', color='gray'); ax.set_title(f'Total incident cost = {24*budget:,}')
            if col == 0: ax.set_ylabel('Physical-null rejection' if row == 0 else 'Ellis rejection power')
            if row == 1: ax.set_xlabel('Validation share of budget after reference cost (%)')
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=2, fontsize=9)
    fig.suptitle(f'Keeping channel-wise validation versus compressing it to a drift bound\n{design["repetitions"]} paired trials per point; pointwise 95% Wilson intervals')
    fig.tight_layout(rect=(0, .075, 1, .94)); save(fig, '01-fixed-budget-validation')

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for ax, group in zip(axes, ['comparator_null', 'ellis_power']):
        rr = [r for r in rates if r['B'] == 1600000 and r['fraction'] == .35 and r['group'] == group]
        for method, label, color in zip(METHODS, LABELS, COLORS):
            errorcurve(ax, rr, method, 'residual', label, color, 100)
        ax.axvline(.3, ls=':', color=COLORS[3]); ax.axhline(.05, ls=':', color='gray')
        ax.set(title='Physical-null rejection' if group == 'comparator_null' else 'Ellis rejection power',
               xlabel='Actual validation-to-science change (percentage points)', ylim=(-.025, 1.025))
        ax.grid(alpha=.2)
    axes[0].set_ylabel('Rejection fraction')
    fig.legend(*axes[0].get_legend_handles_labels(), loc='lower center', ncol=2, fontsize=9)
    fig.suptitle('Validation only helps when its connection to the science response is justified')
    fig.tight_layout(rect=(0, .12, 1, .93)); save(fig, '02-validation-mismatch-stress')

    chosen = next(r for r in rates if r['B'] == 1600000 and r['fraction'] == .35 and r['panel'] == 'allocation')
    with gzip.open(OUT/chosen['file'], 'rt') as stream:
        records = [r for r in csv.DictReader(stream) if r['group'] == 'ellis_power']
    values = np.array([float(r['drift_high']) for r in records])
    # Fixed edges, including an overflow bin; not fitted to the sampled extrema.
    edges = np.r_[np.linspace(0, 2, 81), 100.]
    counts, edges = np.histogram(100*values, bins=edges)
    hist = []
    for i, count in enumerate(counts):
        lo, hi = wilson(int(count), len(values))
        hist.append(dict(left_pp=edges[i], right_pp=edges[i+1], count=int(count),
                         fraction=count/len(values), wilson_low=float(lo), wilson_high=float(hi)))
    csv_write(OUT/'bound-histogram.csv', hist)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    centers = (edges[:-1]+edges[1:])/2
    heights = np.array([r['fraction'] for r in hist])
    axes[0].bar(centers, heights, width=np.diff(edges)*.95, color='#76a6bc')
    axes[0].errorbar(centers, heights,
        yerr=[heights-np.array([r['wilson_low'] for r in hist]), np.array([r['wilson_high'] for r in hist])-heights],
        fmt='none', color='#284c62', capsize=2)
    axes[0].set(xlabel='Inferred maximum reference-to-validation drift (pp)', ylabel='Fraction of trials',
                title='Distribution of the simultaneous upper bound', xlim=(.6, 1.6))
    for budget in design['budgets_per_channel']:
        rr = sorted([r for r in coverage if r['panel'] == 'allocation' and r['B'] == budget and
                     r['group'] == 'ellis_power' and r['V'] > 0], key=lambda r: r['V'])
        axes[1].plot([r['V'] for r in rr], [100*r['upper_median'] for r in rr], 'o-', label=f'Total cost {24*budget:,}')
    axes[1].axhline(100*design['reference_to_validation_drift'], ls='--', color='gray', label='Actual maximum drift')
    axes[1].set(xscale='log', xlabel='Validation probes per class/channel', ylabel='Median upper bound (pp)',
                title='More validation narrows the bound, but retains drift')
    axes[1].legend(fontsize=8); axes[1].grid(alpha=.2)
    fig.suptitle('Synthetic validation data; histogram bars have pointwise Wilson intervals')
    fig.tight_layout(rect=(0, 0, 1, .94)); save(fig, '03-inferred-drift-bounds')


def hashes():
    paths = sorted([*OUT.glob('*.csv'), *OUT.glob('*.gz'), OUT/'design.json', OUT/'summary.json',
                    *FIG.glob('*.png'), *FIG.glob('*.svg')])
    json_write(OUT/'sha256.json', {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--repetitions', type=int, default=500)
    parser.add_argument('--workers', type=int, default=4); args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    previous = json.loads((ROOT/'results/study-007/design.json').read_text())
    design = {key: previous[key] for key in ['K', 'p_ellis', 'q_comparator', 'direction', 'theta', 'parameter_domain']}
    design.update(author='Biswajit Jana', date='2026-11-22', repetitions=args.repetitions, alpha=.05,
        interval_count=120, tail_probability=.05/240, M=30000, guard_eta=.003,
        reference_to_validation_drift=previous['matching_drift'], reference_error=.02,
        budgets_per_channel=[400000, 1600000], validation_fractions=FRACTIONS,
        seed_rule='2026112200 + case index; draw order Y, Rf, Rg, Vf, Vg; null group then Ellis group',
        data_type='synthetic aggregate binomial counts; no observational measurements',
        allocation='N + 2M + 2V = B per channel; 24 channels; fraction = 2V/(B-2M)',
        direction_description='Same prespecified population direction as Study 007; null sign +1, Ellis sign -1',
        source_hash_normalization='Python sources: CRLF converted to LF before SHA256; data: raw bytes',
        environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                         matplotlib=matplotlib.__version__, mpmath=mpmath.__version__),
        source_hashes={name: hashlib.sha256((ROOT/name).read_bytes().replace(b'\r\n', b'\n') if name.endswith('.py')
                      else (ROOT/name).read_bytes()).hexdigest() for name in
                      ['results/study-007/design.json', 'scripts/calibrated_identification.py',
                       'scripts/calibration_drift.py', 'scripts/validation_identification.py',
                       'scripts/run_validation_study.py']})
    design['cases'] = tasks(design); json_write(OUT/'design.json', design)
    rates = []; coverage = []; start = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(case, task, design) for task in design['cases']]
        for i, future in enumerate(as_completed(futures), 1):
            rr, cc = future.result(); rates.extend(rr); coverage.extend(cc)
            print(f'case {i}/{len(futures)} B={rr[0]["B"]} V={rr[0]["V"]} residual={rr[0]["residual"]} {time.time()-start:.1f}s', flush=True)
    rates.sort(key=lambda r: (r['case'], r['group'], r['method']))
    coverage.sort(key=lambda r: (r['case'], r['group']))
    csv_write(OUT/'rates.csv', rates); csv_write(OUT/'coverage.csv', coverage)
    summary = dict(author='Biswajit Jana', experiments=len(design['cases'])*2*args.repetitions,
        aggregate_count_slots=len(design['cases'])*2*args.repetitions*120,
        decision_records=len(rates)*args.repetitions,
        zero_budget_structural_count_slots=sum(c['V'] == 0 for c in design['cases'])*2*args.repetitions*48,
        undetermined_records=sum(r['undetermined'] for r in rates),
        maximum_valid_null_rate=max(r['rate'] for r in rates if r['group'] == 'comparator_null' and r['assumed_transfer_valid']),
        minimum_drift_interval_coverage=min(r['coverage'] for r in coverage))
    json_write(OUT/'summary.json', summary); plots(rates, coverage, design); hashes()
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
