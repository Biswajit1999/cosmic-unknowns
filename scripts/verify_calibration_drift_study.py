"""Audit Study 007 counts, decisions and population bounds. Author: Biswajit Jana."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
from mpmath import iv
from calibration_drift import signed_drift, drift_intervals
from calibrated_identification import continuum_test
from run_calibration_drift_study import evaluate
from transmission_inference import wilson

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/study-007'


def main():
    design = json.loads((OUT / 'design.json').read_text())
    for name, digest in {**design['source_hashes'], **json.loads((OUT / 'sha256.json').read_text())}.items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    certificate = json.loads((OUT / 'population-certificate.json').read_text())
    iv.dps = 70
    gap = iv.mpf(certificate['excluded_gap'])
    for p, lo, hi in zip(certificate['p'], certificate['qlo'], certificate['qhi']):
        lower, upper = iv.mpf(p)-gap, iv.mpf(p)+gap
        assert lo == 0 or bool((iv.mpf(lo)-lower).b <= 0)
        assert hi == 1 or bool((iv.mpf(hi)-upper).a >= 0)
    leaves = sorted(certificate['leaves'])
    assert leaves[0][0] == .1 and leaves[-1][1] == 30.
    for previous, current in zip(leaves, leaves[1:]):
        assert previous[1] == current[0]
    def bound(b, i, q):
        z = iv.mpf(b)*iv.mpf(design['K'][i])
        Q = iv.mpf(q)
        return 2*iv.log((iv.exp(z)-iv.exp(-z))/2)-iv.log(Q/(1-Q))
    for bl, br, i, j in leaves:
        lower = iv.mpf(-16) if i == 0 else bound(bl, i-1, certificate['qhi'][i-1])
        upper = iv.mpf(25) if j == 0 else bound(br, j-1, certificate['qlo'][j-1])
        assert bool((lower-upper).a > 0)
    upper_gap = design['gap']+1e-10
    C, b = [iv.exp(iv.mpf(v)) for v in design['theta']]
    for k, p in zip(design['K'], design['p_ellis']):
        z = b*iv.mpf(k)
        s = (iv.exp(z)-iv.exp(-z))/2
        difference = s*s/(s*s+C)-iv.mpf(p)
        assert bool((difference-iv.mpf(upper_gap)).b <= 0)
        assert bool((difference+iv.mpf(upper_gap)).a >= 0)
    iv.dps = 35
    k, p, q, u = [np.array(design[name]) for name in ['K', 'p_ellis', 'q_comparator', 'direction']]
    f, g = signed_drift(.02, u, design['matching_drift'])
    residual = float(np.max(np.abs((1-f)*q+g*(1-q)-(.02+.96*p))))
    assert residual < 1e-14
    fq,gq=signed_drift(.02,u,design['matching_drift']/2,1)
    fp,gp=signed_drift(.02,u,design['matching_drift']/2,-1)
    midpoint_residual=float(np.max(np.abs((1-fq)*q+gq*(1-q)-((1-fp)*p+gp*(1-p)))))
    assert midpoint_residual < 1e-14
    with (OUT/'population-boundary.csv').open() as stream:
        boundary=list(csv.DictReader(stream))
    assert len(boundary)==101
    for row in boundary:
        d=float(row['actual_drift'])
        assert abs(float(row['minimum_allowance_for_overlap'])+d-design['matching_drift'])<1e-15
        assert abs(float(row['population_margin_matched'])-(design['matching_drift']-2*d))<1e-15
    with (OUT / 'rates.csv').open() as stream:
        rates = list(csv.DictReader(stream))
    experiments = decisions = independent = valid_null_rejections = valid_null_decisions = 0
    for index in range(len(design['budgets'])*len(design['drift_grid'])):
        rr = [r for r in rates if int(r['case']) == index]
        n, m, d = int(rr[0]['N']), int(rr[0]['M']), float(rr[0]['actual_drift'])
        methods = [(f'fixed_{e:g}', e) for e in design['allowances']]+[('matched_bound', d)]
        rng = np.random.default_rng(int(rr[0]['seed']))
        with gzip.open(OUT / rr[0]['file'], 'rt') as stream:
            records = list(csv.DictReader(stream))
        for group, truth, sign in [('comparator_null', q, 1), ('ellis_power', p, -1)]:
            subset = [r for r in records if r['group'] == group]
            assert len(subset) == design['repetitions']
            f, g = signed_drift(.02, u, d, sign)
            for rep, record in enumerate(subset):
                y, fn, fp = [np.array([int(record[f'{prefix}{i}']) for i in range(24)]) for prefix in ['y', 'fn', 'fp']]
                assert int(record['replicate']) == rep
                assert np.array_equal(y, rng.binomial(n, (1-f)*truth+g*(1-truth)))
                assert np.array_equal(fn, rng.binomial(m, .02, size=24))
                assert np.array_equal(fp, rng.binomial(m, .02, size=24))
                statuses = [record[method+'_status'] == 'excluded' for method, eps in sorted(methods, key=lambda v: v[1])]
                assert all(a >= b for a, b in zip(statuses, statuses[1:]))
                experiments += 1
                if rep % 50 == 0:
                    replay = evaluate(y, n, fn, fp, m, k, [eps for method, eps in methods])
                    for method, eps in methods:
                        assert replay[eps]['status'] == record[method+'_status']
                        assert replay[eps]['nodes'] == int(record[method+'_nodes'])
                        decisions += 1
                    # Independently invert and solve the matched bound without the nesting shortcut.
                    fresh = continuum_test(k, *drift_intervals(y, n, fn, fp, m, d))
                    assert (fresh['status'] == 'excluded') == (record['matched_bound_status'] == 'excluded')
                    independent += 1
            for method, eps in methods:
                rate = next(r for r in rr if r['group'] == group and r['method'] == method)
                rejected = sum(r[method+'_status'] == 'excluded' for r in subset)
                assert rejected == int(rate['rejected'])
                assert rejected/len(subset) == float(rate['rate'])
                assert int(rate['bound_contains_truth']) == int(eps >= d)
                assert sum(int(r[method+'_nodes']) for r in subset) == int(rate['nodes'])
                assert sum(r[method+'_status'] == 'undetermined' for r in subset) == int(rate['undetermined'])
                np.testing.assert_allclose(wilson(rejected, len(subset)), [float(rate['wilson_low']), float(rate['wilson_high'])], atol=1e-14)
                if group == 'comparator_null' and eps >= d:
                    valid_null_rejections += rejected
                    valid_null_decisions += len(subset)
    summary = json.loads((OUT / 'summary.json').read_text())
    assert experiments == summary['experiments']
    result = dict(experiments_replayed=experiments, count_entries_replayed=experiments*72,
                  paired_decision_records_replayed=decisions, independent_matched_bound_replays=independent,
                  population_certificate_leaves_verified_at_70_digits=len(leaves),
                  excluded_gap=certificate['excluded_gap'], verified_feasible_gap=upper_gap,
                  exact_exchange_float_residual=residual, artifact_hashes='all verified',
                  midpoint_exchange_float_residual=midpoint_residual,population_boundary_rows_checked=len(boundary),
                  valid_bound_null_rejections=valid_null_rejections, valid_bound_null_decision_records=valid_null_decisions,
                  note='Decision records are paired, not independent; certificates concern saved numerical spectra.')
    (OUT / 'verification.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
