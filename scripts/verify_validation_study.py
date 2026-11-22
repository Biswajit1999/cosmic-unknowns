"""Replay Study 008 and verify interval implications. Author: Biswajit Jana."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
from calibration_drift import signed_drift
from validation_identification import evaluate, validation_regions
from transmission_inference import wilson

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/study-008'
METHODS = ['stale_reference', 'scalar_transfer', 'direct_validation', 'guarded_validation']


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    design = json.loads((OUT/'design.json').read_text())
    for name, digest in json.loads((OUT/'sha256.json').read_text()).items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    for name, digest in design['source_hashes'].items():
        data = (ROOT/name).read_bytes()
        if name.endswith('.py'):
            data = data.replace(b'\r\n', b'\n')
        assert hashlib.sha256(data).hexdigest() == digest, name
    rates = read(OUT/'rates.csv'); coverage = read(OUT/'coverage.csv')
    k, p, q, u = [np.array(design[key]) for key in ['K', 'p_ellis', 'q_comparator', 'direction']]
    d = design['reference_to_validation_drift']
    experiments = audited = nesting_checks = valid_null_rejects = valid_null_records = 0
    for case in design['cases']:
        index, n, m, v = [case[key] for key in ['case', 'N', 'M', 'V']]
        assert n+2*m+2*v == case['B']
        assert abs(2*v/(case['B']-2*m)-case['fraction']) < 1e-14
        rng = np.random.default_rng(2026112200+index)
        with gzip.open(OUT/f'case-{index:02d}.csv.gz', 'rt') as stream:
            records = list(csv.DictReader(stream))
        assert len(records) == 2*design['repetitions']
        for group, truth, sign in [('comparator_null', q, 1), ('ellis_power', p, -1)]:
            fval, gval = signed_drift(.02, u, d, sign)
            fs, gs = signed_drift(.02, u, d+case['residual'], sign)
            subset = [r for r in records if r['group'] == group]
            covered = 0
            for rep, record in enumerate(subset):
                arrays = [np.array([int(record[f'{prefix}{i}']) for i in range(24)])
                          for prefix in ['y', 'rf', 'rg', 'vf', 'vg']]
                expected = [rng.binomial(n, (1-fs)*truth+gs*(1-truth)),
                            rng.binomial(m, .02, 24), rng.binomial(m, .02, 24),
                            rng.binomial(v, fval), rng.binomial(v, gval)]
                assert int(record['replicate']) == rep
                for actual, target in zip(arrays, expected):
                    assert np.array_equal(actual, target)
                covered += float(record['drift_low']) <= d <= float(record['drift_high'])
                experiments += 1
                if rep % 50 == 0:
                    y, rf, rg, vf, vg = arrays
                    result, bounds = evaluate(y, n, rf, rg, m, vf, vg, v, k, design['guard_eta'])
                    np.testing.assert_allclose(bounds, [float(record['drift_low']), float(record['drift_high'])], atol=1e-14)
                    for method in METHODS:
                        assert result[method]['status'] == record[method+'_status']
                        assert result[method]['nodes'] == int(record[method+'_nodes'])
                        audited += 1
                    regions, _ = validation_regions(y, n, rf, rg, m, vf, vg, v, design['guard_eta'])
                    direct = regions['direct_validation']
                    for larger in ['scalar_transfer', 'guarded_validation']:
                        wide = regions[larger]
                        assert np.all(direct[0] >= wide[0]-1e-10)
                        assert np.all(direct[1] <= wide[1]+1e-10)
                        nesting_checks += 1
            cr = next(r for r in coverage if int(r['case']) == index and r['group'] == group)
            assert covered == int(cr['covered'])
            np.testing.assert_allclose(wilson(covered, len(subset)), [float(cr['wilson_low']), float(cr['wilson_high'])], atol=1e-14)
            for key in ['lower', 'upper']:
                field = 'drift_low' if key == 'lower' else 'drift_high'
                assert np.median([float(r[field]) for r in subset]) == float(cr[key+'_median'])
            for method in METHODS:
                rate = next(r for r in rates if int(r['case']) == index and r['group'] == group and r['method'] == method)
                rejected = sum(r[method+'_status'] == 'excluded' for r in subset)
                assert rejected == int(rate['rejected'])
                assert float(rate['rate']) == rejected/len(subset)
                assert int(rate['nodes']) == sum(int(r[method+'_nodes']) for r in subset)
                assert int(rate['undetermined']) == sum(r[method+'_status'] == 'undetermined' for r in subset)
                np.testing.assert_allclose(wilson(rejected, len(subset)), [float(rate['wilson_low']), float(rate['wilson_high'])], atol=1e-14)
                valid = method != 'stale_reference' and (case['residual'] <= design['guard_eta'] if method == 'guarded_validation' else case['residual'] == 0)
                assert int(rate['assumed_transfer_valid']) == valid
                if valid and group == 'comparator_null':
                    valid_null_rejects += rejected; valid_null_records += len(subset)
    chosen = next(r for r in rates if int(r['B']) == 1600000 and float(r['fraction']) == .35 and r['panel'] == 'allocation')
    with gzip.open(OUT/chosen['file'], 'rt') as stream:
        values = [100*float(r['drift_high']) for r in csv.DictReader(stream) if r['group'] == 'ellis_power']
    hist = read(OUT/'bound-histogram.csv')
    edges = [float(r['left_pp']) for r in hist]+[float(hist[-1]['right_pp'])]
    counts, _ = np.histogram(values, edges)
    assert counts.sum() == design['repetitions']
    for row, count in zip(hist, counts):
        assert count == int(row['count'])
        np.testing.assert_allclose(wilson(count, len(values)), [float(row['wilson_low']), float(row['wilson_high'])], atol=1e-14)
    summary = json.loads((OUT/'summary.json').read_text())
    assert experiments == summary['experiments']
    result = dict(experiments_replayed=experiments, count_slots_replayed=experiments*120,
                  solver_decisions_replayed=audited, numerical_region_nesting_checks=nesting_checks,
                  rate_rows_checked=len(rates), drift_coverage_rows_checked=len(coverage),
                  histogram_bins_checked=len(hist), valid_null_rejections=valid_null_rejects,
                  valid_null_paired_decision_records=valid_null_records, hashes='all verified',
                  note='Replay checks implementation consistency; analytical coverage is conditional on response transfer.')
    (OUT/'verification.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
