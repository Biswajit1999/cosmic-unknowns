"""Artifact, ODE replay and fit checks. Author: Biswajit Jana."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from nuisance_identification import NullSpectrum, BOUNDS, STAGES
from conformal_scattering import scatter
from massive_conformal import scatter_massive

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT/'results/study-011'


def read(name):
    with (OUT/(name+'.csv')).open() as stream: return list(csv.DictReader(stream))


def main():
    design = json.loads((OUT/'design.json').read_text())
    for name, digest in json.loads((OUT/'sha256.json').read_text()).items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    for name, digest in design['source_hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes().replace(b'\r\n', b'\n')).hexdigest() == digest, name
    for name, digest in design['parent_hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    templates = read('templates'); model = NullSpectrum(templates)
    checks = read('heldout'); fits = read('fits'); predictions = read('predictions')
    for r in checks:
        v = float(model.functions[int(r['ell'])](float(r['q'])))
        assert abs(v-float(r['interpolated_T'])) < 1e-15
        assert abs(v-float(r['T'])-float(r['error'])) < 1e-15
    with (ROOT/'results/study-010/spectra.csv').open() as stream: parent = list(csv.DictReader(stream))
    truth = {(float(r['beta']), float(r['k']), int(r['ell'])):float(r['T']) for r in parent if float(r['mass']) == .3}
    for r in predictions:
        pars = [float(r[k]) for k in ['scale', 'mass', 'false_positive', 'false_negative']]
        k, ell = float(r['k']), int(r['ell'])
        v = model.predict([k], [ell], pars)[0]
        assert abs(v-float(r['interpolated_probability'])) < 1e-14
        assert float(r['truth']) == truth[(float(r['beta_truth']), k, ell)]
        direct = pars[2]+(1-pars[2]-pars[3])*float(r['direct_T'])
        assert abs(direct-float(r['direct_probability'])) < 1e-15
        assert abs(direct-float(r['truth'])-float(r['direct_residual'])) < 1e-15
        assert abs(direct-v-float(r['interpolation_difference'])) < 1e-14
    replayed = 0
    for r in predictions[::16]:
        q = np.sqrt((float(r['scale'])*float(r['k']))**2-float(r['mass'])**2)
        v = scatter(q, ell=int(r['ell']), form='reduced')['T']
        assert abs(v-float(r['direct_T'])) < 1e-13
        replayed += 1
    for row in fits:
        rr = [r for r in predictions if all(r[k] == row[k] for k in ['beta_truth', 'observables', 'stage'])]
        residual = np.array([float(r['direct_residual']) for r in rr])
        assert abs(float(row['direct_rms'])-np.sqrt(np.mean(residual**2))) < 1e-15
        assert abs(float(row['direct_maximum_residual'])-max(abs(residual))) < 1e-15
        for name, (low, high) in BOUNDS.items(): assert low <= float(row[name]) <= high
    for beta in ['-0.2', '0.2']:
        for obs in ['ell0', 'joint']:
            rr = [next(r for r in fits if r['beta_truth'] == beta and r['observables'] == obs and r['stage'] == s) for s in STAGES]
            assert np.all(np.diff([float(r['interpolated_rms']) for r in rr]) <= 1e-10)
    # Refit the most flexible two-mode cases from all saved deterministic starts.
    refits = 0
    for beta in ['-0.2', '0.2']:
        rr = [r for r in predictions if r['beta_truth'] == beta and r['observables'] == 'joint' and r['stage'] == 'mass_scale_response']
        k = np.array([float(r['k']) for r in rr]); ell = np.array([int(r['ell']) for r in rr]); target = np.array([float(r['truth']) for r in rr])
        parameters, residual, starts = model.fit(k, ell, target, 'mass_scale_response')
        saved = next(r for r in fits if r['beta_truth'] == beta and r['observables'] == 'joint' and r['stage'] == 'mass_scale_response')
        assert abs(np.sqrt(np.mean(residual**2))-float(saved['interpolated_rms'])) < 1e-12
        refits += 1
    parity = read('parity'); derivatives = read('derivatives')
    for r in parity[::6]:
        v = scatter_massive(float(r['k']), float(r['beta']), float(r['mass']), int(r['ell']))['T']
        assert abs(v-float(r['T'])) < 1e-13
    for r in derivatives:
        h = float(r['step'])
        db = (float(r['T_beta_plus'])-float(r['T_beta_minus']))/(2*h)
        dm = (float(r['T_mass_plus'])-float(r['T_mass_minus']))/(2*h)
        assert abs(db-float(r['dT_dbeta'])) < 1e-15
        assert abs(dm-float(r['dT_dmass'])) < 1e-15
        assert abs(db+np.pi*.3*dm/2-float(r['identity_residual'])) < 1e-15
    summary = json.loads((OUT/'summary.json').read_text())
    assert len(templates) == summary['template_integrations'] == 2002
    assert len(predictions) == summary['direct_fit_integrations'] == 504
    assert len(checks) == summary['heldout_integrations'] == 242
    assert len(fits) == summary['fitted_cases'] == 16
    assert len(parity) == summary['parity_integrations'] == 66
    assert 4*len(derivatives) == summary['derivative_integrations'] == 64
    assert max(abs(float(r['error'])) for r in checks) == summary['maximum_heldout_interpolation_error']
    assert max(abs(float(r['interpolation_difference'])) for r in predictions) == summary['maximum_fit_interpolation_error']
    result = dict(hashes='all verified', prediction_rows_checked=len(predictions), direct_fit_replays=replayed,
        heldout_rows_checked=len(checks), flexible_joint_refits=refits, parity_replays=len(parity[::6]),
        derivative_rows_checked=len(derivatives), nested_fit_rms='nonincreasing in saved fits',
        note='Refits reproduce local optimization witnesses; not certified global minima or confidence intervals.')
    (OUT/'verification.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(result))


if __name__ == '__main__': main()
