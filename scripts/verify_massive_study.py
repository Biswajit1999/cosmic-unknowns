"""Independent checks of Study 010 artifacts. Author: Biswajit Jana."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from mpmath import mp
from massive_conformal import scatter_massive, tail_data, scalar_source_obstruction
from conformal_scattering import scatter, stress
from verify_conformal_study import tensor_check

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT/'results/study-010'


def read(name):
    with (OUT/(name+'.csv')).open() as stream:
        return list(csv.DictReader(stream))


def amplitude(row, name):
    return complex(float(row[name+'_real']), float(row[name+'_imag']))


def main():
    for path, digest in json.loads((OUT/'sha256.json').read_text()).items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest, path
    for path, digest in json.loads((OUT/'design.json').read_text())['source_hashes'].items():
        assert hashlib.sha256((ROOT/path).read_bytes().replace(b'\r\n', b'\n')).hexdigest() == digest, path
    tensor_residuals = []; series_checks = 0
    with mp.workdps(70):
        for beta in ['-0.2', '0', '0.2']:
            b = mp.mpf(beta)
            for location in ['-4', '-0.3', '0', '2']:
                x = mp.mpf(location); f = 1+b*(mp.atan(x)-mp.pi/2)
                values = tensor_check(x, mp.mpf('1.1'), b)
                error = abs(values[1]+values[3]+2*b*x/(f**3*(1+x*x)**2))
                assert error < mp.mpf('1e-60')
                tensor_residuals.append(str(error))
            for side in [-1, 1]:
                for ell in [0, 1]:
                    f0 = mp.mpf(1) if side == 1 else 1-mp.pi*b
                    exact = lambda z: ell*(ell+1)*z*z/(1+z*z)+z**4/(1+z*z)**2+mp.mpf('.09')*(f0-side*b*mp.atan(z))**2
                    coefficients = mp.taylor(exact, 0, 14)
                    saved = tail_data(.8, float(b), .3, ell, side, order=12)[2]
                    np.testing.assert_allclose(saved, [float(v) for v in coefficients], atol=3e-14, rtol=3e-14)
                    series_checks += 1
    spectra = read('spectra'); audits = read('convergence')
    for row in spectra:
        t, tf, r = [amplitude(row, name) for name in ['t', 't_flux', 'r']]
        k, b, m = [float(row[name]) for name in ['k', 'beta', 'mass']]
        qr = np.sqrt(k*k-m*m); ql = np.sqrt(k*k-m*m*(1-np.pi*b)**2)
        assert abs(qr-float(row['q_right'])) < 2e-15
        assert abs(ql-float(row['q_left'])) < 2e-15
        assert abs(tf-np.sqrt(ql/qr)*t) < 2e-15
        assert abs(abs(tf)**2-float(row['T'])) < 2e-15
        assert abs(abs(r)**2-float(row['R'])) < 2e-15
        assert abs(abs(abs(tf)**2+abs(r)**2-1)-float(row['flux_error'])) < 3e-15
    replay = 0
    for row in spectra[::27]:
        fresh = scatter_massive(float(row['k']), float(row['beta']), float(row['mass']), int(row['ell']))
        for name in ['t_flux', 'r']:
            assert abs(fresh[name]-amplitude(row, name)) < 1e-13
        replay += 1
    baselines = {}
    for row in audits:
        key = (float(row['k']), float(row['beta']), int(row['ell']))
        if key not in baselines:
            baselines[key] = scatter_massive(key[0], key[1], .3, key[2])
        base = baselines[key]
        assert abs(abs(amplitude(row, 't_flux')-base['t_flux'])-float(row['delta_t_flux'])) < 1e-13
        if row['delta_r']:
            assert abs(abs(amplitude(row, 'r')-base['r'])-float(row['delta_r'])) < 1e-13
        assert abs(abs(float(row['T'])-base['T'])-float(row['delta_transmission_probability'])) < 1e-13
    limits = []
    for ell in [0, 1]:
        for k in [.54, .85, 1.6]:
            for mass in [0., .3]:
                a = scatter_massive(k, 0., mass, ell)
                b = scatter(np.sqrt(k*k-mass*mass), ell=ell, form='reduced')
                limits.append(max(abs(a['r']-b['r']), abs(a['t']-b['t'])))
    assert max(limits) < 3e-8
    profiles = read('source_profiles')
    for row in profiles:
        x, b = float(row['X']), float(row['beta'])
        rho, pr, pt = stress(x, b)
        assert abs(rho+pt-float(row['rho_plus_pt'])) < 1e-14
        assert abs(scalar_source_obstruction(x, b)-float(row['analytic_obstruction'])) < 1e-14
    for row in read('embedding_meridian'):
        x = float(row['X'])
        assert abs(float(row['radius'])-np.sqrt(1+x*x)) < 1e-15
        assert abs(float(row['embedding_z'])-np.arcsinh(x)) < 1e-15
        assert abs((x/np.sqrt(1+x*x))**2+(1/np.sqrt(1+x*x))**2-1) < 5e-16
    for row in read('thresholds'):
        assert abs(float(row['left_threshold'])-.3*(1-np.pi*float(row['beta']))) < 1e-15
    summary = json.loads((OUT/'summary.json').read_text())
    assert len(spectra) == summary['spectral_integrations'] == 972
    assert len(audits) == summary['audit_comparisons'] == 108
    for key, col in [('maximum_flux_error', 'flux_error'), ('maximum_tail_current_error', 'tail_current_error')]:
        assert max(float(r[col]) for r in spectra) == summary[key]
    assert max(max(float(r['delta_t_flux']), float(r['delta_r'] or 0)) for r in audits) == summary['maximum_audit_amplitude_difference']
    assert max(float(r['delta_t_flux']) for r in audits if r['check'] == 'reverse_incidence') == summary['maximum_reciprocity_difference']
    keyed = {(float(r['k']), float(r['beta']), float(r['mass']), int(r['ell'])): r for r in spectra}
    gaps = []; massless_gaps = []
    for key, row in keyed.items():
        if key[1] == 0: continue
        delta = float(row['T'])-float(keyed[(key[0], 0., key[2], key[3])]['T'])
        if key[2] == .3: gaps.append(dict(ell=key[3], beta=key[1], k=key[0], delta_T=delta))
        else: massless_gaps.append(abs(delta))
    assert max(gaps, key=lambda r: abs(r['delta_T'])) == summary['maximum_transmission_probability_gap']
    assert max(abs(float(r['rho_plus_pt'])-float(r['analytic_obstruction'])) for r in profiles) == summary['source_identity_error']
    result = dict(hashes='all verified', spectral_rows_checked=len(spectra), spectral_replays=replay,
        audit_rows_checked=len(audits), audit_baselines_replayed=len(baselines),
        independent_tensor_checks=len(tensor_residuals), tensor_identity_residuals=tensor_residuals,
        high_precision_digits=70, independently_differentiated_tail_series=series_checks,
        ultrastatic_limit_checks=len(limits), maximum_ultrastatic_amplitude_difference=max(limits),
        maximum_massless_probability_difference=max(massless_gaps),
        source_rows_checked=len(profiles), embedding_points_checked=161, threshold_points_checked=121,
        note='Deterministic point checks and convergence tests; no interval-certified global error bound.')
    (OUT/'verification.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'tensor_identity_residuals'}))


if __name__ == '__main__':
    main()
