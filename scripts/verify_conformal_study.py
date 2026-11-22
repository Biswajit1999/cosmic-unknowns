"""Independent tensor and artifact checks. Author: Biswajit Jana."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from mpmath import mp
from conformal_scattering import geometry, stress, scatter, common_potential

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT/'results/study-009'


def connection(x, theta, beta):
    """Christoffel symbols directly from the four-dimensional diagonal metric."""
    f = 1+beta*(mp.atan(x)-mp.pi/2); fp = beta/(1+x*x)
    r = f*mp.sqrt(1+x*x); rp = fp*mp.sqrt(1+x*x)+f*x/mp.sqrt(1+x*x)
    sine, cosine = mp.sin(theta), mp.cos(theta)
    diagonal = [-f*f, f*f, r*r, r*r*sine*sine]
    derivatives = [[mp.mpf(0)]*4 for _ in range(4)]
    derivatives[1] = [-2*f*fp, 2*f*fp, 2*r*rp, 2*r*rp*sine*sine]
    derivatives[2][3] = 2*r*r*sine*cosine
    gamma = [[[((derivatives[b][a] if c == a else 0)+(derivatives[c][a] if b == a else 0)
                 -(derivatives[a][b] if b == c else 0))/(2*diagonal[a])
                for c in range(4)] for b in range(4)] for a in range(4)]
    return diagonal, gamma


def tensor_check(x, theta, beta):
    diagonal, gamma = connection(x, theta, beta)
    def derivative(axis, up, first, second):
        if axis == 1:
            return mp.diff(lambda value: connection(value, theta, beta)[1][up][first][second], x)
        if axis == 2:
            return mp.diff(lambda value: connection(x, value, beta)[1][up][first][second], theta)
        return mp.mpf(0)
    ricci = []
    for a in range(4):
        value = sum(derivative(c, c, a, a)-derivative(a, c, a, c) for c in range(4))
        value += sum(gamma[c][c][d]*gamma[d][a][a]-gamma[c][a][d]*gamma[d][a][c]
                     for c in range(4) for d in range(4))
        ricci.append(value)
    scalar = sum(ricci[a]/diagonal[a] for a in range(4))
    mixed = [ricci[a]/diagonal[a]-scalar/2 for a in range(4)]
    return scalar, -mixed[0], mixed[1], mixed[2], mixed[3]


def read(name):
    with (OUT/name).open() as stream:
        return list(csv.DictReader(stream))


def main():
    design = json.loads((OUT/'design.json').read_text())
    for name, digest in json.loads((OUT/'sha256.json').read_text()).items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    for name, digest in design['source_hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes().replace(b'\r\n', b'\n')).hexdigest() == digest, name
    tensor_rows = []; max_identity = mp.mpf(0)
    with mp.workdps(70):
        for beta in ['-0.2', '0', '0.2']:
            b = mp.mpf(beta)
            for location in ['-4', '-0.3', '0', '2']:
                x = mp.mpf(location)
                values = tensor_check(x, mp.mpf('1.1'), b)
                f = 1+b*(mp.atan(x)-mp.pi/2)
                expected = -2/(f*f*(1+x*x)**2)
                assert abs(values[0]-expected) < mp.mpf('1e-60')
                assert abs(values[3]-values[4]) < mp.mpf('1e-60')
                np.testing.assert_allclose([float(v) for v in values[1:4]], stress(float(x), float(b)), atol=2e-12, rtol=2e-12)
                tensor_rows.append(dict(beta=float(b), X=float(x), ricci=float(values[0]),
                    rho=float(values[1]), pr=float(values[2]), pt=float(values[3]),
                    curvature_identity_residual=str(abs(values[0]-expected))))
                radius = lambda z: (1+b*(mp.atan(z)-mp.pi/2))*mp.sqrt(1+z*z)
                r = radius(x); rpp = mp.diff(radius, x, 2)
                for ell in [0, 1, 2, 5]:
                    for xi in [mp.mpf(0), mp.mpf(1)/6, mp.mpf('.4')]:
                        exact = ell*(ell+1)*f*f/r**2+rpp/r+xi*f*f*values[0]
                        reference = ell*(ell+1)/(1+x*x)+(1-2*xi)/(1+x*x)**2
                        residual = abs(exact-reference)
                        max_identity = max(max_identity, residual)
                        assert residual < mp.mpf('1e-60')
    spectra = read('spectra.csv'); replayed = 0
    for row in spectra[::18]:
        fresh = scatter(float(row['k']), int(row['ell']), float(row['beta']))
        for name in ['r', 't']:
            saved = complex(float(row[name+'_real']), float(row[name+'_imag']))
            assert abs(fresh[name]-saved) < 1e-13
        replayed += 1
    for row in spectra:
        t = complex(float(row['t_real']), float(row['t_imag'])); r = complex(float(row['r_real']), float(row['r_imag']))
        assert abs(abs(t)**2-float(row['T'])) < 2e-15
        assert abs(abs(r)**2-float(row['R'])) < 2e-15
        assert abs(abs(abs(t)**2+abs(r)**2-1)-float(row['flux_error'])) < 2e-15
    audits = read('convergence.csv'); baselines = {}
    for row in audits:
        key = (float(row['k']), int(row['ell']), float(row['beta']))
        if key not in baselines:
            baselines[key] = scatter(*key)
        base = baselines[key]
        t = complex(float(row['t_real']), float(row['t_imag']))
        r = complex(float(row['r_real']), float(row['r_imag']))
        assert abs(abs(t-base['t'])-float(row['delta_t'])) < 1e-13
        assert abs(abs(r-base['r'])-float(row['delta_r'])) < 1e-13
        assert abs(abs(abs(t)**2-base['T'])-float(row['delta_transmission_probability'])) < 1e-13
        assert abs(abs(abs(t)**2+abs(r)**2-1)-float(row['flux_error'])) < 3e-15
    summary = json.loads((OUT/'summary.json').read_text())
    assert len(spectra) == summary['spectral_integrations'] == 729
    assert len(audits) == summary['convergence_comparisons'] == 165
    assert abs(max(float(r['flux_error']) for r in spectra)-summary['max_flux_error']) < 1e-20
    for plane, field in [(False, 'maximum_nonplane_convergence_difference'),
                         (True, 'maximum_plane_boundary_difference')]:
        value = max(max(float(r['delta_r']), float(r['delta_t'])) for r in audits
                    if (r['check'] == 'plane_boundary') == plane)
        assert abs(value-summary[field]) < 1e-20
    keyed = {(float(r['k']), int(r['ell']), float(r['beta'])): r for r in spectra}
    differences = []
    for key, row in keyed.items():
        if key[2] == 0:
            continue
        base = keyed[(key[0], key[1], 0.)]
        for name in ['r', 't']:
            differences.append(abs(complex(float(row[name+'_real']), float(row[name+'_imag']))
                                   -complex(float(base[name+'_real']), float(base[name+'_imag']))))
    assert abs(max(differences)-summary['maximum_geometry_amplitude_difference']) < 1e-20
    profiles = read('profiles.csv')
    for row in profiles:
        x, beta = float(row['X']), float(row['beta'])
        f, _, _, r, _, _ = geometry(x, beta)
        assert abs(float(row['radius'])-r) < 1e-14
        for ell in range(3):
            assert abs(float(row[f'U{ell}'])-common_potential(x, ell)) < 1e-13
        assert abs(float(row['massive_delta_U'])-.09*(f*f-1)) < 1e-14
    throats = read('throats.csv')
    for row in throats:
        b, x = float(row['beta']), float(row['X'])
        f, fp, fpp, r, rp, rpp = geometry(x, b)
        assert abs(rp) < 1e-12 and rpp > 0 and float(row['radial_nec']) < 0
        assert abs(float(row['radial_nec'])+2/(f*f*(1+x*x)**2)) < 1e-12
    result = dict(spectral_rows_checked=len(spectra), original_field_integrations_replayed=replayed,
        convergence_rows_checked=len(audits), audit_baseline_integrations_replayed=len(baselines),
        geometry_profile_rows_checked=len(profiles), throat_rows_checked=len(throats),
        independent_four_dimensional_tensor_checks=tensor_rows, high_precision_digits=70,
        maximum_high_precision_potential_identity_residual=str(max_identity), hashes='all verified',
        note='Point checks support the displayed analytic identities; not interval-certified full-domain ODE error bounds.')
    (OUT/'verification.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'independent_four_dimensional_tensor_checks'}))


if __name__ == '__main__':
    main()
