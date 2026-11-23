"""Replay stored Study 013 results and geometric constraints. Author: Biswajit Jana."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from magnetic_energy import MagneticEnergy, scale_diagnostics

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT/'results/study-013'


def read(name):
    with (OUT/name).open(encoding='utf-8') as stream:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(stream)]


def main():
    summary = json.loads((OUT/'summary.json').read_text(encoding='utf-8'))
    for path, expected in summary['source_hashes'].items():
        assert hashlib.sha256((ROOT/path).read_bytes().replace(b'\r\n', b'\n')).hexdigest() == expected, path
    hashes = json.loads((OUT/'sha256.json').read_text(encoding='utf-8'))
    for path, expected in hashes.items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == expected, path
    model = MagneticEnergy(1024, 1024)
    rows = read('equilibrium.csv'); max_replay = 0.; opt_difference = []
    for i, row in enumerate(rows):
        replay = model.minimum(row['s'])
        difference = abs(replay['y']/row['y']-1); max_replay = max(max_replay, difference)
        assert difference < 1e-9
        assert abs(replay['energy']-row['energy']) < 1e-10
        assert row['y'] > max(1., 1/row['s']) and row['energy'] < 0 and row['curvature'] > 0
        assert abs(row['stationarity']) < 1e-10
        if i % 20 == 0:
            opt = minimize_scalar(lambda logy: model.energy(np.exp(logy), row['s']),
                bounds=(np.log(.5), np.log(max(100., 10/row['s']))), method='bounded', options={'xatol': 1e-12})
            difference = abs(np.exp(opt.x)/replay['y']-1)
            assert opt.success and difference < 2e-7
            opt_difference.append(difference)
    coords = read('field_lines.csv'); flux_residual = 0.
    for row in coords:
        z = row['z_over_d']; r = row['rho_over_d']
        if abs(row['eta']) < .999:
            flux = (z+.5)/np.hypot(z+.5,r)-(z-.5)/np.hypot(z-.5,r)
            flux_residual = max(flux_residual, abs(flux-row['nu']))
    assert flux_residual < 1e-11
    cut = read('tail_omission.csv')
    for row in cut:
        assert row['y_full'] <= row['y_cut'] <= row['bound_upper']
    scales = read('scale_diagnostics.csv')
    for row in scales:
        baseline = min(rows, key=lambda v: abs(v['s']-row['s']))
        pred = scale_diagnostics(row['q'], row['g'], row['s'], baseline['y'], baseline['energy'])
        for k, v in pred.items(): assert abs(v/row[k]-1) < 1e-12
    audit = read('field_length_audit.csv')
    assert max(r['relative_difference'] for r in audit) < 1e-9
    record = dict(author='Biswajit Jana', date='2026-11-23', artifact_hashes_checked=len(hashes),
        equilibrium_replays=len(rows), replay_order=[1024,1024], maximum_relative_y_replay=max_replay,
        independent_energy_minimizations=len(opt_difference), maximum_minimizer_difference=max(opt_difference),
        field_coordinate_rows=len(coords), maximum_field_label_residual=flux_residual,
        omission_bound_checks=len(cut), scale_rows_replayed=len(scales),
        independent_angle_integrals=len(audit), conclusion='Reduced-model checks passed; no full stability or observation claim')
    (OUT/'verification.json').write_text(json.dumps(record, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__': main()
