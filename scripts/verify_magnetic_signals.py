"""Replay Study 014 numerical and analytic checks. Author: Biswajit Jana."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from magnetic_signals import gaussian_average,gaussian_fourier,mouth_probability,eventual_fraction,collection_pulses

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results/study-014'


def read(name):
    with (OUT/name).open(encoding='utf-8') as stream:
        return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(stream)]


def main():
    summary=json.loads((OUT/'summary.json').read_text(encoding='utf-8'))
    for path,h in summary['source_hashes'].items():
        assert hashlib.sha256((ROOT/path).read_bytes().replace(b'\r\n',b'\n')).hexdigest()==h,path
    hashes=json.loads((OUT/'sha256.json').read_text(encoding='utf-8'))
    for path,h in hashes.items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    maximum=0.; replays=0
    for row in read('bandwidth.csv')[::10]:
        value,error=gaussian_average(row['x'],row['bandwidth'],tolerance=2e-12)
        maximum=max(maximum,abs(value-row['average']));replays+=1
        assert abs(value-row['average'])<1e-9
    for row in read('bandwidth_audit.csv'):
        independent=gaussian_fourier(row['x'],row['bandwidth'],row['detuning'])
        assert abs(independent-row['quadrature'])<1e-10
    for row in read('matching_audit.csv'):
        assert row['matching_error']<1e-11 and row['flux_error']<1e-12
    for row in read('time_audit.csv'):assert row['absolute_difference']<1e-10
    for row in read('compact_pulse_audit.csv'):
        assert abs(row['envelope_norm']-1)<1e-12
        assert row['absolute_difference']<1e-10
        if row['duration_over_ell']<2*np.pi:
            assert abs(row['direct_energy']-row['separated_prediction'])<1e-12
    for row in read('collection.csv'):
        p=mouth_probability(row['x']);R=1-p
        # Alternative expression, cancellation safe for the sampled x range.
        expected=p*p*(1-R**(2*row['pulses']))/(1-R*R)
        assert abs(expected-row['incident_energy_fraction'])<1e-12
    for row in read('scenario_audit.csv'):
        x=row['x'];assert abs(x/row['binding_energy_times_re']-.1)<1e-15
        assert abs(row['n90']/collection_pulses(x)-1)<1e-13
        assert row['sigma_t_at_hwhm_over_d']>1e9
        assert row['sigma_t_over_d']<.02
        assert row['bandwidth_over_carrier']<.001
    record=dict(author='Biswajit Jana',date='2026-11-25',artifact_hashes_checked=len(hashes),
        tighter_band_replays=replays,maximum_band_replay_difference=maximum,
        fourier_checks=len(read('bandwidth_audit.csv')),matching_checks=len(read('matching_audit.csv')),
        time_integral_checks=len(read('time_audit.csv')),scenario_checks=len(read('scenario_audit.csv')),
        compact_pulse_checks=len(read('compact_pulse_audit.csv')),
        conclusion='Finite-signal approximation checks passed; no observation or full validity certificate')
    (OUT/'verification.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
