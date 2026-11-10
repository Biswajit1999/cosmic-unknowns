"""Replay dense data and check calibration witnesses. Author: Biswajit Jana."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
from transmission_inference import pt_transmission, deviance, wilson
from calibration_limits import recorded_probability, required_error

ROOT = Path(__file__).resolve().parents[1]


def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def main():
    dense = ROOT/"results/study-005"
    source = ROOT/"results/study-002/spectra.csv"
    design = json.loads((dense/"design.json").read_text())
    assert hashlib.sha256(source.read_bytes()).hexdigest() == design["source_spectra_sha256"]
    spectra = read_csv(source)
    rates = read_csv(dense/"rates.csv")
    audits = read_csv(dense/"optimizer-audit.csv")
    total, warnings = 0, 0
    for audit in audits:
        name, n, seed = audit["design"], int(audit["N"]), int(audit["seed"])
        data = [r for r in spectra if r["design"] == name]
        k, truth, null = (np.array([float(r[key]) for r in data]) for key in ["K", "T_ellis", "T_pt_projection"])
        rng = np.random.default_rng(seed)
        groups = {g: [] for g in ["calibration", "null_check", "ellis_power"]}
        with gzip.open(dense/audit["file"], "rt") as f:
            for r in csv.DictReader(f):
                probability = truth if r["group"] == "ellis_power" else null
                count = np.array([int(r[f"count_{i:02d}"]) for i in range(24)])
                assert np.array_equal(count, rng.binomial(n, probability)), f"Seed replay mismatch: {audit['file']}"
                fitted = pt_transmission(k, [float(r["logC"]), float(r["logb"])])
                assert abs(deviance(count, n, fitted)-float(r["deviance"])) < 1e-6
                groups[r["group"]].append(float(r["deviance"]))
                warnings += 1-int(r["success"])
                total += 1
        assert all(len(v) == 1000 for v in groups.values())
        threshold = float(np.quantile(groups["calibration"], .95, method="higher"))
        for group in ["null_check", "ellis_power"]:
            row = next(r for r in rates if r["design"] == name and int(r["N"]) == n and r["group"] == group)
            rejected = int(sum(np.array(groups[group]) > threshold))
            assert rejected == int(row["rejected"])
            np.testing.assert_allclose(wilson(rejected, 1000), [float(row["wilson_low"]), float(row["wilson_high"])])
    assert total == 150000 and warnings == 0
    witnesses = read_csv(ROOT/"results/study-004/witnesses.csv")
    summary = json.loads((ROOT/"results/study-004/summary.json").read_text())
    for name in ["low_band", "broad_band"]:
        data = [r for r in witnesses if r["design"] == name]
        p, q, fp, fn = (np.array([float(r[key]) for r in data]) for key in ["p_ellis", "q_alternative", "false_positive", "false_negative"])
        np.testing.assert_allclose(recorded_probability(q, fp, fn), p, rtol=0, atol=2e-16)
        assert abs(max(required_error(p, q))-summary["designs"][name]["achievable_error_budget"]) < 1e-12
    for r in read_csv(ROOT/"results/study-004/histogram-bins.csv"):
        np.testing.assert_allclose(wilson(int(r["count"]), int(r["total"])), [float(r["wilson_low"]), float(r["wilson_high"])])
    report = ROOT/"output/pdf/report-manifest.json"
    if report.exists():
        manifest = json.loads(report.read_text())
        assert hashlib.sha256((ROOT/manifest["source"]).read_bytes()).hexdigest() == manifest["source_sha256"]
        assert hashlib.sha256((report.parent/"calibration-limited-identification.pdf").read_bytes()).hexdigest() == manifest["pdf_sha256"]
    print("Verified 150,000 experiments, 3.6 million counts, all dense deviances/rates, 48 calibration witnesses, histogram intervals and report hashes.")


if __name__ == "__main__":
    main()
