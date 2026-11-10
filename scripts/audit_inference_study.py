"""Independent optimizer/bound sensitivity audit. Author: Biswajit Jana."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from transmission_inference import pt_logit, pt_transmission, deviance

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/"results/study-002"


def main():
    with (OUT/"trials.csv").open() as f:
        trials = list(csv.DictReader(f))
    with (OUT/"spectra.csv").open() as f:
        spectra = list(csv.DictReader(f))
    configs = json.loads((OUT/"design.json").read_text())
    k_by_design = {d: np.array([float(r["K"]) for r in spectra if r["design"] == d])
                   for d in configs["designs"]}
    audit, values = [], {}
    for row in trials:
        d, n, group, rep = row["design"], int(row["N"]), row["group"], int(row["repetition"])
        original = float(row["deviance"])
        revised = original
        # Every flagged/boundary fit, plus a fixed sample of ordinary fits.
        if row["optimizer_success"] == "0" or row["boundary"] == "1" or rep < 5:
            k = k_by_design[d]
            y = np.array([int(row[f"count_{i:02d}"]) for i in range(24)])
            theta = np.array([float(row["logC"]), float(row["logb"])])
            def loss(t):
                eta = pt_logit(k, t)
                return float(np.sum(np.logaddexp(0, eta)-y/n*eta))
            fit = minimize(loss, theta, method="Nelder-Mead",
                bounds=[(-25., 40.), (np.log(.01), np.log(60.))],
                options={"xatol": 1e-10, "fatol": 1e-13, "maxiter": 4000})
            candidate = deviance(y, n, pt_transmission(k, fit.x))
            revised = min(original, candidate)
            audit.append(dict(design=d, N=n, group=group, repetition=rep,
                original_deviance=original, alternative_deviance=candidate,
                deviance_improvement=original-revised, alternative_success=int(fit.success),
                original_non_success=int(row["optimizer_success"] == "0"),
                original_boundary=int(row["boundary"] == "1")))
        values.setdefault((d, n, group), [[], []])
        values[(d, n, group)][0].append(original)
        values[(d, n, group)][1].append(revised)
    comparisons = []
    for d in configs["designs"]:
        for n in configs["trials_per_channel"]:
            cal, altcal = values[(d, n, "calibration")]
            t = float(np.quantile(cal, .95, method="higher"))
            at = float(np.quantile(altcal, .95, method="higher"))
            for g in ("null_check", "ellis_power"):
                x, ax = map(np.array, values[(d, n, g)])
                comparisons.append(dict(design=d, N=n, group=g,
                    original_rejected=int(sum(x > t)), alternative_rejected=int(sum(ax > at)),
                    changed_decisions=int(sum((x > t) != (ax > at))),
                    original_threshold=t, alternative_threshold=at))
    with (OUT/"independent-optimizer-audit.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(audit[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(audit)
    summary = {"author": "Biswajit Jana", "algorithm": "bounded Nelder-Mead, wider bounds",
        "logC_bounds": [-25, 40], "b_bounds": [.01, 60], "audited_fits": len(audit),
        "original_non_success_audited": sum(r["original_non_success"] for r in audit),
        "alternative_non_success": sum(1-r["alternative_success"] for r in audit),
        "max_deviance_improvement": max(r["deviance_improvement"] for r in audit),
        "max_non_success_deviance_improvement": max(r["deviance_improvement"] for r in audit if r["original_non_success"]),
        "changed_decisions": sum(r["changed_decisions"] for r in comparisons),
        "comparisons": comparisons,
        "limitation": "Selected audit; not a global optimizer proof or uniform test calibration."}
    (OUT/"independent-audit.json").write_text(json.dumps(summary, indent=2)+"\n", newline="\n")
    paths = [p for folder in (OUT, ROOT/"figures/study-002") for p in sorted(folder.iterdir())
             if p.suffix in (".json", ".csv", ".png", ".svg") and p.name != "sha256.json"]
    hashes = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (OUT/"sha256.json").write_text(json.dumps(hashes, indent=2)+"\n", newline="\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "comparisons"}, indent=2))


if __name__ == "__main__":
    main()
