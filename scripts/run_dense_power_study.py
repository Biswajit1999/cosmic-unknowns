"""Dense, independently sampled power curve. Author: Biswajit Jana."""
import csv
import gzip
import hashlib
import io
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from transmission_inference import fit_pt, pt_transmission, pt_logit, deviance, wilson, BOUNDS

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/"results/study-005"
FIG = ROOT/"figures/study-005"
COUNTS = np.unique(np.rint(np.geomspace(100, 100000, 25)).astype(int)).tolist()
REPS = 1000


def case(args):
    design, n, seed, k, truth, theta, null = args
    rng = np.random.default_rng(seed)
    bootstrap = np.random.default_rng(seed+1000000)
    filename = f"{design}-N{n:06d}.csv.gz"
    groups = {}
    warnings, bounds, fallbacks, improvements = 0, 0, 0, []
    with (OUT/filename).open("wb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as stream:
                writer = csv.writer(stream, lineterminator="\n")
                writer.writerow(["group", "repetition", "deviance", "logC", "logb", "success", "boundary", *[f"count_{i:02d}" for i in range(24)]])
                for group, probability in [("calibration", null), ("null_check", null), ("ellis_power", truth)]:
                    values = []
                    for rep in range(REPS):
                        y = rng.binomial(n, probability)
                        fitted, ok, bound, gradient = fit_pt(k, y/n, start=theta)
                        initial = deviance(y, n, pt_transmission(k, fitted))
                        if not ok or (gradient > 1e-5 and not bound) or rep < 2:
                            t, success, active, grad = fit_pt(k, y/n, start=fitted, multistart=True)
                            if deviance(y, n, pt_transmission(k, t)) <= initial:
                                fitted, ok, bound = t, success, active
                            fallbacks += 1
                        if not ok:
                            def loss(t):
                                eta = pt_logit(k, t)
                                return float(np.sum(np.logaddexp(0, eta)-y/n*eta))
                            alt = minimize(loss, fitted, method="Nelder-Mead", bounds=BOUNDS,
                                options={"xatol": 1e-9, "fatol": 1e-13, "maxiter": 3000})
                            # A successful independent fit within a small likelihood
                            # tolerance resolves a line-search termination warning.
                            if alt.success and loss(alt.x) <= loss(fitted)+1e-12:
                                fitted, ok = alt.x, True
                                bound = any(abs(fitted[i]-edge) < 1e-5 for i, bb in enumerate(BOUNDS) for edge in bb)
                        final = deviance(y, n, pt_transmission(k, fitted))
                        improvements.append(initial-final)
                        warnings += int(not ok)
                        bounds += int(bound)
                        values.append(final)
                        writer.writerow([group, rep, final, *fitted, int(ok), int(bound), *y])
                    groups[group] = np.array(values)
    threshold = float(np.quantile(groups["calibration"], .95, method="higher"))
    threshold_draws = np.quantile(bootstrap.choice(groups["calibration"], (1000, REPS)), .95, axis=1, method="higher")
    rows = []
    for group in ["null_check", "ellis_power"]:
        values = groups[group]
        rejected = int(sum(values > threshold))
        lo, hi = wilson(rejected, REPS)
        sensitivity = np.mean(values[:, None] > threshold_draws[None, :], axis=0)
        rows.append(dict(design=design, N=n, total_incident=24*n, seed=seed, group=group,
            rejected=rejected, repetitions=REPS, rate=rejected/REPS, wilson_low=float(lo), wilson_high=float(hi),
            threshold=threshold, threshold_low=float(np.quantile(threshold_draws, .025)), threshold_high=float(np.quantile(threshold_draws, .975)),
            rate_from_threshold_low=float(np.quantile(sensitivity, .025)), rate_from_threshold_high=float(np.quantile(sensitivity, .975))))
    return rows, {"design": design, "N": n, "seed": seed, "file": filename,
        "unresolved_optimizer_warnings": warnings, "boundary_fits": bounds, "fallback_count": fallbacks,
        "max_deviance_improvement": max(improvements)}


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    source = ROOT/"results/study-002"
    with (source/"spectra.csv").open() as f:
        spectra = list(csv.DictReader(f))
    old = json.loads((source/"summary.json").read_text())
    config = {"author": "Biswajit Jana", "date": "2026-11-10", "N_values": COUNTS,
        "channels": 24, "calibration_repetitions": REPS, "null_check_repetitions": REPS,
        "power_repetitions": REPS, "designs": ["low_band", "broad_band"],
        "base_seed": 20261120, "seed_rule": "base + design_index*1000 + budget_index; bootstrap uses case_seed+1000000",
        "independent_cases": True, "workers": 4, "classification": "synthetic ideal-detector experiment",
        "source_spectra_sha256": hashlib.sha256((source/"spectra.csv").read_bytes()).hexdigest(),
        "fit_bounds": BOUNDS, "threshold_calibration": "pointwise at noiseless comparator projection",
        "interval_scope": "pointwise Wilson over evaluations, conditional on estimated threshold; bootstrap sensitivity separate"}
    (OUT/"design.json").write_text(json.dumps(config, indent=2)+"\n", newline="\n")
    jobs = []
    for d, design in enumerate(config["designs"]):
        s = [r for r in spectra if r["design"] == design]
        k, p, q = (np.array([float(r[key]) for r in s]) for key in ["K", "T_ellis", "T_pt_projection"])
        theta = np.array(old["design_results"][design]["pt_logC_logb"])
        for j, n in enumerate(COUNTS):
            jobs.append((design, n, config["base_seed"]+d*1000+j, k, p, theta, q))
    rates, audits = [], []
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(case, job) for job in jobs]
        for future in as_completed(futures):
            rr, audit = future.result()
            rates.extend(rr)
            audits.append(audit)
            print(f"{len(audits):02d}/50 {audit['design']} N={audit['N']}: power={rr[-1]['rate']:.3f}", flush=True)
    rates.sort(key=lambda r: (r["design"], r["N"], r["group"]))
    audits.sort(key=lambda r: (r["design"], r["N"]))
    for filename, rows in [("rates.csv", rates), ("optimizer-audit.csv", audits)]:
        with (OUT/filename).open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
        "axes.spines.top": False, "axes.spines.right": False, "svg.hashsalt": "dense-power-005"})
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for design, color in [("low_band", "#ad6226"), ("broad_band", "#16798a")]:
        for col, group in enumerate(["ellis_power", "null_check"]):
            r = [v for v in rates if v["design"] == design and v["group"] == group]
            x, y, lo, hi = (np.array([v[key] for v in r]) for key in ["N", "rate", "wilson_low", "wilson_high"])
            axes[col].errorbar(x, y, yerr=[y-lo, hi-y], fmt="o-", markersize=3.5, linewidth=1.2,
                elinewidth=.8, capsize=2, color=color, label=design.replace("_", " "))
            axes[col].set_xscale("log")
            axes[col].set_xlabel("Incident quanta per channel (24 channels)")
            axes[col].set_ylabel("Fraction rejecting the fitted sech² family")
            axes[col].grid(True, alpha=.15)
    axes[0].set(title="Detection power under Ellis truth", ylim=(-.02, 1.035))
    axes[1].set(title="Independent false-positive check", ylim=(0, .12))
    axes[1].axhline(.05, color="#717780", linestyle="--", linewidth=1, label="5% calibration target")
    axes[0].legend(fontsize=9, loc="center left")
    axes[1].legend(fontsize=8)
    fig.suptitle("25 independently simulated budgets · 1,000 evaluations per point · 95% Wilson intervals", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIG/"01-dense-power.png", dpi=220, bbox_inches="tight")
    path = FIG/"01-dense-power.svg"
    fig.savefig(path, bbox_inches="tight", metadata={"Date": None, "Creator": "Biswajit Jana"})
    path.write_text("\n".join(x.rstrip() for x in path.read_text(encoding="utf-8").splitlines())+"\n", encoding="utf-8", newline="\n")
    plt.close(fig)
    summary = {"author": "Biswajit Jana", "experiments": len(jobs)*REPS*3, "channel_counts": len(jobs)*REPS*3*24,
        "unresolved_optimizer_warnings": sum(r["unresolved_optimizer_warnings"] for r in audits),
        "boundary_fits": sum(r["boundary_fits"] for r in audits),
        "max_deviance_improvement": max(r["max_deviance_improvement"] for r in audits),
        "first_tested_broad_budget_lower_Wilson_at_least_90pct": next((r["N"] for r in rates if r["design"] == "broad_band" and r["group"] == "ellis_power" and r["wilson_low"] >= .9), None)}
    (OUT/"summary.json").write_text(json.dumps(summary, indent=2)+"\n", newline="\n")
    paths = [p for folder in (OUT, FIG) for p in sorted(folder.iterdir()) if p.suffix in (".json", ".csv", ".gz", ".png", ".svg") and p.name != "sha256.json"]
    hashes = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (OUT/"sha256.json").write_text(json.dumps(hashes, indent=2)+"\n", newline="\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
