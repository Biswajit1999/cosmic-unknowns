"""Reproduce ideal counting experiment and uncertainty plots. Author: Biswajit Jana."""
import csv
import hashlib
import json
from pathlib import Path
import platform
import time
import numpy as np
import scipy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import xlogy
from wormhole_scattering import scatter
from transmission_inference import (wilson, fit_pt, pt_transmission, physical_parameters,
                                    deviance, numerical_pt, pt_parameters, BOUNDS)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/"results/study-002"
FIG = ROOT/"figures/study-002"
CONFIG = {
    "author": "Biswajit Jana", "classification": "synthetic controlled counting experiment",
    "date": "2026-11-10", "seed": 20261110, "channels": 24,
    "designs": {"low_band": [0.03, 0.2], "broad_band": [0.05, 1.5]},
    "trials_per_channel": [100, 1000, 10000, 100000],
    "calibration_repetitions": 500, "evaluation_repetitions": 500,
    "nominal_false_positive": 0.05, "threshold_bootstraps": 1000,
    "truth": {"xi": 0., "ell": 0, "scale": 1., "domain": 320., "rtol": 2e-11, "atol": 2e-13},
    "pt_fit_bounds_logC_logb": BOUNDS,
    "assumptions": ["independent incident scalar quanta", "known incident counts and wave numbers",
                    "perfect transmission classification", "no absorption or background counts",
                    "fixed massless minimal spherical probe", "no astronomical observations"]}


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False)+"\n", encoding="utf-8", newline="\n")


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def save_figure(fig, name):
    fig.savefig(FIG/(name+".png"), dpi=190, bbox_inches="tight")
    fig.savefig(FIG/(name+".svg"), bbox_inches="tight", metadata={"Date": None, "Creator": "Biswajit Jana"})
    svg = FIG/(name+".svg")
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text(encoding="utf-8").splitlines())+"\n",
                   encoding="utf-8", newline="\n")
    plt.close(fig)


def main():
    started = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    # The specification is written before any counts are drawn.
    write_json(OUT/"design.json", CONFIG)
    rng = np.random.default_rng(CONFIG["seed"])
    # Independent stream for displayed examples: no selection from power trials.
    example_rng = np.random.default_rng(CONFIG["seed"]+1)
    threshold_rng = np.random.default_rng(CONFIG["seed"]+2)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.hashsalt": "cosmic-unknowns-study-002"})
    colors = {"low_band": "#b46120", "broad_band": "#157c91"}
    spectra, numerical, trials, rates, examples, design_summary, optimizer_audit = [], [], [], [], [], {}, []
    cached = {}
    for design, limits in CONFIG["designs"].items():
        k = np.geomspace(*limits, CONFIG["channels"])
        truth = np.array([scatter(v, domain=320, rtol=2e-11, atol=2e-13).T for v in k])
        theta, success, boundary, grad = fit_pt(k, truth, multistart=True)
        if not success or boundary:
            raise RuntimeError("Noise-free comparator fit failed or reached a boundary")
        null = pt_transmission(k, theta)
        height, alpha = physical_parameters(theta)
        kl = xlogy(truth, truth/null)+xlogy(1-truth, (1-truth)/(1-null))
        z = np.exp(theta[1])*k
        jac = np.column_stack([-np.ones_like(k), 2*z/np.tanh(z)])
        fisher = jac.T@((null*(1-null))[:, None]*jac)
        design_summary[design] = {"pt_logC_logb": theta.tolist(), "height": height, "alpha": alpha,
            "sum_Bernoulli_KL": float(kl.sum()), "max_absolute_probability_mismatch": float(max(abs(truth-null))),
            "Fisher_condition_number": float(np.linalg.cond(fisher)), "Fisher_per_trial": fisher.tolist()}
        for i in range(len(k)):
            spectra.append(dict(design=design, channel=i, K=k[i], T_ellis=truth[i],
                                T_pt_projection=null[i], KL_ellis_to_pt=kl[i]))
        for i in range(0, len(k), 3):
            t160 = scatter(k[i], domain=160, rtol=2e-11, atol=2e-13).T
            t640 = scatter(k[i], domain=640, rtol=2e-12, atol=2e-14).T
            tight = scatter(k[i], domain=320, rtol=2e-12, atol=2e-14).T
            numerical.append(dict(design=design, K=k[i], T160=t160, T320=truth[i], T640=t640,
                domain_change_160_320=abs(t160-truth[i]), domain_change_320_640=abs(truth[i]-t640),
                tolerance_change=abs(tight-truth[i]),
                binomial_sd_at_N100000=float(np.sqrt(truth[i]*(1-truth[i])/100000))))
        cached[design] = (k, truth, theta, null, kl)
        print(f"Computed {design}: KL sum={kl.sum():.6g}, max mismatch={max(abs(truth-null)):.6g}", flush=True)
        for n in CONFIG["trials_per_channel"]:
            deviances = {}
            for group, p in [("calibration", null), ("null_check", null), ("ellis_power", truth)]:
                reps = CONFIG["calibration_repetitions"] if group == "calibration" else CONFIG["evaluation_repetitions"]
                ds = []
                for rep in range(reps):
                    y = rng.binomial(n, p)
                    fitted, ok, bound, gradient = fit_pt(k, y/n, start=theta)
                    # All flagged optimizations and a fixed subset get independent starts.
                    audit = rep < 5 or not ok or bound or gradient > 1e-5
                    initial_d = deviance(y, n, pt_transmission(k, fitted))
                    if audit:
                        other, other_ok, other_bound, other_gradient = fit_pt(k, y/n, start=fitted, multistart=True)
                        other_d = deviance(y, n, pt_transmission(k, other))
                        optimizer_audit.append(dict(design=design, N=n, group=group, repetition=rep,
                            initial_deviance=initial_d, multistart_deviance=other_d,
                            improvement=initial_d-other_d))
                        if other_d <= initial_d:
                            fitted, ok, bound, gradient = other, other_ok, other_bound, other_gradient
                    d = deviance(y, n, pt_transmission(k, fitted))
                    ds.append(d)
                    row = dict(design=design, N=n, group=group, repetition=rep, deviance=d,
                        logC=fitted[0], logb=fitted[1], optimizer_success=int(ok), boundary=int(bound), gradient=gradient)
                    row.update({f"count_{j:02d}": int(count) for j, count in enumerate(y)})
                    trials.append(row)
                deviances[group] = np.array(ds)
            cal = deviances["calibration"]
            threshold = float(np.quantile(cal, 0.95, method="higher"))
            thresholds = np.quantile(threshold_rng.choice(cal, (CONFIG["threshold_bootstraps"], len(cal))),
                                      0.95, axis=1, method="higher")
            for group in ("null_check", "ellis_power"):
                ds = deviances[group]
                rejected = int(np.count_nonzero(ds > threshold))
                low, high = wilson(rejected, len(ds))
                varying_rates = np.mean(ds[:, None] > thresholds[None, :], axis=0)
                rates.append(dict(design=design, N=n, total_incident=24*n, group=group,
                    rejected=rejected, repetitions=len(ds), rate=rejected/len(ds), wilson_low=float(low),
                    wilson_high=float(high), threshold=threshold,
                    threshold_bootstrap_low=float(np.quantile(thresholds, .025)),
                    threshold_bootstrap_high=float(np.quantile(thresholds, .975)),
                    rate_from_threshold_low=float(np.quantile(varying_rates, .025)),
                    rate_from_threshold_high=float(np.quantile(varying_rates, .975))))
            print(f"  N={n}: calibrated threshold={threshold:.3f}; power={rates[-1]['rate']:.3f}", flush=True)
        n = 10000
        y = example_rng.binomial(n, truth)
        fit, _, _, _ = fit_pt(k, y/n, start=theta, multistart=True)
        low, high = wilson(y, n)
        for i in range(len(k)):
            examples.append(dict(design=design, N=n, channel=i, K=k[i], transmitted=int(y[i]),
                reflected=int(n-y[i]), fraction=y[i]/n, wilson_low=low[i], wilson_high=high[i],
                T_ellis=truth[i], T_pt_fit=pt_transmission(k[i], fit)))

    check_rows = []
    for h, a in [(0.1, 1.), (0.25, 1.), (1., 1.), (0.8, 1.2)]:
        for k0 in [0.05, 0.3, 1., 1.5]:
            analytic = float(pt_transmission(k0, pt_parameters(h, a)))
            numeric = numerical_pt(k0, h, a)
            check_rows.append(dict(height=h, alpha=a, K=k0, analytic=analytic, numerical=numeric,
                                   absolute_error=abs(analytic-numeric)))
    for name, data in [("spectra", spectra), ("numerical-audit", numerical), ("trials", trials),
                       ("rates", rates), ("example-counts", examples), ("pt-validation", check_rows),
                       ("optimizer-audit", optimizer_audit)]:
        write_csv(OUT/(name+".csv"), data)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex="col", gridspec_kw={"height_ratios": [2, 1]})
    for col, design in enumerate(CONFIG["designs"]):
        data = [r for r in examples if r["design"] == design]
        x, y, lo, hi, p, q = (np.array([r[key] for r in data]) for key in
                             ["K", "fraction", "wilson_low", "wilson_high", "T_ellis", "T_pt_fit"])
        axes[0, col].plot(x, p, color="#157c91", label="Ellis generating probability")
        axes[0, col].plot(x, q, "--", color="#b46120", label="Fitted sech² barrier")
        axes[0, col].errorbar(x, y, yerr=[y-lo, hi-y], fmt=".", color="#303642", capsize=2,
                             label="Synthetic counts: 95% Wilson")
        axes[1, col].errorbar(x, y-q, yerr=[y-lo, hi-y], fmt=".", color="#303642", capsize=2)
        axes[1, col].plot(x, p-q, color="#157c91", label="Expected residual")
        axes[1, col].axhline(0, color="gray", lw=0.8)
        axes[0, col].set_title(design.replace("_", " ").title()+" · 10,000 incidents/channel")
        axes[0, col].set_ylabel("Transmission fraction")
        axes[1, col].set_ylabel("Fraction − fitted probability")
        axes[1, col].set_xlabel("Dimensionless wave number K")
        for ax in axes[:, col]:
            ax.set_xscale("log")
        axes[0, col].legend(fontsize=8)
    fig.suptitle("Simulated counting experiment · pointwise 95% intervals · no observational data", fontsize=12)
    fig.tight_layout()
    save_figure(fig, "01-counts-and-residuals")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.7))
    for design in CONFIG["designs"]:
        for col, group in enumerate(("ellis_power", "null_check")):
            data = [r for r in rates if r["design"] == design and r["group"] == group]
            n, y, low, high = (np.array([r[key] for r in data]) for key in ["N", "rate", "wilson_low", "wilson_high"])
            axes[col].errorbar(n, y, yerr=[y-low, high-y], fmt="o-", capsize=4,
                color=colors[design], label=design.replace("_", " "))
            axes[col].set_xscale("log")
            axes[col].set_xlabel("Incident quanta per channel (24 channels)")
            axes[col].set_ylabel("Fraction rejecting the fitted sech² family")
            axes[col].legend()
    axes[0].set_title("Detection power under Ellis truth")
    axes[0].set_ylim(-.03, 1.04)
    axes[1].set_title("Independent null calibration check")
    axes[1].axhline(.05, color="gray", linestyle="--", label="Target 5%")
    axes[1].set_ylim(0, .16)
    fig.suptitle("500 evaluation trials per point · 95% Wilson Monte Carlo intervals", fontsize=12)
    fig.tight_layout()
    save_figure(fig, "02-power-and-false-positives")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.7))
    x = np.linspace(-5, 5, 501)
    axes[0].plot(x, (1+x*x)**-2, color="#303642", lw=2, label="Ellis, minimal s wave")
    for design, (k, p, theta, q, kl) in cached.items():
        h, a = physical_parameters(theta)
        axes[0].plot(x, h/np.cosh(a*x)**2, "--", color=colors[design], label=design.replace("_", " ")+" fit")
        axes[1].plot(k, 2*10000*kl, "o-", color=colors[design], label=design.replace("_", " "))
    axes[0].set(xlabel="Dimensionless position X", ylabel="Effective potential U(X)", title="Different barriers can have similar spectra")
    axes[1].set(xlabel="Dimensionless wave number K", ylabel="2 N × Bernoulli KL per channel", title="Expected mismatch information at N = 10,000")
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    for ax in axes:
        ax.legend(fontsize=8)
    fig.tight_layout()
    save_figure(fig, "03-potential-and-information")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.7))
    for design in CONFIG["designs"]:
        data = [r for r in numerical if r["design"] == design]
        k = np.array([r["K"] for r in data])
        for key, style, label in [("domain_change_320_640", "o-", "domain change"),
                                   ("tolerance_change", "s:", "tolerance change")]:
            y = np.array([r[key] for r in data])
            axes[0].loglog(k, y, style, color=colors[design], label=design.replace("_", " ")+" "+label)
        ratio = [r["domain_change_320_640"]/r["binomial_sd_at_N100000"] for r in data]
        axes[1].loglog(k, ratio, "o-", color=colors[design], label=design.replace("_", " "))
    axes[0].set(xlabel="K", ylabel="Absolute change in T", title="Numerical changes, not statistical error bars")
    axes[1].set(xlabel="K", ylabel="Domain change / counting standard deviation", title="Numerical scale vs N = 100,000 shot noise")
    for ax in axes:
        ax.legend(fontsize=8)
    fig.tight_layout()
    save_figure(fig, "04-numerical-error-audit")

    summary = {"author": CONFIG["author"], "classification": CONFIG["classification"],
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        "design_results": design_summary, "total_synthetic_experiments": len(trials),
        "total_synthetic_channel_counts": len(trials)*24,
        "optimizer_non_success": sum(1-r["optimizer_success"] for r in trials),
        "boundary_fit_count": sum(r["boundary"] for r in trials),
        "max_multistart_deviance_improvement": max(r["improvement"] for r in optimizer_audit),
        "max_pt_ode_absolute_error": max(r["absolute_error"] for r in check_rows),
        "max_domain_absolute_change": max(r["domain_change_320_640"] for r in numerical),
        "max_tolerance_absolute_change": max(r["tolerance_change"] for r in numerical),
        "max_domain_change_in_counting_sd": max(r["domain_change_320_640"]/r["binomial_sd_at_N100000"] for r in numerical),
        "elapsed_seconds": time.time()-started}
    write_json(OUT/"summary.json", summary)
    paths = [p for folder in (OUT, FIG) for p in sorted(folder.iterdir())
             if p.suffix in (".json", ".csv", ".png", ".svg") and p.name != "sha256.json"]
    for path in paths:
        if path.suffix != ".png":
            path.write_bytes(path.read_bytes().replace(b"\r\n", b"\n"))
    hashes = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    write_json(OUT/"sha256.json", hashes)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
