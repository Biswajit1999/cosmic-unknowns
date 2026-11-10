"""Calibration sensitivity and deviance distributions. Author: Biswajit Jana."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from calibration_limits import (required_error, matching_errors, recorded_probability,
    fit_error_budget, global_search_check, profile_kl, bernoulli_kl)
from transmission_inference import pt_transmission, pt_parameters, wilson
from wormhole_scattering import scatter

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/"results/study-004"
FIG = ROOT/"figures/study-004"


def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def csv_out(name, rows):
    with (OUT/name).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def save(fig, name):
    fig.savefig(FIG/(name+".png"), dpi=180, bbox_inches="tight")
    svg = FIG/(name+".svg")
    fig.savefig(svg, bbox_inches="tight", metadata={"Date": None, "Creator": "Biswajit Jana"})
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text(encoding="utf-8").splitlines())+"\n", encoding="utf-8", newline="\n")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
        "axes.spines.top": False, "axes.spines.right": False, "svg.hashsalt": "calibration-study-004"})
    source = ROOT/"results/study-002"
    spectra, trials, rates = (read_csv(source/name) for name in ["spectra.csv", "trials.csv", "rates.csv"])
    prior = json.loads((source/"summary.json").read_text())
    details, profile, witnesses, histogram_data = {}, [], [], []
    colors = {"low_band": "#ad6226", "broad_band": "#16798a"}
    for design in colors:
        data = [r for r in spectra if r["design"] == design]
        k, p = (np.array([float(r[key]) for r in data]) for key in ["K", "T_ellis"])
        nominal = np.array(prior["design_results"][design]["pt_logC_logb"])
        starts = [nominal, pt_parameters(.5, .8), pt_parameters(1., 1.), pt_parameters(1.5, 1.5)]
        theta, delta, attempts = fit_error_budget(k, p, starts)
        independent = global_search_check(k, p, 20261115)
        if independent["delta"] < delta:
            theta, delta = np.array(independent["theta"]), independent["delta"]
        q = pt_transmission(k, theta)
        fp, fn = matching_errors(p, q)
        adjusted = recorded_probability(q, fp, fn)
        # Direct recomputation on every channel, no interpolation.
        refined = np.array([scatter(v, domain=640, rtol=2e-12, atol=2e-14).T for v in k])
        refined_delta = float(max(required_error(refined, q)))
        details[design] = {"theta_logC_logb": theta.tolist(), "achievable_error_budget": delta,
            "classification": "constructive feasible witness; numerical minimum not certified",
            "nominal_projection_error_budget": float(max(required_error(p, pt_transmission(k, nominal)))),
            "max_recorded_probability_difference": float(max(abs(adjusted-p))),
            "refined_probability_max_change": float(max(abs(refined-p))),
            "refined_witness_error_budget": refined_delta,
            "all_SLSQP_attempts": attempts, "independent_differential_evolution": independent}
        for i in range(len(k)):
            witnesses.append(dict(design=design, channel=i, K=k[i], p_ellis=p[i], q_alternative=q[i],
                false_positive=fp[i], false_negative=fn[i], recorded_alternative=adjusted[i],
                refined_p_ellis=refined[i], per_channel_required_error=max(fp[i], fn[i])))
        # Profiles are dimensionless KL; they are not detection probabilities.
        for fraction in np.linspace(0, 1.15, 47):
            d = float(fraction*delta)
            t, loss, ok = profile_kl(k, p, d, [nominal, theta])
            profile.append(dict(design=design, delta=d, fraction_of_witness_budget=fraction,
                KL_sum=loss, expected_fixed_distribution_logLR_at_N10000=2*10000*loss,
                fit_success=int(ok), logC=t[0], logb=t[1]))
        print(f"{design}: feasible delta={delta:.10g}; matched p error={max(abs(adjusted-p)):.3g}", flush=True)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.6))
    for row, design in enumerate(colors):
        for col, n in enumerate([1000, 10000]):
            ax = axes[row, col]
            null = np.array([float(r["deviance"]) for r in trials if r["design"] == design and int(r["N"]) == n and r["group"] == "null_check"])
            truth = np.array([float(r["deviance"]) for r in trials if r["design"] == design and int(r["N"]) == n and r["group"] == "ellis_power"])
            threshold = float(next(r["threshold"] for r in rates if r["design"] == design and int(r["N"]) == n))
            # Common fixed edges for bin-frequency intervals; no extrema-derived bins.
            edges = np.arange(0., 205., 5.)
            assert min(null.min(), truth.min()) >= edges[0] and max(null.max(), truth.max()) < edges[-1]
            for group, values, color in [("Comparator null", null, "#677283"), ("Ellis truth", truth, colors[design])]:
                count, _ = np.histogram(values, edges)
                ax.stairs(count/len(values), edges, color=color, linewidth=2, label=group)
                lo, hi = wilson(count, len(values))
                mid = (edges[:-1]+edges[1:])/2
                ax.errorbar(mid, count/len(values), yerr=[count/len(values)-lo, hi-count/len(values)], fmt="none", ecolor=color, alpha=.4, capsize=1)
                for i in range(len(count)):
                    histogram_data.append(dict(design=design, N=n, group=group, bin_left=edges[i], bin_right=edges[i+1],
                        count=int(count[i]), total=len(values), probability=count[i]/len(values), wilson_low=lo[i], wilson_high=hi[i]))
            ax.axvline(threshold, color="#bc4153", linestyle="--", lw=1.4, label="Independent calibration threshold")
            ax.set(title=f"{design.replace('_', ' ').title()} · N = {n:,}/channel", xlabel="Refitted binomial deviance", ylabel="Fraction per histogram bin")
            ax.set_xlim(0, 200 if design == "broad_band" and n == 10000 else 80)
            ax.legend(fontsize=7)
    fig.suptitle("Ideal-detector experiments · 500 trials/group · pointwise 95% bin-frequency intervals", fontsize=12)
    fig.tight_layout()
    save(fig, "01-deviance-histograms")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for design in colors:
        data = [r for r in profile if r["design"] == design]
        ax = axes[0 if design == "low_band" else 1]
        ax.plot([100*r["delta"] for r in data], [r["expected_fixed_distribution_logLR_at_N10000"] for r in data], color=colors[design], lw=2)
        ax.axvline(100*details[design]["achievable_error_budget"], color="#bc4153", linestyle="--", label="Constructed exact-match budget")
        ax.set(title=design.replace("_", " ").title(), xlabel="Allowed error rate per class and channel (%)", ylabel="2 N × profiled KL sum at N = 10,000")
        ax.legend(fontsize=8)
    fig.suptitle("Calibration flexibility can erase population-level discrimination", fontsize=12)
    fig.tight_layout()
    save(fig, "02-calibration-profile")

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex="col")
    for col, design in enumerate(colors):
        data = [r for r in witnesses if r["design"] == design]
        k = [r["K"] for r in data]
        axes[0, col].plot(k, [r["p_ellis"] for r in data], color="#16798a", label="Ellis, ideal detector")
        axes[0, col].plot(k, [r["q_alternative"] for r in data], "--", color="#ad6226", label="Alternative, ideal detector")
        axes[0, col].plot(k, [r["recorded_alternative"] for r in data], ".", color="#1e2530", label="Alternative + constructed errors")
        axes[1, col].plot(k, [100*r["false_positive"] for r in data], "o-", color="#ad6226", markersize=3, label="False positive")
        axes[1, col].plot(k, [100*r["false_negative"] for r in data], "s-", color="#16798a", markersize=3, label="False negative")
        axes[0, col].set(title=design.replace("_", " ").title(), ylabel="Recorded transmission probability")
        axes[1, col].set(xlabel="Dimensionless wave number K", ylabel="Constructed classification error (%)")
        for ax in axes[:, col]:
            ax.set_xscale("log")
            ax.legend(fontsize=7)
    fig.suptitle("An explicit observational-equivalence witness · assumed errors, not measured calibration", fontsize=12)
    fig.tight_layout()
    save(fig, "03-equivalent-probabilities")
    csv_out("witnesses.csv", witnesses)
    csv_out("kl-profiles.csv", profile)
    csv_out("histogram-bins.csv", histogram_data)
    summary = {"author": "Biswajit Jana", "date": "2026-11-10", "data_type": "analytical and synthetic-model analysis",
        "source_spectra_sha256": hashlib.sha256((source/"spectra.csv").read_bytes()).hexdigest(),
        "source_trials_sha256": hashlib.sha256((source/"trials.csv").read_bytes()).hexdigest(),
        "model": "p_recorded=(1-fn_i)q_i+fp_i(1-q_i), independent errors in [0,delta] per channel",
        "designs": details, "profile_non_success": sum(1-r["fit_success"] for r in profile),
        "interpretation": "Existence of a probability-matching comparator in an explicitly enlarged nuisance family; not physical wormhole existence.",
        "novelty": "Application and reproducible numerical witnesses are candidate contributions; no first-ever claim."}
    (OUT/"summary.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8", newline="\n")
    paths = [p for folder in (OUT, FIG) for p in sorted(folder.iterdir()) if p.suffix in (".csv", ".json", ".png", ".svg") and p.name != "sha256.json"]
    hashes = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (OUT/"sha256.json").write_text(json.dumps(hashes, indent=2)+"\n", newline="\n")
    print(json.dumps({k: {x:y for x,y in v.items() if x != "all_SLSQP_attempts"} for k,v in details.items()}, indent=2))


if __name__ == "__main__":
    main()
