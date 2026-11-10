"""Reproducible numerical study and figures. Author: Biswajit Jana."""
from pathlib import Path
import csv
import json
import platform
import time
import numpy as np
import scipy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from wormhole_scattering import (scatter, potential, ground_energy,
                                 halfline_reflection, low_energy_coefficient)
from artifact_manifest import write_manifest

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/"results"/"study-001"
FIG = ROOT/"figures"/"study-001"
COLORS = ["#165b7a", "#d17b30", "#418571", "#8a548e", "#b94b4b"]


def write_csv(name, rows):
    with (OUT/name).open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def save(fig, name):
    fig.text(0.99, 0.012, "Biswajit Jana | prescribed geometry | numerical model data",
             ha="right", fontsize=8, color="#62676e")
    for ext in ("png", "svg"):
        path = FIG/f"{name}.{ext}"
        extra = {"metadata": {"Date": None}} if ext == "svg" else {}
        fig.savefig(path, dpi=190, facecolor="white", **extra)
        if ext == "svg":
            path.write_text("\n".join(s.rstrip() for s in path.read_text(encoding="utf-8").splitlines())+"\n", encoding="utf-8")
    plt.close(fig)


def main():
    start = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.alpha": 0.18,
                         "svg.hashsalt": "cosmic-unknowns-study-001"})
    kappas = np.geomspace(0.003, 3, 105)
    xis = [0.0, 1/6, 0.4, 0.49, 0.5]
    spectra = []
    for xi in xis:
        for k in kappas:
            s = scatter(float(k), xi=xi)
            spectra.append({"xi": xi, "K": float(k), "T": s.T, "R": s.R,
                            "t_real": s.transmission.real, "t_imag": s.transmission.imag,
                            "r_real": s.reflection.real, "r_imag": s.reflection.imag,
                            "flux_error": s.flux_error, "domain": s.domain})
        print(f"Spectrum xi={xi:.6g} complete", flush=True)
    write_csv("spectra.csv", spectra)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.9))
    fig.suptitle("Signal transmission depends on the field, not just the throat", fontsize=15)
    x = np.linspace(-4, 4, 500)
    for xi, color in zip(xis, COLORS):
        data = [r for r in spectra if r["xi"] == xi]
        label = f"xi = {xi:.4g}"
        axes[0].plot(x, potential(x, xi), color=color, label=label)
        axes[1].semilogx(kappas, [r["T"] for r in data], color=color, label=label)
    axes[0].set(xlabel="Proper distance X = x/a", ylabel="Effective potential a² V", title="Spherical scalar mode")
    axes[1].set(xlabel="Dimensionless wave number K = ka", ylabel="Transmitted flux fraction T",
                ylim=(-0.02, 1.04), title="Massless probe; outgoing boundaries at X = ±160")
    axes[1].legend(fontsize=9)
    fig.tight_layout(rect=(0, 0.05, 1, 0.93))
    save(fig, "01-transmission")

    # Domain, integration tolerance, and analytic low-energy comparisons are independent.
    convergence = []
    for k in (0.003, 0.03, 0.3, 2.0):
        reference = scatter(k, domain=640, rtol=2e-12, atol=2e-14)
        for domain in (40, 80, 160, 320):
            s = scatter(k, domain=domain)
            convergence.append({"K": k, "domain": domain, "T": s.T,
                                "reference_T": reference.T,
                                "relative_T_error": abs(s.T/reference.T-1),
                                "complex_t_error": abs(s.transmission-reference.transmission),
                                "flux_error": s.flux_error})
    write_csv("domain-convergence.csv", convergence)
    tolerance_rows = []
    for k in (0.03, 0.3, 2.0):
        ref = scatter(k, domain=160, rtol=2e-12, atol=2e-14)
        for rtol in (2e-7, 2e-9, 2e-10):
            s = scatter(k, domain=160, rtol=rtol, atol=rtol/100)
            tolerance_rows.append({"K": k, "rtol": rtol, "relative_T_error": abs(s.T/ref.T-1),
                                   "complex_t_error": abs(s.transmission-ref.transmission)})
    write_csv("tolerance-convergence.csv", tolerance_rows)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.9))
    fig.suptitle("Two independent checks on the scattering calculation", fontsize=15)
    for xi, color in zip(xis[:3], COLORS):
        rows = [r for r in spectra if r["xi"] == xi and r["K"] < 0.12]
        B = low_energy_coefficient(xi)
        axes[0].semilogx([r["K"] for r in rows], [r["T"]/(4*r["K"]**2/B**2) for r in rows],
                         label=f"xi = {xi:.4g}", color=color)
    axes[0].axhline(1, color="black", ls="--", lw=1, label="Analytical limit")
    axes[0].set(xlabel="K = ka", ylabel="Numerical T / leading low-energy T", title="K → 0 benchmark")
    axes[0].legend(fontsize=9)
    for k, color in zip((0.003, 0.03, 0.3, 2.0), COLORS):
        rows = [r for r in convergence if r["K"] == k]
        axes[1].loglog([r["domain"] for r in rows], [max(r["relative_T_error"], 1e-15) for r in rows],
                       "o-", color=color, label=f"K = {k:g}")
    axes[1].set(xlabel="Boundary distance L/a", ylabel="Relative transmission difference",
                title="Compared with L/a = 640 and tighter tolerance")
    axes[1].legend(fontsize=9)
    fig.tight_layout(rect=(0, 0.05, 1, 0.93))
    save(fig, "02-validation")
    print("Scattering validation complete", flush=True)

    # Weak-barrier double scaling: K and q=1-2xi are reduced together.
    scaling = []
    ratios = np.geomspace(0.06, 8, 44)
    for q in (1.0, 0.2, 0.05, 0.01):
        for ratio in ratios:
            k = q*ratio
            s = scatter(k, xi=(1-q)/2)
            limit = ratio**2/(ratio**2+(np.pi/4)**2)
            scaling.append({"q": q, "K": k, "K_over_q": ratio,
                            "T": s.T, "weak_barrier_limit": limit,
                            "absolute_error": abs(s.T-limit)})
    write_csv("critical-scaling.csv", scaling)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.9))
    fig.suptitle("Near transparency, frequency and curvature coupling must be considered together", fontsize=14)
    for q, color in zip((1, .2, .05, .01), COLORS):
        rows = [r for r in scaling if r["q"] == q]
        axes[0].semilogx(ratios, [r["T"] for r in rows], color=color, label=f"q = {q:g}")
        axes[1].loglog(ratios, [max(r["absolute_error"], 1e-12) for r in rows], color=color)
    axes[0].semilogx(ratios, ratios**2/(ratios**2+(np.pi/4)**2), "k--", label="Weak-barrier limit")
    axes[0].set(xlabel="K / q, where q = 1 − 2xi", ylabel="T", ylim=(0, 1.02))
    axes[0].legend(fontsize=9)
    axes[1].set(xlabel="K / q", ylabel="Absolute difference from limit")
    fig.tight_layout(rect=(0, 0.05, 1, 0.93))
    save(fig, "03-critical-scaling")

    stability = []
    for xi in (.505, .51, .52, .55, .6, .7, .8, 1.0):
        g = 2*xi-1
        estimate = np.pi*g/4
        domain = max(80.0, 14/estimate)
        coarse, _, _ = ground_energy(xi, domain=domain, spacing=.04)
        fine, dx, n = ground_energy(xi, domain=domain, spacing=.02)
        halfbox, _, _ = ground_energy(xi, domain=domain/2, spacing=.02)
        extrapolated = (4*fine-coarse)/3  # leading O(dx^2) Richardson removal
        gamma = np.sqrt(-extrapolated)
        row = {"xi": xi, "g": g, "domain": domain, "spacing": dx, "grid_points": n,
               "energy_coarse": coarse, "energy_fine": fine, "energy_richardson": extrapolated,
               "gamma": float(gamma), "leading_gamma": estimate,
               "second_order_gamma": estimate-3*np.pi*g*g/16,
               "relative_grid_change": abs(fine/coarse-1),
               "relative_box_change": abs(fine/halfbox-1),
               "critical_dimensionless_mass": float(gamma)}
        stability.append(row)
        print(f"Stability xi={xi:g}: gamma={gamma:.8g}", flush=True)
    write_csv("stability.csv", stability)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.9))
    fig.suptitle("Probe-field instability: a separate question from transmission", fontsize=15)
    g = np.array([r["g"] for r in stability])
    gamma = np.array([r["gamma"] for r in stability])
    axes[0].loglog(g, gamma, "o-", color=COLORS[0], label="Numerical + grid extrapolation")
    axes[0].loglog(g, np.pi*g/4, "--", color=COLORS[1], label="Leading weak-coupling term")
    axes[0].loglog(g[g<=.2], (np.pi*g/4-3*np.pi*g*g/16)[g<=.2], ":", color=COLORS[2],
                   label="Second-order benchmark (g ≤ 0.2)")
    axes[0].set(xlabel="g = 2xi − 1 > 0", ylabel="Growth rate gamma = sigma a/c")
    axes[0].legend(fontsize=8)
    xx = np.r_[.5, [r["xi"] for r in stability]]
    yy = np.r_[0, gamma]
    axes[1].plot(xx, yy, "o-", color=COLORS[0])
    axes[1].fill_between(xx, 0, yy, color=COLORS[4], alpha=.18, label="Growing spherical probe mode")
    axes[1].fill_between(xx, yy, .7, color=COLORS[2], alpha=.12, label="No negative probe eigenvalue")
    axes[1].set(xlabel="Curvature coupling xi", ylabel="Dimensionless scalar mass mu = m c a / hbar",
                ylim=(0, .7), title="Mass shifts the probe spectrum; background stability untested")
    axes[1].legend(fontsize=8, loc="upper left")
    fig.tight_layout(rect=(0, .05, 1, .93))
    save(fig, "04-stability")

    equivalence = []
    for k in np.geomspace(.02, 2, 32):
        s = scatter(float(k))
        rload, z = halfline_reflection(float(k))
        equivalence.append({"K": k, "full_R": s.R, "load_R": abs(rload)**2,
                            "complex_r_difference": abs(rload-s.reflection),
                            "impedance_real": z.real, "impedance_imag": z.imag})
    write_csv("response-equivalence.csv", equivalence)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.9))
    fig.suptitle("One-sided reflection alone does not uniquely identify the interior", fontsize=15)
    ek = [r["K"] for r in equivalence]
    axes[0].semilogx(ek, [r["full_R"] for r in equivalence], color=COLORS[0], label="Two-sided geometry")
    axes[0].semilogx(ek, [r["load_R"] for r in equivalence], "o", ms=3, mfc="none", color=COLORS[1],
                     label="Identical exterior + exact frequency-dependent load")
    axes[0].set(xlabel="K = ka", ylabel="Reflected flux fraction R")
    axes[0].legend(fontsize=8)
    axes[1].loglog(ek, [max(r["complex_r_difference"], 1e-16) for r in equivalence], color=COLORS[2])
    axes[1].set(xlabel="K = ka", ylabel="Absolute complex reflection difference",
                title="Mathematical equivalence; no material substitute constructed")
    fig.tight_layout(rect=(0, .05, 1, .93))
    save(fig, "05-identifiability")

    # Finite-bandwidth energy fractions from the computed spectra (not detector data).
    packets = []
    for center in (.03, .1, .3, 1.0):
        width = .2*center
        kk = np.linspace(max(.003, center-5*width), min(3, center+5*width), 2001)
        weight = kk**2*np.exp(-.5*((kk-center)/width)**2)
        for xi in xis:
            rows = [r for r in spectra if r["xi"] == xi]
            tt = np.interp(np.log(kk), np.log(kappas), [r["T"] for r in rows])
            fraction = np.trapezoid(weight*tt, kk)/np.trapezoid(weight, kk)
            packets.append({"center_K": center, "spectral_width_K": width, "xi": xi,
                            "energy_fraction": float(fraction),
                            "integration_K_min": float(kk[0]), "integration_K_max": float(kk[-1])})
    write_csv("finite-bandwidth.csv", packets)

    low = [r for r in spectra if r["xi"] == 0 and r["K"] < .025]
    slope, intercept = np.polyfit(np.log([r["K"] for r in low]), np.log([r["T"] for r in low]), 1)
    selected = []
    for xi in (0, 1/6, .49, .5):
        for k in (.01, .1, 1.0):
            s = scatter(k, xi=xi, domain=320)
            selected.append({"xi": xi, "K": k, "T": s.T})
    write_csv("selected-values.csv", selected)
    summary = {"author": "Biswajit Jana", "study_date": "2026-11-10",
               "data_class": "Deterministic numerical model data; no astronomical observations",
               "model": "Fixed Ellis geometry; real scalar probe; spherical sector ell=0",
               "settings": {"domain": 160, "rtol": 2e-10, "atol": 2e-12,
                            "frequency_count": len(kappas), "couplings": xis},
               "max_flux_error": max(r["flux_error"] for r in spectra),
               "max_domain_relative_T_error_at_160": max(r["relative_T_error"] for r in convergence if r["domain"]==160),
               "max_load_equivalence_error": max(r["complex_r_difference"] for r in equivalence),
               "minimal_low_K_fitted_exponent": float(slope),
               "minimal_low_K_fitted_coefficient": float(np.exp(intercept)),
               "minimal_exact_leading_coefficient": float(4/np.pi**2),
               "max_relative_eigenvalue_grid_change": max(r["relative_grid_change"] for r in stability),
               "max_relative_eigenvalue_box_change": max(r["relative_box_change"] for r in stability),
               "selected_values": selected, "stability": stability,
               "environment": {"python": platform.python_version(), "numpy": np.__version__,
                                "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
               "elapsed_seconds": time.perf_counter()-start}
    (OUT/"summary.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    write_manifest()
    print(json.dumps({k:summary[k] for k in ("max_flux_error", "max_domain_relative_T_error_at_160",
                     "minimal_low_K_fitted_exponent", "max_load_equivalence_error", "elapsed_seconds")}, indent=2), flush=True)


if __name__ == "__main__":
    main()
