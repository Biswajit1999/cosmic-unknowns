"""Reproduce a known static wormhole geometry; no physical existence claim.

Author: Biswajit Jana
"""
from pathlib import Path
import csv
import json
import platform
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
C_LIGHT = 299792458.0  # m/s, exact SI definition
G_NEWTON = 6.67430e-11  # m^3 kg^-1 s^-2, fixed nominal value for illustration
STRESS_FACTOR = C_LIGHT**4 / (8 * np.pi * G_NEWTON)


def radius(proper_distance, throat_radius):
    if throat_radius <= 0:
        raise ValueError("The throat radius must be positive.")
    return np.sqrt(np.asarray(proper_distance)**2 + throat_radius**2)


def stresses(areal_radius, throat_radius):
    """Return orthonormal energy density, radial and tangential pressures in SI."""
    r = np.asarray(areal_radius, dtype=float)
    if throat_radius <= 0 or np.any(r < throat_radius):
        raise ValueError("Require r >= r0 > 0.")
    magnitude = STRESS_FACTOR * throat_radius**2 / r**4
    return -magnitude, -magnitude, magnitude


def main():
    figures = ROOT / "figures"
    results = ROOT / "results"
    figures.mkdir(exist_ok=True)
    results.mkdir(exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "svg.hashsalt": "cosmic-unknowns-baseline-v1",
                         "axes.spines.top": False, "axes.spines.right": False})
    fig = plt.figure(figsize=(13, 6.8), facecolor="#faf9f6")
    left = fig.add_subplot(121, projection="3d", facecolor="#faf9f6")
    right = fig.add_subplot(122, facecolor="#faf9f6")
    z = np.linspace(-1.3, 1.3, 90)
    angle = np.linspace(0, 2*np.pi, 110)
    zz, aa = np.meshgrid(z, angle)
    rr = np.cosh(zz)
    left.plot_surface(rr*np.cos(aa), rr*np.sin(aa), zz,
                      color="#277f85", alpha=0.70, linewidth=0, rcount=65, ccount=65)
    left.plot(np.cos(angle), np.sin(angle), np.zeros_like(angle),
              color="#bd542e", linewidth=3)
    left.set(xlabel=r"$X/r_0$", ylabel=r"$Y/r_0$", zlabel=r"$z/r_0$",
             title="A known theoretical throat")
    left.set_box_aspect((1, 1, 0.8))
    left.view_init(elev=18, azim=35)
    left.text2D(0.05, -0.04, "Two-dimensional spatial slice; z is an embedding coordinate.",
                transform=left.transAxes, fontsize=9)
    u = np.linspace(-4, 4, 1001)
    nec = -2 / (1 + u*u)**2
    right.axhline(0, color="#555555", linewidth=1)
    right.plot(u, nec, color="#bd542e", linewidth=2.8)
    right.fill_between(u, nec, 0, color="#bd542e", alpha=0.12)
    right.annotate("Throat: negative radial null energy", xy=(0, -2),
                   xytext=(-3.7, -1.25), fontsize=10,
                   arrowprops={"arrowstyle": "->", "color": "#555555"})
    right.set(xlabel=r"Proper distance $l/r_0$",
              ylabel=r"$(\epsilon+p_r)\,/\,[c^4/(8\pi G r_0^2)]$",
              ylim=(-2.25, 0.25), title="The matter this geometry requires")
    right.grid(alpha=0.15)
    fig.suptitle("A wormhole geometry is only the beginning", fontsize=19, y=0.97)
    fig.text(0.5, 0.89, r"Baseline: $\Phi=0$, $b(r)=r_0^2/r$. Known model, reproduced here.",
             ha="center", fontsize=12)
    fig.text(0.5, 0.05, "Calculated illustration • No observation, realizable matter source, or stability proof",
             ha="center", fontsize=11, color="#555555")
    fig.subplots_adjust(left=0.04, right=0.97, top=0.79, bottom=0.21, wspace=0.25)
    for ext in ("png", "svg"):
        output = figures / f"wormhole-baseline.{ext}"
        options = {"metadata": {"Date": None}} if ext == "svg" else {}
        fig.savefig(output, dpi=180, facecolor=fig.get_facecolor(), **options)
        if ext == "svg":
            output.write_text("\n".join(line.rstrip() for line in output.read_text(encoding="utf-8").splitlines())+"\n", encoding="utf-8")
    plt.close(fig)
    rows = []
    for r0 in (1e-3, 1, 1e3, 1e6):
        energy, radial, tangential = stresses(r0, r0)
        rows.append({"throat_radius_m": r0, "energy_density_J_m3": float(energy),
                     "radial_pressure_Pa": float(radial),
                     "tangential_pressure_Pa": float(tangential),
                     "radial_null_energy_J_m3": float(energy+radial)})
    with (results / "throat-scales.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    record = {"classification": "Known static geometry reproduced; not new physics",
              "model": "Phi=0, b(r)=r0^2/r; no matter action specified",
              "constants": {"c_m_s": C_LIGHT, "G_m3_kg_s2": G_NEWTON},
              "assumptions": ["Einstein equations, Lambda=0", "static spherical symmetry",
                              "prescribed stress-energy, not a quantum-state construction"],
              "not_computed": ["stability", "formation", "quantum energy inequalities",
                               "signal transmission", "observational inference"],
              "environment": {"python": platform.python_version(), "numpy": np.__version__,
                              "matplotlib": matplotlib.__version__},
              "scales": rows}
    (results / "baseline.json").write_text(json.dumps(record, indent=2)+"\n", encoding="utf-8")
    print("Created figures/wormhole-baseline.png and .svg; results/baseline.json and throat-scales.csv")


if __name__ == "__main__":
    main()
