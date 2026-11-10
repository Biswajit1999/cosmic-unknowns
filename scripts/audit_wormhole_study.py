"""Additional checks not inferred from flux conservation. Author: Biswajit Jana."""
import csv
import json
from pathlib import Path
import numpy as np
from scipy.integrate import simpson
from wormhole_scattering import scatter, ground_energy
from artifact_manifest import write_manifest

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/"results"/"study-001"


def main():
    rows = []
    for xi, k in ((1/6, .003), (.4, .003), (.49, .003), (.49, .1), (.495, .0006), (.495, .08)):
        base = scatter(k, xi=xi, domain=160)
        extended = scatter(k, xi=xi, domain=640, rtol=2e-12, atol=2e-14)
        rows.append({"xi": xi, "K": k, "relative_T_change": abs(base.T/extended.T-1),
                     "complex_t_change": abs(base.transmission-extended.transmission)})
    packet_rows = []
    with (OUT/"finite-bandwidth.csv").open(encoding="utf-8") as f:
        packets = [{k:float(v) for k,v in row.items()} for row in csv.DictReader(f)]
    for center, xi in ((.1, 0), (.3, 0), (.03, .49)):
        original = next(r for r in packets if r["center_K"] == center and r["xi"] == xi)
        width = .2*center
        values = []
        for n in (41, 81):
            kgrid = np.linspace(original["integration_K_min"], original["integration_K_max"], n)
            trans = np.array([scatter(float(k), xi=xi).T for k in kgrid])
            weight = kgrid*kgrid*np.exp(-.5*((kgrid-center)/width)**2)
            values.append(float(simpson(weight*trans, x=kgrid)/simpson(weight, x=kgrid)))
        packet_rows.append({"center_K": center, "xi": xi, "direct_41": values[0], "direct_81": values[1],
                            "interpolated": original["energy_fraction"],
                            "relative_interpolation_difference": abs(original["energy_fraction"]/values[1]-1),
                            "relative_quadrature_difference": abs(values[0]/values[1]-1)})
    # Verify expected O(dx^2) eigenvalue convergence before using extrapolation.
    eigen_rows = []
    for xi in (.55, 1.0):
        e = [ground_energy(xi, domain=120, spacing=h)[0] for h in (.08, .04, .02)]
        ratio = (e[0]-e[1])/(e[1]-e[2])
        eigen_rows.append({"xi": xi, "energies": e, "successive_difference_ratio": ratio})
    checks = {"extended_domain": rows, "packet_quadrature": packet_rows,
              "eigenvalue_order": eigen_rows,
              "max_extended_domain_relative_change": max(r["relative_T_change"] for r in rows),
              "max_packet_relative_interpolation_difference": max(r["relative_interpolation_difference"] for r in packet_rows)}
    assert checks["max_extended_domain_relative_change"] < 1e-4
    assert checks["max_packet_relative_interpolation_difference"] < .01
    assert all(3.9 < r["successive_difference_ratio"] < 4.1 for r in eigen_rows)
    (OUT/"additional-audit.json").write_text(json.dumps(checks, indent=2)+"\n", encoding="utf-8")
    write_manifest()
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
