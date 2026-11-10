# Numerical data: study 001

**Author: Biswajit Jana.** All data here are generated from a prescribed theoretical model. They are not telescope observations, detector measurements, or samples from a statistically inferred wormhole population.

Dimensionless conventions: `K = k*a`, `X = x/a`, `xi` is scalar curvature coupling, `ell=0`, `T` and `R` are flux fractions, `gamma = growth_rate*a/c`, and `mu = m*c*a/hbar`.

| File | Contents |
|---|---|
| `spectra.csv` | 105 frequencies × five couplings; flux probabilities, complex amplitudes, boundary distance, and conservation error. |
| `selected-values.csv` | Directly recomputed report values using boundary distance 320. |
| `domain-convergence.csv` | Four frequencies and four boundaries compared with boundary distance 640 and tighter tolerance. |
| `tolerance-convergence.csv` | Three frequencies at three integration tolerances, fixed boundary distance. |
| `critical-scaling.csv` | 44 ratios K/q × four barrier strengths; exact numerical T and weak-barrier approximation. |
| `stability.csv` | Eight couplings; two grid spacings, enlarged/half-size box comparison, extrapolated ground energy and growth rate. |
| `response-equivalence.csv` | 32 frequencies comparing the two-sided geometry and an exact frequency-dependent load. |
| `finite-bandwidth.csv` | All-time energy fractions for four Gaussian spectral centers × five couplings; interpolated main spectrum. |
| `additional-audit.json` | Additional boundary checks, direct packet quadrature, and independent three-spacing eigenvalue convergence. |
| `summary.json` | Run settings, dependency versions, diagnostics and selected output. |
| `sha256.json` | Hashes of generated numerical and figure artifacts. Documentation files are excluded. |

The final digits of CSV entries are retained for reproducibility, not asserted physical precision. Tiny probabilities above one near the exactly transparent case reflect floating-point/integration errors and are not manually clipped. Boundary and resolution diagnostics apply to the sampled cases, not a certified global bound.

The bound-state Richardson estimate uses the leading second-order spacing error. Near-threshold modes require large domains; a positive eigenvalue in a small Dirichlet box would not prove stability on an infinite line. The massive-probe threshold is a spectral shift of the same lowest mode, not a separate nonlinear simulation.

`finite-bandwidth.csv` uses Gaussian spectral power, energy weight K², an explicitly recorded finite integration range, and interpolation in log(K). Direct quadrature comparisons are in `additional-audit.json`. These are not finite-observation-time calculations.

Recreate the numerical products with `scripts/run_wormhole_study.py`, then run `scripts/audit_wormhole_study.py`. Both commands overwrite the generated numerical products; no external data are needed.
