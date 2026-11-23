# Study 012: explicit action and instability audit

**Author: Biswajit Jana · 23 November 2026**

[Report](../../reports/014-explicit-action-and-instability.md) · [Script](../../scripts/scalar_tensor_completion.py)

This is an application of an invertible transformation of the known Ellis phantom action. It supports the background mathematically but retains the phantom kinetic sign and a mapped instability. No healthy source or observed wormhole is claimed.

| File | Contents |
|---|---|
| `couplings.csv` | 1,803 rows: three β values and 601 X values. `phi` = atan(X), `F` is the conformal factor, `curvature_coupling` = h = F⁻², `kinetic_coupling` = z = 2F⁻²+6β²F⁻⁴, and `einstein_kinetic` = z/h−1.5(h_phi/h)² = +2. All dimensionless. |
| `instability.csv` | Twelve finite-difference eigenproblems: half-domain, grid spacing, interior-point count, lowest discrete eigenvalue, growth rate and infinity-norm eigen-equation residual. |
| `unstable_mode.csv` | 4,799 interior points of the L = 30, spacing 0.0125 calculation: coordinate X, normalized mode and published regular master potential −3/(1+X²)². Endpoint values are zero. |
| `verification.json` | 35 background point checks at 70-digit precision, 140 diagonal metric residuals, 35 scalar residuals, kinetic coefficients, numerical stability results, source hashes and interpretation. Long residual strings preserve high-precision output. |
| `sha256.json` | Byte hashes for all CSV, verification output and two PNG/SVG figure pairs. |

The action convention is h R + z(grad phi)² with signature (−,+,+,+). A positive Einstein-frame scalar coefficient is the phantom sign. The source action depends on β; changing β does not exhibit a family in one fixed Jordan action. Background checks use independent metric derivatives and action-coupling derivatives. The ghost-sign identity is analytic; finite checks support implementation rather than proving a global theorem.

The master potential is reproduced from González, Guzmán and Sarbach, [arXiv:0806.0608v2, equation 32 at gamma_1 = 0](https://arxiv.org/html/0806.0608v2). This is a coupled background perturbation operator, not an extra probe field. The coordinate-time growth rate is dimensionless; physical growth is gamma c/a. Richardson extrapolation estimates the continuum rate from second-order finite differences and is not a rigorous error interval. Box and resolution dependence are reported separately. The variational trial function in the report yields an exact negative Rayleigh quotient, independently supporting a negative eigenvalue.

```powershell
.venv\Scripts\python.exe scripts/scalar_tensor_completion.py
.venv\Scripts\python.exe -m unittest discover -s tests -v
```

The script derives/checks the local equations, saves the coefficients and eigenproblem results, and renders `01-action-and-ghost` and `02-coupled-instability`. No observations, Monte Carlo experiments, nonlinear dynamics or quantum corrections are included. All 66 repository tests pass in the recorded project environment. Source hashes normalize CRLF to LF; data hashes are byte exact.

AI assistance was used in implementation and analysis. Biswajit Jana remains the sole listed author; earlier scientific work retains attribution. Independent specialist review is pending.
