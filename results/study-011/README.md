# Study 011 data and reproduction

**Author: Biswajit Jana · 23 November 2026**

[Report](../../reports/013-nuisance-confounding.md) · [Fit implementation](../../scripts/nuisance_identification.py) · [Runner](../../scripts/run_nuisance_study.py) · [Verifier](../../scripts/verify_nuisance_study.py)

All probabilities are noise-free model outputs. Fits are bounded numerical witnesses, not certified global minima, likelihood ratios, confidence intervals or astronomical detections.

| File | Meaning |
|---|---|
| `templates.csv` | 2,002 massless β = 0 ODE evaluations: `q` is the dimensionless wave number, `ell` the angular mode, `T` flux transmission, `flux_error` the conservation diagnostic. |
| `heldout.csv` | 242 independent midpoint ODE values, interpolated T and signed interpolation error. These finite checks are not a full-domain interpolation error bound. |
| `fits.csv` | 16 fit cases. `scale,mass,false_positive,false_negative` are s, ν₀, b, c. `interpolated_*` describe the optimization; `direct_*` use fresh ODE predictions. Maximum residual is absolute probability difference, not percent. |
| `optimizer_starts.csv` | Every deterministic start, success flag, evaluation count and sum of squared residuals. A success flag is not a global optimality certificate. |
| `predictions.csv` | 504 directly evaluated fitted channels. Includes the input truth, surrogate probability, raw direct T, response-adjusted direct probability and signed residual. `interpolation_difference` compares direct and interpolated recorded probabilities. |
| `parity.csv` | 66 spectra: two modes, three K values and 11 signed centered asymmetries. `lam` is λ; centered mass is fixed at 0.3 while `beta,mass` change. |
| `derivatives.csv` | 16 central-difference records, each containing four independently integrated spectra. The two step sizes are 0.001 and 0.0005. Stores raw T values and the derivative identity check. |
| `design.json` | Bounds, shared-response assumption, grid, environment, input byte hash and normalized source hashes. |
| `summary.json` | Counts, checked extrema and complete fit results. |
| `sha256.json` | Artifact hashes, including the three PNG/SVG figure pairs. |
| `verification.json` | Replay and reconstruction results; not itself part of the manifest. |

The parent input is `results/study-010/spectra.csv`, recorded by byte hash. Every fourth K is selected, giving 21 frequencies per mode. The truth has ν = 0.3, β = ±0.2, identity detector response and unit frequency scale. The null has β = 0 and p = b + (1−b−c) E_ell(sqrt((sK)²−ν₀²)). Bounds s ∈ [0.95,1.05], ν₀ ∈ [0.10,0.44], b,c ∈ [0,0.03] are assumptions for this numerical test, not measured calibrations. All bounds keep the null's sampled channels open and inside the template domain.

The objective is unweighted squared error. The response parameters are shared across all frequencies and modes, and b+c remains below one. No gain in observing precision or p-value can be inferred without specifying a noise and calibration model. Finite differences and convergence discrepancies are not statistical error bars.

```powershell
.venv\Scripts\python.exe scripts/run_nuisance_study.py
.venv\Scripts\python.exe scripts/verify_nuisance_study.py
.venv\Scripts\python.exe -m unittest discover -s tests -v
```

The runner performs 2,878 ODE integrations plus local fits, using four workers. The verifier replays 32 fitted predictions and 11 centered cases, checks all 504 predictions and 242 held-out discrepancies, reconstructs all finite differences, and refits both flexible joint cases. The complete fit and parameter records remain available even when plots show only selected summaries. Python source hashes normalize CRLF to LF; data hashes are byte exact. Floating-point and image output may differ by platform.

The three figure pairs are `01-nuisance-residuals`, `02-fitted-spectra` and `03-local-confounding`. Their plotted values come from these files. AI assistance was used for implementation and checks; Biswajit Jana is the sole listed author. Independent specialist review remains pending.
