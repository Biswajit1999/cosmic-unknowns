# Study 010: data, figures and numerical uncertainty

**Author: Biswajit Jana · 22 November 2026**

[Report](../../reports/012-massive-probes-and-source-obstruction.md) · [Solver](../../scripts/massive_conformal.py) · [Runner](../../scripts/run_massive_study.py) · [Verifier](../../scripts/verify_massive_study.py)

These are deterministic calculations on prescribed geometries, not astronomical measurements. The mass parameter is ν = μa, with μ the inverse Compton length. No particle species, actual throat scale, instrument noise or empirical detection significance is assigned.

## Artifact dictionary

| File | Contents and units |
|---|---|
| `spectra.csv` | 972 rows: 81 frequencies × 3 β × 2 masses × 2 modes. `k` = ωa/c, `mass` = μa, `ell` = angular index. `q_right,q_left` are asymptotic dimensionless wave numbers. `r_*`, `t_*`, `t_flux_*` store real/imaginary amplitudes, the last with the velocity normalization. `R,T` are reflected/transmitted flux fractions. `flux_error` = abs(R+T−1); `tail_current_error` is the largest absolute relative outgoing-current error at the boundaries. `domain` = L, `order` = asymptotic series order, `form` identifies the equation, `incidence` the incident end, `evaluations` the solver function count. |
| `convergence.csv` | 108 alternatives against 18 audit baselines. `check` identifies the change. `delta_t_flux`, `delta_r` and `delta_transmission_probability` compare with the baseline. `delta_r` is blank for reverse incidence because reciprocity does not require identical reflection phases. Other columns describe and store the alternative solution. |
| `source_profiles.csv` | 1,803 rows: 601 positions per β. `rho_plus_pt` comes from the required Einstein stress; `analytic_obstruction` independently evaluates −2βX/[F³(1+X²)²]. Both are in units multiplied by 8πGa²/c⁴. |
| `thresholds.csv` | 121 β points; `right_threshold` = ν and `left_threshold` = ν(1−πβ) at ν = 0.3. Analytical channel-opening curves, not near-threshold numerical scattering. |
| `embedding_meridian.csv` | 161 X points, `radius` = sqrt(1+X²), `embedding_z` = asinh(X), all in units of a. Rotated through 91 azimuths for the β = 0 equatorial spatial illustration. The auxiliary embedding height is not a physical extra dimension. |
| `design.json` | Parameter grid, tolerances, phase convention, environment, embedding settings and source hashes normalized from CRLF to LF. |
| `summary.json` | Sampled maxima and counts. `delta_T` inside the maximum-gap record is an absolute probability difference, not a relative percentage. |
| `sha256.json` | Byte hashes for data, design, summary and figure pairs. Scoped `.gitattributes` protects their line endings. |
| `verification.json` | Replay results, independent tensor and Taylor-series checks, limiting-case diagnostics and manifest status. The generated audit output is not itself included in the data manifest. |

All production spectra have right incidence. Thus T = (q_left/q_right)|t|². Reverse-incidence audit rows use the inverse ratio; T = |t_flux|² in either case. The phase basis at each outward radius y is exp[i(qy−η ln(2qy))] times a finite inverse-radius series. This specifies constant-phase conventions and does not describe an independently calibrated source/detector experiment.

## Numerical versus statistical uncertainty

The largest sampled massive difference from β = 0 is 0.0679397421 in T at β = −0.2, K = 0.67, ℓ = 0, μa = 0.3. The maximum complex-amplitude audit difference is 4.082 × 10⁻⁹; production flux imbalance is at most 1.248 × 10⁻⁹. Neither is a confidence interval or a rigorous bound over all possible parameters. The physical difference and the numerical diagnostics are different quantities.

The six audit settings are double boundary radius, tail orders 8 and 16, tolerances reduced by a factor of 100, original-field integration, and reverse incidence. There are 18 audit problems. Pure massless production rows use a β-independent reduced potential, so their exact numerical overlap is not an independent verification of Study 009's analytic identity. That earlier study contains original-field checks.

No statistical error bars are invented for deterministic curves. The audit figure plots discrepancies, with values below 10⁻¹⁶ placed at that display floor. Figure generation uses only saved arrays and their stated transformations. A finite-grid maximum is not a continuous optimum, and no inverse problem with unknown scale, mass or detector response has yet been solved.

## Reproduce

From the repository root, after installing `requirements-lock.txt`:

```powershell
.venv\Scripts\python.exe scripts/run_massive_study.py
.venv\Scripts\python.exe scripts/verify_massive_study.py
.venv\Scripts\python.exe -m unittest discover -s tests -v
```

The runner uses four processes and computes 972 production integrations plus 126 audit integrations (18 baselines and 108 alternatives). The verifier checks all recorded production flux quantities and audit differences, replays 36 production runs and 18 audit baselines, compares 12 ultrastatic limits, constructs 12 four-dimensional tensors at 70 digits, checks 12 Taylor coefficient arrays, and verifies geometry/threshold files. It does not reintegrate every alternative. Hashes detect artifact changes, not scientific correctness by themselves. Numerical and rendering details may vary with platform.

## Figure provenance

All assets in [figures/study-010](../../figures/study-010) are generated with Matplotlib and saved in PNG and SVG:

- `00-research-overview`: β = 0 spatial embedding and ℓ = 0 flux spectra, with equations and interpretation.
- `01-massive-discrimination`: both modes, both masses, and probability differences from β = 0.
- `02-thresholds-and-matter`: channel-opening curves and the single-static-scalar source obstruction.
- `03-numerical-audit`: convergence diagnostics and the need for velocity-weighted transmission.

The sole listed author is Biswajit Jana. AI-assisted implementation and verification are disclosed in the report. Literature sources remain credited, and independent specialist review remains pending.
