# Study 009 data and numerical uncertainty

**Author: Biswajit Jana · 22 November 2026**

[Scientific report](../../reports/011-conformal-scattering-nonuniqueness.md) · [Source](../../scripts/conformal_scattering.py) · [Run script](../../scripts/run_conformal_study.py) · [Verifier](../../scripts/verify_conformal_study.py)

All values are deterministic outputs for prescribed metrics. There are no astronomical measurements, random seeds, counting trials or sampling confidence intervals in this study. Production runs integrate the original field φ separately in each metric, rather than simply reusing the common reduced potential.

## Files and conventions

| File | Rows / contents | Columns and interpretation |
|---|---|---|
| `spectra.csv` | 729 integrations: 81 frequencies × 3 modes × 3 geometries | `k` is K = ωa/c; `ell` is the angular integer; `beta` is the conformal parameter; `xi` is curvature coupling (zero in this sweep). `r_real,r_imag,t_real,t_imag` give complex radiative amplitudes in the shared X convention. `R,T` are their squared moduli. `flux_error` = abs(R+T−1), `domain` is L, `order` is asymptotic tail order, `form=radial` means original φ equation, `evaluations` is the ODE function-call count. |
| `profiles.csv` | 3,603 rows: 1,201 X points × 3 geometries | `X` is dimensionless tortoise coordinate, `F` is lapse square root, `radius` is r/a, `ricci` is a² Ricci scalar. `rho,pr,pt,rho_plus_pr` are required Einstein stresses multiplied by 8πGa²/c⁴ (rho denotes energy density). `U0,U1,U2` are dimensionless minimally coupled radial potentials. `massive_delta_U` = 0.09(F²−1), a potential difference only. |
| `throats.csv` | 101 β points from −0.2 to 0.2 | `X,radius,lapse_sqrt,ricci` evaluated at the radius minimum. `left_lapse_sqrt` = 1−πβ. `radial_nec` is the scaled rho+pr at the throat. `flare_out` = d²(r/a)/dλ² at the throat, with dλ = F dX. |
| `convergence.csv` | 165 comparisons against 30 production-setting audit baselines | `check` identifies doubled domain, tail orders 8/16, tighter tolerance, reduced equation or plane boundaries. `delta_r,delta_t` are absolute complex-amplitude differences; `delta_transmission_probability` is abs(ΔT). Remaining fields record the alternative solution's boundary, equation, amplitudes and flux imbalance. The plane comparison calls the earlier solver, with the same finite L. |
| `design.json` | Design, environment and normalized source hashes | Records parameter grid, conventions, tolerances and tested software versions. Python sources are hashed after CRLF→LF normalization. |
| `summary.json` | Derived extrema and sample counts | Each numerical extremum is reconstructed by the verifier. |
| `sha256.json` | SHA-256 values for CSV, design, summary and figures | Data/figure byte hashes; nested `.gitattributes` prevents text-line-ending conversion for hashed artifacts. |
| `verification.json` | Independent audit output | Direct four-dimensional Christoffel/Ricci checks, 70-digit potential identities, replay counts, profile/throat checks and manifest status. This generated audit log is not itself in the data manifest. |

The amplitudes describe u = (r/a)φ with time normalized at right infinity. Fixed-areal-radius phases or finite-radius φ amplitudes require additional conversions; see the report. There are no units of physical seconds or metres without an independently chosen a. The illustrative massive term uses μa = 0.3; massive scattering is not computed.

## What numerical uncertainty means here

The maximum sampled geometry-to-geometry amplitude disagreement is 4.93 × 10⁻¹⁰. The larger 1.02 × 10⁻⁸ maximum from non-plane convergence comparisons is recorded separately. Neither is a statistical error bar or a rigorous universal upper bound. Audits vary boundary location, tail order, ODE tolerance and equation representation at five frequencies, three modes and two β values. They do not certify every frequency or parameter outside the sampled range.

The plane-boundary comparisons can differ by 0.19825 in complex amplitude while retaining good flux balance. Conservation is a necessary diagnostic, not a phase-accuracy certificate. Figure 4 displays both effects. At extremely small reflected/transmitted amplitudes, absolute agreement does not guarantee useful relative phase accuracy. Reflection phase is hidden where R ≤ 10⁻¹⁰.

## Reproduce

From the repository root, using the locked environment:

```powershell
.venv\Scripts\python.exe scripts/run_conformal_study.py
.venv\Scripts\python.exe scripts/verify_conformal_study.py
.venv\Scripts\python.exe -m unittest discover -s tests -v
```

The runner uses four worker processes, performs 729 production integrations plus 195 audit integrations (30 baselines and 165 alternatives), and creates four PNG/SVG figure pairs. The verifier replays 41 production integrations and 30 audit baselines; it checks the saved differences for all 165 alternatives but does not reintegrate every alternative. It independently constructs 12 four-dimensional tensors and checks 144 high-precision potential identities. Floating-point output and image rendering can vary across environments; regeneration refreshes the artifact manifest. The manifest detects changes, not correctness by itself.

The four figure pairs are `01-distinct-geometries`, `02-multimode-complex-scattering`, `03-noninvariant-observables`, and `04-numerical-audit` under [figures/study-009](../../figures/study-009). All plots use saved numerical columns and stated transformations. No generated decorative image is used as scientific evidence.

## Provenance

Biswajit Jana is the sole listed author. AI-assisted implementation and verification were used; primary-source context is credited in the report. These artifacts support a conditional fixed-background result and have not received independent specialist review.
