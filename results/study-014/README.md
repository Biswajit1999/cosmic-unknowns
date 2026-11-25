# Study 014: finite magnetic scalar signals

**Author: Biswajit Jana · 25 November 2026**

These are deterministic calculations using a published low-frequency matching response, not observations or a full spacetime simulation. [Report and equations](../../reports/016-magnetic-finite-signals.md) · [Primary-source record](../../references/magnetic-signals-2026-11-24.md).

Units: ℏ = c = 1. x = ωr_e; z = ωℓ; b is the standard deviation of incident **radial energy** in z; h is the resonance half-width. Packet intensity has time standard deviation ℓ/(2b). The positive factor R_c = 1−p is not the complex reflected amplitude. Probability/energy ratios in the tables are dimensionless; none is an astronomical detection probability.

| File | Contents |
|---|---|
| `spectra.csv` | 3003 coarse monochromatic samples at three x values. This grid illustrates the comb and is deliberately **not** used for integrals or linewidth measurements. |
| `resolved_peak.csv` | 1001 samples in detuning/h at x = 0.02, resolving the plotted peak locally. |
| `bandwidth.csv` | 303 resonance-centered Gaussian averages, numerical error estimates, single-peak Lorentzian comparison, duration and normalized bandwidth. Lorentzian values are not valid over many resonances. |
| `bandwidth_audit.csv` | 45 independent transformed-integral/Fourier comparisons, including detuning 0.5 at the wrapped cell edge. |
| `matching_audit.csv` | 24 direct linear matching checks and flux residuals. |
| `collection.csv` | Complete-pulse accumulated incident energy and its ratio to the eventual transmitted energy; arrival times refer to the last pulse **center**. |
| `delay_scaling.csv` | 121 carrier frequencies: linewidth, first and eventual fractions, continuous and rounded 90% pulse counts, nominal collection time, resonant packet duration. |
| `pulse_profile.csv` | 10001 times for x = 0.05, b = 10, six amplitude terms; coherent energy flux. This illustrative response benchmark is not independently certified as a viable physical parameter set. |
| `time_audit.csv` | Three direct coherent time integrals versus finite double-sum Gaussian overlaps, with 40 pulses. |
| `compact_pulse_audit.csv` | Six compact cos-squared envelope cases, with durations 0.4, 1.4 and 8 in ell units and quadratic chirps 0 or 7; direct integration versus autocorrelation prediction for 12 output pulses. |
| `parameter_scan.csv` | 111 formal fixed-g flux values, joined to Study 013 s = 0.01 equilibrium; carrier quantum at 10% of reduced binding energy. Not a validity-certified physical region. |
| `scenario_audit.csv` | Four charge scenarios at b = 1000; geometry, energy and localization ratios; variable-frequency versus frozen first/eventual/collected energy. |
| `summary.json` | Settings, headline values, source hashes normalized to LF. |
| `verification.json` | Tighter-tolerance and independent saved-data checks; generated separately by the verifier. |
| `sha256.json` | Byte hashes for twelve CSV files, summary and six figure artifacts. The verification record is regenerated, not self-hashed. |

Three figures each have PNG and SVG versions: [bandwidth](../../figures/study-014/01-bandwidth.png), [collection](../../figures/study-014/02-collection-time.png), [experiment constraints](../../figures/study-014/03-experiment-constraints.png). All plotted curves are calculated from the saved values or explicit analytic comparison formulas in the plotting script.

The code freezes the frequency-dependent mouth coefficients only for the packet approximation, and audits this choice in the listed parameter scenarios. Off-resonance quadrature must resolve both wrapped Gaussian pieces; an internal quadrature error estimate is not sufficient validation. Narrow Gaussian energy spectra assume carrier frequency much larger than bandwidth. Gaussian time tails mean nominal pulse centers are not hard arrival boundaries or exact finite-width collection quantiles.

Run `scripts/run_magnetic_signals.py` then `scripts/verify_magnetic_signals.py`. The former depends on `magnetic_signals.py` and Study 013's `magnetic_energy.py`; existing Study 013 output files need not be regenerated. See the report for unresolved background, energy, ambient-return and quantum-state limitations. Independent specialist review remains pending.
