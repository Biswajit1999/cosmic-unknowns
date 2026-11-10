# Study 002 data guide

**Author: Biswajit Jana**

All files here describe a **synthetic ideal counting experiment**. None contains observed wormholes, telescope measurements, detector strain or digitized literature data. The word “count” always means a pseudorandom binomial draw from an explicitly calculated probability. The physical generating model remains the prescribed Ellis probe problem.

Run `scripts/run_inference_study.py` from the repository environment to regenerate the dataset and four figures. The script writes `design.json` before generating counts. The settings were chosen for a controlled comparison, not registered with an external preregistration service. It uses direct ODE probabilities at 24 logarithmically spaced channels per design, with 48 generating samples in total.

| File | Contents and units |
|---|---|
| `design.json` | Model assumptions, dimensionless frequency ranges, count budgets, fit bounds, repetition counts and seeds. |
| `spectra.csv` | `design`, channel index, dimensionless `K`, exact-model numerical probability `T_ellis`, noiseless fitted `T_pt_projection`, Bernoulli KL divergence in natural-log units. |
| `example-counts.csv` | A separate, unselected example per design at N=10,000. Transmitted and reflected counts sum to N; `fraction` is transmitted/N; `wilson_low/high` are pointwise 95% score limits; fitted and generating probabilities are included. |
| `trials.csv` | Every calibration, independent null-check and Ellis-power experiment. N is incident counts **per channel**. `count_00`–`count_23` are transmitted counts, in spectrum channel order. Includes binomial deviance, fitted log C/log b, optimizer success, bound flag and maximum absolute unprojected gradient of the per-N negative log likelihood. |
| `rates.csv` | Rejection counts and fractions, conditional 95% Wilson Monte Carlo intervals, empirical rejection thresholds, bootstrap threshold sensitivities and induced rejection-rate sensitivities. Total incident count is 24N. |
| `numerical-audit.csv` | Every-third-channel ODE checks. Differences in transmission probability from changing domain and tolerance, and exact binomial standard deviation of the sample fraction at N=100,000. These differences are diagnostics, not statistical confidence intervals. |
| `pt-validation.csv` | Exact sech-squared formula versus an independent ODE calculation, covering both analytic branches and their junction. |
| `optimizer-audit.csv` | Additional-start checks for the first five repetitions of every case and all flagged fits. `improvement` is initial minus multistart deviance. No trials are discarded for poor fit. |
| `independent-optimizer-audit.csv`, `independent-audit.json` | Separate derivative-free refits with wider bounds, covering every original non-success and boundary fit plus a fixed ordinary-fit subset. Reports changes in deviance, thresholds and evaluation decisions. Run `scripts/audit_inference_study.py` after the main script to regenerate. |
| `summary.json` | Model projections, Fisher matrices per incident trial, diagnostics, environment versions and run time. |
| `sha256.json` | Hashes of saved CSV/JSON and scientific PNG/SVG outputs. Run time changes on regeneration, so the summary hash is not intended to be identical across runs. |

Two designs × four count budgets × (500 calibration + 500 null-check + 500 Ellis-power) repetitions give **12,000 experiments and 288,000 channel counts**. Each repeated experiment has a fresh comparator fit. Calibration is pointwise at each design's noiseless comparator projection; it is not a universal composite-null significance guarantee.

NumPy's default generator uses seed 20261110 for main counts, 20261111 for the plotted examples, and 20261112 for threshold bootstraps. The streams are separate. Raw counts are saved to support audits even if library changes affect regeneration. The script's ordered loops specify consumption of each stream. Versions are in `summary.json` and the repository lock file.

Wilson intervals on the power plot are conditional on a threshold estimated from 500 separate calibration experiments. The bootstrap columns report additional threshold sensitivity separately; they are not a combined confidence interval. Bounds and non-success flags must be considered when interpreting optimizer performance, especially in the weakly identified low band.

The conceptual hero lives outside this dataset and is not a scientific figure. See the [derivation](../../notes/003-counting-inference-derivation.md), [report](../../reports/005-transmission-identifiability.md), and [literature extension](../../references/study-002-literature.md).
