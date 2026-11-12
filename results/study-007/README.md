# Study 007 data: calibration transfer drift

**Author: Biswajit Jana. All count data are synthetic.**

[Report](../../reports/009-calibration-transfer-drift.md) · [Literature audit](../../references/calibration-transfer-2026-11-12.md)

- `design.json`: spectra, comparator, optimization attempts, drift direction, seeds, domain and source hashes.
- `case-00.csv.gz` through `case-19.csv.gz`: 1,000 experiments each, with group, replicate, 24 science counts `y0..23`, 24 calibration false negatives `fn0..23`, 24 calibration false positives `fp0..23`, method statuses and node counts. M is per class per channel.
- `rates.csv`: 320 case/group/method rows with rejection totals, pointwise 95% Wilson endpoints, actual/assumed drift, validity flag, seed and filename. Eight methods share each experiment.
- `population-certificate.json`: complete exclusion partition for the lower spectral sup-norm distance bound.
- `population-boundary.csv`: 101 equation evaluations for the candidate overlap boundary, not extra trials.
- `summary.json`: aggregate counts and selected results.
- `sha256.json`: generated data/figure hashes; upstream hashes are in design.json.
- `verification.json`: full count replay, sampled decisions, 70-digit checks and hash outcomes. This verifier output is not self-hashed.

Only `excluded` counts as rejection. Other statuses are `compatible`, `undetermined` and `retained_by_nesting`. The latter conservatively inherits non-rejection for larger allowances. `matched_bound` receives the scenario's true maximum drift; it does not estimate drift.

Probability 0.003 means 0.3 percentage points, not 0.3% of the reference rate. Both reference error rates are 2%. Patterns are fixed before sampling and intentionally difficult, not typical detector models.

From the repository root with pinned requirements installed:

```powershell
.\.venv\Scripts\python.exe scripts/run_calibration_drift_study.py --repetitions 500 --workers 4
.\.venv\Scripts\python.exe scripts/verify_calibration_drift_study.py
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
```

PNG/SVG figures are in `figures/study-007`. The power/null plot has sampling error bars; heatmap intervals are in rates.csv. Population exchange and boundary figures are deterministic and have no sampling error bars. Model and numerical-spectrum uncertainties are separate from Monte Carlo error.
