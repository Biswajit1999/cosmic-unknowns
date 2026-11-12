# Study 006: independent calibration and continuum exclusion

**Author: Biswajit Jana | 12 November 2026**

[Research report and proofs](../../reports/008-joint-calibration-identification.md).

All data are synthetic. Seventy designs (5 science budgets × 7 calibration budgets × 2 true error settings) each contain 500 Ellis and 500 comparator trials: 70,000 experiments. The raw counts are in 70 deterministic gzip CSV shards. Each row contains the group, replicate, solver status/node count, 24 science counts, 24 false-negative calibration counts and 24 false-positive calibration counts. M = 0 entries contain structural zeros, not calibration measurements.

## Files

- `design.json`: probability inputs, domain, seeds, software versions and source hashes.
- `rates.csv`: all 140 design/group rates, pointwise Wilson intervals, unresolved searches and incident-count costs.
- `population-brackets.csv`: full-domain exclusions and original witness allowances; see the verifier's slightly widened interval-checked feasible upper bounds.
- `exclusion-certificates.json`: complete branch partitions and the constraining channels used to exclude each branch.
- `verification.json`: seed-replay counts, solver replay count, 70-digit population checks and feasible witness checks.
- `summary.json`: experiment counts and aggregate diagnostics.
- `sha256.json`: immutable generated data/figure hashes. Verification output and prose documentation are not self-hashed.

## Interpretation

N is incident quanta per science channel. M is known-state calibration injections per class per channel: calibration costs 48M for 24 channels. Error rates describe classification, not simple detection inefficiency. The assumed calibration/science response equality is essential. The 1% generating error is not fixed to 1% when testing the null; it must be learned through intervals.

The mathematical test uses simultaneous Clopper-Pearson regions with per-tail probability 0.05/144. Its exact-arithmetic size bound is uniform over the stated parameter family for each experiment. The numerical implementation uses outward-padded SciPy beta quantiles and response bounds, followed by mpmath interval checks; inverse-beta evaluation is not formally validated. Any incomplete search is non-rejecting.

There are no empirical error bars on deterministic population exclusion bounds. Their limitations are numerical and model-based. Rate plots have Monte Carlo Wilson intervals; the heatmap's intervals are in `rates.csv`.

## Reproduce

```powershell
.\.venv\Scripts\python.exe scripts/run_joint_calibration_study.py --repetitions 500 --workers 4
.\.venv\Scripts\python.exe scripts/verify_joint_calibration_study.py
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
```

The new dependency is mpmath 1.3.0. The verifier replays all 5.04 million stored count entries, recomputes all aggregate rates/intervals, checks every population exclusion leaf at 70 digits and reruns a deterministic 1,400-decision sample. It does not rerun all 70,000 solver decisions a second time; each production rejection was interval-checked during generation.
