# Study 005: dense power experiment

**Author: Biswajit Jana**

This extension replaces the four-point presentation with actual simulations at 25 logarithmically spaced integer budgets between 100 and 100,000 incident quanta per channel. Both frequency designs retain 24 channels. Every design/budget has 1,000 calibration, 1,000 independent null-check and 1,000 independent Ellis-power trials: 150,000 experiments and 3.6 million channel counts in total.

Run `scripts/run_dense_power_study.py`. Four local worker processes handle independent cases. No AI subagents generate or select results. Case-specific seeds make output independent of scheduling. Each experiment refits the alternative; lines between plotted markers only guide the eye.

| File | Meaning |
|---|---|
| `design.json` | Count grid, model source hash, seeds, replication counts, fit bounds and interval definitions. |
| `low_band-N*.csv.gz`, `broad_band-N*.csv.gz` | Standard gzip-compressed CSV, 3,000 trials each. Columns include group, repetition, fitted deviance, log C/log b, optimizer status, boundary flag and transmitted counts in channels 00-23. N and design are in the filename and audit table. |
| `rates.csv` | Rejection counts, conditional Wilson intervals, thresholds and separately bootstrapped threshold sensitivities. |
| `optimizer-audit.csv` | Seeds, data-shard names, unresolved warnings, boundary fits, fallback counts and maximum changes in deviance. |
| `summary.json` | Experiment totals and fit audit summary. Any first tested budget exceeding a criterion describes only this grid and pointwise intervals, not a certified optimal budget. |
| `sha256.json` | Hashes of all shards, summaries and figures. |

The generator reads the 48 direct-ODE probabilities from study 002. It creates new counting trials; it does not relabel old trials or interpolate rejection rates. This is an ideal-detector simulation. The bounded-error ambiguity in study 004 remains relevant even when ideal-detector counts are large.

Intervals are pointwise and conditional on estimated thresholds. Bootstrap sensitivity is saved separately, not combined into the Wilson bars. False-positive calibration is pointwise at the chosen comparator projection and does not control the complete composite family uniformly. Inspect warning and boundary fields before using a run. Compressed CSV can be read with Python's `gzip.open(..., 'rt')` and `csv.DictReader` or ordinary data-analysis software.
