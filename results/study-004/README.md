# Study 004: bounded detector-error analysis

**Author: Biswajit Jana**

This dataset uses theoretical probabilities and the earlier synthetic counts. It contains no astronomical observations or measured detector calibration. Run `scripts/run_calibration_study.py`.

- `witnesses.csv`: every physical comparator probability, target probability, constructed false-positive/false-negative rate and resulting recorded probability; includes a refined Ellis probability for each channel.
- `kl-profiles.csv`: 47 allowed-error budgets per design; numerically profiled summed Bernoulli KL and the fixed-distribution expected twice-log likelihood ratio at N=10,000. These are not p-values or power estimates.
- `histogram-bins.csv`: fixed bins 0, 5, ..., 200 for initial-study deviance distributions; 500 trials per group and pointwise Wilson frequency intervals. Zero-count bin lower limits are exactly zero. Bin frequencies are correlated, and bars are not simultaneous bands.
- `summary.json`: source hashes, feasible witness parameters, all constrained-fit starts, an independent differential-evolution check, numerical-refinement sensitivity and claim limitations.
- `sha256.json`: data and figure hashes.

Each channel has separate false-positive and false-negative nuisance parameters bounded by delta. Recorded transmission is `(1-fn)*q + fp*(1-q)`. These errors describe classification flips, not merely missed detections. The per-channel freedom is an assumption; no smoothness, common response or external calibration data are imposed.

The located budgets are constructive upper bounds on the smallest matching allowance, not certified global optima. At the saved baseline probabilities, the algebraic construction matches exactly to stored precision. Recomputing the Ellis probability changes the matching allowance, especially at small probabilities; the refined values are recorded rather than hidden. The [formal report](../../reports/007-calibration-limited-identification.md) contains the proof and limitations.
