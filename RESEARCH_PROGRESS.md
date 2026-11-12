# Research progress: 12 November 2026

**Author: Biswajit Jana**

Latest completed work: [Study 006: Independent calibration and continuum exclusion](reports/008-joint-calibration-identification.md).

- Defined a joint science/calibration experiment and proved the exact-arithmetic false-positive bound for its confidence-region exclusion rule.
- Implemented interval-checked exclusion over the complete stated comparator domain, with non-rejecting computational fallback.
- Completed 70,000 synthetic experiments and saved all count records, uncertainty estimates and replay checks.
- Narrowed the population matching allowance for the saved finite-channel spectra using verified exclusion branches and a feasible upper witness.
- Reviewed current primary versions of recent MMP transmission/echo and GW231123 papers; documented what has and has not been read or reproduced.

The solved subproblem is conditional statistical identification, not physical wormhole existence. The numerical beta quantiles remain outside the validated-arithmetic certificate. Calibration transfer, physical alternatives and independent review remain essential.

## Next study specification

**Primary question:** how much calibration-to-science response drift can be tolerated while retaining valid null control and useful power?

Prespecify a bound epsilon on each science/reference rate difference. Enlarge reference confidence intervals by epsilon, clip to [0,1], and reuse the continuum exclusion. Simulate both bounded drift and deliberate bound violations. The former must retain the confidence-set argument; the latter should reveal failure rather than be hidden. Compare the science/calibration allocation at fixed incident-count cost, with a separate pilot if optimizing allocations. Do not call a future simulated outcome achieved before computing it.

**Parallel scientific direction for later work:** reproduce finite-time MMP scattering from the current paper versions, then convolve frequency response and impose dimensional observation-time and probe-energy limits. Do not reuse Ellis probabilities for that geometry.

[Data guide](results/study-006/README.md) · [Current literature review](references/research-frontier-2026-11-12.md) · [Original formal report](reports/007-calibration-limited-identification.md)
