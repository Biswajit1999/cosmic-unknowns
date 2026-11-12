# Calibration transfer: primary-source audit

**Author: Biswajit Jana. Checked 12 November 2026.**

This focused update motivates [Study 007](../reports/009-calibration-transfer-drift.md) and supplements the [wormhole frontier review](research-frontier-2026-11-12.md). Reading depth is publisher/repository metadata and abstracts, not independent reproduction of derivations or data analyses. This is not an exhaustive review or a priority claim.

| Primary source and version | What the inspected abstract supports | Consequence here |
|---|---|---|
| Hido, Yamamoto, Chen, Sawada & Miyoki, *Physical Review D* **114**, 062010, 29 September 2026, [DOI 10.1103/vvwz-8b6t](https://journals.aps.org/prd/abstract/10.1103/vvwz-8b6t) | Calibration-line tracking and broadband reference measurements can disagree; a correction and its propagated uncertainty are evaluated. | A precise reference does not itself establish the science response. |
| Baka et al., [LIGO-T2500295-v6](https://dcc.ligo.org/LIGO-T2500295/public), *Correcting misspecification of calibration uncertainties in gravitational-wave data analysis with efficient reweighting* | A mismatch in calibration-uncertainty conventions was corrected with reweighting checked against rerun analyses. The abstract reports no significant change to prior conclusions. | Trace calibration conventions through inference code; this is not evidence that existing GW conclusions are generally unreliable. |
| Virgo Collaboration, [arXiv:2607.19273v1](https://arxiv.org/abs/2607.19273v1), 21 July 2026 | O4 strain reconstruction includes frequency-dependent uncertainties and offline updates; AnalysisReady data and uncertainties are released through GWOSC. | An observational extension needs the applicable response and uncertainty products. |

These papers concern interferometer strain response, including amplitude and phase. Study 007 concerns hypothetical known-class false-negative and false-positive probabilities. Their numerical uncertainties cannot be substituted for this project's drift bound. No observational data from these sources enter Study 007.

Binomial confidence intervals, simultaneous coverage, nuisance-set enlargement and overlap of observation models are established ideas. The application evaluated here is a reproducible benchmark and a narrowly verified population bracket for saved Ellis/comparator spectra. Establishing originality would require a broader systematic review and independent specialist assessment.
