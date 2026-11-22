# Validation and observational distinguishability

**Biswajit Jana · Primary-source review checked 22 November 2026**

This update supports [Study 008](../reports/010-validation-restores-identification.md). It separates calibration practice from idealized wormhole observables; neither literature stream supplies observational data for the synthetic counting study.

## Calibration must retain time and frequency information

The current public [LVK O4 calibration entry, LIGO-T2500288-v6](https://dcc.ligo.org/LIGO-T2500288/public), includes separate uncertainty products for observing periods. Its [version-6 README](https://dcc.ligo.org/public/0202/T2500288/006/README.txt) describes time- and frequency-dependent complex response uncertainties and recommends matching the calibration product to the event time. Its seven-column LIGO format records frequency, median magnitude/phase error and 16th/84th-percentile boundaries. Those are not simultaneous 95% intervals and cannot be substituted for this project's binomial confidence sets. Reading depth: landing page and complete README; numerical archives were not analyzed. One O4b prose date is internally inconsistent with its interval; no calculation here uses it.

Wade et al., [arXiv:2207.00621v2](https://arxiv.org/abs/2207.00621v2), revised 10 January 2023, *Classical and Quantum Gravity* **40**, 035001, develop corrections for time-varying calibration filters. Section VI of the [full text](https://arxiv.org/html/2207.00621v2) studies source-inference consequences using simulated binary-neutron-star signals with deliberately imposed response offsets. The text distinguishes these synthetic errors from realized observing-run errors. Reading depth: abstract, introduction and selected Sections V–VI passages; filter implementation and parameter-inference runs were not reproduced.

The methodological implication for Study 008 is to retain response structure measured by validation, rather than collapse it to one magnitude. This is an inference motivating the benchmark, not a claim that the study reproduces an interferometer. Strain amplitude and phase errors are different quantities from Bernoulli classification errors.

## A different observable still needs source and geometry controls

Zhong, Rebouças, Otoniel and Guo, [arXiv:2610.27144v2](https://arxiv.org/abs/2610.27144v2), revised 2 November 2026, study optical appearance in asymmetric thin-shell wormhole models. Their [inclination discussion and conclusion](https://arxiv.org/html/2610.27144v2) distinguish geometric ray accessibility from brightness determined by emission and motion. Some sampled images approach their black-hole controls; other emission/viewing configurations reveal differences. The authors explicitly restrict these statements to their idealized models rather than establish observational distinguishability. Reading depth: abstract, selected inclination discussion and conclusion; ray tracing and junction calculations were not reproduced.

For future work, a second observable should be tested against physical alternatives with the same source and response assumptions. A different-looking ideal image alone cannot resolve observational identification. This paper uses a different geometry from the Ellis scattering model; its parameters and images have not been imported into the present simulation.

## Scope of this review

The review is targeted, not exhaustive. It establishes relevant precedents and concrete modeling requirements, not originality or priority. The [12 November audit](calibration-transfer-2026-11-12.md) remains applicable. No paper reviewed here establishes a physical wormhole detection or validates this project's hypothetical probe experiment.
