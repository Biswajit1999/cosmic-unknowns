# Literature extension: what a transmission measurement can establish

**Author: Biswajit Jana · Search date: 10 November 2026**

I extend the [earlier audit](literature-audit-2026-11-10.md) with sources relevant to a restricted inference experiment. This is a targeted primary-source review, not an exhaustive or systematic review. I separate papers that motivate a question from formulas that I actually inspect and validate.

| Source | Reading depth | Role in this study |
|---|---|---|
| Boonserm & Visser, *Quasi-normal frequencies: Key analytic results*, [arXiv:1005.4483v3](https://arxiv.org/pdf/1005.4483) | Section 3.8, Eqs. 3.105–3.110 inspected directly in the PDF. | Exact sech-squared transmission comparator. I translate its normalization and verify it against an independent ODE solver. I do not compute its quasinormal modes. |
| Wilson (1927), *Probable Inference, the Law of Succession, and Statistical Inference*, [primary article scan](https://jhanley.biostat.mcgill.ca/c607/ch08/wilson_jasa_1927.pdf), [doi:10.2307/2276774](https://doi.org/10.2307/2276774) | Score-interval derivation on pp. 209–210 inspected. | Binomial count and Monte Carlo rate intervals, with their frequentist interpretation. No claim of exact finite-sample coverage. |
| Azad, Loran & Mostafazadeh (2020), [doi:10.1140/epjc/s10052-020-08668-3](https://doi.org/10.1140/epjc/s10052-020-08668-3) | Relevant derivation inspected in study 001; publication and abstract rechecked in this search. | Minimal Ellis low-frequency coefficient. I combine it with the comparator expansion to exhibit a leading-order degeneracy; this is not a claim to have discovered inverse-scattering nonuniqueness. |
| Damour & Solodukhin (2007), *Wormholes as Black Hole Foils*, [arXiv:0704.2667](https://arxiv.org/abs/0704.2667) | Abstract and publication metadata inspected. | Historical motivation: distinct global geometries can share selected external signatures. I do not adopt its abstract's observational conclusions as universal theorems. |
| Lima Junior, Benone & Crispino (2022), *Scalar scattering by black holes and wormholes*, [arXiv:2211.09886](https://arxiv.org/abs/2211.09886) | Abstract and publication metadata inspected. | Partial-wave and angular-scattering observables offer a broader comparison than my single-mode transmission probability. Their interpolating geometry is not my generating model. |
| De Simone, De Falco & Capozziello (2025), *Can wormholes mirror the quasinormal mode spectrum of Schwarzschild black holes?*, [doi:10.1103/PhysRevD.111.064021](https://doi.org/10.1103/PhysRevD.111.064021) | Publisher abstract and metadata inspected; full text not used. | A recent example of matching a subset of modes through a restricted fit. The paper uses a near-throat parametrization and WKB estimates; my computation does not reproduce that result. |

The last three entries motivate caution about uniqueness; they do not validate my numerical values. The validation comes from the specified likelihood, analytic comparator, independent ODE checks and saved Monte Carlo records.

Searches included “Poschl Teller barrier transmission probability”, “wormhole scalar scattering inverse problem”, “inverse scattering transmission probability phase”, and “Wilson 1927 probable inference”. I followed the exact-barrier search to the original arXiv PDF and the interval search to a scan of Wilson's article. A phaseless inverse-scattering PDF returned an access error; I have not treated its inaccessible contents as reviewed evidence. No new observational catalog was imported, and no published figure was digitized.

## Research position and remaining gap

My question is narrower than proving physical wormholes: given two specified wave-operator families, what frequency coverage and ideal counting precision distinguish them? This produces a reproducible experimental-design result and exposes a low-frequency degeneracy. It is insufficient to distinguish all viable compact objects or to solve an unrestricted inverse problem. A stronger project would derive a competing object's wave operator from its matter and geometry, include multiple modes and phase information, then connect the forward model to an actual instrument and source.

The sech-squared barrier is already a standard exact scattering model. The binomial likelihood, Wilson intervals and Monte Carlo calibration are standard statistical constructions. The particular numerical comparison here is a working analysis; I make no priority claim without a broader citation and specialist review.
