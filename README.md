# Cosmic Unknowns

### What could distinguish a wormhole geometry—and what could support it?

**Author: Biswajit Jana · Version 0.7.0 · 22 November 2026**

I investigate what wave signals can tell us about a prescribed wormhole, which ambiguities survive perfect measurements, and whether the required geometry can be supported by a specified matter model. This is a public working research project with reproducible equations, numerical data and scientific figures. Independent specialist review remains pending.

![Calculated Ellis spatial embedding beside massless and massive scalar transmission spectra](figures/study-010/00-research-overview.png)

*The surface is a calculated equatorial spatial slice of the Ellis metric (β = 0), not a telescope view. The curves are generated model data. [Geometry and image derivation](reports/012-massive-probes-and-source-obstruction.md#6-scientific-image-what-the-readme-surface-represents) · [Data and figure provenance](results/study-010/README.md).*

**[Latest report](reports/012-massive-probes-and-source-obstruction.md)** · **[Data](results/study-010/README.md)** · **[Literature](references/massive-probes-2026-11-22.md)** · **[Research progress](RESEARCH_PROGRESS.md)** · **[Authorship](AUTHORS.md)**

## What the latest study establishes

| Question | Result | Scope |
|---|---|---|
| Can different geometries share all massless scalar modes? | Yes: the conformal family in Study 009 has an identical radial operator. | Fixed curvature coupling and stated phase/field normalization. |
| Does a known mass change the comparison? | Yes: the sampled massive spectra separate by up to **6.79 percentage points** from β = 0. | Fixed mass-scale product, frequency calibration and prescribed metric. |
| Can one static, minimally coupled radial scalar support the nonzero-β metrics? | No: the required stress violates its necessary identity ρ + p_t = 0. | The explicitly stated P(χ,Z) matter class; not a universal matter exclusion. |
| Have physical wormholes been demonstrated? | No. A viable source, coupled stability and observational evidence remain unresolved. | No claim of detection, construction, quantum-gravity solution or certified novelty. |

## A measurable difference within the model

The family is

$$
\frac{ds^2}{a^2}=F_\beta^2[-d\tau^2+dX^2+(1+X^2)d\Omega^2],
\qquad F_\beta=1+\beta(\arctan X-\pi/2).
$$

A minimally coupled scalar with inverse Compton length μ sees

$$
U_\ell=\frac{\ell(\ell+1)}{1+X^2}+\frac{1}{(1+X^2)^2}+(\mu a)^2F_\beta^2.
$$

The last term explains why a mass can reveal information absent from the massless response. I computed **972 spectra** across three geometries, two angular modes and two masses. At μa = 0.3, β = −0.2 and K = ωa/c = 0.67, the ℓ = 0 transmitted flux is **0.06794 lower** than the β = 0 reference. This is the largest sampled difference from that reference, not a global optimum or a detection significance.

![Massive and massless flux spectra with differences from the reference geometry](figures/study-010/01-massive-discrimination.png)

*Solid curves: massive scalar. Dashed curves: massless scalar. Lower panels show the geometry-dependent massive difference. No observation noise is simulated. [Full equations and limits](reports/012-massive-probes-and-source-obstruction.md) · [Raw spectra](results/study-010/spectra.csv).*

Massive scalar wormhole scattering is established in [Azad, Loran and Mostafazadeh (2020)](https://arxiv.org/abs/2010.15023). My extension evaluates this repository's conformal family, including unequal endpoint speeds, long-range phases and independent numerical checks. It does not claim that massive probes are a new idea.

## A signal is only part of the problem

For a single static radial scalar with a minimally coupled action P(χ,Z), the stress tensor necessarily obeys ρ + p_t = 0. The constructed metrics instead require

$$
\frac{8\pi Ga^2}{c^4}(\rho+p_t)=-\frac{2\beta X}{F_\beta^3(1+X^2)^2}.
$$

This excludes that sole-source model whenever β ≠ 0. Other actions require their own background equations and stability analysis. The massive field used to probe the metric is not assumed to support it.

![Channel-opening thresholds and the required stress that excludes a single static minimal scalar](figures/study-010/02-thresholds-and-matter.png)

*Left: analytical channel thresholds; production frequencies lie above both. Right: required stress combination, with zero demanded by the excluded scalar class. [Source derivation](reports/012-massive-probes-and-source-obstruction.md#5-a-restriction-on-supporting-matter).*

## How the calculation was checked

The solver includes logarithmic asymptotic phases and converts amplitude transmission to physical flux using T = (q_left/q_right)|t|². I tested boundary distance, expansion order, integration tolerance, the original field equation, reciprocity and known limiting cases.

| Verification | Recorded result |
|---|---:|
| Alternative numerical comparisons | 108 |
| Largest sampled complex-amplitude audit difference | 4.082 × 10⁻⁹ |
| Largest production flux imbalance | 1.248 × 10⁻⁹ or less |
| Independent four-dimensional tensor checks | 12 at 70-digit precision |
| Independent asymptotic coefficient checks | 12 Taylor expansions |
| Repository unit tests | 57 passed |

[See the numerical audit figure](figures/study-010/03-numerical-audit.png) and [verification record](results/study-010/verification.json). These are deterministic diagnostics, not confidence intervals or rigorous full-domain error bounds. Earlier counting studies retain their statistical error bars; this study does not invent error bars for exact model curves.

## Research library

| Research question | Report | Data / reading record |
|---|---|---|
| Massive probes and a supporting-matter restriction | [Study 010](reports/012-massive-probes-and-source-obstruction.md) | [Data](results/study-010/README.md), [literature](references/massive-probes-2026-11-22.md) |
| Different metrics, identical massless scattering | [Study 009](reports/011-conformal-scattering-nonuniqueness.md) | [Data](results/study-009/README.md), [literature](references/conformal-identification-2026-11-22.md) |
| Independent calibration and fixed measurement cost | [Study 008](reports/010-validation-restores-identification.md) | [Data and intervals](results/study-008/README.md) |
| Calibration-transfer drift and exact overlap | [Study 007](reports/009-calibration-transfer-drift.md) | [Data](results/study-007/README.md) |
| Joint statistical confidence construction | [Study 006](reports/008-joint-calibration-identification.md) | [Data](results/study-006/README.md) |
| Original formal identification report | [Manuscript](reports/007-calibration-limited-identification.md), [PDF](output/pdf/calibration-limited-identification.pdf) | [Originality audit](references/novelty-audit-2026-11-10.md) |
| An observer's hypothetical camera view | [Ray tracing and assumptions](reports/006-what-a-wormhole-could-look-like.md) | [Image provenance](figures/concept/README.md) |
| Scalar transmission and probe stability | [Original study](reports/003-transmission-stability-study.md) | [Data](results/study-001/README.md), [errata and literature](references/literature-audit-2026-11-10.md) |

## Reproduce the latest study

From the repository root on Windows:

```powershell
py -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements-lock.txt
.venv\Scripts\python.exe scripts/run_massive_study.py
.venv\Scripts\python.exe scripts/verify_massive_study.py
.venv\Scripts\python.exe -m unittest discover -s tests -v
```

On Linux/macOS, create the environment with `python3 -m venv .venv` and use `.venv/bin/python`. The tested environment is recorded in [design.json](results/study-010/design.json). The runner saves CSV data, PNG/SVG figures and hashes. [Detailed column definitions and numerical limitations](results/study-010/README.md).

## The next discovery tests

1. **Does the spectral difference survive unknown parameters?** Allow uncertain mass, scale, frequency calibration and detector response before claiming identification.
2. **Which explicit matter action survives the source test?** Verify every background equation, then analyze the coupled metric and matter perturbations.
3. **What could actually be measured?** Specify a physical source and detector for the proposed probe, or a calibrated clock/distance observable. A scalar proxy is not automatically a neutrino or photon experiment.

These are open tasks. Increasing data density alone does not resolve the physical assumptions. [Current progress and remaining checks](RESEARCH_PROGRESS.md).

All scientific figures are calculated from saved models; no astronomical wormhole dataset is claimed. I credit established results and disclose computational/AI assistance in [AUTHORS.md](AUTHORS.md) and the reports. Public availability does not imply peer review or journal submission.

<details>
<summary>Earlier study summaries and archived reproduction instructions</summary>

The following sections retain the earlier research record. The current results and next steps are described above.

## Study 009: geometric nonuniqueness

**[Read Study 009](reports/011-conformal-scattering-nonuniqueness.md)** · **[Data and numerical audit](results/study-009/README.md)** · **[Research progress](RESEARCH_PROGRESS.md)**

I construct a continuous family of distinct regular wormhole metrics with exactly the same massless-scalar radial potentials for every angular mode and fixed curvature coupling. The mechanism is an established harmonic conformal rescaling. Throat areas, curvature and static clock ratios differ, while complex scattering amplitudes coincide in the shared tortoise-coordinate and radiative-field convention. Independently calibrated areal distances and clocks supply additional information.

The study includes **729 original-field integrations**, **165 numerical comparisons**, a direct four-dimensional tensor audit and **four scientific figures**. Across the sampled geometries, complex amplitudes differ by at most **4.93 × 10⁻¹⁰**; separate convergence comparisons reach **1.02 × 10⁻⁸**. These deterministic differences are not statistical confidence limits. All metrics still require negative throat null-energy density; no viable matter source, stability or observed wormhole is established.

![Distinct geometries with common massless scattering](figures/study-009/01-distinct-geometries.png)

*Generated model data, not astronomical observations. [Complex scattering and phase](figures/study-009/02-multimode-complex-scattering.png) · [Clock and massive-probe discriminants](figures/study-009/03-noninvariant-observables.png) · [Boundary-error audit](figures/study-009/04-numerical-audit.png) · [Literature and reading depth](references/conformal-identification-2026-11-22.md).*

## Previous research: independent validation

**[Read Study 008](reports/010-validation-restores-identification.md)** · **[Data and uncertainty guide](results/study-008/README.md)** · **[Research progress](RESEARCH_PROGRESS.md)**

In 29,000 synthetic experiments, retaining the measured response in each channel can recover discrimination that is lost when validation is reduced to one maximum-drift bound. The study compares allocations at fixed incident cost and tests failures when validation does not transfer to science measurements. At one tested allocation, direct validation rejects the barrier family in 500/500 Ellis trials and 0/500 comparator trials. These are conditional model-discrimination results, not wormhole detections.

![Fixed-cost validation power and physical-null rejection](figures/study-008/01-fixed-budget-validation.png)

*Pointwise 95% Wilson intervals over 500 trials per point; all count data are synthetic. Methods share each trial. [Previous drift limit](reports/009-calibration-transfer-drift.md) · [Joint confidence construction](reports/008-joint-calibration-identification.md) · [Current literature context](references/validation-and-observability-2026-11-22.md).*

## Original formal scientific report

**[Read the PDF](output/pdf/calibration-limited-identification.pdf)** · **[Read the source manuscript](reports/007-calibration-limited-identification.md)** · **[Originality audit](references/novelty-audit-2026-11-10.md)**

**[Shareable research poster and LinkedIn caption](outreach/linkedin/README.md)** — a portrait PNG, editable SVG and first-person summary drawn from the archived data.

My original report develops a bounded detector-error criterion and constructs competing models with identical recorded probabilities. In the selected finite-channel benchmark, error allowances of about **0.01224% in the low band** and **3.0943% in the broad band** suffice for explicit matches. These are achievable numerical witnesses under a flexible classification-error model, not measured detector tolerances or certified global minima. They show why more counts alone cannot resolve an exact observation-model ambiguity.

The report contains an analytical proof, equations, likelihood histograms, error bars, calibration profiles, numerical audits and references. Known results are credited. I make no claim that the entire project is unprecedented or that physical wormhole existence is proved.

I investigate foundational questions through explicit mathematical models, literature review and reproducible computation. My first study asks how scalar waves propagate through a prescribed wormhole, when the added field becomes unstable, and what an observer can infer from one-sided reflection.

**[Read the research report](reports/003-transmission-stability-study.md)** · **[Read the derivation](notes/002-scalar-transmission-derivation.md)** · **[Inspect the literature audit](references/literature-audit-2026-11-10.md)**

## New: uncertainty and model distinguishability

I compare the Ellis transmission spectrum with a fitted positive sech-squared barrier. The dense extension uses **25 independently simulated budgets**, **150,000 synthetic experiments** and **3.6 million channel counts**. Each budget/design includes 1,000 calibration, 1,000 false-positive and 1,000 detection trials. The earlier four-budget study remains archived with its own data and report.

In the dense study, broad-band rejection rises from **28.0% at 1,000** incidents per channel to **92.1% at 3,162** and **98.2% at 4,217**. At 10,000, all 1,000 trials reject the comparator, with a conditional 95% Wilson interval of **99.62-100%**. These are ideal-detector model-discrimination results, not astronomical detections. Independent false-positive estimates and threshold sensitivity are reported alongside power.

**[Read study 002](reports/005-transmission-identifiability.md)** · **[Inference equations](notes/003-counting-inference-derivation.md)** · **[Data and uncertainty guide](results/study-002/README.md)** · **[Literature extension](references/study-002-literature.md)**

![Dense detection power with Monte Carlo error bars](figures/study-005/01-dense-power.png)

*Error bars are pointwise 95% Wilson intervals over repeated simulations, conditional on an independently estimated threshold. Additional threshold uncertainty is reported separately. All counting data are synthetic.*

**[Dense simulation data](results/study-005/README.md)** · **[Detector-error witnesses](results/study-004/README.md)**

![Calculated scalar transmission](figures/study-001/01-transmission.png)

## What I have calculated

**[What could a wormhole look like?](reports/006-what-a-wormhole-could-look-like.md)** — a calculated null-ray camera view, legends, finite-observer equations and synthetic-sky provenance.

- 525 spectral samples across five curvature couplings, with a published low-frequency benchmark.
- 176 additional samples testing a weak-barrier scaling approximation near transparency.
- Eight probe-instability growth rates, including a second-order analytical comparison, grid refinement and finite-box checks.
- An exact boundary-response equivalence tested at 32 frequencies.
- Finite-bandwidth energy fractions, with direct-quadrature checks on selected cases.
- Five new scientific figures in PNG and SVG, numerical CSV files, parameter records and artifact hashes.

The minimally coupled spherical mode transmits approximately **0.43%** of its flux at dimensionless wave number **ka = 0.1**, and **77.4%** at **ka = 1**. These are generated model results, not astronomical measurements. The main spectral sweep has a maximum flux imbalance of approximately **3.03 × 10⁻⁹**. Separate domain and resolution diagnostics are reported rather than inferred from that conservation check.

## What these results establish

This is a computational reproduction and analysis within a fixed Ellis geometry. I do not claim a new existence theorem, a stable material source, an observed wormhole, or a solution to the broader quantum-gravity questions. The phantom-supported background has a separate known instability that my additional probe-field calculation does not remove. Specialist review is pending.

## Project map

| Material | Link |
|---|---|
| Main study, results and limitations | [Research report](reports/003-transmission-stability-study.md) |
| Full radial derivation and analytical checks | [Mathematical note](notes/002-scalar-transmission-derivation.md) |
| Sources, errata, assumptions and reading depth | [Literature audit](references/literature-audit-2026-11-10.md) |
| Numerical files and column meanings | [Data guide](results/study-001/README.md) |
| Other foundational questions and mathematical targets | [Research agenda](reports/004-foundational-questions.md) |
| Original shortlist and initial proposal | [Landscape](reports/001-research-landscape.md), [brief](reports/002-wormhole-research-brief.md) |
| Geometry and required stress-energy | [Initial baseline](notes/001-wormhole-baseline.md) |
| Authorship and methodological provenance | [Author](AUTHORS.md), [standards](RESEARCH_STANDARDS.md) |

## Reproduce the work

I used Python 3.14.8; Python 3.11 or newer is the intended baseline. `requirements-lock.txt` records the tested environment, while `requirements.txt` gives minimum dependencies. From the repository root on Windows:

```powershell
py -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe -m unittest discover -s tests -v
.venv\Scripts\python.exe scripts/wormhole_baseline.py
.venv\Scripts\python.exe scripts/run_wormhole_study.py
.venv\Scripts\python.exe scripts/audit_wormhole_study.py
.venv\Scripts\python.exe scripts/run_inference_study.py
.venv\Scripts\python.exe scripts/audit_inference_study.py
.venv\Scripts\python.exe scripts/render_wormhole_camera.py
.venv\Scripts\python.exe scripts/run_calibration_study.py
.venv\Scripts\python.exe scripts/run_dense_power_study.py
.venv\Scripts\python.exe scripts/verify_research_artifacts.py
.venv\Scripts\python.exe scripts/verify_extended_studies.py
```

On Linux or macOS, use `python3 -m venv .venv` and `.venv/bin/python` in the subsequent commands. Study 001 is deterministic; study 002 uses the saved seeds 20261110–20261112 and saves every count. Neither imports astronomical observations. Floating-point details and execution time may vary by platform. The first main study took approximately 99 seconds and the second approximately 80 seconds in the recorded environment; supplementary audits add computation.

## Next research step

PDF build instructions and optional dependencies are in [the report build guide](output/pdf/README.md). All manuscript narrative lists Biswajit Jana as its sole author; scientific sources retain their original attribution.

The next physical step is to specify and verify supporting matter, then assess coupled stability. The next observability step is a calibrated clock/distance experiment or a known-mass probe with correct long-range scattering boundaries. Study 009 shows why more massless angular modes and counts alone cannot remove its exact operator ambiguity. Source and detector nuisance parameters remain part of any observational test.

The repository became public on 10 November 2026. Public availability does not imply peer review or journal submission. Code and manuscript licensing has not yet been selected; source citations preserve credit for prior work. The local working directory is `gravity-frontiers`, while the GitHub repository is [Biswajit1999/cosmic-unknowns](https://github.com/Biswajit1999/cosmic-unknowns).

</details>
