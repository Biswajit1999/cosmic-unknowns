# Literature audit: wormhole transmission, support and stability

**Author: Biswajit Jana**
**Search and inspection date: 10 November 2026**

I separate mathematical constructions, probe dynamics, coupled stability and observations. This is a focused review with a documented search trail, not a claim of exhaustive coverage of all wormhole publications.

## Evidence matrix

| Reference | What I use it for | Scope and reading depth |
|---|---|---|
| Morris & Thorne (1988), [doi:10.1119/1.15620](https://doi.org/10.1119/1.15620) | Static traversable-throat framework. | Publisher summary screened previously; my baseline equations are derived locally. Does not supply a demonstrated material source or formation history. |
| Friedman, Schleich & Witt (1993), [Topological Censorship](https://arxiv.org/abs/gr-qc/9305017), with [1995 erratum](https://doi.org/10.1103/PhysRevLett.75.1872) | Why no-go claims need explicit global and energy assumptions. | Abstract and erratum metadata screened. I do not deploy the theorem as a universal prohibition in this calculation. |
| González, Guzmán & Sarbach (2009), [Linear stability](https://arxiv.org/abs/0806.0608v2) | Separate the coupled phantom-supported background instability from an added scalar probe. | Abstract inspected. I rely on the published result; I have not reproduced their coupled perturbation operator. |
| González, Guzmán & Sarbach (2009), [Nonlinear evolution](https://arxiv.org/abs/0806.1370v2) | Background evolution need not preserve a stationary throat. | Abstract inspected. No nonlinear spacetime evolution is implemented here. |
| Taylor (2014; corrected), [arXiv:1404.7210v3](https://arxiv.org/html/1404.7210v3) | Arbitrary curvature coupling and scalar transmission/stability. | Sections III and version history inspected. Numerical WKB plots were not digitized. |
| Taylor (2017), [Erratum](https://doi.org/10.1103/PhysRevD.95.109904) | Corrected curvature and radial potential. | Relevant equations read directly in the three-page PDF. I use the corrected expression and derive the Ellis specialization independently. |
| Azad, Loran & Mostafazadeh (2020), [Transmission of low-energy scalar waves](https://doi.org/10.1140/epjc/s10052-020-08668-3) | Published analytical benchmark for minimal coupling. | Sections 2 and 5 inspected, especially Eqs. 13 and 51–56. I reproduce the low-energy coefficient; I do not claim to reproduce every result in the paper. |
| Gao, Jafferis & Wall (2017), [Double-trace construction](https://arxiv.org/abs/1608.05687v3) | Controlled example in a coupled-boundary holographic setting. | Abstract and version metadata inspected. Its geometry and boundary interaction are not those of my numerical model. |
| Maldacena, Milekhin & Popov, [Four-dimensional wormholes](https://arxiv.org/abs/1807.04726v3) | Magnetic, fermion-supported alternative to a prescribed phantom geometry. | Abstract and validity-related version notes inspected. Its stress-energy, magnetic charges and scales have not been implemented here. |
| Kontou (2024), [Quantum energy inequalities](https://arxiv.org/html/2405.05963v2) | Restrictions on quantum negative-energy support. | Selected energy-bound sections and discussion inspected. I do not apply a bound outside its field, state, averaging and spacetime assumptions. |
| Freivogel, Fumagalli & Tomašević (2026), [How traversable is a traversable wormhole?](https://arxiv.org/html/2606.12528v2) | Probe-dependent propagation in a magnetic wormhole and open questions about boundary conditions and omitted interactions. | Abstract and discussion inspected. The second mouth and the field content matter. My Ellis calculation is not a reproduction of their magnetic transmission spectrum. |
| Guvendi, Mustafa & Furtado (2026), [Exact scalarization thresholds and onset laws](https://arxiv.org/html/2611.00541v1) | Explicit recent comparator for threshold and second-order onset calculations. | Sections II–VI and IX inspected. Preprint; I independently check only the first spherical threshold and the displayed weak-coupling growth benchmark, not its entire threshold tower or nonlinear results. |
| Jafferis et al. (2022), [Quantum-processor experiment](https://doi.org/10.1038/s41586-022-05424-3), [author explanation](https://inqnet.caltech.edu/wormhole2022/) | Prevent conflation of quantum simulation with an astronomical tunnel. | Author explanation inspected previously. I make no claim that this experiment creates a macroscopic spacetime passage. |
| Kobrin et al. (2023), [Comment](https://arxiv.org/abs/2302.07897) | Critical counterpart to the processor interpretation. | Metadata screened only. Detailed assessment of the experiment and its replies is outside this study. |

## What changed after the review?

I rejected “discovering the critical scalar coupling” as a novelty claim: earlier corrected work already identifies the instability boundary, and a recent preprint develops it further. I also rejected “showing that wormhole signals can be suppressed” as new: low-energy scattering and probe dependence already have a literature.

My contribution in this repository is a transparent computational reproduction, an explicit approximation audit, and a single set of conventions connecting transmission, instability, mass thresholds and response ambiguity. I do not presently claim a new theorem or a first-in-literature result.

The immediate research hypotheses are more specific:

1. A proposed observational signature may lose discriminatory power after physically admissible boundary responses are allowed. My exact-load example establishes only a mathematical ambiguity; testing realistic alternatives remains open.
2. A controlled perturbation to a supporting field may change both propagation and background viability. A fixed-background calculation cannot answer that coupled question.
3. Near a zero-energy resonance, finite bandwidth can invalidate a naive low-frequency estimate. My double-scaling and packet calculations quantify this within the chosen model; a novelty claim would require a broader prior-art search.

## Search log and limitations

Queries included “Ellis wormhole scalar wave scattering transmission”, “nonminimal coupling scalar instability 1/2”, “Transmission of low-energy scalar waves”, “wormhole quantum energy inequalities”, and “topological censorship”. I followed the low-energy paper to its analytical formulas, the Taylor paper to its erratum, and the magnetic-wormhole discussion to its stated limitations. I also checked the recent scalarization preprint against the same operator.

This search is limited by indexing, accessible text and the selected model. It is not a complete citation census. Search-result dates were not used as publication dates; identifiers, version histories and publisher records were preferred. An abstract inspection is recorded as such rather than upgraded to full verification.

## Data provenance

The numerical CSV files are generated by my scripts, not extracted from telescope observations or third-party plots. Analytical curves are evaluated from the displayed formulas. No paper figures or full texts are redistributed. Model parameters are deliberately chosen numerical experiments, not fitted values of a measured wormhole.

This distinction is essential: there is no established observational wormhole dataset that could turn this calculation into an empirical detection. A future observational analysis must specify an actual instrument, source population, noise model, and alternatives before importing astronomical data.
