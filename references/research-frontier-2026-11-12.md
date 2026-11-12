# Research frontier and decision log

**Biswajit Jana | 12 November 2026**

I searched current primary literature to decide what to compute next. This is a targeted review, not an exhaustive priority search. It does not establish that a contribution has no predecessor. No observation of a physical wormhole is established by the papers or calculations below.

## Sources inspected and implications

| Source and verified version | Reading depth this round | Consequence for the project |
|---|---|---|
| Freivogel, Fumagalli & Tomašević, *How traversable is a traversable wormhole?*, [2606.12528v2](https://arxiv.org/html/2606.12528v2), 14 Sep 2026 | Version history, abstract, introduction, model/setup passages and discussion/future directions. Full matching derivation has not been independently reproduced. | MMP scalar transmission depends on frequency and collection time; charged massless fermions probe a different channel. A future extension must distinguish probe type, finite bandwidth, return time and backreaction. It is a separate physical model from Ellis. |
| Mondal & Sadhukhan, *Echoes of the Maldacena-Milekhin-Popov traversable wormhole*, [2511.22671v2](https://arxiv.org/html/2511.22671v2), 16 Sep 2026 | Version history and abstract; selected cavity-mode discussion and table inspected. Numerical integration not reproduced. | v2 corrects a shortcut from v1 and distinguishes mouth modes from extremely delayed cavity dynamics. Do not turn a resonance calculation into a detectable-echo forecast without dimensional timescales. |
| Lai, Lan, Wang & Piao, *Testing the wormhole echo hypothesis for GW231123*, [2602.01615v2](https://arxiv.org/abs/2602.01615v2), 1 Oct 2026 | Version history and revised abstract verified; full HTML obtained but not exhaustively checked. | The positive Bayes factor is waveform/prior sensitive and not unusual in targeted BBH-null ensembles. The authors conclude no robust wormhole-echo evidence. A favorable statistic alone is insufficient. |
| Clopper & Pearson (1934), [publisher record](https://doi.org/10.1093/biomet/26.4.404) | Bibliography verified; the standard interval construction is derived and tested in the local study. Original full paper not read this round. | Established binomial intervals plus the union bound provide a transparent conservative benchmark. They are not new statistical inventions. |
| Berger & Boos (1994), [publisher record](https://doi.org/10.1080/01621459.1994.10476836) | Abstract and bibliography, not full derivation. | Confidence-set treatment of nuisance parameters has prior literature; a later likelihood-based study should compare approaches without claiming the general idea as new. |
| [mpmath interval documentation](https://mpmath.org/doc/current/contexts.html), version 1.3.0 installed | Interval-context contract, supported elementary functions, endpoint behavior and experimental-status note | Recheck branch exclusions with interval arithmetic. Keep exact statistical coverage separate from numerical inverse-beta limitations. |

The first two sources were initially surfaced under older titles or versions; final scientific decisions use the revised primary records. Search snippets and secondary summaries were discovery aids only. In particular, an earlier favorable summary of GW231123 is not a substitute for the October revision.

## Search scope

Queries included wormhole inverse scattering/transmission distinguishability and calibration (2025-2026), wormhole mimickers and inverse echoes, current MMP probe transmission, GW231123 echo hypotheses, unknown classification-transition identifiability, Clopper-Pearson confidence construction and nuisance-parameter confidence sets. I checked arXiv version histories rather than treating crawl dates as publication dates. Some search results were irrelevant, unverified preprints or third-party mirrors and were not used as evidence.

## Work completed as a result

[Study 006](../reports/008-joint-calibration-identification.md) adds known-class calibration counts, a simultaneous-confidence exclusion argument, a continuum parameter search with interval-checked branch exclusions, 70,000 simulated experiments, and numerical bounds on the earlier matching allowance. These are candidate application-level contributions with explicit assumptions, not a proof of physical wormhole existence or a first-ever priority claim.

## Next questions, ranked by falsifiability

1. **Does calibration transfer?** Inject a prespecified drift between reference and science responses. Compare shared-response, bounded-drift and unmodeled-drift tests. Success means retained null control, not simply a larger discovery fraction.
2. **Can calibration cost be reduced?** Use a separate pilot to allocate science, false-negative and false-positive reference counts. Compare equal allocation against the derived local variance rule with a fixed total cost and held-out trials.
3. **Does phase break the intersection?** Derive the complex-amplitude experiment including phase-reference uncertainty. Search for joint magnitude/phase alternatives before simulating apparently decisive observations.
4. **Are proposed resonances resolvable?** For a reproduced MMP calculation, convolve with frequency response and impose finite observation time. Report dimensionless and physical timescales separately and test the probe-energy restrictions.
5. **What survives a physical alternative?** Replace the generic barrier with a fully specified comparator metric and perturbation equations. An astrophysical claim requires that extension before fitting real events.

These questions have explicit failure outcomes. A result that restores ambiguity or rules out a proposed measurement is useful research, even if it does not support a wormhole interpretation.
