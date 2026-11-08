# Project W01: physical viability and signals of theoretical wormholes

**Status:** proposed research direction. Novelty unverified. No detection claim. No stability calculation or quantum matter construction has yet been performed.

## Central question

For a specified wormhole model, can its supporting matter, persistence, and propagation of signals be mutually consistent within the approximation used? If so, which predicted signatures could be distinguished from alternatives?

This is broader than our first calculation. We will narrow it to a single model and one controlled change after reproducing its published result.

## Literature anchors

1. **Morris and Thorne (1988):** a framework for describing static traversable geometries and their physical requirements. Use it to learn what the geometry requires; it is not a formation mechanism. [Paper](https://doi.org/10.1119/1.15620)
2. **Gao, Jafferis and Wall:** a construction using a coupling between two boundaries in an anti-de Sitter setting. It is not automatically a model of a naturally occurring astronomical tunnel. [Paper](https://arxiv.org/abs/1608.05687)
3. **Maldacena, Milekhin and Popov:** a four-dimensional magnetic construction supported by charged fermions, with a restricted regime of validity. [Paper](https://arxiv.org/abs/1807.04726)
4. **Kontou (2024):** quantum energy inequalities restrict the negative-energy distributions that quantum fields can supply. Their field, spacetime, and averaging assumptions matter. [Review](https://arxiv.org/abs/2405.05963)
5. **Freivogel, Fumagalli and Tomašević (2026):** scattering and transmission in the magnetic construction. The discussion explicitly identifies omitted interactions and boundary-condition questions. [Paper, version 2](https://arxiv.org/html/2606.12528v2)

The sources have been screened at the levels listed in the reading log. Full derivation-level verification remains a milestone.

## Candidate questions, to be checked for prior work

**A. Robustness to a specified boundary condition.** In a reproduced scattering model, how do finite observation time and a partially reflecting far boundary change resonances and recovered flux? Start from the paper's stated boundary condition. Introducing reflectivity may require a physically different setup; justify it rather than silently changing the problem.

**B. Robustness to a small mass or allowed interaction.** Does a claimed high-transmission channel survive a controlled perturbation of the field theory? A change to the same fields supporting the throat may also change the background. A calculation on a fixed metric must be labeled a test-field result and cannot establish self-consistent wormhole survival.

**C. Distinguishability from other compact objects.** In a specified simulated experiment, can a wormhole model be separated from a horizonless object with a reflecting surface after fitting nuisance parameters? An echo or resonance alone is not a unique wormhole signature. This question requires forward models for the alternatives, not merely a visually unusual plot.

**Recommendation:** reproduce the transmission baseline before selecting A, B, or C. These are candidate extensions; this scouting exercise has not established that they are unpublished.

## Mathematical route

A common radial wave reduction has the schematic form

$$
\frac{d^2\psi}{dx^2}+[\omega^2-V_{\mathrm{eff}}(x)]\psi=0.
$$

Here $x$ is a model-specific radial coordinate, $\omega$ a frequency, and $V_{\mathrm{eff}}$ a potential derived from the chosen geometry and field equations. This expression uses units $c=1$. It is a starting form, not a potential derived for the magnetic construction in this repository.

For a real, conservative scattering problem, with equal asymptotic wave numbers and compatible flux normalization, $R+T=1$. With unequal wave numbers, transmission includes a flux factor; with absorption or coupled channels, account for the missing or redistributed flux explicitly. Probability conservation and numerical convergence should be checked before interpreting resonances.

Plan a sequence of increasingly demanding checks:

- Recover a known analytical scattering limit and a published curve.
- Refine the grid, integration tolerances, domain size, and boundary placement.
- Vary one assumption at a time and retain the unmodified baseline.
- If observational inference becomes appropriate, specify a detector response, noise model, selection effects, and competing physical explanations.
- Ask an independent specialist to check the proposed novelty and physical interpretation before submission.

## What the present baseline establishes

The [worked note](../notes/001-wormhole-baseline.md) computes a known static geometry and its required stress-energy. The script checks a local conservation identity and geometric regularity in proper radial distance. It does **not** establish that quantum fields realize this stress-energy, that the geometry is stable, that it forms naturally, or that it obeys a relevant quantum energy inequality.

The zero-redshift toy geometry is also **not** the magnetic wormhole used in the 2026 transmission paper. Keeping those models distinct prevents an illustrative calculation from being mistaken for a reproduction of that result.

## Deliverables for the next research cycle

1. An assumption-by-assumption literature matrix and a record of later citing papers.
2. A reproduced baseline with dimension and convergence checks.
3. One written conjecture with an explicit way to reject it.
4. A numerical or analytical result, including failures and null results.
5. A report suitable for specialist criticism, followed by a public repository if desired.

The current repository supplies the initial map and an educational baseline. It does not yet meet these later research milestones.
