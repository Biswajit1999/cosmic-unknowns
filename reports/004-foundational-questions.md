# A mathematical agenda for the remaining foundational questions

**Author: Biswajit Jana**
10 November 2026 · Research agenda, not completed solutions

I retain the broader questions that motivated this project, while treating each as a distinct research program. The completed calculation is [Study 001](003-transmission-stability-study.md). The mathematical targets below specify what a subsequent result would need to establish; none is resolved by the scalar-transmission calculation.

## 1. Can physically supported wormholes exist?

The self-consistent target is a coupled solution, schematically

$$G_{\mu\nu}+\Lambda g_{\mu\nu}
=\frac{8\pi G}{c^4}\left(T^{\rm classical}_{\mu\nu}+\langle T_{\mu\nu}\rangle_{\rm ren}\right),$$

together with the matter equations, an admissible quantum state, appropriate energy bounds, regularity conditions and a perturbation analysis. Choosing a metric and solving backward for a stress tensor does not establish that the source can be realized. My first study does not solve this semiclassical system. Relevant starting points are the [four-dimensional magnetic construction](https://arxiv.org/abs/1807.04726v3) and the [quantum energy inequality review](https://arxiv.org/abs/2405.05963v2).

**Next answerable task:** reproduce one source-supported solution in its stated approximation and identify an independent stability criterion. **Rejection condition:** the candidate leaves its approximation regime or develops a growing physical mode before it can support the proposed signal experiment.

## 2. Does spacetime emerge from entanglement?

In static holographic settings, the leading Ryu–Takayanagi relation is

$$S(A)=\frac{\mathrm{Area}(\gamma_A)}{4G_N},$$

in units $\hbar=c=1$, where $\gamma_A$ is the appropriate bulk minimal surface anchored to the boundary region. This is a model-specific relation, with quantum corrections beyond leading order; it is not a proof that our cosmological spacetime emerges by the same mechanism. [Ryu and Takayanagi](https://arxiv.org/abs/hep-th/0603001)

**Next answerable task:** compute entanglement and reconstruction properties in an explicit finite tensor-network model, then change the network while holding selected entropies fixed. **Rejection condition:** a claimed unique geometry is not identifiable from the chosen information, or depends on an assumed geometric dictionary. A toy-model calculation would establish a property of the model, not the origin of spacetime in nature.

## 3. How is black-hole information recovered?

A useful semiclassical target is the generalized entropy expression

$$S(R)=\min\operatorname*{ext}_{I}\left[
\frac{\mathrm{Area}(\partial I)}{4G_N}+S_{\rm bulk}(R\cup I)\right],$$

again in units $\hbar=c=1$. It organizes entropy calculations with possible island regions $I$. Entropy accounting and operational recovery of a message are different questions. [Almheiri and collaborators](https://arxiv.org/abs/2006.06872)

**Next answerable task:** reproduce an explicit entropy transition in a specified model and separate the assumptions behind it from a proposed decoding procedure. **Rejection condition:** a supposed recovery argument provides only an entropy curve and no justified relation to recoverable information. I have not performed that calculation in this repository.

## 4. Was there a physical phase before the hot Big Bang?

I would compare explicit dynamical models through their primordial perturbation spectra. A broad class of single-effective-mode perturbation problems takes the form

$$v_k''+\left[c_s^2(\eta)k^2-\frac{z''(\eta)}{z(\eta)}\right]v_k=0.$$

Here primes denote conformal-time derivatives; $z$ and $c_s$ must be derived from the particular theory. This equation cannot be transferred unchanged into every bounce or quantum-gravity proposal. The observational target is a model-dependent prediction for primordial correlations, carried through the appropriate transfer functions and compared with data. [Planck's inflation analysis](https://arxiv.org/abs/1807.06211) is an observational starting point, not a universal test of every pre-Big-Bang scenario.

**Next answerable task:** specify one model with a controlled background and perturbation prescription, reproduce its spectrum, and compare it with an agreed baseline after fitting common nuisance parameters. **Rejection condition:** the purported distinctive feature can be reproduced by baseline freedom, or the perturbation calculation crosses an uncontrolled regime. No prior phase has been inferred here.

## 5. What is the global topology of space?

For a compact flat example with a spatial identification length $L_i$, allowed wave numbers can take the discrete form $k_i=2\pi n_i/L_i$. The predicted covariance of sky multipoles can then differ from the statistically isotropic form:

$$\langle a_{\ell m}a^*_{\ell' m'}\rangle
\neq C_\ell\delta_{\ell\ell'}\delta_{mm'}.$$

The covariance depends on the topology, orientation, scale and transfer physics. Not all topologies are flat tori, and a null search cannot exclude a topological scale beyond the observable region. [Planck's topology study](https://doi.org/10.1051/0004-6361/201525829)

**Next answerable task:** perform injection-and-recovery tests on synthetic skies before applying a specific statistic to public CMB maps. **Rejection condition:** a signal fails after masks, anisotropic noise, foreground residuals or the trial-factor penalty are included. No CMB maps have yet been analyzed here.

## 6. What replaces a classical singularity?

I would first distinguish bounded curvature from geodesic completeness. For example, finiteness of the Kretschmann scalar

$$K_{\rm curv}=R_{\alpha\beta\gamma\delta}R^{\alpha\beta\gamma\delta}$$

does not by itself establish that every relevant causal geodesic has a physically admissible continuation. A quantum proposal may require a different criterion involving unitary evolution, observables and the classical limit. [Discussion of singularity-resolution criteria](https://doi.org/10.1088/1361-6382/acb752)

**Next answerable task:** compare several explicit resolution criteria in one fully specified model. **Rejection condition:** a claimed resolution rests only on one finite invariant while evolution remains ill-defined. A complete quantum theory of gravity is not claimed.

## 7. Does physics forbid usable time machines?

An averaged null energy expression often appears in conditional restrictions:

$$\int_\gamma T_{\mu\nu}k^\mu k^\nu\,d\lambda\geq0.$$

Its meaning depends on which null curves, affine parametrization, spacetime and quantum assumptions are permitted. It is not an unconditional chronology theorem. I would need to track those hypotheses alongside the causal structure and backreaction. [Topological censorship and chronology protection review](https://arxiv.org/abs/0801.0735)

**Next answerable task:** audit one precise no-go theorem and one proposed exception, assumption by assumption. **Rejection condition:** the exception quietly violates an assumption that was retained in its claimed interpretation, or the proposed source fails before a closed causal path becomes usable.

## Work sequence

I will develop the wormhole program first because it now has a checked mathematical and numerical baseline. Each additional program needs a separate literature audit and a completed reproduction before a new claim. A publishable contribution may be a bound, a counterexample under explicit conditions, an identified degeneracy or a robust null result; it need not solve the entire foundational question.
