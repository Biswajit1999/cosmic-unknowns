# Transmission inference with explicit counting uncertainty

**Author: Biswajit Jana · 10 November 2026**

I study an ideal controlled experiment, not an astronomical observation. My input is a fixed, massless, minimally coupled scalar spherical mode on the Ellis background. All lengths are measured relative to a reference throat radius; $X=x/a$ and $K=ka$. The generating potential is

$$
U_E(X)=\frac{1}{(1+X^2)^2},\qquad
[-\partial_X^2+U_E(X)]u=K^2u.
$$

I hold the throat scale, curvature coupling and mode fixed. This is a deliberately restricted identifiability calculation. I do not fit a population of unknown wormholes or infer a physical throat radius from the data.

## 1. A different conservative wave barrier

I use the positive sech-squared barrier

$$
U_P(X;H,\alpha)=H\operatorname{sech}^2(\alpha X),\qquad H,\alpha>0.
$$

It is a specified local wave operator. Positivity implies a nonnegative quadratic form for this isolated operator, but says nothing about a material realization or spacetime stability. I do not identify it with a black-hole metric or a complete alternative astrophysical object.

The exact scattering probability is

$$
T_P(K)=\frac{\sinh^2(bK)}{\sinh^2(bK)+C},\qquad b=\frac{\pi}{\alpha},
$$

$$
C=\begin{cases}
\cos^2\!\left[\frac\pi2\sqrt{1-4H/\alpha^2}\right],&4H\leq\alpha^2,\\
\cosh^2\!\left[\frac\pi2\sqrt{4H/\alpha^2-1}\right],&4H>\alpha^2.
\end{cases}
$$

This is the sech-squared result in [Boonserm & Visser, Section 3.8, Eq. 3.108](https://arxiv.org/pdf/1005.4483), translated to a radial operator with unit coefficient on $-\partial_X^2$. The two branches join at $C=1$. I test both branches and their junction against independent numerical integration.

## 2. A low-frequency degeneracy

For the minimal Ellis mode, the published and independently reproduced limit is

$$
T_E(K)\sim \frac{4K^2}{\pi^2}.
$$

See [Azad, Loran & Mostafazadeh](https://doi.org/10.1140/epjc/s10052-020-08668-3) and my [earlier derivation](002-scalar-transmission-derivation.md). Expanding the comparator gives

$$
T_P(K)=\frac{b^2}{C}K^2+
\left(\frac{b^4}{3C}-\frac{b^4}{C^2}\right)K^4+O(K^6).
$$

Therefore **every choice $C=\pi^2b^2/4$ matches the leading Ellis coefficient**. Matching this leading coefficient cannot determine either the barrier shape or the geometry. I do not assume that the full Ellis expansion is an even analytic Taylor series: its algebraic tail requires separate analysis beyond the leading term.

A useful exact relation for the comparator is

$$
\log\frac{T_i}{1-T_i}-\log\frac{T_j}{1-T_j}
=2\log\frac{\sinh(bK_i)}{\sinh(bK_j)}.
$$

The amplitude parameter $C$ cancels. Several frequency pairs can test a common value of $b$, whereas the small-$K$ limit approaches $2\log(K_i/K_j)$ independently of $b$. In practice I use the count likelihood rather than a noisy logit transform, because observed zero and full counts make empirical logits singular.

## 3. Ideal observation model

I define $N$ incident quanta in each of 24 monochromatic channels. Each outcome is classified perfectly as transmitted or reflected. Conditional on the model,

$$
Y_i\sim\operatorname{Binomial}(N,p_i),\qquad p_i=T(K_i),
$$

independently across channels. This uses scattering probabilities to define a thought experiment; it is not a proposed available wormhole apparatus. The source normalization, frequency calibration and detection efficiency are known, and there is no absorption, contamination, wave-packet mixing or background noise. Introducing these effects requires a new likelihood.

The exact count variance is $Np_i(1-p_i)$. For an observed fraction $\hat p=Y/N$, I display the pointwise 95% Wilson interval

$$
\frac{\hat p+z^2/(2N)}{1+z^2/N}
\ \pm\
\frac{z\sqrt{\hat p(1-\hat p)/N+z^2/(4N^2)}}{1+z^2/N},
\quad z=\Phi^{-1}(0.975).
$$

The expression follows inversion of the binomial score criterion in [Wilson (1927), pp. 209–210](https://jhanley.biostat.mcgill.ca/c607/ch08/wilson_jasa_1927.pdf). These are approximate-coverage frequentist intervals, not posterior probabilities or simultaneous coverage of all channels. They remain meaningful at zero counts, unlike a zero-width plug-in normal interval.

## 4. Likelihood, nuisance fitting and information

Let $\theta=(\log C,\log b)$ and

$$
\eta_i=2\log\sinh(e^{\theta_2}K_i)-\theta_1,\qquad p_i=(1+e^{-\eta_i})^{-1}.
$$

I minimize, up to data-dependent constants and the common factor $N$,

$$
L(\theta)=\sum_i\{\log(1+e^{\eta_i})-\hat p_i\eta_i\}.
$$

The stable implementation uses `logaddexp` and `expm1`. The analytic derivative is

$$
\nabla_\theta L=\sum_i(p_i-\hat p_i)
\begin{pmatrix}-1\\2z_i\coth z_i\end{pmatrix},\qquad z_i=e^{\theta_2}K_i.
$$

The Fisher matrix for $N$ trials per channel is

$$
I=N\sum_i p_i(1-p_i)v_iv_i^\top,\qquad
v_i=(-1,2z_i\coth z_i)^\top.
$$

At low frequency, $v_i\to(-1,2)^\top$ in every channel: the leading matrix has rank one. Only $2\log b-\log C$ is identified to that order. I save the finite-design Fisher matrices and their condition numbers in the chosen log-parameter coordinates; condition numbers depend on parameterization and are not coordinate-free discovery metrics.

I fit within $-16\leq\log C\leq25$ and $0.1\leq b\leq30$. These are explicit computational restrictions. A boundary solution is recorded, not silently dropped. Selected fits and all flagged fits receive additional starting points. This is an optimizer audit, not a proof of a unique global likelihood maximum.

## 5. Goodness of fit and repeated-experiment power

My test statistic is binomial deviance relative to the saturated channel model,

$$
D=2\sum_i\left[Y_i\log\frac{Y_i}{N\hat q_i}+
(N-Y_i)\log\frac{N-Y_i}{N(1-\hat q_i)}\right],
$$

where $\hat q_i$ comes from refitting the two comparator parameters to each experiment, and $0\log0=0$. I do not assign a universal chi-square reference law.

For each frequency design I first fit the noiseless Ellis probabilities to get a comparator projection $q^*$. For each $N$, 500 independent experiments from $q^*$ define an empirical 95th-percentile threshold (the `higher` sample quantile). A separate 500 comparator experiments assess false positives at this fixed null point. Another 500 Ellis experiments assess power. Each of those 1,500 experiments is refitted. This is a **pointwise simulation-calibrated test**: it does not control type-I error uniformly over the entire composite comparator family. For real data I would need an appropriate composite-null calibration and systematic-error model.

The displayed rejection-rate error bars are Wilson intervals over the 500 evaluation trials, conditional on the selected threshold. I separately bootstrap the 500 calibration deviances 1,000 times and save the resulting threshold and rate sensitivities. These sensitivity ranges are not combined confidence intervals and do not remove finite-calibration uncertainty.

For fixed probabilities, the expected log-likelihood ratio is

$$
\mathbb E_E\!\left[2\log\frac{\mathcal L_E}{\mathcal L_{q^*}}\right]
=2N\sum_i d_{\rm B}(p_i\Vert q_i^*),
$$

$$
d_{\rm B}(p\Vert q)=p\log(p/q)+(1-p)\log[(1-p)/(1-q)].
$$

The exact identity compares fixed distributions. Refitting changes the statistic's sampling distribution, so I do not substitute this expectation for a p-value or for measured power. The per-channel terms show where these selected models disagree.

## 6. Numerical and physical limits

The generating transmission is directly integrated at every channel; no interpolated surrogate supplies the count probabilities. On a fixed every-third-channel subset I compare domains 160, 320 and 640 and two tolerance settings. The 640 calculation also tightens tolerances, so its difference from the baseline is a combined diagnostic. I separately hold the domain at 320 for the tolerance check. These differences are convergence diagnostics, not certified error bounds or independent Gaussian errors.

Rejecting this comparator would reject its specified potential family under the observation assumptions. It would not uniquely identify a wormhole, prove a viable stress-energy source, determine a cosmic formation mechanism or repair the coupled instability of the phantom-supported Ellis solution. The [exact one-sided load equivalence](002-scalar-transmission-derivation.md) also survives this restricted comparison.
