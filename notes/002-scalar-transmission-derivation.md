# Scalar transmission and instability on a prescribed wormhole

**Author: Biswajit Jana**
**Date: 10 November 2026**
**Classification: independent derivation and numerical reproduction of established model physics.**

I use the Ellis geometry as a tractable laboratory. I do not identify it with a detected astronomical object, or with the magnetically supported construction discussed elsewhere in this project. This note fixes the conventions behind my calculation; prior work is documented in the [literature audit](../references/literature-audit-2026-11-10.md).

## 1. Geometry, field and units

With signature $(-,+,+,+)$, throat radius $a>0$ and proper radial coordinate $x\in\mathbb R$, I write

$$ds^2=-c^2dt^2+dx^2+r(x)^2d\Omega^2,\qquad r(x)=\sqrt{x^2+a^2}.$$

I use the curvature convention in which

$$\mathcal R=-4\frac{r''}{r}+2\frac{1-r'^2}{r^2}
=-\frac{2a^2}{(x^2+a^2)^2}.$$

The probe field satisfies

$$[\Box-\mu_s^2-\xi\mathcal R]\Phi=0,\qquad \mu_s=mc/\hbar.$$

$\xi$ is the dimensionless curvature coupling; $m$ is the scalar particle mass. This is not the phantom field that can support the background. In four spacetime dimensions the usual massless conformal coupling is $\xi=1/6$; the special value $1/2$ appearing below is a property of this radial problem, not the conformal coupling.

Substitute $\Phi=e^{-i\omega t}Y_{\ell m}(\theta,\phi)u(x)/r(x)$. Since

$$\frac1{r^2}\partial_x\!\left[r^2\partial_x(u/r)\right]
=\frac1r\left[u''-\frac{r''}{r}u\right],$$

the radial equation becomes

$$u''+\left[\frac{\omega^2}{c^2}-\mu_s^2-
\frac{\ell(\ell+1)}{x^2+a^2}-\frac{(1-2\xi)a^2}{(x^2+a^2)^2}\right]u=0.$$

This agrees with the corrected curvature potential in [Taylor's erratum](https://doi.org/10.1103/PhysRevD.95.109904). The minimally coupled case can be benchmarked against [Azad, Loran and Mostafazadeh](https://doi.org/10.1140/epjc/s10052-020-08668-3).

Set $X=x/a$, $\Omega=\omega a/c$, $\mu=\mu_sa$ and $K=\sqrt{\Omega^2-\mu^2}$. The dimensionless operator is

$$H_\ell=-\frac{d^2}{dX^2}+U_\ell(X),\qquad
U_\ell=\frac{\ell(\ell+1)}{1+X^2}+\frac{1-2\xi}{(1+X^2)^2},$$

with $H_\ell u=K^2u$. Propagating channels require real $K>0$. Below the mass threshold I cannot use an incident plane-wave flux normalization. All computed transmission curves use $\ell=0$.

## 2. Flux, amplitude and boundary normalization

For a real potential the radial Wronskian current

$$j=\frac{u^*u'-uu'^*}{2i}$$

is conserved. With left incidence,

$$u\sim e^{iKX}+r_s e^{-iKX}\ (X\to-\infty),\qquad
u\sim t_s e^{iKX}\ (X\to+\infty),$$

the identical asymptotic wave numbers give $\mathcal R=|r_s|^2$, $\mathcal T=|t_s|^2$ and $\mathcal R+\mathcal T=1$. Here $\mathcal R$ in the flux equation denotes reflectance, whereas the curvature scalar above is contextually distinct. The code uses `R` for reflectance and never stores curvature under that name.

I integrate the right-outgoing solution inward from $X=L$. At $X=-L$,

$$A_{\rm in}=\tfrac12e^{iKL}[u+u'/(iK)],\qquad
A_{\rm out}=\tfrac12e^{-iKL}[u-u'/(iK)].$$

Then $t_s=1/A_{\rm in}$ and $r_s=A_{\rm out}/A_{\rm in}$. Finite $L$ truncates a nonzero potential tail. An apparently perfect flux balance therefore cannot, on its own, establish agreement with the infinite-domain problem. I vary $L$ and the integration tolerance independently.

## 3. An exact zero-energy check

For $\ell=0$, let $q=1-2\xi$ and set $X=\tan z$, $u=\sec z\,f(z)$. Direct differentiation gives

$$u''=(1+X^2)^{-3/2}(f''+f),\qquad f''+2\xi f=0$$

at $K=0$. For $\xi>0$, $s=\sqrt{2\xi}$ and a right-normalized zero-energy solution is

$$u_+(X)=\sqrt{1+X^2}\frac{\sin[s(\pi/2-\arctan X)]}{s}.$$

It approaches one at the right end. At the left end its growing coefficient is $B(\xi)=\sin(\pi\sqrt{2\xi})/\sqrt{2\xi}$. The limit at $\xi=0$ is $B=\pi$; for $\xi<0$ analytic continuation gives $B=\sinh(\pi\sqrt{-2\xi})/\sqrt{-2\xi}$. Away from zeros of $B$,

$$t_s\sim-\frac{2iK}{B},\qquad \mathcal T\sim\frac{4K^2}{B^2}.$$

In particular, minimal coupling gives $\mathcal T\sim4K^2/\pi^2$. This is a reproduction of a published low-energy benchmark, not an original prediction. At a zero of $B$, this expansion is invalid. It cannot be used to claim divergent transmission.

## 4. The transparency threshold and noncommuting limits

At $\xi=1/2$ and $\ell=0$, $U_0$ vanishes identically. Consequently $t_s=1$ and $\mathcal T=1$ for every $K>0$. This statement is exact for this probe equation. It does not describe electromagnetic waves, all angular modes, or a generic wormhole.

For any fixed $q>0$, by contrast, $\mathcal T\to0$ as $K\to0$. Hence

$$\lim_{q\to0^+}\lim_{K\to0^+}\mathcal T=0,\qquad
\lim_{K\to0^+}\lim_{q\to0^+}\mathcal T=1.$$

I test the crossover by shrinking $K$ and $q$ together. The weak barrier has area $\lambda=\int U_0dX=\pi q/2$. Replacing it by a delta barrier is justified only as the long-wavelength, weak-barrier asymptotic approximation, not as an identity at finite $q$. A delta barrier gives

$$\mathcal T_{\rm weak}=\frac{4K^2}{4K^2+\lambda^2}
=\frac{(K/q)^2}{(K/q)^2+(\pi/4)^2}.$$

My numerical sweep checks the approach to this curve, rather than assuming the approximation is valid at arbitrary frequency. I claim no novelty for this weak-scattering mechanism.

## 5. Why high transmission does not establish stability

For $\xi\leq1/2$, the massless scalar quadratic form is nonnegative: $\langle u,Hu\rangle=\int(|u'|^2+U_\ell|u|^2)dX\geq0$. There is no negative probe eigenvalue. At the critical spherical value there is a zero-energy resonance, not a normalizable growing mode; this statement is about exponential mode instability, not all possible long-time bounds.

For $\xi>1/2$ write $g=2\xi-1>0$, so $U_0=-gW$ with $W=(1+X^2)^{-2}$. The normalized trial function $u_\eta=\sqrt\eta\,e^{-\eta|X|}$ yields

$$\langle H\rangle=\eta^2-g\eta\int e^{-2\eta|X|}W(X)dX
=\eta^2-\frac{\pi g}{2}\eta+o(\eta).$$

For sufficiently small positive $\eta$ this is negative, proving a negative eigenvalue. If $E_0=-\gamma^2$, the massless field has a mode growing as $e^{\gamma ct/a}$. This reproduces the known onset at $\xi=1/2$.

I compare numerical growth rates with the second-order weak-attraction benchmark

$$\gamma\approx\frac{\pi g}{4}-\frac{3\pi g^2}{16}.$$

The integrals behind these two terms are $\int W=\pi/2$ and $\iint W(X)|X-Y|W(Y)dXdY=3\pi/4$. I retain only those displayed orders; I do not infer the size or analytic form of all higher-order remainders from a finite set of numerical points. The recent [scalarization preprint](https://arxiv.org/abs/2611.00541) provides a useful explicit comparator.

For a massive probe, $\Omega^2=E+\mu^2$. The lowest mode therefore grows when $\mu<\gamma$ and ceases to have negative squared frequency when $\mu\geq\gamma$. This simple spectral shift is exact within the prescribed ultrastatic linear model. It is not a mechanism stabilizing the wormhole's supporting matter.

The standard phantom-supported Ellis spacetime has a separate coupled metric–matter instability; see [González, Guzmán and Sarbach](https://arxiv.org/abs/0806.0608). Removing a negative eigenvalue from my additional probe equation does not remove that instability.

## 6. One-sided observations and an exact ambiguity

For each real frequency let $u_+$ be the solution outgoing on the right. Define $Z(K)=u_+'(0)/u_+(0)$ when the denominator is nonzero. I can remove the right half of the calculation and impose $u'(0)=Z(K)u(0)$ on the remaining half-line. Uniqueness of the linear ordinary differential equation guarantees the same left-side reflection amplitude, up to normalization.

Thus a one-sided measurement of this single linear channel determines an exterior response, not automatically the topology of its hidden completion. My numerical comparison checks that identity. I have not constructed a material surface realizing this load, established its engineering feasibility, or shown that a fixed frequency-independent Robin parameter would suffice. Additional channels, physical priors and independent observations could break narrower degeneracies.

## 7. Finite-bandwidth signals

For a massless complex scalar packet with spectral amplitude $A(K)$, integrated energy flux is proportional to $\int_0^\infty K^2|A(K)|^2dK$. I use Gaussian spectral power and compute

$$F_T=\frac{\int K^2|A(K)|^2\mathcal T(K)dK}{\int K^2|A(K)|^2dK}.$$

This is an all-time transmitted energy fraction. It is not a finite observation-time probability, group delay, detector strain, or photon count. In particular, no claim about realistic instrument detectability follows from this calculation.
