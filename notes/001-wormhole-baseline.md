# A checked starting calculation

**Author: Biswajit Jana**

**Classification: reproduction of a known theoretical geometry, not a new solution.** I use the Morris–Thorne form and the familiar zero-redshift choice $b(r)=r_0^2/r$, often associated with the Ellis wormhole geometry. The matter model is not specified here. [Geometric framework](https://doi.org/10.1119/1.15620)

## Geometry and conventions

Take signature $(-,+,+,+)$ and a static, spherically symmetric metric:

$$
ds^2=-e^{2\Phi(r)}c^2dt^2+\frac{dr^2}{1-b(r)/r}+r^2d\Omega^2.
$$

The areal radius is $r$, the throat radius is $r_0>0$, and $d\Omega^2$ is the metric of the unit two-sphere. Set $\Phi=0$ and $b=r_0^2/r$. At the throat, $b(r_0)=r_0$ and $b'(r_0)=-1<1$, meeting the static flare-out criterion.

The divergent radial coefficient at $r=r_0$ is removed by proper distance $l$:

$$
r(l)=\sqrt{l^2+r_0^2},\qquad
ds^2=-c^2dt^2+dl^2+(l^2+r_0^2)d\Omega^2.
$$

Here $l$ runs through zero from one side to the other. This coordinate check is not a dynamical stability proof.

## Required stress-energy

Let $\epsilon$ be local energy density in J/m³, $p_r$ radial pressure, and $p_t$ tangential pressure, both in Pa = J/m³. They are components in an orthonormal static frame. With zero cosmological constant, the static field equations give

$$
\epsilon=\frac{c^4}{8\pi G}\frac{b'}{r^2},\qquad
p_r=\frac{c^4}{8\pi G}\left[-\frac{b}{r^3}+
\frac{2}{r}\left(1-\frac{b}{r}\right)\Phi'\right].
$$

Write $C=c^4/(8\pi G)$. Substituting my chosen functions gives

$$
\epsilon=p_r=-\frac{Cr_0^2}{r^4},\qquad
p_t=+\frac{Cr_0^2}{r^4}.
$$

For a radial null vector with orthonormal components proportional to $(1,1,0,0)$, the null energy condition requires $\epsilon+p_r\geq0$. Instead,

$$
\epsilon+p_r=-\frac{2Cr_0^2}{r^4}<0,
\qquad
(\epsilon+p_r)_{r_0}=-\frac{c^4}{4\pi G r_0^2}.
$$

This is the obstacle exposed by the baseline. It is a requirement on the material source for this prescribed metric, not a supply of realizable negative energy. Ordinary negative gravitational potential energy is not the same quantity.

More generally, a smooth static throat in this metric family with regular redshift function satisfies

$$
(\epsilon+p_r)_{r_0}=\frac{c^4}{8\pi G r_0^2}[b'(r_0)-1]<0.
$$

This conclusion assumes Einstein gravity and the specified static throat conditions. It is not a theorem excluding every wormhole in every theory.

## Independent checks

The stress conservation equation for this anisotropic static source is

$$
\frac{dp_r}{dr}=-(\epsilon+p_r)\Phi'+\frac{2(p_t-p_r)}{r}.
$$

For my choice, both sides equal $4Cr_0^2/r^5$. The test suite compares a numerical derivative against this expression. Other checks confirm throat regularity, flare-out, and inverse-square scaling of the throat stress with radius.

For illustration, the $t=\mathrm{constant}$, equatorial spatial slice can be embedded in Euclidean three-space with

$$
r=r_0\cosh(z/r_0),\qquad l=r_0\sinh(z/r_0).
$$

The vertical coordinate $z$ is an embedding coordinate, not a measured extra spatial dimension. The picture shows only a two-dimensional slice of the spatial geometry. The associated plot uses $u=l/r_0$, for which $(\epsilon+p_r)/(C/r_0^2)=-2/(1+u^2)^2$.

## Interpretation limits

Quantum fields can exhibit negative renormalized energy, but that does not license an arbitrary static negative-energy profile. One must specify the field, state, renormalization prescription, and applicable averaged energy bounds, then solve for backreaction. The baseline has done none of those steps. See [Kontou's review](https://arxiv.org/abs/2405.05963).

No stability, natural formation mechanism, astronomical abundance, or observational identification follows from these equations. Those are research questions beyond this initial calculation.
