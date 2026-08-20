---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.1
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# 9. Characterizations of Mean Square Differentiability and Mean Square Continuity

```{eval-rst}
.. index::
   single: initial value theorem
   single: Wold kernel
   single: mean square differentiability; kernel criterion for
```

{doc}`02_mean_square_continuity_differentiability` characterized mean square continuity and
differentiability through the behavior of the autocovariance $R(\tau)$ near $\tau = 0$.

We now restate those conditions directly in terms of the Wold moving average kernel $p(\tau)$
of {doc}`08_spectral_densities`. We shall apply that test repeatedly to the linear stochastic
differential equations of {doc}`11_linear_sde`.



(thm-msd-representation)=
```{eval-rst}
.. index::
   single: mean square differentiability; kernel criterion
   single: Wold kernel; behaviour at the origin
   single: absolute continuity; of the kernel
```

**Theorem 11.** Let $x(t)$ be a linearly indeterministic covariance stationary process with moving average representation

$$
x(t) = \int_{0}^{\infty} p(\tau) w(t-\tau) \, d\tau
$$

or

$$
x(t) = \tilde P(D) w(t).
$$

Then the derivative $Dx(t)$, understood for the moment as a generalized process, has
representation

```{math}
:label: eq-9-1
Dx(t) = \int_{0}^{\infty} p'(\tau) w(t-\tau) \, d\tau + p(0)w(t)
```

or, in operator form,

$$
Dx(t) = D\tilde P(D)\, w(t).
$$

The impulse $p(0)\,\delta(\tau)$ that {eq}`eq-9-1` carries is what makes the qualification
necessary: when $p(0) \neq 0$ the right side is a generalized process, and the mean square
derivative does not exist. The conditions below say exactly when the impulse is absent.

The Laplace symbol of $D\tilde P(D)$ is $s\tilde P(s) = [\,s\tilde P(s) - p(0)\,] + p(0)$, in which
the bracketed term is the Laplace transform of $p'$ and the constant $p(0)$ is the transform of
the impulse $p(0)\,\delta(\tau)$; this is exactly the split displayed in {eq}`eq-9-1`.

Assume the kernel $p$ is absolutely continuous, so that its derivative $Dp = p'$ exists almost everywhere and $p(\tau) = p(0) + \int_0^\tau p'(u)\, du$. Then, by {eq}`eq-9-1`, the mean square derivative exists as an ordinary stochastic process if and only if (a) $p(0) = 0$, and (b) $\int_{0}^{\infty} |Dp(s)|^2 \, ds < +\infty$.

**Proof.** Write the moving average as  

$$
x(t) = \int_{-\infty}^{t} p(t-s)\, w(s) \, ds, \qquad p(u) = 0 \ \text{ for } u < 0 .
$$

Here $t$ appears both in the integrand and in the upper limit of integration, so differentiating
with respect to $t$ by Leibniz's rule produces one term from each:

$$
Dx(t) = \underbrace{p(0)\, w(t)}_{\text{upper limit}} + \underbrace{\int_{-\infty}^{t} p'(t-s)\, w(s)\, ds}_{\text{integrand}} = p(0)\, w(t) + \int_{0}^{\infty} p'(\tau)\, w(t-\tau)\, d\tau ,
$$

which is {eq}`eq-9-1`. The boundary term $p(0)w(t)$ is a white noise; since white noise is a
generalized process with infinite instantaneous variance ({doc}`04_physical_realizability`),
$Dx(t)$ is an ordinary, finite-variance process only if $p(0) = 0$. When $p(0) = 0$, the
remaining moving average has variance $\int_{0}^{\infty} p'(\tau)^2\, d\tau$, so it is well
defined as a mean square process precisely when $\int_{0}^{\infty} |p'(\tau)|^2\, d\tau < \infty$,
which are conditions (a) and (b). Conversely, when (a) and (b) hold, the absolute continuity of $p$
makes the difference quotients $[\,p(\cdot+\epsilon) - p(\cdot)\,]/\epsilon$ converge in $L^2$ to
$p'$ (extend $p = 0$ on $(-\infty, 0)$; this extension is unbroken at the origin because
$p(0) = 0$). Hence $[\,x(t+\epsilon) - x(t)\,]/\epsilon$ converges in mean square to
$\int_0^\infty p'(\tau)\, w(t-\tau)\, d\tau$, and the mean square derivative exists. This
establishes the equivalence.

The necessity of $p(0) = 0$ can also be seen directly from the autocovariance. With
$E\, w(s) w(s') = \delta(s-s')$, the moving average gives, for $\tau \geq 0$,

$$
R(\tau) = E\, x(t)\, x(t-\tau) = \int_{0}^{\infty} p(b)\, p(b+\tau)\, db .
$$

Differentiating and letting $\tau \downarrow 0$,

$$
R'(0^+) = \int_{0}^{\infty} p(b)\, p'(b)\, db = \tfrac{1}{2}\big[\, p(b)^2 \,\big]_{0}^{\infty} = -\tfrac{1}{2}\, p(0)^2 .
$$

Because $R$ is even, $R'(0^-) = +\tfrac{1}{2}\, p(0)^2$, so $R'$ has a jump of size $-p(0)^2$ at
the origin unless $p(0) = 0$. Mean square differentiability requires $R$ to be twice
differentiable at $\tau = 0$, hence $R'$ continuous there; this forces $p(0) = 0$, in agreement
with condition (a).



A generalization of Theorem 11 to higher order derivatives follows by applying it repeatedly to the successive derivatives of {eq}`eq-9-1`.

(thm-msd-higher-order)=
```{eval-rst}
.. index::
   single: mean square differentiability; higher order
   single: kernel; vanishing derivatives at the origin
```

**Theorem 12.** Assume $p$ is $n-1$ times differentiable with $D^{n-1} p$ absolutely continuous, so that $D^n p$ exists almost everywhere (for $n = 1$ this is the hypothesis of Theorem 11). Then $D^n x(t)$ exists as an ordinary stochastic process if and only if $p(0) = Dp(0) = \ldots = D^{n-1} p(0) = 0$ and

$$
\int_{0}^{\infty} |D^j p(s)|^2 \, ds < +\infty
$$

for $j = 0,\ 1,\, \ldots\, n$.

**Proof.** Argue by induction on $n$, applying {ref}`Theorem 11 <thm-msd-representation>` at each step. Suppose $p(0) = \ldots = D^{n-2} p(0) = 0$, so that the lower-order derivatives carry no impulse and $D^{n-1} x(t)$ has Wold kernel $D^{n-1} p$. Theorem 11 applied to $D^{n-1} x$ then says that $D^n x$ exists as an ordinary process if and only if $D^{n-1} p(0) = 0$ and $\int_0^\infty |D^n p(s)|^2\, ds < \infty$. Collecting these requirements across the induction gives the vanishing conditions $D^j p(0) = 0$ for $j = 0, \ldots, n-1$, together with $D^n p \in L^2$. The intermediate integrability conditions, $\int_0^\infty |D^j p|^2\, ds < \infty$ for $0 < j < n$, then follow from $p \in L^2$ and $D^n p \in L^2$ by interpolation.

We also have:

**Theorem 13.** If $x(t)$ is mean square differentiable and linearly indeterministic with $x(t) = \tilde P(D) w(t)$ being its Wold representation, then

$$
Dx(t) = D \tilde P(D) w(t)
$$

or

$$
Dx(t) = \int_{0}^{\infty} p'(\tau) w(t-\tau) \, d\tau
$$

is a Wold representation for $Dx(t)$.

**Proof.** Since $x$ is mean square differentiable, $p(0) = 0$, so $D\tilde P(D)$ carries no impulse and $Dx(t) = \int_0^\infty p'(\tau)\, w(t-\tau)\, d\tau$. If $\tilde P(s)$ has no zeroes in the open right half of the complex plane, then neither does $s \tilde P(s)$. The extra factor $s$ contributes one zero at the origin, on the imaginary axis rather than inside the right half plane, so $s\tilde P(s)$ meets the condition of the spectral factorization theorem of {doc}`08_spectral_densities` and $w$ is fundamental for $Dx$.

That zero at the origin is the boundary case flagged there. It costs something. The spectral
density of $Dx$ vanishes at $\omega = 0$, so the inverse filter
$1/[\,s\tilde P(s)\,]$ has a pole on the axis of integration and is not square integrable.
Recovering $w(t)$ from $[Dx(v),\, v \leq t]$ therefore requires a limit of square integrable
filters rather than a single one; the two closed linear spaces coincide, but $w(t)$ is not the
output of a square integrable filter applied to $Dx$. This is the continuous time counterpart of
a discrete time moving average with a unit root. Nothing later in the book turns on the
distinction, because the processes of {doc}`11_linear_sde` are differentiated only when
$n - 1 - m \geq 1$, and their derivatives are used as inputs rather than as observables.

We close with mean square continuity, which, in contrast to differentiability, imposes no
condition on $p$ beyond the square integrability already assumed in the Wold representation.

```{eval-rst}
.. index::
   single: mean square continuity; kernel criterion
   single: square integrable kernel
```

**Theorem 14.** Let $x(t)$ be a linearly indeterministic covariance stationary stochastic process with Wold moving average representation $x(t) = \tilde P(D) w(t) = \int_{0}^{\infty} p(\tau) w(t-\tau) \, d\tau$, where $\int_{0}^{\infty} p(\tau)^2 \, d\tau < \infty$. Then $x(t)$ is mean square continuous.

**Proof.** Extending $p(s) = 0$ for $s < 0$ and using $E\, w(t-s) w(t-s') = \delta(s-s')$,

$$
E\, [(x(t+\epsilon) - x(t))^2] = \int_{-\epsilon}^{\infty} |p(s+\epsilon) - p(s)|^2 \, ds = \big\| p(\cdot + \epsilon) - p \big\|_{L^2}^2 .
$$

Translation is continuous in $L^2$: $\| p(\cdot+\epsilon) - p \|_{L^2} \to 0$ as $\epsilon \to 0$ for every $p \in L^2$. Hence the right-hand side tends to $0$, and $x(t)$ is mean square continuous; equivalently, $R(\tau) = \int_{0}^{\infty} p(b)\, p(b+\tau)\, db$ is continuous at $\tau = 0$. Mean square continuity does *not* require $p$ to be pointwise continuous. A moving average against a kernel with jumps, such as the rectangular $p = \mathbf 1_{[0,1]}$, is mean square continuous all the same. Discontinuous kernels of just this kind are central to Marcet's analysis of temporal aggregation in {doc}`23_temporal_aggregation_streamlined`: a one-sided kernel with a jump at the origin, $p(0) \neq 0$, yields a process that is mean square continuous yet fails condition (a) of {ref}`Theorem 11 <thm-msd-representation>`. It is *not* mean square differentiable, hence locally unpredictable ({doc}`13_locally_unpredictable`). It is exactly these discontinuities that generate the aggregation distortions studied there.

Together with {ref}`Theorem 11 <thm-msd-representation>` and {ref}`Theorem 12 <thm-msd-higher-order>`, the following characterization can provide a useful way of testing for mean square differentiability of various orders.

**Criterion (Initial value theorem).** Let $\tilde P(s) = \int_{0}^{\infty} e^{-st} p(t) \, dt$ be the Laplace transform of $p(t)$. Then

$$
p(0) = \lim_{s \to \infty}\, s \tilde P(s).
$$

**Proof.** (For case where $p'(t)$ exists)

$$
\int_{0}^{\infty} p'(t) e^{-st} \, dt = s \tilde P(s) - p(0)
$$

As $s \to \infty$, the integral on the left approaches zero.

This initial value theorem is exactly the tool used in {doc}`11_linear_sde` to count how many
times the solution of an $n$-th order linear stochastic differential equation is mean square
differentiable, and again in {doc}`13_locally_unpredictable`, where $p(0) \neq 0$ certifies
that a process is locally unpredictable.

## Exercises

```{code-cell} ipython3
import numpy as np
from scipy.integrate import quad
```

```{exercise-start}
:label: msd_ex1
```

**Three routes to $p(0)$.** Theorem 11 gives the condition $p(0)=0$ for mean square
differentiability; the initial value theorem evaluates $p(0)$ from the transform; and the
autocovariance argument in the proof shows that $R'(0^+) = -\tfrac12 p(0)^2$. Check that the
three agree for the kernels

$$
p_a(\tau) = b\, e^{-a\tau}, \qquad p_b(\tau) = \tau\, e^{-a\tau},
$$

with $a = 1$, $b = 0.7$. For each, compute

(i) $p(0)$ directly;
(ii) $\lim_{s\to\infty} s\tilde P(s)$, using $\tilde P_a(s) = b/(s+a)$ and
$\tilde P_b(s) = (s+a)^{-2}$;
(iii) $R'(0^+)$ by numerically differentiating
$R(\tau) = \int_0^\infty p(u)p(u+\tau)\,du$, and compare with $-\tfrac12 p(0)^2$.

Which of the two processes is mean square differentiable? Confirm that $p_a$ is the
Ornstein–Uhlenbeck kernel of {doc}`07_wiener_driven_sde`, whose autocovariance
$R(\tau) = \frac{b^2}{2a}e^{-a|\tau|}$ has the kink discussed there.

```{exercise-end}
```

```{solution-start} msd_ex1
:class: dropdown
```

```{code-cell} ipython3
a, b = 1.0, 0.7
kernels = {
    "p_a = b e^{-a t}":  (lambda t: b*np.exp(-a*t),  lambda s: b/(s+a)),
    "p_b = t e^{-a t}":  (lambda t: t*np.exp(-a*t),  lambda s: 1.0/(s+a)**2),
}

for name, (p, P) in kernels.items():
    R    = lambda tau: quad(lambda u: p(u)*p(u+tau), 0, 200, limit=400)[0]
    h    = 1e-6
    Rp0  = (R(h) - R(0.0))/h                      # one-sided derivative R'(0+)
    print(f"{name}")
    print(f"   (i)   p(0)                = {p(0.0):+.8f}")
    print(f"   (ii)  s P(s) at s=1e8     = {1e8*P(1e8):+.8f}")
    print(f"   (iii) R'(0+)              = {Rp0:+.8f}   vs  -p(0)^2/2 = {-0.5*p(0.0)**2:+.8f}")
    print(f"   mean square differentiable: {abs(p(0.0)) < 1e-12}\n")
```

All three routes agree. For $p_a$ the kernel does not vanish at the origin, $R'$ jumps there by
$-p(0)^2$, and the process is not mean square differentiable. This is the Ornstein–Uhlenbeck case, with
$-\tfrac12 p(0)^2 = -b^2/2 = -aR(0)$ as the kink in $R(\tau) = \frac{b^2}{2a}e^{-a|\tau|}$
requires. For $p_b$ the kernel vanishes at the origin, $R'(0^+) = 0$, $R$ is smooth there, and
the process is once mean square differentiable.

```{solution-end}
```
