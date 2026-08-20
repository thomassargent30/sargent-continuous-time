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

# 12. Linear Least Squares Prediction

```{eval-rst}
.. index::
   single: Wiener-Kolmogorov prediction formula
   single: annihilation operator
   single: geometric distributed lead
   single: present value
   single: prediction; as projection
```

```{eval-rst}
.. index::
   single: prediction; Wiener-Kolmogorov formula
   single: projection; linear least squares
   single: forecast; optimal linear
   single: prediction error
```

## (a) The Wiener–Kolmogorov Formula

The Wold moving average representation is useful for representing the linear
least squares $u$-step ahead prediction for a purely linearly indeterministic
process. Let $x(t)$ be a purely linearly indeterministic process with Wold
moving average representation

$$
x(t) = \int^\infty_0 p(s)\, w(t-s)\, ds
$$ (eq-12-wold)

where $\int^\infty_0 p(s)^2\, ds < +\infty$ and $w(t)$ is a fundamental white
noise for $x(t)$. Recall that the property that $w(t)$ is a fundamental white
noise for $x(t)$ means equality of the linear spaces

$$
\begin{aligned}
H_x(-\infty,\, t) \equiv\ & [\, y(t) : y(t) = \int^\infty_0 b(s)\, x(t-s)\, ds, \\
&\ \text{ for any }\ b(s) \in L_2\,[0,\, \infty) \,]
\end{aligned}
$$

and

$$
\begin{aligned}
H_w(-\infty,\, t) \equiv\ & [\, z(t) : z(t) = \int^\infty_0 h(s)\, w(t-s)\, ds \\
&\ \text{ for any }\ h(s) \in L_2\,[0,\, \infty) \,]
\end{aligned}
$$

where $L_2\,[0,\, \infty)$ is the space of square integrable functions, i.e.,
functions $b(s)$ for $0 \leq s < \infty$ such that $\int^\infty_0 b(s)^2\, ds <
+\infty$. The equality of these spaces means that lagged $x$'s contain the
same amount of information as lagged $w$'s.

Since {eq}`eq-12-wold` holds for all $t$, we have

$$
x(t + u) = \int^\infty_{s=-u} p(s+u)\, w(t-s)\, ds,\ \text{ for }\ u \geq 0.
$$

Using the identity of the linear spaces $H_x(-\infty,\,t)$ and $H_w(-\infty,\,t)$, we have that

$$
E\, [\, x(t+u) \mid x(v),\, v \leq t \,] = \int^\infty_{s=0} p(s+u)\, w(t-s)\, ds.
$$ (eq-12-wk)

Equation {eq}`eq-12-wk` is the continuous time Wiener–Kolmogorov formula. It rests on the
equality of the spaces spanned by past $x$'s and past fundamental innovations $w$'s established
in Wold's theorem (Theorem 10 of {doc}`08_spectral_densities`);
{doc}`15_kalman_filter_spectral_factorization` gives the equivalent state-space form, in which
$w$ becomes the Kalman innovations process and this same forecast is computed recursively.

Using operational calculus, the formula can be expressed as

$$
E_t\, x(t+u) = [\, \tilde P(D)\, e^{Du} \,]_+\, w(t)
$$ (eq-12-wkop)

where $[\, \tilde P(s)\, e^{su} \,]_+$ is the time function formed by taking
the inverse Laplace transform of $\tilde P(s)\, e^{su}$, and then multiplying
it by the Heaviside unit step function (i.e., setting values of the time
function for $t < 0$ equal to zero, while leaving values of the function for $t
\geq 0$ unaltered). The operator $[\,\cdot\,]_+$ is known as the *annihilation
operator*. Note that $e^{uD}$ is the operator that shifts a time function
*ahead* $u$ units. Reading property 4 (Delay) of Table 2 in the direction of an
advance rather than a delay, $e^{su}\, \tilde P(s)$ is the Laplace transform of
the advanced function $p(t+u)$; annihilating its negative
part then leaves the kernel $p(s+u),\ s \geq 0$, of {eq}`eq-12-wk`. This
$e^{uD}$ convention is the one used throughout the book, in
{doc}`19_prediction_formulas_continuous_time` and in
{doc}`20_aggregation_inverse_optimal_predictor`.

As a check, for $\tilde P(s) = 1/(a+s)$ the inverse transform of
$e^{su}/(a+s)$ is $e^{-a(t+u)}$ for $t + u \geq 0$; annihilating $t<0$ and
transforming back gives $e^{-au}/(a+s)$, which is {eq}`eq-12-ar1pred`.

As an example of the use of formula {eq}`eq-12-wk`, let $x(t)$ be governed by the
first order stochastic differential equation

$$
(D+a)\, x(t) = w(t),\qquad a > 0
$$

so that $p(\tau) = e^{-a\tau}$. Then formula {eq}`eq-12-wk` gives

$$
\begin{aligned}
E_t\, x(t+u) &= \int^\infty_0 e^{-a(s+u)}\, w(t-s)\, ds \\
&= e^{-au} \int^\infty_0 e^{-as}\, w(t-s)\, ds = e^{-au}\ \frac{1}{D + a}\, w(t)
\end{aligned}
$$

or

$$
E_t\, x(t+u) = e^{-au}\, x(t).
$$ (eq-12-ar1pred)

```{eval-rst}
.. index::
   single: geometric distributed lead; prediction of
   single: present value; of a forecast
   single: discounting; in prediction formulas
```

## (b) A Formula for Predicting "Geometric Distributed Leads"

Such geometric distributed leads are the present values that appear in every asset-pricing
equation, permanent-income model, and quadratic-adjustment-cost Euler equation; they are also
the continuous-time counterpart of the discounted expected sums that define the optimal
feedforward decision rules of {doc}`16_faster_methods_recursive_linear_models`. Evaluating them
is the central computational step in solving continuous-time rational expectations models.

In linear rational expectation models, there often appear terms of the form

$$
E_t \int^\infty_0 e^{\rho u}\, x(t+u)\, du \qquad re(\rho) < 0
$$ (eq-12-glead)

where $x(t)$ is a covariance stationary stochastic process. Where $x(t)$ is
governed by the first order Markov process $(D+a)\, x(t) = w(t)$, equation {eq}`eq-12-ar1pred`
implies that the linear least squares forecast of the geometric distributed
lead {eq}`eq-12-glead` is given by

$$
\int^\infty_0 e^{\rho u}\, E_t\, x(t+u)\, du = \left( \int^\infty_0 e^{\rho u}\, e^{-au}\, du \right) x(t) = \frac{1}{a - \rho}\, x(t).
$$

An approach to the evaluation of {eq}`eq-12-glead` which readily generalizes to $x(t)$'s
governed by higher order linear differential equations is as follows. Denote
the geometric distributed lead to be forecast as

$$
\begin{aligned}
x(t)^{\ast} &= \int^\infty_0 e^{\rho u}\, x(t+u)\, du \\
x(t)^{\ast} &= \left( \frac{-1}{\rho + D} \right)\ \left( \frac{1}{a+D} \right)\, w(t)
\end{aligned}
$$ (eq-12-xstar)

where for $re(\rho) < 0$, $-1/(\rho + s)$ is the (two-sided) Laplace transform
of the time function equal to $e^{-\rho u}$ for $u \leq 0$ and $0$ for $u > 0$.

Obtaining a partial fraction representation of the right side of {eq}`eq-12-xstar` gives

$$
x(t)^{\ast} = \frac{1}{a-\rho}\ \left[ \left( \frac{-1}{\rho + D} \right) + \left( \frac{1}{a+D} \right) \right]\ w(t)
$$

or

$$
x(t)^{\ast} = \left( \frac{1}{a-\rho} \right)\ \left[ -\int^\infty_0 e^{\rho s}\, w(t+s)\, ds + \int^\infty_0 e^{-as}\, w(t-s)\, ds \right]
$$

It then follows that

$$
\begin{aligned}
E_t\, x(t)^{\ast} &= \left( \frac{1}{a-\rho} \right) \int^\infty_0 e^{-as}\, w(t-s)\, ds = \left( \frac{1}{a-\rho} \right)\ \frac{1}{a + D}\ w(t) \\
&= \left( \frac{1}{a - \rho} \right) x(t).
\end{aligned}
$$ (eq-12-Estar)

This approach generalizes readily as follows. Represent equation {eq}`eq-12-Estar` as

$$
\begin{aligned}
E_t\, x(t)^{\ast} &= \left( \frac{1}{a - \rho} \right)\ \left( \frac{1}{a+D} \right)\, w(t) \\
E_t\, x(t)^{\ast} &= \left[ \frac{-\tilde P(D) + \tilde P(-\rho)}{D+\rho} \right]\, w(t)
\end{aligned}
$$ (eq-12-Estarop)

where $\tilde P(D) = 1/(a+D)$. As it happens, Equation {eq}`eq-12-Estarop` holds for *any*
$\tilde P(D)$, where $\tilde P(s)$ is the Laplace transform of a
squared summable function $p(\tau)$ concentrated on $\tau \in [0,\, \infty)$.
Thus, where

$$
x(t) = \int^\infty_0 p(\tau)\, w(t-\tau)\, d\tau
$$

we claim that the generalization of {eq}`eq-12-Estarop` is

$$
E_t \int^\infty_0 e^{\rho s}\, x(t+s)\, ds = \left[ \frac{-\tilde P(D) + \tilde P(-\rho)}{D+\rho} \right]\, w(t).
$$ (eq-12-gen)

The general formula {eq}`eq-12-gen`, valid for any rational $\tilde P(D)$, is established in
{doc}`19_prediction_formulas_continuous_time`.

## Exercises

There is a direct route to the kernel of a geometric distributed lead that provides an
independent check on {eq}`eq-12-gen`. Writing $x(t+s) = \int_0^\infty p(\tau)w(t+s-\tau)d\tau$
and discarding the terms dated after $t$ leaves
$E_t\, x(t+s) = \int_0^\infty p(\tau+s)\, w(t-\tau)\, d\tau$, which is {eq}`eq-12-wk`.
Multiplying by $e^{\rho s}$ and integrating,

```{math}
:label: eq-12-gker
E_t \int_0^\infty e^{\rho s} x(t+s)\, ds = \int_0^\infty g(\tau)\, w(t-\tau)\, d\tau,
\qquad
g(\tau) = \int_0^\infty e^{\rho u}\, p(\tau + u)\, du .
```

Formula {eq}`eq-12-gker` computes the same object as {eq}`eq-12-gen`, but in the time domain
and by quadrature rather than by partial fractions.

```{code-cell} ipython3
import numpy as np
from scipy.integrate import quad
```

```{exercise-start}
:label: pred_ex1
```

Take the second-order process $(D+a_1)(D+a_2)x(t) = w(t)$ with $a_1 = 0.6$, $a_2 = 1.7$, whose
kernel is $p(\tau) = (e^{-a_1\tau} - e^{-a_2\tau})/(a_2 - a_1)$, and let $\rho = -0.4$.

(a) Show by partial fractions that {eq}`eq-12-gen` gives

$$
g(\tau) = \frac{1}{a_2 - a_1}\left[\frac{e^{-a_1\tau}}{a_1 - \rho} - \frac{e^{-a_2\tau}}{a_2 - \rho}\right],
$$

and verify this closed form against the quadrature {eq}`eq-12-gker`.

(b) Confirm that $g(0) = \tilde P(-\rho) = 1/[(a_1-\rho)(a_2-\rho)]$, as the initial value
theorem of {doc}`09_characterizations_ms_differentiability` requires. Note that
$g(0) \neq 0$ even though $p(0) = 0$. {doc}`13_locally_unpredictable` exploits that fact.

(c) Recover the first-order case by letting $a_2 \to \infty$. Note that the kernel itself
vanishes in that limit, so the process must be rescaled: $a_2\, p(\tau) \to e^{-a_1\tau}$, and
correspondingly $a_2\, g(0) \to 1/(a_1 - \rho)$, which is the coefficient in
{eq}`eq-12-Estar`.

```{exercise-end}
```

```{solution-start} pred_ex1
:class: dropdown
```

```{code-cell} ipython3
a1, a2, rho = 0.6, 1.7, -0.4
p     = lambda t: (np.exp(-a1*t) - np.exp(-a2*t))/(a2 - a1)
g_cf  = lambda t: (np.exp(-a1*t)/(a1-rho) - np.exp(-a2*t)/(a2-rho))/(a2 - a1)
Ptil  = lambda s: 1.0/((s+a1)*(s+a2))

print("  tau     g closed form      g by quadrature")
for tau in [0.0, 0.5, 1.5, 4.0]:
    num, _ = quad(lambda u: np.exp(rho*u)*p(tau+u), 0, 200, limit=400)
    print(f"{tau:5.1f}   {g_cf(tau):15.10f}   {num:15.10f}")

print(f"\n(b) g(0)      = {g_cf(0.0):.10f}")
print(f"    P(-rho)   = {Ptil(-rho):.10f}")
print(f"    p(0)      = {p(0.0):.10f}   <- vanishes, but g(0) does not")
```

```{code-cell} ipython3
# (c) as a2 -> infinity, rescaled by a2, the process becomes first order with decay a1
print(f"target 1/(a1-rho) = {1/(a1-rho):.6f}")
for A2 in [10.0, 100.0, 1000.0, 1e5]:
    g0 = 1.0/((a1-rho)*(A2-rho))          # = g(0) with a2 replaced by A2
    print(f"  a2={A2:9.0f}:  a2*g(0) = {A2*g0:.6f}")
```

The quadrature reproduces the partial-fraction formula to ten digits, and $g(0)$ equals
$\tilde P(-\rho)$ exactly. Part (c) shows the first-order limit: the second exponential decays
ever faster and contributes ever less, so after rescaling by $a_2$ the present value collapses
to $x(t)/(a_1-\rho)$, which is {eq}`eq-12-Estar`.

```{solution-end}
```
