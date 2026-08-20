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

```{eval-rst}
.. index::
   single: linear stochastic differential equation; definition
   single: rational spectral density; and linear SDEs
   single: differential operator D
   single: characteristic roots
```

# 11. Linear Stochastic Differential Equations

```{eval-rst}
.. index::
   single: linear stochastic differential equation
   single: rational spectral density
   single: partial fraction expansion
   single: analytic function
```

We now consider a class of higher order stochastic differential equations that are formed
by driving a linear, constant coefficient $n^{th}$ order differential equation with a white
noise. Letting $D^j$ be the $j^{th}$ mean square derivative operator, we consider the
equation

```{math}
:label: eq-11-a1
\theta(D)\, x(t) = w(t)
```

or

```{math}
:label: eq-11-a2
(D^n + \theta_{n-1}\, D^{n-1} + \ldots + \theta_0)\, x(t) = w(t)
```

where

$$
E w(t)\, w(t-\tau) = \delta(\tau)
$$

and

$$
\begin{aligned}
\theta(D) &= (D^n + \theta_{n-1}\, D^{n-1} + \ldots + \theta_0) \\
&= (D - \lambda_1)\ (D - \lambda_2) \ldots (D - \lambda_n)
\end{aligned}
$$

We assume that the real parts of $\lambda_j$, which we denote $re(\lambda_j)$, are less
than zero for $j = 1, \ldots, n$. We assume that the $\lambda_j$'s, which are the zeroes of
$\theta(s)$, are distinct.

By writing {eq}`eq-11-a1` we intend to denote that $x(t)$ is the stochastic process given
by

```{math}
:label: eq-11-b
x(t) = \int^\infty_0 p(\tau)\, w(t - \tau)\, d\tau = \tilde P(D)\, w(t)
```

where

$$
\tilde P(D) \equiv \theta(D)^{-1} =\ \frac{1}{(D - \lambda_1)\, \cdots\, (D - \lambda_n)}.
$$

```{note}
The white noise $w(t)$ is not an ordinary function. Chapters 4 and 6 construct it as a
generalized process with the $\delta$-function autocovariance
$E\,w(t)\,w(t-\tau)=\delta(\tau)$. Its instantaneous variance is infinite, and it means
something only underneath an integral sign. The moving average {eq}`eq-11-b` is such an
integral. Read it as a mean square (Wiener) stochastic integral against the increments of the
Wiener process $W$ of Chapter 6. Write $dW(t) = w(t)\,dt$. Then

$$
x(t) = \int_0^\infty p(\tau)\, w(t-\tau)\, d\tau = \int_0^\infty p(\tau)\, dW(t-\tau),
$$

The right side is the $L^2$ limit of the approximating sums
$\sum_i p(\tau_i)\,[\,W(t-\tau_i) - W(t-\tau_{i+1})\,]$ built from the independent Gaussian
increments of $W$. The limit exists when $\int_0^\infty p(\tau)^2\,d\tau < \infty$. The Itô
isometry then gives $E\,x(t)^2 = \int_0^\infty p(\tau)^2\,d\tau$ and
$R(s) = \int_0^\infty p(\tau)\,p(\tau+s)\,d\tau$. Integration converts a generalized noise into
an ordinary process with finite variance, covariance stationary and mean square continuous. The
kernel is a causal linear filter that weights the recent past of the innovations by $p$. Square
integrability holds here because $p(\tau)=\sum_j g_j e^{\lambda_j\tau}$ decays exponentially.
Each factor $\frac{1}{D-\lambda_j}\,w(t) = \int_0^\infty e^{\lambda_j u}\,w(t-u)\,du$ is a
one-sided exponential smoother, the Ornstein–Uhlenbeck building block of Chapter 7. The
stationary solution sums these smoothers.

The kernel $p(\tau)$ is deterministic. It does not depend on the path of $W$, so how we form
the stochastic integral does not matter. The Itô and Stratonovich integrals differ by a
quadratic variation correction that appears when the integrand depends on the Brownian motion
against which it is integrated, as the $g(x)\,dW$ terms of Chapter 7 do. A deterministic
integrand kills that correction. Here the Itô integral, the Stratonovich integral, and the mean
square integral coincide.
```

```{eval-rst}
.. index::
   single: sum of exponentials kernel
   single: Ito isometry
   single: Stratonovich integral
```

A partial fraction representation of $\tilde P(D)$ is

$$
\tilde P(D) = \sum_{j=1}^{n}\ \frac{g_j}{D - \lambda_j}
$$

where

$$
g_j = \lim_{s\to \lambda_j}\ \tilde P(s)\ (s - \lambda_j).
$$

Now $1/(s-\lambda_j)$ for $re(\lambda_j) < 0$ is the Laplace transform of the function
$[e^{\lambda_j t},\, t \geq 0$; 0 for $t < 0]$. We express this transform pair relationship
as

$$
\frac{1}{s - \lambda_j}\ \leftrightarrow\ e^{\lambda_j\, t}.
$$

This relationship justifies the operational calculus

$$
\frac{1}{D - \lambda_j}\ w(t) = \int^\infty_0 e^{\lambda_j u}\ w(t-u)\, du
$$

for $re(\lambda_j) < 0$.

Thus, we can express {eq}`eq-11-b` as

$$
x(t) = \sum_{j=1}^{n}\ \frac{g_j}{D - \lambda_j}\ w(t)
$$

or

$$
x(t) = \sum_{j=1}^{n} g_j \int^\infty_0 e^{\lambda_j u}\ w(t - u)\, du.
$$

```{eval-rst}
.. index::
   single: stationary solution; of a linear SDE
   single: mean square differentiability; the count n-1-m
```

This equation is said to represent the *stationary solution* of the stochastic differential
equation. These rational-spectral-density processes are precisely the ones that admit a
finite-dimensional state-space realization; {doc}`15_kalman_filter_spectral_factorization`
recasts their Wold representation and spectral factorization as the time-domain output of the
Kalman–Bucy filter.

We now show that the solution is mean square differentiable $n - 1$ times, but not $n$
times, by applying the initial value theorem of {doc}`09_characterizations_ms_differentiability`,
$p(0) = \lim_{s\to\infty} s\tilde P(s)$. The Laplace transform of $p(\tau)$ is $\tilde P(s) = 1/[(s - \lambda_1) \ldots
(s - \lambda_n)]$, that of $p'(\tau)$ is $s \tilde P(s) = s/[(s - \lambda_1) \ldots
(s - \lambda_n)]$, and so on. For $x(t)$ to be mean square differentiable, we require that

$$
p(0) = \lim_{s\to \infty}\ s \tilde P(s) = \lim_{s\to \infty}\ \frac{s}{(s-\lambda_1) \ldots (s - \lambda_n)} = 0,
$$

which obtains so long as $n \geq 2$. For $x(t)$ to be $j$ times mean square
differentiable, we require that

$$
p^{(j - 1)}(0) = \lim_{s \to \infty}\ s^j\, \tilde P(s) = \lim_{s \to \infty}\ \frac{s^j}{(s - \lambda_1) \ldots (s - \lambda_n)} = 0
$$

which holds for $j < n$, but is violated for $j \geq n$.

Therefore, where $x(t)$ is governed by the solution of the $n^{th}$ order stochastic
differential equation, $x(t)$ is mean square differentiable up to orders $(n-1)$, but is
not differentiable $n$ times. (The reader can verify that in the present case of a rational
$p(\tau)$ function, $\int^\infty_0 |D^j p(\tau)|^2\, d\tau < + \infty$ for all $j$, so that
the second condition for mean square differentiability is satisfied.)

Next, we consider a broader class of linear stochastic differential equations generated by
driving an $n^{th}$ order constant coefficient linear differential equation with a linear
combination of a white noise and $m \leq n - 1$ generalized derivatives of this same white
noise. The differential equation is

$$
\theta(D)\, x(t) = \Psi(D)\, w(t)
$$

where

$$
\begin{aligned}
\theta(D) &= (D^n + \theta_{n-1}\, D^{n-1}+\ldots + \theta_0) \\
&= (D - \lambda_1)\ (D - \lambda_2) \ldots (D - \lambda_n),\ re(\lambda_j)< 0,\ j = 1,\ldots, n \\
\Psi(D) &= \Psi_0 + \Psi_1 D + \ldots + \Psi_m D^m,\ \Psi_0 = 1,\ m \leq n-1. \\
E w(t) &\, w(t - \tau) = \delta(\tau).
\end{aligned}
$$

In this equation, $w(t)$ is a fundamental white noise for $x(t)$ if and only if the zeroes
of $\Psi(s)$ lie in the left hand of the complex plane. By the stationary solution of this
equation we mean

$$
x(t) = \tilde P(D)\, w(t) = \int^\infty_0 p(\tau)\, w(t-\tau)\, d\tau
$$

where

$$
\tilde P(s) = \theta(s)^{-1}\, \Psi(s).
$$

A partial fraction representation of $\tilde P(s)$ is

$$
\tilde P(s) = \sum_{j=1}^{n}\ \frac{g_j}{s - \lambda_j}
$$

where $g_j = \lim_{s \to \lambda_j}\ \tilde P(s)\ (s-\lambda_j)$. Therefore $p(\tau) =
\sum_{j=1}^{n}\ g_j\ e^{\lambda_j \tau}$. Proceeding as above, once again we see that the
solution has representation

$$
x(t) = \sum_{j=1}^{n}\ g_j \int^\infty_0 e^{\lambda_j u}\, w(t-u)\, du.
$$

Evidently, $p(\tau) = \sum_{j=1}^{n} g_j\, e^{\lambda_j \tau}$ is continuous in $\tau$.
Therefore $x_t$ is mean square continuous. To test for mean square differentiability, we
must evaluate

$$
p(0) = \lim_{s\to \infty}\, s \tilde P(s) = \lim_{s\to \infty}\ \frac{\Psi_0 s + \Psi_1 s^2 + \ldots + \Psi_m s^{m+1}}{\theta_0 + \theta_1 s + \ldots + s^n}
$$

which equals zero only if $n > m + 1$. Proceeding in this way, it can be verified that
$x(t)$ is mean square differentiable $n - 1 - m$ times, but no more.

```{note}
The kernels $p(\tau) = \sum_j g_j\, e^{\lambda_j \tau}$ that rational linear stochastic
differential equations deliver are analytic. The processes of this chapter are therefore mean
square continuous, and their differentiability depends only on how many of the leading values
$p(0), p'(0), \ldots$ vanish. That is the count $n - 1 - m$ just obtained, an application of the
criteria of {doc}`09_characterizations_ms_differentiability`.
{doc}`23_temporal_aggregation_streamlined` leaves this smooth class. Marcet posits one-sided
kernels that are discontinuous, the leading case being a jump at the origin, $p(0) \neq 0$. Such
a process remains mean square continuous, because Theorem 14 of
{doc}`09_characterizations_ms_differentiability` asks only that the kernel be square integrable,
not that it be continuous. It is not mean square differentiable, so it is locally unpredictable
in the sense of {doc}`13_locally_unpredictable`. Those discontinuities, and the jumps they force
on the sampled discrete-time kernel at every integer, drive the aliasing and the manufactured or
destroyed Granger causality that Chapter 23 analyzes.
```

The following theorem pulls together a number of our earlier results.

```{eval-rst}
.. index::
   single: analytic function; rigidity of
   single: identity theorem
   single: degenerate process
```

**Theorem 15.** Let $x(t)$ be a covariance stationary, zero mean, linearly indeterministic
stochastic process, with Wold representation

$$
x(t) = \int^\infty_0 p(\tau)\, w(t - \tau)\, d\tau,
$$

and covariogram $E x(t)\, x(t-\tau) = R(\tau)$. Suppose that $R(\tau)$ *and* the Wold kernel
$p(\tau)$ are analytic functions, $R$ having a convergent Taylor series representation

```{math}
:label: eq-11-star
R(\tau) = \sum_{n=0}^{\infty} R^{(n)}(0)\ \frac{\tau^n}{n!}
```

where $R^{(n)}(\tau) \equiv \frac{d^n}{d\tau^n}\ R(\tau)$. Then it follows that $x(t)$ is
mean square differentiable $k$ times, for every $k \geq 1$ (the letter $m$ is reserved in this
chapter for the degree of $\Psi$). Indeed, it follows
from {eq}`eq-11-star` that

$$
R^{(2k)}(\tau) = \sum_{n=2k}^{\infty}\, R^{(n)}(0)\ \frac{\tau^{(n-2k)}}{(n - 2k)!}\, ,
$$

so that $R^{(2k)}(\tau)$ exists and itself has a convergent Taylor series representation.
Furthermore, $p(\tau)$ is by hypothesis analytic and so has a convergent
Taylor series representation

$$
p(\tau) = \sum_{n=0}^{\infty} p^{(n)}(0)\ \frac{\tau^n}{n!}\, .
$$

Since $x(t)$ is mean square differentiable an infinite number of times, it follows that
$p^{(n)}(0) = 0$ for all $n \geq 0$, and therefore that $p(\tau) = 0$ for all $\tau \geq 0$.
Therefore $x(t) = 0$ for all $t$.

The theorem states that the only covariance stationary linearly indeterministic process
with analytic autocovariance function $R(\tau)$ and analytic Wold integrating kernel
$p(\tau)$ is the trivial process $x(t) = 0$.

```{note}
A real-valued function $R$ of a real variable is analytic on an interval if around every point
$\tau_0$ of that interval it agrees with a convergent power series

$$
R(\tau) = \sum_{n=0}^{\infty} \frac{R^{(n)}(\tau_0)}{n!}\,(\tau - \tau_0)^n,
$$

having a strictly positive radius of convergence. Equivalently, $R$ is the restriction to the
real axis of a function holomorphic on an open neighborhood of that interval in the complex
plane $\mathbb{C}$. A theorem of complex analysis says that a function complex-differentiable
once on an open set is infinitely complex-differentiable there and equals its own Taylor series.
For complex functions, "differentiable," "infinitely differentiable," and "equal to a convergent
power series" say the same thing. For real functions they do not.

Analytic functions are rigid. The identity theorem says that an analytic function on a connected
domain is determined by its values on any subset having a limit point. If all its derivatives
vanish at a single point, or if it vanishes on an interval however short, then it vanishes
identically.

That rigidity drives Theorem 15, and it asks more than smoothness. A function can be infinitely
differentiable ($C^\infty$) without being analytic. Take $f(\tau) = e^{-1/\tau^2}$ with
$f(0)=0$. Every $f^{(n)}(0)$ vanishes, yet $f$ is not identically zero, so $f$ is not analytic
at the origin. Smoothness of $p(\tau)$ would therefore permit a nonzero kernel all of whose
derivatives vanish at $\tau = 0$. Analyticity of $R$, inherited by $p$, upgrades
$p^{(n)}(0) = 0$ for all $n$ into $p \equiv 0$, and so into $x(t) \equiv 0$.
```

## Exercises

```{code-cell} ipython3
import numpy as np
```

```{exercise-start}
:label: lsde_ex1
```

**The differentiability count.** The text shows that $\theta(D)x = \psi(D)w$ delivers a process
that is mean square differentiable exactly $n - 1 - m$ times, where $n = \deg\theta$ and
$m = \deg\psi$. Verify this two ways for the third-order system

$$
(D^3 + .6D^2 + .4D + .2)\, x(t) = w(t), \qquad n = 3,\ m = 0,
$$

which is the numerical example of {doc}`18_time_aggregation_var`, Table 1.

(a) Compute the roots $\lambda_j$ of $\theta$, form the residues
$g_j = \lim_{s\to\lambda_j}(s-\lambda_j)\tilde P(s)$, and build
$p(\tau) = \sum_j g_j e^{\lambda_j\tau}$. Confirm numerically that
$p(0) = p'(0) = 0$ while $p''(0) \neq 0$, so that the process is twice but not three times mean
square differentiable, i.e. $n-1-m = 2$.

(b) Confirm the same count with the initial value theorem of
{doc}`09_characterizations_ms_differentiability`, by evaluating $s^j \tilde P(s)$ at large
real $s$ for $j = 1, 2, 3$ and seeing which limits vanish.

(c) Now put a numerator of degree $m = 2$ on the same denominator, say
$\psi(s) = 1 + s + s^2$. Show that $p(0) \neq 0$: the process is differentiable
$n - 1 - m = 0$ times, and is therefore locally unpredictable in the sense of
{doc}`13_locally_unpredictable`.

```{exercise-end}
```

```{solution-start} lsde_ex1
:class: dropdown
```

```{code-cell} ipython3
theta = np.array([1, .6, .4, .2])          # s^3 + .6 s^2 + .4 s + .2
lam = np.roots(theta)
print("roots of theta:", np.round(lam, 6))

def residues(lam, psi=np.array([1.0])):
    """g_j = psi(lam_j) / prod_{k != j} (lam_j - lam_k)."""
    g = []
    for j, lj in enumerate(lam):
        denom = np.prod([lj - lk for k, lk in enumerate(lam) if k != j])
        g.append(np.polyval(psi, lj) / denom)
    return np.array(g)

g = residues(lam)
# derivatives of p at 0:  p^(r)(0) = sum_j g_j lam_j^r
for r in range(4):
    val = np.sum(g * lam**r)
    print(f"p^({r})(0) = {val.real:+.10f}  (imag {val.imag:+.1e})")
```

```{code-cell} ipython3
# (b) initial value theorem: p^{(j-1)}(0) = lim_{s->inf} s^j P(s)
Ptil = lambda s: 1.0/np.polyval(theta, s)
for j in [1, 2, 3]:
    print(f"j={j}:  s^{j} P(s) at s=1e6 -> {1e6**j * Ptil(1e6):.6e}")
```

```{code-cell} ipython3
# (c) numerator of degree m = 2 on the same denominator
psi = np.array([1.0, 1.0, 1.0])            # s^2 + s + 1
g2 = residues(lam, psi)
print(f"with m=2:  p(0) = {np.sum(g2).real:+.6f}   (n-1-m = 0 derivatives)")
print(f"check by initial value theorem: s*P(s) at s=1e6 -> "
      f"{1e6*np.polyval(psi,1e6)/np.polyval(theta,1e6):.6f}")
```

Both routes agree. With $\psi = 1$ the first two derivatives of $p$ vanish at the origin and the
third does not, so $x$ is twice mean square differentiable. That is why
{doc}`18_time_aggregation_var` finds the sampled moving-average coefficients of this process so
badly behaved. The kernel is *too smooth* near zero for the discrete innovation to resemble the
continuous one. Raising the numerator degree to $m = n-1 = 2$ moves the process to the other
side of the dividing line, where $p(0)\neq 0$ and local unpredictability sets in.

```{solution-end}
```
