# 11. Linear Stochastic Differential Equations

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
**What integrating the white noise against $p(\tau)$ does.** The white noise $w(t)$ is not an
ordinary function. As in Chapters 4 and 6, it is a *generalized* process: it has the
$\delta$-function autocovariance $E\,w(t)\,w(t-\tau)=\delta(\tau)$, hence infinite instantaneous
variance, and it acquires a well-defined meaning only underneath an integral sign. The
moving-average integral {eq}`eq-11-b` is exactly such an object, and it is most precisely read
as a **mean square (Wiener) stochastic integral** against the increments of the Wiener process
$W$ of Chapter 6. Writing $dW(t) = w(t)\,dt$,

$$
x(t) = \int_0^\infty p(\tau)\, w(t-\tau)\, d\tau = \int_0^\infty p(\tau)\, dW(t-\tau),
$$

is defined as the mean square (that is, $L^2$) limit of the approximating sums
$\sum_i p(\tau_i)\,[\,W(t-\tau_i) - W(t-\tau_{i+1})\,]$ built from the *independent* Gaussian
increments of $W$. Provided $\int_0^\infty p(\tau)^2\,d\tau < \infty$ — which holds here because
$p(\tau)=\sum_j g_j e^{\lambda_j\tau}$ decays exponentially — this limit exists, and the integral
converts the unruly generalized noise into an ordinary, finite-variance, covariance stationary,
mean square continuous process. By the Itô isometry its variance is
$E\,x(t)^2 = \int_0^\infty p(\tau)^2\,d\tau$ and its autocovariance is
$R(s) = \int_0^\infty p(\tau)\,p(\tau+s)\,d\tau$. In short, the kernel acts as a *causal linear
filter* that **colors** — smooths — the white noise, weighting the recent past of the
innovations by $p$. Each elementary factor
$\frac{1}{D-\lambda_j}\,w(t) = \int_0^\infty e^{\lambda_j u}\,w(t-u)\,du$ is a one-sided
exponential smoother — the Ornstein–Uhlenbeck building block of Chapter 7 — and the stationary
solution is a sum of these smoothers.

**Itô versus Stratonovich.** Because the kernel $p(\tau)$ is *deterministic* — it does not depend
on the path of $W$ — there is no ambiguity in how the stochastic integral is formed. The Itô and
Stratonovich integrals differ only by the quadratic-variation correction term that appears when
the integrand depends on the very Brownian motion against which it is integrated (as in the
nonlinear $g(x)\,dW$ terms of Chapter 7); for a deterministic integrand that correction term
vanishes. Here, therefore, the Itô integral, the Stratonovich integral, and the ordinary mean
square integral all coincide, and {eq}`eq-11-b` may be read in whichever of these senses the
reader prefers.
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

This equation is said to represent the *stationary solution* of the stochastic differential
equation. These rational-spectral-density processes are precisely the ones that admit a
finite-dimensional state-space realization; {doc}`13_kalman_filter_spectral_factorization`
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
\Psi(D) &= 1 + \Psi_1 D + \ldots + \Psi_m D^m,\ m \leq n-1. \\
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
**Smooth kernels here; discontinuous kernels in Chapter 23.** The kernels
$p(\tau) = \sum_j g_j\, e^{\lambda_j \tau}$ delivered by rational linear stochastic differential
equations are not merely continuous but *analytic*, so the processes of this chapter are
automatically mean square continuous, and their differentiability is decided entirely by how many
of the leading values $p(0), p'(0), \ldots$ vanish — the count $n - 1 - m$ just obtained, an
application of the criteria of {doc}`09_characterizations_ms_differentiability`.
{doc}`23_temporal_aggregation_streamlined` deliberately leaves this smooth class: Marcet posits
one-sided kernels that are *discontinuous*, the leading case being a jump at the origin,
$p(0) \neq 0$. Such a process is still mean square continuous — Theorem 14 of
{doc}`09_characterizations_ms_differentiability` requires only square integrability of the
kernel, not its continuity — yet it is no longer mean square differentiable, and so is *locally
unpredictable* in the sense of {doc}`15_locally_unpredictable`. It is exactly these
discontinuities, and the jumps they force on the sampled discrete-time kernel at every integer,
that drive the aliasing and the manufactured-or-destroyed Granger causality analyzed there.
```

The following theorem pulls together a number of our earlier results.

**Theorem 15.** Let $x(t)$ be a covariance stationary, zero mean, linearly indeterministic
stochastic process, with Wold representation

$$
x(t) = \int^\infty_0 p(\tau)\, w(t - \tau)\, d\tau,
$$

and covariogram $E x(t)\, x(t-\tau) = R(\tau)$. Suppose that $R(\tau)$ is an analytic
function, having a convergent Taylor series representation

```{math}
:label: eq-11-star
R(\tau) = \sum_{n=0}^{\infty} R^{(n)}(0)\ \frac{\tau^n}{n!}
```

where $R^{(n)}(\tau) \equiv \frac{d^n}{d\tau^n}\ R(\tau)$. Then it follows that $x(t)$ is
mean square differentiable $m$ times, for any $1 \leq m \leq \infty$. Indeed, it follows
from {eq}`eq-11-star` that

$$
R^{(2m)}(\tau) = \sum_{n=2m}^{\infty}\, R^{(n)}(0)\ \frac{\tau^{(n-2m)}}{(n - 2m)!}\, ,
$$

so that $R^{(2m)}(\tau)$ exists and itself has a convergent Taylor series representation.
Furthermore, under the assumed conditions, $p(\tau)$ is analytic and has a convergent
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
**What "analytic" means.** A real-valued function $R$ of a real variable is *(real-)analytic*
on an interval if, around every point $\tau_0$ of that interval, it agrees with a convergent
power series,

$$
R(\tau) = \sum_{n=0}^{\infty} \frac{R^{(n)}(\tau_0)}{n!}\,(\tau - \tau_0)^n,
$$

the series having a strictly positive radius of convergence. Equivalently — and this is the
precise sense the term carries in *complex analysis* — $R$ is the restriction to the real axis
of a function that is **holomorphic** (complex-differentiable) on an open neighborhood of that
interval in the complex plane $\mathbb{C}$. A central theorem of complex analysis is that a
function complex-differentiable even *once* on an open set is automatically infinitely
complex-differentiable there and equals its own Taylor series; so for complex functions
"differentiable," "infinitely differentiable," and "equal to a convergent power series" all
coincide — in sharp contrast to the real case.

Analytic functions are **rigid**. By the identity theorem, an analytic function on a connected
domain is completely determined by its values on any subset that has a limit point; in
particular, if all of its derivatives vanish at a single point — or if it vanishes on any
interval, however small — then it vanishes *identically*.

This rigidity is precisely what powers Theorem 15, and it is strictly stronger than mere
smoothness. A function can be infinitely differentiable ($C^\infty$) without being analytic:
the standard example $f(\tau) = e^{-1/\tau^2}$ (with $f(0)=0$) has $f^{(n)}(0)=0$ for every $n$,
yet is not identically zero, so it is *not* analytic at the origin. Mere smoothness of $p(\tau)$
would therefore permit a nonzero kernel all of whose derivatives vanish at $\tau = 0$; it is the
assumed *analyticity* of $R$ — inherited by $p$ — that upgrades $p^{(n)}(0) = 0$ for all $n$ into
$p \equiv 0$, and hence into the degenerate conclusion $x(t) \equiv 0$.
```
