# 13. Locally Unpredictable Stochastic Processes

Lack of mean square differentiability implies that a process is "locally
unpredictable" or is "locally like a martingale." Sims used the concept of
local unpredictability in his work on asset prices, interest rates, and
consumption. We use the following definition:

**Definition 10.** A stochastic process $x(t)$ with finite second moments is said to
be *locally unpredictable* if

$$
\lim_{\delta \to 0}\ \frac{E_t (x(t+\delta) - E_t x(t+\delta))^2}{E_t (x(t+\delta) - x(t))^2}\ = 1.
$$

This definition makes precise the sense in which $x(t)$ is locally a martingale,
i.e., the sense in which for small $\delta$,

$$
E_t x(t+\delta) \simeq x(t).
$$

We have the following theorem.

**Theorem 16.** Let $x(t)$ be a linearly indeterministic covariance stationary
stochastic process with Wold representation

$$
x(t) = \int^\infty_0 p(\tau)\, w(t-\tau)\, d\tau
$$

where $w(t)$ is a fundamental white noise for the $x(t)$ process, and
$Ew(t) w(t-\tau) = \delta(\tau)$. Assume that $p(\tau)$ is continuous and
continuously differentiable at least twice. Then if $x(t)$ is not mean square
differentiable, $x(t)$ is locally unpredictable.

**Proof:** From the Wiener-Kolmogorov prediction formula, we have that

$$
\frac{E_t (x(t+\delta) - E_t x(t+\delta))^2}{E_t (x(t+\delta) - x(t))^2}
$$

$$
= \frac{\int^\infty_0 p(s)^2\, ds}{\int^\infty_0 p(s)^2\, ds + \int^\infty_0 (p(s+\delta) - p(s))^2\, ds}
$$

Taking the limit as $\delta \to 0$ gives, after applying l'Hospital's rule,

$$
\lim_{\delta \to 0}\ \frac{p(\delta)^2}{p(\delta)^2 + 2 \int^\infty_0 (p(s+\delta) - p(s))\, p'(s+\delta)\, ds}
$$

If $p(0) \neq 0$, this limit is given by

$$
\frac{p(0)^2}{p(0)^2}\ = 1.
$$

Upon noting that $p(0) = 0$ if and only if $x(t)$ is mean square differentiable,
we have the desired result that if $x(t)$ is not mean square differentiable,
then $x(t)$ is locally unpredictable.

Using the preceding theorem and our formula {eq}`eq-12-gen` for geometric distributed leads, it is straightforward to
establish that if $x(t)$ is a covariance stationary stochastic process with Wold
representation

$$
x(t) = \int^\infty_0 p(s)\, w(t-s)\, ds = \tilde P (D)\, w(t),
$$ (eq-13-wold)

then for any $\rho < 0$, the geometric sum of future expected $x$'s,

$$
\begin{aligned}
x(t)\ast &\equiv E_t \int^\infty_0 e^{\rho s}\, x(t+s)\, ds = \left[\frac{- \tilde P (D) + \tilde P(-\rho)}{D+\rho}\right]\, w(t) \\
&\equiv \tilde G(D)\, w(t) \equiv \int^\infty_0 g(s)\, w(t-s)\, ds
\end{aligned}
$$

is locally unpredictable, even if $x(t)$ is mean square differentiable. To show
this, we note that

$$
\begin{aligned}
g(0) &= \lim_{s\to \infty}\, s \tilde G(s) \\
&= \lim_{s\to \infty}\, \left[ \frac{s\tilde P(s)}{s+\rho} + \frac{s \tilde P(-\rho)}{s+\rho} \right] \\
&= \tilde P (-\rho) \neq 0.
\end{aligned}
$$

(Here we are using that $\lim_{s\to \infty} s\tilde P (s) = 0$ by the assumption
of mean square differentiability of $x(t)$.) We know that $\tilde P (-\rho) \neq 0$
because $\tilde P(s)$ has no zeroes in the right half plane, by the assumption
that {eq}`eq-13-wold` is a Wold representation.
