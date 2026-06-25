# 15. Locally Unpredictable Stochastic Processes

Lack of mean square differentiability implies that a process is "locally
unpredictable" or is "locally like a martingale." Sims used the concept of
local unpredictability in his work on asset prices, interest rates, and
consumption. The white innovation $\eta$ of the state-space innovations representation in
{doc}`13_kalman_filter_spectral_factorization` is the leading example: it is locally
unpredictable by construction. We use the following definition:

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
$Ew(t) w(t-\tau) = \delta(\tau)$. Assume that $p(\tau)$ is twice continuously
differentiable, with $p'$ and $p''$ square integrable. Then if $x(t)$ is not mean square
differentiable, $x(t)$ is locally unpredictable.

**Proof:** From the Wiener-Kolmogorov prediction formula, we have that

$$
\frac{E_t (x(t+\delta) - E_t x(t+\delta))^2}{E_t (x(t+\delta) - x(t))^2}
$$

$$
= \frac{\int^\delta_0 p(s)^2\, ds}{\int^\delta_0 p(s)^2\, ds + \int^\infty_0 (p(s+\delta) - p(s))^2\, ds}
$$

Taking the limit as $\delta \to 0$ gives, after applying l'Hospital's rule,

$$
\lim_{\delta \to 0}\ \frac{p(\delta)^2}{p(\delta)^2 + 2 \int^\infty_0 (p(s+\delta) - p(s))\, p'(s+\delta)\, ds}
$$

If $p(0) \neq 0$, this limit is given by

$$
\frac{p(0)^2}{p(0)^2}\ = 1.
$$

Upon noting that $p(0) = 0$ if and only if $x(t)$ is mean square differentiable
(by {doc}`09_characterizations_ms_differentiability`: with $p'$ square integrable, mean square
differentiability reduces to the single condition $p(0) = \lim_{s\to\infty} s\tilde P(s) = 0$),
we have the desired result that if $x(t)$ is not mean square differentiable, then $x(t)$ is
locally unpredictable.

For the rational linear stochastic differential equations of {doc}`11_linear_sde`, this
dividing line is explicit. There $x(t)$ is mean square differentiable $n - 1 - m$ times, where
$n$ and $m$ are the degrees of the denominator operator $\theta(D)$ and the numerator operator
$\psi(D)$; by the initial value theorem, $p(0) = \lim_{s\to\infty} s\tilde P(s) \neq 0$ exactly
when $m = n - 1$. The locally unpredictable members of this family are therefore precisely those
with $n - 1 - m = 0$ — a numerator just one degree below the denominator — while every smoother
member ($m < n - 1$) is locally predictable.

Using the preceding theorem and our formula {eq}`eq-12-gen` for geometric distributed leads from
{doc}`12_prediction`, it is straightforward to
establish that if $x(t)$ is a covariance stationary stochastic process with Wold
representation

$$
x(t) = \int^\infty_0 p(s)\, w(t-s)\, ds = \tilde P (D)\, w(t),
$$ (eq-15-wold)

then for any $\rho < 0$, the geometric sum of future expected $x$'s,

$$
\begin{aligned}
x(t)\ast &\equiv E_t \int^\infty_0 e^{\rho s}\, x(t+s)\, ds = \left[\frac{- \tilde P (D) + \tilde P(-\rho)}{D+\rho}\right]\, w(t) \\
&\equiv \tilde G(D)\, w(t) \equiv \int^\infty_0 g(s)\, w(t-s)\, ds
\end{aligned}
$$

is locally unpredictable, even if $x(t)$ is mean square differentiable. To show
this, we apply the initial value theorem of {doc}`09_characterizations_ms_differentiability`
to the kernel $g$:

$$
\begin{aligned}
g(0) &= \lim_{s\to \infty}\, s \tilde G(s) \\
&= \lim_{s\to \infty}\, \left[ \frac{-s\tilde P(s)}{s+\rho} + \frac{s \tilde P(-\rho)}{s+\rho} \right] \\
&= \tilde P (-\rho) \neq 0.
\end{aligned}
$$

(Here we are using that $\lim_{s\to \infty} s\tilde P (s) = 0$ by the assumption
of mean square differentiability of $x(t)$.) We know that $\tilde P (-\rho) \neq 0$
because $\tilde P(s)$ has no zeroes in the right half plane, by the assumption
that {eq}`eq-15-wold` is a Wold representation.
