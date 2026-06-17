# 9. Characterizations of Mean Square Differentiability and Mean Square Continuity

The following theorem is useful.

**Theorem.** Let $x(t)$ be a linearly indeterministic covariance stationary process with moving average representation

$$
x(t) = \int_{0}^{\infty} p(\tau) w(t-\tau) \, d\tau
$$

or

$$
x(t) = \tilde P(D) w(t).
$$

Then the mean square derivative $Dx(t)$ has representation

```{math}
:label: eq-9-1
Dx(t) = \int_{0}^{\infty} p'(\tau) w(t-\tau) \, d\tau + p(0)w(t) \qquad (+)
```

or

$$
Dx(t) = Dp(t) w(t) + p(0) w(t).
$$

It follows from {eq}`eq-9-1` that the mean square derivative exists as an ordinary stochastic process only if, (a) $p(0) = 0$, and (b) $\int_{0}^{\infty} |Dp(s)| \, ds < +\infty$.

Informal proof of $(+)$:

Applying Leibniz's rule to differentiate with respect to $t$ the expression

$$
x(t) = \int_{0}^{\infty} p(t-s) w(s) \, ds.
$$

We can establish the preceding theorem rigorously as follows. From $x(t) = \int_{-\infty}^{t} p(t-s) w(s) \, ds$, we have that for $\tau > 0$,

$$
Ex(t) x(t-\tau) = \int_{-\infty}^{t} p(t-s) p(t-\tau-s) \, ds.
$$

As established above, in order for $x(t)$ to be mean square differentiable, it is necessary that $\partial R(t, \, t-\tau)/ \partial t \, \partial \tau$ exist, and that it equals zero for $\tau = 0$.

(Details to be filled in)

The generalization of this theorem to higher order derivatives follows by repeated differentiation of (—):

**Theorem.** $D^n x(t)$ exists as an ordinary stochastic process only if $p(0) = Dp(0) = \ldots = D^{n-1} p(0) = 0$ and

$$
\int_{0}^{\infty} |D^j p(s)|^2 \, ds = 0
$$

for $j = 0,\ 1,\, \ldots\, n$.

We also have:

**Theorem.** If $x(t)$ is mean square differentiable and linearly indeterministic with $x(t) = \tilde P(D) w(t)$ being its Wold representation, then

$$
Dx(t) = D \tilde P(D) w(t)
$$

or

$$
Dx(t) = \int_{0}^{\infty} p'(\tau) w(t-\tau) \, d\tau
$$

is a Wold representation for $Dx(t)$.

**Proof.** Note that if $\tilde P(s)$ has no zeroes in the right half of the complex plane, then neither does $s \tilde P(s)$. This is sufficient for the linear space spanned by $[Dx(v),\, v \leq t]$ to equal the linear space spanned by $[w(v),\, v \leq t]$.

The following proposition characterizes mean square continuity.

**Theorem.** Let $x(t)$ be a linearly indeterministic covariance stationary stochastic process with Wold moving average representation $x(t) = \tilde P(D) w(t) = \int_{0}^{\infty} p(\tau) w(t-\tau) \, d\tau$. Then $x(t)$ is mean square continuous if and only if $p(\tau)$ is continuous almost everywhere.

**Proof.** We must verify that $\lim_{\epsilon \to 0} E\, [(x(t+\epsilon) - x(t))^2] = 0$

This limit equals

$$
\begin{aligned}
\lim_{\epsilon \to 0}\ &\left[ \int_{-\epsilon}^{\infty} (p(s+\epsilon) - p(s)) w(t-s) \, ds\right]^2 \\
&= \lim_{\epsilon \to 0} \int_{-\epsilon}^{\infty} |p(s+\epsilon) - p(s)|^2 \, ds \\
&= 0
\end{aligned}
$$

if and only if $p(s)$ is continuous almost everywhere.

Together with results (—), the following characterization can provide a useful way of testing for mean square differentiability of various orders.

**Criterion (Initial value theorem).** Let $\tilde P(s) = \int_{0}^{\infty} e^{-st} p(t) \, dt$ be the Laplace transform of $p(t)$. Then

$$
p(0) = \lim_{s \to \infty}\, s \tilde P(s).
$$

**Proof.** (For case where $p'(t)$ exists)

$$
\int_{0}^{\infty} p'(t) e^{-st} \, dt = s \tilde P(s) - p(0)
$$

As $s \to \infty$, the integral on the left approaches zero.
