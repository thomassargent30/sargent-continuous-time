# 2. Mean Square Continuity and Differentiability of a Stochastic Process

```{eval-rst}
.. index::
   single: mean square continuity
   single: mean square differentiability
   single: mean square derivative
   single: Cauchy criterion
   single: analytic function; autocovariance
   pair: smoothness; predictability
```

Having characterized a process through its autocovariance function $R(\tau)$ in
{doc}`01_covariance_stationary_processes`, we now ask how the smoothness of $R$ at the origin
translates into the continuity and differentiability of $x(t)$ itself. We
begin by defining two notions of stochastic continuity. First, a stronger notion of
continuity than we shall need is "continuity almost everywhere."

```{eval-rst}
.. index::
   single: continuity; almost sure
   single: sample path; continuity of
```

**Definition 3.** The stochastic process $x_t = x(t,w)$ is said to be continuous almost
everywhere at the point $t$ if

$$
\text{Prob}\, \left\{w: \qquad \lim_{\varepsilon \to 0}\, x(t + \varepsilon, w) = x(t, w) \right\} = 1.
$$

Definition 3 states that the realizations of $x(t,w)$ that are continuous at $t$ have
probability 1.

The weaker concept of continuity that we shall use is "mean square continuity."

```{eval-rst}
.. index::
   single: mean square continuity; definition
   single: mean square convergence
```

**Definition 4.** A stochastic process $x_t = x(t,w)$ is said to be *mean square
continuous* at the point $t$ if

$$
\lim_{\varepsilon \to 0} \int \big(x(t + \varepsilon, w) - x(t, w)\big)^2 dP(w) = 0.
$$

We immediately have the following theorem:

**Theorem 1.** The process $x(t,w)$ is mean square continuous at $t$ if and only if
$R(t_1, t_2)$ is continuous in $t_1$ and $t_2$ at $t_1 = t_2 = t$.

**Proof.** Note that

$$
E\big(x(t+\tau) - x(t)\big)^2 = E x(t + \tau)^2 + E x(t)^2 - 2 E x(t + \tau) x(t)
$$

or

$$
E\big(x(t + \tau) - x(t)\big)^2 = R(t + \tau, t + \tau) + R(t,t) - 2R(t + \tau, t).
$$

Taking limits of both sides as $\tau \to 0$ proves the theorem.

Next we have the definition:

**Definition 5.** A stochastic process $x_t = x(t,w)$ is said to be *mean square
continuous* if it is mean square continuous at each point $t \in T$.

We immediately have:

**Theorem 2.** If a stochastic process $x(t,w)$ is mean square continuous, then
$E x(t) = \mu(t)$ is continuous.

**Proof.** For any random variable $z$, $E z^2 =$ variance $z + (Ez)^2 \geq (Ez)^2$. It
follows that

$$
E\,\big\{x(t + \tau, w) - x(t, w)\big\}^2 \geq \big\{E\,\big(x(t + \tau, w) - x(t, w)\big)\big\}^2.
$$

Taking limits as $\tau \to 0$ proves the theorem.

We also have the following theorem:

```{eval-rst}
.. index::
   single: mean square continuity; criterion for
   single: autocovariance function; continuity at the origin
```

**Theorem 3.** A covariance stationary stochastic process $x(t,w)$ is mean square
continuous if and only if its autocorrelation function $R(\tau) = E x(t) x(t-\tau)$ is
continuous for $\tau = 0$.

**Proof.** Theorem 3 is implied by the proof of theorem 1.

Next we turn to a concept of stochastic differentiation. It would be convenient to have a
concept of differentiation that rationalized the following interchange of orders of
integration and differentiation:

$$
\begin{aligned}
\frac{\partial^2}{\partial t_1\,\partial t_2}\ \ R(t_1, t_2) &= \frac{\partial^2}{\partial t_1\,\partial t_2}\ \int x(t_1, w)\, x(t_2, w)\, dP(w) \\
&= \int\ \frac{d}{dt}\ x(t_1, w)\ \frac{d}{dt}\ x(t_2, w)\, dP(w) \\
&= E\,\big\{\tfrac{d}{dt}\ x(t_1)\ \tfrac{d}{dt}\ x(t_2)\big\}.
\end{aligned}
$$

or

$$
\frac{\partial^2}{\partial t_1\,\partial t_2}\ \ R(t_1, t_2) = E\,\big\{\tfrac{d}{dt}\ x(t_1)\ \tfrac{d}{dt}\ x(t_2)\big\}.
$$

We shall see that the concept of mean square differentiability has this property.

```{eval-rst}
.. index::
   single: mean square derivative; definition
   single: mean square differentiability; definition
```

**Definition 6.** A stochastic process $x(t, w)$ is said to have a *mean square derivative*
$x'(t)$ at $t$ if there is a random variable $x'(t)$ such that

$$
\lim_{\epsilon \to 0}\ E\, \biggl\{ \Big( \frac{x(t + \epsilon) - x(t)}{\epsilon}\ - x'(t) \Big)^2\biggr\} = 0.
$$

If $x(t, w)$ has a mean square derivative, we say that $x(t, w)$ is *mean square
differentiable*.

We have the following Cauchy criterion for the mean square differentiability of $x(t,w)$
at $t$:

**Cauchy criterion.** A process $x(t,w)$ is mean square differentiable at $t$ if and only
if

$$
\lim_{\epsilon_1 \to 0,\ \epsilon_2 \to 0}\ E\biggl\{ \Big[ \frac{x(t + \epsilon_1) - x(t)}{\epsilon_1}\ - \ \frac{x(t + \epsilon_2) - x(t)}{\epsilon_2}\Big]^2\biggr\} = 0.
$$

We have the following theorem:

```{eval-rst}
.. index::
   single: autocovariance function; differentiability of
   single: mean square differentiability; criterion in terms of the autocovariance
   pair: smoothness; autocovariance
```

**Theorem 4.** Let $x(t,w)$ be a covariance stationary process with autocorrelation
function $R(\tau)$. If $x'(t)$ exists in the mean square sense, then $R''(0)$ exists.

**Proof.** First notice that

```{math}
:label: eq-2-4
\begin{aligned}
E\, \biggl[ \Big( \frac{x(t + \epsilon) - x(t)}{\epsilon}\Big)^2 \biggr]\ &= \frac{2R(0) - 2R(\epsilon)}{\epsilon^2} \\
&=\ - \ \Big[ \frac{R(\epsilon) - 2R(0) + R(-\epsilon)}{\epsilon^2}\Big].
\end{aligned}
```

Taking limits as $\epsilon \to 0$, it follows that

```{math}
:label: eq-2-5
\lim_{\epsilon \to 0}\ E\, \Big( \frac{x(t + \epsilon) - x(t)}{\epsilon}\Big)^2\ = \ -R''(0),
```

for the limit on the left exists by the assumption of mean square differentiability.

We also have a converse of the above theorem:

**Theorem 5.** Let $x(t,w)$ be a covariance stationary process with autocorrelation
function $R(\tau)$. If $R''(0)$ exists, then $x(t,w)$ is mean square differentiable.

Theorems 4 and 5 together say that, for a covariance stationary process, mean square
differentiability is *equivalent* to the existence of $R''(0)$. Nothing more than the second
derivative *at the origin* is required; the proof below uses only that.

**Proof.** To apply the Cauchy criterion, we shall need to evaluate

$$
\begin{aligned}
E\, &\big\{ \frac{x(t + \epsilon_1) - x(t)}{\epsilon_1}\ \cdot\ \frac{x(t + \epsilon_2) - x(t)}{\epsilon_2}\big\} \\
&= \frac{R(\epsilon_1 - \epsilon_2) - R(\epsilon_1) - R(-\epsilon_2) + R(0)}{\epsilon_1\, \epsilon_2} \\
&= \frac{ \frac{R(\epsilon_1 - \epsilon_2) - R(\epsilon_1)}{\epsilon_2} \ + \ \frac{R(0) - R(-\epsilon_2)}{\epsilon_2} }{\epsilon_1}.
\end{aligned}
$$

Taking limits first as $\epsilon_2 \to 0$, then as $\epsilon_1 \to 0$ gives

```{math}
:label: eq-2-6
\begin{aligned}
\lim_{\epsilon_1,\, \epsilon_2 \to 0}\ &E\big\{ \frac{x(t + \epsilon_1) - x(t)}{\epsilon_1}\ \cdot\ \frac{x(t + \epsilon_2) - x(t)}{\epsilon_2}\big\} = \\
\lim_{\epsilon_1 \to 0}\ & -\ \frac{R'(\epsilon_1) - R'(0)}{\epsilon_1}\ =\ - R''(0).
\end{aligned}
```

Using {eq}`eq-2-5` and {eq}`eq-2-6`, we find that

$$
\begin{aligned}
\lim_{\epsilon_1,\, \epsilon_2 \to 0}\ &E\biggl\{ \Big[ \frac{x(t + \epsilon_1) - x(t)}{\epsilon_1}\ - \ \frac{x(t + \epsilon_2) - x(t)}{\epsilon_2}\Big]^2\biggr\} \\
&= \ - 2R''(0) + 2R''(0) = 0.
\end{aligned}
$$

Thus, if $R''(0)$ exists, then $x(t)$ is mean square differentiable.

Since $R(\tau) = R(-\tau)$, it follows from the fact that $R''(0)$ exists for a mean square
differentiable process that $R'(0) = 0$.

For a nonstationary stochastic process, the counterpart of the two preceding theorems is
the following theorem, which we present without proof

**Theorem 6.** A nonstationary stochastic process $x(t, w)$ is mean square differentiable
if

$$
\frac{\partial^2 R(t_1, t_2)}{\partial t_1\,\partial t_2}
$$

exists for $t_1 = t_2$.

For a mean square differentiable process $x(t, w)$, the following interchange of order of
integration (expectation) and differentiation is appropriate:

$$
\begin{aligned}
E x'(t) &= E\ \lim_{\epsilon \to 0}\ \frac{x(t + \epsilon) - x(t)}{\epsilon} \\
&= \ \lim_{\epsilon \to 0}\ \frac{E x(t + \epsilon) - E x(t)}{\epsilon} \ = \ \frac{d}{dt}\ E x(t).
\end{aligned}
$$

Thus,

$$
E x'(t) = \frac{d}{dt}\ \mu(t).
$$

Thus, we obtain the mean function of $x'(t)$ by once differentiating the mean function
$\mu(t)$ of $x(t)$.

We can derive the autocorrelation function of $x'(t)$ by twice differentiating the
autocorrelation function of $x(t)$. To establish this, we need some additional notation. We
define the autocorrelation function of $x(t)$ as

$$
E x(t_1)\, x(t_2) = R_{xx}(t_1, t_2).
$$

We also define

$$
E x(t_1)\, x'(t_2) = R_{xx'}(t_1, t_2),
$$

and

$$
E x'(t_1)\, x'(t_2) = R_{x'x'}(t_1, t_2).
$$

First, we shall show that

```{math}
:label: eq-2-7
\frac{\partial R_{xx}(t_1, t_2)}{\partial t_2} = R_{xx'}(t_1, t_2).
```

To show this, note that

$$
\begin{aligned}
E &\ \big\{x(t_1) \frac{x(t_2 + \epsilon) - x(t_2)}{\epsilon} \big\} \\
&= \frac{R_{xx}(t_1, t_2 + \epsilon) - R_{xx}(t_1, t_2)}{\epsilon}
\end{aligned}
$$

Taking limits as $\epsilon \to 0$ gives the desired results {eq}`eq-2-7`. Similarly,

$$
R_{x'x'}(t_1, t_2) = \lim\, E\, \big\{ \frac{x(t_1 + \epsilon) - x(t_1)}{\epsilon}\ x'(t_2) \big\}\ = \frac{\partial R_{xx'}(t_1, t_2)}{\partial t_1}
$$

Therefore, we have proved the desired result, which we state in the following theorem:

**Theorem 7.** If the stochastic process $x(t,w)$ has mean square derivative $x'(t)$, then
the autocorrelation function of $x'(t)$ is given by

$$
R_{x'x'}(t_1, t_2) = \frac{\partial^2 R(t_1, t_2)}{\partial t_1\,\partial t_2}.
$$

For covariance stationary processes, we have the immediate corollary.

**Corollary 7C.** If the covariance stationary stochastic process $x(t, w)$ is mean square
differentiable, the autocorrelation function of $x'(t)$ is given by

$$
R_{x'}(\tau) = - R''(\tau)
$$

where $R(\tau)$ is the autocorrelation function of $x(t)$.

By successively applying the preceding reasoning to $x'(t)$ and each of its derivatives in
an evident way, we can prove the following theorem:

**Theorem 8.** Let $x(t, w)$ be a stochastic process with autocorrelation function
$R(t_1, t_2)$. The process is $n$ times mean square differentiable if

$$
\frac{\partial^{2n} R(t_1, t_2)}{\partial t_1^n \partial t_2^n}
$$

exists. The autocorrelation of the $n^{th}$ mean square derivative process $x^{(n)}(t)$
equals $\partial^{2n} R(t_1, t_2)/ \partial t_1^n \partial t_2^n$.

As an example, consider a stochastic process for which

$$
R(\tau) = e^{\lambda |\tau|},\ \lambda < 0
$$

This process is mean square continuous, but not mean square differentiable. (Why?) Next,
consider a process for which

$$
R(\tau) = k_1 e^{\lambda_1 |\tau|} + k_2 e^{\lambda_2 |\tau|};\ \lambda_1,\, \lambda_2\, < 0.
$$

The process is mean square continuous, but is mean square differentiable only if
$\lambda_1\, k_1 + \lambda_2\, k_2 = 0$. (Why?)

The following construction demonstrates a link between the existence of mean square
derivatives of arbitrarily high orders, and the predictability of a series. Let $x(t, w)$
be a covariance stationary process with autocorrelation function $R(\tau)$. We say that
$R(\tau)$ is *analytic* if its derivatives of all orders exist for all $\tau$, and if
$R(\tau)$ has the Taylor (Maclaurin) series representation:

$$
R(\tau) = \sum_{n=0}^{\infty} \ R^{(n)}\ (0)\ \frac{\tau^n}{n!}
$$ (eq-2-analytic)

We note that if $R(\tau)$ is analytic, then for all integer $n > 0$, the $n^{th}$ mean
square derivative $x^{(n)}\,(t)$ exists; we can now state the following theorem.

```{eval-rst}
.. index::
   single: Taylor series; of an autocovariance
```

**Theorem 9.** Let $x(t, w)$ be a covariance stationary stochastic process with analytic
autocorrelation function $R(\tau)$. Then $x(t)$ can be expanded in a Taylor series, i.e.,

$$
x(t + \tau) = \sum_{n=0}^{\infty} x^{(n)}(t) \frac{\tau^n}{n!},\ \text{ for all } \tau > 0.
$$

**Proof.** We have to show that

$$
\begin{aligned}
E\, \big\{ x(t + \tau) - \hat x(t + \tau) \big\}^2 &= E\, \big\{ [ x(t + \tau) - \hat x(t + \tau)] x(t + \tau) \big\} \\
&- E\, \big\{ [x(t + \tau) - \hat x(t + \tau) ] \hat x(t + \tau) \big\} = 0,
\end{aligned}
$$ (eq-2-plus)

where

$$
\hat x(t + \tau) = \sum_{n=0}^{\infty} x^{(n)}(t)\ \frac{\tau^n}{n!}.
$$ (eq-2-xhat)

From the analytic nature of $R(\tau)$ it follows from {eq}`eq-2-analytic` that

$$
R^{(m)}(\tau) = \sum_{n=m}^{\infty} R^{(n)}(0)\ \frac{\tau^{n-m}}{(n-m)!}\, ,\ \text{ for } m \geq 1.
$$ (eq-2-Rm)

It also follows from a Taylor series of $R(\tau + \lambda)$ about $\lambda = 0$ that

$$
R(0) = \sum_{n=0}^{\infty} R^{(n)}(\tau)\ \frac{(-\tau)^n}{n!}.
$$ (eq-2-R0)

Substituting the right side of {eq}`eq-2-xhat` into {eq}`eq-2-plus`, noting by the
reasoning that led to Theorem 8 that $E x^{(n)}(t) x^{(m)}(t - \tau) = (-1)^m R^{(n + m)}(\tau)$,
and using {eq}`eq-2-Rm` and {eq}`eq-2-R0` to evaluate the two terms in braces in
{eq}`eq-2-plus` gives the desired results.

The preceding states that if $R(\tau)$ is analytic, then the stochastic process $x(t)$ is
differentiable an arbitrarily large number of times, and that $x(t)$ is perfectly
forecastable arbitrarily far into the future from values of $x$ and its mean square
derivatives at time $t$. In our work, we shall usually want to deal with stochastic
processes that are only imperfectly forecastable from knowledge of the past. This means
that we shall usually deal with processes for which the autocorrelation function $R(\tau)$
is not analytic. The differentiability criteria of this chapter take a sharper, more usable
form once a process is written in its moving-average (Wold) representation:
{doc}`09_characterizations_ms_differentiability` shows that $x(t)$ is mean square
differentiable precisely when the moving-average kernel satisfies $p(0) = 0$, and
{doc}`13_locally_unpredictable` ties the failure of that condition to local unpredictability.
