# 1. Covariance Stationary Stochastic Processes

A stochastic process is a collection of random variables indexed by a variable $t$
belonging to a set $T$. For a continuous time stochastic process $t \in \mathbb{R}$, or
some interval in $\mathbb{R}$, where $\mathbb{R}$ is the set of real numbers. It is
imagined that there is an underlying probability space, which is characterized by a
3-tuple $(\Omega, \mathcal{F}, P)$, where $\Omega$ is a collection of points,
$\mathcal{F}$ a collection of sets in $\Omega$, and $P$ a probability measure. The sets
in $\mathcal{F}$ are collections of points for which the probability can be computed with
the probability measure $P$. For example, for a set $A \in \mathcal{F}$, $P(A)$ is a real
number in the interval $[0, 1]$. A set $A \in \mathcal{F}$ is called an event. A random
variable is a "measurable function" of $w \in \Omega$. This means that a random variable
$x(w)$ is a function from $\Omega$ to the real line for which the probability is well
defined for all events of the form $A = [w : x(w) \in [a, b]]$ for all intervals with
$a \leq b$.

A stochastic process $x_t$ is a collection of such measurable functions $x_t = x(t,w)$,
for $t \in T$. Hereafter, we set $T = \mathbb{R}$, so that we are studying a continuous
time stochastic process. We begin by assuming that $x(t,w)$ is a scalar random variable,
i.e., $x(t,w)$ is a "univariate stochastic process." In this construction, for each
"drawing" of $w \in \Omega$, one draws an entire function of time $x(t,w)$. Each such
drawing of a function of time $x_t = x(t,w)$ is called a "realization of the stochastic
process," or a "sample path of the stochastic process." We might desire to characterize
the average behavior of these functions $x_t = x(t,w)$ both across realizations or $w$'s,
and across time for a *given* single realization or drawing of $w$. For example,
corresponding to the "mean", we could compute $\int x(t,w)\,dP(w)$, which is the average
across realizations, at a given time $t$; alternatively, we could compute
$\lim_{T \to \infty} (2T)^{-1} \int^T_{-T} x(t,w)\,dt$, the average across $t$ during a
single realization $w$. In these pages we spend most of our time characterizing the
behavior of the stochastic process $x_t = x(t,w)$ across realizations $w$. Later, we
shall see that under some regularity conditions on the stochastic process $x(t,w)$, the
behavior of averages across $w$'s is closely reflected in behavior of averages across
time for a single $w$. These regularity conditions are conditions for "ergodicity," and
must be imposed to acquire a practical theory of estimation. For now, it is important to
realize that averages across time within a single realization are logically distinct from
averages over realizations at a point in time.

We define the following three functions which measure population averages across
realizations:

```{math}
:label: eq-1
\begin{aligned}
\mu(t) &= \int x(t,w)\,dP(w) \qquad t \in \mathbb{R} \\
&= E x(t,w) \\[4pt]
R(t_1,t_2) &= \int x(t_1,w)\,x(t_2,w)\,dP(w); \ t_1, t_2 \in \mathbb{R} \\
&= E x(t_1,w)\,x(t_2,w) \\
C(t_1,t_2) &= \int \big(x(t_1,w) - \mu(t_1)\big)\,\big(x(t_2,w) - \mu(t_2)\big)\,dP(w), \\
&= E \big(x(t_1,w) - \mu(t_1)\big)\,\big(x(t_2,w) - \mu(t_2)\big)\ \ t_1, t_2 \in \mathbb{R}
\end{aligned}
```

We assume that the integrals on the right side exist, which means that the stochastic
process possesses finite first and second moments. The function $\mu(t)$ is called the
*mean function*, while $R(t_1,t_2)$ is called the *autocorrelation function*, and
$C(t_1,t_2)$ is called the *autocovariance function*. Often we shall be studying a
stochastic process for which $\mu(t) = 0$, in which case $C(t_1,t_2) = R(t_1,t_2)$.

```{note}
**A warning about the word "autocorrelation."** Following an older usage common in the
engineering literature on which these notes draw, $R$ is here the *second moment* function
$E\,x(t_1)x(t_2)$ — it is not normalized to lie in $[-1,1]$, and it is not centred unless
$\mu = 0$. Much of modern statistics reserves "autocorrelation" for the normalized quantity
$C(\tau)/C(0)$. Since we shall work almost exclusively with zero-mean processes, where $R$ and
$C$ coincide, no harm follows; but the reader coming from a statistics text should read $R$ as
the autocovariance. {doc}`08_spectral_densities` calls the same object the *covariogram*, and
it is $R$ — not any normalized version of it — whose Fourier transform is the spectral density.
```

In the following definition, we introduce the notion of stationarity.

**Definition 1.** A stochastic process is said to be *covariance stationary* (or
*second-order stationary* or *wide-sense stationary*) if

```{math}
:label: eq-2
\begin{aligned}
\mu(t) &= \mu\ \text{ for all } t \in \mathbb{R} \\
R(t_1,t_2) &= R(t_1 + \varepsilon, t_2 + \varepsilon)\ \text{ for all } t_1, t_2 \in \mathbb{R}\ \text{ and all } \varepsilon \text{ in } \mathbb{R}
\end{aligned}
```

According to {eq}`eq-2`, the mean $\mu(t)$ is a constant function of time, $\mu$, and the
covariance between $x(t_1)$ and $x(t_2)$ equals the covariance between $x(t_1 +
\varepsilon)$ and $x(t_2 + \varepsilon)$, so that the covariance between $x(t)$'s at two
distinct dates depends only on the difference between the two dates. Setting $\varepsilon
= -t_2$, and letting $t_1 - t_2 = \tau$, we have $R(t_1, t_2) = R(t_1 - t_2, 0) = R(t_1,
t_1 - (t_1 - t_2))$ or

$$
R(t_1, t_1 - \tau) = R(\tau, 0).
$$

Adopting an inconsistent but convenient notation, for covariance stationary processes, we
shall denote

$$
E x(t)\, x(t - \tau) = R(\tau)
$$

where $R(\tau)$ is understood as $R(\tau, 0)$. Similarly, we shall also use the notation
$C(\tau, 0) = C(\tau)$.

It is possible to characterize the first and second moments of a covariance stationary
stochastic process by its mean $\mu$ and its autocovariance function $C(\tau)$. It is
useful to know the properties that a function $C(\tau)$ must have in order to be the
autocovariance function of *some* covariance stationary stochastic process. To be the
autocovariance function of a covariance stationary stochastic process, it is necessary and
sufficient that $C(\tau)$ be a positive semidefinite function. There are several
equivalent ways to express the concept of positive semidefiniteness of a function. One way
is as follows. Let $x(t)$ be a stochastic process with mean zero and autocovariance
function $C(\tau)$. Then consider forming a random variable

$$
y_h = \int^\infty_{-\infty} h(s)\, x(t - s)\, ds
$$

where $h$ is a function that is a member of the class $L^2(-\infty, \infty)$ of square
integrable functions, i.e., $\int^\infty_{-\infty} h(s)^2\, ds < +\infty$. Evidently, all
such random variables $y_h$ must have nonnegative variance. This is the requirement that

```{math}
:label: eq-3
\int^\infty_{-\infty} \int^\infty_{-\infty} h(s)\, h(\tau)\, C(\tau - s)\, ds\, d\tau \geq 0
\ \text{ for all } h(s) \in L^2(-\infty, \infty).
```

Condition {eq}`eq-3` is the condition that $C(\tau)$ is positive semidefinite. That
$C(\tau)$ be positive semidefinite is a necessary and sufficient condition for $C(\tau)$
to be an autocovariance function for some covariance stationary stochastic process.

Another way to characterize positive semidefiniteness of $C(\tau)$ is as follows. Let
$t_1, \ldots, t_n$ be any finite set of real numbers. Then use $C(\tau)$ to form the
$(n \times n)$ covariance matrix for the vector $(x(t_1) \ldots x(t_n))$. For every finite
$n$, and every collection $(t_1, \ldots, t_n)$, the covariance matrix must be positive
semidefinite if the function $C(\tau)$ is to be positive semidefinite.

Yet another way to characterize positive semidefiniteness of $C(\tau)$ can be given using
concepts to be introduced in {doc}`08_spectral_densities`. It is true that $C(\tau)$ is positive
semidefinite if and only if the function

$$
S(w) = \int^\infty_{-\infty} C(\tau)\, e^{-i w \tau}\, d\tau \geq 0\ \text{ for all } w \in (-\infty, \infty).
$$

We now define a more restrictive form of stationary:

**Definition 2.** A stochastic process is said to be *stationary* or *strictly
stationary* if for any finite $n$, any collection of dates $[t_1, \ldots, t_n]$, and any
measurable function $g(x(t_1,w), \ldots, x(t_n,w))$, it is true that

$$
\begin{aligned}
\int\ & g\big(x(t_1,w), \ldots, x(t_n,w)\big)\, dP(w) = \\
\int\ & g\big(x(t_1 + \varepsilon, w), \ldots, x(t_n + \varepsilon, w)\big)\, dP(w),
\ \text{ for every } \varepsilon \text{ in } \mathbb{R}.
\end{aligned}
$$

Notice that strict stationarity implies covariance stationarity, but not *vice versa*.
Throughout this book we work with the weaker, second-moment notion. The next chapter shows
that the smoothness of the autocovariance function $R(\tau)$ at the origin governs the
continuity and differentiability of $x(t)$ itself — the first appearance of a
smoothness-versus-predictability theme that recurs in
{doc}`09_characterizations_ms_differentiability` and {doc}`13_locally_unpredictable`.
