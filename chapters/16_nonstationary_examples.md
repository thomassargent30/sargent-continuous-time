# 16. Examples of Nonstationary Processes

We briefly consider an example of a nonstationary process, for which many of the results above can still be used. We consider a process $x(t)$ that is described by

$$
x(t) = \int^\infty_0 p(\tau)\, w(t-\tau)\, d\tau
$$ (eq-16-ma)

where $w(t)$ is a nonstationary white noise with

$$
E w(t)\, w(t-\tau) = \begin{cases}
\delta(\tau) & \text{for } t \geq -T \\
0 & t < -T,\ T > 0.
\end{cases}
$$

This specification states that $w(t) \equiv 0$ for $t < -T$, so that the system can be imagined to start up at $t = -T$. We also assume that

$$
\int^\infty_0 p(\tau)^2\, d\tau = +\infty.
$$

It follows from {eq}`eq-16-ma` that $x(t)$ can also be expressed as

$$
x(t) = \int^t_{-T} p(t-s)\, w(s)\, ds.
$$

It then follows that for $\tau > 0$

$$
E x(t)\, x(t-\tau) = \int^{t-\tau}_{-T} p(t-s)\, p(t-\tau-s)\, ds.
$$ (eq-16-acov)

It can be verified from {eq}`eq-16-acov`, using the mean square differentiability criterion of {doc}`02_mean_square_continuity_differentiability` — the same $p(0) = 0$ condition developed in {doc}`09_characterizations_ms_differentiability` and {doc}`15_locally_unpredictable` — that $x(t)$ is mean square differentiable if (a) $p(0) = 0$, and (b) $p(s)$ is twice differentiable.[^fn16-msd]

As an example, we take a nonstationary process governed by

$$
x(t) = \frac{\beta + D}{D^2}\, w(t) \qquad \beta > 0
$$ (eq-16-proc)

or

$$
x(t) = \int^{t+T}_0 (1 + \beta s)\, w(t-s)\, ds
$$ (eq-16-procint)

where

$$
1 + \beta t \leftrightarrow \frac{\beta + iw}{(iw)^2}\,.
$$

The operator $\dfrac{\beta + D}{D^2}$ has the rational form $\psi(D)/\theta(D)$ of the linear
stochastic differential equations of {doc}`11_linear_sde`, but with the denominator
$\theta(D) = D^2$ placing its double root at the origin rather than in the open left half plane.
It is exactly this departure from the stability condition $\operatorname{re}(\lambda_j) < 0$
assumed there — the characteristic roots lying on the imaginary axis — that makes $x(t)$
nonstationary, with the integrated kernel $p(\tau) = 1 + \beta\tau$ and
$\int_0^\infty p(\tau)^2\, d\tau = \infty$. The operator calculus and prediction formulas of
{doc}`11_linear_sde` and {doc}`12_prediction` nonetheless continue to apply, as the next steps
show.

Applying the Wiener–Kolmogorov prediction formula of {doc}`12_prediction` to {eq}`eq-16-procint` — here extended to this nonstationary process — we have that $\hat E_t\, x(t+v)$ is given by

$$
\hat E_t x(t+v) = \int^{t+T}_0 \left(1 + \beta (s+v)\right) w(t-s)\, ds, \qquad v > 0
$$ (eq-16-forecast)

Since for $E_t x(t+v)$, the "kernel" $p(s) = 1 + \beta (s+v)$ does not satisfy $p(0) = 0$, the process is not mean square differentiable (with respect to $t$); like the processes of {doc}`15_locally_unpredictable`, it is locally unpredictable — its increments behave locally like those of a martingale. However, $d/dv\, E_t x(t+v)$ *does* exist as a mean square derivative. In particular, note that

$$
\begin{aligned}
\lim_{\epsilon \to 0}\ E\ &\left[ \frac{\int^\infty_0 (1 + \beta (\tau + v + \epsilon)) w(t-\tau)\, d\tau - \int^\infty_0 (1 + \beta (\tau + v)) w(t-\tau)\, d\tau}{\epsilon} \right. \\
&\left. - \beta \int^\infty_0 w(t-\tau)\, d\tau \right]^2 = 0.
\end{aligned}
$$

This shows that

```{math}
:label: eq-16-1
\frac{d}{dv}\ E_t\, x(t+v) = \beta \int^{t+T}_0 w(t-\tau)\, d\tau,
```

which is the expression for $d/dv\, E_t x(t+v)$ that is obtained by differentiating the right side of {eq}`eq-16-forecast` formally.

Equation {eq}`eq-16-1` can be expressed as

$$
\frac{d}{dv}\, E_t\, x(t+v) = \frac{\beta}{D}\, w(t).
$$ (eq-16-dvD)

Solving {eq}`eq-16-proc` formally for $w(t)$, and using the result in {eq}`eq-16-dvD` gives

$$
\begin{aligned}
\frac{d}{dv}\ E_t\, x(t+v) &=\ \frac{\beta}{D}\ \frac{D^2}{\beta + D}\ x(t) \\
&= \frac{\beta}{\beta + D}\ Dx(t)
\end{aligned}
$$

```{math}
:label: eq-16-2
\frac{d}{dv}\ E_t\, x(t+v) = \beta \int^\infty_0 e^{-\beta s}\, Dx(t-s)\, ds
```

or

$$
\frac{d}{dv}\ E_t\, x(t+v) = \beta \int^{t+T}_0 e^{-\beta s}\, Dx(t-s)\, ds
$$

which formally expresses the mean square derivative $d/dv\ E_t x(t+v)$ as a geometric distributed lag of the "derivative" of $x(t)$, it being understood that the derivative of $x(t)$ exists only as a generalized stochastic process. Equation {eq}`eq-16-2` is a version of Cagan's adaptive expectations scheme in continuous time. Notice that the expected rate of change $d/dv\ E_t x(t+v)$ is independent of the value of $v$ at which it is evaluated. This is the counterpart in continuous time of the adaptive expectations scheme that was noted by Muth (1960) in discrete time, namely, that the optimal forecast is independent of horizon. This is a special property of the stochastic process {eq}`eq-16-proc` for which adaptive expectations are optimal or rational. {doc}`21_aggregation_inverse_optimal_predictor` develops exactly this example into a full bivariate money-creation/inflation model — solving the "inverse optimal predictor" problem for Cagan's scheme — and asks what sampling and aggregation over time do to it.

## Notes

[^fn16-msd]: Write $x(t) = \int_{-T}^{t} p(t-s)\, w(s)\, ds$ and split the increment over
    $[t, t+\epsilon]$ into the revision of the existing kernel (over $s \in [-T, t]$) and the
    fresh innovations (over $s \in [t, t+\epsilon]$). Since $E\, w(s) w(s') = \delta(s-s')$ and
    the two pieces occupy disjoint intervals of $s$, they are uncorrelated, so
    $$E\big[(x(t+\epsilon) - x(t))^2\big] = \int_{0}^{t+T} \big[p(\tau+\epsilon) - p(\tau)\big]^2\, d\tau + \int_{0}^{\epsilon} p(\tau)^2\, d\tau .$$
    Dividing by $\epsilon^2$, the first term tends to $\int_0^{t+T} p'(\tau)^2\, d\tau$ when $p$
    is twice continuously differentiable, while the second behaves like $p(0)^2/\epsilon$, which
    diverges unless $p(0) = 0$. Hence the difference quotient $[x(t+\epsilon) - x(t)]/\epsilon$
    converges in mean square — making $x$ mean square differentiable, with
    $Dx(t) = \int_{-T}^{t} p'(t-s)\, w(s)\, ds$ — if and only if $p(0) = 0$, the smoothness of
    $p$ supplying convergence of the first term. This is the nonstationary counterpart of the
    boundary term $p(0)\, w(t)$ in {eq}`eq-9-1` of {doc}`09_characterizations_ms_differentiability`.
