# 14. Examples of Nonstationary Processes

We briefly consider an example of a nonstationary process, for which many of the results above can still be used. We consider a process $x(t)$ that is described by

$$
x(t) = \int^\infty_0 p(\tau)\, w(t-\tau)\, d\tau
$$

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

It follows from (—) that $x(t)$ can also be expressed as

$$
x(t) = \int^t_{-T} p(t-s)\, w(s)\, ds.
$$

It then follows that for $\tau > 0$

$$
E x(t)\, x(t-\tau) = \int^{t-\tau}_{-T} p(t-s)\, p(t-\tau-s)\, ds.
$$

It can be verified from (—), using the criterion in theorem (—), that $x(t)$ is mean square differentiable if (a) $p(0) = 0$, and (b) $p(s)$ is twice differentiable. (Put proof in footnote.)

As an example, we take a nonstationary process governed by

$$
x(t) = \frac{\beta + D}{D^2}\, w(t) \qquad \beta > 0
$$

or

$$
x(t) = \int^{t+T}_0 (1 + \beta s)\, w(t-s)\, ds
$$

where

$$
1 + \beta t \leftrightarrow \frac{\beta + iw}{(iw)^2}\,.
$$

Applying our prediction formula to (—), have that $\hat E_t\, x(t+v)$ is given by

$$
\hat E_t x(t+v) = \int^{t+T}_0 \left(1 + \beta (s+v)\right) w(t-s)\, ds, \qquad v > 0
$$

Since for $E_t x(t+v)$, the "kernel" $p(s) = 1 + \beta (s+v)$ does not satisfy $p(0) = 0$, the process is not mean square differentiable (with respect to $t$). However, $d/dv\, E_t x(t+v)$ *does* exist as a mean square derivative. In particular, note that

$$
\begin{aligned}
\lim_{\epsilon \to 0}\ E\ &\left[ \frac{\int^\infty_0 (1 + \beta (\tau + v + \epsilon)) w(t-\tau)\, d\tau - \int^\infty_0 (1 + \beta (\tau + v)) w(t-\tau)\, d\tau}{\epsilon} \right. \\
&\left. - \beta \int^\infty_0 w(t-\tau)\, d\tau \right]^2 = 0.
\end{aligned}
$$

This shows that

```{math}
:label: eq-14-1
\frac{d}{dv}\ E_t\, x(t+v) = \beta \int^\infty_0 w(t-\tau)\, d\tau,
```

which is the expression for $d/dv\, E_t x(t+v)$ that is obtained by differentiating the right side (—) formally.

Equation {eq}`eq-14-1` can be expressed as

$$
\frac{d}{dv}\, E_t\, x(t+v) = \frac{\beta}{D}\, w(t).
$$

Solving (—) formally for $w(t)$, and using the result in (—) gives

$$
\begin{aligned}
\frac{d}{dv}\ E_t\, x(t+v) &=\ \frac{\beta}{D}\ \frac{D^2}{\beta + D}\ x(t) \\
&= \frac{\beta}{\beta + D}\ Dx(t)
\end{aligned}
$$

```{math}
:label: eq-14-2
\frac{d}{dv}\ E_t\, x(t+v) = \beta \int^\infty_0 e^{-\beta s}\, Dx(t-s)\, ds
```

or

$$
\frac{d}{dv}\ E_t\, x(t+v) = \beta \int^{t+T}_0 e^{-\beta s}\, Dx(t-s)\, ds
$$

which formally expresses the mean square derivative $d/dv\ E_t x(t+v)$ as a geometric distributed lag of the "derivative" of $x(t)$, it being understood that the derivative of $x(t)$ exists only as a generalized stochastic process. Equation {eq}`eq-14-2` is a version of Cagan's adaptive expectations scheme in continuous time. Notice that the expected rate of change $d/dv\ E_t x(t+v)$ is independent of the value of $v$ at which it is evaluated. This is the counterpart in continuous time of the adaptive expectations scheme that was noted by Muth (—) in discrete time, namely, that the optimal forecast is independent of horizon. This is a special property of the stochastic process (—) for which adaptive expectations are optimal or rational.
