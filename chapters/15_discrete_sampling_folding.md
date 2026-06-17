# 15. Discrete Sampling: The Folding Formula

Let $x(t)$ be a continuous time, covariance stationary stochastic process with
autocovariogram $R(\tau)$. Let $x_j = x(j)$ be a record of $x$ sampled at the discrete
points in time $j = 0,\ \pm\, T,\ \pm\, 2T, \ldots,$ where $T > 0$ is the sampling
interval. The autocovariogram of the discrete data can be represented as a generalized
function $R^d(\tau)$ where

$$
R^d(\tau) = \sum_{n=-\infty}^{\infty} R(nT)\, \delta(\tau - nT),
$$

so that $R^d(\tau)$ is a train of delta functions with mass $R(nT)$ at
$\tau = 0,\ \pm\, T,\ \pm\, 2T, \ldots\,$. Equation (—) can also be represented as

$$
R^d(\tau) = R(\tau)\ S_T(\tau)
$$

where

$$
S_T(\tau) = \sum_{n=-\infty}^{\infty} \delta(\tau - nT).
$$

From (—), we can express the spectral density of the discrete sampled $x_t$ as

```{math}
:label: eq-15-1
\begin{aligned}
S^d(w) &= \int_{-\infty}^{\infty} R^d(\tau)\, e^{-iw\tau}\, d\tau \\
&= \sum_{n=-\infty}^{\infty} \int_{-\infty}^{\infty} R(\tau)\, e^{-iw\tau}\, \delta(\tau - nT)\, d\tau \\
&= \sum_{n=-\infty}^{\infty} R(nT)\, e^{-iwnT}
\end{aligned}
```

which is the discrete Fourier transform of the sequence $R(nT),\ n = 0,\ \pm 1, \ldots\,$.
It is shown by Papoulis [Papoulis, pp. 42–49] that the Fourier transform of the train of
delta functions $S_T(\tau) = \sum_{n=-\infty}^{\infty} \delta(\tau - nT)$ is given by

```{math}
:label: eq-15-2
\sum_{n=-\infty}^{\infty} \delta(\tau - nT) \leftrightarrow w_0 \sum_{n=-\infty}^{\infty} \delta(w - nw_0),\ w_o = 2\pi/T.
```

Now notice that

```{math}
:label: eq-15-3
w_0 \sum_{n=-\infty}^{\infty} S(w - nw_0) = S(w) \ast w_0 \sum_{n=-\infty}^{\infty} \delta(w - nw_0).
```

The inverse Fourier transform of the right hand side of {eq}`eq-15-3` is

$$
\frac{1}{2\pi}\ \sum_{n=-\infty}^{\infty} R(\tau)\, \delta(\tau - nT) = \frac{1}{2\pi}\ R^d(\tau).
$$

It follows from (—) and (—) that the spectral density of the discrete process $x_j$
satisfies

$$
S^d(w) = \frac{1}{T} \sum_{n=-\infty}^{\infty} S\!\left(w - n\, \frac{2\pi}{T}\right)
$$

where $S^d(w)$ is the spectral density of the discrete data and $S(w)$ is the spectral
density of the continuous time data. Equation (—) is known as the *folding formula*.
