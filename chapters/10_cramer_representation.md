# 10. The Cramér Representation

The spectral density accomplishes an orthogonal decomposition by frequency of the variance of a covariance stationary stochastic process. The Cramér representation exhibits this decomposition in a general way.

Let $x(t)$ be a covariance stationary, zero mean stochastic process with autocorrelation function $R(\tau)$. Let the spectral density of $x$ be denoted $S(w)$, where

$$
S(w) = \int^\infty_{-\infty} R(\tau) e^{-iw\tau} d\tau.
$$

Associated with $x(t)$, there exists a complex valued random process $Z(w)$ defined by

(i) $Z(0) = 0$

(ii) $Z(-w) = \overline{Z(w)}$, where the bar denotes complex conjugation.

(iii) $Z(w)$ is a process with orthogonal increments, and in particular

$$
E Z'(\lambda) \overline{Z'(\mu)} = S(\mu) \delta(\mu - \lambda),
$$

where the prime denotes the (possibly generalized) derivative of $Z(\lambda)$ with respect to $\lambda$. The process $Z(\lambda)$ is called the "random spectral measure" of the $x(t)$ process. In terms of this random process, the $x(t)$ process has the Cramér representation

$$
x(t) = \frac{1}{\sqrt{2\pi}}\, \int^\infty_{-\infty} e^{i\lambda t} dZ(\lambda)
$$

or

```{math}
:label: eq-10-cramer
x(t) = \frac{1}{\sqrt{2\pi}} \int^\infty_{-\infty} e^{i\lambda t}\, Z'(\lambda) d\lambda.
```

To help motivate this representation, we use (—) to calculate $R(\tau) = Ex(t) x(t-\tau)$,

$$
\begin{aligned}
Ex(t) x(t-\tau) &= E\, \frac{1}{2\pi} \int^\infty_{-\infty} e^{i\lambda t}\, Z' d\lambda \int^\infty_{-\infty} e^{-i\mu (t-\tau)}\, \overline{Z'(\mu)} d\mu \\
&= \frac{1}{2\pi} \int^\infty_{-\infty} e^{-i\mu (t - \tau)} \int^\infty_{-\infty} e^{i\lambda t}\, EZ'(\lambda)\, \overline{Z'(\mu)}\, d\lambda d\mu \\
&= \frac{1}{2\pi} \int^\infty_{-\infty} e^{-i\mu (t - \tau)} \int^\infty_{-\infty} e^{i\lambda t}\, S(\mu) \delta(\mu - \lambda) d\lambda d\mu \\
&= \frac{1}{2\pi} \int^\infty_{-\infty} e^{-i\mu (t - \tau)}\, e^{i\mu t}\, S(\mu) d\mu
\end{aligned}
$$

or

$$
R(\tau) = \frac{1}{2\pi} \int^\infty_{-\infty} e^{i\mu \tau} S(\mu) d\mu.
$$

This is the inversion formula (—) for recovering the autocorrelation function from the spectral density.

There is a sense in which the random spectral measure can be defined formally as the Fourier transform of $x(t)$,

$$
Z'(\lambda,\, w) = \frac{1}{\sqrt{2\pi}}\, \int^\infty_{-\infty} e^{-i\lambda t}\, x(t,\, w) dt,
$$

provided that the integral is interpreted delicately. In (—), we have added the argument $w \in \Omega$ explicitly to emphasize the both $x(t,\, w)$ and $Z'(\lambda,\,w)$ are random processes defined on the same underlying probability space $(\Omega,\, \mathcal{F},\, P)$.

Differentiating {eq}`eq-10-cramer` formally with respect to time, we have that the mean square derivative of $x(t)$, if it exists, has Cramér representation

$$
x'(t) = \frac{1}{\sqrt{2\pi}} \int^\infty_{-\infty} i\, \lambda\, e^{i\lambda t}\, Z'(\lambda) d\lambda.
$$

The derivative of the random spectral measure associated with $x'(t)$ is thus $\lambda Z'(\lambda)$, which being proportional to $Z'(\lambda)$ is itself the derivative of a process with independent increments and obeys

$$
E \lambda Z'(\lambda)\ \overline{\mu Z'(\mu)} = \mu^2 S(\mu) \delta(\mu - \lambda).
$$

Thus the spectral density of $Dx(t)$ is $\mu^2 S(\mu)$. More generally, consider the distributed lag

$$
y(t) = \int^\infty_{-\infty} b(\tau) x(t - \tau) d\tau
$$

where $b(\tau) \in L_2\, (-\infty,\, \infty)$. Then using (—), we have

```{math}
:label: eq-10-1
\begin{aligned}
y(t) &= \frac{1}{\sqrt{2\pi}} \int^\infty_{-\infty} b(\tau) \int^\infty_{-\infty} e^{i\lambda (t-\tau)}\, Z'(\lambda) d\lambda \\
&= \frac{1}{\sqrt{2\pi}} \int^\infty_{-\infty} \int^\infty_{-\infty} b(\tau)\, e^{-i\lambda \tau}\, d\tau\, e^{i\lambda t}\, Z'(\lambda) d\lambda \\
y(t) &= \frac{1}{\sqrt{2\pi}} \int^\infty_{-\infty} e^{i\lambda t}\, B(\lambda)\, Z'(\lambda) d\lambda
\end{aligned}
```

where $B(\lambda)$ is the Fourier transform of $b(t)$, namely

$$
B(\lambda) = \int^\infty_{-\infty} b(\tau) e^{-i\lambda \tau}\, d\tau.
$$

From {eq}`eq-10-1`, it follows that the random spectral measure of $y(t)$ has derivative $B(\lambda) Z'(\lambda)$, and that $y(t)$ has spectral density

$$
B(\lambda)\ \overline{B(\lambda)}\ S(\lambda)
$$

or

```{math}
:label: eq-10-2
S_y(w) = |B(w)|^2 S_x(w)
```

where $S_x(w) \equiv S(w)$ is the spectral density of $x(t)$.

Equation {eq}`eq-10-2` can be used to view from another angle the orthogonal decomposition across frequencies induced by $S_x(w)$. For $0 < c < d$, we create the "window" $B_{cd}(w)$ according to

$$
B_{cd}(w) = \begin{cases} 1 & w \in [c,\, d]\ \text{ or }\ [-d,\,-c] \\ 0 & \text{otherwise} \end{cases}
$$

$$
\text{where }\ d > c > 0.
$$

The band-pass filter $B_{cd}(w)$ is illustrated in {numref}`fig-10-1`.

```{figure} figures/fig-10-1_bandpass_window.png
:name: fig-10-1
:width: 90%
:align: center

Figure 1. The band-pass filter ("window") $B_{cd}(w)$, defined on the frequency axis $w \in (-\infty,\, \infty)$. It equals $1$ on the two symmetric frequency bands $[c,\, d]$ and $[-d,\, -c]$ and equals $0$ everywhere else, with $d > c > 0$.
```

Define the time function

$$
\begin{aligned}
b_{cd}(t) &= \frac{1}{2\pi} \int^\infty_{-\infty} B_{cd}(w) e^{+iwt}\, dw \\
&= \frac{1}{\pi}\ \left[ \frac{\sin\ dt}{t}\ - \ \frac{\sin\ ct}{t}\right]
\end{aligned}
$$

We have that $y_{cd}(t)$, defined by

$$
y_{cd}(t) \equiv \int^\infty_{-\infty} b_{cd}(\tau) x(t-\tau) d\tau,
$$

has the spectral density $S_{cd}$ given by

$$
S_{cd}(w) = \begin{cases} S_x(w) & w \in [c,\, d]\ \text{ or }\ [-d,\, -c] \\ 0 & \text{otherwise} \end{cases}
$$

If we choose an interval parameterized by $0 < e < f$ such that $[e,\,f] \cap [c,\, d] = 0$, it follows from the orthogonal increments property of $Z(\lambda)$ that $y_{ef}(t)$ is orthogonal to $y_{cd}(t-\tau)$ for all $t$. For using the Cramér representation, we have

$$
\begin{aligned}
Ey_{cd}(t) y_{ef}(t-\tau) &= \frac{1}{2\pi} \int^\infty_{-\infty} e^{i\lambda t}\, B_{cd}(\lambda) Z'(\lambda) \\
&\quad \int^\infty_{-\infty} e^{-i\mu (t-\tau)}\, \overline{B_{ef}(\mu)}\, Z'(\mu) d\mu d\lambda \\
&= \frac{1}{2\pi} \int^\infty_{-\infty} e^{-i\mu (t - \tau)}\, \overline{B_{ef}(\mu)} \\
&\quad \int^\infty_{-\infty} B_{cd}(\lambda) S(\mu) \delta(\mu - \lambda) d\lambda d\mu \\
&= \frac{1}{2\pi} \int^\infty_{-\infty} e^{-i\mu (t - \tau)}\, \overline{B_{ef}(\mu)}\ B_{cd}(\mu) S(\mu) d\mu \\
&= 0
\end{aligned}
$$

since $\overline{B_{ef}(\mu)}\ B_{cd}(\mu) = 0$ for all $\mu$.

We also have that $y_{ef}(t)$ has spectral density

$$
S_{ef}(w) = \begin{cases} S_x(w) & w \in [e,\, f]\ \text{ or }\ [-f,\, -e] \\ 0 & \text{otherwise} \end{cases}
$$

By filtering the process $x(t)$ with filters formed by taking disjoint frequency intervals, $[c,\,d]$, $[e,\,f]$, we create processes that are orthogonal at all leads and lags, but whose respective spectral densities equal that of $x(t)$ over the intervals $[c,\,d]$ and $[e,\,f]$.

To motivate from another angle the frequency decomposition induced by the spectral density, it is useful to note that the Cramér representation

$$
x(t) = \frac{1}{\sqrt{2\pi}} \int^\infty_{-\infty} e^{i\lambda t}\, Z'(\lambda) d\lambda
$$

can be also represented in a couple of alternative ways. Let the complex valued random process $Z'(\lambda)$ be represented in the alternative forms

$$
\begin{aligned}
z'(\lambda) &= a(\lambda) + b(\lambda)i \\
z'(\lambda) &= r(\lambda)e^{i\theta(\lambda)}.
\end{aligned}
$$

Then it is straightforward to derive the Cramér representation in the forms

```{math}
:label: eq-10-3
x(t) = \frac{1}{\pi} \int^\infty_0\ a(\lambda)\ \cos\ \lambda t d\lambda \ - \ \frac{1}{\pi} \int^\infty_0 b(\lambda)\ \sin\ \lambda t d\lambda
```

where

$$
E (a(\lambda) + ib(\lambda))\ (a(\mu) - ib(\mu)) = S(\mu) \delta(\mu - \lambda);
$$

and

```{math}
:label: eq-10-4
x(t) = \frac{1}{\pi} \int^\infty_0 r(\lambda)\ \cos\ (\lambda t + \theta(\lambda)) d\lambda
```

where

$$
Er(\lambda) e^{i\theta(\lambda)}\, r(\mu) e^{-i\theta(\mu)} = S(\mu) \delta(\mu - \lambda).
$$

Representation {eq}`eq-10-3` expresses $x(t)$ as a weighted sum of cosine and sine waves, with the weights being random processes $a(\lambda)$ and $b(\lambda)$. Equation {eq}`eq-10-4` represents $x(t)$ as a sum of cosine waves with amplitude $r(\lambda)$ and phase $\theta(\lambda)$ being governed by random processes.

Using representation {eq}`eq-10-4`, it is possible to show that

$$
\begin{aligned}
Ex(t)^2 &= \frac{1}{\pi} \int^\infty_0 (a(\lambda)^2 + b(\lambda)^2) d\lambda \\
&= \frac{1}{2\pi} \int^\infty_{-\infty} S(\lambda) d\lambda.
\end{aligned}
$$
