# 6. The Wiener Process

We shall now describe a class of stochastic processes that can be described by stochastic differential equations driven by a "Wiener Process." We shall approach the Wiener process by regarding it as the limit as $\lambda \to \infty$ of a process $x(t)$ that solves

$$
dx(t) = 1/\sqrt\lambda\, \left(dN_1(t) - dN_2(t)\right)
$$

where $N_1$ and $N_2$ are independent Poisson counters with rates $\lambda/2$. We have already seen that for $\lambda > 0$, this process is a random walk with jumps of size $\pm\, 1/\sqrt\lambda$ occurring at times that arrive at a mean rate of $\lambda/2$ per unit of time. As we drive $\lambda$ toward infinity, the jump size $1/\sqrt\lambda$ goes to zero, while the average number of arrivals or jumps per unit of time goes to infinity. We have seen that the mean function and autocorrelation function of the solution $x(t)$ of (—) are given by

$$
\begin{aligned}
E x(t) &= 0 \\
R(t,\, t + \tau) &= \min\ (t,\, t + \tau)
\end{aligned}
$$

and so are independent of $\lambda$. These functions also characterize the process that results in the limit as $\lambda \to \infty$. It turns out that the sample paths of the limiting process is continuous (almost everywhere) but differentiable nowhere. Roughly speaking, the sample paths are continuous because the jumps in the process become infinitesimal as $\lambda \to \infty$; the sample paths are differentiable nowhere because there are "so many" of these infinitesimal jumps in a small interval of time.

We begin by computing the $p^{th}$ moments of the process governed by

$$
dx(t) =\ \frac{1}{\sqrt\lambda}\ \left(dN_1(t) - dN_2\,(t)\right)
$$

where $N_1$ and $N_2$ are two independent Poisson counters with identical rates $\lambda/2 > 0$. Using rule (—), we have that for $p > 0$

$$
\begin{aligned}
dx^p &= \left[(x + \frac{1}{\sqrt\lambda}\,)^p - x^p\right]\, dN_1 + \left[(x - \frac{1}{\sqrt\lambda})^p - x^p\right]\, dN_2 \\
&= \left( px^{p-1}\ \frac{1}{\sqrt\lambda}\ + \binom{p}{2}\, x^{p-2}\ \left(\frac{1}{\sqrt\lambda}\right)^2\ + \cdots + \left(\frac{1}{\sqrt\lambda}\right)^p\right)\, dN_1 \\
&+ \left( -px^{p-1}\ \frac{1}{\sqrt\lambda}\ + \binom{p}{2}\, x^{p-2}\ \left(\frac{1}{\sqrt\lambda}\right)^2 + \cdots \pm\, \left(\frac{1}{\sqrt\lambda}\right)^p\right)\, dN_2.
\end{aligned}
$$

It follows from rule (—) that

$$
\begin{aligned}
\frac{d}{dt}\ E x^p &= \binom{p}{2}\ \left(\frac{1}{\lambda}\right)\ E x^{p-2}\, \lambda + \binom{p}{4}\ \left(\frac{1}{\lambda}\right)^2\ E x^{p-4}\, \lambda + \cdots + \left(\frac{1}{\sqrt\lambda}\right)^{p-1}\, \lambda \\
p\ &\text{ even, }\ p \geq 2
\end{aligned}
$$

$$
\begin{aligned}
\frac{d}{dt}\ E x^p &= \binom{p}{2}\ \left(\frac{1}{\lambda}\right)\ E x^{p-2}\, \lambda + \binom{p}{4}\ \left(\frac{1}{\lambda}\right)^2\, E x^{p-4}\, \cdot\, \lambda + \ldots + \left(\frac{1}{\sqrt\lambda}\right)^{p-1}\ E x \cdot \lambda \\
&\ p\ \text{ odd}
\end{aligned}
$$

We have calculated that $E x(t) = 0$, and that $E x(t)^2 = t$. Using (—) recursively to calculate odd moments, we find that

$$
E x(t)^p = 0\ \text{ for }\ p\ \text{ odd.}
$$

For even powers of $p$ higher than 2, we drive $\lambda \to \infty$ in (—) and find that

$$
\frac{d}{dt}\ E x^p = \binom{p}{2}\, E x^{p-2},
$$

which has the solution

$$
E x^p = \frac{1}{2}\, \int_0^t p\, (p-1)\ E x^{p-2}\, (s)ds.
$$

Using this equation recursively, starting from $p = 2$ gives

$$
\begin{aligned}
E x(t)^2 &= t \\
E x(t)^4 &= \frac{1}{2} \int_0^t 4 \cdot 3 \cdot sds = 3 t^2 \\
E x(t)^6 &= 15 t^3
\end{aligned}
$$

and so on.

Now for the Gaussian density function

$$
f(x) = (2 \pi \sigma)^{-1/2}\ \exp\ \frac{-x^2}{2 \sigma}\, ,
$$

integration by parts shows that $E x^p = 0$ for $p$ odd, while for $p$ even

$$
E x^p = \sigma^{p/2}\ (p - 1)\ (p - 3)\ \ldots\ 1.\ p\ \text{ even}
$$

It is known that if all of the moments of a process equal those of a Gaussian process, then that process is itself Gaussian. Comparing (—) with (—), we can conclude that the limiting process as $\lambda \to \infty$ is *Gaussian* with mean zero and variance $t$. That is, $x(t)$ has density

$$
f\left(x(t)\right) = (2 \pi t)^{-1/2}\ \exp\ \frac{-x(t)^2}{2 t}\, ,\ t \geq 0.
$$

We shall denote the limiting process as $\lambda \to \infty$ as $W(t)$. It is known as *Brownian motion* or a *Wiener process*.

The Wiener process $W(t)$ is characterized by the following properties.

(i) $W(0) = 0$

(ii) $W(t_1) - W(t_2)$ and $W(t_3) - W(t_4)$ are independent if $[t_1,\, t_2]$ and $[t_3,\, t_4]$ are disjoint intervals.

(iii) $E(W(t) - W(\tau))^2 = t - \tau$.

(iv) Sample paths or realizations of $W(t)$ are continuous with probability 1.

Property (ii) is the *independent increments* property, and is inherited from the Poisson process of which $W(t)$ is the limit.
