---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.1
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# 8. Spectral Densities

```{eval-rst}
.. index::
   single: spectral density
   single: power spectrum
   single: covariogram
   single: Wold decomposition
   single: spectral factorization theorem
   single: fundamental white noise
   single: non-fundamentalness
   single: linearly deterministic process
   single: linearly indeterministic process
   single: periodogram
   single: spectral smoothing
   single: Fourier transform; tables of
```

## (a) General Results

Let $x(t)$ be a covariance stationary stochastic process with covariogram

$$
Ex(t)x(t-\tau) = R(\tau), \tau \in R,
$$

where $R$ is a positive semidefinite function. The *power spectrum* or the *spectral
density* of $x$ is defined as

```{math}
:label: eq-8-1
S(w) = \int_{-\infty}^{\infty} e^{-iw\tau}\, R(\tau)\, d\tau,\ w \in (-\infty,\, \infty).
```

Equation {eq}`eq-8-1` defines $S(w)$ as the Fourier transform of $R(\tau)$. Given $S(w)$,
$R(t)$ can be recovered from the inverse Fourier transform

```{math}
:label: eq-8-2
R(\tau) = \frac{1}{2\pi}\ \int_{-\infty}^{\infty} e^{+iw\tau}\, S(w)\, dw
```

From the fact that $R(\tau)$ is a positive semidefinite function, which implies that
$R(\tau) = R(-\tau)$, definition {eq}`eq-8-1` and the inverse relation {eq}`eq-8-2` imply
that $S(w)$ has the following properties

(i) $S(w) = S(-w)$. (from definition {eq}`eq-8-1`)

(ii) $S(w) \geq 0$. (from definition {eq}`eq-8-1` and positive semidefiniteness of $R(t)$)

(iii) $\frac{1}{2\pi} \int_{-\infty}^{\infty} S(w)\, dw = R(0)\quad$ ({eq}`eq-8-2` evaluated
at $\tau = 0$)

Property (iii) asserts that the spectral density achieves a decomposition of the variance
$R(0)$ by frequency. We shall shortly see that this decomposition of variance is actually
into components that are orthogonal across different frequencies.

For the purposes of using {eq}`eq-8-1` and {eq}`eq-8-2`, it is fortunate that tables of
Fourier transform pairs have been prepared. Table 1 is a small table. The reader can verify
the entries in the table by using {eq}`eq-8-1` or {eq}`eq-8-2`.

Fourier transforms have a number of useful operational properties. In particular, let
$F(w) \leftrightarrow f(t)$ mean that $F(w)$ is the Fourier transform of the "time signal"
$f(t)$, i.e.,

$$
F(w) = \int_{-\infty}^{\infty} f(t)\, e^{-iwt}\, dt.
$$

Then a number of simple operations of $F(w)$ or $f(t)$ can be used to generate other
Fourier transform pairs. Some of these are recorded in table 2.

Property (8) of table 2 can be used to find the spectral density of mean square
derivatives of a given process $x(t)$. We saw above that if $x_t$ is a covariance
stationary stochastic process with autocovariogram $R(\tau)$, then the autocovariogram of
$\frac{d^n}{dt^n}\ x(t)$ is $(-1)^n\ \frac{d^{2n} R(\tau)}{d\tau^{2n}}$. It follows from
property 8 of table 2 that the Fourier transform of
$(-1)^n \ \frac{d^{2n}R(\tau)}{d\tau^{2n}}$ is $w^{2n}F(w)$, where $F(w)$ is the Fourier
transform of $R(\tau)$.

The preceding property can be viewed as a limiting case of a more general property. Let
$x(t)$ be a covariance stationary process with covariogram $R(\tau)$. Define $y(t)$ as

$$
y(t) = \int_{-\infty}^{\infty} h(\tau) x(t-\tau)\, d\tau
$$

where $h(\tau) \in L_2\, [-\infty,\, \infty]$. Then calculating $Ey(t)y(t-s)$ gives

```{math}
:label: eq-8-4
\begin{aligned}
R_y (s) &= Ey(t) y(t-s) = E \int_{-\infty}^{\infty} h (\tau) x (t - \tau)\, d\tau\ \int_{-\infty}^{\infty} h(v) x(t-s-v)\, dv \\
&= \int_{-\infty}^{\infty} h(\tau) \int_{-\infty}^{\infty} R_x \left(\tau - (s + v)\right) h (v)\, dv\, d\tau \\
&= \int_{-\infty}^{\infty} h(\tau) R_x\ \ast\ h(\tau - s)\, d\tau \\
R_y (s) &= h\ \ast \ R_x\ \ast \ h(-s)
\end{aligned}
```

Applying property 6 of table 2, the convolution property, twice to {eq}`eq-8-4`, together with
property 12 (reversal) to handle the $h(-s)$ factor, gives

```{math}
:label: eq-8-5
S_y (w) = h(w) S_x (w) h(-w)
```

where $h(w) = \int_{-\infty}^{\infty} h (\tau) e^{-i w \tau}\, d\tau$ is the Fourier
transform of $h(t)$.

The sense in which our preceding result is a special case of {eq}`eq-8-4` is as follows.
The generalized derivative of $\delta (t),\ \delta' (t)$, is a generalized function defined
by

```{math}
:label: eq-8-6
\int_{-\infty}^{\infty} \delta' (t) g(t)\, dt = - g' (0)
```

where $g' (t)$ is a test function that is differentiable at $t = 0$. More generally, the
$n^{th}$ generalized derivative of $\delta (t)$ is defined by

$$
\int_{-\infty}^{\infty} \delta^{(n)}(t) g(t)\, dt = (-1)^n g^{(n)} (0)
$$

where $g(t)$ is a test function that is $n$ times differentiable at $t$. See Papoulis for
definitions of $\delta (t)$ and its generalized derivatives. Both $\delta (t)$ and
$\delta^{(n)}\, (t)$ can be defined in terms of particular limit points of sequences of
functions $h_n (t) \in L_2 [-\infty, \, \infty]$ for each $n$, which are used as
integrating functions in

$$
\int_{-\infty}^{\infty} h_n (t) g(t)\, dt
\qquad \text{and} \qquad
\int_{-\infty}^{\infty}\, \left( \frac{d^n}{dt^n}\right)\, h_n (t)\, g(t)\, dt .
$$

The Fourier transform of $\delta' (t)$ is seen to be, from the definition

```{math}
:label: eq-8-3
\int_{-\infty}^{\infty} \delta' (t) e^{-iwt}\, dt =\ -\ \frac{d}{dt}\ \ e^{-iwt}\Big|_{t=0} = iw.
```

More generally, the Fourier transform of $\delta^{(n)} (t)$ is given by

$$
\delta^{(n)} (t) \leftrightarrow (iw)^n.
$$

Using this result in conjunction with {eq}`eq-8-5` with $h(t) = \delta^{(n)} (t)$
immediately gives that the spectrum of the $n^{th}$ mean square derivative of a process
$x(t)$ with spectrum $S(w)$ is given by $w^{2n} S(w)$.

We record some of these results in table 3.

## Tables

### Table 1

$$
F(w) = \int_{-\infty}^{\infty} f(t) e^{-iwt}\, dt
$$

| $f(t)$ | $F(w)$ |
| --- | --- |
| $e^{-a\vert t\vert}$ | $\frac{2a}{a^2 + w^2}$, $a > 0$ |
| $u(t)$ | $\frac{1}{iw} + \pi\, \delta(w)$ |
| $u(t)t$ | $\frac{1}{(iw)^2} + i\pi\, \delta'(w)$ |
| $\delta (t)$ | $1$ |
| $\delta' (t)$ | $iw$ |
| $\delta^{(n)}\, (t)$ | $(iw)^n$ |
| $\sum_{n=-\infty}^{\infty} \delta (\tau-nT)$ | $w_o \sum_{n=-\infty}^{\infty} \delta (w-n w_0),\ w_0 = \frac{2\pi}{T}$ |
| $e^{-at}\, u(t)$ | $\frac{1}{a + iw}$, $a > 0$ |
| $\begin{cases} 1 & t \in [0,\,1] \\ 0 & t > 1,\ t < 0 \end{cases}$ | $\left(\frac{1 - e^{-iw}}{iw}\right)$ |
| $\begin{cases} 1 & t \in [0,\,1] \\ -1 & t \in [1,\, 2] \\ 0 & t > 2,\ t < 0 \end{cases}$ | $\frac{(1 - e^{-iw})^2}{iw}$ |
| $\begin{cases} t, & t \in [0,\,1] \\ 2-t, & t \in [1,\,2] \\ 0 & t < 0,\ t > 2 \end{cases}$ | $\frac{(1-e^{-iw})^2}{(iw)^2}$ |
| $\frac{t^{n-1}}{(n-1)!}\, e^{-at}\, u(t)$ | $\frac{1}{(a+iw)^n}$, $a > 0$ |
| $e^{iw_0t}$ | $2\pi \delta (w-w_0)$ |
| $\cos\ w_0 t$ | $\pi\, \left[\delta (w-w_0) + \delta\, (w+w_0)\right].$ |

Note: $u(t) = \begin{cases} 1 & t \geq 0 \\ 0 & t < 0 \end{cases}\qquad$ (Heaviside step
function)

The entries for $u(t)$ and $u(t)t$ are the transforms in the sense of generalized functions;
the $\delta$ terms record that these signals do not decay, and they carry all of the mass at
$w = 0$. In every use we make of these entries — the operational calculus of
{doc}`11_linear_sde` and {doc}`12_prediction`, and the nonstationary examples of
{doc}`14_nonstationary_examples` — the operators $1/(iw)$ and $1/(iw)^2$ act on a noise that has
no mass at the origin, so the $\delta$ terms may be, and are, dropped. Writers who work only
with such operators often tabulate $u(t) \leftrightarrow 1/iw$ for that reason.

$\delta (t) =$ Dirac delta generalized function, defined by
$\int_{-\infty}^{\infty} g(t) \delta(t)\, dt = g(0)$ where $g(t)$ is a "test function" that
is continuous at $t=0$.

### Table 2

$$
f(t) \leftrightarrow F(w)\ \text{ means }\ F(w) = \int_{-\infty}^{\infty} f(t) e^{-iwt}\, dt.
$$

Property:

| Property | Transform pair |
| --- | --- |
| 1. Linearity | $a_1 f_1 (t) + a_2 f_2(t) \leftrightarrow a_1 F_1(w) + a_2 F_2 (w)$; $a_1,\ a_2$ are scalars. |
| 2. Symmetry | $F(t) \leftrightarrow 2\pi f(-w)$ |
| 3. Scaling | $f(at) \leftrightarrow \frac{1}{\vert a\vert}\ F\left(\frac{w}{a}\right),\ \text{ a scalar}$ |
| 4. Delay | $f(t-t_0) \leftrightarrow e^{-iwt_o}\, F(w)$ |
| 5. Modulation | $e^{iw_0t}\, f(t) \leftrightarrow F(w-w_0)$ |
| 6. Convolution | $f_1(t) \ast f_2(t) \leftrightarrow F_1(w)F_2(w)$; where $f_1(t) \ast f_2(t) \equiv \int_{-\infty}^{\infty} f_1(t-\tau) f_2(\tau)\, d\tau$ |
| 7. Multiplication | $f_1(t) f_2 (t) \leftrightarrow \frac{1}{2\pi} F_1 (w) \ast F_2 (w)$; where $F_1(w) \ast F_2(w) \equiv \int_{-\infty}^{\infty} F_1(w-s) F_2(s)\, ds$ |
| 8. Time differentiation | $\frac{d^n}{dt^n}\, f(t) \leftrightarrow (iw)^n\ F(w)$ |
| 9. Time integration | $\int_{-\infty}^{t} f(\tau)\, d\tau \leftrightarrow\ \frac{F(w)}{iw} + \pi\ F(0)\, \delta (w)$ |
| 10. Frequency differentiation | $-itf(t)\ \leftrightarrow\ \frac{dF(w)}{dw}$ |
| 11. Frequency integration | $\frac{f (t)}{-it}\ \leftrightarrow\ \int F(w') dw'$ |
| 12. Reversal | $f(-t) \leftrightarrow F(-w)$ |

### Table 3

| Process | Autocovariogram | Spectrum |
| --- | --- | --- |
| $x_t$ | $R(\tau)$ | $S(w)$ |
| $\int_{-\infty}^{\infty} h(\tau) x(t-\tau)\, d\tau$ | $h\, \ast\, R\, \ast\, h(-s)$ | $\vert h(w)\vert^2 S(w)$ |
| $\frac{d}{dt}\ x(t)$ | $-\ \frac{d^2 R(\tau)}{d\tau^2}$ | $w^2 S(w)$ |
| $\frac{d^n}{dt^n}\, x(t)$ | $(-1)^n\ \frac{d^{2n}\, R(\tau)}{d\tau^{2n}}$ | $w^{2n} S(w)$ |

### Table 4

| $x(t)$ | $R(\tau)$ | $S(w)$ |
| --- | --- | --- |
| $w(t)$ | $\delta (\tau)$ | $1$ |
| $(D - \lambda) x(t) = w(t),\, \lambda < 0$ | $\frac{-1}{2\lambda}\ e^{\lambda \vert\tau\vert}$ | $\frac{1}{\lambda^2 + w^2}$ |
| $\theta (D) x(t) = w(t)$ | | $\frac{1}{\theta (iw)\, \theta (-iw)}$ |
| or | | or |
| $(D - \lambda_1)\ldots (D-\lambda_n)\, x(t) = w(t)$ | $\sum_{j=1}^{n}\, k_j\, e^{\lambda_j \vert\tau\vert}$ | $\frac{1}{\Pi_{j=1}^{n}\, (iw-\lambda_j)\, (-iw-\lambda_j)}$ |
| $re\, (\lambda_j) < 0,\ j=1,\, \ldots,\, n$ | $k_j = \lim_{s \to \lambda_j}\ \frac{(s - \lambda_j)}{\theta (s)\, \theta (-s)}$ | |
| $\theta (D) x(t) = \Psi(D) w(t)$ | $\Sigma\, k_j\, e^{\lambda_j \vert\tau\vert}$ | |
| where $\theta (D) = (D-\lambda_1) \ldots (D-\lambda_n)$ | | |
| $re (\lambda_j) < 0$ | $k_j = \lim_{s\to \lambda_j}\ (s-\lambda_j)\ \frac{\Psi (s)}{\theta (s)}\ \frac{\Psi (-s)}{\theta (-s)},$ | $\frac{\Psi (iw)\, \Psi (-iw)}{\theta (iw)\, \theta(-iw)}$ |
| $\Psi (D) = \Psi_0 + \Psi_1 D + \ldots + \Psi_m \, D^m$ | | |
| $m < n$ | | |

## (b) Wold's Theorem

We now state a version of Wold's decomposition theorem in continuous time.

**Theorem 10.** Let $x(t)$ be a covariance stationary stochastic process with autocovariance
function $R(\tau)$ and spectral density $S(w)$. Then $x(t)$ can be represented as

```{math}
:label: eq-8-a
x(t) = \int_0^{\infty} p(\tau)w (t-\tau)\, d\tau + x^d(t)
```

where $Ex^d (t) \cdot \int_0^{\infty} p(\tau) w(t - s - \tau)\, d\tau = 0$ for all $s$, so that
$x^d (t)$ is orthogonal, at all leads and lags, to the moving average
$\int_0^{\infty} p(\tau) w(\,\cdot\, - \tau)\, d\tau$. In
{eq}`eq-8-a`, $w(t)$ is a white noise with autocovariogram

$$
E w(t) w(t-s) = \delta (s),
$$

and $p(\tau)$ is a square integrable function,

$$
\int_0^{\infty} p(\tau)^2\, d\tau < + \infty.
$$

Furthermore, $w(t)$ is a *fundamental white noise for $x(t)$*, which means that minimum
mean squared error $s$-step ahead errors in forecasting $x(t+s)$ as a square integrable
linear functional of $[x(v), v \leq t]$ can be expressed as an integral of $w(\tau)$ for
$t \leq \tau \leq t + s$; in particular,

$$
x(t+s) - E\, \left[ x (t+s) \mid x(v),\, v \leq t\right] = \int_0^{s} p(\tau) w(t+s-\tau)\, d\tau.
$$

In {eq}`eq-8-a`, $x^d (t)$ is a *linearly deterministic* process that can be forecast
arbitrarily far into the future by a linear function of its own past values, or past values
of $x$.

The spectral density of $x^d (t)$ can be represented as

```{math}
:label: eq-8-b
S_d (w) = \sum_{j=1}^{\infty} a_j \pi \left[\delta (w-w_j) + \delta (w + w_j)\right]
```

where $a_j$ are positive constants and where $w_j,\ j = 1,\ \ldots$ is a countable set of
frequencies. Therefore, the autocovariance function of $x^d (t)$ is given by

```{math}
:label: eq-8-c
R^d (\tau) = \sum_{j=1}^{\infty} a_j\ \cos\, (w_j \tau).
```

It follows from {eq}`eq-8-a`, {eq}`eq-8-b` and the convolution property {eq}`eq-8-5` that the
spectral density of $x(t)$ can be represented as

$$
S (w) = P(w) P(-w) + \sum_{j=1}^{\infty} a_j \pi \left[\delta (w-w_j) + \delta (w + w_j)\right]
$$

where $P(w) = \int_0^{\infty} p(\tau) e^{-iw\tau}$, is the Fourier transform of a one-sided,
square-integrable function. The component $\int_0^{\infty} p(\tau) w(t - \tau)\, d\tau$ is
called the *linearly indeterministic* part of the process $x(t)$.

## (c) The Spectral Factorization Theorem

Included in the statement of Wold's theorem is the spectral factorization theorem for
linearly indeterministic processes. We restate this property separately:

**Spectral factorization theorem.** Let $x(t)$ be a covariance stationary, *linearly
indeterministic* process with spectral density $S(w)$. Then $S(w)$ can be factored as

$$
S(w) = \tilde P (iw) \tilde P(-iw)
$$

where $\tilde P (s) = \int_0^{\infty} p(t) e^{-st}\, dt$, where
$\int_0^{\infty} p(t)^2\, dt < + \infty$, and $\tilde P (s)$ has no zero for $s$ in the
right half of the complex plane. This condition on the zeros of $\tilde P(s)$ is the
condition that the white noise $w(t)$ in
$x(t) = \int_0^{\infty} p(\tau) w(t-\tau)\, d\tau$ is fundamental for $x(t)$. The function
$\tilde P(s)$ is the Laplace transform of $p(\tau)$. The Fourier transform $P(w)$ is related
to $\tilde P(s)$ by $P(w) = \tilde P (iw)$.

This factorization is the structural result on which much of the rest of the book leans. It
furnishes the prediction formulas of {doc}`12_prediction`, acquires a time-domain,
state-space counterpart in the Kalman–Bucy filter and Riccati equation of
{doc}`15_kalman_filter_spectral_factorization` (solving the Riccati equation is the
time-domain algorithm that performs this factorization), and is the object identified from
discretely sampled data in {doc}`21_phillips_continuous_time_estimation`.

Thus, we can represent a linearly indeterministic, covariance stationary process $x_t$ as

$$
x(t) = \int_0^{\infty} p(\tau) w(t-\tau)\, d\tau
$$

or

$$
x(t) = \tilde P (D) w(t)
$$

where $D$ is the time derivative operator, and
$\tilde P (s) = \int_0^{\infty} e^{-st} p (t)\, dt$.

As an example of the construction envisioned in Wold's theorem, consider the strictly
linearly indeterministic process $x(t)$ governed by

$$
\begin{aligned}
x(t) = \frac{D - b}{(D-\lambda_1)\ (D-\lambda_2)}\ v(t),\qquad &b > 0 \\
&re (\lambda_j) < 0,\ j=1,\ 2
\end{aligned}
$$ (eq-8-nonfund)

or

$$
x(t) = \int_0^{\infty} r (\tau) v(t-\tau)\, d\tau
$$

where $v(t)$ is a white noise with

$$
Ev(t)v(t-\tau) = \delta(\tau)
$$

and the Laplace transform of $r(\tau)$ is given by

$$
\tilde R (s) =\ \frac{s - b}{(s - \lambda_1)\ (s-\lambda_2)}.
$$

By using a partial fraction representation of $\tilde R(s)$ we find that

$$
r (\tau) = k_1 e^{\lambda_1 \tau} + k_2\, e^{\lambda_2 \tau},\ \tau \geq 0
$$

where

$$
k_j = \lim_{s\to \lambda_j}\ \frac{s - b}{(s - \lambda_1)\ (s-\lambda_2)}\ (s - \lambda_j).
$$

The spectral density of $x(t)$ is given by

$$
S(w) = \tilde R (iw) \tilde R (-iw).
$$

However, this is not the representation alluded to in the statement of Wold's theorem,
because $\tilde R(s)$ has a zero at $s = b > 0$, which is in the right half plane. This
reflects the fact that the space $H_v(- \infty, \, t)$ is strictly larger than
$H_x(- \infty,\,t)$. This right-half-plane zero is the continuous-time prototype of
*non-fundamentalness* — the driving noise spanning a larger information space than the
observable process — which returns as the central difficulty in interpreting vector
autoregressions in {doc}`18_time_aggregation_var` and as the identification problem of
{doc}`21_phillips_continuous_time_estimation`. To see this heuristically, attempt to invert {eq}`eq-8-nonfund`, and to solve for
$v(t)$ as a function of the $x(t)$ process. This gives

$$
v(t) = \ \frac{(D - \lambda_1)\ (D-\lambda_2)}{(D-b)}\ x(t).
$$

Taking the inverse Laplace transform of
$\frac{(s-\lambda_1)\ (s-\lambda_2)}{(s-b)}$, one obtains an equation of the form

$$
\begin{aligned}
v(t) &= \int_0^{\infty} r_1 (\tau) x (t + \tau)\, d\tau + \int_0^{\infty} r_2 (\tau) Dx(t + \tau)\, d\tau \\
&+ \int_0^{\infty} r_3 (\tau) D^2 x(t + \tau)\, d\tau, \int_0^{\infty} r_j (t)^2\, dt < \infty \\
& \text{for }\ j = 1,\ 2,\ 3
\end{aligned}
$$

This expresses $v(t)$ as a sum of square summable integrals of *future* values of
$x(t),\ Dx(t)$, and $D^2 x(t)$. However, $v(t)$ cannot be expressed in terms of square
summable integrals of *lagged* values of $x$ and its derivatives.

To obtain the Wold representation we note that

$$
\begin{aligned}
S(w) &= \frac{(iw-b)\ (-iw-b)}{(iw-\lambda_1)\ (iw-\lambda_2)\ (-iw-\lambda_1)\ (-iw-\lambda_2)} \\
&= \frac{(iw+b)\ (-iw+b)}{(iw-\lambda_1)\ (iw-\lambda_2)\ (-iw-\lambda_1)\ (-iw-\lambda_2)}.
\end{aligned}
$$

That is, the spectral density remains unaltered if we simply change the sign of the real
part of the zero of the numerator polynomial of $\tilde R (iw)$, so that $-b$ is replaced
by $+b$. It follows that the spectral density can be represented as

$$
S(w) = \tilde P (iw) \tilde P (-iw)
$$

where

$$
\tilde P (iw) = \frac{iw + b}{(iw-\lambda_1)\ (iw-\lambda_2)},
$$

where $\tilde P (s)$ now satisfies the hypotheses required in the statement of Wold's
theorem, in particular that $\tilde P(s)$ have no zeroes in the right half plane, and that
$\tilde P (s)$ be the Laplace transform of a square summable function $p(\tau)$ with support
on $[0,\, \infty)$. It follows that a Wold representation for $x(t)$ is

$$
x(t) = \int_0^{\infty} p(\tau) w(t-\tau)\, d\tau
$$ (eq-8-wold)

where $w(t)$ is a fundamental white noise for $x(t)$, and $p(\tau)$ is the inverse transform
of $\tilde P (iw) = (iw+b)/(iw-\lambda_1)\ (iw-\lambda_2)$. By a partial fractions
representation of $\tilde P (s)$, it follows that

$$
p(\tau) = g_1 e^{\lambda_1 \tau} + g_2 e^{\lambda_2 \tau}\quad \tau \geq 0
$$

where
$g_j = \lim_{s\to \lambda_j}\ (s - \lambda_j)\ \frac{(s + b)}{(s - \lambda_1)\ (s - \lambda_2)}$.

We invite the reader to demonstrate how {eq}`eq-8-wold` can be inverted to express $w(t)$ as a sum of
square summable integrals of past values of $x(t),\ Dx(t)$, and $D^2 x(t)$.

## Exercises

These exercises use the Ornstein–Uhlenbeck process of Chapter 7,

$$
dx(t) = -a\,x(t)\,dt + b\,dW(t), \qquad a, b > 0,
$$

as a running example. From Table 4 (the row $(D-\lambda)x = w$ with $\lambda = -a$, scaled
by $b$), its autocovariance and spectral density are

$$
R(\tau) = \frac{b^2}{2a}\, e^{-a|\tau|}, \qquad
S(w) = \frac{b^2}{a^2 + w^2}.
$$

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt
```

```{exercise-start}
:label: spec_ex1
```

Take $a = 1$, $b = 0.7$.

(a) Plot the spectral density $S(w) = b^2/(a^2 + w^2)$.

(b) Verify property (iii) of the text numerically: the spectral density decomposes the
variance, $\frac{1}{2\pi}\int_{-\infty}^{\infty} S(w)\,dw = R(0) = b^2/(2a)$.

(c) Simulate a long OU path, form its **periodogram** (the sample analogue of the spectral
density), and check that, after smoothing, it tracks $S(w)$.

```{exercise-end}
```

```{solution-start} spec_ex1
:class: dropdown
```

```{code-cell} ipython3
a, b = 1.0, 0.7
S = lambda w: b**2 / (a**2 + w**2)

# (a) the spectral density
w = np.linspace(-8, 8, 400)
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(w, S(w), lw=2)
ax.set_xlabel('$w$'); ax.set_ylabel('$S(w)$')
ax.set_title(r'Spectral density of the OU process, $S(w)=b^2/(a^2+w^2)$')
plt.show()

# (b) variance decomposition
wfine = np.linspace(-1500, 1500, 3_000_001)
lhs = np.trapz(S(wfine), wfine) / (2 * np.pi)
print(f"(1/2π)∫S(w)dw = {lhs:.4f},   R(0) = b^2/2a = {b**2 / (2 * a):.4f}")
```

```{code-cell} ipython3
# (c) periodogram of a simulated OU path vs the theoretical spectrum
rng = np.random.default_rng(5)
dt, T = 0.05, 6000.0
steps = int(T / dt)
x = np.empty(steps); x[0] = 0.0
s = b * np.sqrt(dt)
for k in range(steps - 1):
    x[k + 1] = x[k] - a * x[k] * dt + s * rng.normal()
x = x[int(200 / dt):]            # burn-in
x = x - x.mean()

N = len(x)
freqs = np.fft.rfftfreq(N, d=dt) * 2 * np.pi      # angular frequencies
P = (dt / N) * np.abs(np.fft.rfft(x))**2          # periodogram

# band-average to tame the noise
nb = 200
edges = np.linspace(0, 8, nb + 1)
wc = 0.5 * (edges[1:] + edges[:-1])
Pb = np.array([P[(freqs >= edges[i]) & (freqs < edges[i + 1])].mean() for i in range(nb)])

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(wc, Pb, '.', alpha=0.5, label='periodogram (band-averaged)')
ww = np.linspace(0, 8, 400)
ax.plot(ww, S(ww), 'r-', lw=2, label=r'$S(w)=b^2/(a^2+w^2)$')
ax.set_xlabel('$w$'); ax.set_ylabel('spectral density')
ax.legend()
plt.show()
```

The numerical integral of $S(w)/2\pi$ reproduces $R(0) = b^2/(2a)$, and the smoothed
periodogram of the simulated path follows the theoretical spectrum: most of the variance
sits at low frequencies, with a $1/w^2$ tail. The raw periodogram is very noisy, which is why
we band-average — see {ref}`spec_ex2` for why that step is not a cosmetic one.

```{solution-end}
```

```{exercise-start}
:label: spec_ex2
```

**Why the periodogram must be smoothed.** The band-averaging in part (c) above is not a
refinement of the estimator; without it there is no consistency at all. This exercise makes
that concrete.

Simulate the discretely sampled Ornstein–Uhlenbeck process ($a=1$, $b=0.7$, unit sampling),
which is the AR(1) $x_t = e^{-a}x_{t-1} + \eta_t$ with
$\sigma_\eta^2 = \frac{b^2}{2a}(1-e^{-2a})$, and form the periodogram
$I_N(w) = N^{-1}|\sum_t (x_t - \bar x)e^{-iwt}|^2$.

(a) Fix an interior frequency and, over many independent replications, tabulate the mean and
standard deviation of $I_N(w)$ for $N = 128, 512, 2048$. The *mean* sits close to $S(w)$ at
every $N$ — the periodogram is asymptotically unbiased, and at this replication count the
residual bias is not easily separated from Monte Carlo error. The point of the table is the
other column: the *standard deviation* stays put at roughly $S(w)$ itself, however large $N$
becomes.

(b) Explain the result. $I_N(w)$ is the squared modulus of a *single* Fourier coefficient, and
that coefficient is asymptotically complex Gaussian with variance $S(w)$ — a fixed number of
random quantities, however long the record. Hence $I_N(w)/S(w)$ is asymptotically
$\tfrac12\chi^2_2$, with mean $1$ and variance $1$, forever. Lengthening the record buys more
*frequencies*, not more precision at any one of them.

(c) Now average over the $m$ frequencies nearest $w$ and watch the standard deviation fall like
$m^{-1/2}$. Consistency requires letting $m$ grow with $N$ while the bandwidth $m/N$ shrinks —
the conditions $b_N \to 0$ and $Nb_N \to \infty$ of the spectral-window literature.

```{exercise-end}
```

```{solution-start} spec_ex2
:class: dropdown
```

```{code-cell} ipython3
rng = np.random.default_rng(0)
a, b = 1.0, 0.7
phi = np.exp(-a)
sig2 = (b**2/(2*a))*(1 - phi**2)

def periodograms(N, reps):
    """reps independent periodograms of an AR(1) record of length N."""
    out = np.empty((reps, N//2 + 1))
    for r in range(reps):
        e = rng.normal(0, np.sqrt(sig2), N + 500)
        x = np.zeros(N + 500)
        for t in range(1, N + 500):
            x[t] = phi*x[t-1] + e[t]
        x = x[500:]; x = x - x.mean()
        out[r] = np.abs(np.fft.rfft(x))**2 / N
    return out

print("   N     mean I     true S     s.d. I    s.d./S")
for N in [128, 512, 2048]:
    I = periodograms(N, 600)
    k = N//8                                   # a fixed interior frequency
    w = 2*np.pi*k/N
    Strue = sig2/np.abs(1 - phi*np.exp(-1j*w))**2
    print(f"{N:5d}   {I[:,k].mean():8.5f}   {Strue:8.5f}   {I[:,k].std():8.5f}   "
          f"{I[:,k].std()/Strue:6.3f}")
```

```{code-cell} ipython3
# (c) smoothing over m neighbouring ordinates: the s.d. falls like 1/sqrt(m)
N = 2048
I = periodograms(N, 600)
k = N//8
w = 2*np.pi*k/N
Strue = sig2/np.abs(1 - phi*np.exp(-1j*w))**2
print("   m   mean of smoothed I   s.d.    s.d./S   (1/sqrt(m) for comparison)")
for m in [1, 5, 25, 125]:
    sm = I[:, k-m//2 : k+m//2+1].mean(axis=1)
    print(f"{m:4d}   {sm.mean():14.5f}   {sm.std():7.5f}  {sm.std()/Strue:6.3f}"
          f"      {1/np.sqrt(m):6.3f}")
```

The first table is the point: the mean tracks $S(w)$ but the standard deviation does not fall,
staying at roughly $S(w)$ as $N$ increases sixteenfold. The periodogram is asymptotically
unbiased and permanently noisy. The second table shows smoothing doing the work that
lengthening the record cannot: averaging $m$ neighbouring ordinates cuts the standard deviation
by very nearly $m^{-1/2}$, as it would for $m$ independent observations — which is what
neighbouring periodogram ordinates asymptotically are.

That last clause is where the hidden assumption sits. Asymptotic independence across
frequencies is a statement about *fourth* moments, not second, so nothing in this chapter
establishes it; see {doc}`/appendices/ergodicity`. The same appendix explains why
consistency of the $\hat R(\tau)$ feeding the periodogram is itself a fourth-moment question.

```{solution-end}
```
