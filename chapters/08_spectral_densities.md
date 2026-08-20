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

```{eval-rst}
.. index::
   single: spectral density; definition
   single: Fourier transform; of an autocovariance
   single: Laplace transform
   single: cross-spectral density; definition
```

## (a) General Results

Let $x(t)$ be a covariance stationary stochastic process with covariogram

$$
Ex(t)x(t-\tau) = R(\tau), \qquad \tau \in \mathbb{R},
$$

where $R$ is a positive semidefinite function. *Covariogram* names the same function that
{doc}`01_covariance_stationary_processes` calls the autocorrelation function $R(\tau)$ and, for
a zero mean process, the autocovariance function $C(\tau)$. All three names appear in this book,
and the note in that chapter explains why. The *power spectrum* or the *spectral
density* of $x$ is defined as

```{math}
:label: eq-8-1
S(\omega) = \int_{-\infty}^{\infty} e^{-i\omega\tau}\, R(\tau)\, d\tau,\ \omega \in (-\infty,\, \infty).
```

Equation {eq}`eq-8-1` defines $S(\omega)$ as the Fourier transform of $R(\tau)$. Given $S(\omega)$,
$R(t)$ can be recovered from the inverse Fourier transform

```{math}
:label: eq-8-2
R(\tau) = \frac{1}{2\pi}\ \int_{-\infty}^{\infty} e^{+i\omega\tau}\, S(\omega)\, d\omega
```

From the fact that $R(\tau)$ is a positive semidefinite function, which implies that
$R(\tau) = R(-\tau)$, definition {eq}`eq-8-1` and the inverse relation {eq}`eq-8-2` imply
that $S(\omega)$ has the following properties

(i) $S(\omega) = S(-\omega)$. (from definition {eq}`eq-8-1`)

(ii) $S(\omega) \geq 0$. (from definition {eq}`eq-8-1` and positive semidefiniteness of $R(t)$)

(iii) $\frac{1}{2\pi} \int_{-\infty}^{\infty} S(\omega)\, d\omega = R(0)\quad$ ({eq}`eq-8-2` evaluated
at $\tau = 0$)

Property (iii) asserts that the spectral density achieves a decomposition of the variance
$R(0)$ by frequency. We shall shortly see that this decomposition of variance is actually
into components that are orthogonal across different frequencies.

For the purposes of using {eq}`eq-8-1` and {eq}`eq-8-2`, it is fortunate that tables of
Fourier transform pairs have been prepared. Table 1 is a small table. The reader can verify
the entries in the table by using {eq}`eq-8-1` or {eq}`eq-8-2`.

Fourier transforms have a number of useful operational properties. In particular, let
$F(\omega) \leftrightarrow f(t)$ mean that $F(\omega)$ is the Fourier transform of the "time signal"
$f(t)$, i.e.,

$$
F(\omega) = \int_{-\infty}^{\infty} f(t)\, e^{-i\omega t}\, dt.
$$

Then a number of simple operations of $F(\omega)$ or $f(t)$ can be used to generate other
Fourier transform pairs. Some of these are recorded in table 2.

Property (8) of table 2 can be used to find the spectral density of mean square
derivatives of a given process $x(t)$. We saw above that if $x_t$ is a covariance
stationary stochastic process with autocovariogram $R(\tau)$, then the autocovariogram of
$\frac{d^n}{dt^n}\ x(t)$ is $(-1)^n\ \frac{d^{2n} R(\tau)}{d\tau^{2n}}$. It follows from
property 8 of table 2 that the Fourier transform of
$(-1)^n \ \frac{d^{2n}R(\tau)}{d\tau^{2n}}$ is $\omega^{2n}F(\omega)$, where $F(\omega)$ is the Fourier
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
S_y (\omega) = h(\omega) S_x (\omega) h(-\omega)
```

where $h(\omega) = \int_{-\infty}^{\infty} h (\tau) e^{-i \omega \tau}\, d\tau$ is the Fourier
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
\int_{-\infty}^{\infty} \delta' (t) e^{-i\omega t}\, dt =\ -\ \frac{d}{dt}\ \ e^{-i\omega t}\Big|_{t=0} = i\omega.
```

More generally, the Fourier transform of $\delta^{(n)} (t)$ is given by

$$
\delta^{(n)} (t) \leftrightarrow (i\omega)^n.
$$

Using this result in conjunction with {eq}`eq-8-5` with $h(t) = \delta^{(n)} (t)$
immediately gives that the spectrum of the $n^{th}$ mean square derivative of a process
$x(t)$ with spectrum $S(\omega)$ is given by $\omega^{2n} S(\omega)$.

We record some of these results in table 3.

```{eval-rst}
.. index::
   single: transform tables
   single: Fourier transform; table of
```

## Tables

### Table 1

$$
F(\omega) = \int_{-\infty}^{\infty} f(t) e^{-i\omega t}\, dt
$$

| $f(t)$ | $F(\omega)$ |
| --- | --- |
| $e^{-a\vert t\vert}$ | $\frac{2a}{a^2 + \omega^2}$, $a > 0$ |
| $u(t)$ | $\frac{1}{i\omega} + \pi\, \delta(\omega)$ |
| $u(t)t$ | $\frac{1}{(i\omega)^2} + i\pi\, \delta'(\omega)$ |
| $\delta (t)$ | $1$ |
| $\delta' (t)$ | $i\omega$ |
| $\delta^{(n)}\, (t)$ | $(i\omega)^n$ |
| $\sum_{n=-\infty}^{\infty} \delta (\tau-nT)$ | $\omega_0 \sum_{n=-\infty}^{\infty} \delta (\omega-n \omega_0),\ \omega_0 = \frac{2\pi}{T}$ |
| $e^{-at}\, u(t)$ | $\frac{1}{a + i\omega}$, $a > 0$ |
| $\begin{cases} 1 & t \in [0,\,1] \\ 0 & t > 1,\ t < 0 \end{cases}$ | $\left(\frac{1 - e^{-i\omega}}{i\omega}\right)$ |
| $\begin{cases} 1 & t \in [0,\,1] \\ -1 & t \in [1,\, 2] \\ 0 & t > 2,\ t < 0 \end{cases}$ | $\frac{(1 - e^{-i\omega})^2}{i\omega}$ |
| $\begin{cases} t, & t \in [0,\,1] \\ 2-t, & t \in [1,\,2] \\ 0 & t < 0,\ t > 2 \end{cases}$ | $\frac{(1-e^{-i\omega})^2}{(i\omega)^2}$ |
| $\frac{t^{n-1}}{(n-1)!}\, e^{-at}\, u(t)$ | $\frac{1}{(a+i\omega)^n}$, $a > 0$ |
| $e^{i\omega_0t}$ | $2\pi \delta (\omega-\omega_0)$ |
| $\cos\ \omega_0 t$ | $\pi\, \left[\delta (\omega-\omega_0) + \delta\, (\omega+\omega_0)\right].$ |

Note: $u(t) = \begin{cases} 1 & t \geq 0 \\ 0 & t < 0 \end{cases}\qquad$ (Heaviside step
function)

The entries for $u(t)$ and $u(t)t$ are the transforms in the sense of generalized functions;
the $\delta$ terms record that these signals do not decay, and they carry all of the mass at
$\omega = 0$. In every use we make of these entries, the operators $1/(i\omega)$ and $1/(i\omega)^2$ act on a
noise that has no mass at $\omega = 0$, so the $\delta$ terms may be, and are, dropped. Those uses are the
operational calculus of {doc}`11_linear_sde` and {doc}`12_prediction` and the nonstationary
examples of {doc}`14_nonstationary_examples`. Writers who work only
with such operators often tabulate $u(t) \leftrightarrow 1/i\omega$ for that reason.

$\delta (t) =$ Dirac delta generalized function, defined by
$\int_{-\infty}^{\infty} g(t) \delta(t)\, dt = g(0)$ where $g(t)$ is a "test function" that
is continuous at $t=0$.

### Table 2

$$
f(t) \leftrightarrow F(\omega)\ \text{ means }\ F(\omega) = \int_{-\infty}^{\infty} f(t) e^{-i\omega t}\, dt.
$$

Property:

| Property | Transform pair |
| --- | --- |
| 1. Linearity | $a_1 f_1 (t) + a_2 f_2(t) \leftrightarrow a_1 F_1(\omega) + a_2 F_2 (\omega)$; $a_1,\ a_2$ are scalars. |
| 2. Symmetry | $F(t) \leftrightarrow 2\pi f(-\omega)$ |
| 3. Scaling | $f(at) \leftrightarrow \frac{1}{\vert a\vert}\ F\left(\frac{\omega}{a}\right),\ \text{ a scalar}$ |
| 4. Delay | $f(t-t_0) \leftrightarrow e^{-i\omega t_0}\, F(\omega)$ |
| 5. Modulation | $e^{i\omega_0t}\, f(t) \leftrightarrow F(\omega-\omega_0)$ |
| 6. Convolution | $f_1(t) \ast f_2(t) \leftrightarrow F_1(\omega)F_2(\omega)$; where $f_1(t) \ast f_2(t) \equiv \int_{-\infty}^{\infty} f_1(t-\tau) f_2(\tau)\, d\tau$ |
| 7. Multiplication | $f_1(t) f_2 (t) \leftrightarrow \frac{1}{2\pi} F_1 (\omega) \ast F_2 (\omega)$; where $F_1(\omega) \ast F_2(\omega) \equiv \int_{-\infty}^{\infty} F_1(\omega-s) F_2(s)\, ds$ |
| 8. Time differentiation | $\frac{d^n}{dt^n}\, f(t) \leftrightarrow (i\omega)^n\ F(\omega)$ |
| 9. Time integration | $\int_{-\infty}^{t} f(\tau)\, d\tau \leftrightarrow\ \frac{F(\omega)}{i\omega} + \pi\ F(0)\, \delta (\omega)$ |
| 10. Frequency differentiation | $-itf(t)\ \leftrightarrow\ \frac{dF(\omega)}{d\omega}$ |
| 11. Frequency integration | $\frac{f (t)}{-it}\ \leftrightarrow\ \int F(\omega') d\omega'$ |
| 12. Reversal | $f(-t) \leftrightarrow F(-\omega)$ |

### Table 3

| Process | Autocovariogram | Spectrum |
| --- | --- | --- |
| $x_t$ | $R(\tau)$ | $S(\omega)$ |
| $\int_{-\infty}^{\infty} h(\tau) x(t-\tau)\, d\tau$ | $h\, \ast\, R\, \ast\, h(-s)$ | $\vert h(\omega)\vert^2 S(\omega)$ |
| $\frac{d}{dt}\ x(t)$ | $-\ \frac{d^2 R(\tau)}{d\tau^2}$ | $\omega^2 S(\omega)$ |
| $\frac{d^n}{dt^n}\, x(t)$ | $(-1)^n\ \frac{d^{2n}\, R(\tau)}{d\tau^{2n}}$ | $\omega^{2n} S(\omega)$ |

### Table 4

| $x(t)$ | $R(\tau)$ | $S(\omega)$ |
| --- | --- | --- |
| $w(t)$ | $\delta (\tau)$ | $1$ |
| $(D - \lambda) x(t) = w(t),\, \lambda < 0$ | $\frac{-1}{2\lambda}\ e^{\lambda \vert\tau\vert}$ | $\frac{1}{\lambda^2 + \omega^2}$ |
| $\theta (D) x(t) = w(t)$ | | $\frac{1}{\theta (i\omega)\, \theta (-i\omega)}$ |
| or | | or |
| $(D - \lambda_1)\ldots (D-\lambda_n)\, x(t) = w(t)$ | $\sum_{j=1}^{n}\, k_j\, e^{\lambda_j \vert\tau\vert}$ | $\frac{1}{\Pi_{j=1}^{n}\, (i\omega-\lambda_j)\, (-i\omega-\lambda_j)}$ |
| $re\, (\lambda_j) < 0,\ j=1,\, \ldots,\, n$ | $k_j = \lim_{s \to \lambda_j}\ \frac{(s - \lambda_j)}{\theta (s)\, \theta (-s)}$ | |
| $\theta (D) x(t) = \Psi(D) w(t)$ | $\Sigma\, k_j\, e^{\lambda_j \vert\tau\vert}$ | |
| where $\theta (D) = (D-\lambda_1) \ldots (D-\lambda_n)$ | | |
| $re (\lambda_j) < 0$ | $k_j = \lim_{s\to \lambda_j}\ (s-\lambda_j)\ \frac{\Psi (s)}{\theta (s)}\ \frac{\Psi (-s)}{\theta (-s)},$ | $\frac{\Psi (i\omega)\, \Psi (-i\omega)}{\theta (i\omega)\, \theta(-i\omega)}$ |
| $\Psi (D) = \Psi_0 + \Psi_1 D + \ldots + \Psi_m \, D^m$ | | |
| $m < n$ | | |

```{eval-rst}
.. index::
   single: Wold decomposition; continuous time
   single: Wold's theorem
   single: moving average representation; one-sided
   single: linearly deterministic process; definition
   single: linearly indeterministic process; definition
```

## (b) Wold's Theorem

We now state a version of Wold's decomposition theorem in continuous time.

**Theorem 10.** Let $x(t)$ be a covariance stationary stochastic process with autocovariance
function $R(\tau)$ and spectral density $S(\omega)$. Then $x(t)$ can be represented as

```{math}
:label: eq-8-a
x(t) = \int_0^{\infty} p(\tau)\omega (t-\tau)\, d\tau + x^d(t)
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

We restrict the deterministic component to the harmonic case, in which $x^d(t)$ is a sum of
random amplitude sine and cosine waves at a countable set of fixed frequencies. The general Wold
decomposition permits any linearly deterministic $x^d$, including one whose spectral
distribution is singular but carries no atoms. Every use we make of $x^d$ needs only the
harmonic case: the mean square ergodicity criterion of {doc}`10_cramer_representation` and the
seasonal example of {doc}`17_discrete_sampling_folding` both turn on whether one of the
frequencies below is zero. With that restriction, the spectral density of $x^d (t)$ can be
represented as

```{math}
:label: eq-8-b
S_d (\omega) = \sum_{j=1}^{\infty} a_j \pi \left[\delta (\omega-\omega_j) + \delta (\omega + \omega_j)\right]
```

where $a_j$ are positive constants and where $\omega_j,\ j = 1,\ \ldots$ is a countable set of
frequencies. Therefore, the autocovariance function of $x^d (t)$ is given by

```{math}
:label: eq-8-c
R^d (\tau) = \sum_{j=1}^{\infty} a_j\ \cos\, (\omega_j \tau).
```

It follows from {eq}`eq-8-a`, {eq}`eq-8-b` and the convolution property {eq}`eq-8-5` that the
spectral density of $x(t)$ can be represented as

$$
S (\omega) = P(\omega) P(-\omega) + \sum_{j=1}^{\infty} a_j \pi \left[\delta (\omega-\omega_j) + \delta (\omega + \omega_j)\right]
$$

where $P(\omega) = \int_0^{\infty} p(\tau) e^{-i\omega\tau}$, is the Fourier transform of a one-sided,
square-integrable function. The component $\int_0^{\infty} p(\tau) w(t - \tau)\, d\tau$ is
called the *linearly indeterministic* part of the process $x(t)$.

```{eval-rst}
.. index::
   single: spectral factorization theorem; statement
   single: fundamental white noise; definition
   single: non-fundamentalness; definition
   single: transfer function; definition
   single: rational spectral density; factorization of
   single: minimum phase
   single: right half plane zeros
```

## (c) The Spectral Factorization Theorem

Included in the statement of Wold's theorem is the spectral factorization theorem for
linearly indeterministic processes. We restate this property separately:

**Spectral factorization theorem.** Let $x(t)$ be a covariance stationary, *linearly
indeterministic* process with spectral density $S(\omega)$. Then $S(\omega)$ can be factored as

$$
S(\omega) = \tilde P (i\omega) \tilde P(-i\omega)
$$

where $\tilde P (s) = \int_0^{\infty} p(t) e^{-st}\, dt$, where
$\int_0^{\infty} p(t)^2\, dt < + \infty$, and $\tilde P (s)$ has no zero for $s$ in the
*open* right half of the complex plane. This condition on the zeros of $\tilde P(s)$ is the
condition that the white noise $w(t)$ in
$x(t) = \int_0^{\infty} p(\tau) w(t-\tau)\, d\tau$ is fundamental for $x(t)$.

The word *open* matters. A zero of $\tilde P(s)$ strictly inside the right half plane destroys
fundamentalness, as the example below shows. A zero *on* the imaginary axis is a boundary case:
the spectral density vanishes at that frequency, the inverse filter $1/\tilde P(s)$ fails to be
square integrable there, and $\omega$ recovers $x$ only as a limit. We admit such factors, and
{doc}`09_characterizations_ms_differentiability` uses one when it differentiates a Wold
representation. The function
$\tilde P(s)$ is the Laplace transform of $p(\tau)$. The Fourier transform $P(\omega)$ is related
to $\tilde P(s)$ by $P(\omega) = \tilde P (i\omega)$.

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
S(\omega) = \tilde R (i\omega) \tilde R (-i\omega).
$$

However, this is not the representation alluded to in the statement of Wold's theorem,
because $\tilde R(s)$ has a zero at $s = b > 0$, which is in the right half plane. This
reflects the fact that the space $H_v(- \infty, \, t)$ is strictly larger than
$H_x(- \infty,\,t)$. This right-half-plane zero is the continuous-time prototype of
*non-fundamentalness*, in which the driving noise spans a larger information space than the
observable process. It returns as the central difficulty in interpreting vector autoregressions
in {doc}`18_time_aggregation_var` and as the identification problem of
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
S(\omega) &= \frac{(i\omega-b)\ (-i\omega-b)}{(i\omega-\lambda_1)\ (i\omega-\lambda_2)\ (-i\omega-\lambda_1)\ (-i\omega-\lambda_2)} \\
&= \frac{(i\omega+b)\ (-i\omega+b)}{(i\omega-\lambda_1)\ (i\omega-\lambda_2)\ (-i\omega-\lambda_1)\ (-i\omega-\lambda_2)}.
\end{aligned}
$$

That is, the spectral density remains unaltered if we simply change the sign of the real
part of the zero of the numerator polynomial of $\tilde R (i\omega)$, so that $-b$ is replaced
by $+b$. It follows that the spectral density can be represented as

$$
S(\omega) = \tilde P (i\omega) \tilde P (-i\omega)
$$

where

$$
\tilde P (i\omega) = \frac{i\omega + b}{(i\omega-\lambda_1)\ (i\omega-\lambda_2)},
$$

where $\tilde P (s)$ now satisfies the hypotheses required in the statement of Wold's
theorem, in particular that $\tilde P(s)$ have no zeroes in the right half plane, and that
$\tilde P (s)$ be the Laplace transform of a square summable function $p(\tau)$ with support
on $[0,\, \infty)$. It follows that a Wold representation for $x(t)$ is

$$
x(t) = \int_0^{\infty} p(\tau) w(t-\tau)\, d\tau
$$ (eq-8-wold)

where $w(t)$ is a fundamental white noise for $x(t)$, and $p(\tau)$ is the inverse transform
of $\tilde P (i\omega) = (i\omega+b)/(i\omega-\lambda_1)\ (i\omega-\lambda_2)$. By a partial fractions
representation of $\tilde P (s)$, it follows that

$$
p(\tau) = g_1 e^{\lambda_1 \tau} + g_2 e^{\lambda_2 \tau}\quad \tau \geq 0
$$

where
$g_j = \lim_{s\to \lambda_j}\ (s - \lambda_j)\ \frac{(s + b)}{(s - \lambda_1)\ (s - \lambda_2)}$.

It is a good exercise to invert {eq}`eq-8-wold` and so express $w(t)$ as a sum of square
summable integrals of past values of $x(t),\ Dx(t)$, and $D^2 x(t)$.

```{eval-rst}
.. index::
   single: periodogram; definition
   single: spectral smoothing; bandwidth
```

## Exercises

These exercises use the Ornstein–Uhlenbeck process of Chapter 7,

$$
dx(t) = -a\,x(t)\,dt + b\,dW(t), \qquad a, b > 0,
$$

as a running example. From Table 4 (the row $(D-\lambda)x = \omega$ with $\lambda = -a$, scaled
by $b$), its autocovariance and spectral density are

$$
R(\tau) = \frac{b^2}{2a}\, e^{-a|\tau|}, \qquad
S(\omega) = \frac{b^2}{a^2 + \omega^2}.
$$

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt
```

```{exercise-start}
:label: spec_ex1
```

Take $a = 1$, $b = 0.7$.

(a) Plot the spectral density $S(\omega) = b^2/(a^2 + \omega^2)$.

(b) Verify property (iii) of the text numerically: the spectral density decomposes the
variance, $\frac{1}{2\pi}\int_{-\infty}^{\infty} S(\omega)\,d\omega = R(0) = b^2/(2a)$.

(c) Simulate a long OU path, form its **periodogram** (the sample analogue of the spectral
density), and check that, after smoothing, it tracks $S(\omega)$.

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
ax.set_xlabel(r'$\omega$'); ax.set_ylabel(r'$S(\omega)$')
ax.set_title(r'Spectral density of the OU process, $S(\omega)=b^2/(a^2+\omega^2)$')
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
ax.plot(ww, S(ww), 'r-', lw=2, label=r'$S(\omega)=b^2/(a^2+\omega^2)$')
ax.set_xlabel(r'$\omega$'); ax.set_ylabel('spectral density')
ax.legend()
plt.show()
```

The numerical integral of $S(\omega)/2\pi$ reproduces $R(0) = b^2/(2a)$, and the smoothed
periodogram of the simulated path follows the theoretical spectrum: most of the variance
sits at low frequencies, with a $1/\omega^2$ tail. The raw periodogram is very noisy, which is why
we band-average. {ref}`spec_ex2` shows why that step is not cosmetic.

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
$I_N(\omega) = N^{-1}|\sum_t (x_t - \bar x)e^{-i\omega t}|^2$.

(a) Fix an interior frequency and, over many independent replications, tabulate the mean and
standard deviation of $I_N(\omega)$ for $N = 128, 512, 2048$. The *mean* sits close to $S(\omega)$ at
every $N$. The periodogram is asymptotically unbiased, and at this replication count the
residual bias is not easily separated from Monte Carlo error. The point of the table is the
other column: the *standard deviation* stays put at roughly $S(\omega)$ itself, however large $N$
becomes.

(b) Explain the result. $I_N(\omega)$ is the squared modulus of a *single* Fourier coefficient, and
that coefficient is asymptotically complex Gaussian with variance $S(\omega)$, a fixed number of
random quantities however long the record. Hence $I_N(\omega)/S(\omega)$ is asymptotically
$\tfrac12\chi^2_2$, with mean $1$ and variance $1$, forever. Lengthening the record buys more
*frequencies*, not more precision at any one of them.

(c) Now average over the $m$ frequencies nearest $\omega$ and watch the standard deviation fall like
$m^{-1/2}$. Consistency requires letting $m$ grow with $N$ while the bandwidth $m/N$ shrinks. Those are
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

The first table is the point: the mean tracks $S(\omega)$ but the standard deviation does not fall,
staying at roughly $S(\omega)$ as $N$ increases sixteenfold. The periodogram is asymptotically
unbiased and permanently noisy. In the second table, averaging $m$ neighbouring ordinates cuts
the standard deviation by very nearly $m^{-1/2}$, as it would for $m$ independent observations.
Neighbouring periodogram ordinates are asymptotically independent.

That last clause is where the hidden assumption sits. Asymptotic independence across
frequencies is a statement about *fourth* moments, not second, so nothing in this chapter
establishes it; see {doc}`/appendices/ergodicity`. The same appendix explains why
consistency of the $\hat R(\tau)$ feeding the periodogram is itself a fourth-moment question.

```{solution-end}
```
