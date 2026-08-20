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

```{eval-rst}
.. index::
   single: Cramér representation; statement
   single: random spectral measure; definition
   single: spectral distribution function
   single: frequency domain; representation of a process
```

# 10. The Cramér Representation

```{eval-rst}
.. index::
   single: Cramér representation
   single: random spectral measure
   single: orthogonal increments
   single: band-pass filter
```

The spectral density accomplishes an orthogonal decomposition by frequency of the variance of a covariance stationary stochastic process. The Cramér representation exhibits this decomposition in a general way.

Let $x(t)$ be a covariance stationary, zero mean stochastic process with autocorrelation function $R(\tau)$. Let the spectral density of $x$ be denoted $S(\omega)$, where

$$
S(\omega) = \int^\infty_{-\infty} R(\tau) e^{-i\omega\tau} d\tau.
$$

Associated with $x(t)$, there exists a complex valued random process $Z(\omega)$ defined by

(i) $Z(0) = 0$

(ii) $Z(-\omega) = \overline{Z(\omega)}$, where the bar denotes complex conjugation.

(iii) $Z(\omega)$ is a process with orthogonal increments, and in particular

$$
E Z'(\omega) \overline{Z'(\nu)} = S(\nu) \delta(\nu - \omega),
$$

where the prime denotes the (possibly generalized) derivative of $Z(\omega)$ with respect to $\omega$. The process $Z(\omega)$ is called the "random spectral measure" of the $x(t)$ process. This orthogonal-increments random measure reappears as the foundational object $W$ of the Hansen–Sargent prediction calculus in {doc}`19_prediction_formulas_continuous_time`. In terms of this random process, the $x(t)$ process has the Cramér representation

$$
x(t) = \frac{1}{\sqrt{2\pi}}\, \int^\infty_{-\infty} e^{i\omega t} dZ(\omega)
$$

or

```{math}
:label: eq-10-cramer
x(t) = \frac{1}{\sqrt{2\pi}} \int^\infty_{-\infty} e^{i\omega t}\, Z'(\omega) d\omega.
```

To help motivate this representation, we use {eq}`eq-10-cramer` to calculate $R(\tau) = Ex(t) x(t-\tau)$,

$$
\begin{aligned}
Ex(t) x(t-\tau) &= E\, \frac{1}{2\pi} \int^\infty_{-\infty} e^{i\omega t}\, Z' d\omega \int^\infty_{-\infty} e^{-i\nu (t-\tau)}\, \overline{Z'(\nu)} d\nu \\
&= \frac{1}{2\pi} \int^\infty_{-\infty} e^{-i\nu (t - \tau)} \int^\infty_{-\infty} e^{i\omega t}\, EZ'(\omega)\, \overline{Z'(\nu)}\, d\omega d\nu \\
&= \frac{1}{2\pi} \int^\infty_{-\infty} e^{-i\nu (t - \tau)} \int^\infty_{-\infty} e^{i\omega t}\, S(\nu) \delta(\nu - \omega) d\omega d\nu \\
&= \frac{1}{2\pi} \int^\infty_{-\infty} e^{-i\nu (t - \tau)}\, e^{i\nu t}\, S(\nu) d\nu
\end{aligned}
$$

or

$$
R(\tau) = \frac{1}{2\pi} \int^\infty_{-\infty} e^{i\nu \tau} S(\nu) d\nu.
$$

This is the inversion formula {eq}`eq-8-2` for recovering the autocorrelation function from the spectral density.

There is a sense in which the random spectral measure can be defined formally as the Fourier transform of $x(t)$,

$$
Z'(\omega,\, w) = \frac{1}{\sqrt{2\pi}}\, \int^\infty_{-\infty} e^{-i\omega t}\, x(t,\, w) dt,
$$ (eq-10-Zhat)

provided that the integral is interpreted delicately. In {eq}`eq-10-Zhat`, we have added the argument $w \in \Omega$ explicitly to emphasize that both $x(t,\, w)$ and $Z'(\omega,\, w)$ are random processes defined on the same underlying probability space $(\Omega,\, \mathcal{F},\, P)$.

Differentiating {eq}`eq-10-cramer` formally with respect to time, we have that the mean square derivative of $x(t)$, if it exists, has Cramér representation

$$
x'(t) = \frac{1}{\sqrt{2\pi}} \int^\infty_{-\infty} i\, \omega\, e^{i\omega t}\, Z'(\omega) d\omega.
$$

The derivative of the random spectral measure associated with $x'(t)$ is thus $\omega Z'(\omega)$, which being proportional to $Z'(\omega)$ is itself the derivative of a process with orthogonal increments and obeys

$$
E \omega Z'(\omega)\ \overline{\nu Z'(\nu)} = \nu^2 S(\nu) \delta(\nu - \omega).
$$

Thus the spectral density of $Dx(t)$ is $\nu^2 S(\nu)$. More generally, consider the distributed lag

$$
y(t) = \int^\infty_{-\infty} b(\tau) x(t - \tau) d\tau
$$

where $b(\tau) \in L_2\, (-\infty,\, \infty)$. Then using {eq}`eq-10-cramer`, we have

```{math}
:label: eq-10-1
\begin{aligned}
y(t) &= \frac{1}{\sqrt{2\pi}} \int^\infty_{-\infty} b(\tau) \int^\infty_{-\infty} e^{i\omega (t-\tau)}\, Z'(\omega) d\omega \\
&= \frac{1}{\sqrt{2\pi}} \int^\infty_{-\infty} \int^\infty_{-\infty} b(\tau)\, e^{-i\omega \tau}\, d\tau\, e^{i\omega t}\, Z'(\omega) d\omega \\
y(t) &= \frac{1}{\sqrt{2\pi}} \int^\infty_{-\infty} e^{i\omega t}\, B(\omega)\, Z'(\omega) d\omega
\end{aligned}
```

where $B(\omega)$ is the Fourier transform of $b(t)$, namely

$$
B(\omega) = \int^\infty_{-\infty} b(\tau) e^{-i\omega \tau}\, d\tau.
$$

From {eq}`eq-10-1`, it follows that the random spectral measure of $y(t)$ has derivative $B(\omega) Z'(\omega)$, and that $y(t)$ has spectral density

$$
B(\omega)\ \overline{B(\omega)}\ S(\omega)
$$

or

```{math}
:label: eq-10-2
S_y(\omega) = |B(\omega)|^2 S_x(\omega)
```

where $S_x(\omega) \equiv S(\omega)$ is the spectral density of $x(t)$.

Equation {eq}`eq-10-2` can be used to view from another angle the orthogonal decomposition across frequencies induced by $S_x(\omega)$. For $0 < c < d$, we create the "window" $B_{cd}(\omega)$ according to

$$
B_{cd}(\omega) = \begin{cases} 1 & \omega \in [c,\, d]\ \text{ or }\ [-d,\,-c] \\ 0 & \text{otherwise} \end{cases}
$$

$$
\text{where }\ d > c > 0.
$$

The band-pass filter $B_{cd}(\omega)$ is illustrated in {numref}`fig-10-1`.

```{figure} figures/fig-10-1_bandpass_window.png
:name: fig-10-1
:width: 90%
:align: center

Figure 1. The band-pass filter ("window") $B_{cd}(\omega)$, defined on the frequency axis $\omega \in (-\infty,\, \infty)$. It equals $1$ on the two symmetric frequency bands $[c,\, d]$ and $[-d,\, -c]$ and equals $0$ everywhere else, with $d > c > 0$.
```

Define the time function

$$
\begin{aligned}
b_{cd}(t) &= \frac{1}{2\pi} \int^\infty_{-\infty} B_{cd}(\omega) e^{+i\omega t}\, d\omega \\
&= \frac{1}{\pi}\ \left[ \frac{\sin\ dt}{t}\ - \ \frac{\sin\ ct}{t}\right]
\end{aligned}
$$

We have that $y_{cd}(t)$, defined by

$$
y_{cd}(t) \equiv \int^\infty_{-\infty} b_{cd}(\tau) x(t-\tau) d\tau,
$$

has the spectral density $S_{cd}$ given by

$$
S_{cd}(\omega) = \begin{cases} S_x(\omega) & \omega \in [c,\, d]\ \text{ or }\ [-d,\, -c] \\ 0 & \text{otherwise} \end{cases}
$$

If we choose an interval parameterized by $0 < e < f$ such that $[e,\,f] \cap [c,\, d] = 0$, it follows from the orthogonal increments property of $Z(\omega)$ that $y_{ef}(t)$ is orthogonal to $y_{cd}(t-\tau)$ for all $t$. For using the Cramér representation, we have

$$
\begin{aligned}
Ey_{cd}(t) y_{ef}(t-\tau) &= \frac{1}{2\pi} \int^\infty_{-\infty} e^{i\omega t}\, B_{cd}(\omega) Z'(\omega) \\
&\quad \int^\infty_{-\infty} e^{-i\nu (t-\tau)}\, \overline{B_{ef}(\nu)}\, Z'(\nu) d\nu d\omega \\
&= \frac{1}{2\pi} \int^\infty_{-\infty} e^{-i\nu (t - \tau)}\, \overline{B_{ef}(\nu)} \\
&\quad \int^\infty_{-\infty} B_{cd}(\omega) S(\nu) \delta(\nu - \omega) d\omega d\nu \\
&= \frac{1}{2\pi} \int^\infty_{-\infty} e^{-i\nu (t - \tau)}\, \overline{B_{ef}(\nu)}\ B_{cd}(\nu) S(\nu) d\nu \\
&= 0
\end{aligned}
$$

since $\overline{B_{ef}(\nu)}\ B_{cd}(\nu) = 0$ for all $\nu$.

We also have that $y_{ef}(t)$ has spectral density

$$
S_{ef}(\omega) = \begin{cases} S_x(\omega) & \omega \in [e,\, f]\ \text{ or }\ [-f,\, -e] \\ 0 & \text{otherwise} \end{cases}
$$

By filtering the process $x(t)$ with filters formed by taking disjoint frequency intervals, $[c,\,d]$, $[e,\,f]$, we create processes that are orthogonal at all leads and lags, but whose respective spectral densities equal that of $x(t)$ over the intervals $[c,\,d]$ and $[e,\,f]$.

To motivate from another angle the frequency decomposition induced by the spectral density, it is useful to note that the Cramér representation

$$
x(t) = \frac{1}{\sqrt{2\pi}} \int^\infty_{-\infty} e^{i\omega t}\, Z'(\omega) d\omega
$$

can be also represented in a couple of alternative ways. Let the complex valued random process $Z'(\omega)$ be represented in the alternative forms

$$
\begin{aligned}
z'(\omega) &= a(\omega) + b(\omega)i \\
z'(\omega) &= r(\omega)e^{i\theta(\omega)}.
\end{aligned}
$$

Then it is straightforward to derive the Cramér representation in the forms

```{math}
:label: eq-10-3
x(t) = \sqrt{\frac{2}{\pi}} \int^\infty_0\ a(\omega)\ \cos\ \omega t\, d\omega \ - \ \sqrt{\frac{2}{\pi}} \int^\infty_0 b(\omega)\ \sin\ \omega t\, d\omega
```

where

$$
E (a(\omega) + ib(\omega))\ (a(\nu) - ib(\nu)) = S(\nu) \delta(\nu - \omega);
$$

and

```{math}
:label: eq-10-4
x(t) = \sqrt{\frac{2}{\pi}} \int^\infty_0 r(\omega)\ \cos\ (\omega t + \theta(\omega)) d\omega
```

where

$$
Er(\omega) e^{i\theta(\omega)}\, r(\nu) e^{-i\theta(\nu)} = S(\nu) \delta(\nu - \omega).
$$

The factor $\sqrt{2/\pi}$ arises because the conjugate symmetry (ii) makes the negative
frequencies duplicate the positive ones: $x(t) = (2/\sqrt{2\pi})\,\operatorname{Re}
\int_0^\infty e^{i\omega t} Z'(\omega)\, d\omega$, and $2/\sqrt{2\pi} = \sqrt{2/\pi}$.

Representation {eq}`eq-10-3` expresses $x(t)$ as a weighted sum of cosine and sine waves, with the weights being random processes $a(\omega)$ and $b(\omega)$. Equation {eq}`eq-10-4` represents $x(t)$ as a sum of cosine waves with amplitude $r(\omega)$ and phase $\theta(\omega)$ being governed by random processes.

Using representation {eq}`eq-10-3`, it is possible to show that

$$
Ex(t)^2 = \frac{1}{2\pi} \int^\infty_{-\infty} S(\omega)\, d\omega
= \frac{1}{\pi} \int^\infty_0 S(\omega)\, d\omega ,
$$

the last equality using $S(\omega) = S(-\omega)$. Formally
$E\,[a(\omega)^2 + b(\omega)^2] = E\,|Z'(\omega)|^2 = S(\omega)\,\delta(0)$, so the
"random amplitudes" $a(\omega),\ b(\omega)$ are, like $Z'$ itself, generalized processes:
the statement that has ordinary meaning is the one about *increments*, namely that the
variance contributed by the frequency band $[c,\, d]$ is $\frac{1}{\pi}\int_c^d S(\omega)\,
d\omega$, which is what the band-pass construction above computed.

```{eval-rst}
.. index::
   single: mean square ergodicity
   single: ergodicity; mean square
   single: ergodicity; criterion for
   single: ergodicity; spectral criterion
   single: long-run variance
   single: spectral distribution; atom at zero frequency
   single: random constant
   pair: ergodicity; sample mean
   see: ergodic; ergodicity
```

```{eval-rst}
.. index::
   single: time average; as a filter
   single: mean square ergodicity; spectral criterion
   single: random constant; and ergodicity
```

## The time average as a filter: mean square ergodicity

The machinery just assembled settles a question raised in
{doc}`01_covariance_stationary_processes` and left open there: when does averaging a *single
realization* over time reveal the *ensemble* mean? The answer costs nothing extra, because the
time average is simply one more filter of the kind {eq}`eq-10-2` describes.

Let $x(t)$ be covariance stationary with mean $\mu$, and form the time average over a record of
length $2T$,

$$
\bar x_T = \frac{1}{2T}\int^{T}_{-T} x(t)\, dt .
$$

Call $x(t)$ *mean square ergodic* if $E(\bar x_T - \mu)^2 \to 0$ as $T \to \infty$, so that
$\bar x_T$ converges in mean square to $\mu$.

Now observe that $\bar x_T$ is the output of a filter applied to $x$, with weighting function
$b(\tau) = 1/2T$ on $[-T,T]$ and zero elsewhere. Its transfer function is

$$
B_T(\omega) = \frac{1}{2T}\int^{T}_{-T} e^{-i\omega t}\, dt = \frac{\sin \omega T}{\omega T},
$$

so by {eq}`eq-10-1` the centered time average has the Cramér representation
$\bar x_T - \mu = (2\pi)^{-1/2}\int B_T(\omega) Z'(\omega)\, d\omega$, and by
{eq}`eq-10-2` its variance is

```{math}
:label: eq-10-erg
E(\bar x_T - \mu)^2 = \frac{1}{2\pi}\int^\infty_{-\infty}
\left(\frac{\sin \omega T}{\omega T}\right)^{2} S(\omega)\, d\omega .
```

Equation {eq}`eq-10-erg` answers the question by inspection. The kernel
$(\sin \omega T/\omega T)^2$ never exceeds $1$, equals $1$ at $\omega = 0$, and tends to $0$
at every $\omega \neq 0$ as $T$ grows. So as the record lengthens, the filter annihilates
every frequency except the origin: **time averaging is a band-pass filter whose band collapses
onto zero frequency.** What survives in the limit is whatever variance the process carries
*exactly at* $\omega = 0$.

For processes with an ordinary spectral density that limit is zero, and $x$ is mean square
ergodic. But the linearly deterministic component of Wold's theorem
({doc}`08_spectral_densities`) carries $\delta$-functions,
$S_d(\omega) = \sum_j a_j\pi[\delta(\omega - \lambda_j) + \delta(\omega + \lambda_j)]$, and
if one of the $\lambda_j$ is zero then {eq}`eq-10-erg` retains a term that never vanishes. Hence

> **Criterion.** A covariance stationary process is mean square ergodic **if and only if its
> spectral density carries no $\delta$-function at frequency zero**, equivalently iff its
> linearly deterministic component contains no zero-frequency term.

Failure has exactly one form: a *random constant*. If $x(t) \equiv A$ with $EA = 0$ and
$EA^2 = \sigma^2$, then $R(\tau) = \sigma^2$ for every $\tau$, the spectral density is
$2\pi\sigma^2\delta(\omega)$, and $\bar x_T = A$ for every $T$. The time average never
approaches the ensemble mean, because there is nothing to average. Every purely linearly
indeterministic process is therefore mean square ergodic, since its spectral density
$|\tilde P(i\omega)|^2$ is an ordinary function with no atoms. The book's standing assumption
already delivers the property.

Two consequences follow. First, when $\int |R(\tau)|\,d\tau < \infty$, letting
$T \to \infty$ in the time-domain counterpart of {eq}`eq-10-erg`,

$$
E(\bar x_T - \mu)^2 = \frac{1}{2T}\int^{2T}_{-2T}\Big(1 - \frac{|\tau|}{2T}\Big) R(\tau)\, d\tau,
$$

gives

```{math}
:label: eq-10-longrun
2T \cdot E(\bar x_T - \mu)^2 \;\longrightarrow\; \int^\infty_{-\infty} R(\tau)\, d\tau = S(0) .
```

The variance of the sample mean is asymptotically $S(0)/2T$: **the spectral density at the
origin, divided by the length of the record.** The "long-run variance" of time series
econometrics is nothing but the spectrum at zero frequency.

Second, and less comfortably: this settles the sample *mean* only. Whether the sample
*autocovariances* converge to $R(\tau)$ is a question about *fourth* moments. That is the
property on which the estimation in {doc}`21_phillips_continuous_time_estimation` and
{doc}`22_dimensionality_aliasing_problem` relies, which the second-moment theory of this book does not settle.
{doc}`/appendices/ergodicity` takes that up.

```{eval-rst}
.. index::
   single: band decomposition
   single: frequency band; variance carried by
   pair: sampling; frequency bands
```

## The band decomposition and sampling

This chapter exhibits a process as a superposition of mutually orthogonal frequency bands, each
carrying a definite share of the variance. Keep that picture in mind when reading
{doc}`17_discrete_sampling_folding`. Sampling at interval $T$ does not destroy the bands; it
*wraps* them. Every band is translated by a multiple of $2\pi/T$ onto the observable interval
$[-\pi/T,\, \pi/T]$ and there added to whatever was already present, which is precisely what the
folding formula asserts. Because the bands were orthogonal to begin with, their variances simply
add, and the discrete spectrum at a given frequency is the total variance of all the continuous
bands that fold onto it. Nothing tells the observer how that total was divided among them.

That last sentence is the aliasing problem, and {doc}`22_dimensionality_aliasing_problem` turns
it into a construction: given any continuous spectral density, one can build a *different* one,
supported entirely on high-frequency bands, that folds onto the same discrete spectrum. The
alternative process is manufactured with exactly the band-pass window $B_{cd}(\omega)$ of this
chapter, applied at frequencies above the Nyquist rate. Two processes that could hardly look
less alike in continuous time, one concentrated at low frequencies and the other carrying all of
its power above $\pi$, are then indistinguishable in the sampled data.

## Exercises

```{code-cell} ipython3
import numpy as np
from scipy.integrate import quad
```

```{exercise-start}
:label: cramer_ex1
```

**Variance by frequency band, and what sampling does to it.** Take the Ornstein–Uhlenbeck
process of {doc}`07_wiener_driven_sde` with $a = 1$, $b = 0.7$, whose spectral density is
$S(\omega) = b^2/(a^2+\omega^2)$.

(a) Verify that the band variances add up: partition $[0,\infty)$ into bands and check that
$\sum_{\text{bands}} \frac{1}{\pi}\int_c^d S(\omega)\, d\omega = R(0) = b^2/(2a)$, which is property (iii)
of {doc}`08_spectral_densities`.

(b) Report the share of the variance lying *above* the Nyquist frequency $\pi/T$ for
$T = 0.5$, $1$ and $3$. This is the power that sampling at interval $T$ folds back into the
observable band rather than discarding. That quantity governs how badly the formula of
{doc}`17_discrete_sampling_folding` distorts the spectrum.

(c) Confirm the orthogonality claim of the text numerically in the frequency domain: for
disjoint bands, $\int B_{cd}(\omega)\,B_{ef}(\omega)\,S(\omega)\,d\omega = 0$, so the band-pass filtered series are
uncorrelated at all leads and lags.

```{exercise-end}
```

```{solution-start} cramer_ex1
:class: dropdown
```

```{code-cell} ipython3
a, b = 1.0, 0.7
S = lambda w: b**2/(a**2 + w**2)
band_var = lambda c, d: quad(S, c, d, limit=200)[0]/np.pi     # variance in [c,d] and [-d,-c]

edges = [0, 0.5, 1, 2, 4, 8, 16, np.inf]
tot = 0.0
print("   band            variance     share")
for c, d in zip(edges[:-1], edges[1:]):
    v = band_var(c, d); tot += v
    print(f"  [{c:5.1f},{d:6.1f}]   {v:10.6f}   {v/(b**2/(2*a)):7.3%}")
print(f"\n  total = {tot:.8f},   R(0) = b^2/2a = {b**2/(2*a):.8f}")
```

```{code-cell} ipython3
# (b) share of variance above the Nyquist frequency pi/T
for T in [0.5, 1.0, 3.0]:
    nyq = np.pi/T
    print(f"T={T:4.1f}:  Nyquist = {nyq:5.3f},  share of variance above it = "
          f"{band_var(nyq, np.inf)/(b**2/(2*a)):6.2%}")
```

```{code-cell} ipython3
# (c) disjoint bands are orthogonal
B = lambda w, c, d: 1.0 if c <= w <= d else 0.0
cross = quad(lambda w: B(w,0.5,1.0)*B(w,2.0,4.0)*S(w), 0, 8, limit=400)[0]
print(f"cross term for disjoint bands [0.5,1] and [2,4]: {cross:.3e}")
```

The band variances sum to $R(0)$, and disjoint bands contribute nothing to each other. Part (b)
is the quantitative form of the warning in {doc}`17_discrete_sampling_folding`: at $T = 0.5$
about a tenth of the variance sits above the Nyquist frequency, so the folded spectrum departs
from the continuous one only modestly; by $T = 3$ nearly half of it does, and the discrete
spectrum is visibly inflated. The share falls slowly. The Lorentzian tail decays only as
$\omega^{-2}$, so halving the sampling interval does not come close to halving the aliased power. {doc}`22_dimensionality_aliasing_problem` pushes this to its extreme by
building a process whose variance lies *entirely* above the Nyquist frequency, yet which is
indistinguishable from this one in sampled data.

```{solution-end}
```

```{exercise-start}
:label: cramer_ex2
```

**Mean square ergodicity, and how it fails.** Continue with the Ornstein–Uhlenbeck process,
$R(\tau) = \frac{b^2}{2a}e^{-a|\tau|}$, $S(\omega) = b^2/(a^2+\omega^2)$, $a=1$, $b=0.7$.

(a) Evaluate $E(\bar x_T-\mu)^2$ from the time-domain formula
$\frac{1}{2T}\int_{-2T}^{2T}(1-|\tau|/2T)R(\tau)d\tau$ for $T = 5, 20, 100$, and confirm it
approaches $S(0)/2T$ of {eq}`eq-10-longrun`.

(b) Confirm the same numbers from the frequency-domain formula {eq}`eq-10-erg`, integrating
$(\sin\omega T/\omega T)^2 S(\omega)$. The two routes are the same computation seen from the
two sides of the Cramér representation.

(c) Now break ergodicity. Add a random constant: $y(t) = x(t) + A$ with $A$ independent of $x$,
mean zero, variance $\sigma^2 = 0.5$. Show that $E(\bar y_T - \mu)^2 \to \sigma^2$ rather than
to zero, however long the record. The spectral density of $y$ carries a $\delta$-function at
$\omega = 0$.

```{exercise-end}
```

```{solution-start} cramer_ex2
:class: dropdown
```

```{code-cell} ipython3
import numpy as np
from scipy.integrate import quad

a, b = 1.0, 0.7
R = lambda t: (b**2/(2*a))*np.exp(-a*abs(t))
S = lambda w: b**2/(a**2 + w**2)
S0 = b**2/a**2                                   # = S(0) = integral of R

print("  T      time domain      S(0)/2T       frequency domain")
for T in [5.0, 20.0, 100.0]:
    td = quad(lambda t: (1-abs(t)/(2*T))*R(t), -2*T, 2*T, limit=400)[0]/(2*T)
    fd = quad(lambda w: (np.sin(w*T)/(w*T))**2 * S(w), -np.inf, np.inf, limit=800)[0]/(2*np.pi)
    print(f"{T:6.1f}   {td:.8f}     {S0/(2*T):.8f}    {fd:.8f}")
```

```{code-cell} ipython3
# (c) a random constant added: variance of the sample mean no longer vanishes
sig2 = 0.5
print("  T      Var(mean of y) = Var(mean of x) + sigma^2")
for T in [5.0, 20.0, 100.0, 1000.0]:
    td = quad(lambda t: (1-abs(t)/(2*T))*R(t), -2*T, 2*T, limit=400)[0]/(2*T)
    print(f"{T:7.1f}   {td + sig2:.8f}   -> floor of sigma^2 = {sig2}")
```

The time-domain and frequency-domain evaluations agree to eight digits, and both approach
$S(0)/2T$. Adding the random constant puts an atom at the origin, and the variance of the
sample mean falls to $\sigma^2$ and stops: no record, however long, locates $\mu$. Note that
$A$ is invisible in any *centered* second moment. It shifts $R(\tau)$ by $\sigma^2$ at every
lag, which the sample autocovariances of a demeaned record cannot detect.

```{solution-end}
```
