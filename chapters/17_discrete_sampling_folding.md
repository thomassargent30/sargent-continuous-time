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

# 17. Discrete Sampling: The Folding Formula

We now confront the gap between the continuous time theory developed so far and the discretely
sampled data we actually observe. The first question is what point-in-time sampling does to a
process's spectrum — and the answer, the *folding formula*, is the engine of the aliasing and
identification problems that occupy the remainder of the book.

Let $x(t)$ be a continuous time, covariance stationary stochastic process with
autocovariogram $R(\tau)$. Let $x_j = x(jT)$ be a record of $x$ sampled at the discrete
points in time $t = 0,\ \pm\, T,\ \pm\, 2T, \ldots,$ where $T > 0$ is the sampling
interval. The autocovariogram of the discrete data can be represented as a generalized
function $R^d(\tau)$ where

$$
R^d(\tau) = \sum_{n=-\infty}^{\infty} R(nT)\, \delta(\tau - nT),
$$ (eq-17-Rd)

so that $R^d(\tau)$ is a train of delta functions with mass $R(nT)$ at
$\tau = 0,\ \pm\, T,\ \pm\, 2T, \ldots\,$. Equation {eq}`eq-17-Rd` can also be represented as

$$
R^d(\tau) = R(\tau)\ S_T(\tau)
$$

where

$$
S_T(\tau) = \sum_{n=-\infty}^{\infty} \delta(\tau - nT).
$$

From {eq}`eq-17-Rd`, we can express the spectral density of the discrete sampled $x_t$ as

```{math}
:label: eq-17-1
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
:label: eq-17-2
\sum_{n=-\infty}^{\infty} \delta(\tau - nT) \leftrightarrow w_0 \sum_{n=-\infty}^{\infty} \delta(w - nw_0),\ w_o = 2\pi/T.
```

Now notice that

```{math}
:label: eq-17-3
w_0 \sum_{n=-\infty}^{\infty} S(w - nw_0) = S(w) \ast w_0 \sum_{n=-\infty}^{\infty} \delta(w - nw_0).
```

By the multiplication property (property 7 of Table 2 of
{doc}`08_spectral_densities`), multiplication in the time domain corresponds to convolution in
the frequency domain divided by $2\pi$, so the Fourier transform of
$R^d(\tau) = R(\tau)\, S_T(\tau)$ is

$$
\frac{1}{2\pi}\, S(w) \ast w_0 \sum_{n=-\infty}^{\infty} \delta(w - n w_0)
= \frac{w_0}{2\pi} \sum_{n=-\infty}^{\infty} S(w - n w_0)
= \frac{1}{T} \sum_{n=-\infty}^{\infty} S(w - n w_0),
$$

since $w_0/2\pi = 1/T$, the second equality being {eq}`eq-17-3` read from right to left.

Because {eq}`eq-17-1` identifies this transform as $S^d(w)$, the spectral density of the
discrete process $x_j$ satisfies

$$
S^d(w) = \frac{1}{T} \sum_{n=-\infty}^{\infty} S\!\left(w - n\, \frac{2\pi}{T}\right)
$$ (eq-17-fold)

where $S^d(w)$ is the spectral density of the discrete data and $S(w)$ is the spectral
density of the continuous time data. Equation {eq}`eq-17-fold` is known as the *folding formula*.

The folding formula is best read through the Cramér representation of
{doc}`10_cramer_representation`. There the process was exhibited as a superposition of
mutually orthogonal frequency bands, the band $[c,\, d]$ contributing variance
$\frac{1}{\pi}\int_c^d S(w)\, dw$. Sampling does not destroy those bands: it *wraps* them,
translating each by a multiple of $w_0 = 2\pi/T$ onto the observable interval
$[-\pi/T,\, \pi/T]$, where — because the bands were orthogonal — their variances simply add.
Equation {eq}`eq-17-fold` is that addition. The frequency $\pi/T$, half the sampling rate, is
the *Nyquist frequency*; everything above it is folded down on top of something below it.

The folding formula is therefore the engine of the *aliasing problem*: because infinitely many
continuous time frequencies fold onto the same discrete frequency, the sampled spectrum cannot
by itself recover the continuous one — the sampled data record only the *sum* of the folded
bands, never the division of that sum among them. {doc}`21_phillips_continuous_time_estimation` shows this
is the same phenomenon as the multivalued $\lambda = \log\mu$ in Phillips's estimation problem;
{doc}`22_dimensionality_aliasing_problem` counts how many continuous time models survive the
folding; and {doc}`23_temporal_aggregation_streamlined` takes up the related distortions caused
by time-averaging rather than point sampling.

## Exercises

The folding formula says that sampling a continuous time process at interval $T$ produces a
discrete time process whose spectral density is a sum of copies of the continuous spectrum,
shifted by integer multiples of $2\pi/T$ and superimposed ("folded") onto the band
$[-\pi/T, \pi/T]$. The frequency $\pi/T$ is the **Nyquist frequency**. Power that the
continuous process carries above the Nyquist frequency is not lost; it is *aliased* down
into the observable band, distorting the discrete spectrum.

We illustrate this with the Ornstein–Uhlenbeck process of Chapters 7–8, whose continuous
time spectral density is

$$
S(w) = \frac{b^2}{a^2 + w^2}, \qquad a, b > 0.
$$

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt
```

```{exercise-start}
:label: fold_ex1
```

Take $a = 1$, $b = 0.7$.

(a) Implement the folding formula

$$
S^d(w) = \frac{1}{T} \sum_{n=-\infty}^{\infty} S\!\left(w - n\,\frac{2\pi}{T}\right)
$$

(truncating the sum at $|n| \le n_{\max}$). For a *fast* sampling rate ($T = 0.5$) and a
*slow* one ($T = 3$), plot the continuous spectrum $S(w)$ and the folded discrete spectrum
$S^d(w)$ on the observable band $[0, \pi/T]$. Comment on the aliasing.

(b) As an independent check, the discrete spectral density is also the discrete Fourier
transform of the sampled autocovariance, $S^d(w) = \sum_n R(nT)\,e^{-iwnT}$ with
$R(\tau) = \frac{b^2}{2a}e^{-a|\tau|}$ (equation {eq}`eq-17-1`). Verify that this matches
your folding-formula computation.

```{exercise-end}
```

```{solution-start} fold_ex1
:class: dropdown
```

```{code-cell} ipython3
a, b = 1.0, 0.7
S = lambda w: b**2 / (a**2 + w**2)

def folded_spectrum(w, T, nmax=80):
    """Discrete spectral density via the folding formula."""
    return sum(S(w - n * 2 * np.pi / T) for n in range(-nmax, nmax + 1)) / T

fig, axes = plt.subplots(1, 2, figsize=(13, 4))
for ax, T in zip(axes, [0.5, 3.0]):
    w_nyq = np.pi / T                      # Nyquist frequency
    w = np.linspace(0, w_nyq, 400)
    ax.plot(w, S(w), 'b-', lw=2, label='continuous $S(w)$')
    ax.plot(w, folded_spectrum(w, T), 'r--', lw=2, label='folded $S^d(w)$')
    ax.set_title(f'sampling interval $T={T}$  (Nyquist $\\pi/T={w_nyq:.2f}$)')
    ax.set_xlabel('$w$'); ax.set_ylabel('spectral density')
    ax.legend()
plt.tight_layout()
plt.show()
```

```{code-cell} ipython3
# (b) cross-check: folding formula vs DFT of the sampled autocovariance
def Sd_from_acov(w, T, nmax=4000):
    n = np.arange(-nmax, nmax + 1)
    Rn = (b**2 / (2 * a)) * np.exp(-a * np.abs(n * T))   # R(nT)
    return np.array([(Rn * np.exp(-1j * ww * n * T)).sum().real for ww in w])

T = 1.0
w = np.linspace(0, np.pi / T, 6)
print("folding formula :", np.round(folded_spectrum(w, T), 4))
print("DFT of R(nT)    :", np.round(Sd_from_acov(w, T), 4))
```

With the fast sampling rate ($T = 0.5$, Nyquist $\approx 6.3$) the OU spectrum has already
decayed to nearly zero by the Nyquist frequency, so almost no power is folded back and
$S^d \approx S$ on the band. With the slow rate ($T = 3$, Nyquist $\approx 1.05$) a
substantial part of the continuous spectrum lives above the Nyquist frequency; it is
aliased back into the band and the folded spectrum sits visibly *above* the continuous one.
Part (b) confirms that the folding formula and the DFT of the sampled autocovariance agree
to numerical precision — two routes to the same discrete spectral density.

```{solution-end}
```
