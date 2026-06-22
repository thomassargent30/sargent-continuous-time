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
$R(\tau) = \frac{b^2}{2a}e^{-a|\tau|}$ (equation {eq}`eq-15-1`). Verify that this matches
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
