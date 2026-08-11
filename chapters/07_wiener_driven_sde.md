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

# 7. Stochastic Differential Equations Driven by a Wiener Process

We wish to consider a class of stochastic differential equations driven by a Wiener process,

$$
dx(t) = f\left(x(t),\, t\right) dt + g\left(x(t),\,t\right) dW(t).
$$ (eq-7-sde)

We shall derive rules for manipulating such equations by regarding {eq}`eq-7-sde` as the limit as $\lambda \to \infty$ of

$$
dx(t) = f\left(x(t),\, t\right) dt + \frac{1}{\sqrt\lambda}\ g\left(x(t),\, t\right)\ \left(dN_1(t) - dN_2(t)\right),
$$ (eq-7-approx)

where $N_1$ and $N_2$ are two independent Poisson counters with identical rates $\lambda/2$.

The first rule we seek is the counterpart of rule {eq}`eq-5-rule1`. Where $\Psi(x)$ is a function of $x$, $\Psi(x)$ obeys the differential equation

$$
\begin{aligned}
d\Psi(x) &= \left\langle \frac{\partial\Psi}{\partial x}\, , \ f\left(x(t),\, t\right)\right\rangle dt + 1/2\, \left\langle \frac{\partial^2\Psi}{\partial x^2}\, g,\, g(x,\, t)\right\rangle dt \\
&+ \left\langle \frac{\partial\Psi}{\partial x}\, ,\ g(x,\, t)\right\rangle dW(t).
\end{aligned}
$$

```{note}
**The angle-bracket notation.** Throughout this chapter, $\langle u,\, v\rangle$ denotes the
**Euclidean inner product** (dot product) of two vectors $u$ and $v$ in $\mathbb{R}^n$,

$$
\langle u,\, v\rangle = u^\top v = \sum_{i=1}^{n} u_i\, v_i .
$$

Here $x$ is in general a vector-valued state, $\Psi(x)$ is a scalar-valued function,
$\dfrac{\partial\Psi}{\partial x}$ is its **gradient** (a column vector with entries
$\partial\Psi/\partial x_i$), and $\dfrac{\partial^2\Psi}{\partial x^2}$ is its **Hessian**
(the matrix of second partials $\partial^2\Psi/\partial x_i\,\partial x_j$). With this
reading,

- $\left\langle \dfrac{\partial\Psi}{\partial x}\, ,\ f\right\rangle
  = \sum_i \dfrac{\partial\Psi}{\partial x_i}\, f_i$ pairs the gradient with the drift
  vector $f$; and
- $\left\langle \dfrac{\partial^2\Psi}{\partial x^2}\, g,\ g\right\rangle
  = g^\top \dfrac{\partial^2\Psi}{\partial x^2}\, g
  = \sum_{i,\,j} \dfrac{\partial^2\Psi}{\partial x_i\,\partial x_j}\, g_i\, g_j$
  is the **quadratic form** in the diffusion vector $g$ (first the Hessian acts on $g$,
  then we take the inner product with $g$ again).

When $x$ is scalar, every inner product collapses to ordinary multiplication —
$\langle \Psi', f\rangle = \Psi' f$ and $\langle \Psi'' g,\, g\rangle = \Psi'' g^2$ — which
is the form used in the scalar example at the end of this chapter.
```

To generate this rule, we use {eq}`eq-7-approx` for $\lambda > 0$ and using rule {eq}`eq-5-rule1` to obtain

```{math}
:label: eq-7-1
\begin{aligned}
d\Psi(x) &= \left\langle \frac{\partial\Psi}{\partial x}\, ,\ f(x,\, t)\right\rangle dt \\
&+ \left[ \Psi\left(x + \frac{1}{\sqrt\lambda}\ g(x,\, t)\right) - \Psi(x) \right] dN_1(t) \\
&+ \left[ \Psi\left(x - \frac{1}{\sqrt\lambda}\ g(x,\, t)\right) - \Psi(x) \right] dN_2(t).
\end{aligned}
```

Obtain for $\Psi(x \pm\, \frac{1}{\sqrt\lambda}\ g(x,\, t))$ its Taylor series expansion about $x$,

$$
\begin{aligned}
\Psi\left(x \pm \frac{1}{\sqrt\lambda}\ g\left(x,\, t\right)\right) &= \Psi(x) \pm\ \left\langle \frac{\partial\Psi}{\partial x}\, ,\ \frac{1}{\sqrt\lambda}\ g(x,\, t)\right\rangle \\
&+ \frac{1}{2}\ \left(\frac{1}{\sqrt\lambda}\right)^2\ \left\langle \frac{\partial^2\Psi}{\partial x^2}\ g,\ g(x,\, t)\right\rangle + 0\, \left(\frac{1}{\lambda^{3/2}}\right)
\end{aligned}
$$

where $0\, (\frac{1}{\lambda^{3/2}}) \to 0$ as $\lambda \to \infty$. Substituting the Taylor expansions in {eq}`eq-7-1`, we obtain

```{math}
:label: eq-7-2
\begin{aligned}
d\Psi(x) &= \left\langle \frac{\partial\Psi}{\partial x}\, ,\ f(x,\, t)\right\rangle dt + \left\langle \frac{\partial\Psi}{\partial x}\, ,\ g(x,\, t)\right\rangle\ \frac{1}{\sqrt\lambda}\ (dN_1 - dN_2) \\
&+ \frac{1}{2}\ \left\langle \frac{\partial^2\Psi}{\partial x^2}\, g,\, g(x,\, t)\right\rangle\ \frac{1}{\lambda}\ (dN_1 + dN_2) + 0\, \left(\frac{1}{\lambda^{3/2}}\right)
\end{aligned}
```

Now consider the process $z(t)$ governed by

$$
dz(t) = \frac{1}{\lambda}\ (dN_1 + dN_2)
$$

Applying our rules, we find that

$$
\frac{d}{dt}\ Ez(t) = 1
$$

so that

$$
Ez(t) = t.
$$

We also find that

$$
\frac{d}{dt}\ Ez(t)^2 = 2Ez + 1/\lambda,
$$

so that as $\lambda \to \infty$

$$
\frac{d}{dt}\ Ez(t)^2 = 2Ez(t) = 2t.
$$

Therefore

$$
Ez(t)^2 = t^2.
$$

It follows that

$$
E\left(z(t) - Ez(t)\right)^2 = 0,
$$

which implies that as $\lambda \to \infty$,

$$
z(t) = t
$$

with probability 1.

Thus returning to {eq}`eq-7-2` and taking limits as $\lambda \to \infty$, we obtain

$$
\begin{aligned}
d\Psi(x) &= \left\langle \frac{\partial\Psi}{\partial x}\, ,\ f(x,\, t)\right\rangle dt + 1/2\, \left\langle \frac{\partial^2\Psi}{\partial x^2}\, g,\, g(x,\, t)\right\rangle dt \\
&+ \left\langle \frac{\partial\Psi}{\partial x}\, ,\ g(x,\, t)\right\rangle dW(t).
\end{aligned}
$$

This is known as *Ito's rule* for the stochastic differential equation

$$
dx = f(x,\, t)dt + g(x,\, t)dW(t).
$$

A word about where Itô's rule is used later, since the answer is more selective than one might
expect. The *linear* theory that occupies most of this book — the constant-coefficient models
of {doc}`11_linear_sde`, the prediction formulas of {doc}`12_prediction`, and everything built
on them — integrates white noise against a **deterministic** kernel $p(\tau)$. For such an
integrand the second-order term above vanishes, and, as the note in {doc}`11_linear_sde`
records, the Itô integral, the Stratonovich integral, and the ordinary mean square integral all
coincide. The Itô correction is invisible there precisely because nothing is nonlinear.

The correction does real work in one place: the derivation of the Riccati equation of the
Kalman–Bucy filter in {doc}`15_kalman_filter_spectral_factorization`. There one applies the rule
to the *quadratic* function $\Psi(e) = e\,e^\top$ of the estimation error, and the second-order
term supplies exactly the noise-intensity terms $BB^\top + KRK^\top$ that the Riccati equation
balances against the quadratic term $\Sigma C^\top R^{-1} C\,\Sigma$. Beyond that, the value of
Chapters 5 through 7 for what follows is that they *construct* the white noise the linear theory
takes as given, and show that the Poisson and Wiener cases are two faces of one limiting
argument.

Our second rule is

$$
\frac{d}{dt}\ E\Psi(x) = E\ \left\langle \frac{\partial\Psi}{\partial x}\, ,\ f(x,\, t)\right\rangle + \frac{1}{2}\, E\ \left\langle \frac{\partial^2\Psi}{\partial x^2}\, g,\, g\right\rangle.
$$ (eq-7-rule2)

This rule can be derived by the same limiting process. For $\lambda > 0$, we have

$$
\begin{aligned}
d\Psi(x) &= \left\langle \frac{\partial\Psi}{\partial x}\, , \ f(x,\, t)\right\rangle dt + \left( \Psi(x + \frac{1}{\sqrt\lambda}\ g(x,\, t)\right) \\
&- \Psi(x)\bigr) dN_1 + \left(\Psi(x - \frac{1}{\sqrt\lambda}\ g\left(x,\, t)\right) - \Psi(x)\right)dN_2
\end{aligned}
$$

which implies that

$$
\begin{aligned}
\frac{d}{dt}\ E\Psi(x) &= E\, \left\langle \frac{\partial\Psi}{\partial x}\, ,\ f(x,\, t)\right\rangle + E\, \left[\Psi(x + \frac{1}{\sqrt\lambda}\ g(x,\, t)) \right. \\
&+ \left. \Psi\left( x - \frac{1}{\sqrt\lambda}\ g(x,\, t) \right) - 2\ \Psi(x) \right]\frac{\lambda}{2},
\end{aligned}
$$

the factor $\lambda/2$ being the common rate of each of the two counters.

Now take Taylor series expansions of $\Psi\, (x \pm\, \frac{1}{\sqrt\lambda}\ g(x,\, t))$ about $x$ to get

$$
\begin{aligned}
\frac{d}{dt}\ E\Psi(x) &= E\, \left\langle \frac{\partial\Psi}{\partial x}\, ,\ f(x,\, t)\right\rangle + E\, \left\langle \frac{\partial\Psi}{\partial x}\, ,\ \frac{1}{\sqrt\lambda}\ \left(g(x,\, t) - g(x,\, t)\right)\right\rangle \\
&+ E\, \left\langle \frac{\partial^2\Psi}{\partial x^2}\, g,\ g\right\rangle\ \frac{1}{\lambda} \cdot \frac{\lambda}{2} + 0(\lambda^{-1/2}).
\end{aligned}
$$

The first-order terms cancel between the two counters; the two second-order terms are
identical, so together they contribute $(1/\lambda)\langle \Psi_{xx} g,\, g\rangle$, which the
rate $\lambda/2$ converts into $\tfrac12 \langle \Psi_{xx} g,\, g\rangle$.

Taking the limit as $\lambda \to \infty$ gives the desired result {eq}`eq-7-rule2`.

The next result that we desire is for $\tau > 0$,

$$
\frac{d}{d\tau}\ Ex(t) x(t + \tau) = Ex\, (t) f\left(x(t+\tau),\ t+\tau\right).
$$

To obtain this, we take limits as $\lambda \to \infty$ in the formula

$$
\frac{d}{d\tau}\ Ex(t) x(t+\tau) = Ex\, (t)f\left(x(t+\tau),\ t+\tau\right) + Ex\, (t)\ \left[\frac{1}{\sqrt\lambda}\ \left( \frac{\lambda}{2}\ - \ \frac{\lambda}{2}\right)\right].
$$

As an example of the use of these formulas, we take the linear stochastic differential
equation — the *Ornstein–Uhlenbeck process*, the simplest member of the constant-coefficient
class studied systematically in {doc}`11_linear_sde`,

$$
dx(t) = -ax(t)dt + bdW(t),\ x(0) = 0,\ a,\ b > 0.
$$

Applying our formulas, we find that

$$
\begin{aligned}
dx(t)^2 &= \left(-2ax(t)^2 + b^2\right) dt + 2bx(t)\, dW(t). \\
\frac{d}{dt}\ Ex(t)^2 &= -2a Ex(t)^2 + b^2 \\
\frac{d}{dt}\ Ex(t) &= -a Ex(t) \\
\frac{d}{d\tau}\ Ex(t) x(t + \tau) &= -aEx(t) x(t + \tau),\ \tau > 0.
\end{aligned}
$$

## Exercises

The linear stochastic differential equation

$$
dx(t) = -a\, x(t)\, dt + b\, dW(t), \qquad a, b > 0,
$$

is the **Ornstein–Uhlenbeck process**. The exercises below simulate it and check the
moment formulas just derived. We simulate using the **Euler–Maruyama** scheme: on a grid
of spacing $dt$,

$$
x_{k+1} = x_k - a\, x_k\, dt + b\,\sqrt{dt}\; \varepsilon_k, \qquad \varepsilon_k \sim \mathcal{N}(0, 1)\ \text{i.i.d.},
$$

which simply discretizes $dx = -a x\, dt + b\, dW$ with $dW = \sqrt{dt}\,\varepsilon$.

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt

def ou_paths(a, b, T, dt, n, rng, x0=0.0):
    """Simulate n Ornstein-Uhlenbeck paths on [0, T] via Euler-Maruyama."""
    steps = int(T / dt)
    X = np.empty((n, steps + 1))
    X[:, 0] = x0
    s = b * np.sqrt(dt)
    for k in range(steps):
        X[:, k + 1] = X[:, k] - a * X[:, k] * dt + s * rng.normal(size=n)
    return X
```

```{exercise-start}
:label: ou_ex1
```

Take $a = 1$, $b = 0.7$, and start all paths at $x(0) = 0$.

The text's equations integrate (with $x(0)=0$) to the transient variance

$$
E\,x(t)^2 = \frac{b^2}{2a}\bigl(1 - e^{-2at}\bigr),
$$

which approaches the **stationary variance** $b^2/(2a)$ as $t \to \infty$.

(a) Simulate an ensemble of paths and plot a handful of them together with the
theoretical $\pm$ one-standard-deviation band $\pm\sqrt{E\,x(t)^2}$.

(b) Check that the variance across paths at the final time matches the stationary value
$b^2/(2a)$.

```{exercise-end}
```

```{solution-start} ou_ex1
:class: dropdown
```

```{code-cell} ipython3
rng = np.random.default_rng(3)
a, b = 1.0, 0.7
dt, T, n = 0.02, 20.0, 8000

X = ou_paths(a, b, T, dt, n, rng)
t = np.arange(X.shape[1]) * dt

# (a) sample paths and theoretical std band
sd = np.sqrt((b**2 / (2 * a)) * (1 - np.exp(-2 * a * t)))
fig, ax = plt.subplots(figsize=(10, 4))
for i in range(5):
    ax.plot(t, X[i], lw=0.8)
ax.plot(t, sd, 'k--', lw=2, label=r'theoretical $\pm\sqrt{E\,x(t)^2}$')
ax.plot(t, -sd, 'k--', lw=2)
ax.set_xlabel('$t$'); ax.set_ylabel('$x(t)$')
ax.set_title('Ornstein–Uhlenbeck paths and the growing variance band')
ax.legend()
plt.show()

# (b) stationary variance
print(f"variance at t={T}: simulated {X[:, -1].var():.4f},  theory b^2/2a = {b**2 / (2 * a):.4f}")
```

The paths fan out from $0$ and settle into a stationary band whose width is the
mean-reverting balance between the diffusion $b\,dW$ pushing the process out and the drift
$-a x\,dt$ pulling it back.

```{solution-end}
```

```{exercise-start}
:label: ou_ex2
```

In the stationary regime the text shows the autocovariance decays exponentially,

$$
R(\tau) = R(0)\, e^{-a|\tau|}.
$$

Continuing the simulation from {ref}`ou_ex1`, discard an initial burn-in so the process is
(approximately) stationary, then estimate the normalized autocovariance
$R(\tau)/R(0)$ by averaging $x(t)\,x(t+\tau)$ over both time and paths. Compare it with
$e^{-a\tau}$.

```{exercise-end}
```

```{solution-start} ou_ex2
:class: dropdown
```

```{code-cell} ipython3
burn = int(10.0 / dt)            # drop the first 10 time units
Xs = X[:, burn:]
Xs = Xs - Xs.mean()
R0 = np.mean(Xs * Xs)

lags = np.arange(0, int(4.0 / dt))
acf = np.array([np.mean(Xs[:, :Xs.shape[1] - k] * Xs[:, k:]) for k in lags]) / R0
taus = lags * dt

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(taus, acf, 'o', ms=3, label=r'simulated $R(\tau)/R(0)$')
ax.plot(taus, np.exp(-a * taus), 'r-', lw=2, label=r'$e^{-a\tau}$')
ax.set_xlabel(r'$\tau$'); ax.set_ylabel(r'$R(\tau)/R(0)$')
ax.legend()
plt.show()
```

The estimated autocovariance tracks $e^{-a\tau}$. Note that $R(\tau) = R(0)e^{-a|\tau|}$ has a
kink at $\tau = 0$ — exactly like the $e^{-2\lambda|\tau|}$ of the telegraph wave of Chapter 5,
and for the same reason. The kink means that $R'(0^+) = -aR(0)$ while $R'(0^-) = +aR(0)$, so
$R''(0)$ does not exist; by Theorem 4 of
{doc}`02_mean_square_continuity_differentiability` the Ornstein–Uhlenbeck process is therefore,
like the telegraph wave, mean square continuous but **not** mean square differentiable. In the language of
{doc}`11_linear_sde` this is the case $n = 1$, $m = 0$, giving $n - 1 - m = 0$ derivatives; by
{doc}`13_locally_unpredictable` such a process is locally unpredictable. The two processes
differ in their sample paths — the diffusion is continuous, the telegraph wave jumps — but
their second moments are alike in precisely this respect.

```{solution-end}
```
