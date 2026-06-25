# 12. Linear Least Squares Prediction

## (a) The Wiener–Kolmogorov Formula

The Wold moving average representation is useful for representing the linear
least squares $u$-step ahead prediction for a purely linearly indeterministic
process. Let $x(t)$ be a purely linearly indeterministic process with Wold
moving average representation

$$
x(t) = \int^\infty_0 p(s)\, w(t-s)\, ds
$$ (eq-12-wold)

where $\int^\infty_0 p(s)^2\, ds < +\infty$ and $w(t)$ is a fundamental white
noise for $x(t)$. Recall that the property that $w(t)$ is a fundamental white
noise for $x(t)$ means equality of the linear spaces

$$
\begin{aligned}
H_x(-\infty,\, t) \equiv\ & [\, y(t) : y(t) = \int^\infty_0 b(s)\, x(t-s)\, ds, \\
&\ \text{ for any }\ b(s) \in L_2\,[0,\, \infty) \,]
\end{aligned}
$$

and

$$
\begin{aligned}
H_w(-\infty,\, t) \equiv\ & [\, z(t) : z(t) = \int^\infty_0 h(s)\, w(t-s)\, ds \\
&\ \text{ for any }\ h(s) \in L_2\,[0,\, \infty) \,]
\end{aligned}
$$

where $L_2\,[0,\, \infty)$ is the space of square integrable functions, i.e.,
functions $b(s)$ for $0 \leq s < \infty$ such that $\int^\infty_0 b(s)^2\, ds <
+\infty$. The equality of these spaces means that lagged $x$'s contain the
same amount of information as lagged $w$'s.

Since {eq}`eq-12-wold` holds for all $t$, we have

$$
x(t + u) = \int^\infty_{s=-u} p(s+u)\, w(t-s)\, ds,\ \text{ for }\ u \geq 0.
$$

Using the identity of the linear spaces $H_x(-\infty,\,t)$ and $H_w(-\infty,\,t)$, we have that

$$
E\, [\, x(t+u) \mid x(v),\, v \leq t \,] = \int^\infty_{s=0} p(s+u)\, w(t-s)\, ds.
$$ (eq-12-wk)

Equation {eq}`eq-12-wk` is the continuous time Wiener–Kolmogorov formula. It rests on the
equality of the spaces spanned by past $x$'s and past fundamental innovations $w$'s established
in Wold's theorem (Theorem 10 of {doc}`08_spectral_densities`);
{doc}`13_kalman_filter_spectral_factorization` gives the equivalent state-space form, in which
$w$ becomes the Kalman innovations process and this same forecast is computed recursively.

Using operational calculus, the formula can be expressed as

$$
E_t\, x(t+u) = [\, \tilde P(D)\, e^{-Du} \,]_+\, w(t)
$$ (eq-12-wkop)

where $[\, \tilde P(s)\, e^{-su} \,]_+$ is the time function formed by taking
the inverse Laplace transform of $\tilde P(s)\, e^{-su}$, and then convoluting
it with the Heaviside unit step function (i.e., setting values of the time
function for $t < 0$ equal to zero, while leaving values of the function for $t
\geq 0$ unaltered). The operator $[\,\cdot\,]_+$ is known as the *annihilation
operator*. Note by property 4 (Delay) of Table 2 that $e^{-su}\, \tilde P(s)$ is
the Laplace transform of the function $p(s+u)$.

As an example of the use of formula {eq}`eq-12-wk`, let $x(t)$ be governed by the
first order stochastic differential equation

$$
(D+a)\, x(t) = w(t),\qquad a > 0
$$

so that $p(\tau) = e^{-a\tau}$. Then formula {eq}`eq-12-wk` gives

$$
\begin{aligned}
E_t\, x(t+u) &= \int^\infty_0 e^{-a(s+u)}\, w(t-s)\, ds \\
&= e^{-au} \int^\infty_0 e^{-as}\, w(t-s)\, ds = e^{-au}\ \frac{1}{D + a}\, w(t)
\end{aligned}
$$

or

$$
E_t\, x(t+u) = e^{-au}\, x(t).
$$ (eq-12-ar1pred)

## (b) A Formula for Predicting "Geometric Distributed Leads"

Such geometric distributed leads are the present values that appear in every asset-pricing
equation, permanent-income model, and quadratic-adjustment-cost Euler equation; they are also
the continuous-time counterpart of the discounted expected sums that define the optimal
feedforward decision rules of {doc}`14_faster_methods_recursive_linear_models`. Evaluating them
is the central computational step in solving continuous-time rational expectations models.

In linear rational expectation models, there often appear terms of the form

$$
E_t \int^\infty_0 e^{\rho u}\, x(t+u)\, du \qquad re(\rho) < 0
$$ (eq-12-glead)

where $x(t)$ is a covariance stationary stochastic process. Where $x(t)$ is
governed by the first order Markov process $(D+a)\, x(t) = w(t)$, equation {eq}`eq-12-ar1pred`
implies that the linear least squares forecast of the geometric distributed
lead {eq}`eq-12-glead` is given by

$$
\int^\infty_0 e^{\rho u}\, E_t\, x(t+u)\, du = \left( \int^\infty_0 e^{\rho u}\, e^{-au}\, du \right) x(t) = \frac{1}{a - \rho}\, x(t).
$$

An approach to the evaluation of {eq}`eq-12-glead` which readily generalizes to $x(t)$'s
governed by higher order linear differential equations is as follows. Denote
the geometric distributed lead to be forecast as

$$
\begin{aligned}
x(t)^{\ast} &= \int^\infty_0 e^{\rho u}\, x(t+u)\, du \\
x(t)^{\ast} &= \left( \frac{-1}{\rho + D} \right)\ \left( \frac{1}{a+D} \right)\, w(t)
\end{aligned}
$$ (eq-12-xstar)

where for $re(\rho) < 0,\ -1/\rho + s$ is the Laplace transform of the time
function $e^{-\rho u}$ for $u \leq 0,\ 0$ for $u > 0$.

Obtaining a partial fraction representation of the right side of {eq}`eq-12-xstar` gives

$$
x(t)^{\ast} = \frac{1}{a-\rho}\ \left[ \left( \frac{-1}{\rho + D} \right) + \left( \frac{1}{a+D} \right) \right]\ w(t)
$$

or

$$
x(t)^{\ast} = \left( \frac{1}{a-\rho} \right)\ \left[ -\int^\infty_0 e^{\rho s}\, w(t+s)\, ds + \int^\infty_0 e^{-as}\, w(t-s)\, ds \right]
$$

It then follows that

$$
\begin{aligned}
E_t\, x(t)^{\ast} &= \left( \frac{1}{a-\rho} \right) \int^\infty_0 e^{-as}\, w(t-s)\, ds = \left( \frac{1}{a-\rho} \right)\ \frac{1}{a + D}\ w(t) \\
&= \left( \frac{1}{a - \rho} \right) x(t).
\end{aligned}
$$ (eq-12-Estar)

This approach generalizes readily as follows. Represent equation {eq}`eq-12-Estar` as

$$
\begin{aligned}
E_t\, x(t)^{\ast} &= \left( \frac{1}{a - \rho} \right)\ \left( \frac{1}{a+D} \right)\, w(t) \\
E_t\, x(t)^{\ast} &= \left[ \frac{-\tilde P(D) + \tilde P(-\rho)}{D+\rho} \right]\, w(t)
\end{aligned}
$$ (eq-12-Estarop)

where $\tilde P(D) = 1/(a+D)$. As it happens, Equation {eq}`eq-12-Estarop` holds for *any*
$\tilde P(D)$, where $\tilde P(s)$ is the Laplace transform of a
squared summable function $p(\tau)$ concentrated on $\tau \in [0,\, \infty)$.
Thus, where

$$
x(t) = \int^\infty_0 p(\tau)\, w(t-\tau)\, d\tau
$$

we claim that the generalization of {eq}`eq-12-Estarop` is

$$
E_t \int^\infty_0 e^{\rho s}\, x(t+s)\, ds = \left[ \frac{-\tilde P(D) + \tilde P(-\rho)}{D+\rho} \right]\, w(t).
$$ (eq-12-gen)

The general formula {eq}`eq-12-gen`, valid for any rational $\tilde P(D)$, is established in
{doc}`19_prediction_formulas_continuous_time`.
