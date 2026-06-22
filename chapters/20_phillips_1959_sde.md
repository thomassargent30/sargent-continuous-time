# 20. The Estimation of Parameters in Systems of Stochastic Differential Equations

*by A. W. Phillips (London School of Economics)*

```{note}
This chapter reproduces A. W. Phillips, "The Estimation of Parameters in Systems of
Stochastic Differential Equations," *Biometrika*, **46** (1959), no. 1–2, pp. 67–76. The
original section, equation, and footnote numbering of the paper is retained. The companion
essay {doc}`19_appreciation_phillips` discusses this paper and the literature that grew out
of it.
```

## 1. Summary

In many fields of study there occur systems of lagged dependences which can be described
approximately by systems of linear stationary stochastic differential equations. In this
paper we consider the estimation of parameters in such systems, using an approach similar to
that of Bartlett (1946, 1955) and Quenouille (1957), and discuss conditions for
identification of the separate behaviour relationships of a system. Little is known about the
small-sample properties of estimates in systems of this kind, but sampling experiments are
being carried out and the results will be published later.

## 2. Systems of linear stochastic lagged dependences

A system of linear stochastic lagged dependences between $n$ variables $x_i(t)$ which are
functions of continuous time can be described by the system of convolution integrals

```{math}
:label: eq-20-1
x_i(t) = \sum_{j \neq i} \int_0^\infty w_{ij}(h)\, x_j(t-h)\, dh + \xi_i(t) \quad (i,j = 1, 2, \ldots, n),
```

in which each of the weighting functions $w_{ij}(h)$ specifies the magnitude of the
dependence of the value taken by $x_i$ at time $t$ on the values taken by $x_j$ at times
$t-h$, and the $\xi_i(t)$ are stochastic disturbance terms. The Laplace transform of the
weighting function of a lagged dependence is called the transfer function of the dependence.
If all the transfer functions of the system {eq}`eq-20-1` are rational functions, i.e. if

```{math}
:label: eq-20-2
\int_0^\infty e^{-sh} w_{ij}(h)\, dh = \frac{\beta_{ij}(s)}{\alpha_{ij}(s)},
```

where $\beta_{ij}(s)$ and $\alpha_{ij}(s)$ are polynomials of finite degree in $s$,
{eq}`eq-20-1` can be written as the system of differential equations

```{math}
:label: eq-20-3
x_i(t) = \sum_{j \neq i} \frac{\beta_{ij}(D)}{\alpha_{ij}(D)}\, x_j(t) + \xi_i(t) \quad (i,j = 1, 2, \ldots, n)
```

in which $D$ is the differential operator $d/dt$. For the $i$th equation in {eq}`eq-20-1` is
then a particular integral of the $i$th equation in {eq}`eq-20-3`, and if the system is stable
and the initial conditions are at $t = -\infty$ it is the general solution of the $i$th
equation in {eq}`eq-20-3`.

Let $Q(D)$ be an $n \times n$ matrix whose elements $Q_{ij}(D)$ are given by

```{math}
:label: eq-20-4
Q_{ij}(D) = \begin{cases} -\dfrac{\beta_{ij}(D)}{\alpha_{ij}(D)} & \text{if } j \neq i, \\[2ex] 1 & \text{if } j = i, \end{cases}
```

and let $x(t)$ and $\xi(t)$ be column vectors with elements $x_i(t)$ and $\xi_i(t)$. Then
{eq}`eq-20-3` may be written

```{math}
:label: eq-20-5
Q(D)\, x(t) = \xi(t).
```

In many fields of study a system like {eq}`eq-20-3` may be used to represent the behaviour of
a number of separate units which form an interdependent system. Thus the equations may
represent the behaviour of units of physical equipment, as in control engineering, or of
groups of persons or firms, as in economics, or of animals or organisms, as in ecological or
biological studies. In describing real behaviour systems of this sort it is often reasonable
to assume that the behaviour units cannot respond to change instantaneously. This is
equivalent to the assumption that the weighting functions $w_{ij}(h)$ in {eq}`eq-20-1` do not
contain delta-function impulse terms,[^fn20-1] which in turn implies that the numerator
polynomial $\beta_{ij}(s)$ of each transfer function is of lower degree than the denominator
polynomial $\alpha_{ij}(s)$. In this paper we shall consider only systems in which the
numerator of each transfer function is of lower degree than the denominator. It will be seen
later that this restriction plays an important part in the identification of the separate
behaviour equations of a system.

## 3. Specification of the stochastic disturbances $\xi(t)$

We shall consider only systems which are stable and in which the disturbances $\xi_i(t)$ are
stationary stochastic variables whose covariance functions possess rational two-sided Laplace
transforms. It will be convenient to assume that the disturbances $\xi_i(t)$ are generated
from another set of stationary stochastic variables $\zeta_k(t)$ by the transformation

```{math}
:label: eq-20-6
\xi_i(t) = \sum_{k=1}^n \int_0^\infty r_{ik}(h)\, \zeta_k(t-h)\, dh \quad (i = 1, 2, \ldots, n)
```

in which the functions $r_{ik}(h)$ have Laplace transforms $R_{ik}(s)$ of the form

```{math}
:label: eq-20-7
R_{ik}(s) = \frac{\eta_{ik}(s)}{\gamma_{ik}(s)},
```

where $\eta_{ik}(s)$ and $\gamma_{ik}(s)$ are polynomials in $s$, the degree of
$\eta_{ik}(s)$ being lower than that of $\gamma_{ik}(s)$. Equations {eq}`eq-20-6` can then be
written as the set of differential equations

```{math}
:label: eq-20-8
\xi(t) = R(D)\, \zeta(t),
```

where $R(D)$ is a matrix with elements $R_{ik}(D)$ and $\zeta(t)$ is a column vector with
elements $\zeta_k(t)$. In order that the variables $\xi_i(t)$ be stationary it is necessary
that all the roots of the equations $\gamma_{ik}(s) = 0$ have negative real parts.

We define the following covariance functions:

```{math}
:label: eq-20-9
\begin{aligned}
\vartheta(\tau) &= E\{\xi(t)\, \xi'(t-\tau)\}, \\
\psi(\tau) &= E\{\zeta(t)\, \zeta'(t-\tau)\}, \\
\lambda(\tau) &= E\{\zeta(t)\, \xi'(t-\tau)\}
\end{aligned}
```

and denote their Laplace transforms by $\theta(s)$, $\Psi(s)$ and $\Lambda(s)$. Thus

```{math}
:label: eq-20-10
\begin{aligned}
\theta(s) &= \int_{-\infty}^\infty e^{-s\tau} \vartheta(\tau)\, d\tau, \\
\Psi(s) &= \int_{-\infty}^\infty e^{-s\tau} \psi(\tau)\, d\tau, \\
\Lambda(s) &= \int_{-\infty}^\infty e^{-s\tau} \lambda(\tau)\, d\tau.
\end{aligned}
```

It will now be shown that $\theta(s)$ is given by

```{math}
:label: eq-20-11
\theta(s) = R(s)\, \Psi(s)\, R'(-s).
```

The functions $r_{ik}(h)$ in {eq}`eq-20-6` may be regarded as defined but zero for $h < 0$,
and the lower limit of integration extended to $-\infty$. In matrix notation {eq}`eq-20-6`
then becomes

```{math}
:label: eq-20-12
\xi(t) = \int_{-\infty}^\infty r(h)\, \zeta(t-h)\, dh,
```

$r(h)$ being the matrix with elements $r_{ik}(h)$. Post-multiplying {eq}`eq-20-12` by
$\zeta'(t-\tau)$ and taking expectations we obtain using {eq}`eq-20-9`,

```{math}
:label: eq-20-13
\vartheta(\tau) = \int_{-\infty}^\infty r(h)\, \lambda(\tau-h)\, dh.
```

Taking the Laplace transform of {eq}`eq-20-13` and remembering that the transform of the
convolution of two functions is the product of their transforms[^fn20-2] we have

```{math}
:label: eq-20-14
\theta(s) = R(s)\, \Lambda(s).
```

Also, postmultiplying {eq}`eq-20-12` by $\zeta'(t+\tau)$ and taking expectations we find,
using {eq}`eq-20-9`

```{math}
:label: eq-20-15
\lambda'(\tau) = \int_{-\infty}^\infty r(h)\, \psi'(\tau+h)\, dh.
```

On making the substitution $h_1 = -h$ this becomes

```{math}
:label: eq-20-16
\lambda'(\tau) = \int_{-\infty}^\infty r(-h_1)\, \psi'(\tau-h_1)\, dh_1
```

and it follows that

```{math}
:label: eq-20-17
\Lambda'(s) = R(-s)\, \Psi'(s).
```

Substituting the transpose of {eq}`eq-20-17` in {eq}`eq-20-14` we obtain {eq}`eq-20-11`.

We now postulate that $\Psi(s)$ be identically the unit matrix,

```{math}
:label: eq-20-18
\Psi(s) = I.
```

Substitution of {eq}`eq-20-18` in {eq}`eq-20-11` then gives

```{math}
:label: eq-20-19
\theta(s) = R(s)\, R'(-s).
```

Since the spectral density and cross-spectral density functions of a set of variables can be
obtained from the two-sided Laplace transforms of their covariance functions by replacing the
complex number $s$ by the pure imaginary $i\omega$ (the Laplace transform then reducing to a
Fourier transform), the spectral densities of the variables $\zeta(t)$ are constants for all
frequencies $\omega$ and their cross-spectral densities are identically zero. The matrix
$\psi(\tau)$ of their covariance functions is diagonal, the leading elements being
$\delta$-function impulses. Since this implies that the variances are infinite, the variables
$\zeta(t)$ are not physically realizable.[^fn20-3] However, {eq}`eq-20-6` remains a valid
mathematical model for the generation of the physically realizable variables $\xi_i(t)$, whose
spectral density and cross-spectral density functions, given by $\theta(i\omega)$, are proper
rational functions of $i\omega$ as a result of the restriction that $\eta_{ik}(s)$ in
{eq}`eq-20-7` be of lower degree than $\gamma_{ik}(s)$.

## 4. Covariance functions of the system variables $x(t)$

Substituting {eq}`eq-20-8` in {eq}`eq-20-5` we have

```{math}
:label: eq-20-20
Q(D)\, x(t) = R(D)\, \zeta(t).
```

Let $P(D)$ be a diagonal matrix whose $i$th leading element is the lowest common multiple of
the denominators of the elements in the $i$th row of $Q(D)$ and in the $i$th row of $R(D)$,
i.e. the lowest common multiple of the polynomials $\alpha_{ij}(D)$, $j \neq i$, in
{eq}`eq-20-4` and $\gamma_{ik}(D)$, $k = 1, 2, \ldots, n$ in {eq}`eq-20-7`, and let $F(D)$ and
$G(D)$ be matrices defined by

```{math}
:label: eq-20-21
F(D) = P(D)\, Q(D)
```

and

```{math}
:label: eq-20-22
G(D) = P(D)\, R(D).
```

Then premultiplying {eq}`eq-20-20` by $P(D)$ we obtain

```{math}
:label: eq-20-23
F(D)\, x(t) = G(D)\, \zeta(t).
```

The elements of $F(D)$ and $G(D)$ are polynomials in $D$, the diagonal element in any row of
$F(D)$ being of higher degree than the non-diagonal elements in the same row and also of
higher degree than the elements in the same row of $G(D)$. Thus the determinant $|F(D)|$
cannot vanish identically, and we may write {eq}`eq-20-23` as

```{math}
:label: eq-20-24
x(t) = F^{-1}(D)\, G(D)\, \zeta(t).
```

Let $\phi(\tau)$ be the matrix of covariance functions of the system variables $x(t)$, so that

```{math}
:label: eq-20-25
\phi(\tau) = E\{x(t)\, x'(t-\tau)\}
```

and let $\Phi(s)$ be the Laplace transform of $\phi(\tau)$. Then, by a similar argument to
that used in §3, we have

```{math}
:label: eq-20-26
\Phi(s) = F^{-1}(s)\, G(s)\, \Psi(s)\, G'(-s)\, F'^{-1}(-s)
```

or, since $\Psi(s) = I$,

```{math}
:label: eq-20-27
\Phi(s) = F^{-1}(s)\, G(s)\, G'(-s)\, F'^{-1}(-s).
```

If we define the polynomial $p(s)$ by

```{math}
:label: eq-20-28
p(s) = |F(s)| = \prod_{r=1}^m (s - \lambda_r)
```

and write $F^a(s)$ for the adjoint of $F(s)$, {eq}`eq-20-27` becomes

```{math}
:label: eq-20-29
\Phi(s) = \frac{F^a(s)\, G(s)\, G'(-s)\, F'^a(-s)}{p(s)\, p(-s)}.
```

For simplicity we shall consider in this paper only cases in which the roots $\lambda_r$, $r =
1, 2, \ldots, m$, of the equation $p(s) = 0$ are all distinct, and shall also assume that
$p(s)$ does not contain any factor which is common to all the elements of $F^a(s)$. Expanding
the right-hand side of {eq}`eq-20-29` in partial fractions we then obtain

```{math}
:label: eq-20-30
\Phi(s) = \sum_{r=1}^m \frac{K_r}{s - \lambda_r} + \sum_{r=1}^m \frac{K_r'}{-s - \lambda_r},
```

where

```{math}
:label: eq-20-31
K_r = \left. \frac{(s - \lambda_r)\, F^a(s)\, G(s)\, G'(-s)\, F'^a(-s)}{p(s)\, p(-s)} \right|_{s = \lambda_r}.
```

For a simple root $\lambda_r$ the matrix $F(\lambda_r)$ is necessarily simply degenerate. The
adjoint matrix $F^a(\lambda_r)$ is of unit rank and can be written

```{math}
:label: eq-20-32
F^a(\lambda_r) = k_r \kappa_r,
```

where $k_r$ is a column vector which satisfies the equation

```{math}
:label: eq-20-33
F(\lambda_r)\, k_r = 0
```

and $\kappa_r$ is a row vector which satisfies

```{math}
:label: eq-20-34
\kappa_r\, F(\lambda_r) = 0.
```

(See Frazer, Duncan & Collar (1950) §5.6.) Let $l_r$ be a row vector defined by

```{math}
:label: eq-20-35
l_r = \left. \frac{(s - \lambda_r)\, \kappa_r\, G(s)\, G'(-s)\, F'^a(-s)}{p(s)\, p(-s)} \right|_{s = \lambda_r}.
```

Then substituting {eq}`eq-20-32` in {eq}`eq-20-31` and using {eq}`eq-20-35` we have

```{math}
:label: eq-20-36
K_r = k_r\, l_r.
```

Thus the matrix $K_r$ is of unit rank, with columns proportional to the vector $k_r$. It
follows from {eq}`eq-20-36` and {eq}`eq-20-33` that

```{math}
:label: eq-20-37
F(\lambda_r)\, K_r = 0 \quad (r = 1, 2, \ldots, m),
```

a result which could also have been obtained by premultiplying {eq}`eq-20-31` by $F(s)$ and
putting $s = \lambda_r$.

The matrix $\phi(\tau)$ of the system covariance functions is obtained by taking the inverse
transform of {eq}`eq-20-30`. Formally, {eq}`eq-20-30` does not have a unique inverse. However,
for a stationary system the roots $\lambda_r$ are negative or have negative real parts and
$\phi(\tau) \to 0$ as $\tau \to \pm\infty$. It can readily be verified that the only inverse
which is consistent with these conditions is that obtained by taking the first sum in
{eq}`eq-20-30` to be the one-sided transform of $\phi(\tau)$ for $\tau \geq 0$ and the second
sum to be the one-sided transform for $\tau < 0$. The inverse of the first sum in
{eq}`eq-20-30` then gives

```{math}
:label: eq-20-38
\phi(\tau) = \sum_{r=1}^m K_r\, e^{\lambda_r \tau} \quad (\tau \geq 0)
```

while from the inverse of the second sum we see that $\phi(-\tau) = \phi'(\tau)$, a fact which
is also obvious from {eq}`eq-20-25`.

## 5. Calculation of $\Phi(s)$ from given points on $\phi(\tau)$

Let $\phi_\tau$ be the matrix of discrete functions consisting of points at unit intervals of
$\tau$ on the continuous functions $\phi(\tau)$, i.e. let

```{math}
:label: eq-20-39
\phi_\tau = \phi(\tau) \quad (\tau = \ldots, -2, -1, 0, 1, 2, \ldots),
```

and let $M(z)$ be the one-sided generating function defined by

```{math}
:label: eq-20-40
M(z) = \sum_{\tau=0}^\infty \phi_\tau\, z^\tau.
```

Then from {eq}`eq-20-38`, {eq}`eq-20-39` and {eq}`eq-20-40` we find that

```{math}
:label: eq-20-41
M(z) = \sum_{r=1}^m \frac{K_r}{1 - \mu_r z},
```

where

```{math}
:label: eq-20-42
\mu_r = e^{\lambda_r}.
```

We define the polynomial

```{math}
:label: eq-20-43
\Delta(z) = \prod_{r=1}^m (1 - \mu_r z),
```

and seek matrices $U(z)$ and $V(z)$, whose elements are polynomials in $z$, such that

```{math}
:label: eq-20-44
M(z) = U^{-1}(z)\, V(z).
```

If {eq}`eq-20-44` is to hold we require that

```{math}
:label: eq-20-45
|U(z)| = \Delta(z)
```

and that

```{math}
:label: eq-20-46
U(z) \sum_{r=1}^m \frac{K_r}{1 - \mu_r z} = V(z).
```

Multiplying {eq}`eq-20-46` by $1 - \mu_q z$ we obtain

```{math}
:label: eq-20-47
U(z) \left\{ K_q + (1 - \mu_q z) \sum_{r \neq q} \frac{K_r}{1 - \mu_r z} \right\} = (1 - \mu_q z)\, V(z)
```

and on putting $z = 1/\mu_q$, $q = 1, 2, \ldots, m$, this gives

```{math}
:label: eq-20-48
U\!\left( \frac{1}{\mu_q} \right) K_q = 0 \quad (q = 1, 2, \ldots, m).
```

Since $K_q$ is of unit rank, {eq}`eq-20-48` provides $m$ scalar equations in the coefficients
of the polynomials occurring in each row of $U(z)$. If $U(z)$ is given the form

```{math}
:label: eq-20-49
U(z) = I + U_1 z + U_2 z^2 + \ldots + U_c z^c,
```

where $U_1, U_2, \ldots, U_c$ are $n \times n$ matrices whose elements are constants, each row
of $U(z)$ will contain $cn$ coefficients of powers of $z$ (excluding the constant terms, which
are unity or zero). It can be shown that if $c$ is given a suitable value and $cn - m$
coefficients in each row of $U(z)$ are put equal to zero, equations {eq}`eq-20-48` suffice to
determine the remaining $m$ coefficients in each row of $U(z)$, provided the positions of the
zeros are chosen in such a way that $|U(z)|$ is of degree $m$. It will always be possible to
find positions for the zeros such that this condition is satisfied, and equation
{eq}`eq-20-45` will then also be satisfied. If $m$ is an integral multiple of $n$ it will
usually, but not always, be possible to take $c$ equal to $m/n$, in which case $cn - m = 0$
and none of the coefficients is put equal to zero.

When $U(z)$ has been determined $V(z)$ can be found from {eq}`eq-20-46`. It can be seen from
{eq}`eq-20-48` that $1 - \mu_r z$ is a factor of $U(z)\, K_r$ and it follows from this and
{eq}`eq-20-46` that the degree of $V(z)$ is at least one less than that of $U(z)$. Thus $V(z)$
is of the form

```{math}
:label: eq-20-50
V(z) = V_0 + V_1 z + \ldots + V_{c-1} z^{c-1},
```

where $V_0, V_1, \ldots, V_{c-1}$ are matrices whose elements are constants.

Premultiplying {eq}`eq-20-44` by $U(z)$ and writing the matrices in expanded form, we have

```{math}
:label: eq-20-51
(I + U_1 z + \ldots + U_c z^c)(\phi_0 + \phi_1 z + \ldots) = V_0 + V_1 z + \ldots + V_{c-1} z^{c-1}.
```

Equating the coefficients of $z^c$, $z^{c+1}, \ldots$, in {eq}`eq-20-51` we obtain the
difference equation

```{math}
:label: eq-20-52
\phi_\tau + U_1 \phi_{\tau-1} + \ldots + U_c \phi_{\tau-c} = 0 \quad (\tau \geq c),
```

from which the matrices $U_1, U_2, \ldots, U_c$ can be calculated if $\phi_\tau$ is given for a
sufficient number of values of $\tau$. Equating the constant terms in {eq}`eq-20-51` we have

```{math}
:label: eq-20-53
V_0 = \phi_0,
```

and equating the coefficients of $z, z^2, \ldots, z^{c-1}$ we find

```{math}
:label: eq-20-54
V_l = \phi_l + \sum_{k=1}^l U_k \phi_{l-k} \quad (l = 1, 2, \ldots, c-1),
```

from which the matrices $V_0, V_1, \ldots, V_{c-1}$ can be calculated.

Thus given a sufficient number of points at equal intervals on the system covariance functions
$\phi(\tau)$ it is possible to calculate $U(z)$ and $V(z)$. $M(z)$ is then given by
{eq}`eq-20-44`, and can be expanded in partial fractions to give {eq}`eq-20-41`, from which the
matrices $K_r$ and the scalars $\mu_r$ are obtained immediately. The roots $\lambda_r$ of the
equation $p(s) = 0$ can then be calculated from

```{math}
:label: eq-20-55
\lambda_r = \log \mu_r,
```

which follows from {eq}`eq-20-42`, and the matrix $\Phi(s)$ can then be found from the $K_r$
and $\lambda_r$ using {eq}`eq-20-30`.

## 6. Calculation of $Q(s)$ and $\theta(s)$ from $\Phi(s)$ when all equations are of the same order

Consider an $n$ variable system defined by {eq}`eq-20-20` or {eq}`eq-20-23`, in which all the
scalar equations of $F(D)\, x(t) = G(D)\, \zeta(t)$ are of the same order $u$, and which cannot
be reduced to an equivalent system in which any equation is of lower order. The elements on
the leading diagonal of $F(s)$ are polynomials of degree $u$; but as a result of the
restriction on the degree of the numerator terms of the transfer functions in the matrix
$Q(s)$, discussed in §2 above, the non-diagonal elements are of degree $u - 1$ at most. The
coefficient of $s^u$ in each of the diagonal elements of $F(s)$ can be taken as unity. There
are $u$ other coefficients of powers of $s$ (including the constant term) in each of the
diagonal elements and $u$ coefficients in each of the non-diagonal elements. In calculating
$F(s)$ from $\Phi(s)$, or equivalently from the $K_r$ and $\lambda_r$, $r = 1, 2, \ldots, m$,
there are therefore $nu$ coefficients to be determined in each row of $F(s)$. Since $|F(s)|$ is
of degree $nu$ we have $m = nu$. We obtain from {eq}`eq-20-37` a set of $nu$ non-homogeneous
equations in the coefficients of each row of $F(s)$. It can be shown that the equations in each
set are linearly independent. Thus all the coefficients in $F(s)$ are uniquely determined.

Given $F(s)$ we can find $Q(s)$ from

```{math}
:label: eq-20-56
Q(s) = P^{-1}(s)\, F(s),
```

which follows from {eq}`eq-20-21`. It will be noticed that the elements of the diagonal matrix
$P(s)$ are identical with the elements on the leading diagonal of $F(s)$. The numerator and
denominator polynomials of the elements of $Q(s)$ as calculated from {eq}`eq-20-56` may have
common factors which can be cancelled.

Given $F(s)$ and the $K_r$ and $\lambda_r$, $G(s)\, G'(-s)$ is readily obtained from the
equation

```{math}
:label: eq-20-57
G(s)\, G'(-s) = \sum_{r=1}^m \frac{F(s)\, K_r\, F'(-s)}{s - \lambda_r} + \sum_{r=1}^m \frac{F(s)\, K_r'\, F'(-s)}{-s - \lambda_r}
```

which is derived from {eq}`eq-20-27` and {eq}`eq-20-30`. Using {eq}`eq-20-19` and
{eq}`eq-20-22` we then obtain

```{math}
:label: eq-20-58
\theta(s) = P^{-1}(s)\, G(s)\, G'(-s)\, P'^{-1}(-s).
```

The numerator and denominator polynomials of the elements of $\theta(s)$ as calculated from
{eq}`eq-20-58` may also have common factors which can be cancelled.

## 7. Calculation of $Q(s)$ and $\theta(s)$ from $\Phi(s)$ when all the equations are not of the same order, but the disturbances are not cross-correlated

If all the equations in a system are not of the same order they cannot all be identified
unless further restrictions are placed on the system. Consider a system $F(D)\, x(t) = G(D)\,
\zeta(t)$ in which the equations are not all of the same order, and in which the arrangement of
the equations is such that the first equation is of the lowest order and the order of each
successive equation is at least as high as that of the one preceding it. If {eq}`eq-20-27` is
satisfied by the matrices $F(s)$ and $G(s)$ of this system it is also satisfied by $\mu F(s)$
and $\mu G(s)$ where $\mu$ is an arbitrary matrix of constants. We can add to any row of $F(s)$
a multiple of an earlier row corresponding to an equation of lower order without violating the
condition that the non-diagonal elements in each row be of lower degree than the diagonal
elements; but we cannot add a multiple of any succeeding row without violating this condition.
As a result of this condition, therefore, $\mu$ is restricted to a triangular matrix with
zeros above the leading diagonal. Given $\Phi(s)$ for this system, the first equation is
therefore determinate (and so are any other equations of the same order as the first one) but
the remaining equations of higher order than the first one are not identified without further
restrictions on the system.

If it is known that the disturbances to the system are not cross-correlated, $G(s)\, G'(-s)$
must be diagonal. If $\mu F(s)$ and $\mu G(s)$ are to be admissible solutions, therefore, $\mu$
must satisfy the further restriction that the non-diagonal elements of $\mu G(s)\, G'(-s)\,
\mu'$ be identically zero. Writing $g_{ii}$ for the diagonal elements of $G(s)\, G'(-s)$ the
first part of the array formed by the elements lying above the leading diagonal of the matrix
product $\mu G(s)\, G'(-s)\, \mu'$ is found to be

```{math}
:label: eq-20-array
\begin{array}{lll}
\mu_{11} g_{11} \mu_{21} & \mu_{11} g_{11} \mu_{31} & \mu_{11} g_{11} \mu_{41} \\[1ex]
& \mu_{21} g_{11} \mu_{31} + \mu_{22} g_{22} \mu_{32} & \mu_{21} g_{11} \mu_{41} + \mu_{22} g_{22} \mu_{42} \\[1ex]
& & \mu_{31} g_{11} \mu_{41} + \mu_{32} g_{22} \mu_{42} + \mu_{33} g_{33} \mu_{43}.
\end{array}
```

Since $\mu_{ii}$ cannot be made zero, we see from the first element in the array that if
$g_{11} \neq 0$, $\mu_{21}$ must vanish. From the second row of the array we see that if
$g_{11} \neq 0$ and $g_{22} \neq 0$, $\mu_{31}$ and $\mu_{32}$ must vanish. Similarly, from the
third row, $\mu_{41}$, $\mu_{42}$ and $\mu_{43}$ must vanish if $g_{11}$, $g_{22}$ and $g_{33}$
are non-zero. Continuing in this way it can be shown that all the elements $\mu_{ij}$, $(i >
j)$ must vanish if all the elements $g_{ii}$, $(i = 1, 2, \ldots, n-1)$ are non-zero. The
matrix $\mu$ is then diagonal, and all the equations in the system are identified.

To calculate $F(s)$ for the system we have just been considering, given $\Phi(s)$ we first note
that since the equations of lowest order are identified without restriction of the disturbance
covariance functions, it will be possible to find coefficients in at least one of the rows of
$F(s)$ which satisfy the appropriate set of scalar equations derived from {eq}`eq-20-37` for
all the values of $r$. The number of equations will be greater than the number of non-zero
coefficients in this row, but they will form a consistent set with a sufficient number of
independent equations to determine the non-zero coefficients. In the present case the first
equation in the system is of lowest order, so the coefficients in the first row of $F(s)$ can
be found in this way.

From {eq}`eq-20-27` we have

```{math}
:label: eq-20-59
F(s)\, \Phi(s)\, F'(-s) = G(s)\, G'(-s).
```

On the assumption that $G(s)\, G'(-s)$ is diagonal the first column of the matrix product in
{eq}`eq-20-59` yields the equations

```{math}
:label: eq-20-60
\sum_{j=1}^n \sum_{q=1}^n F_{ij}(s)\, \Phi_{jq}(s)\, F_{1q}(-s) = 0 \quad (i = 2, 3, \ldots, n).
```

If we give $\Phi_{jq}(s)$ and $F_{1q}(-s)$ in {eq}`eq-20-60` their known numerical values,
multiply through by $p(s)\, p(-s)$ and equate the coefficients of each power of $s$ we obtain
for each value of $i$ a consistent set of linear equations in the coefficients of the
polynomials in the $i$th row of $F(s)$, which can be solved for these coefficients. $Q(s)$ can
then be calculated from $F(s)$, and $\theta(s)$ is given by

```{math}
:label: eq-20-61
\theta(s) = Q(s)\, \Phi(s)\, Q'(-s).
```

## 8. Stochastic difference equation for discrete points on $x(t)$

Let $\Gamma(z)$ be the two-sided generating function of the discrete functions $\phi_\tau$
defined in {eq}`eq-20-39`, so that

```{math}
:label: eq-20-62
\Gamma(z) = \sum_{\tau=-\infty}^\infty \phi_\tau\, z^\tau.
```

Since $\phi_{-\tau} = \phi_\tau'$ we have from {eq}`eq-20-40`

```{math}
:label: eq-20-63
\Gamma(z) = M(z) + M'(z^{-1}) - \phi_0.
```

Substituting {eq}`eq-20-44` in {eq}`eq-20-63` we obtain

```{math}
:label: eq-20-64
\Gamma(z) = U^{-1}(z)\, V(z) + V'(z^{-1})\, U'^{-1}(z^{-1}) - \phi_0
```

or

```{math}
:label: eq-20-65
U(z)\, \Gamma(z)\, U'(z^{-1}) = W(z)\, W'(z^{-1}),
```

where $W(z)$ is any matrix which satisfies the equation

```{math}
:label: eq-20-66
W(z)\, W'(z^{-1}) = V(z)\, U'(z^{-1}) + U(z)\, V'(z^{-1}) - U(z)\, \phi_0\, U'(z^{-1}).
```

It follows from {eq}`eq-20-65` that $\Gamma(z)$ is also the generating function of the
covariance functions of discrete variables $x_t$ which are generated by the system of
difference equations

```{math}
:label: eq-20-67
U(E^{-1})\, x_t = W(E^{-1})\, \epsilon_t,
```

from stationary stochastic variables $\epsilon_t$ with covariance functions

```{math}
:label: eq-20-68
E\{\epsilon_t\, \epsilon_{t-\tau}'\} = \begin{cases} I & \text{if } \tau = 0, \\ 0 & \text{if } \tau \neq 0. \end{cases}
```

(The symbol $E^{-1}$ in {eq}`eq-20-67` is the shift operator, defined by $E^{-h} y_t =
y_{t-h}$.)

The discrete variables $x_t$ in {eq}`eq-20-67` may be identified with points at unit time
intervals on the continuous variables $x(t)$ in {eq}`eq-20-20`, and it follows that the
sampling properties of estimates of the continuous system {eq}`eq-20-20` may be studied by
considering the sampling properties of estimates of the discrete system {eq}`eq-20-67`. The
work done by Bartlett (1946, 1955), M. G. Kendall (1945, 1946, 1949), Quenouille (1947, 1957),
Whittle (1952, 1953), Wold (1949) and others on sampling properties and goodness-of-fit tests
for autoregressive, moving average and mixed processes is therefore relevant to the estimation
of parameters in systems of stochastic differential equations. Most of the theoretical results
which have been obtained, however, are valid only for large samples.

Some sampling experiments are being carried out on systems of the type considered in this
paper. The results will be published later.

## References

Bartlett, M. S. (1946). On the theoretical specification and sampling properties of
autocorrelated time series. *J. Roy. Statist. Soc. B*, **8**, 27.

Bartlett, M. S. (1955). *Stochastic Processes*. Cambridge University Press.

Frazer, R. A., Duncan, W. J. & Collar, A. R. (1950). *Elementary Matrices*. Cambridge
University Press.

James, H. M., Nichols, N. B. & Phillips, R. S. (1947). *Theory of Servomechanisms*. New York:
McGraw Hill.

Kendall, M. G. (1945). On the analysis of oscillatory time series. *J. Roy. Statist. Soc.*,
**108**, 93.

Kendall, M. G. (1946). *Contributions to the Study of Oscillatory Time Series*. Cambridge
University Press.

Kendall, M. G. (1949). Tables of autoregressive series. *Biometrika*, **36**, 267.

Laning, J. H. & Battin, R. H. (1956). *Random Processes in Automatic Control*. New York:
McGraw Hill.

Miller, K. S. (1956). *Engineering Mathematics*. London: Constable and Co. Ltd.

Quenouille, M. H. (1947). A large-sample test for the goodness-of-fit of autoregressive
schemes. *J. Roy. Statist. Soc. A*, **110**, 123.

Quenouille, M. H. (1957). *The Analysis of Multiple Time Series*. London: Charles Griffin and
Co. Ltd.

Whittle, P. (1952). Tests of fit in time series. *Biometrika*, **39**, 309.

Whittle, P. (1953). The analysis of multiple stationary time-series. *J. Roy. Statist. Soc.
B*, **15**, 125.

Widder, D. V. (1946). *The Laplace Transform*. Princeton University Press.

Wold, H. (1949). A large-sample test for moving averages. *J. Roy. Statist. Soc.*, **11**,
297.

## Notes

[^fn20-1]: The $\delta$-function or unit impulse function, $\delta(h)$ is defined by the
    equations
    $$\delta(h) = \begin{cases} 0 & \text{if } h \neq 0 \\ \infty & \text{if } h = 0 \end{cases} \quad \text{and} \quad \int_{-\infty}^\infty \delta(h)\, dh = 1.$$
    (See James, Nichols & Phillips (1947), chapter 2, especially §2.4.) The Laplace transform
    of $\delta(h)$ is clearly unity. If $\beta_{ij}(s)$ were not of lower degree than
    $\alpha_{ij}(s)$ the partial fraction expansion of $\beta_{ij}(s)/\alpha_{ij}(s)$ would
    contain a constant term and so its inverse transform $w_{ij}(h)$ would contain a
    $\delta$-function impulse term.

[^fn20-2]: See Widder (1946), p. 258, theorem 16b. The conditions for the theorem are
    satisfied as a result of the stationarity of $\zeta(t)$ and the conditions that the roots
    of the equations $\gamma_{ik}(s) = 0$ have negative real parts.

[^fn20-3]: In engineering literature the variables $\zeta(t)$ are called independent white
    noise sources. See Laning & Battin (1956), pp. 136–44 and Miller (1956), pp. 282–4.
