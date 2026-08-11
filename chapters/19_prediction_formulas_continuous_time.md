# 19. Prediction Formulas for Continuous Time Linear Rational Expectations Models

*by Lars Peter Hansen and Thomas J. Sargent*

In this note we derive optimal prediction formulas to be used in solving continuous time
rational expectations models. In these derivations we employ Laplace transforms in a manner
analogous to the use of $z$ transforms for solving discrete time optimal prediction problems
in Hansen and Sargent (1980a, Appendix A). The formulas are intended to play the same role
for continuous time models that the discrete time formulas for optimal predictions of
geometric distributed leads did in Hansen and Sargent (1980a). Within this book, they generalize
the continuous time geometric-distributed-lead forecast of {doc}`12_prediction` from
first-order Markov forcing to the full class of rational (and nonstationary) processes — making
good the promise of equation {eq}`eq-12-gen` there — and they supply the prediction calculus
that {doc}`14_faster_methods_recursive_linear_models` invokes to evaluate the optimal feedforward
{eq}`eq-14-forward`.

## 1. Convolutions and Prediction

Let $L^1$ and $L^2$ denote the spaces of all real-valued Borel measurable functions $\phi$ on
$\mathbb{R}$ that are absolutely integrable and square integrable, respectively. Let $W$ denote
a random measure defined on $\mathbb{R}$ with increments that are orthogonal and second-moment
stationary. In other words,

```{math}
:label: eq-19-1-1
E\left[ W\{[t_2, t_1)\}^2 \right] = t_2 - t_1 \quad \text{for } t_2 > t_1,
```

and

```{math}
:label: eq-19-1-2
E\left[ W\{[t_4, t_3)\}\, W\{[t_2, t_1)\} \right] = 0 \quad \text{for } t_4 > t_3 > t_2 > t_1.
```

Here $W\{[t_2, t_1)\}$ denotes the increment of $W$ over the interval between $t_1$ and $t_2$;
we follow the original in writing the larger endpoint first, so that {eq}`eq-19-1-1` assigns
that interval a variance equal to its length, and {eq}`eq-19-1-2` says that increments over
disjoint intervals are orthogonal.

Using functions in $L^2$ and the random measure $W$, we construct second-moment stationary
processes as convolutions:

```{math}
:label: eq-19-1-3
x(t) = \int_{-\infty}^{+\infty} \phi(\tau)\, dW(t - \tau).
```

The stochastic integral in {eq}`eq-19-1-3` can be interpreted as the limit point of a
mean-square convergent sequence of random variables (e.g. see Rozanov 1967). Relation
{eq}`eq-19-1-3` gives a convenient mapping between the space $L^2$ of functions and the space
$X$ of stochastic processes. It turns out that inner products on these two spaces coincide.
More precisely, let $\phi_1$ and $\phi_2$ be any two functions in $L^2$. An implication of
{eq}`eq-19-1-1` and {eq}`eq-19-1-2`, is

```{math}
:label: eq-19-1-4
\int_{-\infty}^{+\infty} \phi_1(\tau)\phi_2(\tau)\, d\tau = E\left[ x_1(t)\, x_2(t) \right]
```

where $x_1$ and $x_2$ are given by convolution {eq}`eq-19-1-3` using $\phi_1$ and $\phi_2$
respectively.

Let $L^2_+$ denote the subspace of $L^2$ consisting of all functions that are zero on
$(-\infty, 0)$ and let $L^2_-$ denote the subspace of all functions that are zero on
$[0, \infty)$. Clearly, $L^2_+$ and $L^2_-$ are orthogonal and $L^2 = L^2_- \oplus L^2_+$. Any
$\phi$ in $L^2$ can be decomposed uniquely into the sum of two functions $\phi^+ \in L^2_+$ and
$\phi^- \in L^2_-$ via:[^fn19-1]

```{math}
:label: eq-19-1-5
\begin{aligned}
\phi^{+}(t) &\equiv \begin{cases} \phi(t) & t \geq 0 \\ 0 & t < 0 \end{cases} \\[4pt]
\phi^{-}(t) &\equiv \begin{cases} 0 & t \geq 0 \\ \phi(t) & t < 0. \end{cases}
\end{aligned}
```

To formulate the prediction problems of interest, we use the random measure $W$ to induce a
family of information sets indexed by calendar time. Let $H(t)$ denote the space of random
variables $x(t)$ given by {eq}`eq-19-1-3` for $\phi$'s restricted to be in $L^2_+$. It follows
from {eq}`eq-19-1-4` that since $L^2_+$ is a Hilbert space, so is $H(t)$. Furthermore, the
family of Hilbert spaces $\{H(t)\}$ is increasing in the sense that if $t_2 > t_1$, then
$H(t_2) \supset H(t_1)$. Since $H(t)$ is constructed using the random measure $W$, the least
squares projection operator $P[\cdot \mid H(t)]$ onto the space $H(t)$ is given by

```{math}
:label: eq-19-1-6
P\left[ \int_{-\infty}^{+\infty} \phi(\tau)\, dW(t-\tau) \;\Big|\; H(t) \right]
= \int_{-\infty}^{+\infty} \phi^{+}(\tau)\, dW(t-\tau).
```

Hence the prediction process obtained by taking a process $x \in X$ constructed as a
convolution of $\phi$ and $dW$ and projecting it onto $H(t)$ for each $t$ is a convolution of
$\phi^+$ and $dW$ for $\phi^+$ given in {eq}`eq-19-1-5`.

## 2. Transforms

One convenient way to represent functions in $L^2$ and characterize mapping {eq}`eq-19-1-5`
involves the use of transforms. For instance, Fourier transforms are valuable in
characterizing the second moment properties of processes in $X$. For any $\phi$ in
$L^1 \cap L^2$, the Fourier transform of $\phi$ is defined to be

```{math}
:label: eq-19-2-1
\mathcal{F}t(\phi)(\theta) \equiv \int_{-\infty}^{+\infty} \exp(-i\theta t)\, \phi(t)\, dt.
```

There is a well known extension of $\mathcal{F}t$ from $L^1 \cap L^2$ to $L^2$. Using this
extension, the spectral density function for $x$ generated via {eq}`eq-19-1-3` is just
$\lvert \mathcal{F}t(\phi) \rvert^2$.

To characterize the implied second moment properties of the solutions to prediction problems
of the form {eq}`eq-19-1-6`, we use Laplace transforms. These transforms are defined as
follows. For any $\phi$ in $L^2$ and any $\rho$ in $\mathbb{R}$ we construct a new function
$\exp(-\rho t)\phi$. This new function may or may not be in $L^2$ depending on the value of
$\rho$. Whenever it is in $L^2$, we define the Laplace transform to be:

```{math}
:label: eq-19-2-2
\mathcal{L}p(\phi)(\mathbf{c}) \equiv \mathcal{F}t\left[ \exp(-\rho t)\phi \right](\theta)
```

where $\mathbf{c} \equiv \rho + i\theta$.

The question of interest is the following. Given the Laplace transform $\mathcal{L}p(\phi)$ of a
function $\phi \in L^2$, how can we compute or characterize $\mathcal{L}p(\phi^+)$ where
$\phi^+$ is defined in {eq}`eq-19-1-5`? To answer this question, we first study Laplace
transforms of functions $\phi \in L^2_+$. For any such $\phi$, $\exp(-\rho t)\phi$ is also in
$L^2_+$ as long as $\rho > 0$. Hence the Laplace transform $\mathcal{L}p(\phi)(\mathbf{c})$ is
well defined on the closed right plane $\mathbf{C}^+_0$ where
$\mathbf{C}^+_\delta \equiv \{\mathbf{c} \in \mathbf{C} : \operatorname{real}(\mathbf{c}) \geq \delta\}$.
Moreover, $\mathcal{L}p(\phi)$ is analytic in the interior of $\mathbf{C}^+_0$ (relative to
$\mathbf{C}$). For $\delta > 0$ and $\mathbf{c} \in \mathbf{C}^+_\delta$,

```{math}
:label: eq-19-2-3
\begin{aligned}
\lvert \mathcal{L}p(\phi)(\mathbf{c}) \rvert
&\leq \int_0^\infty \lvert \phi(t) \rvert \exp(-\delta t)\, dt \\
&\leq \left[ \int_0^\infty \lvert \phi(t) \rvert^2\, dt \int_0^\infty \exp(-2\delta t)\, dt \right]^{1/2} \\
&\leq \left[ \int_0^\infty \lvert \phi(t) \rvert^2\, dt \,/\, 2\delta \right]^{1/2}
\end{aligned}
```

where the second inequality is an application of the familiar Cauchy–Schwarz inequality. The
right side of {eq}`eq-19-2-3` gives a uniform bound (in $\mathbf{c}$) on $\mathcal{L}p(\phi)$
over the set $\mathbf{C}_{\delta}^+$. This bound becomes arbitrarily small as $\delta$ tends to
plus infinity.

Consider next functions $\phi \in L^2_-$. For any such $\phi$, $\mathcal{L}p(\phi)(\mathbf{c})$
is always well defined for $\mathbf{c}$ in the left half plane
$\mathbf{C}_0^- \equiv \{\mathbf{c} \in \mathbf{C} : \operatorname{real}(\mathbf{c}) \leq 0\}$,
and $\mathcal{L}p(\phi)$ is analytic in the interior of that domain. Define
$\mathbf{C}_\delta^- \equiv \{\mathbf{c} \in \mathbf{C} : \operatorname{real}(\mathbf{c}) \leq \delta\}$.
Mimicking the previous argument, it can be shown that for any $\delta < 0$, $\mathcal{L}p(\phi)$
is bounded on the domain $\mathbf{C}_\delta^-$ and that the bound can be made arbitrarily small
by driving $\delta$ towards minus infinity.

For general functions $\phi$ in $L^2$, $\mathcal{L}p(\phi)$ may only be defined on the imaginary
axis, i.e. for $\operatorname{real}(\mathbf{c}) = 0$. We are interested in a smaller class of
functions, however. Let $\Phi$ be the set of all functions $\phi \in L^2$ such that
$\mathcal{L}p(\phi^-)$ is analytic in the interior of a region $\mathbf{C}_{\delta}^-$ for some
$\delta > 0$. In this case $\mathcal{L}p(\phi)$ is analytic in the interior of the strip
$\mathbf{C}_{\delta}^- \cap \mathbf{C}_0^+$. Furthermore, for any closed interval
$J \subset (0, \delta)$, $\mathcal{L}p(\phi)$ is bounded on
$\{\mathbf{c} \in \mathbf{C} : \operatorname{real}(\mathbf{c}) \in J\}$. Define $\mathcal{A}$ to
be the collection of all Laplace transforms of functions $\phi \in \Phi$.

The following result gives the decomposition for $a \in \mathcal{A}$ corresponding to the
decomposition $\phi = \phi^+ + \phi^-$.

```{admonition} Lemma
:class: note

For any $a \in \mathcal{A}$ there is a unique decomposition $a = a^+ + a^-$ where

(i) $a^-$ is analytic in the interior of $\mathbf{C}_\delta^-$, uniformly bounded on any closed
half plane $\mathbf{C}_\rho^-$ for $\rho < \delta$ and

$$
\lim_{\rho \to -\infty} \max_{\mathbf{c} \in \mathbf{C}_{\rho}^{-}} \lvert a^{-}(\mathbf{c}) \rvert = 0\,;
$$

(ii) $a^+$ is analytic in the interior of $\mathbf{C}_0^+$, uniformly bounded on any closed half
plane $\mathbf{C}_\rho^+$ for any $\rho > 0$.
```

**Proof.** Functions $a^-$ and $a^+$ satisfying (i) and (ii) are obtained by letting
$a^- = \mathcal{L}p(\phi^-)$ and $a^+ = \mathcal{L}p(\phi^+)$. To show that the decomposition is
unique, we let $a = b^+ + b^-$ be any other decomposition where $b^-$ satisfies (i) and $b^+$
satisfies (ii). Note that

$$
a^+ - b^+ = b^- - a^-
$$

at least in the interior of the strip $\mathbf{C}_{\delta}^- \cap \mathbf{C}_0^+$. Since
$a^+ - b^+$ is analytic in the interior of $\mathbf{C}_0^+$ and $b^- - a^-$ is analytic in the
interior of $\mathbf{C}_{\delta}^-$, $b^- - a^-$ can be extended to be analytic on all of
$\mathbf{C}$. Furthermore, the uniform bounds on $a^+ - b^+$ and $b^- - a^-$ on overlapping half
planes ensure that the extension of $b^- - a^-$ is bounded as well. The only functions that are
bounded and analytic on $\mathbf{C}$ are constant. Since $a^-$ and $b^-$ satisfy (i) and the
extension of $b^- - a^-$ to $\mathbf{C}$ is constant, $a^+ - b^+$ must be identically zero.
$\blacksquare$

Decompositions like that given in the Lemma apply to a much more general collection of analytic
functions than the Laplace transforms of functions in $L^2$. For instance, they also apply to
Laplace transforms of generalized functions (e.g. see Beltrami and Wohlers 1966). However,
these more general decompositions may not be unique. For instance, suppose we ignore the
requirement

$$
\lim_{\rho \to -\infty} \max_{\mathbf{c} \in \mathbf{C}_{\rho}^{-}} \lvert a^{-}(\mathbf{c}) \rvert = 0
$$

in (i) of the Lemma. Then one can always add complex numbers to $a^+$ and subtract the same
numbers from $a^-$ to obtain other decompositions of $a$. If in addition, we ignore the bound
restrictions in (i) and (ii) of the Lemma then one can add functions, such as polynomials, that
are analytic in the entire complex plane to $a^+$ and subtract them from $a^-$ to obtain other
decompositions of $a$. Therefore in applying the Lemma to compute $\mathcal{L}p(\phi^+)$, it is
important to check whether the candidates for $\mathcal{L}p(\phi^+)$ and $\mathcal{L}p(\phi^-)$
satisfy the bounds restrictions in (i) and (ii).

## 3. Examples

We now apply the Lemma to obtain frequency domain characterizations of the solutions to
prediction problems that occur in rational expectations models. These problems all have the
following structure. Let $\psi \in L^2_+$, and define $y(t)$ by the convolution:

$$
y(t) = \int_0^{+\infty} \psi(\tau)\, dW(t - \tau).
$$

Construct a new process by forming a forward-looking convolution using a function
$\gamma \in L^1$:[^fn19-2]

```{math}
:label: eq-19-3-1
x(t) = \int_{-\infty}^{+\infty} \gamma(\tau)\, y(t-\tau)\, d\tau
\equiv \int_{-\infty}^{+\infty} \phi(\tau)\, dW(t-\tau)
```

where $\phi$ is given by the convolution:

$$
\phi(\tau) = \int_{-\infty}^{+\infty} \gamma(s)\, \psi(\tau - s)\, ds.
$$

Applying the well known product representation for Fourier transforms of convolutions, we have
that

$$
\mathcal{F}t(\phi) = \mathcal{F}t(\psi)\, \mathcal{F}t(\gamma).
$$

This same result extends to Laplace transforms on the common domain of $\mathcal{L}p(\phi)$ and
$\mathcal{L}p(\gamma)$. For the examples we consider, there will exist a $\delta > 0$ such that
$\mathcal{L}p(\gamma)$ is defined on the interior of $\mathbf{C}_{\delta}^-$. Hence on the
interior of the strip $\mathbf{C}_{\delta}^- \cap \mathbf{C}_{0}^+$,

```{math}
:label: eq-19-3-2
\mathcal{L}p(\phi) = \mathcal{L}p(\psi)\, \mathcal{L}p(\gamma).
```

We now investigate three related examples.

### Example 1

Suppose that $\gamma$ is given by

```{math}
:label: eq-19-3-3
\gamma(t) = \begin{cases} 0 & t \geq 0 \\ \exp(\delta t) & t < 0 \end{cases}.
```

Then

$$
\mathcal{L}p(\gamma)(\mathbf{c}) = \int_{-\infty}^{0} \exp[(\delta - \mathbf{c})t]\, dt = 1/(\delta - \mathbf{c})
$$

for $\operatorname{real}(\mathbf{c}) < \delta$. Thus

$$
\mathcal{L}p(\phi)(\mathbf{c}) = \mathcal{L}p(\psi)(\mathbf{c})/(\delta - \mathbf{c}).
$$

Note that $\mathcal{L}p(\phi)$ is analytic on $\mathbf{C}$ except possibly at the point $\delta$
where it may have a pole. If $\mathcal{L}p(\psi)(\delta)$ is zero, the singularity at $\delta$ is
removable and $\mathcal{L}p(\phi^+) = \mathcal{L}p(\phi)$. Usually $\mathcal{L}p(\phi)$ will have
a pole at $\delta$, and to compute $\mathcal{L}p(\phi^+)$ we must eliminate this pole. One
candidate for $\mathcal{L}p(\phi^+)$ is

$$
a^{+}(\mathbf{c}) = [\mathcal{L}p(\psi)(\delta) - \mathcal{L}p(\psi)(\mathbf{c})]/(\mathbf{c} - \delta).
$$

Notice that the singularity of $a^+$ at $\delta$ is removable. The corresponding choice of $a^-$
is

$$
a^{-}(\mathbf{c}) = a(\mathbf{c}) - a^{+}(\mathbf{c}) = \mathcal{L}p(\psi)(\delta)/(\delta - \mathbf{c}).
$$

It is straightforward to show that $a^+$ and $a^-$ satisfy the requirements of the Lemma.
Therefore,

```{math}
:label: eq-19-3-4
\mathcal{L}p(\phi^{+}) = [\mathcal{L}p(\psi)(\delta) - \mathcal{L}p(\psi)(\mathbf{c})]/(\mathbf{c} - \delta).
```

Formula {eq}`eq-19-3-4` is the continuous time counterpart to formula (5) in Hansen and Sargent
(1980a). It is the operator that evaluates the geometric distributed lead {eq}`eq-12-glead` of
{doc}`12_prediction` and the optimal feedforward {eq}`eq-14-forward` of
{doc}`14_faster_methods_recursive_linear_models`. Indeed, the $\gamma$ of {eq}`eq-19-3-3`
generates precisely the geometric distributed lead of {doc}`12_prediction`: substituting
$u = -\tau$ in {eq}`eq-19-3-1` gives
$x(t) = \int_0^\infty e^{-\delta u} y(t+u)\, du$, so that $\delta = -\rho$ in the notation
there. With that substitution {eq}`eq-19-3-4` reads
$[-\tilde P(\mathbf{c}) + \tilde P(-\rho)]/(\mathbf{c} + \rho)$, which is exactly
{eq}`eq-12-gen`. {doc}`12_prediction` derived that formula only for a first-order
$\tilde P$ and then asserted it for general $\tilde P$; the Lemma above is what establishes
the general case, and in Example 3 below it is turned into a finite recursion.

### Example 2

More generally, suppose

$$
\mathcal{L}p(\gamma)(\mathbf{c}) = p_n(\mathbf{c})/p_d(\mathbf{c})
$$

where $p_n$ and $p_d$ are finite-order polynomials with real coefficients. To ensure that
$p_n(\mathbf{c})/p_d(\mathbf{c})$ is the Laplace transform of a function in $L^2_-$, we assume
that the order of $p_d$ exceeds the order of $p_n$ and that the zeros of $p_d$ are in the
interior of $\mathbf{C}_0^+$. In this case

$$
\mathcal{L}p(\phi)(\mathbf{c}) = \mathcal{L}p(\psi)(\mathbf{c})\, p_n(\mathbf{c})/p_d(\mathbf{c}),
$$

which has poles in the interior of $\mathbf{C}_0^+$ only at the zeros of $p_d$. Let $a_j$ denote
the principal part of the Laurent series expansion of $\mathcal{L}p(\phi)(\mathbf{c})$ at the
$j^{\text{th}}$ zero of $p_d$. It follows from the partial fractions decomposition of a
meromorphic function that

$$
a^{+} = \mathcal{L}p(\phi) - \sum_{j} a_{j}
$$

is analytic in the interior of $\mathbf{C}_0^+$. Furthermore, the principal parts, $a_j$, are
each sums of reciprocals of first and higher-order polynomials and hence satisfy

$$
\lim_{\rho \to -\infty} \max_{\mathbf{c} \in \mathbf{C}_{\rho}^{-}} \lvert a_{j}(\mathbf{c}) \rvert = 0
$$

for each $j$. By construction,

$$
a^- = \sum_j a_j
$$

satisfies (i) of the Lemma where $\delta$ is the real part of the zero of $p_d$ closest to the
imaginary axis and is bounded on $\mathbf{C}_{\rho}^-$ for any $\rho < \delta$. Therefore, we
have the following generalization of {eq}`eq-19-3-4`:

$$
\mathcal{L}p(\phi^{+}) = \mathcal{L}p(\phi) - \sum_{j} a_{j}.
$$

### Example 3

Suppose that $\gamma$ is given by {eq}`eq-19-3-3`, and $\mathcal{L}p(\psi)$ is a rational
function:

$$
\mathcal{L}p(\psi) = q_n/q_d
$$

where $q_n$ and $q_d$ are polynomials with real coefficients. To guarantee that $q_n/q_d$ is the
Laplace transform of a function in $L^2_+$, we assume that the order of $q_d$ exceeds the order
of $q_n$ and that the zeros of $q_d$ are in the interior of $\mathbf{C}_0^-$. Solution
{eq}`eq-19-3-4` now becomes

$$
\mathcal{L}p(\phi)(\mathbf{c}) = q_n(\mathbf{c})/[q_d(\mathbf{c})(\delta - \mathbf{c})].
$$

From Example 1, we know that

```{math}
:label: eq-19-3-5
\begin{aligned}
\mathcal{L}p(\phi^{+})(\mathbf{c})
&= [q_{n}(\delta)/q_{d}(\delta) - q_{n}(\mathbf{c})/q_{d}(\mathbf{c})]/(\mathbf{c} - \delta) \\
&= [q_{n}(\delta)q_{d}(\mathbf{c}) - q_{n}(\mathbf{c})q_{d}(\delta)]/[q_{d}(\mathbf{c})q_{d}(\delta)(\mathbf{c} - \delta)].
\end{aligned}
```

The right side of {eq}`eq-19-3-5` has a removable singularity at $\delta$ by construction. This
is evident because the polynomial $[q_n(\delta)q_d(\mathbf{c}) - q_n(\mathbf{c})q_d(\delta)]$ has
a zero at $\delta$. Canceling the common factor $(\mathbf{c} - \delta)$ in the numerator and
denominator results in

$$
\mathcal{L}p(\phi^+)(\mathbf{c}) = q_n^+(\mathbf{c})/q_d^+(\mathbf{c})
$$

where

$$
q_d^+(\mathbf{c}) = q_d(\mathbf{c})q_d(\delta)
$$

and $q_n^+$ satisfies

```{math}
:label: eq-19-3-6
q_n^+(\mathbf{c})(\mathbf{c} - \delta) = [q_n(\delta)q_d(\mathbf{c}) - q_n(\mathbf{c})q_d(\delta)].
```

By equating coefficients of the polynomials on both sides of {eq}`eq-19-3-6`, one can construct
a linear system of equations in the coefficients of $q_n^+(\mathbf{c})$. In fact there is a
recursive structure to this equation system that can be exploited as follows. Let $\eta_j$
denote the coefficient on $\mathbf{c}^j$ in
$[q_n(\delta)q_d(\mathbf{c}) - q_n(\mathbf{c})q_d(\delta)]$ and let $\epsilon_j$ denote the
corresponding coefficient in $q_n^+(\mathbf{c})$. Then

$$
-\delta\epsilon_0 = \eta_0
$$

and

$$
\epsilon_{j-1} - \delta \epsilon_j = \eta_j \quad \text{for } j \geq 1
$$

which can be solved recursively beginning with $\epsilon_0$. The solution to this recursion
gives a continuous time counterpart to formulas reported in Hansen and Sargent (1980a, 1981b)
for autoregressive and autoregressive moving-average processes.

## 4. Vector Information Structures

Suppose that $W$ is a $k$-dimensional vector random measure with second moment stationary
increments. We now replace {eq}`eq-19-1-1` and {eq}`eq-19-1-2` with

$$
E\left[ W\{[t_2, t_1)\}\, W\{[t_2, t_1)\}' \right] = (t_2 - t_1)I_k \quad \text{for } t_2 > t_1,
$$

and

$$
E\left[ W\{[t_4, t_3)\}\, W\{[t_2, t_1)\}' \right] = 0 \quad \text{for } t_4 > t_3 > t_2 > t_1
$$

where $I_k$ is a $k$-dimensional identity matrix. Processes in $X$ are now constructed using a
$k$-dimensional vector $\phi$ of functions in $L^2$ via:

$$
x(t) = \int_{-\infty}^{+\infty} \phi(\tau) \cdot dW(t - \tau).
$$

The analyses in Sections 2 and 3 extend by applying the decompositions to each of the $k$
Laplace transforms of entries in $\phi$. In Example 1 formula {eq}`eq-19-3-4` still applies where
$\mathcal{L}p(\psi)$ is the vector of Laplace transforms of entries in $\psi$. The recursions
derived in Example 3 still apply where $q_n$ is now a $k$-dimensional vector of polynomials, each
with orders less than the scalar polynomial $q_d$.

## 5. Nonstationarities

In Section 3, the assumption that $\psi \in L_+^2$ guaranteed that process $y$ is second moment
stationary. Our analysis can be extended to a more general class of processes, however. To
accommodate nonstationarities, it is most convenient to think of the underlying information
process as starting at some initial time, say $t = 0$. Hence we imagine {eq}`eq-19-1-1` and
{eq}`eq-19-1-2` holding for nonnegative values of $t_1$, $t_2$, $t_3$ and $t_4$, and we assume
that the random measure of any interval contained in $(-\infty, 0)$ is zero. This permits
formula {eq}`eq-19-3-1` to be well defined for a much larger class of functions $\psi$. We might
view the process $y$ as being the deviation from a path that is perfectly predictable from time
zero forward. We impose the weaker requirement that $\exp(-\rho t)\psi$ be in $L^2$ for strictly
positive values of $\rho$ which allows for polynomial growth in the second moment of $y$.[^fn19-3]
The calculations in Examples 1 through 3 still apply. In the case of Example 3, to accommodate
polynomial growth we now allow $q_d$ to have zeros on the imaginary axis of the complex plane
$\mathbf{C}$.

## Notes

[^fn19-1]: The uniqueness of this decomposition requires some qualification. Elements of $L^2$
    are only defined up to an equivalence of functions that are equal almost everywhere. Hence
    from the vantage point of $L^2$, the construction of $\phi^+$ and $\phi^-$ at a particular
    point, say $t = 0$, is inconsequential.

[^fn19-2]: We take the right side of equation {eq}`eq-19-3-1` as the definition of $x(t)$.
    Alternatively, for particular classes of $\gamma$ we could define $x(t)$ using finite sum
    approximations for the middle integral.

[^fn19-3]: In Examples 1 and 3 it is also possible to allow for exponential growth in the second
    moments of $y$ as long as $\psi \exp(-\sigma t)$ is in $L^2_+$ for some $\sigma$ satisfying
    $0 < \sigma < \delta$. The transform analysis now applies to the narrower strip
    $\mathbf{C}^-_{\delta} \cap \mathbf{C}^+_{\sigma}$.
