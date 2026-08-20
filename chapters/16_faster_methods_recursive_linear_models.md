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

# 16. Faster Methods for Solving Recursive Linear Models of Dynamic Economies

```{eval-rst}
.. index::
   single: matrix sign algorithm
   single: Hamiltonian matrix
   single: certainty equivalence
   single: feedback and feedforward
   single: Arrow-Debreu prices
   single: social planning problem
```

*An application and extension of the continuous-time linear regulator of {doc}`15_kalman_filter_spectral_factorization`.*

```{note}
This chapter paraphrases and streamlines Lars Peter Hansen, John Heaton, and Thomas J.
Sargent, "Faster Methods for Solving Continuous Time Recursive Linear Models of Dynamic
Economies," Chapter 7 of *Rational Expectations Econometrics* (L. P. Hansen and T. J.
Sargent, eds., Westview Press, Underground Classics in Economics). We present it as an
application of the optimal linear regulator and algebraic Riccati equation of
{doc}`15_kalman_filter_spectral_factorization`.
```

Section (d) of {doc}`15_kalman_filter_spectral_factorization` posed the continuous-time linear
regulator, minimizing $\int_0^\infty (x^\top Q x + u^\top R u)\, dt$ subject to
$\dot x = A x + B u$, and solved it through an algebraic Riccati equation. That regulator computes the
equilibria of a large class of dynamic stochastic economies. Hansen and Sargent's research
program associates with a competitive economy a fictitious social planning problem with a
quadratic objective and linear constraints. The planner's allocation coincides with the
competitive equilibrium allocation, and the gradient of the planner's value function gives the
Arrow–Debreu prices that support it. That gradient is linear in the state, because the value
function is quadratic. Computing the equilibrium therefore reduces to solving an optimal linear
regulator.

This chapter develops that reduction and then extends the basic regulator in three ways that
matter for economics:

1. It allows the forcing processes, the preference and endowment shocks, to be **arbitrary**
   square-integrable functions of time rather than outputs of a finite system of differentiable
   equations. The constraints become *nonautonomous*. The optimal decision rule then
   **separates** into a *feedback* (control) part and a *feedforward* (prediction) part. That is
   the continuous-time form of the certainty-equivalence separation that Lucas and Sargent
   (1981) emphasized. The feedforward part is exactly the geometric
   distributed lead of {doc}`12_prediction`, computed in general by the companion
   {doc}`19_prediction_formulas_continuous_time`.
2. It supplies a fast solver for the Riccati equation, the *matrix sign algorithm*, that avoids
   eigenvalue decompositions and exploits the Hamiltonian structure. It scales to the large
   state vectors that realistic economies require.
3. It connects the solution to the **spectral factorization theorem** of
   {doc}`08_spectral_densities`: the Euler equations of the model are factored exactly as a
   rational spectral density is factored.

```{eval-rst}
.. index::
   single: linear quadratic economy
   single: social planning problem; fictitious
   single: competitive equilibrium; and the planning problem
   single: value function; gradient of
```

## (a) A linear-quadratic economy as a planning problem

A representative household has time-separable preferences over an $n_s$-vector $s(t)$ of
*services*,

```{math}
:label: eq-16-pref
-\tfrac12 \int_0^\infty e^{-\rho t}\, [s(t) - b(t)]^\top [s(t) - b(t)]\, dt, \qquad \rho > 0,
```

where $b(t)$ is an exogenous *bliss* (preference-shock) process and $\rho$ is the subjective
discount rate. Services are produced from a stock of *household capital* $h(t)$ and a vector of
*consumption goods* $c(t)$,

```{math}
:label: eq-16-services
s(t) = \Lambda\, h(t) + \Pi\, c(t),
```

and the two capital stocks, household capital $h(t)$ and physical capital $k(t)$, accumulate
according to linear laws of motion driven by consumption and by investment $i(t)$:

```{math}
:label: eq-16-capital
D h(t) = \Delta_h\, h(t) + \Theta_h\, c(t), \qquad
D k(t) = \Delta_k\, k(t) + \Theta_k\, i(t),
```

with $D$ the time-derivative operator of the earlier chapters. Output is produced from physical
capital and an endowment $f(t)$ and split between consumption and investment through a resource
constraint

```{math}
:label: eq-16-resource
\Phi_c\, c(t) + \Phi_i\, i(t) = \Gamma\, k(t) + f(t).
```

The matrices $\Lambda, \Pi, \Delta_h, \Theta_h, \Delta_k, \Theta_k, \Phi_c, \Phi_i, \Gamma$ are
the deep parameters of preferences and technology; $\Lambda$ and $\Pi$ govern how durable and
how substitutable consumption goods are over time, so the model accommodates habit persistence
and durability. The endowment and bliss shocks are linear images of an $n_z$-vector forcing
process $\hat z(t)$,

```{math}
:label: eq-16-forcing
f(t) = \Xi_f\, \hat z(t), \qquad b(t) = \Xi_b\, \hat z(t).
```

Stacking the capital stocks into a composite state $\hat x(t) = [\,h(t)^\top\ k(t)^\top\,]^\top$
and taking investment as the control $\hat u(t) = i(t)$, the laws of motion
{eq}`eq-16-capital` together with the resource constraint {eq}`eq-16-resource` collapse to a
single **linear state-space system** of exactly the form met in
{doc}`15_kalman_filter_spectral_factorization`,

```{math}
:label: eq-16-statehat
D \hat x(t) = \hat A\, \hat x(t) + B_u\, \hat u(t) + B_z\, \hat z(t), \qquad \hat x(0) = \mu,
```

and the preference criterion {eq}`eq-16-pref`, after substituting out $s(t) - b(t)$, becomes a
quadratic form in $(\hat u, \hat x, \hat z)$,

```{math}
:label: eq-16-objhat
-\tfrac12 \int_0^\infty e^{-\rho t}\,
\begin{bmatrix} \hat u \\ \hat x \\ \hat z \end{bmatrix}^\top
\Omega
\begin{bmatrix} \hat u \\ \hat x \\ \hat z \end{bmatrix} dt,
\qquad \Omega \succeq 0 .
```

The planning problem is to choose $\{\hat u(t)\}$ to maximize {eq}`eq-16-objhat` subject to
{eq}`eq-16-statehat`. This is a discounted optimal linear regulator with an exogenous forcing
term $B_z \hat z$. It is the regulator of Chapter 15 augmented by a forcing process.

```{eval-rst}
.. index::
   single: discounting; removal of
   single: undiscounted transformation
```

## (b) Removing the discount

The discount factor is absorbed into the state by a change of variables. With $\epsilon = \rho/2$,
define the deflated variables

```{math}
:label: eq-16-deflate
x(t) = e^{-\epsilon t}\, \hat x(t), \qquad u(t) = e^{-\epsilon t}\, \hat u(t), \qquad
z(t) = e^{-\epsilon t}\, \hat z(t).
```

Because $D x = -\epsilon x + e^{-\epsilon t} D\hat x$, the state equation {eq}`eq-16-statehat`
becomes the **undiscounted** system

```{math}
:label: eq-16-state
D x(t) = A\, x(t) + B_u\, u(t) + B_z\, z(t), \qquad A \equiv \hat A - \epsilon I, \quad x(0) = \mu,
```

and the criterion {eq}`eq-16-objhat` becomes the undiscounted quadratic form

```{math}
:label: eq-16-obj
-\tfrac12 \int_0^\infty
\begin{bmatrix} u \\ x \\ z \end{bmatrix}^\top \Omega
\begin{bmatrix} u \\ x \\ z \end{bmatrix} dt .
```

Deflation has simply shifted every eigenvalue of the transition matrix left by $\epsilon$;
provided the original capital-depreciation matrices have eigenvalues with real parts below
$\epsilon$, the deflated transition matrix $A$ is a stability matrix and the deflated stocks are
square integrable. We have reduced the economy to a **standard undiscounted optimal linear
regulator with forcing**. It is the problem of Chapter 15, Section (d), with the forcing term
$B_z z$ and a general cross-penalized criterion $\Omega$.

```{eval-rst}
.. index::
   single: invariant subspace; stable
   single: linear regulator; Hamiltonian of
```

## (c) The optimal linear regulator and its Hamiltonian

Partition $\Omega$ conformably with $(u, x, z)$ into blocks $\Omega_{uu}, \Omega_{ux}, \ldots$
The regulator differs from the textbook problem of Chapter 15 in two ways that the economic
application forces upon us. First, the forcing $z(t)$ is allowed to be *any* square-integrable
function. It need **not** solve a differential equation, so the constraints are nonautonomous.
Second, stability of $x$ is imposed as an explicit constraint rather than emerging
automatically, because with general $z$ a stable state path is not guaranteed.

Attaching a co-state (Lagrange multiplier) process $\lambda(t)$ to the state equation
{eq}`eq-16-state` and forming the Lagrangian, the first-order conditions in $u$ and $x$, after
eliminating the control, combine the state and co-state into a single first-order system
governed by a **Hamiltonian matrix** $H$:

```{math}
:label: eq-16-ham
\begin{bmatrix} D x_s(t) \\ D x_c(t) \end{bmatrix}
= H \begin{bmatrix} x_s(t) \\ x_c(t) \end{bmatrix} + K\, z(t),
\qquad
H = \begin{bmatrix} H_{11} & H_{12} \\ H_{21} & H_{22} \end{bmatrix},
```

where $x_c$ is the (transformed) co-state. The blocks are built from $A$, $B_u$ and the
partitions of $\Omega$, for instance $H_{11} = A - B_u \Omega_{uu}^{-1}\Omega_{ux}$ and
$H_{12} = -B_u \Omega_{uu}^{-1} B_u^\top$. $H$ has the defining **Hamiltonian symmetries**

$$
H_{22} = -H_{11}^\top, \qquad H_{12} = H_{12}^\top, \qquad H_{21} = H_{21}^\top .
$$

These symmetries are the structural fingerprint of an optimization problem; they are what make
the eigenvalues of $H$ occur in pairs $\pm\zeta$ symmetric about the imaginary axis, exactly $n$
of them with strictly negative real parts.

**Solving the homogeneous ($z \equiv 0$) problem.** Following Vaughan (1969), take the
eigendecomposition $H = E J E^{-1}$ with the $n$ stable eigenvalues in the leading block of $J$.
Stability of the deflated state requires the solution to lie in the **stable invariant
subspace** of $H$, which forces the co-state to be a linear function of the state,

```{math}
:label: eq-16-Mx
x_c(t) = M_x\, x_s(t), \qquad M_x = E_{21} E_{11}^{-1},
```

where $E_{11}, E_{21}$ are blocks of the stable eigenvectors. The matrix $M_x$ is symmetric
and positive semidefinite, and it is the stabilizing solution of an algebraic Riccati
equation, the same object as the cost-to-go matrix $P$ of the regulator in Chapter 15. It satisfies

```{math}
:label: eq-16-riccati
[\,-M_x \quad I\,]\; H \begin{bmatrix} I \\ M_x \end{bmatrix} = 0 ,
```

which is the Riccati equation {eq}`eq-15-careg` written in Hamiltonian form. The optimal value
of the planning criterion, as a function of the initial condition, is the quadratic
$-\tfrac12\, \mu^\top M_x\, \mu$. Its gradient $-M_x \mu$ delivers the supporting Arrow–Debreu
prices.

**The feedback–feedforward decomposition.** Substituting {eq}`eq-16-Mx` back into the
Hamiltonian system {eq}`eq-16-ham` gives a recursive law of motion for the optimally controlled
state,

```{math}
:label: eq-16-decision
D x_s(t) = \underbrace{(H_{11} + H_{12} M_x)\, x_s(t)}_{\text{feedback (control)}}
\; + \; \underbrace{H_{12}\, w(t) + K_1\, z(t)}_{\text{feedforward (prediction)}} .
```

The **feedback** matrix $H_{11} + H_{12} M_x$ is a stability matrix. It is the closed-loop
transition $A - B_u F$ of the Chapter 15 regulator, with feedback gain $F$ built from $M_x$. The
**feedforward** term carries the influence of the forcing process through an auxiliary process
$w(t)$ that depends on the *expected future* path of $z$. The solution thus splits into a part
that stabilizes the state (control) and a part that forecasts the shocks (prediction). That is
the separation principle in continuous time.

```{eval-rst}
.. index::
   single: matrix sign function
   single: quadratic convergence
   single: Riccati equation; numerical solution of
```

## (d) Faster Riccati solvers: the matrix sign algorithm

Vaughan's eigenvalue method computes $M_x$ from a full eigendecomposition of $H$, which is
expensive and numerically delicate when state vectors are large or eigenvalues nearly coincide.
The "faster methods" of the title replace it with the **matrix sign algorithm** of Roberts
(1971) and Denman and Beavers (1976), which computes $M_x$ by a simple, quadratically convergent
iteration that never forms the eigenvectors.

For any matrix $G$ with no purely imaginary eigenvalues, the iteration

```{math}
:label: eq-16-sign
\mathcal{R}(G) = \tfrac12\big(G + G^{-1}\big)
```

converges, when started at $G_0 = H$, to the matrix sign function $\operatorname{sign}(H)$,
a matrix with eigenvalues $-1$ on the stable invariant subspace and $+1$ on the unstable one.
Partitioning the limit $H^\infty = \operatorname{sign}(H)$ conformably, the stabilizing Riccati
solution is recovered algebraically; equivalently, $M_x$ is the least-squares solution of the
overdetermined system

```{math}
:label: eq-16-signextract
\begin{bmatrix} H^\infty_{12} \\ H^\infty_{22} + I \end{bmatrix} M_x
= - \begin{bmatrix} H^\infty_{11} + I \\ H^\infty_{21} \end{bmatrix}.
```

Anderson (1978) observed that the Hamiltonian structure ($H_{22} = -H_{11}^\top$) is *preserved*
by the iteration {eq}`eq-16-sign`, so the algorithm can be carried out on the half-size blocks,
roughly halving the work. The matrix sign algorithm is the continuous-time analogue of the
doubling algorithms used for discrete-time quadratic control. It solves the *same* algebraic
Riccati equation as `scipy`'s `solve_continuous_are` used in Chapter 15, but by an iteration that
scales gracefully, which is the point of the paper. {ref}`ex-16-1` implements it and checks it
against
the Chapter 15 solver.

**The feedforward as prediction.** It remains to compute the feedforward process $w(t)$. When
the forcing is itself autonomous, with $D z(t) = A_{zz} z(t)$ and $A_{zz}$ stable, the
feedforward collapses to a linear function $w(t) = M_z z(t)$, where $M_z$ solves a Sylvester equation linear
in its entries,

```{math}
:label: eq-16-sylvester
(H_{12} M_x + H_{11})^\top M_z + M_z A_{zz} = K_2 - M_x K_1 ,
```

solvable directly or by a second application of the matrix sign algorithm to a triangular
composite system. For *general* (nonautonomous) $z$, the feedforward is a forward convolution,
a present value of the expected future forcing,

```{math}
:label: eq-16-forward
w(t) = \int_0^\infty \exp\!\big[(H_{11} + H_{12} M_x)^\top \tau\big]\,(K_2 - M_x K_1)\, z(t+\tau)\, d\tau ,
```

in which the *transpose of the stable feedback matrix* governs the decay of the weights on
future shocks. (Equation {eq}`eq-16-forward` is written for a $z$ path known at $t$, the
perfect-foresight, or deterministic, case. When $z$ is stochastic, $z(t+\tau)$ is replaced
throughout by its least squares forecast $\hat E_t\, z(t+\tau)$; certainty equivalence, which
holds because the criterion is quadratic and the constraints linear, is what makes this
substitution legitimate and leaves the feedback matrix $H_{11} + H_{12}M_x$ unchanged.)
Equation {eq}`eq-16-forward` is exactly a **geometric distributed lead** of the
kind introduced in {doc}`12_prediction`: the optimal control today responds to a discounted
integral of expected future shocks. Evaluating it for a general forcing process is the task of
the companion {doc}`19_prediction_formulas_continuous_time`, which supplies the Laplace-transform
prediction calculus that turns {eq}`eq-16-forward` into a finite matrix computation.

```{eval-rst}
.. index::
   single: Euler equations; spectral factorization of
   single: geometric distributed lead; in decision rules
   single: rational spectral density matrix
```

## (e) Spectral factorization in the Euler equations

Take a pure investment-with-adjustment-cost economy: no household capital, so
$s(t) = \Pi c(t)$, with costs of adjusting capital that may penalize higher derivatives of the
capital stock. The first-order conditions can then be written as continuous-time **Euler
equations** in the capital stock. Disentangling the discount scaling, these take the
form

$$
[\,Q(D-\epsilon)\, \Omega_{11}\, P(D-\epsilon)\,]\, \hat y(t) = -\, Q(D - \epsilon)\, \Omega_{12}\, \hat z(t),
$$

where $P(\cdot)$ is a backward derivative operator and $Q(\cdot)$ a forward convolution operator.
The operator on the left is a polynomial in $D$ whose symbol, evaluated on the imaginary axis,
is a **rational spectral density matrix** $F(\theta)$. Solving the Euler equation requires
factoring $F$ into a part with all its zeros in the left half plane and its mirror image. That
applies the **spectral factorization theorem** of {doc}`08_spectral_densities`,

$$
F(\theta) = \hat P(-i\theta)^\top\, V\, \hat P(i\theta) \big/ \big[(-i\theta+\epsilon)(i\theta+\epsilon)\big]^{\ell+1},
$$

with $\hat P$ nonsingular in the left half plane. The stable factor furnishes a one-sided
forward inverse, and a matrix partial-fractions expansion of it yields the decision rule. Thus
the regulator solution of Section (c) and the spectral-factorization route of this section are
two faces of the same computation, the time-domain Riccati equation and the frequency-domain
factorization that Chapter 15 showed to be equivalent.

A particularly clean instance is Heaton's (1989) model of consumer durables, in which household
capital is the only endogenous state and the optimal net investment in the durable stock obeys

```{math}
:label: eq-16-heaton
D k(t) = [\,f(t) - b(t)\,] - \rho \int_0^\infty e^{-\rho\tau}\, [\,f(t+\tau) - b(t+\tau)\,]\, d\tau .
```

Net investment compares the *current* endowment-relative-to-bliss gap with a $\rho$-discounted
**present value of the future** gaps. A permanent-income logic governs it: the durable stock is
accumulated whenever the present is favorable relative to what is expected to come. The present
value on the right is, once again, the geometric distributed lead of {doc}`12_prediction`.
{ref}`ex-16-2` works this rule out for an autoregressive forcing process.

## Exercises

```{code-cell} ipython3
import numpy as np
from scipy.linalg import solve_continuous_are
import matplotlib.pyplot as plt
```

```{exercise-start}
:label: ex-16-1
```

**The matrix sign algorithm for the algebraic Riccati equation.** Implement the iteration
{eq}`eq-16-sign` and use it to solve the regulator algebraic Riccati equation
$A^\top X + X A - X B R^{-1} B^\top X + Q = 0$ for the system

$$
A = \begin{bmatrix} 0 & 1 \\ -2 & -3 \end{bmatrix}, \quad
B = \begin{bmatrix} 0 \\ 1 \end{bmatrix}, \quad
Q = \begin{bmatrix} 1 & 0 \\ 0 & 0 \end{bmatrix}, \quad
R = [\,1\,] .
$$

(a) Form the Hamiltonian $H = \begin{bmatrix} A & -BR^{-1}B^\top \\ -Q & -A^\top \end{bmatrix}$,
iterate {eq}`eq-16-sign` to the matrix sign function, and extract $M_x$ from
{eq}`eq-16-signextract`.

(b) Confirm that $M_x$ equals the solution returned by Chapter 15's `solve_continuous_are`, that
it satisfies the Riccati equation to machine precision, and that the closed loop $A - BR^{-1}B^\top M_x$
is stable. Observe the **quadratic convergence** of the iteration.

```{exercise-end}
```

```{solution-start} ex-16-1
:class: dropdown
```

```{code-cell} ipython3
A = np.array([[0., 1.], [-2., -3.]])
B = np.array([[0.], [1.]])
Q = np.array([[1., 0.], [0., 0.]])
R = np.array([[1.]])
G = B @ np.linalg.inv(R) @ B.T
H = np.block([[A, -G], [-Q, -A.T]])

def matrix_sign(W, tol=1e-13, maxit=50):
    """Roberts (1971) iteration R(G) = 1/2 (G + G^{-1}) -> sign(W)."""
    S = W.astype(float).copy()
    deltas = []
    for _ in range(maxit):
        Snew = 0.5 * (S + np.linalg.inv(S))
        d = np.linalg.norm(Snew - S) / np.linalg.norm(Snew)
        deltas.append(d)
        S = Snew
        if d < tol:
            break
    return S, deltas

Hinf, deltas = matrix_sign(H)
n = A.shape[0]
S11, S12 = Hinf[:n, :n], Hinf[:n, n:]
S21, S22 = Hinf[n:, :n], Hinf[n:, n:]

# Extract M_x as the least-squares solution of (5.6)-(5.8) / eq-16-signextract
lhs = np.vstack([S12, S22 + np.eye(n)])
rhs = -np.vstack([S11 + np.eye(n), S21])
M_x, *_ = np.linalg.lstsq(lhs, rhs, rcond=None)

P = solve_continuous_are(A, B, Q, R)        # Chapter 15 solver
residual = A.T @ M_x + M_x @ A - M_x @ G @ M_x + Q

print("‖M_x - P‖              =", np.linalg.norm(M_x - P))
print("Riccati residual norm  =", np.linalg.norm(residual))
print("closed-loop eigenvalues=", np.linalg.eigvals(A - G @ M_x).round(4))
print("‖ΔS‖ per iteration     =", [f"{d:.1e}" for d in deltas])
```

The sign-algorithm matrix $M_x$ matches `solve_continuous_are` and annihilates the Riccati
residual; the closed loop is stable; and the per-iteration change collapses from $O(10^{-1})$ to
below $10^{-13}$ in a handful of steps. That is the quadratic convergence that makes the matrix
sign algorithm fast.

```{solution-end}
```

```{exercise-start}
:label: ex-16-2
```

**Feedforward as a present value: Heaton's durable-stock rule.** Let the endowment-minus-bliss
gap $g(t) \equiv f(t) - b(t)$ follow the autoregression $D g(t) = -a\, g(t)$ with $a > 0$, so
that $E_t\, g(t+\tau) = e^{-a\tau} g(t)$.

(a) Using the prediction formula of {doc}`12_prediction`, show that the present value in
Heaton's rule {eq}`eq-16-heaton` is $\rho\int_0^\infty e^{-\rho\tau} E_t g(t+\tau)\, d\tau
= \dfrac{\rho}{\rho + a}\, g(t)$, so that net investment is
$D k(t) = \dfrac{a}{\rho + a}\, g(t)$.

(b) Verify the present-value coefficient numerically, and simulate an autoregressive $g$ to
display the implied net-investment path.

```{exercise-end}
```

```{solution-start} ex-16-2
:class: dropdown
```

(a) With $E_t g(t+\tau) = e^{-a\tau} g(t)$, the geometric distributed lead of
{doc}`12_prediction` gives $\int_0^\infty e^{-\rho\tau} e^{-a\tau}\, d\tau = 1/(\rho+a)$, so the
discounted present value is $\rho/(\rho+a)\, g(t)$. Subtracting it from the current gap $g(t)$
leaves $D k(t) = \big[1 - \rho/(\rho+a)\big] g(t) = \dfrac{a}{\rho+a}\, g(t)$. The more
*persistent* the gap (smaller $a$), the more of it is expected to last, the larger the present
value, and the *less* is invested today.

```{code-cell} ipython3
rho, a = 0.05, 0.8

pv_closed = 1.0 / (rho + a)                              # closed form
tau = np.linspace(0, 300, 3_000_001)
pv_num = np.trapz(np.exp(-(rho + a) * tau), tau)         # numerical present value
print(f"present-value coefficient: closed form {pv_closed:.6f}, numerical {pv_num:.6f}")
print(f"net-investment coefficient a/(rho+a) = {a/(rho+a):.6f}")
```

```{code-cell} ipython3
# Simulate an Ornstein-Uhlenbeck (AR(1)) gap g and plot the induced net investment Dk.
rng = np.random.default_rng(1)
dt, T = 0.01, 60.0
N = int(T / dt)
g = np.empty(N); g[0] = 1.0
sig = 0.3
for t in range(N - 1):
    g[t + 1] = g[t] - a * g[t] * dt + sig * np.sqrt(dt) * rng.normal()
Dk = (a / (rho + a)) * g
tt = np.arange(N) * dt

fig, ax = plt.subplots(figsize=(9, 3.5))
ax.plot(tt, g, lw=0.9, label=r'gap $g(t)=f(t)-b(t)$')
ax.plot(tt, Dk, lw=1.3, label=r'net investment $Dk(t)=\frac{a}{\rho+a}\,g(t)$')
ax.axhline(0, color='k', lw=0.5)
ax.set_xlabel('$t$'); ax.legend(); plt.show()
```

Net investment tracks a damped fraction $a/(\rho+a)$ of the current gap: the durable stock is
built up exactly when the present endowment, relative to bliss, exceeds its expected future
present value.

```{solution-end}
```
