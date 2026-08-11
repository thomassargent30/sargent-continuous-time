# 0. Introduction: A Reader's Guide

This book is about the theory of systems of linear stochastic differential equations and
some of their applications in economics. It develops the continuous time machinery —
covariance stationary processes, the Wiener process and Itô calculus, spectral densities
and the spectral factorization theorem, the Cramér representation, linear least squares
prediction, and discrete sampling — and then puts it to work. The "work" is the
construction and interpretation of quantitative economic models in the rational
expectations tradition — the research program of Lars Peter Hansen and Thomas Sargent, in
which the decision rules of optimizing agents are tied, through cross-equation restrictions,
to the stochastic laws of motion of the forces that buffet them. This introductory chapter
is a map. It explains what each chapter does, why the chapters are ordered as they are, and
how the separate pieces assemble into a single set of tools for economic modeling.

## Why continuous time?

Economic agents make decisions, and information arrives, more or less continuously; the
quarterly or monthly spacing of the data we observe is an artifact of measurement, not of
the economy. Modeling in continuous time has three payoffs that recur throughout these
notes.

First, it yields an unusually clean *operator calculus*. Once a process is written as the
output of a differential operator driven by noise, $\theta(D)x(t) = \psi(D)w(t)$, questions
about second moments, prediction, and smoothness become questions about rational functions
of the complex variable $s$ and their Laplace transforms. Much of the book is the
development of this calculus.

Second, continuous time draws a sharp line between the *smoothness* of a process and its
*predictability*. A process is mean square differentiable exactly when its moving-average
kernel vanishes at the origin; it is "locally unpredictable" — locally a martingale —
exactly when that kernel does not. This dichotomy, invisible in a purely discrete-time
treatment, turns out to carry real economic content for asset prices, consumption, and
investment.

Third, precisely because the theory lives in continuous time while the data are discrete,
it lets us see *what sampling and time-averaging do to inference*. Aliasing, aggregation
bias, spurious Granger causality, and the non-fundamentalness of innovations are not
nuisances to be assumed away; they are phenomena the continuous-time viewpoint predicts and
quantifies.

## The object at the center: the Wold representation

Almost everything in these notes organizes around a single object. A zero-mean, covariance
stationary, linearly indeterministic process has a fundamental moving-average (Wold)
representation

$$
x(t) = \int_0^\infty p(\tau)\, w(t-\tau)\, d\tau,
$$

where $w(t)$ is a *fundamental* white noise — one whose current and past values span exactly
the same information as current and past $x$. Two readings of this formula run side by side
through the book.

- **The econometrician's reading.** $p(\tau)$ is an impulse response; $w(t)$ is the
  one-step-ahead forecast error, the "innovation" recovered from the data. Summaries of
  $p(\tau)$ (or its discrete-time counterpart) are what vector autoregressions report.

- **The agent's reading.** $w(t)$ is *news*: the surprise that arrives at $t$. A rational
  agent's forecast of any future variable is a known functional of current and past $w$'s,
  and the agent's decisions are convolutions of $p$ with the kernels implied by the
  agent's optimization problem.

The deep claim of the rational expectations program is that these two readings must be
consistent. The continuous-time machinery developed here is what makes the consistency
operational — and what reveals where it can break down.

## How the chapters fit together

The material falls into five movements.

### I. Foundations: stationarity and the mean square calculus (Chapters 1–2)

[Chapter 1](01_covariance_stationary_processes.md) defines a stochastic process on a
probability space, introduces covariance stationarity, and shows that the autocovariance
function $R(\tau)$ — equivalently, a positive semidefinite function — carries all the
second-moment information we will use. [Chapter 2](02_mean_square_continuity_differentiability.md)
develops the *mean square calculus*: a process is mean square continuous or differentiable
according to the smoothness of $R(\tau)$ at the origin, and the autocovariance of a
derivative is obtained by differentiating $R$. The chapter closes with the first appearance
of the smoothness–predictability theme: if $R$ is *analytic*, the process is infinitely
differentiable and perfectly forecastable — a degenerate case we will deliberately avoid,
because interesting economic series are imperfectly predictable.

### II. Building processes from elementary noise (Chapters 3–7)

Here we construct processes from the ground up, out of two canonical sources of randomness,
both obtained by the *same* limiting device.

[Chapter 3](03_poisson_counting_process.md) introduces the Poisson counting process and its
generalized derivative $dN/dt$, a first example of *white noise*: an object with a
$\delta$-function autocovariance that exists only as a generalized process.
[Chapter 4](04_physical_realizability.md) pauses to make precise what it means for such an
object to be "physically realizable," and reassures us that in economic models generalized
processes appear only under integral signs, so that observable variables remain ordinary.
[Chapter 5](05_poisson_driven_processes.md) develops the calculus of stochastic differential
equations driven by Poisson noise — the jump analogue of Itô's rule — and works the random
telegraph wave and the generalized Poisson process.

[Chapter 6](06_wiener_process.md) then obtains the Wiener process (Brownian motion) as the
limit of scaled differences of two independent Poisson counters as the arrival rate goes to
infinity:
jumps shrink, arrivals proliferate, and a Gaussian, continuous-but-nowhere-differentiable
process emerges. [Chapter 7](07_wiener_driven_sde.md) derives Itô's rule and the
moment-recursion formulas for Wiener-driven stochastic differential equations by the same
limiting argument used for Poisson noise — so the reader sees diffusion and jump models as
two faces of one construction.

### III. Representation theory: spectra, Wold, factorization, Cramér (Chapters 8–10)

This is the structural core. [Chapter 8](08_spectral_densities.md) introduces the spectral
density as the Fourier transform of $R(\tau)$, assembles the operational Fourier/Laplace
toolkit (the tables used everywhere downstream), and states the two theorems on which the
rest of the book leans: the **Wold decomposition** and the **spectral factorization
theorem**, $S(w) = \tilde P(iw)\tilde P(-iw)$ with $\tilde P$ free of right-half-plane zeros.
The worked example in which the numerator has a right-half-plane zero is worth dwelling on:
it is the continuous-time prototype of *non-fundamentalness*, the situation in which the
noise spans a larger space than the observables — exactly the identification problem that
later complicates the interpretation of vector autoregressions.

[Chapter 9](09_characterizations_ms_differentiability.md) ties the smoothness theme to the
Wold kernel: $x(t)$ is mean square differentiable iff $p(0)=0$, and the initial-value
theorem $p(0)=\lim_{s\to\infty}s\tilde P(s)$ turns this into a transform computation.
[Chapter 10](10_cramer_representation.md) gives the Cramér representation, exhibiting the
process as an orthogonal superposition of random sinusoids across frequencies, and uses
band-pass filtering to make the frequency decomposition of variance concrete.

### IV. Linear models, prediction, and their economic payoff (Chapters 11–14)

Now the machinery is aimed at economics. [Chapter 11](11_linear_sde.md) studies the
workhorse class $\theta(D)x(t)=\psi(D)w(t)$ of linear constant-coefficient stochastic
differential equations: rational spectral densities, exponential moving-average kernels, and
a differentiability count ($n-1-m$ derivatives) read straight off the orders of the
operators.

[Chapter 12](12_prediction.md) derives the continuous-time Wiener–Kolmogorov prediction
formula and then the single most important economic object in the book: the forecast of a
*geometric distributed lead*,

$$
E_t \int_0^\infty e^{\rho u}\, x(t+u)\, du,
$$

the present value that appears in every asset-pricing equation, permanent-income model, and
quadratic-adjustment-cost Euler equation.

The next two chapters collect the payoff. [Chapter 13](13_locally_unpredictable.md) shows
that such present-value variables are *locally unpredictable* — locally martingales — even
when the underlying $x$ is smooth, because the present-value operator forces the kernel to be
nonzero at the origin. This is the formal counterpart of the observation, due to Sims and to
Hall, that asset prices and consumption should look like martingales at high frequency, and it
is the terminus of the chain begun in Chapter 2: $R''(0) \to p(0) \to n-1-m \to$ local
unpredictability. [Chapter 14](14_nonstationary_examples.md) extends the prediction calculus to
nonstationary processes and shows that Cagan's adaptive expectations scheme is the *optimal*
forecast for a particular process — a small but telling instance of the rational-expectations
principle that expectations formulas are not free parameters but implications of the law of
motion. [Chapter 20](20_aggregation_inverse_optimal_predictor.md) later develops this very
example into a full account of what aggregation over time does to Cagan's scheme.

### V. State space, filtering, and the computation of equilibria (Chapters 15–16)

The prediction theory of Chapters 11–12 was built in the *frequency domain*, out of Laplace
transforms of moving-average kernels. The same theory has an equivalent *time-domain* form,
and that form is what one computes with.

[Chapter 15](15_kalman_filter_spectral_factorization.md) recasts prediction in
*state-space* terms. For processes with a rational spectral density, the Wold representation
and the spectral factorization theorem of Chapter 8 acquire a finite-dimensional counterpart:
a linear state-space model $(A,B,C)$ together with the **Kalman–Bucy filter**,
whose steady state solves an algebraic Riccati equation and delivers the **innovations
representation** — the state-space realization of the Wold representation, with the Kalman
gain in place of the moving-average kernel. Solving the Riccati equation is the time-domain
algorithm behind spectral factorization. The same chapter introduces the **continuous-time
linear regulator** and its Bellman equation, and the **duality** that makes optimal filtering
and optimal control one Riccati equation read two ways — the estimation-and-control pairing
that recurs throughout the Hansen–Sargent program. The innovation it constructs is also the
object whose recoverability from sampled data is at issue throughout Part II.

[Chapter 16](16_faster_methods_recursive_linear_models.md) (Hansen, Heaton, and Sargent) turns
the linear regulator into a tool for *computing the equilibria* of dynamic linear-quadratic
economies: a fictitious social planning problem whose solution is an optimal linear regulator,
with Arrow–Debreu prices read off the gradient of the value function. It extends the Chapter-15
machinery to nonautonomous forcing — so the decision rule *separates* into a feedback (control)
part and a feedforward (prediction) part — and supplies fast Riccati solvers, the matrix sign
algorithm, that scale to the large state vectors realistic economies demand.

### VI. From continuous time to data: sampling, aggregation, and identification (Chapters 17–23)

The final movement confronts the gap between a continuous-time theory and discretely sampled,
often time-averaged, data. [Chapter 17](17_discrete_sampling_folding.md) derives the
*folding formula*, which shows how sampling aliases high continuous-time frequencies into the
discrete spectrum — the precise sense in which a sampling interval throws information away, and
the engine of the aliasing problem taken up in Chapters 21 and 22.

The companion papers that close the book apply the whole apparatus to the central problem of
inference. In the versions presented here they are retold in the third person, as the present
account of how each paper *uses* — and extends — the machinery built up above.

They fall into two groups. Chapters 18 and 20 ask what sampling does to the *interpretation* of
discrete-time evidence; Chapters 21 and 22 ask whether the continuous-time primitive can be
*recovered* from that evidence at all. Chapter 19 supplies machinery both groups use, and
Chapter 23 generalizes the first group.

[Chapter 18](18_time_aggregation_var.md) (Hansen and Sargent) asks when the discrete-time
moving-average coefficients of a sampled process resemble the continuous-time kernel $p(\tau)$,
and exhibits the systematic ways in which vector autoregressions can mislead.
[Chapter 19](19_prediction_formulas_continuous_time.md) (Hansen and Sargent) supplies the
general Laplace-transform prediction machinery — the analytic decomposition and annihilation
operator — that solves continuous-time rational expectations models, generalizing the
geometric-lead formula of Chapter 12 and computing the feedforward terms of Chapter 16.
[Chapter 20](20_aggregation_inverse_optimal_predictor.md) (Hansen and Sargent) is the natural
sequel to Chapter 14: it constructs the continuous-time money-creation/inflation process for
which Cagan's adaptive expectations are *exactly* optimal — using the prediction calculus of
Chapter 12 — and then samples it, showing how aggregation over time manufactures Granger
causality from money to inflation that is wholly absent in continuous time, and gauging the
bias in Cagan's $\lambda \approx e^{-\beta}$ approximation.

The next two chapters turn to *identification* — recovering the continuous-time
primitive from its discrete shadow. [Chapter 21](21_phillips_continuous_time_estimation.md)
retells A. W. Phillips's landmark 1959 paper, framed by an appreciation by Hansen and Sargent,
as a synthesis of nearly everything the book has built: rational transfer functions and the
spectral factorization theorem (Chapter 8), the sum-of-exponentials autocovariance of a
rational process (Chapter 11), white noise as a physically unrealizable but admissible input
(Chapter 4), and the sampling map $\mu = e^{\lambda}$ that turns the folding formula of
Chapter 17 into the aliasing problem.
[Chapter 22](22_dimensionality_aliasing_problem.md) (Hansen and Sargent) closes the loop on
aliasing: it *counts* the continuous-time models consistent with given discrete data and shows
that the restriction to a rational spectral density (Chapters 8 and 11) collapses the
uncountable ambiguity of the general covariance-stationary case — and even the countable
ambiguity that P. C. B. Phillips (1973) had identified in the first-order Markov case — down to
a *finite* set, often a single model, with fine enough sampling resolving it entirely.

The book closes with [Chapter 23](23_temporal_aggregation_streamlined.md), an account of Albert
Marcet's work, which recasts the comparison of continuous and discrete Wold representations as
an exercise in the $L^2$ projection of Chapter 19: the sampled moving-average kernel is
*literally* the projection of the continuous kernel onto its own integer shifts. Marcet's
central point is that *discontinuities* in the one-sided continuous kernel — the very feature
that $p(0)\neq 0$ creates, and that Chapters 9 and 13 tie to local unpredictability — open up
possibilities for the sampled kernel (inflated leading coefficients, within-interval structure
that becomes invisible, manufactured or destroyed Granger causality) that the smooth,
two-sided distributed lags of Sims and Geweke never had.

## Threads to follow

Several ideas recur and are worth tracking deliberately across chapters.

- **Fundamentalness (invertibility).** Introduced through the spectral factorization theorem
  in [Chapter 8](08_spectral_densities.md), the distinction between fundamental and
  non-fundamental white noise is the continuous-time root of the difficulties in interpreting
  vector autoregressions taken up in [Chapter 18](18_time_aggregation_var.md) and
  [Chapter 23](23_temporal_aggregation_streamlined.md). When agents see news the
  econometrician cannot recover from the observables, the two readings of the Wold
  representation diverge. Two distinct conditions do the work here and are worth keeping
  apart: *positive semidefiniteness* of the spectral density is what permits a factorization
  at all, while the absence of *right-half-plane zeros* is what picks out the fundamental
  factor from among the many. It is the first of these that returns in
  [Chapter 22](22_dimensionality_aliasing_problem.md) — imposed there on the *implied
  continuous-time* spectral density, or equivalently on the intensity matrix $V_k$ — as the
  test deciding which aliased continuous-time models are admissible at all.

- **Smoothness versus predictability.** The chain $R''(0)$ ([Chapter 2](02_mean_square_continuity_differentiability.md))
  $\to$ $p(0)$ ([Chapter 9](09_characterizations_ms_differentiability.md)) $\to$ local
  unpredictability ([Chapter 13](13_locally_unpredictable.md)) is one continuous argument. It
  is what lets us say which economic variables should behave like martingales and which should
  be smooth — and, in [Chapter 23](23_temporal_aggregation_streamlined.md), the same
  kernel-at-the-origin condition governs how badly time aggregation distorts a model.

- **Prediction as an operator.** The forecast is a projection, realized analytically by the
  annihilation operator $[\,\cdot\,]_+$ and the Laplace-transform decomposition. It appears
  first in [Chapter 12](12_prediction.md), takes state-space form via the Kalman filter in
  [Chapter 15](15_kalman_filter_spectral_factorization.md), drives the feedforward decision
  rules of [Chapter 16](16_faster_methods_recursive_linear_models.md), is generalized in
  [Chapter 19](19_prediction_formulas_continuous_time.md), reappears as $L^2$ projection onto
  the span of sampled kernels in [Chapter 23](23_temporal_aggregation_streamlined.md), and
  underlies the optimal-forecast constructions of [Chapters 14](14_nonstationary_examples.md)
  and [20](20_aggregation_inverse_optimal_predictor.md).

- **The continuous–discrete bridge.** [Chapters 17](17_discrete_sampling_folding.md),
  [18](18_time_aggregation_var.md), and [23](23_temporal_aggregation_streamlined.md) are a
  sustained argument that the continuous-time model is the primitive object and the
  discrete-time model a derived, lossy summary of it — and that recovering the primitive from
  the summary requires exactly the cross-equation restrictions a structural rational
  expectations model supplies.

- **Aliasing and identification.** A second strand of that bridge asks whether the lost
  information can be recovered at all. The folding formula ([Chapter 17](17_discrete_sampling_folding.md))
  and Phillips's multivalued $\lambda = \log\mu$ ([Chapter 21](21_phillips_continuous_time_estimation.md))
  pose the aliasing problem; [Chapter 22](22_dimensionality_aliasing_problem.md) counts its
  solutions and shows that the rational-spectral-density restriction makes them finite; and
  [Chapter 20](20_aggregation_inverse_optimal_predictor.md) exhibits the complementary route by
  which a structural rational-expectations model's cross-equation restrictions pin the
  continuous-time model down.

## How to read these notes

The chapters are cumulative, and a first reading is best done in order. Readers who want the
shortest path to the economics can, after [Chapters 1–2](01_covariance_stationary_processes.md),
proceed directly to the representation theory of [Chapter 8](08_spectral_densities.md) and then
to prediction in [Chapters 11–14](11_linear_sde.md), treating the Poisson and Wiener
constructions of [Chapters 3–7](03_poisson_counting_process.md) as the source of the white
noise that drives everything else. Readers primarily interested in econometric inference can
read [Chapter 8](08_spectral_densities.md), then [Chapter 17](17_discrete_sampling_folding.md)
and the companion papers, [Chapters 18–23](18_time_aggregation_var.md) — adding the prediction
calculus of [Chapters 12–14](12_prediction.md) before the identification chapters
[21](21_phillips_continuous_time_estimation.md) and
[22](22_dimensionality_aliasing_problem.md). Readers whose interest is computational can go
from [Chapter 12](12_prediction.md) directly to
[Chapters 15–16](15_kalman_filter_spectral_factorization.md) and thence to
[Chapter 19](19_prediction_formulas_continuous_time.md), which supplies the prediction formulas
those chapters need. The companion papers are largely
self-contained but presuppose the Wold representation, spectral factorization, and the
prediction calculus of the first part.

The prerequisites are modest: probability at the level of expectations and covariances, the
elements of Fourier and Laplace transforms (collected in the tables of
[Chapter 8](08_spectral_densities.md)), and a willingness to manipulate the derivative
operator $D$ as an algebraic symbol. Throughout, $E$ denotes expectation, $\mathbb{R}$ the
real line, $D$ the time-derivative operator, and $\delta(\cdot)$ the Dirac delta generalized
function; a process is written $x_t = x(t,w)$ with $w$ indexing the underlying probability
space.
