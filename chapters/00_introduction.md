# 0. Introduction: A Reader's Guide

These notes develop the theory of continuous time stochastic processes and then put it to
work. The "work" is the construction and interpretation of quantitative economic models in
the rational expectations tradition — the research program of Lars Peter Hansen and Thomas
Sargent, in which the decision rules of optimizing agents are tied, through cross-equation
restrictions, to the stochastic laws of motion of the forces that buffet them. This
introductory chapter is a map. It explains what each chapter does, why the chapters are
ordered as they are, and how the separate pieces assemble into a single set of tools for
economic modeling.

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
limit of scaled, compensated Poisson differences as the arrival rate goes to infinity:
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

### IV. Linear models, prediction, and economics (Chapters 11–14)

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
quadratic-adjustment-cost Euler equation. [Chapter 13](13_locally_unpredictable.md) shows
that such present-value variables are *locally unpredictable* — locally martingales — even
when the underlying $x$ is smooth, because the present-value operator forces the kernel to be
nonzero at the origin. This is the formal counterpart of the observation, due to Sims and to
Hall, that asset prices and consumption should look like martingales at high frequency.
[Chapter 14](14_nonstationary_examples.md) extends the prediction calculus to nonstationary
processes and shows that Cagan's adaptive expectations scheme is the *optimal* forecast for a
particular process — a small but telling instance of the rational-expectations principle that
expectations formulas are not free parameters but implications of the law of motion.

### V. From continuous time to data: sampling and aggregation (Chapters 15–19)

The final movement confronts the gap between a continuous-time theory and discretely sampled,
often time-averaged, data. [Chapter 15](15_discrete_sampling_folding.md) derives the
*folding formula*, which shows how sampling aliases high continuous-time frequencies into the
discrete spectrum — the precise sense in which a sampling interval throws information away.

The four companion papers that close the book apply the whole apparatus to the central
problem of inference. [Chapter 16](16_time_aggregation_var.md) (Hansen and Sargent) asks when
the discrete-time moving-average coefficients of a sampled process resemble the continuous-time
kernel $p(\tau)$, and exhibits the systematic ways in which vector autoregressions can mislead.
[Chapter 17](17_prediction_formulas_continuous_time.md) (Hansen and Sargent) provides the
general Laplace-transform prediction machinery — the analytic decomposition and annihilation
operator — that solves continuous-time rational expectations models, generalizing the
geometric-lead formula of Chapter 12. [Chapter 18](18_temporal_aggregation_streamlined.md)
(Marcet) develops, in the language of $L^2$ projection introduced in Chapter 17, the exact
relationship between the continuous and discrete Wold representations: how sampling and
unit-averaging contaminate impulse responses, distort Granger-causality, and bias the
apparent importance of variables, and in what sense finer sampling recovers the continuous
model. Finally, [Chapter 19](19_phillips_continuous_time_estimation.md) presents A. W. Phillips's
landmark 1959 paper on inferring a continuous-time system from discrete-time data — read as an
application of the book's machinery and framed by an appreciation by Hansen and Sargent — and
through it surveys the whole literature on aggregation over time, placing the preceding chapters
within the program of identifying continuous-time rational expectations models from discrete
data.

## Threads to follow

Several ideas recur and are worth tracking deliberately across chapters.

- **Fundamentalness (invertibility).** Introduced through the spectral factorization theorem
  in [Chapter 8](08_spectral_densities.md), the distinction between fundamental and
  non-fundamental white noise is the continuous-time root of the difficulties in interpreting
  vector autoregressions taken up in [Chapter 16](16_time_aggregation_var.md) and
  [Chapter 18](18_temporal_aggregation_streamlined.md). When agents see news the
  econometrician cannot recover from the observables, the two readings of the Wold
  representation diverge.

- **Smoothness versus predictability.** The chain $R''(0)$ ([Chapter 2](02_mean_square_continuity_differentiability.md))
  $\to$ $p(0)$ ([Chapter 9](09_characterizations_ms_differentiability.md)) $\to$ local
  unpredictability ([Chapter 13](13_locally_unpredictable.md)) is one continuous argument. It
  is what lets us say which economic variables should behave like martingales and which should
  be smooth.

- **Prediction as an operator.** The forecast is a projection, realized analytically by the
  annihilation operator $[\,\cdot\,]_+$ and the Laplace-transform decomposition. It appears
  first in [Chapter 12](12_prediction.md), is generalized in
  [Chapter 17](17_prediction_formulas_continuous_time.md), and reappears as $L^2$ projection
  onto the span of sampled kernels in [Chapter 18](18_temporal_aggregation_streamlined.md).

- **The continuous–discrete bridge.** [Chapters 15](15_discrete_sampling_folding.md),
  [16](16_time_aggregation_var.md), and [18](18_temporal_aggregation_streamlined.md) are a
  sustained argument that the continuous-time model is the primitive object and the
  discrete-time model a derived, lossy summary of it — and that recovering the primitive from
  the summary requires exactly the cross-equation restrictions a structural rational
  expectations model supplies.

## How to read these notes

The chapters are cumulative, and a first reading is best done in order. Readers who want the
shortest path to the economics can, after [Chapters 1–2](01_covariance_stationary_processes.md),
proceed directly to the representation theory of [Chapter 8](08_spectral_densities.md) and then
to prediction in [Chapters 11–13](11_linear_sde.md), treating the Poisson and Wiener
constructions of [Chapters 3–7](03_poisson_counting_process.md) as the source of the white
noise that drives everything else. Readers primarily interested in econometric inference can
read [Chapter 8](08_spectral_densities.md), then [Chapter 15](15_discrete_sampling_folding.md)
and the companion papers, [Chapters 16–18](16_time_aggregation_var.md). The companion papers
are largely self-contained but presuppose the Wold representation, spectral factorization, and
the prediction calculus of the first part.

The prerequisites are modest: probability at the level of expectations and covariances, the
elements of Fourier and Laplace transforms (collected in the tables of
[Chapter 8](08_spectral_densities.md)), and a willingness to manipulate the derivative
operator $D$ as an algebraic symbol. Throughout, $E$ denotes expectation, $\mathbb{R}$ the
real line, $D$ the time-derivative operator, and $\delta(\cdot)$ the Dirac delta generalized
function; a process is written $x_t = x(t,w)$ with $w$ indexing the underlying probability
space.
