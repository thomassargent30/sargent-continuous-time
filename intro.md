#  Continuous Time Stochastic Processes

## Thomas J Sargent

This Jupyter book develops parts of the theory of **continuous time stochastic processes** used in modern macroeconomics and quantitative economics.
It is based on lecture notes that I prepared for a 2nd year PhD class at the University of Minnesota during the 1983-1984 academic year.  As we often did in the economics department at Minnesota in those days, I taught the course not because I knew very much about the subject, but because I wanted to learn about it.

Starting from covariance stationarity and mean square calculus, my 1983 notes built up the Poisson counting process, stochastic differential equations driven by Poisson and Wiener processes (Itô's rule), spectral densities, the Wold decomposition and spectral
factorization theorem, the Cramér representation, linear stochastic differential
equations, linear least squares prediction, locally unpredictable processes, and discrete sampling (the folding formula).

Chapters 1 through 14 and Chapter 17 are revised and improved versions of my 1983 class notes.
Chapters 15 and 16 are new. Chapter 15 recasts the prediction theory of Chapters 11 and 12 in
state space form and derives the Kalman-Bucy filter and the continuous time linear regulator.
Chapter 16 reports work of Hansen, Heaton, and Sargent on computing equilibria of linear
quadratic economies; it sits in Part I because it develops the linear regulator of Chapter 15.

```{note}
The figures in my 1983 notes  were hand-drawn. They have been redrawn here; the
scripts that generate them live alongside the images in `chapters/figures/`.
```

Part II applies the Part I tools to time aggregation, prediction, and the relation between
continuous and discrete time models. Chapters 18 and 19 reprint work of Hansen and Sargent and
carry their byline. Chapters 20 through 23 report and paraphrase four further papers rather than
reprinting them. To all six I add exercises, answers, and Python code.

## Why continuous time

Economic agents decide, and information arrives, more or less continuously. The quarterly or
monthly spacing of our data is an artifact of measurement, not of the economy.

Continuous time buys three things.

First, an operator calculus. Write a process as the output of a differential operator driven by
a white noise, $\theta(D) x(t) = \psi(D) w(t)$. Questions about second moments, prediction, and
smoothness become questions about rational functions of a complex variable $s$. Much of Part I
develops this calculus.

Second, a sharp line between the smoothness of a process and its predictability. The process
$x(t)$ is mean square differentiable when its moving average kernel satisfies $p(0) = 0$; it is
locally unpredictable when $p(0) \neq 0$. A discrete time treatment hides this distinction.
Chapter 13 uses it to say which economic variables should behave like martingales.

Third, an account of what sampling and time averaging do to inference. Aliasing, aggregation
bias, spurious Granger causality, and non-fundamental innovations are consequences that a
continuous time model predicts and quantifies. Part II computes them.

## The Wold representation

Most of what follows organizes around one object. A zero mean, covariance stationary, linearly
indeterministic process $x(t)$ has a fundamental moving average representation

$$
x(t) = \int_0^\infty p(\tau)\, w(t-\tau)\, d\tau,
$$

in which $w(t)$ is fundamental for $x(t)$: current and past $w$'s span the space that current
and past $x$'s span.

An econometrician reads $p(\tau)$ as an impulse response and $w(t)$ as a one step ahead
forecast error. Summaries of $p(\tau)$ are what vector autoregressions report.

An agent inside a model reads $w(t)$ as news. That agent's forecast of a future variable is a
known function of current and past $w$'s, and that agent's decision rules are convolutions of
$p$ with kernels that the agent's optimum problem implies.

### Part I. Continuous Time Stochastic Processes

1. Covariance Stationary Stochastic Processes
2. Mean Square Continuity and Differentiability of a Stochastic Process
3. The Poisson Counting Process
4. The Concept of Physical Realizability
5. Stochastic Processes Driven by a Poisson Counting Process
6. The Wiener Process
7. Stochastic Differential Equations Driven by a Wiener Process
8. Spectral Densities (General Results; Wold's Theorem; the Spectral Factorization Theorem)
9. Characterizations of Mean Square Differentiability and Mean Square Continuity
10. The Cramér Representation
11. Linear Stochastic Differential Equations
12. Linear Least Squares Prediction (Wiener–Kolmogorov; Geometric Distributed Leads)
13. Locally Unpredictable Stochastic Processes
14. Examples of Nonstationary Processes
15. State-Space Models, the Kalman Filter, and Spectral Factorization
16. Faster Methods for Solving Recursive Linear Models of Dynamic Economies (Hansen, Heaton, and Sargent)
17. Discrete Sampling: The Folding Formula

Chapters 1 and 2 define covariance stationarity and the mean square calculus. Chapters 3
through 7 construct white noise twice, once from a Poisson counter and once from a Wiener
process, by one limiting argument. Chapters 8 through 10 state the Wold decomposition, the
spectral factorization theorem, and the Cramér representation. Chapters 11 and 12 solve
$\theta(D)x(t) = \psi(D)w(t)$ and derive the prediction formulas. Chapters 13 and 14 collect
two consequences: which variables behave like martingales, and which process makes Cagan's
adaptive expectations rational. Chapters 15 and 16 rebuild the prediction theory in state space
form and compute with it. Chapter 17 derives the folding formula and opens the way to Part II.

### Part II. Companion Papers: Aggregation, Prediction, Identification

18. Time Aggregation (reprinted from Hansen and Sargent, *Two Difficulties in Interpreting Vector Autoregressions*, §2)
19. Prediction Formulas for Continuous Time Linear Rational Expectations Models (reprinted from a note by Hansen and Sargent)
20. Aggregation Over Time and the Inverse Optimal Predictor Problem for Adaptive Expectations in Continuous Time (Hansen and Sargent, *IER*, 1983)
21. Inferring a Continuous-Time System from Discrete-Time Data: An Appreciation of A. W. Phillips (1959) (presenting Phillips, *Biometrika*, 1959; the framing draws on an appreciation of Phillips by Hansen and Sargent)
22. The Dimensionality of the Aliasing Problem in Models with Rational Spectral Densities (Hansen and Sargent, *Econometrica*, 1983)
23. Temporal Aggregation of Economic Time Series (Marcet, in *Rational Expectations Econometrics*, eds. Hansen and Sargent, Westview, 1991)

Chapters 18 and 20 ask what aggregation over time does to vector autoregressions and to Granger
causality. Chapters 21 and 22 ask whether discretely sampled data can recover a continuous time
model at all. Chapter 19 supplies prediction formulas that both groups use. Chapter 23 states
the general theory of what sampling does to a moving average kernel.

### Appendices

Background material that the chapters in Parts I and II presuppose but do not develop.

A1. Ergodicity and the Consistent Estimation of Second Moments. Conditions under which time
averages computed from a single realization converge to the ensemble moments that this book
characterizes, and what those conditions do and do not buy for estimating autocovariances,
spectra, and cross spectra.

## Threads

Four ideas recur.

**Fundamentalness.** Chapter 8 states the spectral factorization theorem and distinguishes
fundamental from non-fundamental white noise. Two conditions do separate work. Positive
semidefiniteness of the spectral density permits a factorization; absence of right half plane
zeros picks the fundamental factor out of the many. When agents see news that an econometrician
cannot recover from observables, the econometrician's reading of the Wold representation and
the agent's reading diverge. Chapters 18 and 23 measure the divergence. Chapter 22 imposes
positive semidefiniteness on an implied continuous time spectral density and uses it to decide
which aliased models are admissible.

**Smoothness versus predictability.** $R''(0)$ in Chapter 2 becomes $p(0)$ in Chapter 9,
becomes the count $n-1-m$ in Chapter 11, and becomes local unpredictability in Chapter 13. That
is one argument carried across four chapters. The same condition on the kernel at the origin
governs how badly time aggregation distorts a model in Chapter 23.

**Prediction as an operator.** A forecast is a projection. The annihilation operator
$[\,\cdot\,]_+$ and a Laplace transform decomposition compute it. Chapter 12 introduces it,
Chapter 15 gives it state space form, Chapter 16 uses it for feedforward decision rules,
Chapter 19 generalizes it, and Chapter 23 recovers it as an $L^2$ projection onto integer
shifts of a kernel.

**The continuous-discrete bridge.** Chapters 17, 18, and 23 treat a continuous time model as
the primitive object and a discrete time model as a lossy summary of it. Chapters 21 and 22 ask
what recovers the primitive from the summary. Cross equation restrictions of a structural
rational expectations model are one answer; fine enough sampling is another.

## How to read this book

A first reading is best done in order. Three shorter paths work.

For the economics, read Chapters 1 and 2, then Chapter 8, then Chapters 11 through 14, taking
the white noise of Chapters 3 through 7 as given.

For econometric inference, read Chapter 8, then Chapter 17, then Chapters 18 through 23, adding
Chapters 12 through 14 before Chapters 21 and 22.

For computation, read Chapter 12, then Chapters 15 and 16, then Chapter 19.

Prerequisites are modest: probability at the level of expectations and covariances, elements of
Fourier and Laplace transforms, and a willingness to manipulate the derivative operator $D$ as
an algebraic symbol. Chapter 8 tabulates the transforms that later chapters use.

## Notation

Throughout, $E$ denotes the mathematical expectation operator and $\mathbb{R}$ the real
line. A continuous time stochastic process is written $x_t = x(t,w)$, where $w$ indexes
the underlying probability space. $D$ denotes the time-derivative operator and $\delta(\cdot)$
the Dirac delta generalized function.

The letter $w$ does two jobs. Written with an argument, $w(t)$ is a white noise. Written alone,
$w$ is a point of the underlying probability space, as in $x(t,w)$. Frequency is always
$\omega$, and a second frequency, when one is needed, is $\nu$. Spectral densities are
$S(\omega)$, and $\lambda$ is reserved for the roots of a differential operator and for the
arrival rate of a Poisson process.
