# Preliminary Introduction to Continuous Time Stochastic Processes

These notes develop the theory of **continuous time stochastic processes** used in
modern macroeconomics and quantitative economics:

> Thomas J. Sargent, *Preliminary Introduction to Continuous Time Stochastic Processes*.

Starting from covariance stationarity and mean square calculus, the notes build up the
Poisson counting process, stochastic differential equations driven by Poisson and Wiener
processes (Itô's rule), spectral densities, the Wold decomposition and spectral
factorization theorem, the Cramér representation, linear stochastic differential
equations, linear least squares prediction, locally unpredictable processes, and discrete
sampling (the folding formula).

## About this edition

Relative to the original typed manuscript, this MyST/Jupyter Book edition:

- **Corrects typographical and mathematical errors** present in the LaTeX source.
- **Splits the chapter into one Markdown file per section** for easier reading and
  cross-referencing online.
- Renders all mathematics with MathJax so equations can be searched, linked, and copied.

## Contents

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
13. State-Space Models, the Kalman Filter, and Spectral Factorization
14. Faster Methods for Solving Recursive Linear Models of Dynamic Economies (Hansen, Heaton, and Sargent)
15. Locally Unpredictable Stochastic Processes
16. Examples of Nonstationary Processes
17. Discrete Sampling: The Folding Formula

## Companion Papers

The second part collects companion papers that apply the continuous time machinery
developed above to problems of time aggregation, prediction, and the relationship between
continuous and discrete time models. They originally appeared as separate chapters of a related
volume and are reproduced here for convenience.

18. Time Aggregation (Hansen and Sargent, from *Two Difficulties in Interpreting Vector Autoregressions*, §2)
19. Prediction Formulas for Continuous Time Linear Rational Expectations Models (Hansen and Sargent)
20. Inferring a Continuous-Time System from Discrete-Time Data: An Appreciation of A. W. Phillips (1959) (presenting Phillips, *Biometrika*, 1959, with Hansen and Sargent)
21. Aggregation Over Time and the Inverse Optimal Predictor Problem for Adaptive Expectations in Continuous Time (Hansen and Sargent, *IER*, 1983)
22. The Dimensionality of the Aliasing Problem in Models with Rational Spectral Densities (Hansen and Sargent, *Econometrica*, 1983)
23. Temporal Aggregation of Economic Time Series (Marcet)

## Notation

Throughout, $E$ denotes the mathematical expectation operator and $\mathbb{R}$ the real
line. A continuous time stochastic process is written $x_t = x(t,w)$, where $w$ indexes
the underlying probability space. $D$ denotes the time-derivative operator and $\delta(\cdot)$
the Dirac delta generalized function.

```{note}
A number of figures in the original manuscript were hand-drawn and are not yet
reproduced here; their locations are marked in the text with a *(figure)* placeholder.
```
