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
13. Locally Unpredictable Stochastic Processes
14. Examples of Nonstationary Processes
15. Discrete Sampling: The Folding Formula

## Notation

Throughout, $E$ denotes the mathematical expectation operator and $\mathbb{R}$ the real
line. A continuous time stochastic process is written $x_t = x(t,w)$, where $w$ indexes
the underlying probability space. $D$ denotes the time-derivative operator and $\delta(\cdot)$
the Dirac delta generalized function.

```{note}
A number of figures in the original manuscript were hand-drawn and are not yet
reproduced here; their locations are marked in the text with a *(figure)* placeholder.
```
