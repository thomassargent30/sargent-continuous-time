# A1. Ergodicity and the Consistent Estimation of Second Moments

```{note}
This appendix redeems a promise made in {doc}`/chapters/01_covariance_stationary_processes` and
left unkept by the rest of the book. It is background rather than a chapter: nothing in
Chapters 1–23 depends on it, but every claim in those chapters that some quantity is
*estimable* does.
```

## A1.1 The promise this appendix redeems

{doc}`/chapters/01_covariance_stationary_processes` draws a distinction and then defers it.
Having defined the population moments as averages *across realizations*,

$$
\mu(t) = \int x(t,w)\, dP(w), \qquad
C(t_1,t_2) = \int \big(x(t_1,w)-\mu(t_1)\big)\big(x(t_2,w)-\mu(t_2)\big)\, dP(w),
$$

it observes that one could instead average *across time within a single realization*,
$\lim_{T\to\infty}(2T)^{-1}\int_{-T}^{T} x(t,w)\, dt$, and then says that under "conditions for
ergodicity" the two kinds of average agree, and that such conditions "must be imposed to
acquire a practical theory of estimation."

This matters because *every* estimator in the second part of the book — the sampled
autocovariances $\phi_\tau$ that Phillips inverts in
{doc}`/chapters/21_phillips_continuous_time_estimation`, the matrix covariogram
$\Gamma_0, \Gamma_1$ factored in {doc}`/chapters/20_aggregation_inverse_optimal_predictor`, the
discrete spectral density $S^d(\omega)$ of {doc}`/chapters/17_discrete_sampling_folding`, the
vector autoregressions whose interpretation occupies
{doc}`/chapters/18_time_aggregation_var` — is a *time* average, computed from one realization,
and is being asked to converge to an *ensemble* quantity. Ergodicity is the hypothesis that
licenses the substitution.

The one substantive message is worth stating in advance:

```{important}
**Mean square ergodicity is a restriction on second moments. Covariance ergodicity is a
restriction on *fourth* moments.** Everything the book has built characterizes second moments,
which is exactly enough to settle when the sample *mean* is consistent — and not enough to
settle when the sample *autocovariances* are.
```

**Notation.** Frequency is written $\omega$ here, as in Chapters 19, 22 and 23, rather than the
$w$ of Chapters 8, 10 and 17, so as not to collide with the white noise $w(t)$. Write
$\tilde x(t) = x(t)-\mu$ for the centered process and
$C(\tau) = E\,\tilde x(t)\tilde x(t-\tau)$ for the autocovariance (the book's $R(\tau)$ when
$\mu = 0$). When a spectral density exists it is
$S(\omega) = \int C(\tau)e^{-i\omega\tau}d\tau$; when it does not — the case that matters
below — we use the *spectral distribution* $F$, the nondecreasing function of bounded variation
with

$$
C(\tau) = \frac{1}{2\pi}\int_{-\infty}^{\infty} e^{i\omega\tau}\, dF(\omega),
$$

which every covariance stationary process possesses by Bochner's theorem. $S$ exists precisely
when $F$ is absolutely continuous, and then $dF = S\,d\omega$.

## A1.2 Mean square ergodicity

**Definition A1 (mean square ergodicity).** Let $x(t)$ be covariance stationary with mean $\mu$,
and write the time average over a record of length $2T$ as
$\bar x_T = (2T)^{-1}\int_{-T}^{T}x(t)\,dt$. The process is *mean square ergodic*, or *ergodic
in the mean*, if

$$
\lim_{T\to\infty} E\big(\bar x_T - \mu\big)^2 = 0 .
$$

This is the weakest thing one could ask: a single realization, observed long enough, reveals
the ensemble mean. It is a statement in the same $L^2$ sense used throughout the book — the
mean square convergence of {doc}`/chapters/02_mean_square_continuity_differentiability` — not
almost-sure convergence. The mean square version is the one second-moment theory can settle,
and the one that delivers consistency of estimators.

### The variance of the time average

Everything follows from one computation. Since $E\bar x_T = \mu$,

$$
E(\bar x_T-\mu)^2 = \frac{1}{4T^2}\int_{-T}^{T}\!\!\int_{-T}^{T} C(t-s)\, dt\, ds .
$$

The set of $(t,s)$ in the square with $t-s \in d\tau$ has measure $(2T-|\tau|)\,d\tau$ for
$|\tau|\le 2T$, so the double integral collapses:

```{math}
:label: eq-erg-var
E(\bar x_T - \mu)^2 = \frac{1}{2T}\int_{-2T}^{2T}\Big(1 - \frac{|\tau|}{2T}\Big) C(\tau)\, d\tau .
```

Formula {eq}`eq-erg-var` is the whole of first-order ergodic theory. The weight
$1-|\tau|/(2T)$ is the Bartlett triangle: a record of length $2T$ sees the lag $\tau$ only
through the $2T-|\tau|$ pairs of dates that are $\tau$ apart and both inside the window. The
process is mean square ergodic if and only if this Cesàro average of $C$ vanishes — Slutsky's
condition. It is necessary and sufficient but awkward to check; the frequency domain does
better.

### A necessary and sufficient condition

{doc}`/chapters/10_cramer_representation` derives the criterion by treating $\bar x_T$ as a
filter, and the result is worth recording here in the form used below.

```{admonition} Theorem A1
:class: tip
A covariance stationary process is mean square ergodic **if and only if its spectral
distribution has no atom at frequency zero**, i.e. iff $F(\{0\}) = 0$. Indeed
$\lim_{T\to\infty} E(\bar x_T-\mu)^2 = F(\{0\})/2\pi$.
```

Mean square ergodicity therefore fails in exactly one way: when the process contains a *random
constant* — a component of nonzero variance that never moves — for then no amount of time
averaging can average it away.

### Sufficient conditions

Theorem A1 requires knowing $F$ near the origin. The following are what one checks in practice.

**(S1) $C(\tau)\to0$ as $|\tau|\to\infty$.** If the autocovariance dies out at all, the Cesàro
average in {eq}`eq-erg-var` dies out. This is the workhorse condition. It is *not* necessary —
$C$ may oscillate without decaying and still average to zero — but it covers every model in
this book.

**(S2) $\int_{-\infty}^{\infty}|C(\tau)|d\tau < \infty$.** Absolute integrability gives
ergodicity *and* a rate. Letting $T\to\infty$ in {eq}`eq-erg-var` by dominated convergence,

```{math}
:label: eq-erg-rate
2T\cdot E(\bar x_T-\mu)^2 \;\longrightarrow\; \int_{-\infty}^{\infty} C(\tau)\, d\tau = S(0).
```

**The variance of the sample mean is asymptotically $S(0)/2T$** — the spectral density at zero
frequency divided by the length of the record. The "long-run variance" of time series
econometrics is nothing but the spectrum at the origin.

**(S3) $x$ is purely linearly indeterministic.** This is what the book assumes nearly
everywhere, and it *implies* mean square ergodicity, so no separate hypothesis is needed. If
$x(t) = \int_0^\infty p(\tau)w(t-\tau)d\tau$ with $p \in L^2$ is a Wold representation in the
sense of {doc}`/chapters/08_spectral_densities`, then $S(\omega)=|\tilde P(i\omega)|^2$ exists
as an ordinary function, so $F$ is absolutely continuous and has no atom anywhere. Equivalently,
$C(\tau) = \langle p, p_\tau\rangle$ with $p_\tau(u) = p(u+\tau)$, and $p_\tau \to 0$ weakly in
$L^2$, so $C(\tau)\to0$ and (S1) applies.

### Where the failure lives in the book

The canonical counterexample is a single random variable frozen in time: $x(t)\equiv A$ with
$EA=0$, $EA^2=\sigma^2$. Then $C(\tau)=\sigma^2$ for every $\tau$ and $\bar x_T = A$ for every
$T$. Its spectral distribution is one atom, $dF = 2\pi\sigma^2\delta(\omega)d\omega$, sitting
exactly where Theorem A1 forbids.

This is not exotic; it is the *linearly deterministic* component that Wold's theorem isolates.
{doc}`/chapters/08_spectral_densities` gives that component the spectral density
$S_d(\omega) = \sum_j a_j\pi[\delta(\omega-\omega_j)+\delta(\omega+\omega_j)]$ and the
autocovariance $R^d(\tau)=\sum_j a_j\cos(\omega_j\tau)$. Those $\delta$'s are exactly the atoms
of $F$, so Theorem A1 translates into the book's own decomposition:

> A covariance stationary process is mean square ergodic **iff its linearly deterministic
> component contains no zero-frequency term** — iff no $\omega_j = 0$ in $S_d$.

A deterministic component at a *nonzero* frequency, such as a fixed seasonal cycle, does no
harm: it averages out over a long record. Only the frequency-zero atom survives averaging.
{doc}`/chapters/17_discrete_sampling_folding` shows that *sampling can move an atom onto
frequency zero*, which is where this becomes an aliasing question.

## A1.3 Covariance ergodicity

**Definition A2 (covariance ergodicity).** With $x$ covariance stationary and

$$
\hat C_T(\tau) = \frac{1}{2T}\int_{-T}^{T}\tilde x(t+\tau)\,\tilde x(t)\, dt,
$$

the process is *covariance ergodic*, or *ergodic in the second moments*, if for every fixed
$\tau$,

$$
\lim_{T\to\infty}E\big(\hat C_T(\tau)-C(\tau)\big)^2 = 0 .
$$

(In a finite record one integrates over $|t|\le T-|\tau|$; dividing by $2T$ rather than
$2T-|\tau|$ costs a bias of order $|\tau|/T$ but keeps $\hat C_T$ a positive semidefinite
function of $\tau$ — the property {doc}`/chapters/01_covariance_stationary_processes` identifies
as necessary for $\hat C_T$ to be *some* process's autocovariance, and which spectral estimation
needs.)

### The reduction that makes it tractable

Fix $\tau$ and define the *product process* $y_\tau(t) = \tilde x(t+\tau)\tilde x(t)$. If $x$ is
stationary to fourth order then $y_\tau$ is covariance stationary with mean
$Ey_\tau(t) = C(\tau)$, and $\hat C_T(\tau)$ is exactly its time average. Hence:

> **Covariance ergodicity of $x$ is mean square ergodicity of the product process $y_\tau$, for
> every $\tau$.**

Everything above applies verbatim, one lag at a time. By {eq}`eq-erg-var`,

```{math}
:label: eq-erg-cov
E\big(\hat C_T(\tau)-C(\tau)\big)^2
= \frac{1}{2T}\int_{-2T}^{2T}\Big(1-\frac{|s|}{2T}\Big)\, C_{y_\tau}(s)\, ds,
\qquad
C_{y_\tau}(s) = \operatorname{cov}\big(y_\tau(t+s),\, y_\tau(t)\big).
```

### Why fourth moments are unavoidable

Look at what $C_{y_\tau}$ is:

$$
C_{y_\tau}(s) = E\big[\tilde x(t+s+\tau)\,\tilde x(t+s)\,\tilde x(t+\tau)\,\tilde x(t)\big] - C(\tau)^2 .
$$

It is a *fourth* moment of $x$. This is the essential difference between the two concepts, and
why covariance ergodicity cannot be settled by anything in Chapters 1–17: those chapters
characterize a process entirely through $C(\tau)$ and $S(\omega)$, and two processes with
identical second moments can have wildly different fourth moments, one covariance ergodic and
the other not. Something must be assumed beyond covariance stationarity.

### The Gaussian route

If $x$ is *Gaussian*, its fourth moments are determined by its second moments and the question
closes. By Isserlis's theorem, for jointly Gaussian zero-mean variables
$E[X_1X_2X_3X_4] = E[X_1X_2]E[X_3X_4]+E[X_1X_3]E[X_2X_4]+E[X_1X_4]E[X_2X_3]$. Applying it with
$X_1=\tilde x(t+s+\tau)$, $X_2=\tilde x(t+s)$, $X_3=\tilde x(t+\tau)$, $X_4=\tilde x(t)$, the
first pairing reproduces $C(\tau)^2$ and cancels against the subtracted mean, leaving

```{math}
:label: eq-erg-isserlis
C_{y_\tau}(s) = C(s)^2 + C(s+\tau)\, C(s-\tau).
```

```{admonition} Theorem A2 (Gaussian case)
:class: tip
Let $x$ be a stationary Gaussian process. Then $x$ is covariance ergodic if the Cesàro average
of $C(s)^2 + C(s+\tau)C(s-\tau)$ vanishes for every $\tau$ — in particular if
$C(s)\to 0$ as $|s|\to\infty$.
```

So **for a Gaussian process the same condition $C\to0$ that delivers mean square ergodicity also
delivers covariance ergodicity.** If moreover $C\in L^2$ — equivalently, by Parseval,
$\int S(\omega)^2 d\omega<\infty$ — then {eq}`eq-erg-cov` is $O(1/T)$ and

$$
2T\cdot E\big(\hat C_T(\tau)-C(\tau)\big)^2 \longrightarrow
\int_{-\infty}^{\infty}\big[C(s)^2 + C(s+\tau)C(s-\tau)\big]ds .
$$

Every process generated by the linear stochastic differential equations of
{doc}`/chapters/11_linear_sde` with Gaussian white noise satisfies this:
$p(\tau)=\sum_j g_j e^{\lambda_j\tau}$ with $\operatorname{re}\lambda_j<0$ gives an
exponentially decaying $C$, hence $C \in L^1\cap L^2$.

### The general route: fourth cumulants

Without Gaussianity, decompose the fourth moment into cumulants. Define

$$
\kappa_4(u_1,u_2,u_3) = \operatorname{cum}\big\{\tilde x(t+u_1),\tilde x(t+u_2),\tilde x(t+u_3),\tilde x(t)\big\},
$$

independent of $t$ when $x$ is stationary to fourth order, and identically zero when $x$ is
Gaussian. Then {eq}`eq-erg-isserlis` acquires one extra term:

```{math}
:label: eq-erg-cumulant
C_{y_\tau}(s) = C(s)^2 + C(s+\tau)C(s-\tau) + \kappa_4(s+\tau,\, s,\, \tau).
```

```{admonition} Theorem A3 (general case)
:class: tip
A fourth-order stationary process is covariance ergodic if $C(s)\to0$ as $|s|\to\infty$ and, for
each $\tau$, the Cesàro average of $\kappa_4(s+\tau,s,\tau)$ in $s$ vanishes. A convenient
sufficient condition for the latter is absolute integrability of the fourth cumulant function,
$\iiint|\kappa_4(u_1,u_2,u_3)|\,du_1du_2du_3<\infty$.
```

Absolute integrability of cumulants through fourth order is the standard hypothesis in the
theory of spectral estimation (Brillinger, Hannan). Intuitively it says that fourth-order
dependence, like second-order dependence, dies out with separation in time.

## A1.4 Why these are the ingredients for consistent estimation

### Autocovariances: both concepts are needed

In practice $\mu$ is unknown and replaced by $\bar x_T$. Expanding, the demeaned estimator
equals the known-mean estimator of Definition A2 minus $(\bar x_T-\mu)^2$, up to end effects of
order $|\tau|/T$. Therefore

$$
\underbrace{\text{consistency of } \hat C_T(\tau)}_{\text{what we want}}
\;\Longleftarrow\;
\underbrace{\text{covariance ergodicity}}_{\text{fourth moments}}
\;+\;
\underbrace{\text{mean square ergodicity}}_{\text{second moments}} .
$$

The two do different jobs: covariance ergodicity makes the second-moment time average converge,
mean square ergodicity makes the *centering* harmless. Failure of the latter alone wrecks the
former, and it is worth seeing how. Let $x(t)=\mu+A+\nu(t)$ with $A$ of mean zero and variance
$\sigma^2$ and $\nu$ a well-behaved ergodic noise with autocovariance $C_\nu$. The true
autocovariance is $C(\tau)=\sigma^2+C_\nu(\tau)$. But the sample mean converges to $\mu+A$, not
$\mu$, so the centered data are $x(t)-\bar x_T \to \nu(t)-\bar\nu_T$ and

$$
\hat C_T(\tau) \longrightarrow C_\nu(\tau) = C(\tau)-\sigma^2 .
$$

The estimator converges — to the wrong thing, understating the autocovariance by exactly
$\sigma^2$ at *every* lag. Demeaning removes the random constant along with the mean, and with
it all trace of the non-ergodic component. Nothing in the sample reveals the error.

This is the guarantee that {doc}`/chapters/21_phillips_continuous_time_estimation` needs.
Phillips fits a generating function $M(z)=\sum_{\tau\ge0}\phi_\tau z^\tau$ to the *sampled*
autocovariances and inverts $\lambda_r=\log\mu_r$ to recover the continuous-time poles; every
$\phi_\tau$ entering that computation is a time average from one realization. The same is true
of the matrix covariogram factored in
{doc}`/chapters/20_aggregation_inverse_optimal_predictor` and of the $(B_0,W_0)$ pair that
{doc}`/chapters/22_dimensionality_aliasing_problem` calls estimable.

### Spectral densities: why the periodogram is *not* consistent

Define the periodogram of a record of length $2T$,

```{math}
:label: eq-erg-periodogram
I_T(\omega) = \frac{1}{2T}\left|\int_{-T}^{T}\big(x(t)-\bar x_T\big)e^{-i\omega t}\, dt\right|^2 .
```

It is the sample analogue of the definition of $S(\omega)$ in
{doc}`/chapters/08_spectral_densities`, being the Fourier transform of the sample
autocovariances:
$I_T(\omega)=\int_{-2T}^{2T}(1-|\tau|/2T)\hat C_T(\tau)e^{-i\omega\tau}d\tau$. Two facts about
it are in tension.

*It is asymptotically unbiased*: if $\int|C|<\infty$ and $S$ is continuous at $\omega$, then
$EI_T(\omega)\to S(\omega)$. Covariance ergodicity does that work, through the convergence of
the $\hat C_T(\tau)$ inside.

*It is nonetheless inconsistent*: $\operatorname{Var}I_T(\omega)\to S(\omega)^2$ for
$\omega\ne0$, and $2S(0)^2$ at the origin. The reason is structural. $I_T(\omega)$ is the
squared modulus of a *single* Fourier coefficient, itself asymptotically complex Gaussian with
variance $S(\omega)$ — a fixed number of random quantities, however long the record. So
$I_T(\omega)/S(\omega)$ is asymptotically $\tfrac12\chi^2_2$, with mean 1 and variance 1, and
stays that way forever. Lengthening the record buys more *frequencies*, not more precision at
any one of them. The exercises of {doc}`/chapters/08_spectral_densities` verify this
numerically.

### Smoothing restores consistency

Consistency is recovered by averaging the periodogram over a band of frequencies that shrinks
slowly enough to accumulate an increasing number of them:

$$
\hat S_T(\omega) = \int W_T(\omega-\lambda)I_T(\lambda)\, d\lambda,
\qquad W_T(\cdot)=\frac{1}{b_T}W\!\Big(\frac{\cdot}{b_T}\Big),
$$

or equivalently, in lag-window form,
$\hat S_T(\omega)=\int_{-M_T}^{M_T}k(\tau/M_T)\hat C_T(\tau)e^{-i\omega\tau}d\tau$ with $k$ a
lag window ($k(0)=1$, $k$ even, vanishing outside $[-1,1]$). The two are Fourier duals;
$M_T \sim 1/b_T$. Under the cumulant conditions of Theorem A3,

$$
\operatorname{bias}\hat S_T(\omega) = O(b_T) \text{ or } O(b_T^2),
\qquad
\operatorname{Var}\hat S_T(\omega) \approx \frac{S(\omega)^2}{2Tb_T}\int W(u)^2du,
$$

so choosing the bandwidth to satisfy

$$
b_T \to 0 \quad\text{and}\quad T\,b_T \to \infty
$$

drives both to zero. The two requirements are the familiar bias–variance trade-off: the band
must narrow, so that we estimate $S$ at $\omega$ rather than an average over a wide
neighbourhood, but it must narrow *more slowly* than $1/T$ so that the number of ordinates
averaged, roughly $Tb_T$, grows without bound.

Ergodicity now appears in two distinct places. Covariance ergodicity makes each
$\hat C_T(\tau)$ converge, controlling the *bias*. The fourth-cumulant condition controls the
*variance*, by ensuring that periodogram ordinates at distinct frequencies are asymptotically
uncorrelated, so that averaging $Tb_T$ of them reduces variance by a factor $Tb_T$. Without the
fourth-order condition, smoothing need not help at all.

### Cross-spectra, coherence, and a trap

For a vector process write $C_{jk}(\tau)=E\,\tilde x_j(t+\tau)\tilde x_k(t)$ and
$S_{jk}(\omega)=\int C_{jk}(\tau)e^{-i\omega\tau}d\tau$, whose real part is the co-spectrum and
whose negative imaginary part is the quadrature spectrum. The derived quantities are the squared
coherence and phase,

$$
\mathcal{K}^2_{jk}(\omega)=\frac{|S_{jk}(\omega)|^2}{S_{jj}(\omega)S_{kk}(\omega)},
\qquad \varphi_{jk}(\omega)=\arg S_{jk}(\omega).
$$

Everything above carries over. Cross-covariance ergodicity is again mean square ergodicity of a
product process, now $y_\tau(t)=\tilde x_j(t+\tau)\tilde x_k(t)$, and again a fourth-moment
condition; in the Gaussian case Isserlis gives

$$
\operatorname{cov}\big(y_\tau(t+s),y_\tau(t)\big) = C_{jj}(s)C_{kk}(s) + C_{jk}(s+\tau)C_{kj}(s-\tau),
$$

so $C_{jk}(s)\to0$ for all $j,k$ suffices — the direct analogue of Theorem A2.

There is, however, a trap with no univariate counterpart. Define the cross-periodogram from the
finite Fourier transforms $d_j(\omega)=\int_{-T}^{T}\tilde x_j(t)e^{-i\omega t}dt$ as
$I_{jk}=(2T)^{-1}d_j\overline{d_k}$. Then the *unsmoothed* estimate of squared coherence is

$$
\frac{|I_{jk}|^2}{I_{jj}I_{kk}} = \frac{|d_j|^2|d_k|^2}{|d_j|^2|d_k|^2} = 1
$$

**identically**, at every frequency, for every pair of series, whatever the truth. Two
independent white noises yield an estimated squared coherence of exactly one. Smoothing is not
a refinement here; it is the entire content of the estimator. The moral is the previous one,
sharpened: a single frequency ordinate carries a fixed amount of information no matter how long
the record, so consistency at a point in the frequency domain always comes from pooling
neighbouring frequencies — and it is the fourth-order condition that makes the pooling work.

## A1.5 Discrete records: span, not sampling rate

Every criterion above involves $T\to\infty$, where $2T$ is the *span*. Sampling more finely over
a fixed span does not help: letting $\Delta\to0$ with $T$ fixed drives the sample mean to
$\bar x_T$, a random variable with variance $S(0)/2T>0$, not to the constant $\mu$. No amount
of intra-window resolution removes it.

This is a corrective to a natural misreading of
{doc}`/chapters/23_temporal_aggregation_streamlined`. Marcet's Propositions 6 and 7 show that as
the sampling interval shrinks the discrete forecasts converge to the continuous ones and the
discrete impulse responses recover the continuous kernel. Those are statements about
*population* objects. Fine sampling recovers the continuous-time structure; only a long span
estimates it.

Sampling can also destroy ergodicity outright, by folding a deterministic component onto
frequency zero. That mechanism, and the example that exhibits it, are given in
{doc}`/chapters/17_discrete_sampling_folding`.

## A1.6 Summary

| | mean square ergodicity | covariance ergodicity |
|---|---|---|
| **object** | $\bar x_T \to \mu$ | $\hat C_T(\tau)\to C(\tau)$ for each $\tau$ |
| **order of moments involved** | second | **fourth** |
| **exact criterion** | $F$ has no atom at $\omega=0$ | product process $y_\tau$ has no atom at $\omega=0$ |
| **practical sufficient condition** | $C(\tau)\to0$ | Gaussian and $C(\tau)\to0$; in general, plus a fourth-cumulant condition |
| **automatic for** | any purely linearly indeterministic process | Gaussian linear SDEs of Chapter 11 |
| **rate, when available** | $\operatorname{Var}\bar x_T \sim S(0)/2T$ | $O(1/T)$ when $C\in L^2$ |
| **what it licenses** | centering; the long-run variance | sample autocovariances; spectral and cross-spectral estimation |

Three points are worth carrying away.

1. **The two concepts are not variations on one theme.** Mean square ergodicity is settled by
   the theory this book contains; covariance ergodicity is not, and requires either Gaussianity
   or an explicit condition on fourth cumulants. Covariance stationarity does not give it to
   you.

2. **Ergodicity delivers consistent autocovariances, and that is all.** It does *not* deliver a
   consistent spectral density estimator. The periodogram is asymptotically unbiased and
   permanently noisy; consistency at a frequency comes only from smoothing across neighbouring
   frequencies, at a bandwidth satisfying $b_T\to0$ and $Tb_T\to\infty$. For cross-spectra the
   point is sharper: without smoothing the estimated coherence is identically one.

3. **Consistency is a statement about span.** Sampling a fixed window ever more finely estimates
   nothing consistently, and sampling can even destroy ergodicity by folding a deterministic
   component onto frequency zero.

## References

Bartlett, M. S. (1946). On the theoretical specification and sampling properties of
autocorrelated time series. *Journal of the Royal Statistical Society B*, **8**, 27–41.

Brillinger, D. R. (1981). *Time Series: Data Analysis and Theory*. San Francisco: Holden-Day.

Doob, J. L. (1953). *Stochastic Processes*. New York: Wiley.

Grenander, U., and M. Rosenblatt (1957). *Statistical Analysis of Stationary Time Series*. New
York: Wiley.

Hannan, E. J. (1970). *Multiple Time Series*. New York: Wiley.

Ibragimov, I. A., and Yu. A. Rozanov (1978). *Gaussian Random Processes*. New York:
Springer-Verlag.

Parzen, E. (1962). *Stochastic Processes*. San Francisco: Holden-Day.

Priestley, M. B. (1981). *Spectral Analysis and Time Series*. London: Academic Press.

Rozanov, Y. A. (1967). *Stationary Random Processes*. San Francisco: Holden-Day.

Slutsky, E. (1938). Sur les fonctions aléatoires presque périodiques et sur la décomposition
des fonctions aléatoires stationnaires en composantes. *Actualités Scientifiques et
Industrielles*, **738**, 33–55.
