# Resolved missing references (`(—)`) — chapters 01–16

For each `(—)` placeholder this document proposes the specific referent. Line numbers
are the **current** file line numbers (after the frontmatter / figures / exercises added
in earlier sessions). Chapters **01, 04, 11, and 16 contain no placeholders.**

Conventions used below:
- When the target equation already carries a MyST label, I cite it (e.g. `{eq}`eq-3-1``).
- Most displayed equations are **unlabeled** in the source; for these I identify the
  target by its content and give the line(s) where it sits, so you can either add a
  `:label:` and reference it or wire in a number by hand.
- Items that point to a **theorem / section / chapter / citation** rather than an
  equation are flagged as such.

---

## Chapter 02 — `02_mean_square_continuity_differentiability.md`

This chapter's proof of the Taylor-series theorem (L310–346) uses several ad-hoc inline
tags — `(*)`, `(+)`, `(a)`, `(b)` — *and* two `(—)` placeholders. Resolving the `(—)`
requires pinning all of them, so they are listed together.

- **L310 — "Theorem (—)."** *(theorem label)* → This is the next theorem after **Theorem 8**
  (L271). Number it **Theorem 9**.
- **L344 — "Substituting the right side of (—) into (+)"** → `(—)` = the series for
  $\hat x(t+\tau)$, i.e. $\hat x(t+\tau)=\sum_{n=0}^\infty x^{(n)}(t)\tfrac{\tau^n}{n!}$ at
  **L328–330**. *(Here `(+)` = the squared-error identity at L319–324, whose "two terms in
  braces" are referenced at L346.)*
- **L345 — "the reasoning that led to theorem (—)"** *(theorem ref)* → the autocorrelation-of-derivatives
  result, **Theorem 8** (L271–279) together with **Corollary 7C** (L259–266); this is what
  gives $E\,x^{(n)}(t)x^{(m)}(t-\tau)=(-1)^m R^{(n+m)}(\tau)$.

*Related dangling tags worth fixing in the same pass (not `(—)`):* `(*)` at L332 = the
analytic Taylor representation $R(\tau)=\sum R^{(n)}(0)\tau^n/n!$ (L303–305); `(a)` (cited
L345) = the differentiated series $R^{(m)}(\tau)=\sum_{n=m}^\infty R^{(n)}(0)\tfrac{\tau^{n-m}}{(n-m)!}$
(L334–336); `(b)` is tagged at L341; `(+)` = the identity at L319–324.

---

## Chapter 03 — `03_poisson_counting_process.md`

- **L169 — "It follows from (—) that $C(t_1,t_2)=\lambda\min(t_1,t_2)$"** → `(—)` = the
  autocorrelation $R(t_1,t_2)=\lambda\min(t_1,t_2)+\lambda^2 t_1 t_2$ at **L164–166** (with
  $E\,N(t)=\lambda t$, so $C=R-\lambda^2 t_1 t_2$).
- **L175 — "From (—) it follows that $N(t)$ is mean square continuous"** → `(—)` = the
  autocovariance $C(t_1,t_2)=\lambda\min(t_1,t_2)$ at **L171–173** (continuous at
  $t_1=t_2$; invoke Theorem 1 / Theorem 3 of Ch. 2).
- **L189 — "It follows from (—) that $\partial^2 R/\partial t_1\partial t_2=\lambda\delta+\lambda^2$"**
  → `(—)` = $\partial R/\partial t_1=\lambda u(t_2-t_1)+\lambda^2 t_2$ at **L184–186**.
- **L204 — "Equation (—) states that $\partial^2 R/\partial t_1\partial t_2$ does not exist…"**
  → `(—)` = $\partial^2 R/\partial t_1\partial t_2=\lambda\delta(t_2-t_1)+\lambda^2$ at **L191–193**.
- **L311 — "The process $Y(t)$ defined by (—)"** → `(—)` = the shot-noise result
  $Y(t)=\sum_i L(t-t_i)$ in the aligned block at **L303–309** (equivalently the distributed
  lag $Y(t)=\int_0^\infty L(\tau)z(t-\tau)\,d\tau$ at L289–291).
- **L313 — "In Section (—), we shall show how to calculate its … moments"** *(section ref)*
  → forward reference; the first/second moments of filtered white noise / shot noise are
  computed in **Chapter 5** (the generalized-Poisson-process discussion, L483–555) and via
  spectra in **Chapter 8** ({eq}`eq-8-4`–{eq}`eq-8-5`).

---

## Chapter 05 — `05_poisson_driven_processes.md`

Three "rules" are derived in this chapter (all currently unlabeled):
- **Rule 1** (change of variables / Itô for Poisson): $d\Psi(x)=\langle\partial\Psi/\partial x,f\rangle dt+\sum_i(\Psi(x+g_i)-\Psi(x))dN_i$ at **L161–166**.
- **Rule 2** (expectations): $\tfrac{d}{dt}E x=Ef+\sum_i Eg_i\,\lambda_i$ at **L258–260**.
- **Rule 3** (autocovariance): $\tfrac{d}{d\tau}E x(t)x(t+\tau)=\dots$ at **L330–334** (≡ L321–326).

- **L159 — "from definition (—) that $\Psi(x)$ satisfies…"** *(definition ref)* → the
  **Definition** of an Itô-sense solution (L31–43); equivalently {eq}`eq-5-1`.
- **L168 — "Together with definition (—), equation (—) states that"** → definition `(—)` =
  the Itô-sense solution **Definition** (L31–43); equation `(—)` = **Rule 1** (L161–166).
- **L182 — "the rule (—) and its implication…"** → **Rule 1** (L161–166).
- **L256 — "Using this result in (—) and (—), dividing by $\Delta$…"** → first `(—)` = the
  increment expansion $E(x(t+\Delta)-x(t))=E\int f\,ds+\sum_i E\int g_i\,dN_i$ (**L236–241**);
  second `(—)` = the decomposition $\int g_i\,dN_i=E\int g_i(dN_i-\lambda_i ds)+E\int\lambda_i g_i\,ds$ (**L245–250**).
- **L262 — "As an example of the use of (—)…"** → **Rule 2** (L258–260).
- **L274 — "Applying rule (—), we find…"** → **Rule 2** (L258–260).
- **L319 — "which implies, via our rule (—), that"** → **Rule 2** (L258–260).
- **L400 — "Applying (—), we find that"** → **Rule 3** (L330–334).
- **L418 — "Applying formulas (—), we find that"** → **Rule 2** (L258–260).
- **L430 — "We use formula (—) to derive a differential equation for $x^2$"** → **Rule 1** (L161–166).
- **L445 — "Applying formula (—) to the above equation"** → **Rule 2** (L258–260).
- **L457 — "Finally, applying formula (—)…"** → **Rule 3** (L330–334).
- **L501 — "$a(t)$ plays the role of $g(x,t)$ in formula (—). Applying rule (—)…"** → first
  `(—)` = the vector Poisson SDE $dx=f\,dt+\sum_i g_i\,dN_i$ (**L153–155**, ≡ {eq}`eq-5-1`);
  second `(—)` = **Rule 1** (L161–166).
- **L507 — "applying rule (—), we find"** → **Rule 2** (L258–260).
- **L525 — "we have from rule (—) that"** → **Rule 3** (L330–334).
- **L537 — "The derivative $dy/dt$ given by (—)"** → `(—)` = $dy/dt=\sum_i a(\tau_i)\delta(t-\tau_i)$ at **L491–493**.
- **L549 — "By substituting (—) into (—)…"** → first `(—)` = $w(t)=dy/dt=\sum_i a(\tau_i)\delta(t-\tau_i)$ (**L539–541**);
  second `(—)` = $z(t)=\int_0^\infty h(\tau)w(t-\tau)\,d\tau$ (**L545–547**).

---

## Chapter 06 — `06_wiener_process.md`

- **L22 — "the mean and autocorrelation of the solution $x(t)$ of (—)"** → `(—)` = the SDE
  $dx=\tfrac{1}{\sqrt\lambda}(dN_1-dN_2)$ at **L18–20**.
- **L39 — "Using rule (—), we have that for $p>0$"** → **Rule 1** of Chapter 5 (the
  change-of-variables rule, Ch. 5 L161–166).
- **L49 — "It follows from rule (—) that"** → **Rule 2** of Chapter 5 (the expectation rule,
  Ch. 5 L258–260).
- **L65 — "Using (—) recursively to calculate odd moments"** → `(—)` = the moment recursion
  for **$p$ odd**, $\tfrac{d}{dt}Ex^p=\dots$, at **L58–63**.
- **L71 — "we drive $\lambda\to\infty$ in (—)"** → `(—)` = the moment recursion for **$p$ even**
  at **L51–56**.
- **L107 — "Comparing (—) with (—), we can conclude … Gaussian"** → first `(—)` = the computed
  moments $Ex^2=t,\;Ex^4=3t^2,\;Ex^6=15t^3$ at **L85–91**; second `(—)` = the Gaussian
  even-moment formula $Ex^p=\sigma^{p/2}(p-1)(p-3)\cdots1$ at **L103–105**.

---

## Chapter 07 — `07_wiener_driven_sde.md`

- **L22 — "regarding (—) as the limit as $\lambda\to\infty$ of"** → `(—)` = the Wiener-driven
  SDE $dx=f\,dt+g\,dW$ at **L18–20**.
- **L30 — "the counterpart of rule (—)"** → **Rule 1** of Chapter 5 (change-of-variables,
  Ch. 5 L161–166).
- **L39 — "we use (—) for $\lambda>0$ and using rule (—) to obtain"** → first `(—)` = the
  Poisson-approximation SDE $dx=f\,dt+\tfrac{1}{\sqrt\lambda}g(dN_1-dN_2)$ at **L24–26**;
  second rule `(—)` = **Rule 1** of Chapter 5 (Ch. 5 L161–166). *(The result obtained is
  {eq}`eq-7-1`.)*
- **L167 — "Taking the limit as $\lambda\to\infty$ gives the desired result (—)"** → `(—)` =
  "our second rule," $\tfrac{d}{dt}E\Psi(x)=E\langle\partial\Psi/\partial x,f\rangle+\tfrac12 E\langle\partial^2\Psi/\partial x^2 g,g\rangle$, at **L136–138**.

---

## Chapter 08 — `08_spectral_densities.md`

- **L293 — "{eq}`eq-8-a`, {eq}`eq-8-b` and the convolution property (—)"** → `(—)` = the
  filtering/convolution spectral relation {eq}`eq-8-5` ($S_y(w)=h(w)S_x(w)h(-w)$);
  equivalently the convolution property (property 6) of **Table 2**.
- **L388 — "attempt to invert (—), and to solve for $v(t)$"** → `(—)` = the non-fundamental
  representation $x(t)=\dfrac{D-b}{(D-\lambda_1)(D-\lambda_2)}v(t)$ (L330–334), i.e.
  $x(t)=\int_0^\infty r(\tau)v(t-\tau)\,d\tau$ (L338–340).
- **L453 — "demonstrate how (—) can be inverted to express $w(t)$…"** → `(—)` = the Wold
  representation $x(t)=\int_0^\infty p(\tau)w(t-\tau)\,d\tau$ at **L425–427**.

---

## Chapter 09 — `09_characterizations_ms_differentiability.md`

- **L50 — "by repeated differentiation of (—)"** → `(—)` = the first-derivative
  representation {eq}`eq-9-1`, $Dx(t)=\int_0^\infty p'(\tau)w(t-\tau)\,d\tau+p(0)w(t)$
  (also tagged `(+)` at L21).
- **L96 — "Together with results (—), the following characterization…"** *(results ref)* →
  the differentiability conditions established earlier in the chapter: the necessary
  conditions $p(0)=0$, $\int_0^\infty|Dp|\,ds<\infty$ (L30) and the higher-order theorem
  $p(0)=Dp(0)=\cdots=D^{n-1}p(0)=0$ (L52–58).

---

## Chapter 10 — `10_cramer_representation.md`

- **L36 — "we use (—) to calculate $R(\tau)$"** → `(—)` = the Cramér representation
  {eq}`eq-10-cramer`.
- **L53 — "This is the inversion formula (—)"** → the spectral-density inversion formula —
  here re-derived as $R(\tau)=\tfrac{1}{2\pi}\int e^{i\mu\tau}S(\mu)\,d\mu$ (L49–51);
  originally {eq}`eq-8-2` in Chapter 8.
- **L61 — "In (—), we have added the argument $w\in\Omega$ explicitly"** → `(—)` = the
  Fourier-transform definition $Z'(\lambda,w)=\tfrac{1}{\sqrt{2\pi}}\int e^{-i\lambda t}x(t,w)\,dt$
  at **L57–59**.
- **L81 — "Then using (—), we have"** → `(—)` = the Cramér representation {eq}`eq-10-cramer`
  (used to substitute for $x(t-\tau)$).

---

## Chapter 12 — `12_prediction.md`

*(All displayed equations in this chapter are unlabeled.)*

- **L39 — "Since (—) holds for all $t$"** → `(—)` = the Wold MA representation
  $x(t)=\int_0^\infty p(s)w(t-s)\,ds$ at **L10–12**.
- **L51 — "Equation (—) is the continuous time Wiener–Kolmogorov formula"** → `(—)` =
  $E[x(t+u)\mid x(v),v\le t]=\int_0^\infty p(s+u)w(t-s)\,ds$ at **L47–49**.
- **L67 — "As an example of the use of formula (—)"** → the Wiener–Kolmogorov formula
  (L47–49), in operator form $E_t x(t+u)=[\tilde P(D)e^{-Du}]_+ w(t)$ (L55–57).
- **L74 — "Then formula (—) gives"** → same Wiener–Kolmogorov formula (L47–49 / L55–57).
- **L98 — "equation (—) implies that the … forecast … is given by"** → `(—)` = the prediction
  result $E_t x(t+u)=e^{-au}x(t)$ at **L85–87**.
- **L100 — "the geometric distributed lead (—)"** → `(—)` = $E_t\int_0^\infty e^{\rho u}x(t+u)\,du$
  at **L93–95**.
- **L106 — "An approach to the evaluation of (—)"** → `(—)` = the geometric distributed lead
  $E_t\int_0^\infty e^{\rho u}x(t+u)\,du$ at **L93–95**.
- **L120 — "a partial fraction representation of the right side of (—)"** → `(—)` =
  $x(t)^*=\bigl(\tfrac{-1}{\rho+D}\bigr)\bigl(\tfrac{1}{a+D}\bigr)w(t)$ at **L110–115**.
- **L141 — "Represent equation (—) as"** → `(—)` = the forecast $E_t x(t)^*=\bigl(\tfrac{1}{a-\rho}\bigr)\int_0^\infty e^{-as}w(t-s)\,ds=\bigl(\tfrac{1}{a-\rho}\bigr)x(t)$ at **L134–139**.
- **L150 — "Equation (—) holds for *any* $\tilde P(D)$"** → `(—)` = the operator formula
  $E_t x(t)^*=\bigl[\tfrac{-\tilde P(D)+\tilde P(-\rho)}{D+\rho}\bigr]w(t)$ at **L143–148**.
- **L159 — "the generalization of (—) is"** → `(—)` = the same operator formula at **L143–148**.
- **L165 — "This generalization of (—) will be established in chapter (—)."** → first `(—)` =
  $E_t\int_0^\infty e^{\rho s}x(t+s)\,ds=\bigl[\tfrac{-\tilde P(D)+\tilde P(-\rho)}{D+\rho}\bigr]w(t)$
  at **L161–163**; "chapter (—)" *(chapter ref)* = forward reference, most likely
  **Chapter 17** ("Prediction Formulas for Continuous Time Linear Rational Expectations Models").

---

## Chapter 13 — `13_locally_unpredictable.md`

- **L62 — "Using the preceding theorem and our formulas (—)"** → `(—)` = the
  geometric-distributed-lead prediction formula of **Chapter 12**
  ($E_t\int_0^\infty e^{\rho s}x(t+s)\,ds=\bigl[\tfrac{-\tilde P(D)+\tilde P(-\rho)}{D+\rho}\bigr]w(t)$,
  Ch. 12 L161–163), together with the **initial-value theorem** of Chapter 9 (L98–102).
- **L93 — "by the assumption that (—) is a Wold representation"** → `(—)` = the Wold
  representation $x(t)=\int_0^\infty p(s)w(t-s)\,ds=\tilde P(D)w(t)$ at **L66–68**.

---

## Chapter 14 — `14_nonstationary_examples.md`

- **L24 — "It follows from (—) that $x(t)$ can also be expressed as"** → `(—)` = the MA
  representation $x(t)=\int_0^\infty p(\tau)w(t-\tau)\,d\tau$ at **L5–7** (with the start-up
  white-noise spec, L11–16).
- **L36 — "It can be verified from (—), using the criterion in theorem (—)"** → first `(—)` =
  the autocovariance $E x(t)x(t-\tau)=\int_{-T}^{t-\tau}p(t-s)p(t-\tau-s)\,ds$ at **L32–34**
  (≡ the representation $x(t)=\int_{-T}^t p(t-s)w(s)\,ds$, L26–28); "theorem (—)" *(theorem
  ref)* = the mean-square-differentiability criterion — **Theorem 6** of Chapter 2 (the
  nonstationary criterion), reinforced by the $p(0)=0$ conditions of Chapter 9.
- **L56 — "Applying our prediction formula to (—)"** → `(—)` = the process
  $x(t)=\int_0^{t+T}(1+\beta s)w(t-s)\,ds$ at **L46–48** (≡ $x(t)=\tfrac{\beta+D}{D^2}w(t)$, L40–42).
- **L78 — "obtained by differentiating the right side (—) formally"** → `(—)` = the forecast
  $\hat E_t x(t+v)=\int_0^{t+T}(1+\beta(s+v))w(t-s)\,ds$ at **L58–60**.
- **L86 — "Solving (—) formally for $w(t)$, and using the result in (—) gives"** → first `(—)` =
  the process definition $x(t)=\tfrac{\beta+D}{D^2}w(t)$ at **L40–42** (solve $w=\tfrac{D^2}{\beta+D}x$);
  second `(—)` = $\tfrac{d}{dv}E_t x(t+v)=\tfrac{\beta}{D}w(t)$ at **L82–84**.
- **L106 — "noted by Muth (—) in discrete time"** *(citation)* → **Muth (1960)**, "Optimal
  Properties of Exponentially Weighted Forecasts," *JASA* 55, 299–306.
- **L106 — "a special property of the stochastic process (—)"** → `(—)` = the process
  $x(t)=\tfrac{\beta+D}{D^2}w(t)$ at **L40–42** (≡ L46–48).

---

## Chapter 15 — `15_discrete_sampling_folding.md`

- **L27 — "Equation (—) can also be represented as $R^d(\tau)=R(\tau)S_T(\tau)$"** → `(—)` =
  $R^d(\tau)=\sum_n R(nT)\delta(\tau-nT)$ at **L20–23** (current).
- **L39 — "From (—), we can express the spectral density…"** → `(—)` = the same train-of-deltas
  $R^d(\tau)=\sum_n R(nT)\delta(\tau-nT)$ (the definition of $R^d$), leading into {eq}`eq-15-1`.
- **L72 — "It follows from (—) and (—) that the spectral density … satisfies [folding formula]"**
  → first `(—)` = {eq}`eq-15-1` ($S^d(w)=\sum_n R(nT)e^{-iwnT}$); second `(—)` = the
  Poisson-summation pair {eq}`eq-15-2`/{eq}`eq-15-3` (and the inverse-transform line
  $\tfrac{1}{2\pi}\sum_n R(\tau)\delta(\tau-nT)=\tfrac{1}{2\pi}R^d(\tau)$ just below {eq}`eq-15-3`).
- **L80 — "Equation (—) is known as the *folding formula*"** → `(—)` = the immediately
  preceding displayed equation $S^d(w)=\tfrac{1}{T}\sum_n S\!\left(w-n\tfrac{2\pi}{T}\right)$
  (the line above the sentence).
