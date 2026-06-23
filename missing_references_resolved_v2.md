# Missing references — what was done (v2)

This is the follow-up to `missing_references_resolved.md`. Every `(—)` placeholder in
chapters 01–16 has now been **resolved in the source files**. This document reports what
each placeholder became.

**Method.** Where a placeholder pointed at an unlabeled displayed equation, I gave that
equation a MyST label using the dollar-math label syntax `$$ … $$ (label)` and replaced the
`(—)` with a `{eq}`label`` cross-reference. Where it pointed at a theorem, section,
chapter, or citation, I substituted the concrete reference in prose.

**Verification.** A full `jupyter-book build .` succeeds with **zero "undefined label"
warnings**; the only build warnings are the two pre-existing "document isn't included in any
toctree" notices (`18_temporal_aggregation.md` and these report files). All 47 new labels
resolve, including the cross-chapter references.

**Totals.** 0 `(—)` remain in chapters 01–16. 47 new equation labels were introduced:

| Chapter | New labels |
|---|---|
| 02 | `eq-2-analytic`, `eq-2-plus`, `eq-2-xhat`, `eq-2-Rm`, `eq-2-R0` |
| 03 | `eq-3-R`, `eq-3-C`, `eq-3-dR`, `eq-3-d2R`, `eq-3-shot` |
| 05 | `eq-5-rule1`, `eq-5-rule2`, `eq-5-rule3`, `eq-5-Eincr`, `eq-5-decomp`, `eq-5-dydt`, `eq-5-wt`, `eq-5-z` |
| 06 | `eq-6-sde`, `eq-6-even`, `eq-6-odd`, `eq-6-moments`, `eq-6-gauss` |
| 07 | `eq-7-sde`, `eq-7-approx`, `eq-7-rule2` |
| 08 | `eq-8-nonfund`, `eq-8-wold` |
| 10 | `eq-10-Zhat` |
| 12 | `eq-12-wold`, `eq-12-wk`, `eq-12-wkop`, `eq-12-ar1pred`, `eq-12-glead`, `eq-12-xstar`, `eq-12-Estar`, `eq-12-Estarop`, `eq-12-gen` |
| 13 | `eq-13-wold` |
| 14 | `eq-14-ma`, `eq-14-acov`, `eq-14-proc`, `eq-14-procint`, `eq-14-forecast`, `eq-14-dvD` |
| 15 | `eq-15-Rd`, `eq-15-fold` |

(Chapters 09 used the existing `eq-9-1`; chapter 10 also reuses the existing `eq-10-cramer`
and the chapter-8 label `eq-8-2`; chapters 01, 04, 11, 16 had no placeholders.)

---

## Chapter 02
- "Theorem (—)" → **Theorem 9** (numbered after the existing Theorem 8).
- The proof's ad-hoc tags were also wired up: the analytic series is now `eq-2-analytic`
  (was `(*)`), the squared-error identity `eq-2-plus` (was `(+)`), the $\hat x$ series
  `eq-2-xhat`, the differentiated series `eq-2-Rm` (was `(a)`), and `eq-2-R0` (was `(b)`,
  with the trailing `(b)` tag removed).
- The closing sentence now reads: "Substituting the right side of {eq}`eq-2-xhat` into
  {eq}`eq-2-plus`, noting by the reasoning that led to **Theorem 8** … and using
  {eq}`eq-2-Rm` and {eq}`eq-2-R0` to evaluate the two terms in braces in {eq}`eq-2-plus`."

## Chapter 03
- "follows from (—) that $C=\lambda\min$" → `{eq}`eq-3-R``.
- "From (—) … mean square continuous" → `{eq}`eq-3-C``.
- "It follows from (—) that $\partial^2R…$" → `{eq}`eq-3-dR``.
- "Equation (—) states that $\partial^2R…$ does not exist" → `{eq}`eq-3-d2R``.
- "The process $Y(t)$ defined by (—)" → `{eq}`eq-3-shot``.
- "In Section (—), we shall show how to calculate its … moments" → **"In Chapter 8"**
  *(inferred forward reference — see flags below)*.

## Chapter 05
- The three manipulation rules are now `eq-5-rule1` (change of variables), `eq-5-rule2`
  (expectations), `eq-5-rule3` (autocovariance); the increment expansion is `eq-5-Eincr`,
  its decomposition `eq-5-decomp`, and the generalized-Poisson objects `eq-5-dydt`,
  `eq-5-wt`, `eq-5-z`.
- "from definition (—)" / "Together with definition (—)" → prose **"the definition of an
  Itô-sense solution"** (of `{eq}`eq-5-1``), since the solution definition is an admonition,
  not a numbered equation.
- All "rule/formula (—)" → the matching `{eq}`eq-5-rule1/2/3`` per the v1 mapping; the two
  expansion references → `{eq}`eq-5-Eincr`` and `{eq}`eq-5-decomp``; "$g$ … in formula (—)"
  → `{eq}`eq-5-1``; "$dy/dt$ given by (—)" → `{eq}`eq-5-dydt``; "substituting (—) into (—)"
  → `{eq}`eq-5-wt`` into `{eq}`eq-5-z``.

## Chapter 06
- "solution $x(t)$ of (—)" → `{eq}`eq-6-sde``.
- "Using rule (—)" → `{eq}`eq-5-rule1`` *(cross-chapter)*; "It follows from rule (—)" →
  `{eq}`eq-5-rule2`` *(cross-chapter)*.
- "Using (—) recursively … odd moments" → `{eq}`eq-6-odd``; "drive $\lambda\to\infty$ in
  (—)" → `{eq}`eq-6-even``.
- "Comparing (—) with (—)" → `{eq}`eq-6-moments`` with `{eq}`eq-6-gauss``.

## Chapter 07
- "regarding (—) as the limit" → `{eq}`eq-7-sde``.
- "counterpart of rule (—)" → `{eq}`eq-5-rule1`` *(cross-chapter)*.
- "we use (—) … and using rule (—)" → `{eq}`eq-7-approx`` and `{eq}`eq-5-rule1``.
- "the desired result (—)" → `{eq}`eq-7-rule2``.

## Chapter 08
- "the convolution property (—)" → `{eq}`eq-8-5`` (the filtering relation
  $S_y=h(w)S_x(w)h(-w)$).
- "attempt to invert (—)" → `{eq}`eq-8-nonfund`` (the non-fundamental representation).
- "how (—) can be inverted to express $w(t)$" → `{eq}`eq-8-wold``.

## Chapter 09
- "by repeated differentiation of (—)" → `{eq}`eq-9-1``.
- "Together with results (—)" → prose **"Together with the preceding theorems on the
  existence of mean square derivatives"** (the targets are this chapter's unnumbered
  theorems, so a prose reference replaces the placeholder).

## Chapter 10
- "we use (—) to calculate $R(\tau)$" → `{eq}`eq-10-cramer``.
- "the inversion formula (—)" → `{eq}`eq-8-2`` *(cross-chapter — the inversion formula
  defined in Chapter 8)*.
- "In (—), we have added the argument $w\in\Omega$" → `{eq}`eq-10-Zhat``.
- "Then using (—)" → `{eq}`eq-10-cramer``.

## Chapter 12
- New labels for the nine prediction equations (`eq-12-wold`, `eq-12-wk`, `eq-12-wkop`,
  `eq-12-ar1pred`, `eq-12-glead`, `eq-12-xstar`, `eq-12-Estar`, `eq-12-Estarop`,
  `eq-12-gen`).
- "Since (—) holds" → `{eq}`eq-12-wold``; "Equation (—) is the … Wiener–Kolmogorov formula"
  → `{eq}`eq-12-wk``; "use of formula (—)" / "Then formula (—) gives" → `{eq}`eq-12-wk``;
  "equation (—) implies … lead (—)" → `{eq}`eq-12-ar1pred`` and `{eq}`eq-12-glead``;
  "evaluation of (—)" → `{eq}`eq-12-glead``; "right side of (—)" → `{eq}`eq-12-xstar``;
  "Represent equation (—)" → `{eq}`eq-12-Estar``; "Equation (—) holds for any $\tilde P$" /
  "generalization of (—)" → `{eq}`eq-12-Estarop``.
- "This generalization of (—) will be established in chapter (—)" → "This generalization of
  `{eq}`eq-12-gen`` will be established in **Chapter 17**" *(inferred chapter reference — see
  flags)*.

## Chapter 13
- "our formulas (—)" → "our formula `{eq}`eq-12-gen`` for geometric distributed leads"
  *(cross-chapter)*.
- "the assumption that (—) is a Wold representation" → `{eq}`eq-13-wold``.

## Chapter 14
- New labels `eq-14-ma`, `eq-14-acov`, `eq-14-proc`, `eq-14-procint`, `eq-14-forecast`,
  `eq-14-dvD`.
- "It follows from (—)" → `{eq}`eq-14-ma``; "verified from (—), using the criterion in
  theorem (—)" → "verified from `{eq}`eq-14-acov``, using the criterion in **Theorem 6 of
  Chapter 2**"; "prediction formula to (—)" → `{eq}`eq-14-procint``; "right side (—)
  formally" → `{eq}`eq-14-forecast``; "Solving (—) … using the result in (—)" →
  `{eq}`eq-14-proc`` and `{eq}`eq-14-dvD``.
- "noted by Muth (—)" → **"Muth (1960)"**; "the stochastic process (—)" → `{eq}`eq-14-proc``.

## Chapter 15
- New labels `eq-15-Rd` (the delta-train autocovariance) and `eq-15-fold` (the folding
  formula).
- "Equation (—) can also be represented" / "From (—), we can express" → `{eq}`eq-15-Rd``.
- "It follows from (—) and (—)" → `{eq}`eq-15-1`` and `{eq}`eq-15-3``.
- "Equation (—) is known as the folding formula" → `{eq}`eq-15-fold``.

---

## Cross-chapter references introduced
These `{eq}` references point to labels in other chapters (all resolve in the build):
`eq-5-rule1` (used in ch. 6, 7), `eq-5-rule2` (ch. 6), `eq-8-2` (ch. 10), `eq-12-gen`
(ch. 13).

## Inferences you should double-check
These three resolutions were **inferred from context**, not stated in the source — please
confirm they match your intent:

1. **Ch. 3, "In Section (—)"** → rendered as **"In Chapter 8."** The shot-noise first/second
   moments are most naturally computed with the chapter-8 filtering machinery
   ({eq}`eq-8-4`–{eq}`eq-8-5`), but Chapter 5's generalized-Poisson discussion is also a
   candidate.
2. **Ch. 12, "in chapter (—)"** → rendered as **"Chapter 17."** That companion paper
   ("Prediction Formulas for Continuous Time Linear Rational Expectations Models") is where
   the geometric-lead generalization is established in this volume; in the original
   manuscript it may have pointed to a different later chapter.
3. **Ch. 14, "Muth (—)"** → rendered as **"Muth (1960)."** No bibliography entry was added
   (these chapters carry no reference list); if you want a full citation,
   Muth, J. F. (1960), "Optimal Properties of Exponentially Weighted Forecasts," *JASA* 55,
   299–306.

## Style notes
- A few placeholders pointed at objects with no equation number (the Itô-solution
  *definition* in Ch. 5; the unnumbered *theorems* in Ch. 9). These were resolved with prose
  references rather than `{eq}` links, since there is no numbered target to point at. If you
  prefer, those definitions/theorems could be converted to labeled `{prf:definition}` /
  numbered theorems and cross-referenced.
- New labels use descriptive slugs (e.g. `eq-5-rule1`, `eq-12-glead`) rather than sequential
  numbers, to avoid collisions with the chapters' existing `eq-N-k` labels and to keep the
  source readable. The *displayed* equation numbers are still assigned automatically by
  Sphinx in source order.

---

## Addendum — numbered & cross-referenced definitions/theorems

The "style note" follow-up: the two placeholders that pointed at *unnumbered* objects (the
Itô-solution definition in Ch. 5; the differentiability theorems in Ch. 9) have now been
converted to **numbered, anchored** objects and the prose references replaced with real
hyperlinks. `sphinx-proof` is not installed, so rather than `{prf:}` directives I used the
book's existing hand-numbered **bold** style plus MyST `(label)=` target anchors, and
referenced them with the `{ref}`text <label>`` role.

To keep the book's two running sequences gap-free, the intervening unnumbered items were
also numbered:

- **Definitions** now run **1–10** with no gaps: Ch. 1 (1–2), Ch. 2 (3–6), Ch. 4 (7–8),
  Ch. 5 (**9**, the Itô-sense solution, anchor `def-ito-solution`), Ch. 13 (10).
- **Theorems** now run **1–16** with no gaps: Ch. 2 (1–9), Ch. 8 (10), Ch. 9
  (**11** anchor `thm-msd-representation`, **12** anchor `thm-msd-higher-order`, 13, 14),
  Ch. 11 (15), Ch. 13 (16).

Cross-references made live (each now renders as an internal hyperlink):
- Ch. 5: the two references to the Itô-sense solution now link to it, rendering as
  "Definition 9" — in "it follows immediately from Definition 9 …" and "Together with
  Definition 9, equation …".
- Ch. 9: the testing remark now links to the two differentiability theorems, rendering as
  "Together with Theorem 11 and Theorem 12, the following characterization …".

Named theorems/criteria were left as-is (e.g. the Cauchy criterion in Ch. 2, the
"Spectral factorization theorem" in Ch. 8, the "Criterion (Initial value theorem)" in
Ch. 9), consistent with the manuscript's use of named results alongside numbered ones.

The build resolves all three new `{ref}` links to internal hyperlinks ("Definition 9",
"Theorem 11", "Theorem 12") with no warnings.
