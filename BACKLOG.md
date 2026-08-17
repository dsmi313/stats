# Backlog

Everything cut from the earlier 22-section construction plan, plus what has moved
in and out since. Nothing is deleted; each item is here with one line saying why
it left the talk. If any of these earns its way back, an existing session has to
leave to make room.

## This round: two in, two out

The talk stays at nine sessions. Two sessions were added — the real-data run
(session 7) and composition two ways (session 8) — so two had to leave:

- Four intervals for one parameter (old session 5) — merged into session 3, which
  now reads all four intervals off the same rate while keeping the profile curve
  as its figure. One script instead of two, no interval recipe lost.
- Monte Carlo marginalization (old session 7) — cut to the backlog. Its identity,
  that drawing a posterior value and pushing it through the calculation equals
  integrating the nuisance out, is shown by the grid integral in session 4 and is
  the working engine of the coverage study in session 6, so a standalone numerical
  demonstration is more than the talk needs.

Both additions were promotions of material the backlog already accounted for:
composition (Section 7 and Sections 8 to 10 below) and the real-data run (Section
14 below), rather than newly invented sessions.

## Folded into the nine sessions

- Notation glossary — folded into the session-2 skeleton table, where each
  symbol is defined against its real production name.
- Section 2, MLE for the window-count estimator — this is now session 1.
- Section 3, variance and the delta method — the delta interval survives as one
  of the four intervals in session 3; the standalone section is cut.
- Section 4, parametric bootstrap and coverage — the mechanism lives in session 6;
  the standalone section is cut.
- Section 5, nonparametric bootstrap and distributional diagnostic — same, the
  mechanism lives in session 6.
- Section 7, multinomial composition — promoted this round to session 8,
  composition two ways, after starting folded into the skeleton table.
- Sections 8 to 10, the Delomas and Hess likelihood in three parts — the PBT
  multinomial likelihood is now built by hand in session 8 and anchored to
  `PBT_expand_calc_MLE()`; the deep three-part build is still more than the talk
  needs.
- Section 11, SCRAPI ratio estimator — folded into session 1 and the skeleton.
- Section 12, the compound bootstrap — folded into session 6, where composed
  intervals actually get examined.
- Section 16, escapeLGD internals — promoted, not cut: reading `escapeLGD` is now
  the spine of every session's Locate part rather than one late section.
- Section 19, Bayesian sidebar — folded into session 3 (Bayesian credible
  interval) and session 5 (Bayesian GE fit).
- Section 21, year-one synthesis memo — replaced by session 9, talk assembly.

## Cut to backlog

- Section 6, stratification and when pooling lies — a real issue in both tools,
  but the talk does not need a standalone treatment of it to make its point.
- Monte Carlo marginalization — cut this round; see "two in, two out" above.
- Section 13, reading SCOBI source and the lgr2SCRAPI audit — a code-provenance
  audit, not something the talk argues.
- Section 14, end-to-end run on MY2024 — the MY2025 real-data step it was to
  supersede is now session 7 in its own right; a second full run on MY2024 stays
  redundant and cut.
- Section 15, fishCompTools as a parallel implementation — a third
  implementation to compare against is out of scope for a two-tool talk.
- Section 17, FSA broader fisheries bootstrap and stratification idioms — general
  fisheries background, not load-bearing for the argument.
- Section 18, TropFishR MLE in a different fisheries domain — a cross-domain
  detour; interesting, not needed.
- Section 20, STADEM, the state-space answer at the dam — a different modeling
  paradigm; a comparison to it would be its own talk.
- Section 22, the integrated Bayesian routing and composition model — the model
  is already built and in production, so the talk reverse engineers it rather
  than constructing it; its architecture write-up is not needed to give the talk.
- `R/PartI/mit_studios/` — the fish-themed port of the MIT 18.05 studios is
  useful scaffolding but is learning practice, not talk material.
