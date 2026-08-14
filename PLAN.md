# Reverse engineering EASE and SCRAPI2

This repository has one purpose: to get me able to give a 45-minute talk to
IDFG colleagues, working title *"uncertainty in escapement estimation at Lower
Granite: adults and smolts, one framework."* Every session here exists to make
some part of that talk explainable out loud. If a session does not serve the
talk, it is not here; it is in `BACKLOG.md` with one line saying why.

## The one idea the talk is built on

EASE for adults (`escapeLGD`) and SCRAPI2 for smolts (`smoltEASE`) are not two
topics. They are one estimator structure applied at two life stages. The whole
talk is the claim that the same skeleton runs at both stages, that the pieces
line up almost row for row, and that the few places they diverge are the
interesting places. The generic skeleton, in the order the estimators apply it:

1. A total count at the dam (window count for adults, trap count for smolts).
2. A sampled fraction that the count is expanded by (counting-window open
   fraction for adults, trap sample rate for smolts).
3. A detection/observation expansion for fish the count structurally misses
   (nighttime passage for adults, guidance efficiency into the bypass for
   smolts).
4. Composition proportions that split the expanded total into groups (PBT and
   genetic stock, rear type, age).
5. A stage-specific adjustment (fallback and reascension for adults; wild
   fraction split for smolts).
6. Uncertainty propagated to an interval by resampling and posterior draws.

Both production tools are already built and in production. This plan is not a
path toward building them. It is reverse engineering: every session starts from
a toy version I write by hand, then points at the exact production file and
function that does the same job, and says how the two production
implementations handle it and where they part ways.

## Where every description of EASE comes from

Every description of EASE in this repository is derived from the `escapeLGD`
source, never from general knowledge and never from the package README. The
same holds for SCRAPI2 and the `smoltEASE` source. The "Locate" part of each
session names files and functions by their real names in those two repositories.

## Session format

Every session is one self-contained R script in `R/`, written MIT 18.05 style:
a stated objective, a worked example, and a short exercise. Each follows the
Kery and Schaub simulate-build-function loop. Five parts, in this order:

1. **Objective.** One sentence naming what I should be able to explain out loud
   afterward. Not what I should be able to code.
2. **Simulate.** Generate data from known truth. Truth is a variable at the top
   of the script.
3. **Build.** Write the estimator by hand. The likelihood is written out and
   optimized with `optim` or `nlminb`. No package does the estimating work.
4. **Recover.** Show the estimate against known truth, show the interval, plot
   it.
5. **Locate.** Name the exact file and function in **both** `escapeLGD` and
   `smoltEASE` that does this same job in production. In two sentences per side,
   say how each production version differs from the toy, and whether the two
   production implementations handle it the same way.

Each script also writes `docs/sessionNN_explain.md`: how I would explain that
one concept to a colleague in three minutes, in plain English, with no notation.

## The nine sessions

The sequence is centered on likelihood profiling. Sessions 1 and 2 lay the
generic skeleton; sessions 3 and 4 are the profiling core; sessions 5 through 8
are the four ways uncertainty becomes an interval and what those intervals
actually claim; session 9 assembles the talk.

### Session 1 — One count, expanded, by likelihood

**Objective:** explain why a count at a dam divided by the fraction of fish it
saw is a maximum likelihood estimate, and why that single idea is the base of
both the adult and the smolt estimator.

Binomial window count, likelihood written by hand, MLE recovered. This session
already existed under the old plan and is rewritten to the new format. The toy
estimator is written generically enough that both production estimators map
onto it: quantities are named generically (total count, sampled fraction,
composition proportions, expansion, uncertainty sources), not with
smolt-specific names.

**Locate:** `escapeLGD` `expand_wc_binom_night()` in
`R/night_fall_reascend_wc_binom.R`; `smoltEASE` `thetahat()` inside `SCRAPI2()`
in `R/SCRAPI2.R` (`dailypass <- Tally / Ptrue`).

### Session 2 — The shared skeleton

**Objective:** explain, component by component, how the adult and smolt
estimators are the same machine, and name the two or three places they are
genuinely not.

No simulation here. The deliverable is one table, derived from reading both
sources, with one row per estimator component. Columns: the generic component;
how `escapeLGD` implements it for adults; how `smoltEASE` implements it for
smolts; and whether the two differ structurally or only in parameterization.
Every uncertainty source each tool propagates is listed with its mechanism
(bootstrap, posterior draws, delta method, or not propagated at all), and any
source one propagates while the other does not is flagged. Written to
`docs/session02_shared_skeleton.md`. This table is the spine of the talk.

### Session 3 — Profile likelihood

**Objective:** explain how an interval falls out of the shape of the likelihood
curve for one parameter, without any normal approximation.

Plot the log-likelihood curve for one parameter, drop the line where twice the
log-likelihood falls by the chi-squared cutoff, and read the interval off the
curve.

**Locate:** the single-parameter binomial rate inside `escapeLGD`
`nightFall()` (nighttime passage, `p_night`) and the fallback rate `p_fa`; the
per-rate resampling in `smoltEASE` `SCRAPI2()`.

### Session 4 — Profile versus marginalize a nuisance parameter

**Objective:** explain the difference between profiling out a nuisance parameter
and integrating it out, and show they give nearly the same interval here and why.

Fit the two-parameter joint likelihood surface, then on one figure show the
profile trace for the parameter of interest against the marginalized curve.
The production anchor is exact: `escapeLGD`'s fallback likelihood is literally a
two-parameter binomial likelihood in P(fallback) and the nuisance P(reascend |
fallback). Sessions 3 and 4 both point at this one function; keeping them
separate holds each script inside the 80-line and one-minute budget while
staying teachable, which is why the plan is nine sessions rather than folding
the surface into session 3 for eight.

**Locate:** `escapeLGD` `fallback_log_likelihood()` and
`gradient_fallback_log_likelihood()` in `R/fallback_reascend_likelihood.R`,
optimized in `nightFall()`; `smoltEASE` has no two-parameter reascension
likelihood, which is itself a point the talk makes.

### Session 5 — Four intervals for one parameter

**Objective:** explain what each of the four intervals means and why they agree
here and can disagree elsewhere.

One parameter, one plot, four intervals: delta method, profile likelihood,
bootstrap percentile, and Bayesian credible.

**Locate:** bootstrap percentile is the interval both production tools actually
ship (`quantile()` on the bootstrap matrix in `SCRAPI2()` and in
`apply_fallback_rates()`); the profile and delta intervals are what the
by-hand likelihood work implies; the Bayesian credible interval is what
`smoltEASE` `fit_ge_model()` produces for GE.

### Session 6 — Guidance efficiency fit both ways

**Objective:** explain what guidance efficiency is, why it needs its own model,
and what fitting it by maximum likelihood versus Bayesian buys and costs.

Fit the GE relationship both ways on simulated data with known truth: maximum
likelihood in `glmmTMB` and Bayesian by hand. Simulate-only, so the truth stays
known and the two fits can be judged against it; the truth values are drawn to
look like the real MY2025 GE-versus-spill shape (see `data/`) without the script
depending on the data to run.

**Locate:** `smoltEASE` `fit_ge_model()` (JAGS multistate mark-recapture for the
route-selection probability) and `prep_ge_data()` in `R/fit_ge_model.R` and
`R/prep_ge_data.R`; `escapeLGD` has no GE model, because nighttime passage plays
the structurally analogous role and is estimated as a plain binomial rate in
`nightFall()`.

### Session 7 — Monte Carlo marginalization

**Objective:** explain how drawing posterior values and pushing each through the
escapement calculation is the same thing as integrating the nuisance parameter
out, demonstrated numerically.

Find the loop in `smoltEASE` that draws GE and GSI posterior values and pushes
them through the escapement calculation, and show numerically that it equals the
integral from session 4. Do the same for whatever `escapeLGD` draws or resamples.

**Locate:** `smoltEASE` `SCRAPI2()` bootstrap loop, where `ge_day_mat[, b]` and
the per-iteration GSI draw column enter `thetahat()`, fed by
`generate_ge_draws()` and `sim_gsi_draws()`; `escapeLGD` `HNC_expand_unkGSI()`
(GSI posterior draw columns) and the binomial resampling in `nightFall()`.

### Session 8 — What the composed interval actually claims

**Objective:** explain what the SCRAPI2 interval claims to cover and whether a
simulation says it delivers that coverage, and how the adult interval compares.

The SCRAPI2 interval is composed of a nonparametric bootstrap stacked on
posterior draws of GE and GSI. State what that composed interval claims, then
run a simulation study to check whether it has the coverage the claim implies.
Do the same for the adult interval from `escapeLGD` and compare.

**Locate:** `smoltEASE` `SCRAPI2()` CI construction (`quantile(theta.b, ...)`
over a bootstrap that folds in GE and GSI draws); `escapeLGD`
`apply_fallback_rates()` CI construction (`quantile()` over composition
bootstrap times fallback bootstrap).

### Session 9 — Talk assembly

**Objective:** deliver the talk from two or three figures.

Assemble the two or three figures that carry the whole thing: the shared
skeleton, the profiling core, and the coverage result. No new estimator work;
this session selects and finishes figures already produced.

## Rules this repository holds itself to

- No script longer than 80 lines. No script that takes more than a minute to run.
- No new session is added without deleting one. Tangents go to `BACKLOG.md`.
- No banner comments (no `# =====`, no `# -----`). Plain comments only.
- Every heading in every markdown file has real text. No empty or placeholder
  headings.
- No dates are set anywhere in this repository. I add dates myself.

## Repository layout

- `PLAN.md` — this file.
- `README.md` — one-screen orientation.
- `BACKLOG.md` — everything cut from the earlier 22-section plan, each with one
  line saying why it was cut.
- `R/sessionNN_*.R` — one self-contained script per session.
- `docs/sessionNN_explain.md` — the three-minute plain-English explanation each
  session writes.
- `docs/session02_shared_skeleton.md` — the component-by-component comparison
  table.
- `data/` — real MY2025 steelhead inputs, kept as the reference shapes the
  simulations imitate. No session depends on them to run; see `data/README.md`.
</content>
</invoke>
