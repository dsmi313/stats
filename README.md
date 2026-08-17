# stats

Working repository for a 45-minute talk to IDFG colleagues: *"uncertainty in
escapement estimation at Lower Granite: adults and smolts, one framework."*

## What this is

Nine self-contained R sessions that reverse engineer the two production
escapement tools at Lower Granite Dam — `escapeLGD` (EASE, adults) and
`smoltEASE` (SCRAPI2, smolts) — treating them as one estimator structure applied
at two life stages. Most sessions write a toy version by hand, then point at the
exact production file and function that does the same job; session 7 instead runs
production `SCRAPI2()` on the real MY2025 inputs and checks it by hand.

## Where to start

1. Read `PLAN.md`. It states the one idea the talk is built on, the five-part
   session format, and the nine sessions in order.
2. Read `docs/session02_shared_skeleton.md` for the component-by-component
   comparison of the adult and smolt estimators. That table is the spine.
3. Work the sessions in order in `R/`. Run each from the repository root
   (`source("R/session01_count_expansion.R")`); each prints its recovery, writes
   a figure to `figs/`, and writes a plain-English explanation to
   `docs/sessionNN_explain.md`. Session 2 also writes the skeleton table to
   `docs/session02_shared_skeleton.md`. Base R throughout, with two exceptions:
   session 5 uses `glmmTMB` (with a `glm` fallback), and session 7 requires the
   `smoltEASE` package because it runs production `SCRAPI2()` on the real data.

## What is not here

Everything cut from the earlier construction-oriented plan is in `BACKLOG.md`,
each item with one line saying why it was cut. Nothing was deleted.

## Ground rules

Each session teaches one idea and one idea only: simulate, build the estimator by
hand, recover against truth, one figure, Locate. A second figure means a second
session. Single-use helpers are inlined; no defensive input checking, since these
are teaching scripts, not package code. No session over a minute to run. No banner
comments. No placeholder headings. No dates — those get added by hand.
