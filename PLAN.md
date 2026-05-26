# smoltEASE / SCRAPI Learning Plan

This repository is organized around a structured learning and development plan for understanding EASE/SCRAPI estimation using real fisheries examples and R code. The plan follows MIT 18.05-style sessions, with a focus on the statistical foundations, production code, and a Bayesian integrated model for smolt abundance and composition.

## Repository purpose

- Capture the learning plan for each section, with links to materials.
- Host R code and examples for the exercises described in each section.
- Bridge mathematics, simulation, bootstrap, and likelihood-based estimation with real fisheries data and SCRAPI/EASE workflows.
- Document the path from simple binomial estimators to an integrated Bayesian routing and composition model.

## Structure

- `README.md`: Overview and quick navigation.
- `PLAN.md`: Full detailed plan and resources.
- `R/PartI` through `R/PartVIII`: Session folders for evolving R code and analysis.
- `docs/`: Supporting documentation, notes, and reference links.

## Plan outline

### Notation Glossary

Standard symbols used throughout the repository, with EASE-specific names and MIT 18.05 aliases.

### Part I — Building EASE from Scratch

Section 1 — Binomial window count [done]
Section 2 — MLE for the window-count estimator [in progress]
Section 3 — Variance, the delta method, and the joint estimator
Section 4 — Parametric bootstrap and coverage
Section 5 — Nonparametric bootstrap and distributional diagnostic
Section 6 — Stratification: when pooling lies
Section 7 — Multinomial composition

### Part II — The Delomas & Hess Likelihood

Section 8 — PBT tag rates and the naive estimator's bias
Section 9 — D&H likelihood part 1: cell probabilities by hand
Section 10 — D&H likelihood part 2: joint optim and GSI uncertainty

### Part III — SCRAPI Internals

Section 11 — SCRAPI ratio estimator (daily abundance)
Section 12 — The compound bootstrap (bootsmolt)
Section 13 — Reading SCOBI source + lgr2SCRAPI audit
Section 14 — End-to-end run on MY2024

### Part IV — Parallel Implementation

Section 15 — fishCompTools as a parallel implementation

### Part V — Outside Perspectives

Section 16 — escapeLGD internals: reading your own production tool
Section 17 — FSA: broader fisheries bootstrap and stratification idioms
Section 18 — TropFishR: MLE in a different fisheries domain

### Part VI — Bayesian Sidebar and STADEM

Section 19 — Bayesian sidebar: the same problem, different lens
Section 20 — STADEM: the state-space answer at your dam

### Part VII — Synthesis

Section 21 — Year-one synthesis memo

### Part VIII — SCRAPI Innovation: Integrated Bayesian Routing and Composition Model

Section 22 — The integrated model: architecture

## Next steps

1. Add R scripts for each section into the corresponding `R/Part*` folder.
2. Populate examples and simulation code for Sections 1–7.
3. Add audit notes and source-review summaries for Parts II–III.
4. Build and validate the `smoltEASE` integrated model in Part VIII.

## Useful links

The plan includes MIT 18.05 resources, Delomas & Hess (2021), Hance et al. (2024), CSS 2023 Chapter 8, SCOBI, escapeLGD, smoltEASE, and other fisheries modeling references.

---

> For a finer-grained implementation, create one R script per section and add data or notebooks as the analysis progresses.
