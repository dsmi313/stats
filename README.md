# stats

A structured learning and development repository for SCRAPI, EASE, and smolt abundance estimation.

## Goal

This repository captures a planned pathway from MIT 18.05-style statistical foundations through real fisheries examples, SCRAPI/EASE internals, bootstrap and likelihood estimation, and a Bayesian integrated model.

## What is included

- `PLAN.md`: the full detailed learning plan, session structure, and resources.
- `R/PartI` .. `R/PartVIII`: folders for R code, simulations, and analysis for each major section.
- `docs/`: documentation, notes, and reference material.

## How to use this repo

1. Read `PLAN.md` to understand the overall workflow and section goals.
2. Add or expand R scripts inside the appropriate `R/Part*` folders.
3. Use `docs/` for notes, diagnostics, and audit write-ups.

## Suggested next files

- `R/PartI/section01_binomial_window_count.R`
- `R/PartI/section02_mle_window_count.R`
- `R/PartII/section08_pbt_bias.R`
- `R/PartIII/section11_scrapi_ratio_estimator.R`
- `R/PartVIII/section22_integrated_bayesian_model.R`

## Notes

This repository is intended as both a learning notebook and a code base for developing the `smoltEASE` integrated model in a reproducible way.
