# Part I - Building EASE from Scratch

Seven sections plus a fish-themed port of every MIT 18.05 R studio. By the
end of Part I you can simulate every binomial estimator EASE uses, build
CIs for them by parametric and nonparametric bootstrap, and explain why
three strata schemes exist.

## Sections (the learning-plan spine)

Each section file is self-contained and saves its figures under
`docs/figures/PartI/`. Run from the repository root.

| File | Topic | Maps to PLAN.md |
| --- | --- | --- |
| `section01_binomial_window_count.R` | window count W ~ Binomial(a_d, r); a_d_hat = w/r | Section 1 |
| `section02_mle_window_count.R` | binomial MLE, likelihood surface, probability != likelihood | Section 2 |
| `section03_variance_delta_method.R` | variance propagation, delta method, joint estimator | Section 3 |
| `section04_parametric_bootstrap.R` | parametric bootstrap CI + coverage check | Section 4 |
| `section05_nonparametric_bootstrap.R` | nonparametric + stratified resample + chi-square diagnostic | Section 5 |
| `section06_stratification.R` | pooled vs stratified estimators; ">=100 wild fish" rule | Section 6 |
| `section07_multinomial_composition.R` | multinomial MLE; D&H skeleton (tagged + untagged) | Section 7 |

## MIT studio ports

`mit_studios/` contains fish-themed translations of all 10 MIT 18.05
(Spring 2022) R studios. Each MIT studio wrapper function
(`studioN_problem_Ka()`) is renamed `studioN_problem_Ka_fish()` and
operates on the EASE/SCRAPI domain. See `mit_studios/README.md` for the
full studio-to-section map.

Highlights:

- Studio 1 -> PIT-tag collision probability (birthday paradox in fish form)
- Studio 4 -> two trap supervisors sharing shifts (covariance from
  shared signal)
- Studio 5 -> sequential GSI markers updating posterior on stock identity
  (Bayesian version of the iterative replacement in PLAN.md Section 10)
- Studio 6 -> sonar-tagged fish position from Cauchy bank-projections
- Studio 10 -> bootstrap CI coverage on log-normal smolt lengths (the
  realistic stress test referenced in PLAN.md Sections 5 and 14)

## Running

```r
# from the repository root
source("R/PartI/section01_binomial_window_count.R")
source("R/PartI/mit_studios/studio01_pit_tag_collisions.R")
```

Required packages: `ggplot2`, `dplyr`, `purrr`, `tidyr`. Studio 10
optionally uses `matrixStats` for column-wise statistics (falls back to
`apply` if not installed).
