# Part I - Building EASE from Scratch (targeted at SCOBI + escapeLGD)

Seven sections that build the data structures and statistical machinery
used in two production repositories:

- https://github.com/mackerman44/SCOBI  (SCRAPI, thetahat, bootsmolt)
- https://github.com/delomast/escapeLGD (PBT_log_likelihood, fallback,
  nighttime/window-count binomials)

Every section script keeps the production variable names and reproduces
key code blocks verbatim wherever possible. After Part I you can open
either repo and read it without translating notation in your head.

## Sections

| File | Topic | Mirrors |
| --- | --- | --- |
| `section01_binomial_window_count.R` | window count + SCRAPI `pass`/`passdata` table + escapeLGD `wc` tibble | SCRAPI.r:196-247; night_fall_reascend_wc_binom.R:130-150 |
| `section02_mle_window_count.R` | binomial MLE + inverse-SR weighting (mApply 1/SR) | SCRAPI.r:74-126 |
| `section03_variance_delta_method.R` | delta method + verbatim `fallback_log_likelihood`, `optimllh`, `gradient_fallback_log_likelihood`, optim call | fallback_reascend_likelihood.R (full); night_fall_reascend_wc_binom.R:71-75, 110 |
| `section04_parametric_bootstrap.R` | bootsmolt daily-count layer (line-by-line); vectorized escapeLGD equivalent | SCRAPI.r:139-145; night_fall_reascend_wc_binom.R:145-151 |
| `section05_nonparametric_bootstrap.R` | weighted FishWH/FishDat resample, `thetahat_toy()`, Error 1 trigger | SCRAPI.r:128-189 |
| `section06_stratification.R` | Cpattern + Collaps assignment loop + Error 3 trigger | SCRAPI.r:217, 254-269 |
| `section07_multinomial_composition.R` | `softMax`, `PBT_log_likelihood`, `PBT_optimllh`, `PBT_expand_calc_MLE`, `PBT_expand_calc` (accounting), `PBT_breakdown` (bootstrap) | composition_estimation_utils.R:15-110, 440-486 |

See `repo_map.md` for the consolidated line-by-line index.

## MIT studio ports (supplemental)

`mit_studios/` contains fish-themed translations of every MIT 18.05 R
studio (1 through 10). They are pedagogical companions to the seven
sections above and not the primary path. See `mit_studios/README.md` for
the studio-to-section map.

## Recommended workflow

1. Clone the production repos locally:

   ```bash
   git clone https://github.com/mackerman44/SCOBI ../SCOBI
   git clone https://github.com/delomast/escapeLGD ../escapeLGD
   ```

2. Work through one section. Each prints results to the console and
   saves figures under `docs/figures/PartI/`.

3. Open the cited source file in parallel and walk through the lines
   listed in the "Repo pointer" header of the section script.

4. When you can read those lines without surprise, move to the next
   section. By the end of Section 7 you should be able to follow
   `SCRAPI()` end-to-end and recognize every helper called inside
   `escapeLGD::HNC_expand()`.

## Running

```r
# from the repository root
source("R/PartI/section01_binomial_window_count.R")
source("R/PartI/section02_mle_window_count.R")
# ...
```

Required packages: `ggplot2`, `dplyr`, `purrr`, `tibble`, `tidyr`.
