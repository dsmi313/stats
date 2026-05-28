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

Each section is split into MIT's problem/solution/test format:

| Problems file | Solutions file | Test driver |
| --- | --- | --- |
| `section01_binomial_window_count.R` | `section01_binomial_window_count-solutions.R` | `section01_binomial_window_count-test.R` |
| `section02_mle_window_count.R` | `section02_mle_window_count-solutions.R` | `section02_mle_window_count-test.R` |
| `section03_variance_delta_method.R` | `section03_variance_delta_method-solutions.R` | `section03_variance_delta_method-test.R` |
| `section04_parametric_bootstrap.R` | `section04_parametric_bootstrap-solutions.R` | `section04_parametric_bootstrap-test.R` |
| `section05_nonparametric_bootstrap.R` | `section05_nonparametric_bootstrap-solutions.R` | `section05_nonparametric_bootstrap-test.R` |
| `section06_stratification.R` | `section06_stratification-solutions.R` | `section06_stratification-test.R` |
| `section07_multinomial_composition.R` | `section07_multinomial_composition-solutions.R` | `section07_multinomial_composition-test.R` |

Each problems file contains wrapper functions in MIT 18.05 format:

```r
sectionN_problem_Xx_fish <- function(args) {
  cat("\n----------------------------------\n")
  cat("Problem Xx: <title>\n")
  # Arguments described above.
  # Do not change the above code.
  # ********* YOUR CODE HERE ***********
}
```

Problems and what they mirror:

| Section | Problems | Mirrors |
| --- | --- | --- |
| 1 | 1a single-day a_d_hat = w/r; 1b escapeLGD wc tibble; 1c vectorized parametric bootstrap; 1d replicated unbiasedness | `night_fall_reascend_wc_binom.R:130-150` |
| 2 | 2a replicated MLE centering; 2b grid log-likelihood; 2c optim continuous relaxation | `SCRAPI.r:74-126` |
| 3 | 3a layer-1 variance; 3b layer-2 variance; 3c delta method; 3d fallback joint MLE (escapeLGD verbatim); 3e nightFall bootstrap | `fallback_reascend_likelihood.R`; `night_fall_reascend_wc_binom.R:71-110` |
| 4 | 4a bootsmolt daily-loop (verbatim); 4b vectorized escapeLGD; 4c coverage check; 4d plot bootstrap distribution | `SCRAPI.r:139-145`; `night_fall_reascend_wc_binom.R:150` |
| 5 | 5a inverse-SR weighting (Horvitz-Thompson); 5b toy data frames; 5c FishWH resample (verbatim); 5d FishDat resample (verbatim); 5e thetahat_toy; 5f full bootsmolt; 5g Error 1 trigger | `SCRAPI.r:74-126, 128-189` |
| 6 | 6a simulate shifting p_n; 6b Cpattern + pooled vs stratified; 6c Collaps loop (verbatim); 6d Error 3 trigger; 6e parametric bootstrap CI | `SCRAPI.r:217, 254-269` |
| 7 | 7a rmultinom three-stock simulation; 7b multinomial log-likelihood + optim(); 7c joint two-sample log-likelihood (D&H skeleton); 7d nonparametric bootstrap CI | `composition_estimation_utils.R:62-66` |

See `repo_map.md` for the consolidated line-by-line index.

## Section checkpoints (what you should learn)

Use these as a quick self-check before moving to the next section.

### Section 1 checklist
- [ ] Explain why the window-count estimator is \(\hat{a}_d = w/r\).
- [ ] Build a season `wc` tibble with `sWeek`, `truth`, and observed `wc`.
- [ ] Run a vectorized parametric bootstrap and interpret its 95% CI.
- [ ] Demonstrate empirically that \(\hat{a}_d\) is unbiased under repeated simulation.

### Section 2 checklist
- [ ] Write and evaluate the log-likelihood on a parameter grid.
- [ ] Compare grid-search intuition to `optim()` output for the same model.
- [ ] Explain why the continuous-relaxation MLE is centered near truth in simulation.

### Section 3 checklist
- [ ] Separate and compute layer-1 vs layer-2 uncertainty terms.
- [ ] Apply the delta method to propagate variance through a nonlinear estimator.
- [ ] Reproduce the fallback joint-MLE logic used in escapeLGD.
- [ ] Compare analytic variance ideas to bootstrap behavior for night-fall windows.

### Section 4 checklist
- [ ] Implement the bootsmolt daily-loop structure from SCOBI.
- [ ] Implement a vectorized escapeLGD-style bootstrap for season totals.
- [ ] Check bootstrap interval coverage in repeated simulation.
- [ ] Read and explain a bootstrap distribution plot (shape, center, spread).

### Section 5 checklist
- [ ] Explain why summing 1/SR within each Strat × PGrp cell gives an unbiased composition estimate.
- [ ] Build toy FishWH/FishDat-style data frames used by resampling code.
- [ ] Perform nonparametric resampling for FishWH and FishDat (verbatim logic).
- [ ] Compute `thetahat` on resampled data and aggregate bootstrap totals.
- [ ] Trigger and interpret the Section 5 Error 1 path deliberately.

### Section 6 checklist
- [ ] Simulate changing `p_n` over time and diagnose pooled-estimator bias.
- [ ] Construct `Cpattern`-based pooled vs stratified comparisons.
- [ ] Reproduce the `Collaps` loop logic and explain each grouped quantity.
- [ ] Trigger Error 3 and explain why the stratification guard exists.
- [ ] Build a parametric-bootstrap CI under stratified structure.

### Section 7 checklist
- [ ] Simulate a multinomial sample and verify the MLE equals observed proportions.
- [ ] Write the multinomial log-likelihood and confirm `optim()` recovers `counts / sum(counts)`.
- [ ] Explain why the joint log-likelihood for two multinomial samples equals the sum of their individual log-likelihoods.
- [ ] Build a nonparametric bootstrap CI for stock proportions and interpret its width.

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

2. Work through one section. Each prints results to the console.

3. Open the cited source file in parallel and walk through the lines
   listed in the "Repo pointer" header of the section script.

4. When you can read those lines without surprise, move to the next
   section. By the end of Section 7 you should be able to follow
   `SCRAPI()` end-to-end and recognize every helper called inside
   `escapeLGD::HNC_expand()`.

## Running

The workflow mirrors MIT 18.05 studios:

```r
# Option A -- work on the stubs yourself
source("R/PartI/section01_binomial_window_count.R")
source("R/PartI/section01_binomial_window_count-test.R")

# Option B -- run the reference answers
source("R/PartI/section01_binomial_window_count-solutions.R")
source("R/PartI/section01_binomial_window_count-test.R")
```

Required packages: `ggplot2`, `dplyr`, `tibble`. Section 6 also uses `dplyr`.

## Setup helper (recommended)

Use the setup script to install missing dependencies and print package versions:

```r
source("R/PartI/setup_part1.R")
```

If you prefer strict reproducibility, you can also manage a lockfile workflow
with `renv` (for example, `renv::init()` and `renv::snapshot()` after setup).
