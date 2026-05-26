# MIT 18.05 Studios - Fish Ports

Fish-themed translations of each of the ten MIT 18.05 (Spring 2022) R
Studios. The originals live at
https://ocw.mit.edu/courses/18-05-introduction-to-probability-and-statistics-spring-2022/pages/studio-resources/.
Each fish port keeps MIT's wrapper-function structure
(`studioN_problem_K_fish()`) and the same statistical content, but swaps
generic toy scenarios (coin tosses, casino bets, lighthouses, polls) for
the EASE/SCRAPI domain (window counts, PIT detections, GSI markers,
sonar-tag locations, smolt lengths).

## Studio map

| MIT studio | Fish file | EASE/SCRAPI hook |
| --- | --- | --- |
| 1. Birthday matches | `studio01_pit_tag_collisions.R` | PIT-tag bin collisions in subsampling; GSI baseline-drop intuition |
| 2. Binomial distributions | `studio02_binomial_window_count.R` | Ladder window count; payoff = marked-fish sample |
| 3. Histograms | `studio03_histograms_passage_times.R` | Inter-arrival times at the ladder; CLT for averaged delays |
| 4. Covariance and correlation | `studio04_covariance_supervisors.R` | Two trap supervisors sharing shifts; daily detection covariance |
| 5. Discrete Bayesian updating | `studio05_discrete_bayes_stock_id.R` | Posterior on stock identity from sequential GSI markers (Section 10) |
| 6. Discretized continuous Bayes | `studio06_continuous_bayes_antenna_location.R` | Sonar-tag bank-projection -> Cauchy update on fish position |
| 7. Significance testing | `studio07_significance_detection.R` | Testing PIT-detector rate; alpha vs P(H0 | reject) |
| 8. NHST simulations | `studio08_nhst_fish_data.R` | F-stat, z-test, chi-square, ANOVA on stock-vs-stock length data |
| 9. Simulating CIs | `studio09_simulating_cis_lengths.R` | z- and t-CIs for mean smolt length; coverage + Bayesian posterior |
| 10. Bootstrap CIs | `studio10_bootstrap_cis_lengths.R` | Bootstrap CIs on Normal vs Log-Normal length data |

## How the studios feed the seven Part I sections

| Section | Drawn-from studios |
| --- | --- |
| Section 1 - Binomial window count | Studio 2 (primary), Studio 3 (visualization tips) |
| Section 2 - MLE for the window-count estimator | Studio 2, Studio 4 (CLT framing) |
| Section 3 - Variance and the delta method | Studio 3, Studio 4 (variance/covariance machinery) |
| Section 4 - Parametric bootstrap and coverage | Studio 9, Studio 10 |
| Section 5 - Nonparametric bootstrap + diagnostic | Studio 10 (primary), Studio 8 (chi-square diagnostic) |
| Section 6 - Stratification | Studio 9 (coverage), Studio 8 (ANOVA / between-group variance) |
| Section 7 - Multinomial composition | Studio 5 (discrete Bayes on stock id), Studio 7 (rejection regions for composition tests) |

Studios 5, 6, 7 are also the structural primers for Part VI (Bayesian
sidebar) and Section 19 of the learning plan -- they are placed here so
the entire MIT corpus is available without duplicating later in the repo.

## Running the studios

From the repository root:

```r
source("R/PartI/mit_studios/studio01_pit_tag_collisions.R")
source("R/PartI/mit_studios/studio02_binomial_window_count.R")
# ...
```

All figures land under `docs/figures/PartI/`. Each studio file is
self-contained; it sets its own seed and creates the figures directory if
it does not exist.
