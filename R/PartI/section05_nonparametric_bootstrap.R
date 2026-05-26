# Section 5 — Nonparametric bootstrap and distributional diagnostic

# Goals:
# - Resample fish labels nonparametrically
# - Compare parametric and nonparametric CIs
# - Build a diagnostic for overdispersion

library(ggplot2)
library(dplyr)

set.seed(2026)
fish <- tibble(id = 1:30, stock = sample(c("A", "B"), size = 30, replace = TRUE, prob = c(0.6, 0.4)))

# TODO: implement weighted nonparametric resampling and CI construction
# TODO: build a simple overdispersion diagnostic for binomial counts
