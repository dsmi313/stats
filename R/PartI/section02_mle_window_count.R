# Section 2 — MLE for the window-count estimator

# Goals:
# - Simulate repeated binomial samples
# - Show that the MLE for a_d is w / r
# - Plot likelihood as a function of a_d for fixed w and r

library(ggplot2)
library(purrr)

set.seed(2026)
a_d_true <- 500
rr <- 5/6

# Simulate repeated estimates
sim_estimates <- map_dbl(1:1000, ~ rbinom(1, a_d_true, rr) / rr)

# TODO: visualize distribution and compare mean to truth

# TODO: plot binomial likelihood for a fixed observed w
