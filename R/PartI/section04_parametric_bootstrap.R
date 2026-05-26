# Section 4 — Parametric bootstrap and coverage

# Goals:
# - Build a parametric bootstrap CI for a_d and a_t
# - Check coverage with repeated simulated seasons
# - Document the bootstrap workflow

library(ggplot2)
library(dplyr)

set.seed(2026)
a_d_true <- 500
rr <- 5/6
w_obs <- rbinom(1, size = a_d_true, prob = rr)
a_d_hat <- w_obs / rr

# Parametric bootstrap for a_d
B <- 10000
w_star <- rbinom(B, size = round(a_d_hat), prob = rr)
a_d_star <- w_star / rr
ci <- quantile(a_d_star, c(0.025, 0.975))

cat("95% CI for a_d:", ci, "\n")

# TODO: extend to joint bootstrap for a_t and coverage assessment
