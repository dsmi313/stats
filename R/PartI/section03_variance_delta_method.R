# Section 3 — Variance, the delta method, and the joint estimator

# Goals:
# - Simulate nighttime correction uncertainty
# - Compare variance propagation with delta method and bootstrap intuition
# - Visualize right-skew in a_t estimates

library(ggplot2)
library(dplyr)

set.seed(2026)
a_d <- 500
p_n <- 0.25
nreps <- 1000

d_n <- rbinom(nreps, size = a_d, prob = p_n)
p_n_hat <- d_n / a_d

a_t_hat <- a_d / (1 - p_n_hat)

# TODO: compute joint uncertainty when a_d and p_n are both estimated
# TODO: add delta method approximation and compare empirical variance
