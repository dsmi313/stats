# Section 7 — Multinomial composition

# Goals:
# - Simulate multinomial stock composition data
# - Verify that observed proportions are the multinomial MLE
# - Build and optimize a multinomial log-likelihood

library(ggplot2)
library(dplyr)

set.seed(2026)
prob <- c(0.6, 0.3, 0.1)
counts <- rmultinom(1, size = 100, prob = prob)
observed_props <- counts / sum(counts)

cat("Observed proportions:", observed_props, "\n")

# TODO: implement multinomial log-likelihood and optim() verification
# TODO: add joint multinomial skeleton for tagged vs untagged samples
