# Section 1 — Binomial window count

# Goals:
# - Simulate binomial window counts
# - Visualize sampling variability and estimator centering
# - Relate window count estimation to EASE/SCRAPI notation

library(ggplot2)
library(dplyr)

# Example: generate a single season window count
set.seed(2026)
a_d_true <- 500
rr <- 5/6
w <- rbinom(1, size = a_d_true, prob = rr)

# Estimate daytime escapement from window count
a_d_hat <- w / rr

cat("Observed window count:", w, "\n")
cat("Estimated daytime escapement:", a_d_hat, "\n")

# TODO: extend with replicate() simulation and histogram
