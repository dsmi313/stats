# Section 3 - Variance, the delta method, and the joint estimator
# ---------------------------------------------------------------
# The nighttime correction has two layers of uncertainty:
#   1. window count w (Section 1)
#   2. nighttime proportion p_n estimated from PIT-tag day/night counts
# Corrected total escapement is a_t = a_d / (1 - p_n).
# This is a nonlinear transform with a noisy denominator. Its variance
# is given exactly by simulation and approximately by the delta method.
#
# Reading: PLAN.md Section 3, MIT Class 5a (Variance of Discrete RVs),
# Class 4b, Class 5d, R Studio 8.

library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)

set.seed(2026)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# --- 1. Single source of noise: p_n only -----------------------------------
d_a       <- 200L     # PIT-tagged adults available to estimate p_n
p_n_true  <- 0.25     # truth: 25% of adults pass at night
a_d_fixed <- 500L     # daytime escapement held fixed for layer 1
nreps     <- 10000L

d_n_sim <- rbinom(nreps, size = d_a, prob = p_n_true)
p_n_hat <- d_n_sim / d_a
a_t_one <- a_d_fixed / (1 - p_n_hat)

cat("Layer 1 (only p_n noisy):\n")
cat("  mean(a_t) =", round(mean(a_t_one), 1), "\n")
cat("  sd(a_t)   =", round(sd(a_t_one), 1), "\n\n")

# --- 2. Two sources: a_d and p_n both noisy --------------------------------
a_d_true <- 500L
r_sample <- 5/6

w_sim   <- rbinom(nreps, size = a_d_true, prob = r_sample)
a_d_sim <- w_sim / r_sample
a_t_two <- a_d_sim / (1 - p_n_hat)

cat("Layer 2 (a_d and p_n both noisy):\n")
cat("  mean(a_t) =", round(mean(a_t_two), 1), "\n")
cat("  sd(a_t)   =", round(sd(a_t_two), 1), "\n\n")

# --- 3. Delta-method variance approximation --------------------------------
# X = a_d_hat, Y = p_n_hat, g(X, Y) = X / (1 - Y).
# Assuming independence:
#   Var(g) ~ (dg/dX)^2 * Var(X) + (dg/dY)^2 * Var(Y)
var_a_d <- a_d_true * (1 - r_sample) / r_sample^2     # Var(W/r)
var_p_n <- p_n_true * (1 - p_n_true) / d_a            # Var(d_n/d_a)
g_dX <- 1 / (1 - p_n_true)
g_dY <- a_d_true / (1 - p_n_true)^2
var_a_t_delta <- g_dX^2 * var_a_d + g_dY^2 * var_p_n

cat("Delta-method check:\n")
cat("  Var(a_d_hat)             =", round(var_a_d, 2), "\n")
cat("  Var(p_n_hat)             =", round(var_p_n, 5), "\n")
cat("  Var(a_t_hat) [delta]     =", round(var_a_t_delta, 2), "\n")
cat("  sd  (a_t_hat) [delta]    =", round(sqrt(var_a_t_delta), 2), "\n")
cat("  sd  (a_t_hat) [empir]    =", round(sd(a_t_two), 2), "\n\n")

# --- 4. Fallback correction ------------------------------------------------
p_f       <- 0.05
a_t_corr  <- a_t_two * (1 - p_f)
cat("Layer 3 (with fallback p_f = 0.05):\n")
cat("  mean(a_t_corr) =", round(mean(a_t_corr), 1), "\n")
cat("  sd(a_t_corr)   =", round(sd(a_t_corr), 1), "\n\n")

# --- 5. Plot the right-skew of a_t across layers ---------------------------
skew_tbl <- bind_rows(
  tibble(layer = "1: only p_n noisy",        a_t = a_t_one),
  tibble(layer = "2: a_d and p_n noisy",     a_t = a_t_two),
  tibble(layer = "3: + fallback correction", a_t = a_t_corr)
)

p_skew <- ggplot(skew_tbl, aes(a_t, fill = layer)) +
  geom_histogram(bins = 60, alpha = 0.7) +
  facet_wrap(~ layer, ncol = 1, scales = "free_y") +
  geom_vline(xintercept = a_d_true / (1 - p_n_true),
             colour = "black", linetype = "dashed") +
  labs(title = "Section 3 - Layered uncertainty in a_t",
       subtitle = "Dashed line = true expected a_t before fallback",
       x = expression(hat(a)[t]),
       y = "Replicates") +
  theme(legend.position = "none")
ggsave(file.path(plots_dir, "section03_layered_uncertainty.png"), p_skew,
       width = 6, height = 7, dpi = 150)

# --- 6. Coverage of a normal CI vs the truth ------------------------------
# Closed-form variance plus normal CI: how often does it cover?
nseasons <- 500L
hits <- map_lgl(seq_len(nseasons), function(seed) {
  set.seed(seed)
  w   <- rbinom(1, a_d_true, r_sample)
  d_n <- rbinom(1, d_a, p_n_true)
  a_d_h <- w / r_sample
  p_n_h <- d_n / d_a
  if (p_n_h >= 1) return(FALSE)
  var_d  <- a_d_h * (1 - r_sample) / r_sample
  var_p  <- p_n_h * (1 - p_n_h) / d_a
  est    <- a_d_h / (1 - p_n_h)
  se_est <- sqrt((1 / (1 - p_n_h))^2 * var_d +
                 (a_d_h / (1 - p_n_h)^2)^2 * var_p)
  truth  <- a_d_true / (1 - p_n_true)
  abs(est - truth) <= 1.96 * se_est
})
cat("Normal-CI coverage from delta-method SE across",
    nseasons, "seasons =", mean(hits), "\n")
cat("(Section 4 replaces this with bootstrap CIs that handle the skew.)\n")

# Section 3 payoff ----------------------------------------------------------
# The right-skew in layer 2 is exactly why bootsmolt() does not use a
# normal-approximation CI; it samples a_t directly from these distributions
# in Section 12. The delta method is the closed-form benchmark to compare.
