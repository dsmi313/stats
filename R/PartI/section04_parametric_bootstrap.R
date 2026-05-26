# Section 4 - Parametric bootstrap and coverage
# ---------------------------------------------------------------
# MIT R Studio 9 (Simulating CIs) + R Studio 10 (Bootstrap CIs) cast
# into LWG units. A 95% CI is only honest if it covers the truth ~95%
# of the time across simulated seasons.
#
# Reading: PLAN.md Section 4, MIT Class 22 (Confidence Intervals),
# Class 24 (Bootstrap CIs), R Studio 9 and 10.

library(ggplot2)
library(dplyr)
library(purrr)

set.seed(2026)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# --- 1. Parametric bootstrap for a_d ---------------------------------------
boot_a_d <- function(w, r, B = 10000) {
  a_d_hat  <- w / r
  w_star   <- rbinom(B, size = round(a_d_hat), prob = r)
  a_d_star <- w_star / r
  list(estimate = a_d_hat,
       ci       = quantile(a_d_star, c(0.025, 0.975)),
       draws    = a_d_star)
}

a_d_true <- 500L
r_sample <- 5/6
w_obs    <- rbinom(1, a_d_true, r_sample)
boot_d   <- boot_a_d(w_obs, r_sample, B = 10000)

cat("Parametric bootstrap for a_d, one season:\n")
cat("  a_d_hat =", boot_d$estimate, "\n")
cat("  95% CI  =", round(boot_d$ci, 1), "\n\n")

# --- 2. Joint parametric bootstrap for a_t ---------------------------------
# Re-draw both w and d_n from their fitted Binomials, then compute a_t*.
boot_a_t <- function(w, r, d_a, d_n, B = 10000) {
  a_d_hat  <- w / r
  p_n_hat  <- d_n / d_a
  w_star   <- rbinom(B, size = round(a_d_hat), prob = r)
  d_n_star <- rbinom(B, size = d_a,            prob = p_n_hat)
  a_d_star <- w_star  / r
  p_n_star <- d_n_star / d_a
  a_t_star <- a_d_star / (1 - p_n_star)
  list(estimate = a_d_hat / (1 - p_n_hat),
       ci       = quantile(a_t_star, c(0.025, 0.975)),
       draws    = a_t_star)
}

d_a      <- 200L
p_n_true <- 0.25
d_n_obs  <- rbinom(1, d_a, p_n_true)
boot_t   <- boot_a_t(w_obs, r_sample, d_a, d_n_obs, B = 10000)

cat("Joint parametric bootstrap for a_t, one season:\n")
cat("  a_t_hat =", round(boot_t$estimate, 1), "\n")
cat("  95% CI  =", round(boot_t$ci, 1), "\n\n")

# --- 3. Coverage check across many simulated seasons ----------------------
coverage_one <- function(seed, a_d_true, r, d_a, p_n_true, B = 2000) {
  set.seed(seed)
  w   <- rbinom(1, a_d_true, r)
  d_n <- rbinom(1, d_a, p_n_true)
  ci  <- boot_a_t(w, r, d_a, d_n, B = B)$ci
  truth <- a_d_true / (1 - p_n_true)
  (truth >= ci[1]) && (truth <= ci[2])
}

nseasons <- 500L
hits <- map_lgl(seq_len(nseasons),
                ~ coverage_one(seed = .x,
                               a_d_true = a_d_true, r = r_sample,
                               d_a = d_a, p_n_true = p_n_true,
                               B = 1500))
cat("Coverage of the 95% bootstrap CI for a_t across",
    nseasons, "seasons =", mean(hits), "\n\n")

# --- 4. Plot the bootstrap distribution -----------------------------------
p_boot <- ggplot(tibble(a_t = boot_t$draws), aes(a_t)) +
  geom_histogram(bins = 50, fill = "steelblue", colour = "white") +
  geom_vline(xintercept = boot_t$ci,
             colour = "firebrick", linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = a_d_true / (1 - p_n_true),
             colour = "black", linewidth = 1) +
  labs(title = "Section 4 - Parametric bootstrap distribution of a_t",
       subtitle = "Dashed = 95% CI, solid = truth",
       x = expression(hat(a)[t]),
       y = "Bootstrap draws")
ggsave(file.path(plots_dir, "section04_boot_a_t.png"), p_boot,
       width = 6, height = 4, dpi = 150)

# --- 5. CI width as PIT sample size grows ----------------------------------
# More PIT-tagged fish ==> narrower CI on a_t. Confirms the bias-variance
# story before Section 6 asks where to draw the >=100 wild-fish line.
sample_sweep <- map_dfr(c(50, 100, 200, 400, 800), function(da) {
  draws <- boot_a_t(w_obs, r_sample, da,
                    d_n = rbinom(1, da, p_n_true), B = 1500)$draws
  tibble(d_a = da,
         lower = unname(quantile(draws, 0.025)),
         upper = unname(quantile(draws, 0.975)),
         width = upper - lower)
})

cat("CI width as PIT sample size grows:\n")
print(sample_sweep)
cat("\n")

p_sweep <- ggplot(sample_sweep, aes(d_a, width)) +
  geom_line(colour = "steelblue") +
  geom_point(size = 2) +
  labs(title = "Section 4 - 95% CI width shrinks with PIT sample size",
       x = "d_a (PIT-tagged adults)", y = "CI width on a_t")
ggsave(file.path(plots_dir, "section04_ci_width.png"), p_sweep,
       width = 6, height = 4, dpi = 150)

# Section 4 payoff ----------------------------------------------------------
# This is the inner loop of bootsmolt(). Section 12 wraps the same pattern
# around composition resampling and a denominator that itself varies daily.
