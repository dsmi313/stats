#---------------------------------------------------------
# File:   section01_binomial_window_count-solutions.R
# Part I, Section 1 solutions
#---------------------------------------------------------
# Section 1 solutions ----

library(ggplot2)
library(dplyr)
library(tibble)


#--------------------------------------
# Problem 1a: Simulate one day at the ladder window
section1_problem_1a_fish <- function(a_d_true, r_sample) {
  cat("\n----------------------------------\n")
  cat("Problem 1a: Single day window count and a_d_hat\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  w_obs   <- rbinom(1, size = a_d_true, prob = r_sample)
  a_d_hat <- w_obs / r_sample

  cat("a_d_true =", a_d_true, "  r_sample =", r_sample, "\n")
  cat("Observed window count w =", w_obs, "\n")
  cat("a_d_hat = w / r          =", a_d_hat, "\n")
  invisible(list(w = w_obs, a_d_hat = a_d_hat))
}


# Problem 1b: Build escapeLGD-style `wc` tibble for a season
section1_problem_1b_fish <- function(n_weeks, wc_prop, mean_truth) {
  cat("\n----------------------------------\n")
  cat("Problem 1b: Season-long escapeLGD wc tibble (no GE)\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  wc <- tibble(
    sWeek = seq_len(n_weeks),
    truth = round(runif(n_weeks, mean_truth * 0.5, mean_truth * 1.5)),
    wc    = NA_integer_
  )
  wc$wc <- rbinom(n_weeks, size = wc$truth, prob = wc_prop)

  cat("n_weeks =", n_weeks, "  wc_prop =", wc_prop,
      "  mean_truth =", mean_truth, "\n")
  cat("First rows of wc tibble:\n")
  print(head(wc, 6))
  invisible(wc)
}


# Problem 1c: Vectorized parametric bootstrap of season total
section1_problem_1c_fish <- function(wc, wc_prop, boots) {
  cat("\n----------------------------------\n")
  cat("Problem 1c: Vectorized parametric bootstrap of season total\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  # escapeLGD line 150 reproduced:
  boot_mat <- vapply(seq_len(nrow(wc)),
                     function(i) rbinom(boots, wc$wc[i], wc_prop) / wc_prop,
                     numeric(boots))
  season_totals <- rowSums(boot_mat)
  ci <- quantile(season_totals, c(0.025, 0.975))

  cat("boots =", boots, "  wc_prop =", wc_prop, "\n")
  cat("Point estimate (sum(wc / r)) =", sum(wc$wc / wc_prop), "\n")
  cat("95% bootstrap CI             = [", round(ci[1]), ",",
      round(ci[2]), "]\n")
  invisible(list(season_totals = season_totals, ci = ci))
}


# Problem 1d: Replicate to verify a_d_hat = w/r is unbiased
section1_problem_1d_fish <- function(a_d_true, r_sample, nreps) {
  cat("\n----------------------------------\n")
  cat("Problem 1d: Replicated unbiasedness check\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  sims <- vapply(seq_len(nreps),
                 function(i) rbinom(1, a_d_true, r_sample) / r_sample,
                 numeric(1))
  cat("a_d_true =", a_d_true, "  r_sample =", r_sample,
      "  nreps =", nreps, "\n")
  cat("mean(a_d_hat) =", mean(sims),
      "   sd(a_d_hat) =", sd(sims), "\n")

  p <- ggplot(tibble(a_d_hat = sims), aes(a_d_hat)) +
    geom_histogram(bins = 30, fill = "steelblue", colour = "white") +
    geom_vline(xintercept = a_d_true, colour = "firebrick", linewidth = 1) +
    labs(title = "Section 1 - Window-count estimator a_d_hat = w / r",
         x = expression(hat(a)[d]), y = "Replicates")
  invisible(sims)
}


# Narrative walkthrough ----
#
# The wrapper functions above are the correct answers to the problem set.
# What follows explains why they work and shows three ways the model breaks
# when assumptions fail. Run this after you have completed the problem set.

set.seed(2026)

a_d_true   <- 400
r_sample   <- 5 / 6
mean_truth <- 400
n_weeks    <- 20
boots      <- 2000
nreps      <- 5000

# Why a_d_hat = w / r is unbiased ----
#
# Each fish is a Bernoulli trial: seen with probability r_sample, missed
# otherwise. The observed window count W is the sum of those trials, so
# W ~ Binomial(a_d_true, r_sample).
#
# E(W) = a_d_true * r_sample
#
# Solving for a_d_true: a_d_hat = W / r_sample.
# Because E(W / r_sample) = a_d_true * r_sample / r_sample = a_d_true,
# the estimator is unbiased — it is centered on the true value.
#
# The simulation in 1d confirms this numerically. The mean of a_d_hat
# across nreps replicates should be very close to a_d_true = 400.

sims <- vapply(seq_len(nreps),
               function(i) rbinom(1, a_d_true, r_sample) / r_sample,
               numeric(1))

cat("Unbiasedness check\n")
cat("mean(a_d_hat) =", round(mean(sims), 2), "  true =", a_d_true, "\n")
cat("sd(a_d_hat)   =", round(sd(sims), 2),
    "  theoretical =",
    round(sqrt(a_d_true * r_sample * (1 - r_sample)) / r_sample, 2), "\n\n")

# Failure 1 — wrong r ----
#
# The estimator is only unbiased if r_sample matches the true detection rate.
# If the observer is on TikTok, or visibility is poor, the true rate is lower.
#
# E(a_d_hat) = a_d_true * r_true / r_assumed
#
# If r_true < r_assumed, every single estimate is too low by the same factor.
# This is not random noise — it is a systematic undercount that compounds
# across every week of the season.

r_assumed  <- 5 / 6
r_moderate <- 5 / 8   # one hour short of the full window
r_extreme  <- 0.05    # observer on TikTok

N_hat_correct  <- rbinom(nreps, a_d_true, r_assumed)  / r_assumed
N_hat_moderate <- rbinom(nreps, a_d_true, r_moderate) / r_assumed
N_hat_extreme  <- rbinom(nreps, a_d_true, r_extreme)  / r_assumed

cat("Failure 1: wrong assumed r\n")
cat("True a_d_true =", a_d_true, "\n")
cat("correct  (r_true = 5/6): mean =", round(mean(N_hat_correct),  1), "\n")
cat("moderate (r_true = 5/8): mean =", round(mean(N_hat_moderate), 1),
    " -- missing ~", round(a_d_true - mean(N_hat_moderate), 0), "fish/day\n")
cat("TikTok   (r_true = 0.05): mean =", round(mean(N_hat_extreme),  1),
    " -- missing ~", round(a_d_true - mean(N_hat_extreme),  0), "fish/day\n\n")

op <- par(no.readonly = TRUE)
par(mfrow = c(1, 3), mar = c(4, 4, 4, 1))

hist(N_hat_correct,
     breaks = 35, col = "honeydew3", border = "white",
     xlab = "a_d_hat", main = "Correct r = 5/6\nUnbiased")
abline(v = a_d_true, col = "firebrick", lwd = 2)
legend("topleft", legend = "True N", col = "firebrick", lwd = 2, bty = "n")

hist(N_hat_moderate,
     breaks = 35, col = "khaki2", border = "white",
     xlim = c(min(N_hat_moderate) * 0.98, a_d_true * 1.05),
     xlab = "a_d_hat", main = "True r = 5/8, assumed 5/6\nModerate bias")
abline(v = a_d_true,             col = "firebrick", lwd = 2)
abline(v = mean(N_hat_moderate), col = "navy",      lwd = 2, lty = 2)
legend("topleft",
       legend = c("True N", "Mean a_d_hat"),
       col = c("firebrick", "navy"), lwd = 2, lty = c(1, 2), bty = "n")

hist(N_hat_extreme,
     breaks = 20, col = "lightcoral", border = "white",
     xlab = "a_d_hat", main = "True r = 0.05, assumed 5/6\nCatastrophic bias")
abline(v = mean(N_hat_extreme), col = "navy", lwd = 2, lty = 2)
mtext(paste("True N =", a_d_true, "is off this axis to the right"),
      side = 3, line = 0.2, cex = 0.75, col = "firebrick")
legend("topright", legend = "Mean a_d_hat",
       col = "navy", lwd = 2, lty = 2, bty = "n")

par(op)

# Failure 2 — random weekly r, CI too narrow ----
#
# Now r is correct on average but varies week to week — some days sharp,
# some days poor conditions. The point estimate stays approximately unbiased
# because the errors cancel. The CI does not.
#
# EASE builds its bootstrap CI assuming r is fixed. The true variance of W
# when r is random has an extra term EASE never sees:
#
#   Var(W | random r) = N * E(r) * (1 - E(r))   <- what EASE models
#                     + N^2 * Var(r)             <- what EASE misses
#
# The reported CI is too narrow. Managers use that CI to decide whether
# escapement is safely above the goal. When the CI is artificially tight,
# harvests get greenlit that should have been held back.
#
# We demonstrate this by computing empirical coverage: how often does the
# EASE 95% CI actually contain the true season total? It should be 95%.

r_sd      <- 0.08
beta_var  <- r_sd^2
alpha_par <- r_assumed * (r_assumed * (1 - r_assumed) / beta_var - 1)
beta_par  <- (1 - r_assumed) * (r_assumed * (1 - r_assumed) / beta_var - 1)

n_seasons     <- 2000
boots_cov     <- 500
ci_contains_N <- logical(n_seasons)

for (s in seq_len(n_seasons)) {
  r_week     <- rbeta(n_weeks, alpha_par, beta_par)
  W_week     <- rbinom(n_weeks, size = a_d_true, prob = r_week)
  N_hat_week <- round(W_week / r_assumed)
  boot_s     <- vapply(seq_len(boots_cov), function(b) {
    sum(rbinom(n_weeks, size = N_hat_week, prob = r_assumed) / r_assumed)
  }, numeric(1))
  ci_s              <- quantile(boot_s, c(0.025, 0.975))
  true_total        <- a_d_true * n_weeks
  ci_contains_N[s]  <- ci_s[1] <= true_total & true_total <= ci_s[2]
}

cat("Failure 2: random weekly r\n")
cat("Nominal CI coverage:   95%\n")
cat("Empirical CI coverage:", round(mean(ci_contains_N) * 100, 1), "%\n\n")

# Failure 3 — r degrades systematically across the season ----
#
# Early season: good conditions, r near 5/6 (the ceiling — you cannot
# observe more fish than pass during your watch window).
# Late season: fatigue, weather, shorter days push r down toward 0.5.
#
# EASE assumes 5/6 all season. Early weeks are fine. Late weeks are silently
# biased low. No individual week looks alarming. But the undercount
# accumulates into the season total where it matters for management.

r_true_season <- seq(5 / 6, 0.5, length.out = n_weeks)
N_week_true   <- rep(a_d_true, n_weeks)
W_season      <- rbinom(n_weeks, size = N_week_true, prob = r_true_season)
N_hat_wrong   <- W_season / r_assumed

cat("Failure 3: seasonal r degradation (5/6 -> 0.5)\n")
cat("True season total:           ", sum(N_week_true), "\n")
cat("EASE season total (wrong r): ", round(sum(N_hat_wrong), 0), "\n")
cat("Bias:                        ",
    round(sum(N_hat_wrong) - sum(N_week_true), 0), "fish\n\n")

op <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

plot(seq_len(n_weeks), r_true_season,
     type = "b", pch = 16, col = "firebrick", lwd = 2,
     ylim = c(0.4, 1.0),
     xlab = "Statistical week", ylab = "Sampling fraction r",
     main = "True r degrades across season\nEASE assumes constant 5/6")
abline(h = r_assumed, col = "navy", lwd = 2, lty = 2)
legend("topright",
       legend = c("True r", "Assumed r"),
       col = c("firebrick", "navy"), lwd = 2, lty = c(1, 2), bty = "n")

plot(seq_len(n_weeks), cumsum(N_hat_wrong) - cumsum(N_week_true),
     type = "b", pch = 16, col = "firebrick", lwd = 2,
     xlab = "Statistical week",
     ylab = "Cumulative bias (N_hat - N_true)",
     main = "Undercount accumulates\nas season progresses")
abline(h = 0, col = "gray50", lty = 2)

par(op)
