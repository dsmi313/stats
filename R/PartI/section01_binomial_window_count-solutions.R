# File: section01_binomial_window_count.R
# Part I, Section 1 — Binomial window count (adult ladder passage)
#
# The problem:
#   We do not watch the counting window 24 hours a day. We observe a known
#   fraction r of the passage window. Any adult passing outside that window
#   is missed.
#
# The model:
#   Each fish is an independent Bernoulli trial: detected (1) or missed (0).
#   The observed window count W is the sum of those trials, so
#   W ~ Binomial(a_d_true, r_sample), with E(W) = a_d_true * r_sample.
#
# The estimator:
#   Solving E(W) = a_d_true * r_sample for a_d_true gives
#   a_d_hat = W / r_sample.
#   That is the core of every window-count expansion in EASE.
#
# Three ways the model breaks (covered after the problem set walkthrough):
#   1. You assume the wrong r. N_hat is systematically biased — same
#      direction every single day.
#   2. r varies randomly day to day. The CI is too narrow. You look
#      precise when you are not.
#   3. r degrades across the season. The undercount accumulates silently
#      into the season total.
#
# Problem set map:
#   1a -> single day simulation
#   1b -> season-long wc tibble
#   1c -> vectorized parametric bootstrap
#   1d -> replicated unbiasedness check
#
# EASE pointer (escapeLGD/R/night_fall_reascend_wc_binom.R):
#   line 143: wc_binom <- list(wc %>% mutate(wc = round(wc / wc_prop)))
#   line 150: wc_binom[[2]][,i] <- rbinom(boots, wc[i], wc_prop) / wc_prop
#
# GE note:
#   This script is adults only. Guidance efficiency is a smolt concept
#   and does not enter until Part III.

library(tibble)
library(dplyr)

set.seed(2026)

# shared parameters — match what the problem set passes as arguments
a_d_true  <- 400     # true number of adults passing the ladder
r_sample  <- 5 / 6   # sampling fraction (fraction of window observed)
mean_truth <- 400    # mean true weekly passage for the season tibble
n_weeks   <- 20      # number of statistical weeks in the season
boots     <- 2000    # bootstrap iterations

# Problem 1a — single day: fish-level Bernoulli trials ----
#
# Each of the a_d_true adults either passes during the watched window or not.
# rbinom(1, a_d_true, r_sample) collapses all fish-level trials into one draw.
# That is mathematically the same as sum(rbinom(a_d_true, 1, r_sample))
# but much faster.

w         <- rbinom(1, size = a_d_true, prob = r_sample)  # observed count
a_d_hat   <- w / r_sample                                 # EASE point estimate

cat("Problem 1a\n")
cat("a_d_true =", a_d_true, "  r_sample =", round(r_sample, 4), "\n")
cat("Observed window count w =", w, "\n")
cat("a_d_hat = w / r         =", round(a_d_hat, 2), "\n\n")

# Problem 1b — season-long wc tibble ----
#
# EASE works week by week across the season. Each week has:
#   sWeek: statistical week index
#   truth: true adults passing that week (unobserved)
#   wc:    observed window count ~ Binomial(truth, r_sample)
#
# truth is drawn from a uniform band around mean_truth so weeks vary
# realistically. This mirrors the structure escapeLGD expects.

wc <- tibble(
  sWeek = seq_len(n_weeks),
  truth = round(runif(n_weeks, mean_truth * 0.5, mean_truth * 1.5)),
  wc    = NA_integer_
)
wc$wc <- rbinom(n_weeks, size = wc$truth, prob = r_sample)

cat("Problem 1b\n")
cat("n_weeks =", n_weeks, "  wc_prop =", round(r_sample, 4),
    "  mean_truth =", mean_truth, "\n")
print(head(wc, 6))
cat("\n")

# Problem 1c — vectorized parametric bootstrap of the season total ----
#
# For each bootstrap draw b, resample every week independently:
#   rbinom(1, wc[i], r_sample) / r_sample
# This treats the observed count wc[i] as if it were the true count and
# re-applies the sampling + expansion. Summing across weeks gives one
# bootstrap season total. The 2.5% and 97.5% quantiles are the 95% CI.
#
# vapply applies this across all weeks in one call, returning a boots x n_weeks
# matrix. rowSums collapses each bootstrap draw to a season total.

boot_mat <- vapply(
  seq_len(nrow(wc)),
  function(i) rbinom(boots, wc$wc[i], r_sample) / r_sample,
  numeric(boots)
)
season_totals <- rowSums(boot_mat)   # one season total per bootstrap draw
ci            <- quantile(season_totals, c(0.025, 0.975))

cat("Problem 1c\n")
cat("boots =", boots, "  wc_prop =", round(r_sample, 4), "\n")
cat("Point estimate (sum(wc / r)) =", round(sum(wc$wc / r_sample), 0), "\n")
cat("True season total            =", sum(wc$truth), "\n")
cat("95% bootstrap CI             = [", round(ci[1], 0), ",",
    round(ci[2], 0), "]\n\n")

# Problem 1d — replicated unbiasedness check ----
#
# If a_d_hat = w / r_sample is unbiased, then E(a_d_hat) = a_d_true.
# We verify this by simulating nreps independent days and computing
# the mean and sd of a_d_hat across those replicates.
# The histogram should be centered on the red line (true value).

nreps <- 5000
sims  <- vapply(
  seq_len(nreps),
  function(i) rbinom(1, a_d_true, r_sample) / r_sample,
  numeric(1)
)

cat("Problem 1d\n")
cat("a_d_true =", a_d_true, "  r_sample =", round(r_sample, 4),
    "  nreps =", nreps, "\n")
cat("mean(a_d_hat) =", round(mean(sims), 2),
    "   sd(a_d_hat) =", round(sd(sims), 2), "\n")
cat("Theoretical E(a_d_hat) =", a_d_true, "\n")
cat("Theoretical sd         =",
    round(sqrt(a_d_true * r_sample * (1 - r_sample)) / r_sample, 2), "\n\n")

op <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

hist(sims,
     breaks = 35, col = "lightsteelblue", border = "white",
     xlab = expression(hat(a)[d]),
     main = "Problem 1d: a_d_hat = w / r\nshould be centered on true N")
abline(v = a_d_true, col = "firebrick", lwd = 2)
legend("topright", legend = "True a_d", col = "firebrick", lwd = 2, bty = "n")

# also show the bootstrap season totals from 1c for comparison
hist(season_totals,
     breaks = 35, col = "honeydew3", border = "white",
     xlab = "Bootstrap season total",
     main = "Problem 1c: bootstrap CI\naround season total")
abline(v = sum(wc$truth),        col = "firebrick", lwd = 2)
abline(v = sum(wc$wc / r_sample), col = "navy",     lwd = 2, lty = 2)
abline(v = ci,                   col = "gray40",    lwd = 1, lty = 3)
legend("topright",
       legend = c("True total", "Point estimate", "95% CI bounds"),
       col = c("firebrick", "navy", "gray40"),
       lwd = c(2, 2, 1), lty = c(1, 2, 3), bty = "n")

par(op)

# Extension: three failure modes ----
#
# The problems above assumed r_sample is correct. The following sections
# break that assumption in three different ways. Read these after completing
# the problem set.

# Failure 1 — wrong r, fixed ----
#
# EASE assumes r_sample = 5/6. If the true detection rate is lower —
# because the observer drifted off task or conditions were bad —
# every single estimate is biased in the same direction.
#
# E(a_d_hat) = a_d_true * r_true / r_assumed
#
# If r_true < r_assumed, a_d_hat is too low. Every day. The bias compounds
# across the season into a serious undercount of escapement.

r_assumed  <- 5 / 6   # what EASE thinks
r_extreme  <- 0.05    # TikTok observer
r_moderate <- 5 / 8   # one hour short of the full window

W_correct  <- rbinom(nreps, size = a_d_true, prob = r_assumed)
W_moderate <- rbinom(nreps, size = a_d_true, prob = r_moderate)
W_extreme  <- rbinom(nreps, size = a_d_true, prob = r_extreme)

N_hat_correct  <- W_correct  / r_assumed
N_hat_moderate <- W_moderate / r_assumed
N_hat_extreme  <- W_extreme  / r_assumed

cat("Failure 1: wrong assumed r\n")
cat("True a_d_true =", a_d_true, "\n")
cat("Mean a_d_hat (r_true = 5/6, correct):  ",
    round(mean(N_hat_correct), 1), "\n")
cat("Mean a_d_hat (r_true = 5/8, moderate): ",
    round(mean(N_hat_moderate), 1),
    " -- missing ~", round(a_d_true - mean(N_hat_moderate), 0), "fish per day\n")
cat("Mean a_d_hat (r_true = 0.05, TikTok):  ",
    round(mean(N_hat_extreme), 1),
    " -- missing ~", round(a_d_true - mean(N_hat_extreme), 0), "fish per day\n\n")

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
abline(v = a_d_true,               col = "firebrick", lwd = 2)
abline(v = mean(N_hat_moderate),   col = "navy",      lwd = 2, lty = 2)
legend("topleft",
       legend = c("True N", "Mean a_d_hat"),
       col = c("firebrick", "navy"), lwd = 2, lty = c(1, 2), bty = "n")

# true N is ~16x the mean estimate here — can't show on the same axis
hist(N_hat_extreme,
     breaks = 20, col = "lightcoral", border = "white",
     xlab = "a_d_hat", main = "True r = 0.05, assumed 5/6\nCatastrophic bias")
abline(v = mean(N_hat_extreme), col = "navy", lwd = 2, lty = 2)
mtext(paste("True N =", a_d_true, "is off this axis to the right"),
      side = 3, line = 0.2, cex = 0.75, col = "firebrick")
legend("topright", legend = "Mean a_d_hat",
       col = "navy", lwd = 2, lty = 2, bty = "n")

par(op)

# Failure 2 — random daily r, CI too narrow ----
#
# Now r is correct on average but varies randomly week to week.
# The point estimate is approximately unbiased. The CI is not.
#
# EASE bootstrap variance assumes fixed r, capturing only:
#   N * r * (1 - r) / r^2
#
# True variance when r is random adds:
#   N^2 * Var(r) / r^2
#
# EASE misses that second term. The reported CI is too narrow.
#
# Why it matters for management:
#   Managers use the CI to judge whether escapement is safely above goal.
#   If the true CI overlaps the escapement threshold but EASE's doesn't,
#   a harvest gets greenlit that should have been held back.
#   Example:
#     Escapement goal = 3500 fish
#     EASE reports:    a_d_hat = 3800, 95% CI [3650, 3950] -> looks safe
#     True CI would be:              95% CI [3200, 4400] -> lower bound
#                                                           is below goal
#     The fish don't care that the point estimate was unbiased.
#
# Demonstrated by computing empirical CI coverage across many seasons.
# If the model is right, 95% of CIs should contain the true season total.

r_true_mean <- 5 / 6
r_sd        <- 0.08

beta_var  <- r_sd^2
alpha_par <- r_true_mean * (r_true_mean * (1 - r_true_mean) / beta_var - 1)
beta_par  <- (1 - r_true_mean) * (r_true_mean * (1 - r_true_mean) / beta_var - 1)

n_seasons     <- 2000
boots_cov     <- 500
ci_contains_N <- logical(n_seasons)

for (s in seq_len(n_seasons)) {
  r_week  <- rbeta(n_weeks, alpha_par, beta_par)        # true r each week
  W_week  <- rbinom(n_weeks, size = a_d_true, prob = r_week)

  # EASE bootstrap: expand W to point estimates first, then resample.
  # E[Binomial(N_hat_week, r) / r] = N_hat_week, so the bootstrap centers
  # at the EASE point estimate. Using raw W_week instead centers at
  # sum(W_week) ~ r * true total, which is wrong.
  N_hat_week  <- round(W_week / r_assumed)
  boot_season <- vapply(seq_len(boots_cov), function(b) {
    sum(rbinom(n_weeks, size = N_hat_week, prob = r_assumed) / r_assumed)
  }, numeric(1))

  ci_s          <- quantile(boot_season, c(0.025, 0.975))
  true_total    <- a_d_true * n_weeks
  ci_contains_N[s] <- ci_s[1] <= true_total & true_total <= ci_s[2]
}

empirical_coverage <- mean(ci_contains_N)

cat("Failure 2: random weekly r\n")
cat("Nominal CI coverage:   95%\n")
cat("Empirical CI coverage:", round(empirical_coverage * 100, 1), "%\n")
cat("If coverage << 95%, EASE is overconfident about its precision.\n\n")

# Failure 3 — r degrades systematically across the season ----
#
# Early season: good conditions, r close to 5/6 (the ceiling —
# you cannot observe more fish than pass during your watch window).
# Late season: fatigue and weather push r down toward 0.5.
#
# EASE assumes r = 5/6 for every week. Early weeks are fine.
# Late weeks are not — dividing by too large an r makes a_d_hat too low.
# No single week looks wrong enough to trigger an alarm, but the bias
# accumulates silently into the season total.

r_early       <- 5 / 6   # true r at start of season (matches assumption)
r_late        <- 0.5     # true r by end of season
r_true_season <- seq(r_early, r_late, length.out = n_weeks)  # linear decline

N_week_true        <- rep(a_d_true, n_weeks)
W_season           <- rbinom(n_weeks, size = N_week_true, prob = r_true_season)
N_hat_wrong        <- W_season / r_assumed       # EASE estimate (wrong r)
N_hat_correct_seas <- W_season / r_true_season   # what correct r would give

cat("Failure 3: seasonal r degradation (5/6 -> 0.5)\n")
cat("True season total:             ", sum(N_week_true), "\n")
cat("EASE season total (wrong r):   ", round(sum(N_hat_wrong), 0), "\n")
cat("EASE season total (correct r): ", round(sum(N_hat_correct_seas), 0), "\n")
cat("Bias:                          ",
    round(sum(N_hat_wrong) - sum(N_week_true), 0), "fish\n\n")

op <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

plot(
  seq_len(n_weeks), r_true_season,
  type = "b", pch = 16, col = "firebrick", lwd = 2,
  ylim = c(0.4, 1.0),
  xlab = "Statistical week", ylab = "Sampling fraction r",
  main = "True r degrades across season\nEASE assumes constant 5/6"
)
abline(h = r_assumed, col = "navy", lwd = 2, lty = 2)
legend("topright",
       legend = c("True r", "Assumed r"),
       col = c("firebrick", "navy"), lwd = 2, lty = c(1, 2), bty = "n")

cumulative_bias <- cumsum(N_hat_wrong) - cumsum(N_week_true)

plot(
  seq_len(n_weeks), cumulative_bias,
  type = "b", pch = 16, col = "firebrick", lwd = 2,
  xlab = "Statistical week",
  ylab = "Cumulative bias (N_hat - N_true)",
  main = "Undercount accumulates\nas season progresses"
)
abline(h = 0, col = "gray50", lty = 2)

par(op)

# Section summary ----
#
# When the model is right:
#   W ~ Binomial(a_d_true, r_sample), E(W) = a_d_true * r_sample,
#   a_d_hat = W / r_sample is unbiased.
#
# Failure 1 — wrong r, fixed:
#   E(a_d_hat) = a_d_true * r_true / r_assumed. If r_true < r_assumed,
#   every estimate is too low. Systematic, not noise. Compounds all season.
#
# Failure 2 — r varies randomly across weeks:
#   Point estimate approximately unbiased. CI is wrong. EASE misses the
#   N^2 * Var(r) variance term. Coverage drops well below 95%. Managers
#   make harvest decisions with false precision.
#
# Failure 3 — r degrades systematically across the season:
#   Early weeks fine. Late weeks silently biased low. Season total
#   accumulates a meaningful undercount invisible at the weekly level.
#
# Next: Section 2 connects the window-count model to the MLE for a_d_true,
# which is where escapeLGD's optimization step lives.
