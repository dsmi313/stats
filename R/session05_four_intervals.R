# Session 5: four intervals for one parameter
# Objective: explain what each of four intervals means and why they agree here
# and can disagree elsewhere.

# One rate again, so the four methods are compared on identical data. Stands in
# for any single rate in either estimator (a passage rate, a fallback rate).
truth <- list(rate = 0.15, n_trials = 300L)
set.seed(5)
k <- rbinom(1, truth$n_trials, truth$rate); n <- truth$n_trials
p_hat <- k / n

# 1. Delta method (Wald): estimate plus/minus 1.96 standard errors, symmetric.
se <- sqrt(p_hat * (1 - p_hat) / n)
ci_delta <- p_hat + c(-1, 1) * qnorm(0.975) * se

# 2. Profile likelihood: rates within half the chi-squared(1) drop of the peak.
grid <- seq(1e-4, 1 - 1e-4, length.out = 4000)
ll <- dbinom(k, n, grid, log = TRUE)
ci_prof <- range(grid[ll >= max(ll) - qchisq(0.95, 1) / 2])

# 3. Bootstrap percentile: resample the count, re-estimate, take quantiles.
boot <- rbinom(4000, n, p_hat) / n
ci_boot <- quantile(boot, c(0.025, 0.975))

# 4. Bayesian credible: Jeffreys prior Beta(0.5,0.5), conjugate Beta posterior.
ci_bayes <- qbeta(c(0.025, 0.975), k + 0.5, n - k + 0.5)

# Recover. Print all four and draw them as stacked bars against the truth.
mat <- rbind(delta = ci_delta, profile = ci_prof, bootstrap = ci_boot, bayes = ci_bayes)
print(round(cbind(lower = mat[, 1], est = p_hat, upper = mat[, 2]), 3))

png("figs/session05_four_intervals.png", width = 900, height = 600)
plot(NA, xlim = range(mat) + c(-0.02, 0.02), ylim = c(0.5, 4.5), yaxt = "n",
     xlab = "rate", ylab = "", main = "Session 5: four intervals for one rate")
axis(2, at = 4:1, labels = rownames(mat), las = 1)
for (i in 1:4) {
  y <- 5 - i
  segments(mat[i, 1], y, mat[i, 2], y, lwd = 4, col = "steelblue")
  points(p_hat, y, pch = 19)
}
abline(v = truth$rate, col = "firebrick", lwd = 3)
dev.off()

# Exercise. Set rate = 0.02 and n_trials = 60 and rerun. The delta interval runs
# below zero while the other three stay positive and skew right; say out loud why
# only the symmetric method breaks near the boundary.

# Locate.
# escapeLGD and smoltEASE both ship only the bootstrap percentile: quantile() on
#   the bootstrap matrix in apply_fallback_rates() (adults) and quantile(theta.b)
#   in SCRAPI2() (smolts). The profile and delta intervals here are what the
#   by-hand likelihood work implies but neither tool computes.
# The Bayesian credible interval is real production on the smolt side: fit_ge_model()
#   returns a posterior for guidance efficiency, and its central quantiles are a
#   credible interval of exactly this kind. That is the interval session 6 fits.

writeLines(c(
"How I would explain session 5 in three minutes",
"",
"There is more than one way to put a range on an estimate, and four of them show",
"up in this work. The delta method draws a symmetric band using a standard error.",
"The likelihood method reads the range off the fit curve. The bootstrap resamples",
"the data thousands of times and takes the middle 95 percent. The Bayesian method",
"combines the data with a prior and reports the middle of the resulting belief.",
"",
"When the sample is large and the rate is well away from zero or one, all four",
"land in nearly the same place, so it hardly matters which we quote. They only",
"split apart when data are thin or the rate hugs a boundary, and there the",
"symmetric delta band is the first to mislead. The tools ship the bootstrap; the",
"Bayesian one is what the guidance-efficiency model produces."
), "docs/session05_explain.md")
