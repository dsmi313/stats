# Session 3: profile likelihood
# Objective: explain how an interval falls out of the shape of the likelihood
# curve for one parameter, without any normal approximation.

# Truth: a single rate. This stands in for nighttime passage p_night on the
# adult side, or any one binomial rate resampled on the smolt side.
truth <- list(rate = 0.12, n_trials = 400L)
set.seed(3)

# Simulate. Observe k successes out of n. Here: PIT tags detected passing at
# night out of all tags whose passage time was read.
k <- rbinom(1, truth$n_trials, truth$rate)
n <- truth$n_trials

# Build the log-likelihood by hand and walk it across a grid of the rate.
loglik <- function(p) dbinom(k, n, p, log = TRUE)
grid <- seq(0.001, 0.999, length.out = 2000)
ll <- vapply(grid, loglik, numeric(1))
p_hat <- k / n

# Recover. The likelihood-ratio interval is every rate whose log-likelihood sits
# within half the chi-squared(1) cutoff of the peak. Read the endpoints straight
# off the curve; no standard error, no symmetry assumed.
cutoff <- max(ll) - qchisq(0.95, df = 1) / 2
inside <- grid[ll >= cutoff]
ci <- range(inside)
cat(sprintf("truth %.3f | MLE %.3f | 95%% profile interval [%.3f, %.3f]\n",
            truth$rate, p_hat, ci[1], ci[2]))

png("figs/session03_profile_likelihood.png", width = 900, height = 600)
plot(grid, ll, type = "l", lwd = 2, xlim = c(max(0, p_hat - 0.12), p_hat + 0.12),
     ylim = c(cutoff - 3, max(ll) + 0.5), xlab = "rate", ylab = "log-likelihood",
     main = "Session 3: an interval read off the likelihood curve")
abline(h = cutoff, col = "steelblue", lwd = 2, lty = 2)
abline(v = ci, col = "steelblue", lwd = 2, lty = 3)
abline(v = truth$rate, col = "firebrick", lwd = 3)
points(p_hat, max(ll), pch = 19)
dev.off()

# Exercise. Cut n_trials to 40 and rerun. The curve gets broad and visibly
# lopsided, and the interval stops looking symmetric around the MLE; say out
# loud why a normal-approximation interval would misplace the endpoints here.

# Locate.
# escapeLGD (adults): nightFall() in R/night_fall_reascend_wc_binom.R estimates
#   the nighttime passage rate as p_night = nightPass / totalPass, exactly the
#   k/n MLE here, and carries its uncertainty by binomial bootstrap rather than
#   by reading this curve. The single-rate fallback branch (p_fa =
#   numReascend / totalPass) is the same one-parameter binomial.
# smoltEASE (smolts): SCRAPI2() has no standalone rate to profile; each daily
#   rate is resampled inside the bootstrap (rbinom on est_daily). The one-
#   parameter likelihood here is the atom those resamples approximate, so the
#   two production tools land on the same interval by simulation that this
#   session reaches by curve shape.

writeLines(c(
"How I would explain session 3 in three minutes",
"",
"When we estimate a single rate, like the share of fish that slip past at night,",
"we get one best number: detected over total. But a range matters more than the",
"number. The honest way to get the range is to ask which other rates explain the",
"data almost as well as the best one.",
"",
"Plot how well each candidate rate fits, and you get a hump peaking at the best",
"estimate. Draw a horizontal line a fixed drop below the peak and the rates",
"above the line are the ones the data cannot rule out; where the line cuts the",
"hump are the interval endpoints. Nothing here assumes the range is symmetric,",
"which is why it stays honest when the data are thin and the hump is lopsided."
), "docs/session03_explain.md")
