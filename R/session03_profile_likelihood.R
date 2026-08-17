# Session 3: profile likelihood, and four intervals off one curve
# Objective: explain how an interval falls out of the shape of the likelihood
# curve for one rate, and what the four interval recipes each mean on that rate.

# Truth: a single rate. Stands in for nighttime passage p_night on the adult
# side, or any one binomial rate resampled on the smolt side.
truth <- list(rate = 0.12, n_trials = 400L)
set.seed(3)

# Simulate. Observe k successes out of n: PIT tags detected passing at night out
# of all tags whose passage time was read.
k <- rbinom(1, truth$n_trials, truth$rate)
n <- truth$n_trials
p_hat <- k / n

# Build the log-likelihood by hand and walk it across a grid of the rate.
loglik <- function(p) dbinom(k, n, p, log = TRUE)
grid <- seq(1e-4, 1 - 1e-4, length.out = 4000)
ll   <- vapply(grid, loglik, numeric(1))

# Four intervals for the one rate, on identical data.
cutoff  <- max(ll) - qchisq(0.95, df = 1) / 2      # profile: chi-sq(1) drop
ci_prof <- range(grid[ll >= cutoff])
se       <- sqrt(p_hat * (1 - p_hat) / n)           # delta: symmetric Wald band
ci_delta <- p_hat + c(-1, 1) * qnorm(0.975) * se
ci_boot  <- quantile(rbinom(4000, n, p_hat) / n, c(0.025, 0.975))  # bootstrap
ci_bayes <- qbeta(c(0.025, 0.975), k + 0.5, n - k + 0.5)  # Bayes, Jeffreys prior

# Recover. Print the MLE, the profile interval read off the curve, and all four.
mat <- rbind(profile = ci_prof, delta = ci_delta, bootstrap = ci_boot, bayes = ci_bayes)
cat(sprintf("truth %.3f | MLE %.3f | profile interval [%.3f, %.3f]\n",
            truth$rate, p_hat, ci_prof[1], ci_prof[2]))
print(round(cbind(lower = mat[, 1], est = p_hat, upper = mat[, 2]), 3))
gap <- max(abs(mat[, 1] - mean(mat[, 1])), abs(mat[, 2] - mean(mat[, 2])))
cat(sprintf("widest endpoint gap across the four: %.4f\n", gap))

png("figs/session03_profile_likelihood.png", width = 900, height = 600)
plot(grid, ll, type = "l", lwd = 2, xlim = c(max(0, p_hat - 0.12), p_hat + 0.12),
     ylim = c(cutoff - 3, max(ll) + 0.5), xlab = "rate", ylab = "log-likelihood",
     main = "Session 3: an interval read off the likelihood curve")
abline(h = cutoff, col = "steelblue", lwd = 2, lty = 2)
abline(v = ci_prof, col = "steelblue", lwd = 2, lty = 3)
abline(v = truth$rate, col = "firebrick", lwd = 3)
points(p_hat, max(ll), pch = 19)
dev.off()

# Exercise. Set rate = 0.02 and n_trials = 60 and rerun. The curve turns broad and
# lopsided, the delta interval runs below zero while the other three stay positive
# and skew right; say why only the symmetric method breaks near a boundary.

# Locate.
# escapeLGD (adults): nightFall() estimates p_night = nightPass / totalPass, the
#   k/n MLE here; the fallback branch p_fa = numReascend / totalPass is the same
#   one-parameter binomial. Both are carried by binomial bootstrap.
# smoltEASE (smolts): SCRAPI2() has no standalone rate to profile; each daily rate
#   is resampled inside the bootstrap. Of the four, the bootstrap percentile is the
#   one both tools ship: quantile() on the bootstrap matrix in apply_fallback_rates()
#   and quantile(theta.b) in SCRAPI2(). Delta and profile are what the by-hand
#   likelihood implies; the Bayesian credible interval is what fit_ge_model()
#   produces for guidance efficiency, the interval session 5 fits.

dir <- if (gap < 0.01) "land in nearly the same place" else "already pull apart"
writeLines(c(
"How I would explain session 3 in three minutes", "",
sprintf("When we estimate one rate, like the share of fish slipping past at night, we"),
sprintf("get one best number: detected over total is %.3f. A range matters more, and the", p_hat),
"honest way to get one is to ask which other rates explain the data almost as well.",
"Plot that and a hump peaks at the estimate; drop a line a fixed amount below the",
sprintf("peak and where it cuts the hump are the endpoints, here [%.3f, %.3f], with no", ci_prof[1], ci_prof[2]),
"symmetry assumed.", "",
"There is more than one recipe for that range: a symmetric delta band from a standard",
"error, the likelihood read off the curve, a bootstrap of the count, and a Bayesian",
sprintf("blend of prior and data. Here the widest the four disagree at either end is %.4f,", gap),
sprintf("so they %s. They split only when data are thin or the rate hugs a", dir),
"boundary, and there the delta band misleads first. Both tools ship the bootstrap."
), "docs/session03_explain.md")
