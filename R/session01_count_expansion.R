# Session 1: one count, expanded, by likelihood
# Objective: explain out loud why a count at a dam divided by the fraction of
# fish it saw is a maximum likelihood estimate, and why that one idea is the
# base of both the adult (EASE) and smolt (SCRAPI2) estimators.

# Truth, set at the top. Generic names so both life stages map onto them.
#   total_count      truly passing fish (adult window total, or smolt run size)
#   sampled_fraction fraction the count structurally sees (window-open fraction
#                    for adults; sample rate times guidance for smolts)
truth <- list(total_count = 8000L, sampled_fraction = 0.10)
set.seed(1)

# Simulate. We observe only the counted fish; the rest pass unseen.
observed_count <- rbinom(1, size = truth$total_count, prob = truth$sampled_fraction)

# Build the estimator by hand. Given observed count y and known fraction r, the
# unknown total N enters a Binomial(N, r) likelihood. dbinom needs an integer
# size, so relax the binomial coefficient with lgamma and optimise over real N.
neg_loglik <- function(N, y, r) {
  if (N < y) return(1e10)
  logchoose <- lgamma(N + 1) - lgamma(y + 1) - lgamma(N - y + 1)
  -(logchoose + y * log(r) + (N - y) * log(1 - r))
}
fit <- optim(par = observed_count / truth$sampled_fraction, fn = neg_loglik,
             y = observed_count, r = truth$sampled_fraction, method = "Brent",
             lower = observed_count, upper = 3 * observed_count / truth$sampled_fraction)
N_hat <- fit$par

# Recover. The MLE is simply the expansion count/fraction. Interval by the exact
# parametric bootstrap production uses: resample the count, re-expand.
boot_N <- rbinom(2000, round(N_hat), truth$sampled_fraction) / truth$sampled_fraction
ci <- quantile(boot_N, c(0.05, 0.95))
cat(sprintf("truth %d | MLE %.0f | closed form y/r %.0f | 90%% CI [%.0f, %.0f]\n",
            truth$total_count, N_hat, observed_count / truth$sampled_fraction, ci[1], ci[2]))

png("figs/session01_count_expansion.png", width = 900, height = 600)
hist(boot_N, breaks = 40, col = "grey80", border = "white",
     main = "Session 1: expanded count and its bootstrap interval",
     xlab = "estimated total count")
abline(v = truth$total_count, col = "firebrick", lwd = 3)
abline(v = ci, col = "steelblue", lwd = 2, lty = 2)
dev.off()

# Exercise. Drop sampled_fraction to 0.02 and rerun. The point estimate stays
# roughly unbiased but the interval widens sharply; say out loud why the same
# count buys far less certainty when a smaller fraction is seen.

# Locate.
# escapeLGD (adults): expand_wc_binom_night() in
#   R/night_fall_reascend_wc_binom.R expands the window count as
#   round(wc / wc_prop) and bootstraps it as rbinom(boots, wc, wc_prop)/wc_prop,
#   the same expansion and parametric bootstrap used here, with wc_prop as the
#   sampled fraction. The adult count's only filter is the counting window, so
#   there is no guidance term.
# smoltEASE (smolts): thetahat() inside SCRAPI2() in R/SCRAPI2.R forms
#   dailypass <- Tally / Ptrue with Ptrue = SampleRate * GuidanceEfficiency, and
#   bootstraps it as rbinom(1, est_daily, Ptrue). Structurally identical; the
#   smolt sampled fraction is a product of two terms rather than one, which is
#   why guidance efficiency earns its own session. Both stages share this exact
#   expand-and-binomial-bootstrap.

writeLines(c(
  "How I would explain session 1 in three minutes",
  "",
  "We never count every fish at the dam. We see a known slice of them: the",
  "hours the counting window is open, or the fraction the trap sampled. If we",
  "saw a tenth of the fish and counted 800, our best guess for the whole run is",
  "8000. That division is not a rule of thumb; it is the single most likely",
  "total given what we saw, which is what a maximum likelihood estimate means.",
  "",
  "The catch is certainty. The point estimate barely moves when the seen slice",
  "shrinks, but the range around it grows fast, because a thin slice could have",
  "come from many different true totals. Everything else in both the adult and",
  "smolt estimators is layered on top of this one move: count, divide by the",
  "fraction seen, and carry the uncertainty forward."
), "docs/session01_explain.md")
