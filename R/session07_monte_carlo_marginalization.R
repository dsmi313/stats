# Session 7: Monte Carlo marginalization
# Objective: explain how drawing posterior values and pushing each through the
# escapement calculation is the same thing as integrating the nuisance out.

# A daily count expanded to escapement depends on guidance efficiency GE, which
# is not known exactly: the GE model leaves a posterior. Truth here is that
# posterior (a Beta standing in for the fit_ge_model draws) plus a fixed count
# and sample rate.
truth <- list(count = 1200L, rate = 0.05, ge_a = 18, ge_b = 90)  # mean GE ~ 0.17
set.seed(7)

# Escapement as a function of the nuisance GE. GE enters through 1/GE.
escapement <- function(ge) truth$count / (truth$rate * ge)

# Deterministic marginalization: integrate escapement over the GE posterior on a
# fine grid. This is the session-4 move written as an explicit integral.
g  <- seq(0.02, 0.5, length.out = 5000)
w  <- dbeta(g, truth$ge_a, truth$ge_b); w <- w / sum(w)
int_mean <- sum(escapement(g) * w)
Ngrid <- escapement(g); ord <- order(Ngrid); cw <- cumsum(w[ord])
int_ci <- Ngrid[ord][c(which.min(abs(cw - 0.05)), which.min(abs(cw - 0.95)))]

# Monte Carlo marginalization: draw GE from the same posterior, one draw per
# iteration, exactly what generate_ge_draws() hands SCRAPI2, and push each
# through the calculation. No integral is written; the draws do the integrating.
B <- 20000
N_draws <- escapement(rbeta(B, truth$ge_a, truth$ge_b))
mc_mean <- mean(N_draws); mc_ci <- quantile(N_draws, c(0.05, 0.95))

cat(sprintf("integral mean %.0f   MC mean %.0f\n", int_mean, mc_mean))
cat(sprintf("integral 90%% [%.0f, %.0f]   MC 90%% [%.0f, %.0f]\n",
            int_ci[1], int_ci[2], mc_ci[1], mc_ci[2]))

png("figs/session07_monte_carlo_marginalization.png", width = 900, height = 600)
run <- cumsum(N_draws) / seq_len(B)
plot(run, type = "l", lwd = 2, xlab = "Monte Carlo draws", ylab = "mean escapement",
     ylim = int_mean + c(-1, 1) * 400,
     main = "Session 7: drawing GE equals integrating GE out")
abline(h = int_mean, col = "firebrick", lwd = 3)
legend("topright", c("running MC mean", "grid integral"), lwd = c(2, 3),
       col = c("black", "firebrick"), bty = "n")
dev.off()

# Exercise. Replace GE with a second nuisance, a per-fish stock draw, and average
# escapement over both draws at once. Say out loud why stacking two draw sources
# still just integrates over their joint distribution, one iteration at a time.

# Locate.
# smoltEASE (smolts): the SCRAPI2() bootstrap loop is this demonstration in
#   production. Iteration b reads ge_day_mat[, b] (one GE posterior column, fed by
#   generate_ge_draws()) and one gsiDraws column, pushes both through thetahat(),
#   and the spread of the results is the marginal over GE and GSI, computed by
#   drawing rather than integrating.
# escapeLGD (adults): HNC_expand_unkGSI() consumes one GSI posterior column per
#   iteration the same way, marginalizing genetic stock; nightFall() resamples
#   passage and fallback rates with rbinom, which marginalizes those rates. Same
#   Monte Carlo integral, same identity, both life stages.

writeLines(c(
"How I would explain session 7 in three minutes",
"",
"When a number we need, like guidance efficiency, is uncertain, we do not want a",
"single plugged-in value; we want to account for every value it might take,",
"weighted by how likely each is. Written as math that is an integral, and it can",
"look intimidating.",
"",
"There is a shortcut that gives the identical answer. Draw a few thousand values",
"of the uncertain number from its range, run the escapement calculation once for",
"each, and look at the spread of results. The average lands on the integral and",
"the middle band is the honest interval. That is all the bootstrap loop is doing",
"when it reads one guidance-efficiency draw and one stock draw per pass: it is",
"integrating out our uncertainty by drawing, not by calculus."
), "docs/session07_explain.md")
