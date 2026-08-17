# Session 5: guidance efficiency fit both ways
# Objective: explain what guidance efficiency is, why it needs its own model, and
# what fitting it by maximum likelihood versus Bayesian buys and costs.

# Simulate-only, so truth stays known. logit(GE) is linear in standardized spill,
# shaped like the real MY2025 pattern (see data/): GE low, falling as spill rises.
# Few strata, so the two fits can actually disagree. GE is route-selection psi.
truth <- list(alpha = -1.4, beta = -0.7, n_strata = 8L, n_tag = 60L)
set.seed(6)
spill_std <- scale(sort(runif(truth$n_strata, 0, 100)))[, 1]
psi   <- plogis(truth$alpha + truth$beta * spill_std)
nsz   <- rep(truth$n_tag, truth$n_strata)
y     <- rbinom(truth$n_strata, nsz, psi)   # tagged fish taking the bypass route

# Maximum likelihood in glmmTMB (falls back to glm if glmmTMB is absent).
if (requireNamespace("glmmTMB", quietly = TRUE)) {
  m  <- glmmTMB::glmmTMB(cbind(y, nsz - y) ~ spill_std, family = binomial, data = data.frame(y, nsz, spill_std))
  co <- summary(m)$coefficients$cond
} else {
  message("glmmTMB not installed; using glm for the ML fit.")
  co <- summary(glm(cbind(y, nsz - y) ~ spill_std, binomial))$coefficients
}
ml <- data.frame(est = co[, 1], lo = co[, 1] - 1.96 * co[, 2], hi = co[, 1] + 1.96 * co[, 2])

# Bayesian by hand: a short Metropolis sampler over (alpha, beta), weak priors.
logpost <- function(p) sum(dbinom(y, nsz, plogis(p[1] + p[2] * spill_std), log = TRUE)) +
  dnorm(p[1], 0, 5, log = TRUE) + dnorm(p[2], 0, 5, log = TRUE)
niter <- 30000; chain <- matrix(NA_real_, niter, 2); cur <- c(0, 0); lc <- logpost(cur)
for (i in 1:niter) {
  prop <- cur + rnorm(2, 0, c(0.25, 0.3)); lp <- logpost(prop)
  if (log(runif(1)) < lp - lc) { cur <- prop; lc <- lp }
  chain[i, ] <- cur
}
post <- chain[-(1:5000), ]
bayes <- data.frame(est = colMeans(post), lo = apply(post, 2, quantile, 0.025),
                    hi = apply(post, 2, quantile, 0.975))

# Recover. Compare truth against both fits, and draw the GE-vs-spill curves.
cat(sprintf("truth alpha %.2f beta %.2f\n", truth$alpha, truth$beta))
cat(sprintf("ML    alpha %.2f [%.2f, %.2f]  beta %.2f [%.2f, %.2f]\n",
            ml$est[1], ml$lo[1], ml$hi[1], ml$est[2], ml$lo[2], ml$hi[2]))
cat(sprintf("Bayes alpha %.2f [%.2f, %.2f]  beta %.2f [%.2f, %.2f]\n",
            bayes$est[1], bayes$lo[1], bayes$hi[1], bayes$est[2], bayes$lo[2], bayes$hi[2]))

png("figs/session05_ge_fit_both_ways.png", width = 900, height = 600)
xs <- seq(min(spill_std), max(spill_std), length.out = 100)
plot(spill_std, y / nsz, pch = 19, ylim = c(0, 1), xlab = "standardized spill",
     ylab = "guidance efficiency", main = "Session 5: GE fit by ML and Bayesian")
lines(xs, plogis(truth$alpha + truth$beta * xs), col = "firebrick", lwd = 3)
lines(xs, plogis(ml$est[1] + ml$est[2] * xs), col = "black", lwd = 2, lty = 2)
lines(xs, plogis(bayes$est[1] + bayes$est[2] * xs), col = "darkorange", lwd = 2, lty = 3)
legend("topright", c("truth", "ML", "Bayes"), lwd = c(3, 2, 2), lty = c(1, 2, 3), col = c("firebrick", "black", "darkorange"), bty = "n")
dev.off()

# Exercise. Set n_strata = 5L and n_tag = 2L and rerun. With two tags a stratum the
# counts (2,2,1,0,0) separate perfectly along spill: the ML slope runs off to
# beta = -33.7 with SE = 161087 (glmmTMB warns), a useless interval, while the prior
# keeps the Bayesian slope finite near -5.9, 95% [-13.3, -1.7]. Say what that buys.

# Locate.
# smoltEASE (smolts): fit_ge_model() in R/fit_ge_model.R fits this for real as a
#   Bayesian multistate mark-recapture in JAGS, with logit(psi) linear in spill
#   plus stratum process error and nested shrinkage, far richer than the two-
#   parameter logistic here; generate_ge_draws()/prep_ge_data() turn PIT
#   detections into the posterior draws SCRAPI2 consumes.
# escapeLGD (adults): has no GE model, because the adult count has no bypass to
#   route into. The analogue, nighttime passage, is a plain binomial rate in
#   nightFall(), not a fitted curve.

wml <- ml$hi[2] - ml$lo[2]; wbayes <- bayes$hi[2] - bayes$lo[2]
agree <- if (abs(ml$est[2] - bayes$est[2]) < 0.1) "land on essentially the same slope" else "disagree on the slope"
writeLines(c(
"How I would explain session 5 in three minutes",
"",
"Guidance efficiency is the share of smolts taking the bypass route where the trap",
"can sample them. It changes with spill, and under heavy spill it is both low and",
"hard to pin down, so we cannot plug in one number; it needs its own model.",
"",
"We can fit that model two ways. Maximum likelihood finds the single best curve and",
sprintf("a range from its curvature; here it recovers the spill slope %.2f against a truth", ml$est[2]),
sprintf("of %.2f. The Bayesian fit adds a mild prior; here it gives %.2f. On this run the", truth$beta, bayes$est[2]),
sprintf("two %s, with slope intervals %.2f and %.2f wide. The prior earns its keep", agree, wml, wbayes),
"only when data are thin enough to separate; the production smolt model is Bayesian",
"for exactly that thin-data reason."
), "docs/session05_explain.md")
