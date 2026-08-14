# Session 8: what the composed interval actually claims
# Objective: explain what the SCRAPI2 interval claims to cover and whether a
# simulation says it delivers, and how the adult interval compares.

# The SCRAPI2 interval stacks a count bootstrap on posterior draws of GE (and
# GSI). It claims about 90% coverage of true escapement, but only if the
# posteriors it is fed are correct. This study checks that claim three ways.
truth <- list(N = 6000L, rate = 0.05, detect = 0.17)  # detect = GE or night rate
set.seed(8)

# One coverage engine. draw_detect(B) supplies the B detection values the interval
# integrates over: a GE posterior (smolt) or a bootstrapped rate (adult). The
# count bootstrap uses the point detection, the divisor uses the draws, exactly
# as SCRAPI2 uses est_daily against ge_day_mat[, b].
run_cov <- function(draw_detect, nsim = 300, B = 400, alpha = 0.1) {
  N <- truth$N; R <- truth$rate; d0 <- truth$detect
  hit <- logical(nsim)
  for (s in seq_len(nsim)) {
    y    <- rbinom(1, N, R * d0)
    d    <- draw_detect(B)
    dhat <- mean(d)
    yb   <- rbinom(B, round(y / (R * dhat)), R * dhat)
    ci   <- quantile(yb / (R * d), c(alpha / 2, 1 - alpha / 2))
    hit[s] <- ci[1] <= N && N <= ci[2]
  }
  mean(hit)
}

# A Beta posterior with a target mean and modest spread (concentration 200).
beta_post <- function(mean_d, conc = 200) function(B) rbeta(B, mean_d * conc, (1 - mean_d) * conc)

res <- c(
  "smolt GE (calibrated)"   = run_cov(beta_post(truth$detect)),
  "smolt GE (biased +10%)"  = run_cov(beta_post(truth$detect * 1.10)),
  "adult night (bootstrap)" = run_cov(function(B) rbinom(B, 250, truth$detect) / 250))
print(round(res, 3))

png("figs/session08_interval_coverage.png", width = 900, height = 600)
barplot(res, ylim = c(0, 1), col = "steelblue", las = 1, ylab = "coverage",
        main = "Session 8: does the 90% interval cover 90%?")
abline(h = 0.90, col = "firebrick", lwd = 3)
dev.off()

# Exercise. Push the bias to +25% and watch coverage collapse, then set it back to
# 0 and widen the posterior spread instead (concentration 40). Say out loud which
# failure a wider-but-centered posterior fixes and which it does not.

# Locate.
# smoltEASE (smolts): SCRAPI2() builds its CI as quantile(theta.b, c(alph/2,
#   1-alph/2)) over a bootstrap that folds in one GE draw and one GSI draw per
#   iteration. The study shows what that construction quietly assumes: the GE
#   posterior must be centered right, or the interval covers the wrong place.
# escapeLGD (adults): apply_fallback_rates() builds its CI as quantile() over the
#   composition bootstrap multiplied by the fallback bootstrap. Its detection-side
#   input, nighttime passage, is a bootstrapped binomial rate rather than a fitted
#   posterior, so it cannot be miscentered by a model the way GE can; that is the
#   practical difference the coverage bars show.

writeLines(c(
"How I would explain session 8 in three minutes",
"",
"Our escapement range is built by two moves stacked together: resample the catch,",
"and draw guidance efficiency from its model. We call it a 90 percent interval, so",
"the fair test is to invent a world where we know the true answer, run the whole",
"machine hundreds of times, and count how often the interval actually catches the",
"truth. It should be 90 out of 100.",
"",
"When the guidance-efficiency model is centered correctly, it is. When that model",
"is biased even ten percent, the interval still looks tight and confident but",
"lands in the wrong place, and coverage falls well short. The adult interval leans",
"on a resampled rate instead of a fitted model, so it cannot be thrown off that",
"particular way. The lesson: a composed interval is only as honest as the model",
"feeding it, and we should say so out loud when we quote one."
), "docs/session08_explain.md")
