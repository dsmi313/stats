# Session 6: what the composed interval actually claims
# Objective: explain what the SCRAPI2 interval claims to cover and whether a
# simulation says it delivers, and how the adult interval compares.

# The SCRAPI2 interval stacks a count bootstrap on posterior draws of GE (and
# GSI). It claims about 90% coverage of true escapement, but only if the GE model
# it is fed is centered right. To test that fairly, the true GE varies run to run
# from the very distribution a calibrated interval integrates over.
truth <- list(N = 6000L, rate = 0.05, detect = 0.17)  # detect = GE or night rate
set.seed(8)

# One coverage engine. post_draw(B, d_true) supplies the B detection values the
# interval integrates over. The count bootstrap uses the point detection, the
# divisor uses the draws, exactly as SCRAPI2 uses est_daily against ge_day_mat[, b].
run_cov <- function(post_draw, nsim = 500, B = 400, alpha = 0.1, conc = 400) {
  N <- truth$N; R <- truth$rate; hit <- logical(nsim)
  for (s in seq_len(nsim)) {
    d_true <- rbeta(1, truth$detect * conc, (1 - truth$detect) * conc)
    y  <- rbinom(1, N, R * d_true)
    d  <- post_draw(B, d_true)
    dh <- mean(d)
    yb <- rbinom(B, round(y / (R * dh)), R * dh)
    ci <- quantile(yb / (R * d), c(alpha / 2, 1 - alpha / 2))
    hit[s] <- ci[1] <= N && N <= ci[2]
  }
  mean(hit)
}

# A GE posterior fixed at mean m (integrates the population GE spread), and an
# adult-style bootstrap of a night rate estimated from 250 PIT tags.
ge_post <- function(m) function(B, dt) rbeta(B, m * 400, (1 - m) * 400)
res <- c(
  "smolt GE (calibrated)"   = run_cov(ge_post(truth$detect)),
  "smolt GE (biased +10%)"  = run_cov(ge_post(truth$detect * 1.10)),
  "adult night (bootstrap)" = run_cov(function(B, dt) { ph <- rbinom(1, 250, dt) / 250; rbinom(B, 250, ph) / 250 }))
print(round(res, 3))

png("figs/session06_interval_coverage.png", width = 900, height = 600)
barplot(res, ylim = c(0, 1), col = "steelblue", las = 1, ylab = "coverage",
        main = "Session 6: does the 90% interval cover 90%?")
abline(h = 0.90, col = "firebrick", lwd = 3)
dev.off()

# Exercise. Push the bias to 1.20 and watch coverage fall further, then set it
# back to 1.00 and instead widen the true-GE spread (conc = 100) while leaving the
# posterior at conc = 400. Say out loud which failure is a centering problem and
# which is a spread problem, and which one honesty about the model can fix.

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

cal <- res[[1]]; bia <- res[[2]]; adu <- res[[3]]
drop <- 0.90 - bia
writeLines(c(
"How I would explain session 6 in three minutes",
"",
"Our escapement range is built by two moves stacked together: resample the catch,",
"and draw guidance efficiency from its model. We call it a 90 percent interval, so",
"the fair test is to invent a world where we know the true answer, run the whole",
"machine hundreds of times, and count how often the interval actually catches the",
"truth. It should catch it about 90 times in 100.",
"",
"When the guidance-efficiency model is centered correctly, it does: coverage here",
sprintf("is %.0f percent. Bias the model up ten percent and it still looks tight and", 100 * cal),
sprintf("confident but lands off-center, and coverage slips to %.0f percent, short of the", 100 * bia),
sprintf("nominal by about %.0f points. The adult interval leans on a resampled rate rather", 100 * drop),
sprintf("than a fitted model, so it cannot be miscentered that way and holds near %.0f. The", 100 * adu),
"lesson: a composed interval is only as honest as the model feeding it, and we",
"should say so out loud when we quote one."
), "docs/session06_explain.md")
