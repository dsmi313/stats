# Session 9: talk assembly
# Objective: deliver the talk from three figures: the shared skeleton, the
# profiling core, and the coverage result. No new estimator work; this reruns the
# carrying pieces at small size and lays them out as one talk figure.
set.seed(9)

png("figs/session09_talk_assembly.png", width = 1500, height = 520)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), cex = 1.05)

# Panel 1: the shared skeleton. Six shared rows; the divergences in orange.
comp <- c("count", "sampled fraction", "detection expansion", "composition",
          "stage adjustment", "interval")
flag <- c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE)
plot(NA, xlim = c(0, 3), ylim = c(0.5, 6.8), axes = FALSE, xlab = "", ylab = "",
     main = "One skeleton, two life stages")
text(0.05, 6:1, comp, adj = 0, col = ifelse(flag, "darkorange", "black"))
text(2.25, 6.7, "adult smolt", adj = 0, cex = 0.85)
for (i in 1:6) points(c(2.3, 2.7), c(7 - i, 7 - i), pch = 19,
                      col = ifelse(flag[i], "darkorange", "steelblue"))
legend("bottomleft", c("shared", "diverges"), pch = 19,
       col = c("steelblue", "darkorange"), bty = "n", cex = 0.9)

# Panel 2: the profiling core. A likelihood curve and the interval read off it.
k <- 45; n <- 400; grid <- seq(1e-3, 0.3, length.out = 2000)
ll <- dbinom(k, n, grid, log = TRUE)
cut <- max(ll) - qchisq(0.95, 1) / 2; ci <- range(grid[ll >= cut])
plot(grid, ll, type = "l", lwd = 2, xlim = c(0.05, 0.18), ylim = c(cut - 3, max(ll) + 0.5),
     xlab = "rate", ylab = "log-likelihood", main = "An interval from the curve")
abline(h = cut, col = "steelblue", lwd = 2, lty = 2)
abline(v = ci, col = "steelblue", lwd = 2, lty = 3)

# Panel 3: the coverage result. Does the 90% interval cover 90%? The true GE
# varies run to run, so the calibrated interval sits near nominal and the biased
# one under-covers (the session 6 engine, condensed).
cov <- function(mean_d) {
  N <- 6000; R <- 0.05; d0 <- 0.17; h <- logical(500)
  for (s in 1:500) {
    dt <- rbeta(1, d0 * 400, (1 - d0) * 400); y <- rbinom(1, N, R * dt)
    d <- rbeta(400, mean_d * 400, (1 - mean_d) * 400)
    dh <- mean(d); yb <- rbinom(400, round(y / (R * dh)), R * dh)
    ci <- quantile(yb / (R * d), c(.05, .95)); h[s] <- ci[1] <= N && N <= ci[2]
  }
  mean(h)
}
covbars <- c(calibrated = cov(0.17), `biased +20%` = cov(0.204))
barplot(covbars, ylim = c(0, 1),
        col = "steelblue", ylab = "coverage", main = "Does 90% cover 90%?")
abline(h = 0.90, col = "firebrick", lwd = 3)
dev.off()
cat(sprintf("talk figure written; coverage calibrated %.2f vs biased %.2f\n",
            covbars[1], covbars[2]))

# Exercise. Rehearse the talk from these three panels alone, in order, spending
# one minute on each. If a panel needs a sentence the figure does not support,
# that sentence belongs in an earlier session, not on the slide.

# Locate. This session is assembly, so its anchors are the sessions it draws from:
# the skeleton (session 2, read from both sources), the profile (sessions 3 and 4,
# the escapeLGD fallback likelihood), and coverage (session 6, the SCRAPI2 and
# apply_fallback_rates CI construction). The real-data check (session 7) and the
# composition divergence (session 8) stand behind it. No new production code here.

writeLines(c(
"How I would give the talk in three minutes",
"",
"Escapement at Lower Granite is estimated for adults and for smolts, and it looks",
"like two separate methods. It is not. Both count the fish they can see, divide by",
"the fraction they saw, expand for the fish the count structurally misses, split",
"the total into groups, and put an honest range on the answer. One machine, two",
"life stages. Three parts genuinely differ, and naming them is most of the talk.",
"",
"The range is not decoration. It comes from the shape of the likelihood and from",
"drawing our uncertainties through the calculation, which is the same as averaging",
"over what we do not know. And when we quote a ninety percent interval, we owe it",
sprintf("a check: does it cover ninety percent? Panel three says it covers %.0f%% when the",
        100 * covbars[1]),
sprintf("guidance-efficiency model is centered right, and slips to %.0f%% when that model is",
        100 * covbars[2]),
"off by twenty percent, quietly and confidently. That is the whole talk: one",
"framework, honest intervals, and knowing exactly what each interval claims."
), "docs/session09_explain.md")
