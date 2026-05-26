# Studio 10 - Bootstrap CIs for mean, median, sd of fish lengths
# ---------------------------------------------------------------
# Source: https://ocw.mit.edu/courses/18-05-introduction-to-probability-and-statistics-spring-2022/mit18_05_s22_studio10-instructions.pdf
#
# MIT Studio 10 compares percentile vs basic bootstrap CIs for mean,
# median, and standard deviation. Fish version: smolt lengths simulated
# from either Normal or Log-Normal distributions. Fulton-style body
# masses and lengths are commonly log-normal in field data, so the
# log-normal case is the realistic stress test.

set.seed(2026)

# Helper: column-wise statistics via matrixStats if available -----------
col_stat <- function(mat, fun) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    switch(fun,
           mean   = matrixStats::colMeans2(mat),
           median = matrixStats::colMedians(mat),
           sd     = matrixStats::colSds(mat))
  } else {
    apply(mat, 2, get(fun))
  }
}

# Bootstrap CI sweep for a generic statistic ----------------------------
bootstrap_cis <- function(sample, stat_fun, n_boot, confidence) {
  n      <- length(sample)
  boots  <- matrix(sample(sample, size = n * n_boot, replace = TRUE),
                   nrow = n, ncol = n_boot)
  draws  <- col_stat(boots, stat_fun)
  obs    <- get(stat_fun)(sample)
  alpha  <- 1 - confidence
  qs_d   <- quantile(draws - obs, c(1 - alpha/2, alpha/2))
  percentile <- quantile(draws, c(alpha/2, 1 - alpha/2))
  basic      <- obs - qs_d
  list(point      = obs,
       percentile = unname(percentile),
       basic      = unname(basic))
}

ci_error <- function(true_value, ci) {
  !(true_value >= ci[1] && true_value <= ci[2])
}

# Problem 1: type-1 CI error rates on Normal data -----------------------
studio10_problem_1_fish <- function(true_mean, true_sd, n_data,
                                    n_boot, n_trials, confidence) {
  err <- list(mean = list(pct = 0L, basic = 0L),
              median = list(pct = 0L, basic = 0L),
              sd = list(pct = 0L, basic = 0L))
  true_median <- true_mean
  for (t in seq_len(n_trials)) {
    x <- rnorm(n_data, mean = true_mean, sd = true_sd)
    for (stat in c("mean", "median", "sd")) {
      true_v <- switch(stat,
                       mean = true_mean,
                       median = true_median,
                       sd = true_sd)
      cis <- bootstrap_cis(x, stat, n_boot, confidence)
      err[[stat]]$pct   <- err[[stat]]$pct   + ci_error(true_v, cis$percentile)
      err[[stat]]$basic <- err[[stat]]$basic + ci_error(true_v, cis$basic)
    }
  }
  cat("Studio 10 Problem 1: Normal data CI error rates (",
      n_trials, " trials, confidence ", confidence, ")\n", sep = "")
  for (stat in c("mean", "median", "sd"))
    cat(sprintf("  %-7s percentile %.4f   basic %.4f\n",
                stat,
                err[[stat]]$pct   / n_trials,
                err[[stat]]$basic / n_trials))
  cat("\n")
}
studio10_problem_1_fish(true_mean = 100, true_sd = 12,
                        n_data = 30, n_boot = 1000,
                        n_trials = 400, confidence = 0.95)

# Problem 2a: log-normal mean/median/sd from rlnorm parameters ---------
studio10_problem_2a_fish <- function(meanlog, sdlog) {
  mu     <- exp(meanlog + sdlog^2 / 2)
  med    <- exp(meanlog)
  sigma2 <- (exp(sdlog^2) - 1) * exp(2 * meanlog + sdlog^2)
  cat("Studio 10 Problem 2a: log-normal moments\n")
  cat("  meanlog =", meanlog, "   sdlog =", sdlog, "\n")
  cat("  mean   =", round(mu, 3),
      "   median =", round(med, 3),
      "   sd =", round(sqrt(sigma2), 3), "\n\n")
  invisible(list(mean = mu, median = med, sd = sqrt(sigma2)))
}
studio10_problem_2a_fish(meanlog = 4.6, sdlog = 0.3)

# Problem 2b: bootstrap CI error rates on log-normal data --------------
studio10_problem_2b_fish <- function(meanlog, sdlog, n_data,
                                     n_boot, n_trials, confidence) {
  true_vals <- studio10_problem_2a_fish(meanlog, sdlog)
  err <- list(mean = list(pct = 0L, basic = 0L),
              median = list(pct = 0L, basic = 0L),
              sd = list(pct = 0L, basic = 0L))
  for (t in seq_len(n_trials)) {
    x <- rlnorm(n_data, meanlog = meanlog, sdlog = sdlog)
    for (stat in c("mean", "median", "sd")) {
      true_v <- true_vals[[stat]]
      cis <- bootstrap_cis(x, stat, n_boot, confidence)
      err[[stat]]$pct   <- err[[stat]]$pct   + ci_error(true_v, cis$percentile)
      err[[stat]]$basic <- err[[stat]]$basic + ci_error(true_v, cis$basic)
    }
  }
  cat("Studio 10 Problem 2b: log-normal data CI error rates (",
      n_trials, " trials, confidence ", confidence, ")\n", sep = "")
  for (stat in c("mean", "median", "sd"))
    cat(sprintf("  %-7s percentile %.4f   basic %.4f\n",
                stat,
                err[[stat]]$pct   / n_trials,
                err[[stat]]$basic / n_trials))
  cat("Notice that the mean and sd CIs under-cover on right-skewed data,\n")
  cat("while the median CI stays close to nominal. This is exactly the\n")
  cat("under-coverage issue Sections 5 and 14 watch for in real MY2024\n")
  cat("counts when the data are over-dispersed.\n\n")
}
studio10_problem_2b_fish(meanlog = 4.6, sdlog = 0.3,
                          n_data = 30, n_boot = 1000,
                          n_trials = 400, confidence = 0.95)
