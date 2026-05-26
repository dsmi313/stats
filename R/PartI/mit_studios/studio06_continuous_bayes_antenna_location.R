# Studio 6 - Discretized continuous Bayesian updating
# ---------------------------------------------------------------
# Source: https://ocw.mit.edu/courses/18-05-introduction-to-probability-and-statistics-spring-2022/mit18_05_s22_studio6-instructions.pdf
#
# MIT Studio 6 uses Cauchy's Lighthouse to teach continuous Bayesian
# updating on a grid. Fish version: a sonar-tagged adult sturgeon at unknown
# river-km theta emits pings in random directions. A detection station on
# the bank, 1 river-km offshore, hears each ping with bank-side projection
# x = theta + tan(phi), where phi is uniform on (-pi/2, pi/2).
#
# Bank-side projections x_i are Cauchy(theta, 1) distributed. Given 15
# ping projections we use a discretized continuous Bayesian update on a
# Uniform prior for theta.

set.seed(2026)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# Problem 0: simulate averaged Normal samples (CLT preview) ---------------
studio6_problem0_fish <- function() {
  n_per_trial <- 9
  ntrials <- 5000
  raw   <- rnorm(ntrials, mean = 10, sd = 6)
  avgs  <- rowMeans(matrix(rnorm(ntrials * n_per_trial, mean = 10, sd = 6),
                           ncol = n_per_trial))
  png(file.path(plots_dir, "studio06_normal_vs_avg.png"),
      width = 900, height = 400)
  par(mfrow = c(1, 2))
  hist(raw,  breaks = 30, freq = FALSE,
       main = "Raw Norm(10, 36)",   col = "lightyellow")
  hist(avgs, breaks = 30, freq = FALSE,
       main = "Average of 9 draws", col = "lightyellow",
       xlim = range(raw))
  par(mfrow = c(1, 1))
  dev.off()
  cat("Averaging shrinks the spread by sqrt(9) = 3 and pulls the\n")
  cat("distribution toward Normal -- classical CLT.\n\n")
}
studio6_problem0_fish()

# Problem 1: Cauchy vs Normal -- introduce fat tails ---------------------
studio6_problem_1_fish <- function() {
  xs <- seq(-6, 6, length.out = 600)
  png(file.path(plots_dir, "studio06_cauchy_vs_normal.png"),
      width = 600, height = 400)
  plot(xs, dnorm(xs), type = "l", lwd = 2, col = "orange",
       ylim = c(0, 0.45),
       xlab = "x", ylab = "density",
       main = "Studio 6 - Standard Normal (orange) vs Cauchy (blue)")
  lines(xs, dcauchy(xs, location = 0, scale = 1),
        col = "steelblue", lwd = 2)
  legend("topright", legend = c("Normal(0, 1)", "Cauchy(0, 1)"),
         col = c("orange", "steelblue"), lwd = 2, bty = "n")
  dev.off()
  cat("Cauchy density falls off as 1/x^2, far slower than Normal's exp(-x^2/2)\n")
  cat("-- this is why averaging Cauchy samples does NOT narrow the spread.\n\n")
}
studio6_problem_1_fish()

# Problem 2: estimate sonar-tag position via Bayesian updating ------------
studio6_problem_2_fish <- function(true_theta = -2.5,
                                   n_pings    = 15,
                                   theta_min  = -10,
                                   theta_max  =  10,
                                   dtheta     = 0.05) {
  # Simulate ping bank-side projections from Cauchy(theta_true, 1)
  data <- rcauchy(n_pings, location = true_theta, scale = 1)

  # Plot the data
  png(file.path(plots_dir, "studio06_ping_data.png"),
      width = 600, height = 300)
  stripchart(data, method = "stack", pch = 19, col = "steelblue",
             xlim = c(theta_min, theta_max),
             xlab = "bank-side ping projection (river-km)",
             main = "Studio 6 - Sonar pings on the river bank")
  dev.off()

  # Discretized uniform prior on theta
  theta_grid <- seq(theta_min, theta_max, by = dtheta)
  prior <- rep(1, length(theta_grid))
  prior <- prior / sum(prior)

  # Iterate Bayesian update
  posteriors <- matrix(0, nrow = n_pings + 1, ncol = length(theta_grid))
  posteriors[1, ] <- prior
  for (i in seq_len(n_pings)) {
    lik <- dcauchy(data[i], location = theta_grid, scale = 1)
    new <- posteriors[i, ] * lik
    posteriors[i + 1, ] <- new / sum(new)
  }

  map_estimates <- theta_grid[apply(posteriors, 1, which.max)]
  final_post    <- posteriors[nrow(posteriors), ]
  final_map     <- theta_grid[which.max(final_post)]

  # Plot prior + every posterior on one set of axes
  png(file.path(plots_dir, "studio06_posteriors.png"),
      width = 700, height = 400)
  plot(theta_grid, prior, type = "n",
       ylim = c(0, max(posteriors)),
       xlab = expression(theta),
       ylab = "discretized density",
       main = "Studio 6 - Prior and posteriors after each ping")
  cols <- colorRampPalette(c("grey80", "steelblue"))(nrow(posteriors))
  for (i in seq_len(nrow(posteriors))) {
    lines(theta_grid, posteriors[i, ], col = cols[i])
  }
  abline(v = true_theta, col = "firebrick", lty = 2, lwd = 2)
  dev.off()

  # Plot MAP trajectory
  png(file.path(plots_dir, "studio06_map_trajectory.png"),
      width = 600, height = 400)
  plot(0:n_pings, map_estimates, type = "b", pch = 19,
       col = "steelblue", lwd = 1,
       xlab = "ping number", ylab = "MAP estimate of theta",
       main = "Studio 6 - MAP estimate trajectory")
  abline(h = true_theta, col = "firebrick", lty = 2)
  dev.off()

  # Plot final posterior
  png(file.path(plots_dir, "studio06_final_posterior.png"),
      width = 600, height = 400)
  plot(theta_grid, final_post, type = "l",
       col = "steelblue", lwd = 2,
       xlab = expression(theta), ylab = "posterior",
       main = paste0("Studio 6 - Final posterior, MAP = ",
                     round(final_map, 2)))
  abline(v = c(true_theta, final_map),
         col = c("firebrick", "darkgreen"), lty = 2, lwd = 2)
  legend("topleft",
         legend = c(paste0("truth = ",     round(true_theta, 2)),
                    paste0("MAP estimate = ", round(final_map, 2))),
         col = c("firebrick", "darkgreen"), lty = 2, lwd = 2, bty = "n")
  dev.off()

  cat("Sonar-tag position estimate:\n")
  cat("  truth         =", true_theta, "\n")
  cat("  final MAP     =", final_map, "\n")
  cat("  search the bank near river-km", round(final_map, 2),
      "for the obscure path.\n\n")
  invisible(list(theta_grid = theta_grid,
                 posteriors  = posteriors,
                 map         = map_estimates))
}
studio6_problem_2_fish()
