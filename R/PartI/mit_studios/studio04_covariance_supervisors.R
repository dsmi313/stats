# Studio 4 - Covariance and correlation (two trap supervisors)
# ---------------------------------------------------------------
# Source: https://ocw.mit.edu/courses/18-05-introduction-to-probability-and-statistics-spring-2022/mit18_05_s22_studio4-instructions.pdf
#
# MIT Studio 4 uses two casino dealers (Axel, Barto) playing roulette
# together and one playing alone, then measures covariance and correlation
# of their daily winnings.
#
# Fish version: two trap supervisors (A and B) covering the LWG juvenile
# bypass. Each "bet" is one hour of trap coverage; the outcome is +1 if at
# least one PIT-tagged fish is detected that hour, else -1 (relative to a
# baseline expectation). Hours covered jointly contribute to both
# supervisors' daily totals -> covariance.

set.seed(2026)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# Per-hour signed outcome: +1 if a tagged fish is detected (prob p_detect),
# else -1.
hour_outcome <- function(n_hours, p_detect = 18/38) {
  ifelse(rbinom(n_hours, size = 1, prob = p_detect) == 1, 1, -1)
}

# Problem 1a: covariance / correlation of A and B daily totals ------------
studio4_problem_1a_fish <- function(n_together, n_B_alone, ntrials) {
  totals_A <- numeric(ntrials)
  totals_B <- numeric(ntrials)
  for (i in seq_len(ntrials)) {
    shared <- hour_outcome(n_together)
    alone  <- hour_outcome(n_B_alone)
    totals_A[i] <- sum(shared)              # A only works shared hours
    totals_B[i] <- sum(shared) + sum(alone) # B works shared + alone hours
  }
  cat("Supervisor A and B daily totals (n_together =", n_together,
      ", n_B_alone =", n_B_alone, "):\n")
  cat("  E[A] =", round(mean(totals_A), 3),
      "   Var[A] =", round(var(totals_A), 3), "\n")
  cat("  E[B] =", round(mean(totals_B), 3),
      "   Var[B] =", round(var(totals_B), 3), "\n")
  cat("  Cov(A, B)  =", round(cov(totals_A, totals_B), 3), "\n")
  cat("  Cor(A, B)  =", round(cor(totals_A, totals_B), 3), "\n\n")
  invisible(list(A = totals_A, B = totals_B))
}

studio4_problem_1a_fish(n_together = 10, n_B_alone =  0, ntrials = 5000)
studio4_problem_1a_fish(n_together = 10, n_B_alone = 10, ntrials = 5000)
studio4_problem_1a_fish(n_together = 10, n_B_alone = 50, ntrials = 5000)

# Problem 1b: explain the behavior --------------------------------------
studio4_problem_1b_fish <- function() {
  cat("As B works more hours alone:\n")
  cat("  Var(B) grows linearly with n_B_alone (independent hours add up).\n")
  cat("  Cov(A, B) stays at n_together * Var(per-hour outcome), so the\n")
  cat("  correlation drops because the shared portion of B's variance\n")
  cat("  becomes a smaller fraction of the total. This is the same reason\n")
  cat("  Section 6 stratification widens the marginal CI: independent\n")
  cat("  hours add variance without adding shared signal.\n\n")
}
studio4_problem_1b_fish()

# Problem 2: CLT on the day's total signed detections --------------------
studio4_problem_2_fish <- function(n_hours_per_trial, ntrials) {
  totals <- replicate(ntrials, sum(hour_outcome(n_hours_per_trial)))
  p_detect <- 18/38
  per_hour_mean <- 2 * p_detect - 1            # E[+1*p + (-1)*(1-p)]
  per_hour_var  <- 1 - per_hour_mean^2         # Var of +-1 outcome
  mu    <- n_hours_per_trial * per_hour_mean
  sd_x  <- sqrt(n_hours_per_trial * per_hour_var)

  png(file.path(plots_dir, "studio04_clt_supervisor.png"),
      width = 600, height = 400)
  hist(totals, breaks = 30, freq = FALSE,
       col = "lightyellow", border = "steelblue",
       main = "Studio 4 - CLT on supervisor daily total",
       xlab = "signed total detections over the day")
  xs <- seq(min(totals), max(totals), length.out = 200)
  lines(xs, dnorm(xs, mean = mu, sd = sd_x),
        col = "firebrick", lwd = 2)
  dev.off()

  cat("CLT check for daily total (n_hours_per_trial =",
      n_hours_per_trial, "):\n")
  cat("  empirical mean =", round(mean(totals), 3),
      "  (target =", round(mu, 3), ")\n")
  cat("  empirical sd   =", round(sd(totals), 3),
      "  (target =", round(sd_x, 3), ")\n\n")
  invisible(totals)
}
studio4_problem_2_fish(n_hours_per_trial = 100, ntrials = 5000)
