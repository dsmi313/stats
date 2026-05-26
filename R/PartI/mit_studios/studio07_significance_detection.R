# Studio 7 - Significance testing and probability of hypotheses
# ---------------------------------------------------------------
# Source: https://ocw.mit.edu/courses/18-05-introduction-to-probability-and-statistics-spring-2022/mit18_05_s22_studio7-instructions.pdf
#
# MIT Studio 7 illustrates why P(rejection | H0) is NOT P(H0 | rejection).
# Fish version: a PIT detector's true detection probability for a tagged
# fish is either nominal (theta_H0 = 0.5) or higher (theta_HA > 0.5). We
# test a sample of n_passes detection trials and use a right-sided test.
# A secret prior over (H0, HA) lets us measure both conditional and
# joint error rates.

set.seed(2026)

# Helper: rejection region for a right-sided binomial test ----------------
right_sided_rejection <- function(n_tosses, theta_H0, alpha) {
  # qbinom is the smallest k with P(X <= k) >= 1 - alpha.
  k <- qbinom(1 - alpha, size = n_tosses, prob = theta_H0)
  # Step up by one if the actual significance overshoots alpha.
  while (1 - pbinom(k - 1, n_tosses, theta_H0) > alpha) k <- k + 1
  list(reject_if_at_or_above = k,
       actual_significance   = 1 - pbinom(k - 1, n_tosses, theta_H0))
}

# Problem 1: rejection region, actual significance, power ----------------
studio7_problem_1_fish <- function(theta_HA, alpha, n_passes) {
  stopifnot(theta_HA > 0.5)
  rr      <- right_sided_rejection(n_passes, 0.5, alpha)
  power   <- 1 - pbinom(rr$reject_if_at_or_above - 1, n_passes, theta_HA)
  cat("Studio 7 Problem 1 (n_passes =", n_passes,
      ", alpha =", alpha, ", theta_HA =", theta_HA, "):\n")
  cat("  reject H0 if detections >=", rr$reject_if_at_or_above, "\n")
  cat("  actual significance       =", round(rr$actual_significance, 4), "\n")
  cat("  power                     =", round(power, 4), "\n\n")
  invisible(list(rr = rr, power = power))
}
studio7_problem_1_fish(theta_HA = 0.7, alpha = 0.05, n_passes = 20)
studio7_problem_1_fish(theta_HA = 0.6, alpha = 0.05, n_passes = 20)

# Problem 2: omniscient simulation tracking both error directions --------
studio7_problem_2_fish <- function(theta_HA, alpha, n_passes,
                                   n_trials, secret_prior) {
  rr  <- right_sided_rejection(n_passes, 0.5, alpha)
  k_r <- rr$reject_if_at_or_above

  is_H0       <- runif(n_trials) < secret_prior[1]
  thetas      <- ifelse(is_H0, 0.5, theta_HA)
  detections  <- rbinom(n_trials, size = n_passes, prob = thetas)
  rejected    <- detections >= k_r

  n_reject  <- sum(rejected)
  n_type1   <- sum(rejected &  is_H0)
  n_type2   <- sum(!rejected & !is_H0)

  p_rej_H0  <- n_type1 / max(sum(is_H0), 1)
  p_H0_rej  <- n_type1 / max(n_reject, 1)
  p_rej_HA  <- sum(rejected & !is_H0) / max(sum(!is_H0), 1)
  p_HA_rej  <- sum(rejected & !is_H0) / max(n_reject, 1)
  p_rej     <- mean(rejected)

  cat("Studio 7 Problem 2 (n_trials =", n_trials,
      ", prior =", paste(secret_prior, collapse = ", "), "):\n")
  cat("  rejections        =", n_reject, "\n")
  cat("  type 1 errors     =", n_type1,  "\n")
  cat("  type 2 errors     =", n_type2,  "\n")
  cat("  P(reject | H0)    =", round(p_rej_H0, 4), "\n")
  cat("  P(H0 | reject)    =", round(p_H0_rej, 4), "\n")
  cat("  P(reject | HA)    =", round(p_rej_HA, 4), "\n")
  cat("  P(HA | reject)    =", round(p_HA_rej, 4), "\n")
  cat("  P(reject)         =", round(p_rej, 4), "\n\n")
  invisible(list(p_rej_H0 = p_rej_H0, p_H0_rej = p_H0_rej))
}

# Problem 3a: all H0 -> rejections are all type-1 errors -----------------
studio7_problem_3a_fish <- function(theta_HA, alpha, n_passes, n_trials) {
  cat("Problem 3a: experiment runs only nominal detectors (prior = c(1, 0))\n")
  studio7_problem_2_fish(theta_HA, alpha, n_passes,
                         n_trials, secret_prior = c(1.0, 0.0))
  cat("Every rejection here is a false positive, so P(H0 | reject) = 1.\n")
  cat("Significance alpha caps the FALSE-POSITIVE RATE among H0 detectors,\n")
  cat("not the probability of H0 given a rejection.\n\n")
}
studio7_problem_3a_fish(theta_HA = 0.7, alpha = 0.05,
                         n_passes = 20, n_trials = 5000)

# Problem 3b: all HA -> rejections all correct ---------------------------
studio7_problem_3b_fish <- function(theta_HA, alpha, n_passes, n_trials) {
  cat("Problem 3b: experiment runs only upgraded detectors (prior = c(0, 1))\n")
  studio7_problem_2_fish(theta_HA, alpha, n_passes,
                         n_trials, secret_prior = c(0.0, 1.0))
  cat("All detectors are HA, so every rejection is a true positive and\n")
  cat("P(HA | reject) = 1. Compare to alpha and power above.\n\n")
}
studio7_problem_3b_fish(theta_HA = 0.7, alpha = 0.05,
                         n_passes = 20, n_trials = 5000)

# Problem 3c: cite the distinction ---------------------------------------
studio7_problem_3c_fish <- function() {
  cat("Significance is P(reject | H0). It is a property of the test.\n")
  cat("P(H0 | reject) is a posterior; you need a prior to compute it.\n")
  cat("Frequentist SCRAPI/EASE CIs only report the former; the integrated\n")
  cat("Bayesian model in Section 22 reports the latter.\n\n")
}
studio7_problem_3c_fish()

# Problem 3d: the chant ---------------------------------------------------
studio7_problem_3d_fish <- function() {
  for (i in 1:5)
    cat("THE SIGNIFICANCE IS NOT THE PROBABILITY OF AN ERROR GIVEN REJECTION!\n")
  cat("It is the probability of rejection given H0.\n")
  cat("Frequentists don't compute P(Error | rejection).\n\n")
}
studio7_problem_3d_fish()

# Problem 4 (optional): exact Bayesian inversion -------------------------
studio7_problem_4_fish <- function(theta_HA, alpha, n_passes, prior) {
  rr <- right_sided_rejection(n_passes, 0.5, alpha)
  k_r <- rr$reject_if_at_or_above
  P_rej_H0 <- 1 - pbinom(k_r - 1, n_passes, 0.5)
  P_rej_HA <- 1 - pbinom(k_r - 1, n_passes, theta_HA)
  P_H0_rej <- (P_rej_H0 * prior[1]) /
              (P_rej_H0 * prior[1] + P_rej_HA * prior[2])
  P_HA_rej <- 1 - P_H0_rej
  cat("Problem 4 exact Bayesian inversion (prior = ",
      paste(prior, collapse = ", "), "):\n", sep = "")
  cat("  P(H0 | reject) =", round(P_H0_rej, 4), "\n")
  cat("  P(HA | reject) =", round(P_HA_rej, 4), "\n\n")
}
studio7_problem_4_fish(theta_HA = 0.7, alpha = 0.05, n_passes = 20,
                        prior = c(0.9, 0.1))
studio7_problem_4_fish(theta_HA = 0.7, alpha = 0.05, n_passes = 20,
                        prior = c(0.5, 0.5))
