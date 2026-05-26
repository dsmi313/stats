# Studio 9 - Simulating confidence intervals (mean fish length)
# ---------------------------------------------------------------
# Source: https://ocw.mit.edu/courses/18-05-introduction-to-probability-and-statistics-spring-2022/mit18_05_s22_studio9-instructions.pdf
#
# MIT Studio 9 simulates type-1 confidence-interval error rates for z- and
# t-CIs, then contrasts them with the Bayesian posterior probability that
# theta is in the CI. Fish version: theta is the true mean smolt length;
# data are Normal(theta, sigma^2).

set.seed(2026)

# Problem 1a: type-1 error rate of z-CIs ---------------------------------
studio9_problem_1a_fish <- function(theta_vals, theta_prior, sigma,
                                    n_data, confidence, n_trials) {
  z_star <- qnorm(1 - (1 - confidence) / 2)
  errors <- 0L
  last_ci <- NULL
  for (t in seq_len(n_trials)) {
    theta <- sample(theta_vals, size = 1, prob = theta_prior)
    x     <- rnorm(n_data, mean = theta, sd = sigma)
    xbar  <- mean(x)
    se    <- sigma / sqrt(n_data)
    ci    <- xbar + c(-1, 1) * z_star * se
    if (theta < ci[1] || theta > ci[2]) errors <- errors + 1L
    last_ci <- ci
  }
  cat("Studio 9 Problem 1a: z-CIs\n")
  cat("  last CI       =", round(last_ci, 3), "\n")
  cat("  type-1 errors =", errors, "/", n_trials,
      "  rate =", round(errors / n_trials, 4),
      "  (expected 1 - c =", round(1 - confidence, 4), ")\n\n")
}
studio9_problem_1a_fish(theta_vals = c(95, 100, 105),
                        theta_prior = c(0.3, 0.4, 0.3),
                        sigma = 10, n_data = 25,
                        confidence = 0.95, n_trials = 5000)

# Problem 1b: type-1 error rate of t-CIs (sigma unknown) -----------------
studio9_problem_1b_fish <- function(theta_vals, theta_prior, sigma,
                                    n_data, confidence, n_trials) {
  errors <- 0L
  last_ci <- NULL
  t_star  <- qt(1 - (1 - confidence) / 2, df = n_data - 1)
  for (t in seq_len(n_trials)) {
    theta <- sample(theta_vals, size = 1, prob = theta_prior)
    x     <- rnorm(n_data, mean = theta, sd = sigma)
    xbar  <- mean(x)
    se    <- sd(x) / sqrt(n_data)
    ci    <- xbar + c(-1, 1) * t_star * se
    if (theta < ci[1] || theta > ci[2]) errors <- errors + 1L
    last_ci <- ci
  }
  cat("Studio 9 Problem 1b: t-CIs (sigma estimated from data)\n")
  cat("  last CI       =", round(last_ci, 3), "\n")
  cat("  type-1 errors =", errors, "/", n_trials,
      "  rate =", round(errors / n_trials, 4), "\n\n")
}
studio9_problem_1b_fish(theta_vals = c(95, 100, 105),
                        theta_prior = c(0.3, 0.4, 0.3),
                        sigma = 10, n_data = 25,
                        confidence = 0.95, n_trials = 5000)

# Problem 1c: Bayesian posterior probability theta is in the CI ----------
studio9_problem_1c_fish <- function(theta_vals, theta_prior, sigma,
                                    n_data, confidence, xbar) {
  # Likelihood of xbar under each candidate theta (using known sigma).
  lik <- dnorm(xbar, mean = theta_vals, sd = sigma / sqrt(n_data))
  post <- theta_prior * lik
  post <- post / sum(post)

  z_star <- qnorm(1 - (1 - confidence) / 2)
  ci     <- xbar + c(-1, 1) * z_star * sigma / sqrt(n_data)
  in_ci  <- theta_vals >= ci[1] & theta_vals <= ci[2]

  p_prior <- sum(theta_prior[in_ci])
  p_post  <- sum(post[in_ci])

  cat("Studio 9 Problem 1c: prior + posterior probability theta in CI\n")
  cat("  CI          =", round(ci, 3), "\n")
  cat("  prior:      "); print(setNames(theta_prior, theta_vals))
  cat("  posterior:  "); print(round(setNames(post, theta_vals), 4))
  cat("  P(theta in CI | prior)     =", round(p_prior, 4), "\n")
  cat("  P(theta in CI | posterior) =", round(p_post, 4), "\n\n")
}
studio9_problem_1c_fish(theta_vals = c(95, 100, 105),
                        theta_prior = c(0.3, 0.4, 0.3),
                        sigma = 10, n_data = 25,
                        confidence = 0.95, xbar = 99)

# Problem 2 (optional): rule-of-thumb 95% poll CI ------------------------
studio9_problem_2_fish <- function(true_theta, n) {
  votes <- rbinom(1, size = n, prob = true_theta)
  p_hat <- votes / n
  margin <- 1 / sqrt(n)  # rule-of-thumb 95% margin
  cat("Studio 9 Problem 2: poll on supporting a hatchery moratorium\n")
  cat("  n =", n, "   p_hat =", round(p_hat, 3),
      "   margin (~95%) =", round(margin, 3), "\n")
  cat("  rule-of-thumb 95% CI =", round(p_hat - margin, 3),
      "to", round(p_hat + margin, 3), "\n\n")
}
studio9_problem_2_fish(true_theta = 0.55, n = 1000)
