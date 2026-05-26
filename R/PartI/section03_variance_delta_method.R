#---------------------------------------------------------
# File:   section03_variance_delta_method.R
# Part I, Section 3 - Variance, delta method, escapeLGD fallback joint MLE
#
# Repo pointers:
#   escapeLGD/R/fallback_reascend_likelihood.R lines 13-15:
#     fallback_log_likelihood(pf, pre_f, dfr, df, dr, dt) =
#       dbinom(dfr, df, pre_f, log=TRUE) + dbinom(dr, dt, pf*pre_f, log=TRUE)
#   escapeLGD/R/night_fall_reascend_wc_binom.R line 110:
#     rbinom(boots, totalPass[i], p_night[i]) / totalPass[i]
#
# Place every answer inside the wrapper functions below.
#---------------------------------------------------------
# Section 3 ----

#--------------------------------------
# Problem 3a: Layer 1 -- only p_n is noisy
section3_problem_3a_fish <- function(d_a, p_n_true, a_d_fixed, nreps) {
  cat("\n----------------------------------\n")
  cat("Problem 3a: Layer 1 -- only p_n noisy\n")

  # Arguments:
  #   d_a       = PIT-tagged adults available to estimate p_n
  #   p_n_true  = true nighttime proportion
  #   a_d_fixed = daytime escapement (held fixed for this layer)
  #   nreps     = number of replicates
  #
  # Draw d_n ~ Binomial(d_a, p_n_true), compute p_n_hat = d_n / d_a,
  # then a_t = a_d_fixed / (1 - p_n_hat). Report mean and sd.

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 3b: Layer 2 -- a_d and p_n both noisy
section3_problem_3b_fish <- function(a_d_true, r_sample, d_a, p_n_true, nreps) {
  cat("\n----------------------------------\n")
  cat("Problem 3b: Layer 2 -- a_d and p_n both noisy\n")

  # Simulate w ~ Binomial(a_d_true, r_sample) and d_n ~ Binomial(d_a, p_n_true)
  # nreps times. Compute a_t = (w / r_sample) / (1 - d_n / d_a).
  # Report mean and sd.

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 3c: Delta-method variance approximation
section3_problem_3c_fish <- function(a_d_true, r_sample, p_n_true, d_a) {
  cat("\n----------------------------------\n")
  cat("Problem 3c: Delta-method variance for a_t = a_d / (1 - p_n)\n")

  # Independence assumed. Use the closed-form:
  #   Var(a_d) = a_d_true * (1 - r_sample) / r_sample^2
  #   Var(p_n) = p_n_true * (1 - p_n_true) / d_a
  #   Var(a_t) ~ (dg/dX)^2 Var(a_d) + (dg/dY)^2 Var(p_n)
  # where g(X, Y) = X / (1 - Y).

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 3d: escapeLGD fallback joint MLE (verbatim port)
section3_problem_3d_fish <- function(pf_true, pre_f_true, dt) {
  cat("\n----------------------------------\n")
  cat("Problem 3d: escapeLGD fallback joint MLE\n")

  # Arguments:
  #   pf_true    = true P(fallback)
  #   pre_f_true = true P(reascend | fallback)
  #   dt         = total ladder ascensions
  #
  # Simulate df ~ Binomial(dt, pf_true), dfr ~ Binomial(df, pre_f_true),
  # dr ~ Binomial(dt, pf_true * pre_f_true). Build:
  #   fallback_log_likelihood(pf, pre_f, dfr, df, dr, dt) =
  #     dbinom(dfr, df, pre_f, log=TRUE) + dbinom(dr, dt, pf*pre_f, log=TRUE)
  # Then run optim() with L-BFGS-B, fnscale = -1, lower = 1e-7, upper = 1 - 1e-7.

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 3e: escapeLGD nighttime-rate bootstrap
section3_problem_3e_fish <- function(d_a, p_n_true, boots) {
  cat("\n----------------------------------\n")
  cat("Problem 3e: escapeLGD nighttime-rate bootstrap\n")

  # Reproduce escapeLGD line 110:
  #   rbinom(boots, totalPass, p_night) / totalPass

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}
