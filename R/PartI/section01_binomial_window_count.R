#---------------------------------------------------------
# File:   section01_binomial_window_count.R
# Part I, Section 1 - Binomial window count (adults at the ladder)
#
# Window counts are ADULTS at the ladder; estimator is a_d_hat = w / r.
# There is NO guidance efficiency here -- GE is a smolt-bypass concept
# that enters in Part III.
#
# Repo pointer (escapeLGD/R/night_fall_reascend_wc_binom.R):
#   line 143:  wc_binom <- list(wc %>% mutate(wc = round(wc / wc_prop)))
#   line 150:  wc_binom[[2]][,i] <- rbinom(boots, wc[i], wc_prop) / wc_prop
#
# Place every answer inside the wrapper functions below.
# Run section01_binomial_window_count-test.R after sourcing this file.
#---------------------------------------------------------
# Section 1 ----

#--------------------------------------
# Problem 1a: Simulate one day at the ladder window
section1_problem_1a_fish <- function(a_d_true, r_sample) {
  cat("\n----------------------------------\n")
  cat("Problem 1a: Single day window count and a_d_hat\n")

  # Arguments:
  #   a_d_true = true daytime adult escapement
  #   r_sample = sampling fraction r (e.g. 5/6)

  # Simulate w via rbinom, compute a_d_hat = w / r, and print both.

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 1b: Build escapeLGD-style `wc` tibble for a season
section1_problem_1b_fish <- function(n_weeks, wc_prop, mean_truth) {
  cat("\n----------------------------------\n")
  cat("Problem 1b: Season-long escapeLGD wc tibble (no GE)\n")

  # Arguments:
  #   n_weeks    = number of statistical weeks
  #   wc_prop    = sampling fraction r
  #   mean_truth = mean true weekly passage (for simulating truth)
  #
  # Build a tibble with columns sWeek, truth, wc, where wc ~ Binomial(truth,
  # wc_prop). Return it invisibly. No GE column!

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 1c: Vectorized parametric bootstrap of season total (escapeLGD style)
section1_problem_1c_fish <- function(wc, wc_prop, boots) {
  cat("\n----------------------------------\n")
  cat("Problem 1c: Vectorized parametric bootstrap of season total\n")

  # Arguments:
  #   wc      = tibble with columns sWeek, wc (raw window counts per week)
  #   wc_prop = sampling fraction r
  #   boots   = number of bootstrap iterations
  #
  # For each row i, draw rbinom(boots, wc$wc[i], wc_prop) / wc_prop and
  # sum across rows to get bootstrap season totals. Report point estimate
  # and 95% CI.

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 1d: Replicate to verify a_d_hat = w/r is unbiased
section1_problem_1d_fish <- function(a_d_true, r_sample, nreps) {
  cat("\n----------------------------------\n")
  cat("Problem 1d: Replicated unbiasedness check\n")

  # Arguments:
  #   a_d_true = true daytime adult escapement
  #   r_sample = sampling fraction r
  #   nreps    = number of replicates
  #
  # Replicate (w <- rbinom(1, a_d_true, r_sample); a_d_hat <- w / r_sample)
  # nreps times. Report mean and sd of the estimator and draw a histogram.

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}
