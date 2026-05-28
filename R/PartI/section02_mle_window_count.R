#---------------------------------------------------------
# File:   section02_mle_window_count.R
# Part I, Section 2 - Binomial MLE for the window count
#
# Problems 2a-2c: a_d_hat = w / r is the Binomial MLE; verified by
# grid search and continuous relaxation with optim().
# Inverse-SR weighting (formerly 2d) moved to Section 5, problem 5a.
#
# Repo pointer (SCOBI/R/SCRAPI.r):
#   line 75:  dailypass <- passage$Tally / passage$Ptrue    (per-day MLE)
#
# Place every answer inside the wrapper functions below.
# Run section02_mle_window_count-test.R after sourcing this file.
#---------------------------------------------------------
# Section 2 ----

#--------------------------------------
# Problem 2a: Replicate to verify centering of the MLE
section2_problem_2a_fish <- function(a_d_true, r_sample, nreps) {
  cat("\n----------------------------------\n")
  cat("Problem 2a: Replicated MLE centering check\n")

  # Arguments:
  #   a_d_true = true daytime adult escapement
  #   r_sample = sampling fraction r
  #   nreps    = number of replicates
  #
  # Replicate the closed-form MLE nreps times and report mean and sd.

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 2b: Grid evaluation of the log-likelihood for one observed w
section2_problem_2b_fish <- function(w_obs, r_sample, a_d_max) {
  cat("\n----------------------------------\n")
  cat("Problem 2b: Grid evaluation of the binomial log-likelihood\n")

  # Arguments:
  #   w_obs    = one observed window count
  #   r_sample = sampling fraction r
  #   a_d_max  = upper end of the grid for a_d
  #
  # Evaluate dbinom(w_obs, size = a_d, prob = r_sample, log = TRUE) over
  # a_d in seq(max(w_obs, 1), a_d_max). Plot and report the grid MLE.

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 2c: optim() on a continuous relaxation (lgamma)
section2_problem_2c_fish <- function(w_obs, r_sample) {
  cat("\n----------------------------------\n")
  cat("Problem 2c: optim() continuous-relaxation MLE\n")

  # Arguments:
  #   w_obs    = observed window count
  #   r_sample = sampling fraction r
  #
  # dbinom requires integer size; use lgamma to relax over real a_d.
  # Run optim(method = "Brent") on the negative log-likelihood.

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 2d: → moved to Section 5 as problem 5a.
# Inverse-SR weighting belongs alongside thetahat() and the smolt-trap
# machinery it serves.  See section05_nonparametric_bootstrap.R,
# section5_problem_5a_fish().
