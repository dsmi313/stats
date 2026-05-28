#---------------------------------------------------------
# File:   section02_mle_window_count-solutions.R
# Part I, Section 2 solutions
#---------------------------------------------------------
# Section 2 solutions ----

library(ggplot2)
library(tibble)


#--------------------------------------
# Problem 2a: Replicated MLE centering
section2_problem_2a_fish <- function(a_d_true, r_sample, nreps) {
  cat("\n----------------------------------\n")
  cat("Problem 2a: Replicated MLE centering check\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  estimates <- vapply(seq_len(nreps),
                      function(i) rbinom(1, a_d_true, r_sample) / r_sample,
                      numeric(1))
  cat("a_d_true =", a_d_true, "  r_sample =", r_sample,
      "  nreps =", nreps, "\n")
  cat("mean(a_d_hat) =", mean(estimates),
      "   sd(a_d_hat) =", sd(estimates), "\n")
  invisible(estimates)
}


# Problem 2b: Grid evaluation of the log-likelihood
section2_problem_2b_fish <- function(w_obs, r_sample, a_d_max) {
  cat("\n----------------------------------\n")
  cat("Problem 2b: Grid evaluation of the binomial log-likelihood\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  a_d_grid <- seq(max(w_obs, 1L), a_d_max, by = 1L)
  loglik   <- dbinom(w_obs, size = a_d_grid, prob = r_sample, log = TRUE)
  mle_grid <- a_d_grid[which.max(loglik)]
  cat("w_obs =", w_obs, "  r_sample =", r_sample, "\n")
  cat("Grid MLE         =", mle_grid, "\n")
  cat("Closed-form w/r  =", w_obs / r_sample, "\n")

  p <- ggplot(tibble(a_d = a_d_grid, loglik = loglik),
              aes(a_d, loglik)) +
    geom_line(colour = "steelblue", linewidth = 1) +
    geom_vline(xintercept = w_obs / r_sample,
               colour = "firebrick", linewidth = 1, linetype = "dashed") +
    labs(title = "Section 2 - Binomial log-likelihood for a_d",
         x = expression(a[d]), y = expression(log~L(a[d])))
  invisible(list(a_d_grid = a_d_grid, loglik = loglik, mle = mle_grid))
}


# Problem 2c: optim() on a continuous relaxation
section2_problem_2c_fish <- function(w_obs, r_sample) {
  cat("\n----------------------------------\n")
  cat("Problem 2c: optim() continuous-relaxation MLE\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  neg_loglik_relaxed <- function(a_d, w, r) {
    if (a_d < w) return(Inf)
    -(lgamma(a_d + 1) - lgamma(w + 1) - lgamma(a_d - w + 1) +
      w * log(r) + (a_d - w) * log(1 - r))
  }
  fit <- optim(par = w_obs * 1.1, fn = neg_loglik_relaxed,
               w = w_obs, r = r_sample,
               method = "Brent", lower = w_obs, upper = 5000)
  cat("w_obs =", w_obs, "  r_sample =", r_sample, "\n")
  cat("optim MLE        =", round(fit$par, 2), "\n")
  cat("Closed-form w/r  =", w_obs / r_sample, "\n")
  invisible(fit)
}



# Narrative walkthrough ----
# The wrapper functions above ARE the correct answers.
# What follows is plain runnable code that explains why each answer works
# and shows where the MLE machinery breaks in real data.

set.seed(2026)

# ---- shared parameters (match stub argument names) ----
a_d_true <- 500L   # true daytime adult escapement
r_sample <- 5/6    # fraction of ladder width covered by the window
nreps    <- 1000L  # Monte Carlo repetitions for the centering check
a_d_max  <- 700L   # upper bound for the log-likelihood grid

# ---- Problem 2a: the MLE is w / r — prove it by simulation ----
# Each day a_d_true fish pass through the ladder.  The window counts only the
# fraction r_sample of the ladder width, so the window count w is:
#   w ~ Binomial(n = a_d_true, p = r_sample)
#
# The log-likelihood for a_d given observed w is:
#   log L(a_d | w) = log C(a_d, w) + w log(r) + (a_d - w) log(1 - r)
# Differentiating with respect to a_d and setting to zero gives a_d_hat = w / r.
# Here we verify that formula centres on the truth across many replicates.

estimates <- vapply(seq_len(nreps),
                    function(i) rbinom(1, a_d_true, r_sample) / r_sample,
                    numeric(1))
# rbinom(1, a_d_true, r_sample) draws one window count w.
# Dividing by r_sample is the MLE for that draw.
# vapply repeats this nreps times and collects results into a numeric vector.

cat("Problem 2a\n")
cat("  mean(a_d_hat) =", round(mean(estimates), 1),
    "  sd =", round(sd(estimates), 1), "\n")
# mean ≈ 500 confirms the estimator is unbiased.
# sd quantifies the shot-to-shot noise at this sampling rate.

# ---- Problem 2b: see the log-likelihood surface directly ----
# Rather than trusting the algebra, we evaluate log L at every integer a_d
# from w_obs up to a_d_max.  The peak should land exactly at w / r.

w_obs <- rbinom(1, a_d_true, r_sample)   # one window count for this example

# The grid must start at w_obs: if a_d < w_obs the window count would exceed
# the number of fish that passed, which has probability zero.
a_d_grid <- seq(max(w_obs, 1L), a_d_max, by = 1L)

# dbinom evaluates the Binomial probability mass function at each grid point.
# log = TRUE gives log-probabilities, which are numerically stable for small p.
loglik   <- dbinom(w_obs, size = a_d_grid, prob = r_sample, log = TRUE)

mle_grid <- a_d_grid[which.max(loglik)]   # integer with highest log-likelihood

cat("Problem 2b\n")
cat("  w_obs =", w_obs, "\n")
cat("  Grid MLE    =", mle_grid,         "\n")   # nearest integer to w/r
cat("  w / r       =", w_obs / r_sample, "\n")   # analytic MLE
# They agree to within rounding — the surface peaks exactly where the
# algebra predicts.

# ---- Problem 2c: optim() when closed forms do not exist ----
# dbinom requires integer size, so we cannot differentiate through it with
# standard optimisers.  We relax by replacing the factorial ratio with lgamma,
# which extends log(n!) smoothly to all real n >= 0.
#   lgamma(n + 1) == log(n!)  for integer n
#   lgamma(x + 1) is smooth and differentiable for real x > 0

neg_loglik_relaxed <- function(a_d, w, r) {
  if (a_d < w) return(Inf)   # constraint: can't observe more fish than passed
  # Relaxed log C(a_d, w) = lgamma(a_d+1) - lgamma(w+1) - lgamma(a_d-w+1)
  -(lgamma(a_d + 1) - lgamma(w + 1) - lgamma(a_d - w + 1) +
    w * log(r) + (a_d - w) * log(1 - r))
}
# The negative sign turns maximisation into the minimisation that optim() does.

# method = "Brent" is for 1-D problems: efficient and exact on an interval.
fit <- optim(par = w_obs * 1.1, fn = neg_loglik_relaxed,
             w = w_obs, r = r_sample,
             method = "Brent", lower = w_obs, upper = 5000)

cat("Problem 2c\n")
cat("  optim MLE   =", round(fit$par, 2), "\n")
cat("  w / r       =", w_obs / r_sample,  "\n")
# The continuous optimiser recovers the analytic answer.
# This pattern — lgamma relaxation + Brent — is what escapeLGD uses when
# the likelihood has no closed-form derivative.

# ---- Forward pointer ----
# Section 3 combines the a_d estimator (Sections 1-2) with the nighttime-
# proportion estimator p_n_hat to form the total-escapement estimator
#   a_t = a_d / (1 - p_n)
# and derives its variance via the delta method.
# Inverse-SR weighting (formerly problem 2d) is developed as problem 5a in
# Section 5, where it serves as the explicit motivator for thetahat().
