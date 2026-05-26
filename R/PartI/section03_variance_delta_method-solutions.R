#---------------------------------------------------------
# File:   section03_variance_delta_method-solutions.R
# Part I, Section 3 solutions
#---------------------------------------------------------
# Section 3 solutions ----

#--------------------------------------
# Problem 3a: Layer 1 -- only p_n noisy
section3_problem_3a_fish <- function(d_a, p_n_true, a_d_fixed, nreps) {
  cat("\n----------------------------------\n")
  cat("Problem 3a: Layer 1 -- only p_n noisy\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  d_n_sim <- rbinom(nreps, size = d_a, prob = p_n_true)
  p_n_hat <- d_n_sim / d_a
  a_t     <- a_d_fixed / (1 - p_n_hat)
  cat("d_a =", d_a, "  p_n_true =", p_n_true,
      "  a_d_fixed =", a_d_fixed, "  nreps =", nreps, "\n")
  cat("mean(a_t) =", mean(a_t), "   sd(a_t) =", sd(a_t), "\n")
  invisible(a_t)
}


# Problem 3b: Layer 2 -- a_d and p_n both noisy
section3_problem_3b_fish <- function(a_d_true, r_sample, d_a, p_n_true, nreps) {
  cat("\n----------------------------------\n")
  cat("Problem 3b: Layer 2 -- a_d and p_n both noisy\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  w_sim   <- rbinom(nreps, size = a_d_true, prob = r_sample)
  d_n_sim <- rbinom(nreps, size = d_a, prob = p_n_true)
  a_d_sim <- w_sim / r_sample
  p_n_hat <- d_n_sim / d_a
  a_t     <- a_d_sim / (1 - p_n_hat)
  cat("a_d_true =", a_d_true, "  r_sample =", r_sample,
      "  d_a =", d_a, "  p_n_true =", p_n_true,
      "  nreps =", nreps, "\n")
  cat("mean(a_t) =", mean(a_t), "   sd(a_t) =", sd(a_t), "\n")
  invisible(a_t)
}


# Problem 3c: Delta-method variance approximation
section3_problem_3c_fish <- function(a_d_true, r_sample, p_n_true, d_a) {
  cat("\n----------------------------------\n")
  cat("Problem 3c: Delta-method variance for a_t = a_d / (1 - p_n)\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  var_a_d <- a_d_true * (1 - r_sample) / r_sample^2
  var_p_n <- p_n_true * (1 - p_n_true) / d_a
  g_dX <- 1 / (1 - p_n_true)
  g_dY <- a_d_true / (1 - p_n_true)^2
  var_a_t_delta <- g_dX^2 * var_a_d + g_dY^2 * var_p_n
  cat("Var(a_d)             =", var_a_d, "\n")
  cat("Var(p_n)             =", var_p_n, "\n")
  cat("Var(a_t) [delta]     =", var_a_t_delta, "\n")
  cat("sd(a_t)  [delta]     =", sqrt(var_a_t_delta), "\n")
  invisible(list(var_a_d = var_a_d, var_p_n = var_p_n,
                 var_a_t_delta = var_a_t_delta))
}


# Problem 3d: escapeLGD fallback joint MLE
section3_problem_3d_fish <- function(pf_true, pre_f_true, dt) {
  cat("\n----------------------------------\n")
  cat("Problem 3d: escapeLGD fallback joint MLE\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  # Simulate data
  df_count  <- rbinom(1, dt, pf_true)
  dfr_count <- rbinom(1, df_count, pre_f_true)
  dr_count  <- rbinom(1, dt, pf_true * pre_f_true)

  # escapeLGD fallback_reascend_likelihood.R lines 13-15 verbatim:
  fallback_log_likelihood <- function(pf, pre_f, dfr, df, dr, dt) {
    dbinom(dfr, df, pre_f,    log = TRUE) +
    dbinom(dr,  dt, pf*pre_f, log = TRUE)
  }
  optimllh <- function(par, dfr, df, dr, dt) {
    fallback_log_likelihood(par[1], par[2], dfr, df, dr, dt)
  }

  opts <- optim(par = c(0.1, 0.9), fn = optimllh,
                dfr = dfr_count, df = df_count,
                dr  = dr_count,  dt = dt,
                control = list(fnscale = -1, maxit = 1000),
                method = "L-BFGS-B",
                upper = 1 - 1e-7, lower = 1e-7)
  cat("Simulated: dt =", dt, "  df =", df_count,
      "  dfr =", dfr_count, "  dr =", dr_count, "\n")
  cat("Truth: pf =", pf_true, "  pre_f =", pre_f_true, "\n")
  cat("MLE:   pf =", round(opts$par[1], 3),
      "  pre_f =", round(opts$par[2], 3), "\n")
  invisible(opts)
}


# Problem 3e: escapeLGD nighttime-rate bootstrap
section3_problem_3e_fish <- function(d_a, p_n_true, boots) {
  cat("\n----------------------------------\n")
  cat("Problem 3e: escapeLGD nighttime-rate bootstrap\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  # escapeLGD night_fall_reascend_wc_binom.R line 110 verbatim:
  p_night_boot <- rbinom(boots, d_a, p_n_true) / d_a
  cat("d_a =", d_a, "  p_n_true =", p_n_true, "  boots =", boots, "\n")
  cat("mean(p_night_boot) =", mean(p_night_boot),
      "   sd =", sd(p_night_boot), "\n")
  invisible(p_night_boot)
}
