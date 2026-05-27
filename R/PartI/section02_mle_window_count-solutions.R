#---------------------------------------------------------
# File:   section02_mle_window_count-solutions.R
# Part I, Section 2 solutions
#---------------------------------------------------------
set.seed(2026)

# Section 2 solutions ----

library(ggplot2)
library(tibble)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

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
  ggsave(file.path(plots_dir, "section02_loglik.png"), p,
         width = 6, height = 4, dpi = 150)
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


# Problem 2d: Inverse-SR weighting in SCRAPI's thetahat
section2_problem_2d_fish <- function(stocks, strats, SR, n_fish) {
  cat("\n----------------------------------\n")
  cat("Problem 2d: Inverse-SR weighting (SCRAPI thetahat pattern, SMOLT trap)\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  AllPrime <- tibble(
    Strat = sample(strats, size = n_fish, replace = TRUE),
    PGrp  = sample(stocks, size = n_fish, replace = TRUE),
    SR    = SR
  )
  # SCRAPI line 94 reproduced:
  Primarystrata <- tapply(1 / AllPrime$SR,
                          list(factor(AllPrime$Strat, levels = strats),
                               factor(AllPrime$PGrp,  levels = stocks)),
                          sum)
  Primarystrata[is.na(Primarystrata)] <- 0
  Primaryproportions <- prop.table(Primarystrata, margin = 1)

  cat("AllPrime rows =", n_fish, "  SR =", round(SR, 3), "\n")
  cat("Primarystrata (inverse-SR weighted counts):\n")
  print(round(Primarystrata, 1))
  cat("Primaryproportions (row-normalized within stratum):\n")
  print(round(Primaryproportions, 3))
  invisible(list(AllPrime = AllPrime,
                 Primarystrata = Primarystrata,
                 Primaryproportions = Primaryproportions))
}
