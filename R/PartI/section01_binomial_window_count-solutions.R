#---------------------------------------------------------
# File:   section01_binomial_window_count-solutions.R
# Part I, Section 1 solutions
#---------------------------------------------------------
# Section 1 solutions ----

library(ggplot2)
library(dplyr)
library(tibble)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

#--------------------------------------
# Problem 1a: Simulate one day at the ladder window
section1_problem_1a_fish <- function(a_d_true, r_sample) {
  cat("\n----------------------------------\n")
  cat("Problem 1a: Single day window count and a_d_hat\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  w_obs   <- rbinom(1, size = a_d_true, prob = r_sample)
  a_d_hat <- w_obs / r_sample

  cat("a_d_true =", a_d_true, "  r_sample =", r_sample, "\n")
  cat("Observed window count w =", w_obs, "\n")
  cat("a_d_hat = w / r          =", a_d_hat, "\n")
  invisible(list(w = w_obs, a_d_hat = a_d_hat))
}


# Problem 1b: Build escapeLGD-style `wc` tibble for a season
section1_problem_1b_fish <- function(n_weeks, wc_prop, mean_truth) {
  cat("\n----------------------------------\n")
  cat("Problem 1b: Season-long escapeLGD wc tibble (no GE)\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  wc <- tibble(
    sWeek = seq_len(n_weeks),
    truth = round(runif(n_weeks, mean_truth * 0.5, mean_truth * 1.5)),
    wc    = NA_integer_
  )
  wc$wc <- rbinom(n_weeks, size = wc$truth, prob = wc_prop)

  cat("n_weeks =", n_weeks, "  wc_prop =", wc_prop,
      "  mean_truth =", mean_truth, "\n")
  cat("First rows of wc tibble:\n")
  print(head(wc, 6))
  invisible(wc)
}


# Problem 1c: Vectorized parametric bootstrap of season total
section1_problem_1c_fish <- function(wc, wc_prop, boots) {
  cat("\n----------------------------------\n")
  cat("Problem 1c: Vectorized parametric bootstrap of season total\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  # escapeLGD line 150 reproduced:
  boot_mat <- vapply(seq_len(nrow(wc)),
                     function(i) rbinom(boots, wc$wc[i], wc_prop) / wc_prop,
                     numeric(boots))
  season_totals <- rowSums(boot_mat)
  ci <- quantile(season_totals, c(0.025, 0.975))

  cat("boots =", boots, "  wc_prop =", wc_prop, "\n")
  cat("Point estimate (sum(wc / r)) =", sum(wc$wc / wc_prop), "\n")
  cat("95% bootstrap CI             = [", round(ci[1]), ",",
      round(ci[2]), "]\n")
  invisible(list(season_totals = season_totals, ci = ci))
}


# Problem 1d: Replicate to verify a_d_hat = w/r is unbiased
section1_problem_1d_fish <- function(a_d_true, r_sample, nreps) {
  cat("\n----------------------------------\n")
  cat("Problem 1d: Replicated unbiasedness check\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  sims <- vapply(seq_len(nreps),
                 function(i) rbinom(1, a_d_true, r_sample) / r_sample,
                 numeric(1))
  cat("a_d_true =", a_d_true, "  r_sample =", r_sample,
      "  nreps =", nreps, "\n")
  cat("mean(a_d_hat) =", mean(sims),
      "   sd(a_d_hat) =", sd(sims), "\n")

  p <- ggplot(tibble(a_d_hat = sims), aes(a_d_hat)) +
    geom_histogram(bins = 30, fill = "steelblue", colour = "white") +
    geom_vline(xintercept = a_d_true, colour = "firebrick", linewidth = 1) +
    labs(title = "Section 1 - Window-count estimator a_d_hat = w / r",
         x = expression(hat(a)[d]), y = "Replicates")
  ggsave(file.path(plots_dir, "section01_estimator_hist.png"), p,
         width = 6, height = 4, dpi = 150)
  invisible(sims)
}
