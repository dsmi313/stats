#---------------------------------------------------------
# File:   section07_multinomial_composition-solutions.R
# Part I, Section 7 solutions
#---------------------------------------------------------
library(dplyr)

softmax <- function(x) { e <- exp(x - max(x)); e / sum(e) }

# Problem 7a
section7_problem_7a_fish <- function(n = 200, props = c(0.6, 0.3, 0.1)) {
  cat("\n----------------------------------\n")
  cat("Problem 7a: rmultinom — MLE is observed proportions\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  obs <- rmultinom(1, size = n, prob = props)[, 1]
  cat("True props:  ", round(props, 3), "\n")
  cat("MLE (sample):", round(obs / sum(obs), 3), "\n")

  reps <- replicate(1000, {
    x <- rmultinom(1, size = n, prob = props)[, 1]
    x[1] / sum(x)
  })
  hist(reps, breaks = 30, main = "Stock 1 proportion MLE (1,000 sims)",
       xlab = "Estimated proportion", col = "steelblue")
  abline(v = props[1], col = "red", lwd = 2)
  invisible(obs / sum(obs))
}


# Problem 7b
section7_problem_7b_fish <- function(counts = c(120L, 60L, 20L)) {
  cat("\n----------------------------------\n")
  cat("Problem 7b: multinomial log-likelihood; maximize with optim()\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  ll <- function(par, counts) sum(counts * log(softmax(par)))
  K   <- length(counts)
  fit <- optim(rep(0, K), ll, counts = counts,
               control = list(fnscale = -1), method = "BFGS")
  mle_opt  <- softmax(fit$par)
  mle_true <- counts / sum(counts)
  cat("optim MLE:    ", round(mle_opt, 4), "\n")
  cat("counts/n MLE: ", round(mle_true, 4), "\n")
  invisible(mle_opt)
}


# Problem 7c
section7_problem_7c_fish <- function(
    tagged_counts   = c(S1 = 42L, S2 = 35L, S3 = 18L),
    untagged_counts = c(S1 = 80L, S2 = 60L, S3 = 30L, Unassigned = 32L)) {

  cat("\n----------------------------------\n")
  cat("Problem 7c: product of two multinomials — skeleton of Delomas & Hess\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  K <- length(tagged_counts)
  joint_ll <- function(par, tc, uc) {
    p_tag  <- softmax(par[1:K])
    p_full <- softmax(c(par[1:K], par[K + 1]))
    sum(tc * log(p_tag)) + sum(uc * log(p_full))
  }
  fit <- optim(rep(0, K + 1), joint_ll,
               tc = tagged_counts, uc = untagged_counts,
               control = list(fnscale = -1), method = "BFGS")
  p_stocks <- softmax(fit$par[1:K])
  p_full   <- softmax(c(fit$par[1:K], fit$par[K + 1]))
  cat("Stock proportions:\n")
  print(round(setNames(p_stocks, names(tagged_counts)), 4))
  cat("Full composition (incl. Unassigned):\n")
  print(round(setNames(p_full, c(names(tagged_counts), "Unassigned")), 4))
  invisible(list(p_stocks = p_stocks, p_full = p_full))
}


# Problem 7d
section7_problem_7d_fish <- function(counts = c(S1 = 120L, S2 = 60L, S3 = 20L),
                                      boots  = 10000L) {
  cat("\n----------------------------------\n")
  cat("Problem 7d: nonparametric bootstrap CI for stock proportions\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  fish   <- rep(names(counts), counts)
  stocks <- names(counts)
  boot_mat <- matrix(NA_real_, nrow = boots, ncol = length(counts))
  colnames(boot_mat) <- stocks
  for (b in seq_len(boots)) {
    smp <- sample(fish, length(fish), replace = TRUE)
    tbl <- table(factor(smp, levels = stocks))
    boot_mat[b, ] <- as.numeric(tbl) / sum(tbl)
  }
  cis <- apply(boot_mat, 2, quantile, c(0.025, 0.975))
  cat("Bootstrap 95% CIs:\n"); print(round(cis, 3))
  cat("Observed proportions:\n"); print(round(counts / sum(counts), 3))
  invisible(cis)
}
