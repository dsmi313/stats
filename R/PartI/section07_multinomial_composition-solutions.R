#---------------------------------------------------------
# File:   section07_multinomial_composition-solutions.R
# Part I, Section 7 solutions
#---------------------------------------------------------
set.seed(2026)

# Section 7 solutions ----

library(dplyr)
library(tibble)

#--------------------------------------
# Problem 7a: softMax
section7_problem_7a_fish <- function(x) {
  cat("\n----------------------------------\n")
  cat("Problem 7a: softMax with numerical stability\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  pr <- exp(x - max(x))
  pr / sum(pr)
}


# Problem 7b: PBT_log_likelihood (verbatim)
section7_problem_7b_fish <- function(pGroups, pW, nGroups, nUntag, tagRates) {
  cat("\n----------------------------------\n")
  cat("Problem 7b: PBT_log_likelihood\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  y           <- pW + sum((1 - tagRates) * pGroups)
  pGroups_tag <- pGroups * tagRates
  dmultinom(x = c(nGroups, nUntag),
            prob = c(pGroups_tag, y), log = TRUE)
}


# Problem 7c: PBT_optimllh wrapper
section7_problem_7c_fish <- function(par, nGroups, nUntag, tagRates) {
  par_prob <- section7_problem_7a_fish(par)
  pGroups  <- par_prob[1:(length(par_prob) - 1)]
  pW       <- par_prob[length(par_prob)]
  section7_problem_7b_fish(pGroups = pGroups, pW = pW,
                           nGroups = nGroups, nUntag = nUntag,
                           tagRates = tagRates)
}


# Problem 7d: PBT_expand_calc_MLE
section7_problem_7d_fish <- function(nGroups, nUntag, tagRates) {
  cat("\n----------------------------------\n")
  cat("Problem 7d: PBT_expand_calc_MLE (escapeLGD lines 28-51)\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  K <- length(nGroups)
  propsMLE <- optim(par = rep(1, K + 1),
                    fn  = section7_problem_7c_fish,
                    nGroups = nGroups, nUntag = nUntag,
                    tagRates = tagRates$tagRate,
                    control = list(fnscale = -1, maxit = 1000),
                    method = "BFGS")
  if (propsMLE$convergence != 0) warning("Convergence error in PBT_expand_calc_MLE")
  par_prob <- section7_problem_7a_fish(propsMLE$par)
  res <- tibble(group = c(tagRates$group, "Unassigned"), prop = par_prob)
  cat("MLE composition:\n"); print(res)
  invisible(res)
}


# Problem 7e: PBT_expand_calc (TotEx accounting)
section7_problem_7e_fish <- function(nGroups, nUntag, tagRates) {
  cat("\n----------------------------------\n")
  cat("Problem 7e: PBT_expand_calc (TotEx accounting)\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  pr <- tibble(group   = tagRates$group,
               count   = as.numeric(nGroups),
               tagRate = tagRates$tagRate) |>
    mutate(expand = count / tagRate,
           diff   = expand - count)
  pr <- pr |> bind_rows(tibble(group = "Unassigned",
                               count = nUntag, tagRate = 1,
                               expand = 0, diff = 0))
  if (pr$count[pr$group == "Unassigned"] < sum(pr$diff)) {
    pr$diff   <- pr$count[pr$group == "Unassigned"] *
                 (pr$diff / sum(pr$diff))
    pr$expand <- pr$count + pr$diff
    pr$expand[pr$group == "Unassigned"] <- 0
  } else {
    pr$expand[pr$group == "Unassigned"] <-
      pr$expand[pr$group == "Unassigned"] - sum(pr$diff)
  }
  res <- pr |> mutate(prop = expand / sum(expand)) |> select(group, prop)
  cat("TotEx composition:\n"); print(res)
  invisible(res)
}


# Problem 7f: PBT_breakdown nonparametric bootstrap
section7_problem_7f_fish <- function(nGroups, nUntag, tagRates, boots) {
  cat("\n----------------------------------\n")
  cat("Problem 7f: PBT_breakdown -- nonparametric bootstrap of TotEx\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  values <- c(rep(names(nGroups), nGroups), rep("Unassigned", nUntag))
  out <- matrix(NA_real_, nrow = boots,
                ncol = length(nGroups) + 1L,
                dimnames = list(NULL, c(names(nGroups), "Unassigned")))
  for (b in seq_len(boots)) {
    smp <- sample(values, length(values), replace = TRUE)
    counts <- table(factor(smp, levels = c(names(nGroups), "Unassigned")))
    fit <- section7_problem_7e_fish(counts[names(nGroups)],
                                    counts["Unassigned"], tagRates)
    out[b, fit$group] <- fit$prop
  }
  cis <- apply(out, 2, quantile, c(0.025, 0.975), na.rm = TRUE)
  cat("Bootstrap 95% CIs:\n")
  print(round(cis, 3))
  invisible(out)
}
