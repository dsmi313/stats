#---------------------------------------------------------
# File:   section07_multinomial_composition.R
# Part I, Section 7 - escapeLGD's PBT composition likelihood
#
# Repo pointer (escapeLGD/R/composition_estimation_utils.R):
#   lines 15-18:   softMax(x)
#   lines 62-66:   PBT_log_likelihood(pGroups, pW, nGroups, nUntag, tagRates)
#                    y <- pW + sum((1 - tagRates) * pGroups)
#                    pGroups_tag <- pGroups * tagRates
#                    dmultinom(x = c(nGroups, nUntag),
#                              prob = c(pGroups_tag, y), log = TRUE)
#   lines 77-84:   PBT_optimllh(par, nGroups, nUntag, tagRates)
#   lines 28-51:   PBT_expand_calc_MLE(values, tagRates)
#   lines 440-455: PBT_expand_calc (accounting / TotEx method)
#   lines 466-486: PBT_breakdown (nonparametric bootstrap)
#---------------------------------------------------------
# Section 7 ----

#--------------------------------------
# Problem 7a: softMax with numerical stability (escapeLGD lines 15-18)
section7_problem_7a_fish <- function(x) {
  cat("\n----------------------------------\n")
  cat("Problem 7a: softMax with numerical stability\n")

  # Reproduce escapeLGD lines 15-18:
  #   pr <- exp(x - max(x))
  #   pr / sum(pr)

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 7b: PBT_log_likelihood (escapeLGD lines 62-66 verbatim)
section7_problem_7b_fish <- function(pGroups, pW, nGroups, nUntag, tagRates) {
  cat("\n----------------------------------\n")
  cat("Problem 7b: PBT_log_likelihood\n")

  # Reproduce escapeLGD lines 62-66 verbatim:
  #   y           <- pW + sum((1 - tagRates) * pGroups)
  #   pGroups_tag <- pGroups * tagRates
  #   dmultinom(x = c(nGroups, nUntag),
  #             prob = c(pGroups_tag, y), log = TRUE)

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 7c: PBT_optimllh wrapper (escapeLGD lines 77-84)
section7_problem_7c_fish <- function(par, nGroups, nUntag, tagRates) {
  cat("\n----------------------------------\n")
  cat("Problem 7c: PBT_optimllh -- softmax wrapper for optim()\n")

  # par_prob <- softMax(par)
  # pGroups  <- par_prob[1:(length(par_prob) - 1)]
  # pW       <- par_prob[length(par_prob)]
  # PBT_log_likelihood(pGroups, pW, nGroups, nUntag, tagRates)
  # (Use the functions you implemented in 7a and 7b.)

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 7d: PBT_expand_calc_MLE -- run optim() to recover proportions
section7_problem_7d_fish <- function(nGroups, nUntag, tagRates) {
  cat("\n----------------------------------\n")
  cat("Problem 7d: PBT_expand_calc_MLE (escapeLGD lines 28-51)\n")

  # K <- length(nGroups)
  # optim(par = rep(1, K + 1),
  #       fn = section7_problem_7c_fish,
  #       nGroups = nGroups, nUntag = nUntag, tagRates = tagRates$tagRate,
  #       control = list(fnscale = -1, maxit = 1000),
  #       method = "BFGS")
  # Then apply softMax to optim's par and return a tibble with
  # group = c(tagRates$group, "Unassigned") and prop = par_prob.

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 7e: PBT_expand_calc -- accounting/TotEx method (escapeLGD lines 440-455)
section7_problem_7e_fish <- function(nGroups, nUntag, tagRates) {
  cat("\n----------------------------------\n")
  cat("Problem 7e: PBT_expand_calc (TotEx accounting)\n")

  # For each hatchery group i: expand[i] = nGroups[i] / tagRates[i].
  # Subtract the overage (sum(expand - nGroups)) from nUntag. Normalize.

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 7f: PBT_breakdown -- nonparametric bootstrap (escapeLGD lines 466-486)
section7_problem_7f_fish <- function(nGroups, nUntag, tagRates, boots) {
  cat("\n----------------------------------\n")
  cat("Problem 7f: PBT_breakdown -- nonparametric bootstrap of TotEx\n")

  # Build the value vector values <- c(rep(names(nGroups), nGroups),
  #                                     rep("Unassigned", nUntag)).
  # For each bootstrap iteration: sample(values, replace = TRUE), recompute
  # counts, call section7_problem_7e_fish. Return a matrix of bootstrap
  # proportions (rows = iterations).

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}
