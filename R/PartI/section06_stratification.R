#---------------------------------------------------------
# File:   section06_stratification.R
# Part I, Section 6 - Cpattern, Collaps assignment, SCRAPI Error 3
#
# Repo pointers (SCOBI/R/SCRAPI.r):
#   line 217:      Cpattern <- unique(cbind(pass$Week, pass$Collaps))
#   lines 254-259: Collaps assignment loop (Error 3 trigger)
#---------------------------------------------------------
# Section 6 ----

#--------------------------------------
# Problem 6a: Simulate a season with shifting p_n
section6_problem_6a_fish <- function(n_weeks, days_per_wk,
                                     p_n_early, p_n_late,
                                     a_d_per_day = 100L,
                                     d_a_per_day = 30L) {
  cat("\n----------------------------------\n")
  cat("Problem 6a: Simulate a season with shifting p_n\n")

  # Build a tibble with columns Week, Collaps (EARLY/LATE),
  # p_n_true, a_d, d_n, p_n_hat, a_t_true. Return it.

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 6b: Cpattern, pooled vs stratified estimator
section6_problem_6b_fish <- function(season) {
  cat("\n----------------------------------\n")
  cat("Problem 6b: Cpattern + pooled vs stratified estimator\n")

  # Build Cpattern = unique(cbind(season$Week, season$Collaps)).
  # Compute pooled p_n_pooled = sum(d_n) / (nrow*d_a), a_t_pooled.
  # Compute stratified p_n_hat per Collaps, a_t per Collaps, sum.
  # Report both and the truth_total.

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 6c: SCRAPI's Collaps assignment loop (verbatim port)
section6_problem_6c_fish <- function(pass, All) {
  cat("\n----------------------------------\n")
  cat("Problem 6c: Collaps assignment by date match (SCRAPI lines 254-259)\n")

  # Reproduce SCRAPI lines 254-259 verbatim:
  #   nAll <- nrow(All)
  #   All$Collaps <- character(nAll)
  #   AllDates    <- unique(All$CollectionDate)
  #   for (d in AllDates) {
  #     CollStrat <- pass$Collaps[which(pass$SampleEndDate == d)]
  #     All$Collaps[which(All$CollectionDate == d)] <- CollStrat
  #   }
  # Return the updated All tibble.

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 6d: Trigger SCRAPI Error 3 with a duplicated FPC date
section6_problem_6d_fish <- function(pass, All) {
  cat("\n----------------------------------\n")
  cat("Problem 6d: Trigger SCRAPI Error 3 with duplicated FPC date\n")

  # Duplicate the first row of pass, then rerun the Collaps loop. Catch the
  # warning / error from R recycling (CollStrat has length 2).

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 6e: Parametric-bootstrap CI: pooled vs stratified
section6_problem_6e_fish <- function(season, B) {
  cat("\n----------------------------------\n")
  cat("Problem 6e: Parametric-bootstrap CI for pooled vs stratified\n")

  # For each of B iterations, draw d_n_b ~ Binomial(d_a, p_n_hat). Compute
  # pooled and stratified a_t. Report both 95% CIs and CI widths.

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}
