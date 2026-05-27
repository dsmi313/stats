#---------------------------------------------------------
# File:   section06_stratification-solutions.R
# Part I, Section 6 solutions
#---------------------------------------------------------
set.seed(2026)

# Section 6 solutions ----

library(dplyr)
library(tibble)

#--------------------------------------
# Problem 6a: Simulate a season with shifting p_n
section6_problem_6a_fish <- function(n_weeks, days_per_wk,
                                     p_n_early, p_n_late,
                                     a_d_per_day = 100L,
                                     d_a_per_day = 30L) {
  cat("\n----------------------------------\n")
  cat("Problem 6a: Simulate a season with shifting p_n\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  ndays  <- n_weeks * days_per_wk
  season <- tibble(
    Week    = rep(seq_len(n_weeks), each = days_per_wk),
    Collaps = rep(c("EARLY", "LATE"), each = ndays / 2),
    p_n_true = ifelse(rep(c("EARLY", "LATE"), each = ndays / 2) == "EARLY",
                      p_n_early, p_n_late)
  )
  season <- season |>
    mutate(a_d      = a_d_per_day,
           d_n      = rbinom(n(), size = d_a_per_day, prob = p_n_true),
           p_n_hat  = d_n / d_a_per_day,
           a_t_true = a_d / (1 - p_n_true))
  cat("n_weeks =", n_weeks, "  ndays =", ndays,
      "  truth =", sum(season$a_t_true), "\n")
  invisible(season)
}


# Problem 6b: Cpattern + pooled vs stratified estimator
section6_problem_6b_fish <- function(season) {
  cat("\n----------------------------------\n")
  cat("Problem 6b: Cpattern + pooled vs stratified estimator\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  Cpattern <- unique(cbind(season$Week, season$Collaps))
  cat("Cpattern (SCRAPI line 217):\n"); print(Cpattern)

  d_a_per_day <- 30L
  p_n_pooled  <- sum(season$d_n) / (nrow(season) * d_a_per_day)
  a_t_pooled  <- sum(season$a_d) / (1 - p_n_pooled)

  strat_summary <- season |>
    group_by(Collaps) |>
    summarise(p_n_hat = sum(d_n) / (n() * d_a_per_day),
              a_d_sum = sum(a_d), .groups = "drop") |>
    mutate(a_t = a_d_sum / (1 - p_n_hat))
  a_t_stratified <- sum(strat_summary$a_t)

  truth_total <- sum(season$a_t_true)
  cat("Truth        =", round(truth_total), "\n")
  cat("Pooled       =", round(a_t_pooled),
      "   bias =", round(a_t_pooled - truth_total), "\n")
  cat("Stratified   =", round(a_t_stratified),
      "   bias =", round(a_t_stratified - truth_total), "\n")
  invisible(list(Cpattern = Cpattern,
                 a_t_pooled = a_t_pooled,
                 a_t_stratified = a_t_stratified,
                 truth_total = truth_total))
}


# Problem 6c: SCRAPI's Collaps assignment loop
section6_problem_6c_fish <- function(pass, All) {
  cat("\n----------------------------------\n")
  cat("Problem 6c: Collaps assignment by date match (SCRAPI lines 254-259)\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  nAll <- nrow(All)
  All$Collaps <- character(nAll)
  AllDates    <- unique(All$CollectionDate)
  for (d in AllDates) {
    CollStrat <- pass$Collaps[which(pass$SampleEndDate == d)]
    All$Collaps[which(All$CollectionDate == d)] <- CollStrat
  }
  cat("Collaps assigned to", nAll, "fish records.\n")
  print(head(All, 5))
  invisible(All)
}


# Problem 6d: Trigger SCRAPI Error 3
section6_problem_6d_fish <- function(pass, All) {
  cat("\n----------------------------------\n")
  cat("Problem 6d: Trigger SCRAPI Error 3 with duplicated FPC date\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  pass_bad <- bind_rows(pass, pass[1, ])
  AllDates <- unique(All$CollectionDate)
  caught <- tryCatch({
    for (d in AllDates) {
      CollStrat <- pass_bad$Collaps[which(pass_bad$SampleEndDate == d)]
      All$Collaps[which(All$CollectionDate == d)] <- CollStrat
    }
    "no error fired"
  }, warning = function(w) conditionMessage(w),
     error   = function(e) conditionMessage(e))
  cat("Caught:", caught, "\n")
  cat("This is SCRAPI Error 3 (Collaps recycling).\n")
  cat("Fix: 1:1 date correspondence between trap and FPC data.\n")
  invisible(caught)
}


# Problem 6e: Parametric-bootstrap CI: pooled vs stratified
section6_problem_6e_fish <- function(season, B) {
  cat("\n----------------------------------\n")
  cat("Problem 6e: Parametric-bootstrap CI for pooled vs stratified\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  d_a <- 30L
  boot_season <- function(stratified) {
    estimates <- numeric(B)
    if (stratified) {
      p_strata <- season |>
        group_by(Collaps) |>
        summarise(p_n_hat = sum(d_n) / (n() * d_a), .groups = "drop")
      for (b in seq_len(B)) {
        smpl <- season |>
          left_join(p_strata, by = "Collaps", suffix = c("", ".s")) |>
          mutate(d_n_b = rbinom(n(), size = d_a, prob = p_n_hat.s))
        strat_b <- smpl |>
          group_by(Collaps) |>
          summarise(p_n_b   = sum(d_n_b) / (n() * d_a),
                    a_d_sum = sum(a_d), .groups = "drop")
        estimates[b] <- sum(strat_b$a_d_sum / (1 - strat_b$p_n_b))
      }
    } else {
      p_pool <- sum(season$d_n) / (nrow(season) * d_a)
      for (b in seq_len(B)) {
        d_n_b <- rbinom(nrow(season), size = d_a, prob = p_pool)
        p_n_b <- sum(d_n_b) / (nrow(season) * d_a)
        estimates[b] <- sum(season$a_d) / (1 - p_n_b)
      }
    }
    quantile(estimates, c(0.025, 0.975))
  }
  ci_pool   <- boot_season(FALSE)
  ci_strat  <- boot_season(TRUE)
  cat("B =", B, "\n")
  cat("Pooled       CI = [", round(ci_pool[1]),  ",",
      round(ci_pool[2]),  "]   width =", round(diff(ci_pool)),  "\n")
  cat("Stratified   CI = [", round(ci_strat[1]), ",",
      round(ci_strat[2]), "]   width =", round(diff(ci_strat)), "\n")
  invisible(list(ci_pool = ci_pool, ci_strat = ci_strat))
}
