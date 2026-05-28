#---------------------------------------------------------
# File:   section06_stratification-solutions.R
# Part I, Section 6 solutions
#---------------------------------------------------------
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


# Narrative walkthrough ----
# The wrapper functions above ARE the correct answers.
# What follows is plain runnable code that explains why pooling p_n across
# strata introduces bias, traces the Collaps assignment loop, and shows
# exactly how SCRAPI Error 3 fires.

library(tibble)
library(dplyr)
set.seed(2026)

# ---- shared parameters (match stub argument names and test values) ----
n_weeks    <- 8L      # total weeks in the season
days_per_wk <- 7L     # days sampled per week
p_n_early  <- 0.10    # nighttime fraction in early season (low — fish run in daylight)
p_n_late   <- 0.40    # nighttime fraction in late season (fish shift to night)
a_d_per_day <- 100L   # fixed daytime passage per day
d_a_per_day <- 30L    # PIT-tagged adults available each day to estimate p_n
B          <- 1500L   # bootstrap iterations for 6e

# ---- Problem 6a: simulate a season with a structural shift in p_n ----
# The season is split into EARLY and LATE strata.  p_n is stable within each
# stratum but very different between them — the classic stratification scenario.

ndays  <- n_weeks * days_per_wk
season <- tibble(
  Week    = rep(seq_len(n_weeks), each = days_per_wk),
  # First half of the season = EARLY, second half = LATE.
  Collaps = rep(c("EARLY", "LATE"), each = ndays / 2),
  p_n_true = ifelse(rep(c("EARLY", "LATE"), each = ndays / 2) == "EARLY",
                    p_n_early, p_n_late)
)
season <- season |>
  mutate(
    a_d     = a_d_per_day,
    # d_n ~ Binomial(d_a_per_day, p_n_true): count of tagged fish moving at night.
    d_n     = rbinom(n(), size = d_a_per_day, prob = p_n_true),
    p_n_hat = d_n / d_a_per_day,          # MLE for p_n on that day
    a_t_true = a_d / (1 - p_n_true)       # total escapement if p_n were known exactly
  )
cat("Problem 6a\n")
cat("  Season truth =", round(sum(season$a_t_true)), "\n")

# ---- Problem 6b: Cpattern — the map from weeks to strata ----
# SCRAPI line 217 extracts a unique Week → Collaps lookup table.
# Every week belongs to exactly one stratum; Cpattern makes that explicit.

Cpattern <- unique(cbind(season$Week, season$Collaps))
cat("Problem 6b: Cpattern\n")
print(Cpattern)

# Pooled estimator: compute a single p_n_pooled by ignoring the EARLY/LATE split.
# All d_n observations are mixed together.
p_n_pooled <- sum(season$d_n) / (nrow(season) * d_a_per_day)
a_t_pooled <- sum(season$a_d) / (1 - p_n_pooled)
# Because LATE has much higher p_n, the pooled estimate is pulled upward,
# causing the denominator (1 - p_n_pooled) to be too small and a_t to be inflated.

# Stratified estimator: estimate p_n separately within each Collaps group,
# compute a_t for each stratum, then add them up.
strat_summary <- season |>
  group_by(Collaps) |>
  summarise(p_n_hat = sum(d_n) / (n() * d_a_per_day),
            a_d_sum = sum(a_d), .groups = "drop") |>
  mutate(a_t = a_d_sum / (1 - p_n_hat))
a_t_stratified <- sum(strat_summary$a_t)
truth_total    <- sum(season$a_t_true)

cat("  Truth      =", round(truth_total), "\n")
cat("  Pooled     =", round(a_t_pooled),
    "  bias =", round(a_t_pooled - truth_total), "\n")
cat("  Stratified =", round(a_t_stratified),
    "  bias =", round(a_t_stratified - truth_total), "\n")
# Pooled bias is systematically positive when p_n is higher in the strata with
# more fish (Jensen's inequality strikes again through the 1/(1-p_n) nonlinearity).

# ---- Build pass and All for problems 6c and 6d ----
# These mirror what the test driver constructs.

nFish      <- 200L
fish_weeks <- sample(season$Week, size = nFish, replace = TRUE)
All <- tibble(
  CollectionDate = format(as.Date("2024-04-01") + (fish_weeks - 1L) * 7L,
                          "%m/%d/%Y"),
  GenStock       = sample(c("LOSALM", "CHMBLN", "IMNAHA"),
                          size = nFish, replace = TRUE)
)
# All: individual fish PIT/genetic records with a collection date.

pass <- season |>
  mutate(SampleEndDate = format(as.Date("2024-04-01") + (Week - 1L) * 7L,
                                "%m/%d/%Y")) |>
  group_by(Week, Collaps, SampleEndDate) |>
  summarise(SampleCount = sum(a_d), .groups = "drop")
# pass: one row per week summarising the FPC (fish passage counter) trap data.
# SampleEndDate is used to join fish records back to their weekly Collaps stratum.

# ---- Problem 6c: Collaps assignment loop (SCRAPI lines 254-259) ----
# Each fish in All needs to know which Collaps stratum it belongs to.
# SCRAPI matches on date: for each unique CollectionDate in All, look up the
# corresponding Collaps in pass via SampleEndDate.

nAll     <- nrow(All)
All$Collaps <- character(nAll)    # pre-allocate the new column
AllDates <- unique(All$CollectionDate)
for (d in AllDates) {
  CollStrat <- pass$Collaps[which(pass$SampleEndDate == d)]
  # which() finds the row(s) in pass whose SampleEndDate matches date d.
  # CollStrat should be exactly one value — one stratum per date.
  All$Collaps[which(All$CollectionDate == d)] <- CollStrat
  # Assign that stratum to every fish collected on date d.
}
cat("Problem 6c: Collaps assigned to", nAll, "fish records\n")
cat("  Unique Collaps values assigned:", unique(All$Collaps), "\n")

# ---- Problem 6d: trigger SCRAPI Error 3 ----
# Error 3 fires when a single CollectionDate matches MORE THAN ONE row in pass
# (a duplicated FPC date).  CollStrat then has length > 1, and the assignment
# to All$Collaps recycles or errors because the lengths don't match.

pass_bad <- bind_rows(pass, pass[1, ])   # duplicate the first FPC date
caught <- tryCatch({
  for (d in AllDates) {
    CollStrat <- pass_bad$Collaps[which(pass_bad$SampleEndDate == d)]
    # On the duplicated date, CollStrat has length 2 instead of 1.
    All$Collaps[which(All$CollectionDate == d)] <- CollStrat
    # R tries to assign 2 values to a vector of length 1 — triggers a warning.
  }
  "no error fired"
}, warning = function(w) conditionMessage(w),
   error   = function(e) conditionMessage(e))
cat("Problem 6d caught:", caught, "\n")
cat("  Fix: ensure 1:1 date correspondence between trap (All) and FPC (pass).\n")

# ---- Problem 6e: bootstrap CI — pooled vs stratified ----
# Parametric bootstrap propagates p_n uncertainty into a_t uncertainty.
# Under the pooled approach, a single p_n_pooled is resampled each iteration.
# Under the stratified approach, each stratum gets its own resampled p_n_b.
# The CI comparison reveals two things:
#   1. The pooled CI is wider (extra between-stratum variance masquerades as noise).
#   2. The pooled CI is shifted upward (the bias from 6b persists).

d_a <- d_a_per_day

# Pooled bootstrap
p_pool     <- sum(season$d_n) / (nrow(season) * d_a)
pool_ests  <- numeric(B)
for (b in seq_len(B)) {
  d_n_b <- rbinom(nrow(season), size = d_a, prob = p_pool)
  p_n_b <- sum(d_n_b) / (nrow(season) * d_a)
  pool_ests[b] <- sum(season$a_d) / (1 - p_n_b)
}
ci_pool <- quantile(pool_ests, c(0.025, 0.975))

# Stratified bootstrap
p_strata <- season |>
  group_by(Collaps) |>
  summarise(p_n_hat = sum(d_n) / (n() * d_a), .groups = "drop")
strat_ests <- numeric(B)
for (b in seq_len(B)) {
  smpl <- season |>
    left_join(p_strata, by = "Collaps", suffix = c("", ".s")) |>
    mutate(d_n_b = rbinom(n(), size = d_a, prob = p_n_hat.s))
  strat_b <- smpl |>
    group_by(Collaps) |>
    summarise(p_n_b   = sum(d_n_b) / (n() * d_a),
              a_d_sum = sum(a_d), .groups = "drop")
  strat_ests[b] <- sum(strat_b$a_d_sum / (1 - strat_b$p_n_b))
}
ci_strat <- quantile(strat_ests, c(0.025, 0.975))

cat("Problem 6e (B =", B, ")\n")
cat("  Truth      =", round(truth_total), "\n")
cat("  Pooled CI  = [", round(ci_pool[1]),  ",", round(ci_pool[2]),
    "]  width =", round(diff(ci_pool)), "\n")
cat("  Stratified CI = [", round(ci_strat[1]), ",", round(ci_strat[2]),
    "]  width =", round(diff(ci_strat)), "\n")

# ---- Extension: how much bias grows with the p_n contrast ----
# The pooling bias is roughly proportional to the variance of p_n across strata.
# Here we sweep the contrast while keeping the mean p_n fixed and measure bias.

mean_pn   <- (p_n_early + p_n_late) / 2
contrasts <- c(0.00, 0.10, 0.20, 0.30, 0.40)   # p_n_late - p_n_early
cat("\nExtension: pooling bias vs p_n contrast (mean p_n fixed at", mean_pn, ")\n")
for (delta in contrasts) {
  pn_e <- mean_pn - delta / 2
  pn_l <- mean_pn + delta / 2
  if (pn_e < 0 || pn_l > 1) next
  truth_c  <- ndays / 2 * (a_d_per_day / (1 - pn_e)) +
              ndays / 2 * (a_d_per_day / (1 - pn_l))
  pn_pool  <- mean_pn   # pooled estimate is just the mean when groups are equal size
  pooled_c <- ndays * a_d_per_day / (1 - pn_pool)
  cat("  delta =", delta, "  truth =", round(truth_c),
      "  pooled =", round(pooled_c),
      "  bias =", round(pooled_c - truth_c), "\n")
}
# Bias grows nonlinearly with delta because 1/(1-p_n) is convex.
# Even a modest contrast of 0.20 (typical in LGD data) produces meaningful
# upward bias — justifying SCRAPI's mandatory Collaps stratification.

# ---- Forward pointer ----
# Section 7 replaces the two-outcome (Wild/HNC) model with a full multinomial,
# generalising the stock-composition estimator to K >= 2 stocks and building
# the joint-likelihood skeleton of the Delomas & Hess PBT estimator.
