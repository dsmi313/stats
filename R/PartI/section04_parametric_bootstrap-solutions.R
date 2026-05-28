#---------------------------------------------------------
# File:   section04_parametric_bootstrap-solutions.R
# Part I, Section 4 solutions
#---------------------------------------------------------
# Section 4 solutions ----

library(ggplot2)
library(tibble)


#--------------------------------------
# Problem 4a: bootsmolt's daily-count step (SCRAPI verbatim)
section4_problem_4a_fish <- function(LGDdaily, B) {
  cat("\n----------------------------------\n")
  cat("Problem 4a: bootsmolt daily-count step (SCRAPI lines 139-145)\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  dailypass <- round(LGDdaily$Tally / LGDdaily$Ptrue)
  ndays     <- nrow(LGDdaily)
  out       <- matrix(NA_real_, nrow = B, ncol = ndays)
  for (b in seq_len(B)) {
    cntstar <- numeric(ndays)
    for (i in seq_len(ndays)) {
      if (dailypass[i] != 0) {
        cntstar[i] <- rbinom(1, dailypass[i], LGDdaily$Ptrue[i])
      }
    }
    dailyStar <- data.frame(Stratum = LGDdaily$Stratum,
                            Tally   = cntstar,
                            Ptrue   = LGDdaily$Ptrue)
    out[b, ] <- dailyStar$Tally / dailyStar$Ptrue
  }
  cat("B =", B, "  ndays =", ndays, "  point total =",
      round(sum(LGDdaily$Tally / LGDdaily$Ptrue)), "\n")
  invisible(out)
}


# Problem 4b: Vectorized escapeLGD-style daily-count bootstrap
section4_problem_4b_fish <- function(LGDdaily, B) {
  cat("\n----------------------------------\n")
  cat("Problem 4b: Vectorized escapeLGD daily-count bootstrap\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  dailypass <- round(LGDdaily$Tally / LGDdaily$Ptrue)
  out <- sapply(seq_len(nrow(LGDdaily)), function(i) {
    rbinom(B, dailypass[i], LGDdaily$Ptrue[i]) / LGDdaily$Ptrue[i]
  })
  cat("B =", B, "  ndays =", nrow(LGDdaily), "\n")
  invisible(out)
}


# Problem 4c: Coverage check across simulated seasons
section4_problem_4c_fish <- function(ndays, true_pass_per_day, Ptrue,
                                     B, nseasons) {
  cat("\n----------------------------------\n")
  cat("Problem 4c: 95% bootstrap-CI coverage across simulated seasons\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  hits <- vapply(seq_len(nseasons), function(s) {
    set.seed(s)
    Tally <- rbinom(ndays, true_pass_per_day, Ptrue)
    LD <- data.frame(Stratum = 1L, Tally = Tally, Ptrue = Ptrue)
    draws  <- section4_problem_4b_fish(LD, B = B)
    totals <- rowSums(draws)
    ci     <- quantile(totals, c(0.025, 0.975))
    truth  <- ndays * true_pass_per_day
    (truth >= ci[1]) && (truth <= ci[2])
  }, logical(1))
  cat("ndays =", ndays, "  true_pass_per_day =", true_pass_per_day,
      "  Ptrue =", Ptrue, "  B =", B,
      "  nseasons =", nseasons, "\n")
  cat("Coverage of 95% CI =", mean(hits), "\n")
  invisible(hits)
}


# Problem 4d: Plot the bootstrap distribution of the season total
section4_problem_4d_fish <- function(LGDdaily, B) {
  cat("\n----------------------------------\n")
  cat("Problem 4d: Plot bootstrap distribution of season total\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  draws  <- section4_problem_4b_fish(LGDdaily, B = B)
  totals <- rowSums(draws)
  ci     <- quantile(totals, c(0.025, 0.975))

  p <- ggplot(tibble(total = totals), aes(total)) +
    geom_histogram(bins = 60, fill = "steelblue", colour = "white") +
    geom_vline(xintercept = ci, colour = "firebrick",
               linewidth = 1, linetype = "dashed") +
    labs(title = "Section 4 - Bootstrap season totals (bootsmolt inner loop)",
         subtitle = "Dashed lines = 95% CI",
         x = "season total smolt passage", y = "Bootstrap draws")
  invisible(totals)
}


# Narrative walkthrough ----
# The wrapper functions above ARE the correct answers.
# What follows is plain runnable code that explains why parametric bootstrap
# works for smolt passage, compares the slow SCRAPI loop to the fast vectorised
# form, and shows where coverage collapses.

set.seed(2026)

# ---- shared parameters (match stub argument names and test values) ----
Ptrue          <- 5/6 * 0.45   # SampleRate * GuidanceEfficiency for smolt trap
ndays          <- 60L          # days in the toy season
true_pass_per_day <- 200L      # true daily smolt passage
B              <- 5000L        # bootstrap iterations
nseasons       <- 100L         # seasons for the coverage simulation

# Build a 60-day toy passage data frame matching the test driver layout.
LGDdaily <- data.frame(
  Stratum = rep(seq_len(6L), each = 10L),   # six weekly strata, 10 days each
  Tally   = rbinom(ndays, size = true_pass_per_day, prob = Ptrue),
  Ptrue   = Ptrue
)
# Tally   = fish actually counted in the trap (the raw sensor reading).
# Ptrue   = combined probability that a passing fish reaches the counter.
# The point estimate for daily passage is Tally / Ptrue (expansion by 1/p).

# ---- Problem 4a: SCRAPI bootsmolt daily-count loop (slow, verbatim) ----
# SCRAPI lines 139-145 bootstrap the Tally by resampling from its implied
# Binomial distribution.  The idea:
#   dailypass = round(Tally / Ptrue)   ← infer the true passage that would
#                                         have produced this Tally
#   cntstar   = rbinom(1, dailypass, Ptrue)   ← resample a new Tally
#   expanded  = cntstar / Ptrue               ← expand back to passage scale
# Repeat for each day, then sum across days for a bootstrap season total.

dailypass <- round(LGDdaily$Tally / LGDdaily$Ptrue)
# round() is SCRAPI's deliberate choice: dailypass must be an integer because
# it is the size argument to rbinom().

ndays_act <- nrow(LGDdaily)
out_a <- matrix(NA_real_, nrow = 500L, ncol = ndays_act)   # B=500 for speed demo
for (b in seq_len(500L)) {
  cntstar <- numeric(ndays_act)
  for (i in seq_len(ndays_act)) {
    # Guard: if the inferred passage is zero, cntstar stays zero (no fish to resample).
    if (dailypass[i] != 0)
      cntstar[i] <- rbinom(1, dailypass[i], LGDdaily$Ptrue[i])
  }
  dailyStar <- data.frame(Stratum = LGDdaily$Stratum,
                           Tally   = cntstar,
                           Ptrue   = LGDdaily$Ptrue)
  out_a[b, ] <- dailyStar$Tally / dailyStar$Ptrue   # expanded daily estimates
}
cat("Problem 4a (SCRAPI loop, B=500)\n")
cat("  Point total =", round(sum(LGDdaily$Tally / LGDdaily$Ptrue)), "\n")
cat("  Boot CI = [", round(quantile(rowSums(out_a), 0.025)),
    ",", round(quantile(rowSums(out_a), 0.975)), "]\n")

# ---- Problem 4b: vectorised escapeLGD equivalent ----
# The nested loop above is O(B × ndays).  sapply collapses the inner loop by
# drawing B samples for each day all at once using rbinom's vectorised form.

out_b <- sapply(seq_len(nrow(LGDdaily)), function(i) {
  rbinom(B, dailypass[i], LGDdaily$Ptrue[i]) / LGDdaily$Ptrue[i]
})
# out_b is B × ndays.  Each row is one bootstrap season's daily expanded estimates.
# rowSums gives B season totals.

season_totals <- rowSums(out_b)
ci_b <- quantile(season_totals, c(0.025, 0.975))
cat("Problem 4b (vectorised, B=", B, ")\n", sep = "")
cat("  Boot CI = [", round(ci_b[1]), ",", round(ci_b[2]), "]\n")
# Should agree with 4a; vectorised version is ~10× faster.

# ---- Problem 4c: coverage across simulated seasons ----
# A 95% CI should contain the truth 95% of the time.
# We simulate nseasons independent seasons, build a CI for each, and count
# how often the true season total falls inside.

hits <- vapply(seq_len(nseasons), function(s) {
  set.seed(s)
  Tally_s <- rbinom(ndays, true_pass_per_day, Ptrue)
  LD_s    <- data.frame(Stratum = 1L, Tally = Tally_s, Ptrue = Ptrue)
  pass_s  <- round(LD_s$Tally / LD_s$Ptrue)
  draws   <- sapply(seq_len(nrow(LD_s)), function(i)
    rbinom(1000L, pass_s[i], LD_s$Ptrue[i]) / LD_s$Ptrue[i])
  ci_s    <- quantile(rowSums(draws), c(0.025, 0.975))
  truth_s <- ndays * true_pass_per_day
  (truth_s >= ci_s[1]) && (truth_s <= ci_s[2])
}, logical(1))

cat("Problem 4c\n")
cat("  Coverage =", mean(hits),
    " (target 0.95, expect 0.92-0.97 with nseasons =", nseasons, ")\n")
# Coverage close to 0.95 validates the parametric bootstrap for this setting.
# It can drop below nominal when Ptrue is very small (sparse counts) — see extension.

# ---- Problem 4d: distribution of one season's bootstrap totals ----
# The histogram of rowSums(out_b) shows the sampling distribution of the
# season-total estimator.  The shape should be roughly normal by CLT since
# we are summing 60 independent daily estimates.

totals_4d <- rowSums(out_b)
ci_4d     <- quantile(totals_4d, c(0.025, 0.975))
cat("Problem 4d\n")
cat("  Bootstrap season-total CI = [", round(ci_4d[1]),
    ",", round(ci_4d[2]), "]\n")
cat("  Point estimate =", round(sum(LGDdaily$Tally / LGDdaily$Ptrue)), "\n")

# ---- Extension: coverage collapses when Ptrue is very small ----
# When Ptrue is low, Tally can be zero on some days even when fish did pass.
# rbinom(1, dailypass=0, p) always returns 0, so those days contribute
# nothing to the bootstrap — the CI systematically undercovers the truth.

coverage_by_p <- vapply(c(0.10, 0.20, 0.40, 0.75), function(p) {
  hs <- vapply(seq_len(50L), function(s) {
    set.seed(s)
    tly <- rbinom(20L, 100L, p)
    dp  <- round(tly / p)
    drs <- sapply(seq_along(dp), function(i)
      rbinom(500L, max(dp[i], 0L), p) / p)
    ci  <- quantile(rowSums(drs), c(0.025, 0.975))
    truth <- 20L * 100L
    (truth >= ci[1]) && (truth <= ci[2])
  }, logical(1))
  mean(hs)
}, numeric(1))

cat("\nExtension: CI coverage vs Ptrue\n")
for (k in seq_along(c(0.10, 0.20, 0.40, 0.75)))
  cat("  Ptrue =", c(0.10, 0.20, 0.40, 0.75)[k],
      "  coverage =", round(coverage_by_p[k], 2), "\n")
# Coverage degrades as Ptrue shrinks because zero-Tally days are
# informative but the bootstrap treats them as no-fish days.
# escapeLGD addresses this via the guidance-efficiency correction — a
# separate multiplier that adjusts for days when the trap was not fishing.

# ---- Forward pointer ----
# Section 5 layers the nonparametric (FishWH and FishDat) resamples on top
# of this daily-count bootstrap to propagate stock-composition uncertainty
# through the full bootsmolt estimator.
