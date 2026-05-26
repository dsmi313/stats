#---------------------------------------------------------
# File:   section04_parametric_bootstrap.R
# Part I, Section 4 - Parametric bootstrap: bootsmolt's daily-count layer
#
# Repo pointer (SCOBI/R/SCRAPI.r):
#   lines 139-145: cntstar[i] <- rbinom(1, dailypass[i], LGDdaily$Ptrue[i])
# Repo pointer (escapeLGD/R/night_fall_reascend_wc_binom.R):
#   line 150: wc_binom[[2]][,i] <- rbinom(boots, wc[i], wc_prop) / wc_prop
#
# This is SCRAPI smolt-trap territory: Ptrue = SampleRate * GuidanceEfficiency
# is the combined detection probability for the smolt bypass system.
#---------------------------------------------------------
# Section 4 ----

#--------------------------------------
# Problem 4a: Implement bootsmolt's daily-count step (SCRAPI verbatim)
section4_problem_4a_fish <- function(LGDdaily, B) {
  cat("\n----------------------------------\n")
  cat("Problem 4a: bootsmolt daily-count step (SCRAPI lines 139-145)\n")

  # Arguments:
  #   LGDdaily = data.frame(Stratum, Tally, Ptrue) -- smolt-trap layout
  #   B        = number of bootstrap iterations
  #
  # Mirror SCRAPI lines 139-145 EXACTLY:
  #   dailypass <- round(LGDdaily$Tally / LGDdaily$Ptrue)
  #   cntstar = numeric(ndays)
  #   for (i in 1:ndays)
  #     if (dailypass[i] != 0) cntstar[i] <- rbinom(1, dailypass[i], LGDdaily$Ptrue[i])
  # Then build dailyStar <- data.frame(Stratum, Tally = cntstar, Ptrue).
  # Repeat B times and return a B x ndays matrix of expanded estimates
  # (cntstar / Ptrue).

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 4b: Vectorized escapeLGD-style daily-count bootstrap
section4_problem_4b_fish <- function(LGDdaily, B) {
  cat("\n----------------------------------\n")
  cat("Problem 4b: Vectorized escapeLGD daily-count bootstrap\n")

  # Reproduce escapeLGD line 150 for each day i:
  #   rbinom(B, dailypass[i], LGDdaily$Ptrue[i]) / LGDdaily$Ptrue[i]
  # Return a B x ndays matrix.

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 4c: Coverage check across simulated seasons
section4_problem_4c_fish <- function(ndays, true_pass_per_day, Ptrue,
                                     B, nseasons) {
  cat("\n----------------------------------\n")
  cat("Problem 4c: 95% bootstrap-CI coverage across simulated seasons\n")

  # Simulate nseasons seasons. For each: generate Tally, build LGDdaily,
  # run section4_problem_4b_fish, take quantile(rowSums, c(.025, .975)),
  # check whether the true season total (ndays * true_pass_per_day) is in.

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 4d: Plot the bootstrap distribution of the season total
section4_problem_4d_fish <- function(LGDdaily, B) {
  cat("\n----------------------------------\n")
  cat("Problem 4d: Plot bootstrap distribution of season total\n")

  # Run section4_problem_4b_fish, sum each row, plot a histogram of the
  # season totals with the 95% CI marked.

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}
