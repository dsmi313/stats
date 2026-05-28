#---------------------------------------------------------
# File:   section05_nonparametric_bootstrap.R
# Part I, Section 5 - bootsmolt's weighted stratum resample + thetahat
#
# Repo pointers (SCOBI/R/SCRAPI.r):
#   line 94:   Primarystrata <- mApply(1/Fish$SR, list(Strat, PGrp), sum)
#   line 131:  theta.b <- matrix(numeric(p*B), ncol = p)
#   lines 148-157: weighted resample of FishWH by stratum (True weights)
#   lines 159-169: weighted resample of FishDat by stratum (SR weights)
#   lines 74-126:  thetahat function (signature you will mirror)
#   line 175:      theta.b[b, ] <- c(eststar[[1]], t(eststar[[4]]), ...)
#
# Inverse-SR weighting (problem 5a) is the Horvitz-Thompson estimator at
# the heart of thetahat() -- build it first, then assemble the full estimator.
#
# Place every answer inside the wrapper functions below.
# Run section05_nonparametric_bootstrap-test.R after sourcing this file.
#---------------------------------------------------------
# Section 5 ----

#--------------------------------------
# Problem 5a: Inverse-SR weighting -- the Horvitz-Thompson step inside thetahat()
section5_problem_5a_fish <- function(stocks, strats, SR, n_fish) {
  cat("\n----------------------------------\n")
  cat("Problem 5a: Inverse-SR weighting (SCRAPI thetahat pattern, SMOLT trap)\n")

  # Arguments:
  #   stocks = character vector of stock names (PGrp values)
  #   strats = character vector of stratum names (Strat values)
  #   SR     = combined smolt-trap detection probability (SampleRate * GuidanceEfficiency)
  #   n_fish = number of sampled fish
  #
  # Build an AllPrime tibble with columns Strat, PGrp, SR (one row per
  # fish). Use tapply(1/SR, list(Strat, PGrp), sum) to reproduce SCRAPI
  # line 94, then row-normalize to get within-stratum proportions.
  # This inverse-weighting step recurs inside thetahat() (problem 5e).

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 5b: Build toy RearData, AllPrime, passdata matching SCRAPI columns
section5_problem_5b_fish <- function(stocks, strats, nfish) {
  cat("\n----------------------------------\n")
  cat("Problem 5b: Build toy SCRAPI-shaped data frames\n")

  # Arguments:
  #   stocks = character vector of stock names (PGrp values)
  #   strats = character vector of stratum names (Strat/Stratum values)
  #   nfish  = number of fish to simulate
  #
  # Build:
  #   AllPrime = tibble(Strat, PGrp, SR)
  #   RearData = tibble(Rear, Stratum, True)
  #   passdata = tibble(Stratum, Tally, Ptrue)
  # Return all three in a list.

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 5c: Verbatim port of SCRAPI's FishWH weighted resample
section5_problem_5c_fish <- function(FishWH, strats) {
  cat("\n----------------------------------\n")
  cat("Problem 5c: Weighted FishWH resample (SCRAPI lines 148-157)\n")

  # Mirror SCRAPI.r:148-157 verbatim:
  #   H <- 0
  #   for (h in strats) {
  #     justwk <- FishWH[FishWH$Stratum == h, ]
  #     nwk    <- nrow(justwk)
  #     i      <- sample.int(nwk, replace = TRUE, prob = unlist(justwk$True))
  #     wkstar <- justwk[i, ]
  #     if (H == 0) { WHstar <- wkstar; H <- 1 }
  #     else        { WHstar <- rbind(WHstar, wkstar) }
  #   }

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 5d: Verbatim port of SCRAPI's FishDat weighted resample
section5_problem_5d_fish <- function(FishDat, strats) {
  cat("\n----------------------------------\n")
  cat("Problem 5d: Weighted FishDat resample (SCRAPI lines 159-169)\n")

  # Same loop as 5c but Strat column and SR weights (not Stratum/True).

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 5e: Toy thetahat() matching SCRAPI's signature
section5_problem_5e_fish <- function(passage, RearDat, Fish, Pgrps, strats) {
  cat("\n----------------------------------\n")
  cat("Problem 5e: thetahat_toy() (SCRAPI lines 74-126)\n")

  # Inputs (matching SCRAPI lines 74-87):
  #   passage = data.frame(Stratum, Tally, Ptrue)
  #   RearDat = data.frame(Rear, Stratum, True)
  #   Fish    = data.frame(Strat, PGrp, SR)
  # Output: list(TotalWild, WildStrata, Primaryproportions, Primaryests)
  # Use:
  #   dailypass <- passage$Tally / passage$Ptrue
  #   bystrata  <- tapply(dailypass, passage$Stratum, sum)
  #   HNCWstrat <- tapply(1/RearDat$True, list(RearDat$Stratum, RearDat$Rear), sum)
  #   Primarystrata <- tapply(1/Fish$SR, list(Fish$Strat, Fish$PGrp), sum)
  # The last line is problem 5a applied inside thetahat().

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 5f: Full bootsmolt loop on toy data
section5_problem_5f_fish <- function(passdata, RearData, AllPrime,
                                     stocks, strats, B) {
  cat("\n----------------------------------\n")
  cat("Problem 5f: Full bootsmolt loop (SCRAPI lines 128-189)\n")

  # Build theta.b matrix with p = 1 + length(stocks) columns. For b = 1 use
  # the real data; for b > 1 use 4a's daily-count bootstrap + 5c/5d
  # weighted resamples + 5e's thetahat_toy(). Stack TotalWild and
  # Primaryests into theta.b[b, ]. Return CIs from quantile().

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 5g: Deliberately trigger Error 1 (theta.b dimension mismatch)
section5_problem_5g_fish <- function() {
  cat("\n----------------------------------\n")
  cat("Problem 5g: Trigger SCRAPI Error 1 with a rare stock\n")

  # Add a 4th rare stock (e.g. LOCLWR) with only 2 fish in 1 stratum, then
  # run bootstrap iterations where the rare stock can disappear entirely
  # from the resample. Show recycling error / warning.

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}
