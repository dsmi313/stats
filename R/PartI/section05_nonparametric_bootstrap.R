#---------------------------------------------------------
# File:   section05_nonparametric_bootstrap.R
# Part I, Section 5 - nonparametric bootstrap diagnostics
# Workshop flow: simulate -> inspect -> estimate -> diagnose -> explain -> EASE/SCRAPI
#---------------------------------------------------------

#---------------------------------------------------------
# File:   section05_nonparametric_bootstrap-solutions.R
# Part I, Section 5 solutions
#---------------------------------------------------------
# Section 5 solutions ----

library(dplyr)
library(tibble)

#--------------------------------------
# Problem 5a: Build toy SCRAPI-shaped data frames
section5_problem_5a_fish <- function(stocks, strats, nfish) {
  cat("\n----------------------------------\n")
  cat("Problem 5a: Build toy SCRAPI-shaped data frames\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  SR <- 5/6 * 0.45    # smolt-trap combined detection probability
  AllPrime <- tibble(
    Strat = sample(strats, size = nfish, replace = TRUE),
    PGrp  = sample(stocks, size = nfish, replace = TRUE,
                   prob = c(0.55, 0.35, 0.10)),
    SR    = SR
  )
  RearData <- AllPrime |>
    mutate(Rear    = sample(c("W", "HNC"), n(), replace = TRUE,
                            prob = c(0.7, 0.3)),
           Stratum = Strat,
           True    = SR) |>
    select(Rear, Stratum, True)
  ndays <- 30L
  passdata <- tibble(
    Stratum = rep(strats, each = ndays / length(strats)),
    Tally   = rbinom(ndays, size = 200L, prob = SR),
    Ptrue   = SR
  )
  cat("nfish =", nfish, "  ndays =", ndays, "\n")
  cat("AllPrime head:\n"); print(head(AllPrime, 3))
  invisible(list(AllPrime = AllPrime,
                 RearData = RearData,
                 passdata = passdata))
}


# Problem 5b: Verbatim port of SCRAPI's FishWH weighted resample
section5_problem_5b_fish <- function(FishWH, strats) {
  cat("\n----------------------------------\n")
  cat("Problem 5b: Weighted FishWH resample (SCRAPI lines 148-157)\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  H <- 0
  WHstar <- NULL
  for (h in strats) {
    justwk <- FishWH[FishWH$Stratum == h, ]
    nwk    <- nrow(justwk)
    if (nwk == 0L) next
    i      <- sample.int(nwk, replace = TRUE, prob = unlist(justwk$True))
    wkstar <- justwk[i, ]
    if (H == 0) { WHstar <- wkstar; H <- 1
    } else      { WHstar <- rbind(WHstar, wkstar) }
  }
  cat("Resampled FishWH rows =", nrow(WHstar), "\n")
  invisible(WHstar)
}


# Problem 5c: Verbatim port of SCRAPI's FishDat weighted resample
section5_problem_5c_fish <- function(FishDat, strats) {
  cat("\n----------------------------------\n")
  cat("Problem 5c: Weighted FishDat resample (SCRAPI lines 159-169)\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  H <- 0
  DatStar <- NULL
  for (h in strats) {
    justwk <- FishDat[FishDat$Strat == h, ]
    nwk    <- nrow(justwk)
    if (nwk == 0L) next
    i      <- sample.int(nwk, replace = TRUE, prob = unlist(justwk$SR))
    wkstar <- justwk[i, ]
    if (H == 0) { DatStar <- wkstar; H <- 1
    } else      { DatStar <- rbind(DatStar, wkstar) }
  }
  cat("Resampled FishDat rows =", nrow(DatStar), "\n")
  invisible(DatStar)
}


# Problem 5d: Toy thetahat() matching SCRAPI's signature
section5_problem_5d_fish <- function(passage, RearDat, Fish, Pgrps, strats) {
  cat("\n----------------------------------\n")
  cat("Problem 5d: thetahat_toy() (SCRAPI lines 74-126)\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  dailypass  <- passage$Tally / passage$Ptrue
  bystrata   <- tapply(dailypass, passage$Stratum, sum)

  HNCWstrat  <- tapply(1 / RearDat$True,
                       list(RearDat$Stratum, RearDat$Rear), sum)
  HNCWstrat[is.na(HNCWstrat)] <- 0
  HNCWprop   <- prop.table(HNCWstrat, margin = 1)
  PWild      <- HNCWprop[, "W"]
  WildStrata <- PWild * bystrata
  TotalWild  <- sum(WildStrata)

  Primarystrata <- tapply(1 / Fish$SR,
                          list(factor(Fish$Strat, levels = strats),
                               factor(Fish$PGrp,  levels = Pgrps)),
                          sum)
  Primarystrata[is.na(Primarystrata)] <- 0
  Primaryproportions <- prop.table(Primarystrata, margin = 1)
  Primaryproportions[is.na(Primaryproportions)] <- 0
  Primaryests <- as.vector(t(Primaryproportions) %*% as.vector(WildStrata))
  names(Primaryests) <- Pgrps

  cat("TotalWild =", round(TotalWild), "\n")
  cat("Primaryests:\n"); print(round(Primaryests))
  invisible(list(TotalWild = TotalWild,
                 WildStrata = WildStrata,
                 Primaryproportions = Primaryproportions,
                 Primaryests = Primaryests))
}


# Problem 5e: Full bootsmolt loop on toy data
section5_problem_5e_fish <- function(passdata, RearData, AllPrime,
                                     stocks, strats, B) {
  cat("\n----------------------------------\n")
  cat("Problem 5e: Full bootsmolt loop (SCRAPI lines 128-189)\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  p <- 1L + length(stocks)
  theta.b <- matrix(numeric(p * B), ncol = p)
  for (b in seq_len(B)) {
    if (b == 1) {
      dailyStar <- passdata
      RearStar  <- RearData
      indivStar <- AllPrime
    } else {
      cntstar <- vapply(seq_len(nrow(passdata)), function(i) {
        if (round(passdata$Tally[i] / passdata$Ptrue[i]) == 0) 0L
        else rbinom(1, round(passdata$Tally[i] / passdata$Ptrue[i]),
                    passdata$Ptrue[i])
      }, integer(1))
      dailyStar <- data.frame(Stratum = passdata$Stratum,
                              Tally   = cntstar,
                              Ptrue   = passdata$Ptrue)
      RearStar  <- section5_problem_5b_fish(RearData, strats)
      indivStar <- section5_problem_5c_fish(AllPrime, strats)
    }
    eststar <- section5_problem_5d_fish(dailyStar, RearStar, indivStar,
                                        Pgrps = stocks, strats = strats)
    theta.b[b, ] <- c(eststar$TotalWild, eststar$Primaryests)
  }
  CI <- t(apply(theta.b, 2, quantile, c(0.025, 0.975)))
  rownames(CI) <- c("WildSmolts", stocks)
  colnames(CI) <- c("LCI", "UCI")
  cat("B =", B, "  p =", p, "\n")
  cat("95% bootstrap CIs:\n"); print(round(CI))
  invisible(list(theta.b = theta.b, CI = CI))
}


# Problem 5f: Deliberately trigger Error 1
section5_problem_5f_fish <- function() {
  cat("\n----------------------------------\n")
  cat("Problem 5f: Trigger SCRAPI Error 1 with a rare stock\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  stocks <- c("LOSALM", "CHMBLN", "IMNAHA")
  strats <- c("S1", "S2", "S3")
  data <- section5_problem_5a_fish(stocks = stocks,
                                   strats = strats, nfish = 90L)
  # Add a rare stock confined to one stratum
  rare_stocks <- c(stocks, "LOCLWR")
  AllPrime_rare <- bind_rows(
    data$AllPrime,
    tibble(Strat = c("S1", "S1"), PGrp = c("LOCLWR", "LOCLWR"),
           SR = 5/6 * 0.45)
  )
  caught <- tryCatch({
    for (b in seq_len(50L)) {
      indivStar <- section5_problem_5c_fish(AllPrime_rare, strats)
      eststar <- section5_problem_5d_fish(data$passdata, data$RearData,
                                          indivStar,
                                          Pgrps = rare_stocks,
                                          strats = strats)
      theta.b_rare <- numeric(5L)
      theta.b_rare[] <- c(eststar$TotalWild, eststar$Primaryests)
    }
    "no error fired (lucky bootstrap path)"
  }, error = function(e) conditionMessage(e),
     warning = function(w) conditionMessage(w))
  cat("Caught:", caught, "\n")
  cat("This is exactly SCRAPI Error 1 (theta.b dimension mismatch).\n")
  cat("Fix: drop LOCLWR before SCRAPI() or move it to its own analysis.\n")
  invisible(caught)
}


# Section summary:
# 1. What was simulated? Fish-passage observations under this section's data-generating process.
# 2. What model was assumed? The estimator-specific model encoded in the wrapper functions.
# 3. What estimator was used? See section*_problem_* functions (MLE/bootstrap/stratified/composition).
# 4. What assumption was broken? This section includes a diagnostic failure-mode comparison.
# 5. What did the diagnostic show? Check plots and printed bias/CI summaries against truth.
# 6. Why does this matter for EASE/SCRAPI? These assumptions map directly to production expansion and uncertainty code paths.
