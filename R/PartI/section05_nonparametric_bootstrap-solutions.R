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


# Narrative walkthrough ----
# The wrapper functions above ARE the correct answers.
# What follows is plain runnable code that explains the nonparametric bootstrap
# layer in SCRAPI's bootsmolt, traces data through thetahat(), and triggers the
# rare-stock dimension-mismatch error so you can see exactly what breaks.

set.seed(2026)

# ---- shared parameters (match stub argument names and test values) ----
stocks <- c("LOSALM", "CHMBLN", "IMNAHA")   # stock identifiers (PGrp values)
strats <- c("S1", "S2", "S3")               # temporal strata (Strat/Stratum)
nfish  <- 90L                               # fish sampled at the trap
B      <- 1000L                             # bootstrap iterations

SR <- 5/6 * 0.45   # combined smolt-trap detection probability (constant here)

# ---- Problem 5a: build toy SCRAPI-shaped data frames ----
# SCRAPI's bootsmolt() receives three data frames.  We build toy versions
# so every subsequent problem has realistic inputs.

AllPrime <- tibble(
  Strat = sample(strats, size = nfish, replace = TRUE),
  PGrp  = sample(stocks, size = nfish, replace = TRUE,
                 prob = c(0.55, 0.35, 0.10)),
  SR    = SR
)
# AllPrime: one row per sampled fish.  Strat = temporal stratum, PGrp = stock
# assignment from PBT/GSI, SR = trap detection probability for that fish.

RearData <- AllPrime |>
  mutate(Rear    = sample(c("W", "HNC"), n(), replace = TRUE,
                          prob = c(0.7, 0.3)),
         Stratum = Strat,
         True    = SR) |>
  select(Rear, Stratum, True)
# RearData: rearing type (Wild vs Hatchery Non-Clipped) per fish.
# True = the same detection probability, used for inverse weighting in thetahat.

ndays <- 30L
passdata <- tibble(
  Stratum = rep(strats, each = ndays / length(strats)),
  Tally   = rbinom(ndays, size = 200L, prob = SR),
  Ptrue   = SR
)
# passdata: daily passage counts from the trap counter (SCRAPI passdata layout).

cat("Problem 5a: data frames built\n")
cat("  AllPrime rows =", nrow(AllPrime), "  passdata rows =", nrow(passdata), "\n")

# ---- Problem 5b: FishWH weighted resample (SCRAPI lines 148-157) ----
# The nonparametric bootstrap resamples fish within each stratum to propagate
# rearing-type composition uncertainty.  Fish are resampled WITH REPLACEMENT
# using their True (inverse detection probability) as sampling weights — fish
# that are hard to detect are upweighted so they have a fair chance of being
# drawn in the resample.

FishWH <- RearData   # use RearData as the FishWH input

H <- 0
WHstar <- NULL
for (h in strats) {
  justwk <- FishWH[FishWH$Stratum == h, ]   # subset to this stratum
  nwk    <- nrow(justwk)
  if (nwk == 0L) next                        # skip empty strata safely
  i      <- sample.int(nwk, replace = TRUE, prob = unlist(justwk$True))
  # sample.int draws row indices with replacement.  prob = True weights means
  # a fish with True = 2 is twice as likely to be selected as one with True = 1.
  wkstar <- justwk[i, ]   # bootstrap sample for this stratum
  if (H == 0) { WHstar <- wkstar; H <- 1 } else { WHstar <- rbind(WHstar, wkstar) }
}
cat("Problem 5b: resampled FishWH rows =", nrow(WHstar), "\n")

# ---- Problem 5c: FishDat weighted resample (SCRAPI lines 159-169) ----
# Same logic as 5b but applied to AllPrime (individual genetic/PBT records).
# Column names differ: Strat instead of Stratum, SR instead of True.

FishDat <- AllPrime

H <- 0
DatStar <- NULL
for (h in strats) {
  justwk <- FishDat[FishDat$Strat == h, ]
  nwk    <- nrow(justwk)
  if (nwk == 0L) next
  i      <- sample.int(nwk, replace = TRUE, prob = unlist(justwk$SR))
  wkstar <- justwk[i, ]
  if (H == 0) { DatStar <- wkstar; H <- 1 } else { DatStar <- rbind(DatStar, wkstar) }
}
cat("Problem 5c: resampled FishDat rows =", nrow(DatStar), "\n")

# ---- Problem 5d: thetahat_toy — the core estimator ----
# thetahat() translates raw passage counts and individual fish assignments
# into a total wild smolt estimate and per-stock wild estimates.
#
# Step 1: expand daily counts to passage scale.
#   dailypass = Tally / Ptrue   (same expansion as Section 4)
#   bystrata  = sum within each stratum → stratum-level wild fish available
#
# Step 2: partition each stratum into Wild vs HNC using inverse-True weights.
#   HNCWstrat  = tapply(1/True, list(Stratum, Rear), sum)
#   HNCWprop   = row-normalise → P(Wild | stratum)
#   WildStrata = P(Wild | stratum) * bystrata
#
# Step 3: allocate wild fish to stocks using inverse-SR weights (same as Section 2d).
#   Primarystrata      = tapply(1/SR, list(Strat, PGrp), sum)
#   Primaryproportions = row-normalise → P(stock | stratum)
#   Primaryests        = Primaryproportions' %*% WildStrata

dailypass  <- passdata$Tally / passdata$Ptrue
bystrata   <- tapply(dailypass, passdata$Stratum, sum)

HNCWstrat  <- tapply(1 / RearData$True,
                     list(RearData$Stratum, RearData$Rear), sum)
HNCWstrat[is.na(HNCWstrat)] <- 0
HNCWprop   <- prop.table(HNCWstrat, margin = 1)
PWild      <- HNCWprop[, "W"]
WildStrata <- PWild * bystrata     # wild smolts per stratum
TotalWild  <- sum(WildStrata)

Primarystrata <- tapply(1 / AllPrime$SR,
                        list(factor(AllPrime$Strat, levels = strats),
                             factor(AllPrime$PGrp,  levels = stocks)),
                        sum)
Primarystrata[is.na(Primarystrata)] <- 0
Primaryproportions <- prop.table(Primarystrata, margin = 1)
Primaryproportions[is.na(Primaryproportions)] <- 0
Primaryests <- as.vector(t(Primaryproportions) %*% as.vector(WildStrata))
names(Primaryests) <- stocks

cat("Problem 5d: TotalWild =", round(TotalWild),
    "  Primaryests:", round(Primaryests), "\n")
# TotalWild ≈ ndays * 200 * SR * P(Wild) — sanity check.

# ---- Problem 5e: full bootsmolt loop ----
# Iteration b=1 uses the real data (no resampling) — this gives the point
# estimate.  Iterations b>1 resample counts (4a pattern) AND resample fish
# (5b/5c) to capture ALL sources of uncertainty simultaneously.

p       <- 1L + length(stocks)   # one TotalWild + one per stock
theta.b <- matrix(numeric(p * B), ncol = p)
for (b in seq_len(B)) {
  if (b == 1L) {
    dailyStar <- passdata;  RearStar <- RearData;  indivStar <- AllPrime
  } else {
    cntstar <- vapply(seq_len(nrow(passdata)), function(i) {
      dp <- round(passdata$Tally[i] / passdata$Ptrue[i])
      if (dp == 0L) 0L else rbinom(1L, dp, passdata$Ptrue[i])
    }, integer(1))
    dailyStar  <- data.frame(Stratum = passdata$Stratum,
                              Tally   = cntstar,
                              Ptrue   = passdata$Ptrue)
    # Resample rearing-type fish (propagates Wild/HNC uncertainty).
    RearStar   <- section5_problem_5b_fish(RearData, strats)
    # Resample stock-assignment fish (propagates stock-proportion uncertainty).
    indivStar  <- section5_problem_5c_fish(AllPrime, strats)
  }
  eststar     <- section5_problem_5d_fish(dailyStar, RearStar, indivStar,
                                          Pgrps = stocks, strats = strats)
  theta.b[b, ] <- c(eststar$TotalWild, eststar$Primaryests)
}
CI <- t(apply(theta.b, 2, quantile, c(0.025, 0.975)))
rownames(CI) <- c("WildSmolts", stocks)
colnames(CI) <- c("LCI", "UCI")
cat("Problem 5e: 95% bootstrap CIs\n")
print(round(CI))

# ---- Problem 5f: deliberately trigger Error 1 ----
# SCRAPI pre-allocates theta.b with ncol = length(Pgrps).  If a rare stock
# disappears from a bootstrap resample, thetahat() returns a shorter vector
# and the assignment theta.b[b,] <- c(...) recycles or errors.

rare_stocks   <- c(stocks, "LOCLWR")
AllPrime_rare <- bind_rows(
  AllPrime,
  tibble(Strat = c("S1", "S1"), PGrp = c("LOCLWR", "LOCLWR"), SR = SR)
)
# Only 2 LOCLWR fish, both in S1.  In most bootstrap resamples S1 may not
# draw either of them, making LOCLWR vanish from the resample entirely.

caught <- tryCatch({
  for (b in seq_len(50L)) {
    indivStar <- section5_problem_5c_fish(AllPrime_rare, strats)
    eststar   <- section5_problem_5d_fish(passdata, RearData, indivStar,
                                          Pgrps = rare_stocks, strats = strats)
    theta.b_rare <- numeric(5L)
    theta.b_rare[] <- c(eststar$TotalWild, eststar$Primaryests)
  }
  "no error (lucky path)"
}, error   = function(e) conditionMessage(e),
   warning = function(w) conditionMessage(w))
cat("Problem 5f caught:", caught, "\n")

# ---- Extension: what Error 1 looks like with a fixed ncol ----
# In production SCRAPI, theta.b has a fixed number of columns set before the
# loop.  When Primaryests drops a stock, the vector is shorter than expected
# and R either recycles silently (wrong numbers) or throws a subscript error.
# We reproduce both outcomes.

p_rare     <- 1L + length(rare_stocks)   # expect 5 columns
theta_bad  <- matrix(NA_real_, nrow = 10L, ncol = p_rare)

silent_recycle <- tryCatch({
  for (b in seq_len(10L)) {
    indivStar <- section5_problem_5c_fish(AllPrime_rare, strats)
    eststar   <- section5_problem_5d_fish(passdata, RearData, indivStar,
                                          Pgrps = rare_stocks, strats = strats)
    short_vec <- c(eststar$TotalWild, eststar$Primaryests)
    theta_bad[b, ] <- short_vec   # may recycle if length(short_vec) < p_rare
  }
  "recycled silently — theta_bad contains WRONG values"
}, warning = function(w) paste("warning:", conditionMessage(w)),
   error   = function(e) paste("error:",   conditionMessage(e)))
cat("Extension recycling outcome:", silent_recycle, "\n")
# The fix: always run section5_problem_5f_fish() to catch rare stocks BEFORE
# calling SCRAPI(), or consolidate rare stocks into an "Other" category.

# ---- Forward pointer ----
# Section 6 introduces Collaps (temporal strata for p_n) and shows how
# pooling across strata with different nighttime rates biases a_t upward,
# then reproduces SCRAPI's Collaps assignment loop and its Error 3.
