#---------------------------------------------------------
# File:   section05_nonparametric_bootstrap-solutions.R
# Part I, Section 5 solutions
#---------------------------------------------------------
# Section 5 solutions ----

library(dplyr)
library(tibble)

#--------------------------------------
# Problem 5a: Inverse-SR weighting -- the Horvitz-Thompson step inside thetahat()
section5_problem_5a_fish <- function(stocks, strats, SR, n_fish) {
  cat("\n----------------------------------\n")
  cat("Problem 5a: Inverse-SR weighting (SCRAPI thetahat pattern, SMOLT trap)\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  AllPrime <- tibble(
    Strat = sample(strats, size = n_fish, replace = TRUE),
    PGrp  = sample(stocks, size = n_fish, replace = TRUE),
    SR    = SR
  )
  # SCRAPI line 94 reproduced:
  Primarystrata <- tapply(1 / AllPrime$SR,
                          list(factor(AllPrime$Strat, levels = strats),
                               factor(AllPrime$PGrp,  levels = stocks)),
                          sum)
  Primarystrata[is.na(Primarystrata)] <- 0
  Primaryproportions <- prop.table(Primarystrata, margin = 1)

  cat("AllPrime rows =", n_fish, "  SR =", round(SR, 3), "\n")
  cat("Primarystrata (inverse-SR weighted counts):\n")
  print(round(Primarystrata, 1))
  cat("Primaryproportions (row-normalized within stratum):\n")
  print(round(Primaryproportions, 3))
  invisible(list(AllPrime = AllPrime,
                 Primarystrata = Primarystrata,
                 Primaryproportions = Primaryproportions))
}


# Problem 5b: Build toy SCRAPI-shaped data frames
section5_problem_5b_fish <- function(stocks, strats, nfish) {
  cat("\n----------------------------------\n")
  cat("Problem 5b: Build toy SCRAPI-shaped data frames\n")

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


# Problem 5c: Verbatim port of SCRAPI's FishWH weighted resample
section5_problem_5c_fish <- function(FishWH, strats) {
  cat("\n----------------------------------\n")
  cat("Problem 5c: Weighted FishWH resample (SCRAPI lines 148-157)\n")

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


# Problem 5d: Verbatim port of SCRAPI's FishDat weighted resample
section5_problem_5d_fish <- function(FishDat, strats) {
  cat("\n----------------------------------\n")
  cat("Problem 5d: Weighted FishDat resample (SCRAPI lines 159-169)\n")

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


# Problem 5e: Toy thetahat() matching SCRAPI's signature
section5_problem_5e_fish <- function(passage, RearDat, Fish, Pgrps, strats) {
  cat("\n----------------------------------\n")
  cat("Problem 5e: thetahat_toy() (SCRAPI lines 74-126)\n")

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


# Problem 5f: Full bootsmolt loop on toy data
section5_problem_5f_fish <- function(passdata, RearData, AllPrime,
                                     stocks, strats, B) {
  cat("\n----------------------------------\n")
  cat("Problem 5f: Full bootsmolt loop (SCRAPI lines 128-189)\n")

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
      RearStar  <- section5_problem_5c_fish(RearData, strats)
      indivStar <- section5_problem_5d_fish(AllPrime, strats)
    }
    eststar <- section5_problem_5e_fish(dailyStar, RearStar, indivStar,
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


# Problem 5g: Deliberately trigger Error 1
section5_problem_5g_fish <- function() {
  cat("\n----------------------------------\n")
  cat("Problem 5g: Trigger SCRAPI Error 1 with a rare stock\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  stocks <- c("LOSALM", "CHMBLN", "IMNAHA")
  strats <- c("S1", "S2", "S3")
  data <- section5_problem_5b_fish(stocks = stocks,
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
      indivStar <- section5_problem_5d_fish(AllPrime_rare, strats)
      eststar <- section5_problem_5e_fish(data$passdata, data$RearData,
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

# ---- Problem 5a: inverse-SR weighting — the Horvitz-Thompson motivator ----
# Before assembling thetahat(), understand the weighting step at its core.
# At the smolt bypass, SR = SampleRate * GuidanceEfficiency is the probability
# that a passing fish reaches the counter.  A fish caught when SR = 0.20 represents
# 1/0.20 = 5 fish in the total run; one caught when SR = 0.40 represents 2.5 fish.
# Summing 1/SR within each Strat × PGrp cell gives the Horvitz-Thompson estimator
# for the number of fish in that cell — unbiased regardless of how SR varies.

n_fish <- nfish   # same value, problem 5a uses n_fish as argument name

AllPrime_5a <- tibble(
  Strat = sample(strats, size = n_fish, replace = TRUE),
  PGrp  = sample(stocks, size = n_fish, replace = TRUE),
  SR    = SR
)
# One row per sampled fish.  SCRAPI calls this AllPrime (the full individual-fish table).

# SCRAPI line 94: tapply sums 1/SR within every Strat × PGrp cell.
Primarystrata_5a <- tapply(1 / AllPrime_5a$SR,
                           list(factor(AllPrime_5a$Strat, levels = strats),
                                factor(AllPrime_5a$PGrp,  levels = stocks)),
                           sum)
Primarystrata_5a[is.na(Primarystrata_5a)] <- 0

# Row-normalize: convert weighted counts into within-stratum stock proportions.
Primaryproportions_5a <- prop.table(Primarystrata_5a, margin = 1)

cat("Problem 5a: Primaryproportions (each row sums to 1.0):\n")
print(round(Primaryproportions_5a, 3))
# When SR is constant, 1/SR cancels in the ratio and proportions equal raw counts.
# When SR varies by stratum (common in real data), the weighting corrects for
# unequal detectability.  This exact tapply(1/SR, ...) call appears verbatim
# inside thetahat() as problem 5e — now you know why it's there.

# ---- Problem 5b: build toy SCRAPI-shaped data frames ----
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

cat("Problem 5b: data frames built\n")
cat("  AllPrime rows =", nrow(AllPrime), "  passdata rows =", nrow(passdata), "\n")

# ---- Problem 5c: FishWH weighted resample (SCRAPI lines 148-157) ----
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
cat("Problem 5c: resampled FishWH rows =", nrow(WHstar), "\n")

# ---- Problem 5d: FishDat weighted resample (SCRAPI lines 159-169) ----
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
cat("Problem 5d: resampled FishDat rows =", nrow(DatStar), "\n")

# ---- Problem 5e: thetahat_toy — the core estimator ----
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
# Step 3: allocate wild fish to stocks using inverse-SR weights (same as problem 5a).
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

cat("Problem 5e: TotalWild =", round(TotalWild),
    "  Primaryests:", round(Primaryests), "\n")
# TotalWild ≈ ndays * 200 * SR * P(Wild) — sanity check.

# ---- Problem 5f: full bootsmolt loop ----
# Iteration b=1 uses the real data (no resampling) — this gives the point
# estimate.  Iterations b>1 resample counts (4a pattern) AND resample fish
# (5c/5d) to capture ALL sources of uncertainty simultaneously.

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
    RearStar   <- section5_problem_5c_fish(RearData, strats)
    # Resample stock-assignment fish (propagates stock-proportion uncertainty).
    indivStar  <- section5_problem_5d_fish(AllPrime, strats)
  }
  eststar     <- section5_problem_5e_fish(dailyStar, RearStar, indivStar,
                                          Pgrps = stocks, strats = strats)
  theta.b[b, ] <- c(eststar$TotalWild, eststar$Primaryests)
}
CI <- t(apply(theta.b, 2, quantile, c(0.025, 0.975)))
rownames(CI) <- c("WildSmolts", stocks)
colnames(CI) <- c("LCI", "UCI")
cat("Problem 5f: 95% bootstrap CIs\n")
print(round(CI))

# ---- Problem 5g: deliberately trigger Error 1 ----
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
    indivStar <- section5_problem_5d_fish(AllPrime_rare, strats)
    eststar   <- section5_problem_5e_fish(passdata, RearData, indivStar,
                                          Pgrps = rare_stocks, strats = strats)
    theta.b_rare <- numeric(5L)
    theta.b_rare[] <- c(eststar$TotalWild, eststar$Primaryests)
  }
  "no error (lucky path)"
}, error   = function(e) conditionMessage(e),
   warning = function(w) conditionMessage(w))
cat("Problem 5g caught:", caught, "\n")

# ---- Extension: what Error 1 looks like with a fixed ncol ----
# In production SCRAPI, theta.b has a fixed number of columns set before the
# loop.  When Primaryests drops a stock, the vector is shorter than expected
# and R either recycles silently (wrong numbers) or throws a subscript error.
# We reproduce both outcomes.

p_rare     <- 1L + length(rare_stocks)   # expect 5 columns
theta_bad  <- matrix(NA_real_, nrow = 10L, ncol = p_rare)

silent_recycle <- tryCatch({
  for (b in seq_len(10L)) {
    indivStar <- section5_problem_5d_fish(AllPrime_rare, strats)
    eststar   <- section5_problem_5e_fish(passdata, RearData, indivStar,
                                          Pgrps = rare_stocks, strats = strats)
    short_vec <- c(eststar$TotalWild, eststar$Primaryests)
    theta_bad[b, ] <- short_vec   # may recycle if length(short_vec) < p_rare
  }
  "recycled silently — theta_bad contains WRONG values"
}, warning = function(w) paste("warning:", conditionMessage(w)),
   error   = function(e) paste("error:",   conditionMessage(e)))
cat("Extension recycling outcome:", silent_recycle, "\n")
# The fix: always run section5_problem_5g_fish() to catch rare stocks BEFORE
# calling SCRAPI(), or consolidate rare stocks into an "Other" category.

# ---- Forward pointer ----
# Section 6 introduces Collaps (temporal strata for p_n) and shows how
# pooling across strata with different nighttime rates biases a_t upward,
# then reproduces SCRAPI's Collaps assignment loop and its Error 3.
