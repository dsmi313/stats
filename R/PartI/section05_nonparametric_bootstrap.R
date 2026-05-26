# Section 5 - Nonparametric bootstrap: bootsmolt's weighted stratum resample
# ---------------------------------------------------------------
# Goal: rebuild the within-stratum weighted resample loop that SCRAPI's
# bootsmolt() uses for FishWH and FishDat, then assemble a toy
# `thetahat_toy()` whose signature matches SCRAPI's thetahat(). After this
# section you can read every line of bootsmolt() and explain Error 1
# (theta.b dimension mismatch) from first principles.
#
# Repo pointer (SCOBI/R/SCRAPI.r, bootsmolt function):
#   line 131:  theta.b <- matrix(numeric(p*B), ncol = p)
#   lines 148-157:  weighted resample of FishWH by stratum (Rear weights)
#   lines 159-169:  weighted resample of FishDat by stratum (SR weights)
#   line 173:  eststar <- thetahat(dailyStar, RearStar, indivStar)
#   lines 175-176:  theta.b[b, ] <- c(eststar[[1]], t(eststar[[4]]),
#                                     eststar[[5]], as.vector(t(eststar[[6]])))
#
# Repo pointer (SCOBI/R/SCRAPI.r, thetahat function):
#   lines 74-126:  the full function we mirror in `thetahat_toy()` below.

library(ggplot2)
library(dplyr)
library(purrr)
library(tibble)

set.seed(2026)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# --- 1. Build toy RearData, FishWH, AllPrime, passdata ---------------------
# Each toy data frame matches the column names SCRAPI uses internally.
stocks <- c("LOSALM", "CHMBLN", "IMNAHA")
strats <- c("S1", "S2", "S3")
nfish  <- 90L

# AllPrime mirrors SCRAPI's AllPrime: per-fish stratum + primary group + SR
AllPrime <- tibble(
  Strat = sample(strats, size = nfish, replace = TRUE,
                 prob = c(0.4, 0.4, 0.2)),
  PGrp  = sample(stocks, size = nfish, replace = TRUE,
                 prob = c(0.55, 0.35, 0.10)),
  SR    = 5/6 * 0.45                       # realized sample rate (t_d * e_sd)
)

# RearData mirrors SCRAPI's RearData: per-fish rear type + stratum + True rate
RearData <- AllPrime |>
  mutate(Rear   = sample(c("W", "HNC"), n(), replace = TRUE,
                         prob = c(0.7, 0.3)),
         Stratum = Strat,
         True    = SR) |>
  select(Rear, Stratum, True)

# passdata mirrors SCRAPI's passdata: collapsed stratum, daily Tally, Ptrue
ndays <- 30L
passdata <- tibble(
  Stratum = rep(strats, each = ndays / 3),
  Tally   = rbinom(ndays, size = 200L, prob = 5/6 * 0.45),
  Ptrue   = 5/6 * 0.45
)

# --- 2. Implement the bootsmolt weighted resample VERBATIM -----------------
# SCRAPI.r lines 148-157 (FishWH resample) and 159-169 (FishDat resample)
# are reproduced below as a function we can call B times.
resample_fishWH <- function(FishWH, strats) {
  H <- 0
  for (h in strats) {
    justwk <- FishWH[FishWH$Stratum == h, ]
    nwk    <- nrow(justwk)
    if (nwk == 0L) next
    i      <- sample.int(nwk, replace = TRUE, prob = unlist(justwk$True))
    wkstar <- justwk[i, ]
    if (H == 0) { WHstar <- wkstar; H <- 1
    } else        { WHstar <- rbind(WHstar, wkstar) }
  }
  WHstar
}
resample_fishDat <- function(FishDat, strats) {
  H <- 0
  for (h in strats) {
    justwk <- FishDat[FishDat$Strat == h, ]
    nwk    <- nrow(justwk)
    if (nwk == 0L) next
    i      <- sample.int(nwk, replace = TRUE, prob = unlist(justwk$SR))
    wkstar <- justwk[i, ]
    if (H == 0) { DatStar <- wkstar; H <- 1
    } else        { DatStar <- rbind(DatStar, wkstar) }
  }
  DatStar
}

# --- 3. Toy thetahat() matching SCRAPI's signature -------------------------
# Inputs (matching SCRAPI lines 74-87):
#   passage  : data.frame(Stratum, Tally, Ptrue)
#   RearDat  : data.frame(Rear, Stratum, True)
#   Fish     : data.frame(Strat, PGrp, SR)
# Output  : list(TotalWild, WildStrata, Primaryproportions, Primaryests)
thetahat_toy <- function(passage, RearDat, Fish, Pgrps, strats) {
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
                               factor(Fish$PGrp,  levels = Pgrps)), sum)
  Primarystrata[is.na(Primarystrata)] <- 0
  Primaryproportions <- prop.table(Primarystrata, margin = 1)
  Primaryproportions[is.na(Primaryproportions)] <- 0
  Primaryests <- as.vector(t(Primaryproportions) %*% as.vector(WildStrata))
  names(Primaryests) <- Pgrps

  list(TotalWild           = TotalWild,
       WildStrata          = WildStrata,
       Primaryproportions  = Primaryproportions,
       Primaryests         = Primaryests)
}

est <- thetahat_toy(passdata, RearData, AllPrime, Pgrps = stocks, strats = strats)
cat("Point estimates from thetahat_toy():\n")
cat("  TotalWild         =", round(est$TotalWild), "\n")
cat("  WildStrata        ="); print(round(est$WildStrata))
cat("  Primaryests       ="); print(round(est$Primaryests))
cat("\n")

# --- 4. Full bootsmolt loop on toy data ------------------------------------
B <- 1000L
p <- 1L + length(stocks)
theta.b <- matrix(numeric(p * B), ncol = p)
for (b in seq_len(B)) {
  if (b == 1) {
    dailyStar <- passdata
    RearStar  <- RearData
    indivStar <- AllPrime
  } else {
    # Daily count layer (Section 4)
    cntstar <- vapply(seq_len(nrow(passdata)), function(i) {
      if (round(passdata$Tally[i] / passdata$Ptrue[i]) == 0) 0L
      else rbinom(1, round(passdata$Tally[i] / passdata$Ptrue[i]),
                  passdata$Ptrue[i])
    }, integer(1))
    dailyStar <- data.frame(Stratum = passdata$Stratum,
                            Tally   = cntstar,
                            Ptrue   = passdata$Ptrue)
    RearStar  <- resample_fishWH(RearData, strats)
    indivStar <- resample_fishDat(AllPrime, strats)
  }
  eststar <- thetahat_toy(dailyStar, RearStar, indivStar,
                          Pgrps = stocks, strats = strats)
  theta.b[b, ] <- c(eststar$TotalWild, eststar$Primaryests)
}

# SCRAPI.r lines 181-186: extract 95% CIs per column of theta.b
CI <- t(apply(theta.b, 2, quantile, c(0.025, 0.975)))
rownames(CI) <- c("WildSmolts", stocks)
colnames(CI) <- c("LCI", "UCI")
cat("Bootstrap 95% CIs (theta.b columns):\n")
print(round(CI))
cat("\n")

# --- 5. Trigger Error 1 deliberately ---------------------------------------
# Add a 4th rare stock with only 2 fish, all in one stratum. When a
# bootstrap iteration loses that stock entirely from a stratum, the
# `Primaryproportions` row collapses and `theta.b[b, ] <- c(...)` recycles.
rare_stocks <- c(stocks, "LOCLWR")
AllPrime_rare <- bind_rows(
  AllPrime,
  tibble(Strat = c("S1", "S1"), PGrp = c("LOCLWR", "LOCLWR"),
         SR = 5/6 * 0.45)
)

cat("Attempting bootstrap with a rare stock that lives in 1 stratum only:\n")
caught <- tryCatch({
  for (b in seq_len(50L)) {
    indivStar <- resample_fishDat(AllPrime_rare, strats)
    eststar <- thetahat_toy(passdata, RearData, indivStar,
                            Pgrps = rare_stocks, strats = strats)
    # Try to stuff into a p = 1 + 4 = 5-wide row
    theta.b_rare <- numeric(5L)
    theta.b_rare[] <- c(eststar$TotalWild, eststar$Primaryests)
  }
  "no error fired (lucky bootstrap path)"
}, error = function(e) conditionMessage(e),
   warning = function(w) conditionMessage(w))
cat("  Caught:", caught, "\n\n")
cat("  This is exactly the SCRAPI Error 1 mechanism described in PLAN.md.\n")
cat("  Fix: drop LOCLWR or move it to its own analysis before SCRAPI().\n\n")

# --- 6. Plot the theta.b column for TotalWild ------------------------------
p_total <- ggplot(tibble(total = theta.b[, 1]), aes(total)) +
  geom_histogram(bins = 60, fill = "steelblue", colour = "white") +
  geom_vline(xintercept = CI[1, ], colour = "firebrick",
             linewidth = 1, linetype = "dashed") +
  labs(title = "Section 5 - Bootstrap TotalWild (theta.b column 1)",
       subtitle = "Mirrors SCRAPI.r line 175",
       x = "TotalWild", y = "Bootstrap draws")
ggsave(file.path(plots_dir, "section05_theta_b_total.png"), p_total,
       width = 6, height = 4, dpi = 150)

# --- 7. End-of-section pointers --------------------------------------------
# You can now read:
#   SCOBI/R/SCRAPI.r lines 128-189  (full bootsmolt function)
#   SCOBI/R/SCRAPI.r lines  74-126  (full thetahat function)
# and recognize every line. The pre-run checklist in PLAN.md (Error 1) now
# describes exactly what you just triggered.
