# Section 2 - MLE and the inverse-sample-rate weighting in SCRAPI
# ---------------------------------------------------------------
# Goal: connect the Binomial MLE  a_d_hat = w / r  to the inverse-
# sample-rate weighting that SCRAPI's thetahat() uses everywhere.
#
# Repo pointer (SCOBI/R/SCRAPI.r, thetahat function):
#   line 75:  dailypass <- passage$Tally / passage$Ptrue
#                                ^ this IS the per-day MLE
#   line 78:  HNCWstrat <- mApply(1/RearDat$True, list(RearDat$Stratum,
#                                                      RearDat$Rear), sum)
#   line 94:  Primarystrata <- mApply(1/Fish$SR, list(Fish$Strat,
#                                                     Fish$PGrp), sum)
#                                ^ each fish counts as 1/SR (Hansen-Hurwitz)
#
# Reading sample: SCOBI/R/SCRAPI.r lines 74-126 (the full thetahat function)

library(ggplot2)
library(dplyr)
library(purrr)
library(tibble)

set.seed(2026)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# --- 1. The closed-form MLE -------------------------------------------------
a_d_true <- 500L
r_sample <- 5/6

w_obs    <- rbinom(1, a_d_true, r_sample)
a_d_hat  <- w_obs / r_sample
cat("w =", w_obs, "  a_d_hat = w / r =", a_d_hat,
    "  (this is SCRAPI line 75)\n\n")

# --- 2. Verify via grid evaluation of the log-likelihood --------------------
a_d_grid <- seq(max(w_obs, 1), 700, by = 1L)
loglik   <- dbinom(w_obs, size = a_d_grid, prob = r_sample, log = TRUE)
mle_grid <- a_d_grid[which.max(loglik)]
cat("Grid MLE =", mle_grid, "   closed-form w/r =", w_obs / r_sample, "\n\n")

p_lik <- ggplot(tibble(a_d = a_d_grid, loglik = loglik),
                aes(a_d, loglik)) +
  geom_line(colour = "steelblue", linewidth = 1) +
  geom_vline(xintercept = w_obs / r_sample,
             colour = "firebrick", linewidth = 1, linetype = "dashed") +
  labs(title = "Section 2 - Binomial log-likelihood for a_d",
       subtitle = paste0("w = ", w_obs, ", r = 5/6; peak at w/r"),
       x = expression(a[d]), y = expression(log~L(a[d])))
ggsave(file.path(plots_dir, "section02_loglik.png"), p_lik,
       width = 6, height = 4, dpi = 150)

# --- 3. optim() on a continuous relaxation (lgamma) -------------------------
neg_loglik_relaxed <- function(a_d, w, r) {
  if (a_d < w) return(Inf)
  -(lgamma(a_d + 1) - lgamma(w + 1) - lgamma(a_d - w + 1) +
    w * log(r) + (a_d - w) * log(1 - r))
}
fit <- optim(par = w_obs * 1.1, fn = neg_loglik_relaxed,
             w = w_obs, r = r_sample,
             method = "Brent", lower = w_obs, upper = 5000)
cat("optim continuous-relaxation MLE =", round(fit$par, 2), "\n\n")

# --- 4. Transition: from adult window count to SCRAPI's smolt trap ---------
# Steps 1-3 above were about ADULT window counts: a_d_hat = w / r.
# SCRAPI is the SMOLT estimator -- different fish, different filter.
# A smolt at the bypass passes TWO filters:
#   (a) trap rate    t_d  (was the gate open?)
#   (b) guidance eff e_sd (did the fish route into the bypass?)
# Combined detection probability p_sd = t_d * e_sd. This `SR` IS that
# smolt-trap probability. The same binomial MLE applies: each sampled
# fish counts as 1/SR fish in the run.
#
# SCRAPI's "AllPrime" data frame has one row per sampled SMOLT:
#   columns: Strat (collapsed week), PGrp (e.g. GenStock), [SGrp,] SR
# It weighs each fish by 1/SR (Hansen-Hurwitz / Horvitz-Thompson form).
SR <- 5/6 * 0.45   # smolt-trap combined detection p_sd = t_d * e_sd
AllPrime <- tibble(
  Strat = sample(c("S1", "S2", "S3"), size = 30,
                 replace = TRUE, prob = c(0.4, 0.4, 0.2)),
  PGrp  = sample(c("LOSALM", "CHMBLN", "IMNAHA"),
                 size = 30, replace = TRUE,
                 prob = c(0.5, 0.35, 0.15)),
  SR    = SR
)

# Reproduce SCRAPI line 94 verbatim with base R `tapply` standing in for plyr::mApply
Primarystrata <- tapply(1 / AllPrime$SR,
                        list(Strat = AllPrime$Strat,
                             PGrp  = AllPrime$PGrp),
                        sum)
Primarystrata[is.na(Primarystrata)] <- 0
cat("Primarystrata (inverse-SR weighted counts per stratum x stock):\n")
print(round(Primarystrata, 1))
cat("\n")

# SCRAPI line 96: column-normalize to get per-stock proportions
Primaryproportions <- prop.table(Primarystrata, margin = 1)
cat("Primaryproportions (row-normalized within stratum):\n")
print(round(Primaryproportions, 3))
cat("\n")

# --- 5. Show that 1/SR weighting equals MLE expansion -----------------------
# If we observed n fish in a sample of N true passing fish, each fish has
# sampling probability SR. The MLE for N is sum(1/SR) over sampled fish,
# which is exactly the inverse-sample-rate sum used in SCRAPI.
n_fish_sampled <- nrow(AllPrime)
N_hat          <- sum(1 / AllPrime$SR)
cat("Inverse-SR expansion: ", n_fish_sampled, "fish sampled at SR =",
    round(SR, 3), "  ->  N_hat =", round(N_hat, 1), "\n")
cat("Equivalent to MLE for a binomial detection process.\n\n")

# --- 6. Replicate to confirm the estimator is unbiased ----------------------
nreps <- 1000L
estimates <- map_dbl(seq_len(nreps),
                     ~ rbinom(1, a_d_true, r_sample) / r_sample)
cat("Across", nreps, "replicates: mean(a_d_hat) =",
    round(mean(estimates), 1), "  (target", a_d_true, ")\n")

# --- 7. End-of-section pointers ---------------------------------------------
# Once Section 2 makes sense you should be able to read every line of
#   SCOBI/R/SCRAPI.r:74-126 (thetahat function)
# without surprise. In particular:
#   - dailypass / Ptrue is the daily MLE
#   - mApply(1/SR, ...) is the inverse-sample-rate (Hansen-Hurwitz) sum
#   - prop.table(Freqs, margin = ...) turns weighted counts into proportions
