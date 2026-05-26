# Section 4 - Parametric bootstrap: the daily-count layer of bootsmolt()
# ---------------------------------------------------------------
# Goal: implement the EXACT daily-count bootstrap inner loop from
# SCRAPI's bootsmolt() and compare it to the vectorized analogue in
# escapeLGD::expand_wc_binom_night().
#
# Repo pointer (SCOBI/R/SCRAPI.r, bootsmolt function):
#   lines 139-145:
#     dailypass <- round(LGDdaily$Tally / LGDdaily$Ptrue)
#     cntstar = numeric(ndays)
#     for (i in 1:ndays) {
#       if (dailypass[i] != 0) cntstar[i] <- rbinom(1, dailypass[i],
#                                                    LGDdaily$Ptrue[i])
#     }
#     dailyStar <- data.frame(Stratum = LGDdaily$Stratum,
#                              Tally   = cntstar,
#                              Ptrue   = LGDdaily$Ptrue)
#
# Repo pointer (escapeLGD/R/night_fall_reascend_wc_binom.R):
#   lines 145-151:
#     wc_binom[[2]][,i] <- rbinom(boots, wc_binom[[1]]$wc[i], wc_prop) / wc_prop
#   (vectorized: one rbinom call returns all bootstrap draws at once)

library(ggplot2)
library(dplyr)
library(purrr)
library(tibble)

set.seed(2026)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# --- 1. Build a SCRAPI-style LGDdaily data frame ----------------------------
ndays <- 60L
LGDdaily <- tibble(
  Stratum = rep(seq_len(6), each = 10L),
  Tally   = rbinom(ndays, size = 200L, prob = (5/6) * 0.45),
  Ptrue   = 5/6 * 0.45     # = SampleRate x GuidanceEfficiency
)
cat("First 5 rows of LGDdaily (matches SCRAPI passdata layout):\n")
print(head(LGDdaily, 5))
cat("\n")

# --- 2. Implement bootsmolt's daily-count step VERBATIM --------------------
B <- 5000L
bootsmolt_daily <- function(LGDdaily, B) {
  dailypass <- round(LGDdaily$Tally / LGDdaily$Ptrue)
  ndays     <- nrow(LGDdaily)
  out <- matrix(NA_real_, nrow = B, ncol = ndays)
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
  out
}
boot_draws <- bootsmolt_daily(LGDdaily, B)
season_totals <- rowSums(boot_draws)
ci <- quantile(season_totals, c(0.025, 0.975))
cat("Season-total parametric bootstrap (SCRAPI inner loop):\n")
cat("  point      =", round(sum(LGDdaily$Tally / LGDdaily$Ptrue)), "\n")
cat("  95% CI     = [", round(ci[1]), ",", round(ci[2]), "]\n")
cat("  CI width   =", round(diff(ci)), "\n\n")

# --- 3. The escapeLGD vectorized equivalent --------------------------------
# escapeLGD does the same thing in one rbinom call per day. Faster but
# functionally identical.
bootsmolt_daily_vec <- function(LGDdaily, B) {
  dailypass <- round(LGDdaily$Tally / LGDdaily$Ptrue)
  sapply(seq_len(nrow(LGDdaily)), function(i) {
    rbinom(B, dailypass[i], LGDdaily$Ptrue[i]) / LGDdaily$Ptrue[i]
  })
}
t_loop <- system.time({
  res_loop <- bootsmolt_daily(LGDdaily, B = 500L)
})[["elapsed"]]
t_vec  <- system.time({
  res_vec  <- bootsmolt_daily_vec(LGDdaily, B = 500L)
})[["elapsed"]]
cat("Timing on B = 500 bootstraps:\n")
cat("  SCRAPI inner loop (line 142-144):", round(t_loop, 2), "s\n")
cat("  escapeLGD vectorized (line 150):  ", round(t_vec,  2), "s\n\n")

# --- 4. Coverage check ------------------------------------------------------
# Sanity: simulate 500 seasons, compute a 95% bootstrap CI each time, check
# how often the truth lies inside.
ndays_cov <- 20L; Ptrue <- 5/6 * 0.45; true_pass_per_day <- 200L
nseasons  <- 200L
hits <- vapply(seq_len(nseasons), function(s) {
  set.seed(s)
  Tally <- rbinom(ndays_cov, true_pass_per_day, Ptrue)
  LD <- tibble(Stratum = 1L, Tally = Tally, Ptrue = Ptrue)
  draws <- bootsmolt_daily_vec(LD, B = 1000L)
  totals <- rowSums(draws)
  ci  <- quantile(totals, c(0.025, 0.975))
  truth <- ndays_cov * true_pass_per_day
  (truth >= ci[1]) && (truth <= ci[2])
}, logical(1))
cat("95% CI coverage across", nseasons, "seasons =", mean(hits), "\n\n")

# --- 5. Plot the season-total bootstrap distribution -----------------------
p_boot <- ggplot(tibble(total = season_totals), aes(total)) +
  geom_histogram(bins = 60, fill = "steelblue", colour = "white") +
  geom_vline(xintercept = ci, colour = "firebrick",
             linewidth = 1, linetype = "dashed") +
  labs(title = "Section 4 - Bootstrap season totals (bootsmolt inner loop)",
       subtitle = "Dashed lines = 95% CI",
       x = "season total smolt passage", y = "Bootstrap draws")
ggsave(file.path(plots_dir, "section04_boot_season.png"), p_boot,
       width = 6, height = 4, dpi = 150)

# --- 6. End-of-section pointers --------------------------------------------
# You can now read SCRAPI's daily-count bootstrap block at:
#   SCOBI/R/SCRAPI.r lines 128-145 (bootsmolt setup + daily layer)
# Then move to Section 5 to build the within-stratum weighted resample of
# FishWH and FishDat (lines 148-169).
