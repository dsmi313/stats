# Section 1 - Binomial window count
# ---------------------------------------------------------------
# Goal: simulate the window count and reproduce the EXACT passage-data
# layout that SCRAPI() and escapeLGD::expand_wc_binom_night() consume.
# Every column name below matches the production source.
#
# Repo pointer (SCOBI/R/SCRAPI.r):
#   line 227:  pass$true      <- pass[, samrate] * pass[, guidance]
#   line 228:  pass$estimated <- pass[, tally]   / pass$true
#   line 232:  passdata <- data.frame(Stratum = ..., Tally = ..., Ptrue = ...)
#
# Repo pointer (escapeLGD/R/night_fall_reascend_wc_binom.R):
#   lines 143-150:  wc_binom <- list(wc %>% mutate(wc = round(wc / wc_prop)));
#                   wc_binom[[2]][,i] <- rbinom(boots, wc[i], wc_prop) / wc_prop
#
# Run from the repository root.

library(ggplot2)
library(dplyr)
library(purrr)
library(tibble)

set.seed(2026)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# --- 1. One day at the ladder window ----------------------------------------
# a_d = true daytime adults; r = sampling fraction; w = observed window count.
a_d_true <- 500L
r_sample <- 5/6

w_obs   <- rbinom(1, size = a_d_true, prob = r_sample)
a_d_hat <- w_obs / r_sample
cat("Single day: w =", w_obs, "  a_d_hat = w/r =", a_d_hat, "\n\n")

# --- 2. Build a season-long `pass` table with SCRAPI's exact column names ---
# SCRAPI's `passageData` input uses these columns:
#   Week, SampleEndDate, SampleCount, SampleRate, GuidanceEfficiency, Collapse
# We simulate a 12-week MY2024-ish season for clarity.
n_weeks <- 12L
days_per_week <- 7L
ndays <- n_weeks * days_per_week

pass <- tibble(
  Week               = rep(seq_len(n_weeks), each = days_per_week),
  SampleEndDate      = format(seq.Date(as.Date("2024-04-01"),
                                       by = "day", length.out = ndays),
                              "%m/%d/%Y"),
  SampleRate         = 5/6,                          # sampling fraction r
  GuidanceEfficiency = 0.45,                         # toy e_sd
  SampleCount        = rbinom(ndays, size = 200L,
                              prob = (5/6) * 0.45),  # daily trap count
  Collapse           = rep(seq_len(n_weeks / 2),     # collapse two weeks
                           each = 2 * days_per_week) # into one stratum
)

# Lines 227-228 of SCRAPI.r reproduced verbatim:
pass$true      <- pass$SampleRate * pass$GuidanceEfficiency
pass$estimated <- pass$SampleCount / pass$true

cat("First 4 rows of the SCRAPI-style pass table:\n")
print(head(pass, 4))
cat("\n")

# Line 232 of SCRAPI.r reproduced verbatim:
passdata <- data.frame(Stratum = pass$Collapse,
                       Tally   = pass$SampleCount,
                       Ptrue   = pass$true)

# Seasonal sums by collapsed stratum -- matches SCRAPI lines 241-244:
passcollaps <- tapply(pass$estimated, pass$Collapse, sum)
cat("Estimate of total smolts by collapsed stratum (SCRAPI line 242):\n")
print(round(passcollaps))
cat("\nTotal smolts =", round(sum(passcollaps)), "\n\n")

# --- 3. Build an escapeLGD-style `wc` tibble --------------------------------
# escapeLGD::expand_wc_binom_night() takes a tibble with two columns:
#   sWeek, wc  (sWeek = statistical week, wc = raw window count)
# escapeLGD then expands by wc_prop:  wc <- round(wc / wc_prop)
wc_prop <- 5/6
wc <- pass |>
  group_by(Week) |>
  summarise(wc = sum(rbinom(n(), size = 500L, prob = wc_prop)),
            .groups = "drop") |>
  rename(sWeek = Week)
cat("First rows of the escapeLGD-style wc tibble:\n")
print(head(wc, 4))
cat("\n")

# escapeLGD line 143:
wc_expanded <- wc |> mutate(wc = round(wc / wc_prop))
cat("After expansion by wc_prop = 5/6 (escapeLGD line 143):\n")
print(head(wc_expanded, 4))
cat("\n")

# --- 4. Replicated simulation to confirm a_d_hat = w/r is unbiased ----------
nreps <- 1000L
sims  <- map_dbl(seq_len(nreps),
                 ~ rbinom(1, a_d_true, r_sample) / r_sample)
cat("Across", nreps, "simulated days:\n")
cat("  mean(a_d_hat) =", round(mean(sims), 1),
    "  (target", a_d_true, ")\n")
cat("  sd(a_d_hat)   =", round(sd(sims), 1), "\n\n")

# --- 5. Plot the simulated estimator distribution ---------------------------
p_hist <- ggplot(tibble(a_d_hat = sims), aes(a_d_hat)) +
  geom_histogram(bins = 30, fill = "steelblue", colour = "white") +
  geom_vline(xintercept = a_d_true, colour = "firebrick", linewidth = 1) +
  labs(title = "Section 1 - Window-count estimator a_d_hat = w/r",
       subtitle = paste0("a_d = ", a_d_true, ", r = 5/6, ", nreps, " days"),
       x = expression(hat(a)[d]), y = "Replicates")
ggsave(file.path(plots_dir, "section01_estimator_hist.png"), p_hist,
       width = 6, height = 4, dpi = 150)

# --- 6. Plot the exact pmf of W ---------------------------------------------
pmf_tbl <- tibble(k = 380:460,
                  prob = dbinom(380:460, size = a_d_true, prob = r_sample))
p_pmf <- ggplot(pmf_tbl, aes(k, prob)) +
  geom_col(fill = "steelblue") +
  labs(title = "Section 1 - Exact pmf of W ~ Binomial(500, 5/6)",
       x = "window count w", y = "P(W = w)")
ggsave(file.path(plots_dir, "section01_window_pmf.png"), p_pmf,
       width = 6, height = 4, dpi = 150)

# --- 7. End-of-section: what to read next in the production repos -----------
# After this section you can open and follow these source blocks line by line:
#   SCOBI/R/SCRAPI.r       lines 196-247  (passage data prep)
#   escapeLGD/R/night_fall_reascend_wc_binom.R lines 130-150 (wc expansion)
