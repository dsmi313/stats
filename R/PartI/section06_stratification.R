# Section 6 - Stratification, Collaps assignment, and SCRAPI Error 3
# ---------------------------------------------------------------
# Goal: build SCRAPI's collapsing pattern from raw weeks and reproduce the
# date-by-date Collaps assignment that triggers Error 3 when trap and FPC
# date tables disagree.
#
# Repo pointers (SCOBI/R/SCRAPI.r):
#   line 217:  Cpattern <- unique(cbind(pass[, PASSstrat],
#                                       pass[, PASScollaps]))
#   lines 254-259 (the "All$Collaps" assignment loop that fires Error 3):
#     All$Collaps <- numeric(nAll)
#     AllDates    <- unique(All[, FISHdate])
#     for (d in AllDates) {
#       CollStrat <- pass[which(pass[, PASSdate] == d), PASScollaps]
#       All$Collaps[which(All[, FISHdate] == d)] <- CollStrat
#     }

library(ggplot2)
library(dplyr)
library(purrr)
library(tibble)

set.seed(2026)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# --- 1. Build a season with p_n shifting halfway through -------------------
# Mirrors PLAN.md Section 6: pooling lies when p_n is time-varying.
n_weeks      <- 8L
days_per_wk  <- 7L
ndays        <- n_weeks * days_per_wk
p_n_early    <- 0.10
p_n_late     <- 0.40
a_d_per_day  <- 100L
d_a_per_day  <- 30L

season <- tibble(
  Week    = rep(seq_len(n_weeks), each = days_per_wk),
  Collaps = rep(c("EARLY", "LATE"),                # 2-stratum collapse
                each = ndays / 2),
  p_n_true = ifelse(Collaps == "EARLY", p_n_early, p_n_late),
  a_d      = a_d_per_day,
  d_n      = rbinom(ndays, d_a_per_day, p_n_true),
  p_n_hat  = d_n / d_a_per_day,
  a_t_true = a_d / (1 - p_n_true)
)
truth_total <- sum(season$a_t_true)

# SCRAPI line 217 reproduced verbatim:
Cpattern <- unique(cbind(season$Week, season$Collaps))
cat("Cpattern (SCRAPI line 217):\n")
print(Cpattern); cat("\n")

# --- 2. Pooled vs stratified estimator -------------------------------------
p_n_pooled  <- sum(season$d_n) / (nrow(season) * d_a_per_day)
a_t_pooled  <- sum(season$a_d) / (1 - p_n_pooled)

strat_summary <- season |>
  group_by(Collaps) |>
  summarise(p_n_hat = sum(d_n) / (n() * d_a_per_day),
            a_d_sum = sum(a_d), .groups = "drop") |>
  mutate(a_t = a_d_sum / (1 - p_n_hat))
a_t_stratified <- sum(strat_summary$a_t)

cat("Truth        =", round(truth_total), "\n")
cat("Pooled       =", round(a_t_pooled),
    "   bias =", round(a_t_pooled - truth_total), "\n")
cat("Stratified   =", round(a_t_stratified),
    "   bias =", round(a_t_stratified - truth_total), "\n\n")

# --- 3. Reproduce the SCRAPI Collaps assignment loop -----------------------
# SCRAPI lines 254-259 - we replicate the loop on toy trap data.
nFish <- 200L
fish_weeks <- sample(season$Week, size = nFish, replace = TRUE)
All <- tibble(
  CollectionDate = format(as.Date("2024-04-01") + (fish_weeks - 1L) * 7L,
                          "%m/%d/%Y"),
  GenStock       = sample(c("LOSALM", "CHMBLN", "IMNAHA"),
                          size = nFish, replace = TRUE)
)
pass <- season |>
  mutate(SampleEndDate = format(as.Date("2024-04-01") + (Week - 1L) * 7L,
                                "%m/%d/%Y")) |>
  group_by(Week, Collaps, SampleEndDate) |>
  summarise(SampleCount = sum(a_d), .groups = "drop")

# Verbatim port of SCRAPI lines 254-259:
nAll <- nrow(All)
All$Collaps <- character(nAll)
AllDates    <- unique(All$CollectionDate)
for (d in AllDates) {
  CollStrat <- pass$Collaps[which(pass$SampleEndDate == d)]
  All$Collaps[which(All$CollectionDate == d)] <- CollStrat
}
cat("First 5 rows of All with Collaps assigned (mirrors SCRAPI line 258):\n")
print(head(All, 5)); cat("\n")

# --- 4. Trigger Error 3 deliberately ---------------------------------------
# Error 3 fires when a trap date matches multiple FPC rows so CollStrat has
# length > 1. We duplicate one row of `pass` and rerun the loop.
pass_bad <- bind_rows(pass, pass[1, ])
cat("Attempting Collaps assignment with a duplicated FPC date:\n")
caught <- tryCatch({
  for (d in AllDates) {
    CollStrat <- pass_bad$Collaps[which(pass_bad$SampleEndDate == d)]
    All$Collaps[which(All$CollectionDate == d)] <- CollStrat
  }
  "no error fired"
}, warning = function(w) conditionMessage(w),
   error   = function(e) conditionMessage(e))
cat("  Caught:", caught, "\n")
cat("  This is the SCRAPI Error 3 mechanism (Collaps recycling).\n")
cat("  Fix: ensure 1:1 date correspondence between trap data and FPC.\n\n")

# --- 5. Parametric-bootstrap CI for pooled vs stratified -------------------
boot_season <- function(season, B, stratified) {
  d_a <- 30L
  estimates <- numeric(B)
  if (stratified) {
    p_strata <- season |>
      group_by(Collaps) |>
      summarise(p_n_hat = sum(d_n) / (n() * d_a), .groups = "drop")
    for (b in seq_len(B)) {
      smpl <- season |>
        left_join(p_strata, by = "Collaps", suffix = c("", ".s")) |>
        mutate(d_n_b = rbinom(n(), size = d_a, prob = p_n_hat.s))
      strat_b <- smpl |>
        group_by(Collaps) |>
        summarise(p_n_b   = sum(d_n_b) / (n() * d_a),
                  a_d_sum = sum(a_d), .groups = "drop")
      estimates[b] <- sum(strat_b$a_d_sum / (1 - strat_b$p_n_b))
    }
  } else {
    p_pool <- sum(season$d_n) / (nrow(season) * d_a)
    for (b in seq_len(B)) {
      d_n_b <- rbinom(nrow(season), size = d_a, prob = p_pool)
      p_n_b <- sum(d_n_b) / (nrow(season) * d_a)
      estimates[b] <- sum(season$a_d) / (1 - p_n_b)
    }
  }
  list(ci = quantile(estimates, c(0.025, 0.975)), draws = estimates)
}

set.seed(2026)
bp <- boot_season(season, B = 1500L, stratified = FALSE)
bs <- boot_season(season, B = 1500L, stratified = TRUE)
cat("Parametric bootstrap CIs (B = 1500):\n")
cat("  pooled      = [", round(bp$ci[1]), ",", round(bp$ci[2]), "]",
    "   width =", round(diff(bp$ci)), "\n")
cat("  stratified  = [", round(bs$ci[1]), ",", round(bs$ci[2]), "]",
    "   width =", round(diff(bs$ci)), "\n\n")

# --- 6. Visualize pooled bias vs stratified spread -------------------------
draws_tbl <- bind_rows(
  tibble(method = "pooled",     a_t = bp$draws),
  tibble(method = "stratified", a_t = bs$draws)
)
p_strat <- ggplot(draws_tbl, aes(a_t, fill = method)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = truth_total, colour = "black", linewidth = 1) +
  labs(title = "Section 6 - Pooled vs stratified estimator",
       subtitle = "Solid line = season truth",
       x = expression(hat(a)[t]^{season}),
       y = "Replicates", fill = NULL)
ggsave(file.path(plots_dir, "section06_pool_vs_strat.png"), p_strat,
       width = 6, height = 4, dpi = 150)

# --- 7. End-of-section pointers --------------------------------------------
# You can now read:
#   SCOBI/R/SCRAPI.r lines 215-225 (collapsing pattern, Cpattern)
#   SCOBI/R/SCRAPI.r lines 254-269 (Collaps and true-rate assignment loops)
# Both of these are where Error 3 originates. The pre-run checklist in
# PLAN.md (no extra FPC rows; 1:1 date correspondence) is now self-evident.
