# Section 6 - Stratification: when pooling lies
# ---------------------------------------------------------------
# If a parameter shifts across the run, a pooled estimate is biased.
# Stratifying removes the bias at the cost of wider CIs. This justifies
# the three EASE strata schemes and the ">=100 wild fish" rule from SCOBI.
#
# Reading: PLAN.md Section 6, MIT Class 17a (Frequentist School),
# Class 22 (Confidence Intervals), R Studio 8 and 9.

library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)

set.seed(2026)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# --- 1. Simulate a season with shifting p_n --------------------------------
# Early half: low nighttime passage. Late half: high.
n_days       <- 60L
p_n_early    <- 0.10
p_n_late     <- 0.40
a_d_per_day  <- 100L     # daytime escapement per day (known for clarity)
d_a_per_day  <- 30L      # PIT-tagged adults per day in the proportion sample

season <- tibble(day = seq_len(n_days),
                 stratum  = ifelse(day <= n_days / 2, "early", "late"),
                 p_n_true = ifelse(stratum == "early", p_n_early, p_n_late))
season <- season |>
  mutate(a_d      = a_d_per_day,
         d_n      = rbinom(n(), size = d_a_per_day, prob = p_n_true),
         p_n_hat  = d_n / d_a_per_day,
         a_t_true = a_d / (1 - p_n_true))

truth_total <- sum(season$a_t_true)

# --- 2. Pooled estimate (biased) -------------------------------------------
p_n_pooled  <- sum(season$d_n) / (nrow(season) * d_a_per_day)
a_t_pooled  <- sum(season$a_d) / (1 - p_n_pooled)

# --- 3. Stratified estimate (unbiased) -------------------------------------
strat_summary <- season |>
  group_by(stratum) |>
  summarise(p_n_hat = sum(d_n) / (n() * d_a_per_day),
            a_d_sum = sum(a_d),
            .groups = "drop") |>
  mutate(a_t = a_d_sum / (1 - p_n_hat))
a_t_stratified <- sum(strat_summary$a_t)

cat("Pooled vs stratified vs truth:\n")
cat("  truth         =", round(truth_total), "\n")
cat("  pooled  a_t   =", round(a_t_pooled),
    "   bias =", round(a_t_pooled - truth_total), "\n")
cat("  stratified    =", round(a_t_stratified),
    "   bias =", round(a_t_stratified - truth_total), "\n\n")

# --- 4. CI width via parametric bootstrap: pooled vs stratified -----------
boot_season <- function(season, B, stratified) {
  d_a <- 30L
  estimates <- numeric(B)
  if (stratified) {
    p_strata <- season |>
      group_by(stratum) |>
      summarise(p_n_hat = sum(d_n) / (n() * d_a), .groups = "drop")
    for (b in seq_len(B)) {
      smpl <- season |>
        left_join(p_strata, by = "stratum") |>
        mutate(d_n_b = rbinom(n(), size = d_a, prob = p_n_hat))
      strat_b <- smpl |>
        group_by(stratum) |>
        summarise(p_n_b   = sum(d_n_b) / (n() * d_a),
                  a_d_sum = sum(a_d),
                  .groups = "drop")
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
boot_pool  <- boot_season(season, B = 1500, stratified = FALSE)
boot_strat <- boot_season(season, B = 1500, stratified = TRUE)

cat("Parametric bootstrap CIs (B = 1500):\n")
cat("  pooled      CI =", round(boot_pool$ci),
    "  width =", round(diff(boot_pool$ci)), "\n")
cat("  stratified  CI =", round(boot_strat$ci),
    "  width =", round(diff(boot_strat$ci)), "\n\n")

# --- 5. Stratum sample-size sweep -----------------------------------------
# How wide is the stratified CI as the late stratum PIT sample shrinks?
sweep_sizes <- c(5L, 10L, 20L, 50L)
set.seed(2026)
sweep <- map_dfr(sweep_sizes, function(d_a_late) {
  s2 <- season |>
    mutate(d_a_eff = ifelse(stratum == "late", d_a_late, 30L),
           d_n     = rbinom(n(), size = d_a_eff, prob = p_n_true),
           p_n_hat = d_n / d_a_eff)
  draws <- numeric(1000L)
  p_strata <- s2 |>
    group_by(stratum) |>
    summarise(p_n_hat = sum(d_n) / sum(d_a_eff), .groups = "drop")
  for (b in seq_along(draws)) {
    sb <- s2 |>
      left_join(p_strata, by = "stratum", suffix = c("", ".s")) |>
      mutate(d_n_b = rbinom(n(), size = d_a_eff, prob = p_n_hat.s))
    strat_b <- sb |>
      group_by(stratum) |>
      summarise(p_n_b   = sum(d_n_b) / sum(d_a_eff),
                a_d_sum = sum(a_d),
                .groups = "drop")
    draws[b] <- sum(strat_b$a_d_sum / (1 - strat_b$p_n_b))
  }
  ci <- quantile(draws, c(0.025, 0.975))
  tibble(d_a_late = d_a_late,
         lower    = unname(ci[1]),
         upper    = unname(ci[2]),
         width    = unname(diff(ci)))
})

cat("Stratified CI width as the late-stratum PIT sample shrinks:\n")
print(sweep); cat("\n")

p_sweep <- ggplot(sweep, aes(d_a_late, width)) +
  geom_line(colour = "steelblue") +
  geom_point(size = 2) +
  labs(title = "Section 6 - Stratified CI width vs late-stratum sample",
       x = "d_a in late stratum (per day)",
       y = "Bootstrap 95% CI width")
ggsave(file.path(plots_dir, "section06_strat_sweep.png"), p_sweep,
       width = 6, height = 4, dpi = 150)

# --- 6. Plot the pooled bias and stratified spread ------------------------
draws_tbl <- bind_rows(
  tibble(method = "pooled",     a_t = boot_pool$draws),
  tibble(method = "stratified", a_t = boot_strat$draws)
)

p_strat <- ggplot(draws_tbl, aes(a_t, fill = method)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = truth_total,
             colour = "black", linewidth = 1) +
  labs(title = "Section 6 - Pooling vs stratification",
       subtitle = "Solid line = season truth",
       x = expression(hat(a)[t]^{season}),
       y = "Replicates", fill = NULL)
ggsave(file.path(plots_dir, "section06_pool_vs_strat.png"), p_strat,
       width = 6, height = 4, dpi = 150)

# Section 6 payoff ----------------------------------------------------------
# The ">=100 wild fish per stratum" rule used in SCOBI is derived from this
# bias-variance tradeoff. After Section 6 you can defend the rule from
# first principles instead of citing the SOP.
