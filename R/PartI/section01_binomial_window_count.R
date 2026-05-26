# Section 1 - Binomial window count
# ---------------------------------------------------------------
# MIT 18.05 R Studio 2 (Binomial Distributions) translated into the
# Lower Granite ladder window. Every coin flip in MIT Studio 2 is a
# fish passing the ladder during the sampled portion of the hour.
#
# EASE/SCRAPI mapping:
#   a_d       total daytime adults that actually passed (truth)
#   r         window sampling fraction (e.g. 50 min / 60 min = 5/6)
#   w         observed window count = rbinom(1, a_d, r)
#   a_d_hat   estimator = w / r
#
# Reading: PLAN.md Section 1, MIT Class 4a (Discrete RVs), R Studio 2.
# Run from the repository root: source("R/PartI/section01_binomial_window_count.R")

library(ggplot2)
library(dplyr)
library(purrr)

set.seed(2026)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# --- 1. One season, one window count ----------------------------------------
a_d_true <- 500L      # truth: 500 fish passed the ladder window today
r_sample <- 5/6       # we sample 50 minutes out of every hour

w_obs   <- rbinom(1, size = a_d_true, prob = r_sample)
a_d_hat <- w_obs / r_sample

cat("Single season:\n")
cat("  truth a_d   =", a_d_true, "\n")
cat("  window w    =", w_obs, "\n")
cat("  a_d_hat     =", a_d_hat, "\n\n")

# --- 2. MIT Studio 2, Problem 1a (fish version) -----------------------------
# Estimate P(W = k) and P(W <= k) for a single hour by simulation.
studio2_problem_1a_fish <- function(n_fish, r_sample, k, ntrials) {
  ws <- rbinom(ntrials, size = n_fish, prob = r_sample)
  list(p_eq = mean(ws == k),
       p_le = mean(ws <= k))
}

# --- 3. MIT Studio 2, Problem 1b (fish version) -----------------------------
# Exact P(W = k) using choose() rather than dbinom() so the formula stays
# visible. This is the binomial pmf in its primitive form.
studio2_problem_1b_fish <- function(n_fish, r_sample, k) {
  choose(n_fish, k) * r_sample^k * (1 - r_sample)^(n_fish - k)
}

est_60   <- studio2_problem_1a_fish(n_fish = 60, r_sample = r_sample,
                                    k = 50, ntrials = 1e5)
exact_60 <- studio2_problem_1b_fish(n_fish = 60, r_sample = r_sample, k = 50)
dbinom_60 <- dbinom(50, size = 60, prob = r_sample)

cat("Single-hour pmf check (60 fish, k = 50):\n")
cat("  simulated P(W = 50) =", round(est_60$p_eq, 4),
    "  P(W <= 50) =", round(est_60$p_le, 4), "\n")
cat("  exact     P(W = 50) =", round(exact_60, 4),
    "  (dbinom = ", round(dbinom_60, 4), ")\n\n")

# --- 4. MIT Studio 2, Problem 2 (fish version) ------------------------------
# A reward-and-penalty rule for the trap supervisor: for each marked fish in
# a daily sample of size 10, payoff is k^2 - 7k. Decide whether the rule
# nets the supervisor money when the mark rate is 0.6.
fish_payoff <- function(k) k^2 - 7 * k

payoff_curve <- tibble(k = 0:10, payoff = fish_payoff(0:10))
p_payoff <- ggplot(payoff_curve, aes(k, payoff)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  labs(title = "Section 1 - Daily marked-fish payoff curve",
       subtitle = "Mark rate p = 0.6, ntosses = 10 marked-or-not draws",
       x = "k = marked fish in daily draw", y = "payoff ($)")
ggsave(file.path(plots_dir, "section01_payoff.png"), p_payoff,
       width = 6, height = 4, dpi = 150)

mark_rate <- 0.6
n_draws   <- 10
exact_value <- sum(fish_payoff(0:n_draws) *
                   dbinom(0:n_draws, size = n_draws, prob = mark_rate))
sim_value   <- mean(fish_payoff(rbinom(1e5, n_draws, mark_rate)))
cat("Marked-fish payoff (Problem 2 in fish form):\n")
cat("  exact expected payoff      =", round(exact_value, 3), "\n")
cat("  simulated  average payoff  =", round(sim_value, 3), "\n")
cat("  decision: ", ifelse(exact_value > 0, "good bet", "bad bet"), "\n\n")

# --- 5. Replicate 1000 seasons and verify a_d_hat centers on truth ---------
nreps <- 1000L
sims  <- map_dbl(seq_len(nreps),
                 ~ rbinom(1, a_d_true, r_sample) / r_sample)

cat("Replicated season estimator:\n")
cat("  mean(a_d_hat) =", round(mean(sims), 1),
    "  (target", a_d_true, ")\n")
cat("  sd(a_d_hat)   =", round(sd(sims), 1), "\n\n")

# --- 6. Histogram of the simulated estimator distribution ------------------
p_hist <- ggplot(tibble(a_d_hat = sims), aes(a_d_hat)) +
  geom_histogram(bins = 30, fill = "steelblue", colour = "white") +
  geom_vline(xintercept = a_d_true,
             colour = "firebrick", linewidth = 1) +
  labs(title = "Section 1 - Binomial window-count estimator",
       subtitle = paste0("a_d = ", a_d_true, ", r = 5/6, ", nreps, " seasons"),
       x = expression(hat(a)[d]),
       y = "Replicates")
ggsave(file.path(plots_dir, "section01_estimator_hist.png"), p_hist,
       width = 6, height = 4, dpi = 150)

# --- 7. Exact pmf for the daily window count -------------------------------
pmf_tbl <- tibble(k = 380:460,
                  prob = dbinom(380:460, size = a_d_true, prob = r_sample))
p_pmf <- ggplot(pmf_tbl, aes(k, prob)) +
  geom_col(fill = "steelblue") +
  labs(title = "Section 1 - Exact pmf of W ~ Binomial(500, 5/6)",
       x = "window count w", y = "P(W = w)")
ggsave(file.path(plots_dir, "section01_window_pmf.png"), p_pmf,
       width = 6, height = 4, dpi = 150)

# Section 1 payoff -----------------------------------------------------------
# Every later EASE estimator that has the form (count / known proportion)
# is a binomial MLE. This file is the template for that pattern.
