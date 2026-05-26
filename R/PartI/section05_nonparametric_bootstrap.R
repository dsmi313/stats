# Section 5 - Nonparametric bootstrap and distributional diagnostic
# ---------------------------------------------------------------
# Parametric bootstrap (Section 4) assumes a model. Nonparametric
# bootstrap resamples the actual data. EASE uses nonparametric for
# composition uncertainty. This section also builds the dispersion
# diagnostic that will be run on real MY2024 counts in Section 14.
#
# Reading: PLAN.md Section 5, MIT Class 23a/23b, Class 24, R Studio 10.

library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)

set.seed(2026)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# --- 1. Toy trap sample with stock labels ----------------------------------
# PTAGIS-style 6-letter stock codes.
stocks <- c("LOSALM", "CHMBLN", "IMNAHA")
true_p <- c(LOSALM = 0.50, CHMBLN = 0.35, IMNAHA = 0.15)
n_trap <- 30L
fish   <- sample(stocks, size = n_trap, replace = TRUE, prob = true_p)
cat("Observed trap composition:\n"); print(table(fish)); cat("\n")

# --- 2. Nonparametric bootstrap of stock proportions -----------------------
np_boot_props <- function(labels, B = 5000) {
  n     <- length(labels)
  uniq  <- sort(unique(labels))
  draws <- replicate(B, {
    boot <- sample(labels, size = n, replace = TRUE)
    vapply(uniq, function(s) mean(boot == s), numeric(1))
  })
  cis <- apply(draws, 1, quantile, c(0.025, 0.975))
  list(point = vapply(uniq, function(s) mean(labels == s), numeric(1)),
       cis   = cis,
       draws = draws)
}

np <- np_boot_props(fish, B = 5000)
cat("Nonparametric 95% CIs for stock proportions:\n")
print(round(rbind(point = np$point, np$cis), 3)); cat("\n")

# --- 3. Parametric (multinomial) bootstrap for comparison -----------------
param_boot_props <- function(labels, B = 5000) {
  uniq  <- sort(unique(labels))
  p_hat <- vapply(uniq, function(s) mean(labels == s), numeric(1))
  draws <- rmultinom(B, size = length(labels), prob = p_hat) / length(labels)
  cis   <- apply(draws, 1, quantile, c(0.025, 0.975))
  list(point = p_hat, cis = cis)
}

pp <- param_boot_props(fish, B = 5000)
cat("Parametric (multinomial) 95% CIs for the same data:\n")
print(round(rbind(point = pp$point, pp$cis), 3)); cat("\n")

# --- 4. Stratified nonparametric bootstrap ---------------------------------
# Two time strata; resample within each; combine using stratum weights.
strat_boot <- function(labels, strata, weights, B = 2000) {
  uniq_s <- sort(unique(labels))
  draws  <- replicate(B, {
    pooled <- map_dfr(unique(strata), function(stratum) {
      idx  <- which(strata == stratum)
      boot <- sample(labels[idx], size = length(idx), replace = TRUE)
      tibble(stratum = stratum,
             stock   = uniq_s,
             p       = vapply(uniq_s, function(s) mean(boot == s), numeric(1)),
             w       = weights[[stratum]])
    })
    pooled |>
      group_by(stock) |>
      summarise(p = sum(p * w) / sum(w), .groups = "drop") |>
      pull(p)
  })
  cis <- apply(draws, 1, quantile, c(0.025, 0.975))
  list(point = rowMeans(draws), cis = cis)
}

stratum  <- c(rep("early", 15), rep("late", 15))
weights  <- list(early = 0.6, late = 0.4)
strat    <- strat_boot(fish, stratum, weights, B = 1500)
cat("Stratified nonparametric 95% CIs (early/late weights 0.6/0.4):\n")
print(round(rbind(point = strat$point, strat$cis), 3)); cat("\n")

# --- 5. Edge case: rare stock in a small stratum ---------------------------
# Replace the late stratum with 5 fish, all LOSALM. Watch IMNAHA CI collapse.
small_labels <- c(fish[1:15], rep("LOSALM", 5))
small_strata <- c(rep("early", 15), rep("late", 5))
small_weights <- list(early = 0.75, late = 0.25)
strat_small  <- strat_boot(small_labels, small_strata, small_weights, B = 1500)
cat("Edge case (late stratum has 5 fish, all LOSALM):\n")
print(round(rbind(point = strat_small$point, strat_small$cis), 3))
cat("Notice IMNAHA CI is degenerate -> this is the same failure mode\n")
cat("that triggers the SCRAPI theta.b dimension error (Diagnostic Error 1).\n\n")

# --- 6. Distributional diagnostic: compare Binomial to overdispersed -------
ndays      <- 30L
trials     <- 500L
p_true     <- 5/6
mu         <- trials * p_true
bin_counts <- rbinom(ndays, size = trials, prob = p_true)

target_var <- 4 * trials * p_true * (1 - p_true)  # 4x the Binomial variance
nb_size    <- mu^2 / (target_var - mu)
nb_counts  <- rnbinom(ndays, mu = mu, size = nb_size)

diag_tbl <- tibble(
  source = rep(c("Binomial", "Overdispersed (NB)"), each = ndays),
  count  = c(bin_counts, nb_counts)
)
cat("Mean and variance by source:\n")
print(diag_tbl |>
        group_by(source) |>
        summarise(mean = mean(count), var = var(count), .groups = "drop"))
cat("\n")

# --- 7. Chi-squared GoF diagnostic (reusable) -----------------------------
# Use as a screening diagnostic, not a strict hypothesis test.
bin_diagnostic <- function(counts, n_trials, p) {
  rng        <- min(counts):max(counts)
  observed   <- as.numeric(table(factor(counts, levels = rng)))
  expected_p <- dbinom(rng, size = n_trials, prob = p)
  fit <- chisq.test(x = observed,
                    p = expected_p / sum(expected_p),
                    rescale.p = TRUE,
                    simulate.p.value = TRUE, B = 2000)
  list(chi = unname(fit$statistic), p_value = fit$p.value)
}

diag_bin <- bin_diagnostic(bin_counts, trials, p_true)
diag_nb  <- bin_diagnostic(nb_counts,  trials, p_true)
cat("Chi-squared GoF against Binomial(500, 5/6):\n")
cat("  Binomial sample      : chi =", round(diag_bin$chi, 2),
    "  p =", round(diag_bin$p_value, 3), "\n")
cat("  Overdispersed sample : chi =", round(diag_nb$chi, 2),
    "  p =", round(diag_nb$p_value, 3), "\n\n")

# --- 8. Plot dispersion of the two samples --------------------------------
p_disp <- ggplot(diag_tbl, aes(count, fill = source)) +
  geom_histogram(bins = 25, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = mu, colour = "black", linetype = "dashed") +
  labs(title = "Section 5 - Dispersion diagnostic",
       subtitle = "Same mean; NB has fatter tails",
       x = "daily count", y = "Replicates", fill = NULL)
ggsave(file.path(plots_dir, "section05_dispersion.png"), p_disp,
       width = 6, height = 4, dpi = 150)

# Section 5 payoff ----------------------------------------------------------
# Stratified nonparametric resampling is the skeleton of bootsmolt().
# bin_diagnostic() is the function used on real MY2024 w and p_n samples
# in Section 14 to flag any overdispersion for the Section 21 memo.
