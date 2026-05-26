# Section 1 - Binomial window count (adults at the ladder)
# ---------------------------------------------------------------
# Window counts are ADULTS passing the ladder window. The estimator is
#     a_d_hat = w / r
# where r is the sampling fraction (e.g. 50 min / 60 min = 5/6).
# There is NO guidance efficiency here -- GE is a smolt-bypass concept
# (it shows up in SCRAPI's smolt pass table in Part III, not Part I).
#
# Repo pointer (escapeLGD/R/night_fall_reascend_wc_binom.R):
#   line 143:  wc_binom <- list(wc %>% mutate(wc = round(wc / wc_prop)))
#   line 150:  wc_binom[[2]][,i] <- rbinom(boots, wc[i], wc_prop) / wc_prop
#
# Where SCRAPI's GE comes in (smolts, not adults):
#   SCOBI/R/SCRAPI.r line 227:  pass$true <- SampleRate * GuidanceEfficiency
# That formula is the smolt-bypass combined detection probability
# (trap rate x guidance efficiency). It is the right model for the smolt
# trap, NOT for adult window counts at the ladder. Part III revisits it.

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
cat("Single day:\n")
cat("  truth a_d  =", a_d_true, "\n")
cat("  window w   =", w_obs, "\n")
cat("  a_d_hat    = w / r =", a_d_hat, "\n\n")

# --- 2. Season-long adult window counts in escapeLGD's `wc` layout ----------
# escapeLGD::expand_wc_binom_night() takes:
#   wc       tibble with columns (sWeek, wc) -- raw window count per week
#   wc_prop  scalar -- the sampling fraction r (one number, not a column)
# No guidance efficiency anywhere. Adults at the ladder pass one filter
# (the window sampling fraction), not two.
n_weeks <- 12L
wc_prop <- 5/6                # the only sampling fraction relevant to adults

wc <- tibble(
  sWeek = seq_len(n_weeks),
  # weekly truth varies; each week we observe a Binomial(truth, wc_prop)
  truth = round(runif(n_weeks, 2500, 5000)),
  wc    = NA_integer_
)
wc$wc <- rbinom(n_weeks, size = wc$truth, prob = wc_prop)

cat("escapeLGD-style wc tibble (adult window counts, no GE):\n")
print(head(wc |> select(sWeek, wc), 6))
cat("\n")

# escapeLGD line 143 reproduced verbatim:
wc_expanded <- wc |> mutate(wc_expanded = round(wc / wc_prop))
cat("After expansion by wc_prop = 5/6 (escapeLGD line 143):\n")
print(head(wc_expanded |> select(sWeek, wc, wc_expanded, truth), 6))
cat("\n")
cat("Season expanded total (escapeLGD method) =",
    sum(wc_expanded$wc_expanded), "\n")
cat("Season truth                              =",
    sum(wc$truth), "\n\n")

# --- 3. Vectorized parametric bootstrap, escapeLGD style --------------------
# escapeLGD line 150 reproduced verbatim:
#   wc_binom[[2]][,i] <- rbinom(boots, wc[i], wc_prop) / wc_prop
boots <- 5000L
boot_mat <- vapply(seq_len(nrow(wc)),
                   function(i) rbinom(boots, wc$wc[i], wc_prop) / wc_prop,
                   numeric(boots))
season_totals <- rowSums(boot_mat)
ci <- quantile(season_totals, c(0.025, 0.975))
cat("escapeLGD-style parametric bootstrap (no GE involved):\n")
cat("  point estimate (sum(wc / r)) =", sum(wc$wc / wc_prop), "\n")
cat("  95% CI                        = [", round(ci[1]), ",",
    round(ci[2]), "]\n\n")

# --- 4. Replicated simulation: a_d_hat = w/r is unbiased --------------------
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
  labs(title = "Section 1 - Window-count estimator a_d_hat = w / r",
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

# --- 7. End-of-section: what to read next -----------------------------------
# After this section you can read end-to-end:
#   escapeLGD/R/night_fall_reascend_wc_binom.R lines 130-150
#     (expand_wc_binom_night: adult window count expansion + bootstrap)
#
# DO NOT confuse this with SCRAPI's pass table. SCRAPI is the SMOLT
# estimator; its `pass$true = SampleRate * GuidanceEfficiency` is the
# combined detection probability for the smolt bypass system. Part III
# covers that. For adult window counts (this section), there is just r.
