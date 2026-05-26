# Section 3 - Variance, the delta method, and escapeLGD::nightFall()
# ---------------------------------------------------------------
# Goal: build up to the exact fallback joint likelihood used in
# escapeLGD/R/fallback_reascend_likelihood.R and the nighttime/fallback
# bootstrap loop in escapeLGD/R/night_fall_reascend_wc_binom.R.
#
# Repo pointers (escapeLGD/R/fallback_reascend_likelihood.R):
#   line 13:  fallback_log_likelihood <- function(pf, pre_f, dfr, df, dr, dt)
#   line 14:  return( dbinom(dfr, df, pre_f, log=TRUE) +
#                     dbinom(dr,  dt, pf*pre_f, log=TRUE) )
#
# Repo pointers (escapeLGD/R/night_fall_reascend_wc_binom.R):
#   lines 71-75:  optim(par = c(.1, .9), fn = optimllh,
#                       gr  = gradient_fallback_log_likelihood,
#                       dfr = ..., df = ..., dr = ..., dt = ...,
#                       control = list(fnscale = -1, maxit = 1000),
#                       method = "L-BFGS-B", upper = 1 - 1e-7, lower = 1e-7)
#   line 110:  rbinom(boots, totalPass[i], p_night[i]) / totalPass[i]

library(ggplot2)
library(dplyr)
library(purrr)
library(tibble)

set.seed(2026)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# --- 1. Nighttime correction: one source of noise ---------------------------
d_a       <- 200L     # PIT-tagged adults observed
p_n_true  <- 0.25     # truth: nighttime proportion
a_d_fixed <- 500L
nreps     <- 10000L

d_n_sim <- rbinom(nreps, size = d_a, prob = p_n_true)
p_n_hat <- d_n_sim / d_a
a_t_one <- a_d_fixed / (1 - p_n_hat)
cat("Layer 1 (only p_n noisy):  mean =", round(mean(a_t_one), 1),
    "  sd =", round(sd(a_t_one), 1), "\n\n")

# --- 2. Two sources of noise: a_d and p_n both estimated --------------------
a_d_true <- 500L; r_sample <- 5/6
w_sim   <- rbinom(nreps, size = a_d_true, prob = r_sample)
a_d_sim <- w_sim / r_sample
a_t_two <- a_d_sim / (1 - p_n_hat)
cat("Layer 2 (a_d and p_n both noisy):  mean =",
    round(mean(a_t_two), 1), "  sd =", round(sd(a_t_two), 1), "\n\n")

# --- 3. Delta-method variance approximation --------------------------------
# X = a_d_hat, Y = p_n_hat, g(X, Y) = X / (1 - Y).
# Independence assumed (data from disjoint subsamples).
var_a_d <- a_d_true * (1 - r_sample) / r_sample^2
var_p_n <- p_n_true * (1 - p_n_true) / d_a
g_dX <- 1 / (1 - p_n_true)
g_dY <- a_d_true / (1 - p_n_true)^2
var_a_t_delta <- g_dX^2 * var_a_d + g_dY^2 * var_p_n
cat("Delta-method:  Var(a_t) ~", round(var_a_t_delta, 2),
    "  sd(a_t)_delta =", round(sqrt(var_a_t_delta), 2),
    "  sd(a_t)_empir =", round(sd(a_t_two), 2), "\n\n")

# --- 4. The escapeLGD fallback joint likelihood (verbatim) -----------------
# Verbatim from escapeLGD/R/fallback_reascend_likelihood.R lines 13-15
fallback_log_likelihood <- function(pf, pre_f, dfr, df, dr, dt) {
  dbinom(dfr, df, pre_f,    log = TRUE) +
  dbinom(dr,  dt, pf*pre_f, log = TRUE)
}

# Verbatim wrapper that escapeLGD passes to optim() (line 61 of the file)
optimllh <- function(par, dfr, df, dr, dt) {
  pf    <- par[1]
  pre_f <- par[2]
  fallback_log_likelihood(pf, pre_f, dfr, df, dr, dt)
}

# Analytical gradient (lines 26-50 of fallback_reascend_likelihood.R)
gradient_fallback_log_likelihood <- function(par, dfr, df, dr, dt) {
  pf <- par[1]; pre_f <- par[2]
  c(
    (dr / pf) - (((dt - dr) * pre_f) / (1 - (pf * pre_f))),
    (dfr / pre_f) - ((df - dfr) / (1 - pre_f)) + (dr / pre_f) -
      (((dt - dr) * pf) / (1 - (pf * pre_f)))
  )
}

# --- 5. Simulate fallback data and run the joint MLE ----------------------
# Truth: 5% of ascending fish fall back. Of those, 60% later reascend.
pf_true    <- 0.05
pre_f_true <- 0.60
dt <- 2000L
dfr_count <- df_count <- dr_count <- 0L
df_count  <- rbinom(1, dt, pf_true)                # spillway-detected fallbacks
dfr_count <- rbinom(1, df_count, pre_f_true)        # of those, how many reascend
dr_count  <- rbinom(1, dt, pf_true * pre_f_true)    # ladder-observed reascensions

cat("Simulated counts: dt =", dt, "  df =", df_count,
    "  dfr =", dfr_count, "  dr =", dr_count, "\n")

# Verbatim optim call (escapeLGD lines 71-75)
opts <- optim(par = c(0.1, 0.9), fn = optimllh,
              gr  = gradient_fallback_log_likelihood,
              dfr = dfr_count, df = df_count,
              dr  = dr_count,  dt = dt,
              control = list(fnscale = -1, maxit = 1000),
              method = "L-BFGS-B",
              upper = 1 - 1e-7, lower = 1e-7)
cat("Joint MLE:  pf =", round(opts$par[1], 3),
    "  pre_f =", round(opts$par[2], 3),
    "  (truth pf =", pf_true,
    ", pre_f =", pre_f_true, ")\n\n")

# --- 6. Bootstrap nighttime passage rate the escapeLGD way -----------------
# escapeLGD/R/night_fall_reascend_wc_binom.R line 110 reproduced verbatim:
#   rbinom(boots, totalPass[i], p_night[i]) / totalPass[i]
boots <- 2000L
totalPass <- d_a
p_night_boot <- rbinom(boots, totalPass, p_n_true) / totalPass

cat("Bootstrap p_night rates (escapeLGD line 110):\n")
cat("  mean =", round(mean(p_night_boot), 4),
    "  sd =", round(sd(p_night_boot), 4), "\n\n")

# --- 7. Visualize the layered uncertainty ----------------------------------
skew_tbl <- bind_rows(
  tibble(layer = "1: only p_n noisy",        a_t = a_t_one),
  tibble(layer = "2: a_d and p_n noisy",     a_t = a_t_two),
  tibble(layer = "3: + fallback correction", a_t = a_t_two * (1 - pf_true))
)
p_skew <- ggplot(skew_tbl, aes(a_t, fill = layer)) +
  geom_histogram(bins = 60, alpha = 0.7) +
  facet_wrap(~ layer, ncol = 1, scales = "free_y") +
  geom_vline(xintercept = a_d_true / (1 - p_n_true),
             colour = "black", linetype = "dashed") +
  labs(title = "Section 3 - Layered uncertainty in a_t",
       subtitle = "Dashed line = true expected a_t before fallback",
       x = expression(hat(a)[t]),
       y = "Replicates") +
  theme(legend.position = "none")
ggsave(file.path(plots_dir, "section03_layered_uncertainty.png"), p_skew,
       width = 6, height = 7, dpi = 150)

# --- 8. End-of-section pointers --------------------------------------------
# You can now read:
#   escapeLGD/R/fallback_reascend_likelihood.R (full file, 71 lines)
#   escapeLGD/R/night_fall_reascend_wc_binom.R lines 22-116 (nightFall function)
# Compare line by line to the optim() call you just ran.
