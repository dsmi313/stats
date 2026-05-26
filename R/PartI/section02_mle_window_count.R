# Section 2 - MLE for the window-count estimator
# ---------------------------------------------------------------
# MIT 18.05 Class 10b (Maximum Likelihood Estimates) + R Studio 2.
# We treat the daytime escapement a_d as the parameter to estimate
# from a single observed window count w. The Binomial MLE
# a_d_hat = w / r is verified three ways:
#   (a) replicate() simulation centering on the truth,
#   (b) grid evaluation of the log-likelihood,
#   (c) optim() on a continuous relaxation via lgamma.
#
# Reading: PLAN.md Section 2, MIT Class 10a/10b, R Studio 2.
# StatQuest: "In Statistics, Probability is not Likelihood" - the
# probability-vs-likelihood plot below illustrates that distinction.

library(ggplot2)
library(dplyr)
library(purrr)

set.seed(2026)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

a_d_true <- 500L
r_sample <- 5/6

# --- 1. Repeated estimates: a_d_hat centers on the truth -------------------
nreps <- 1000L
estimates <- map_dbl(seq_len(nreps),
                     ~ rbinom(1, a_d_true, r_sample) / r_sample)

cat("Replicated MLE check:\n")
cat("  mean(a_d_hat) =", round(mean(estimates), 1),
    "  (target", a_d_true, ")\n")
cat("  sd(a_d_hat)   =", round(sd(estimates), 1), "\n\n")

# --- 2. Log-likelihood for a single observed w -----------------------------
# Given one observed w, view L(a_d) = P(W = w | a_d, r) as a function of a_d.
w_obs    <- rbinom(1, a_d_true, r_sample)
a_d_grid <- seq(max(w_obs, 1), 700, by = 1L)
loglik   <- dbinom(w_obs, size = a_d_grid, prob = r_sample, log = TRUE)
mle_grid <- a_d_grid[which.max(loglik)]

cat("Observed window count w =", w_obs, "\n")
cat("  closed-form MLE w/r  =", w_obs / r_sample, "\n")
cat("  grid MLE             =", mle_grid, "\n\n")

p_lik <- ggplot(tibble(a_d = a_d_grid, loglik = loglik),
                aes(a_d, loglik)) +
  geom_line(colour = "steelblue", linewidth = 1) +
  geom_vline(xintercept = w_obs / r_sample,
             colour = "firebrick", linewidth = 1, linetype = "dashed") +
  labs(title = "Section 2 - Binomial log-likelihood for a_d",
       subtitle = paste0("w = ", w_obs, ", r = 5/6; peak at w/r"),
       x = expression(a[d]),
       y = expression(log~L(a[d])))
ggsave(file.path(plots_dir, "section02_loglik.png"), p_lik,
       width = 6, height = 4, dpi = 150)

# --- 3. optim() on a continuous relaxation ---------------------------------
# dbinom requires integer size; use lgamma to relax over real-valued a_d.
neg_loglik_relaxed <- function(a_d, w, r) {
  if (a_d < w) return(Inf)
  -(lgamma(a_d + 1) - lgamma(w + 1) - lgamma(a_d - w + 1) +
    w * log(r) + (a_d - w) * log(1 - r))
}

fit <- optim(par = w_obs * 1.1,
             fn  = neg_loglik_relaxed,
             w = w_obs, r = r_sample,
             method = "Brent", lower = w_obs, upper = 5000)
cat("optim continuous-relaxation MLE =", round(fit$par, 2),
    "  (target w/r =", w_obs / r_sample, ")\n\n")

# --- 4. Probability is not likelihood --------------------------------------
# Probability view: fix a_d, sweep observable w (the pmf).
# Likelihood view: fix w, sweep parameter a_d.
pmf_tbl <- tibble(w = 380:460,
                  density = dbinom(380:460, size = a_d_true, prob = r_sample),
                  view = "Probability: P(W | a_d = 500)")
lik_tbl <- tibble(a_d = a_d_grid,
                  density = exp(loglik - max(loglik)),
                  view = "Likelihood: L(a_d | w obs)")

# Each panel has its own x-axis (w vs a_d) so the two views stay separate.
panel_tbl <- bind_rows(
  pmf_tbl |> rename(x = w),
  lik_tbl |> rename(x = a_d)
)

p_pl <- ggplot(panel_tbl, aes(x, density)) +
  geom_line(colour = "steelblue", linewidth = 1) +
  facet_wrap(~ view, scales = "free", ncol = 1) +
  labs(title = "Section 2 - Probability is not Likelihood",
       subtitle = "Same dbinom() viewed two ways",
       x = NULL, y = "density / relative likelihood")
ggsave(file.path(plots_dir, "section02_prob_vs_lik.png"), p_pl,
       width = 6, height = 6, dpi = 150)

# --- 5. Bias and SE summary ------------------------------------------------
true_var <- a_d_true * (1 - r_sample) / r_sample^2
cat("Sanity check on Var(a_d_hat):\n")
cat("  closed form Var(W/r)  =", round(true_var, 2),
    "  -> sd =", round(sqrt(true_var), 2), "\n")
cat("  empirical sd          =", round(sd(estimates), 2), "\n")

# Section 2 payoff ----------------------------------------------------------
# The same likelihood machinery scales to the nighttime PIT proportion,
# trap detection probability, and PBT tag rates in Sections 3 and 7-10.
