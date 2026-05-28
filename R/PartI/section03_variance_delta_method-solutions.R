#---------------------------------------------------------
# File:   section03_variance_delta_method-solutions.R
# Part I, Section 3 solutions
#---------------------------------------------------------
# Section 3 solutions ----

#--------------------------------------
# Problem 3a: Layer 1 -- only p_n noisy
section3_problem_3a_fish <- function(d_a, p_n_true, a_d_fixed, nreps) {
  cat("\n----------------------------------\n")
  cat("Problem 3a: Layer 1 -- only p_n noisy\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  d_n_sim <- rbinom(nreps, size = d_a, prob = p_n_true)
  p_n_hat <- d_n_sim / d_a
  a_t     <- a_d_fixed / (1 - p_n_hat)
  cat("d_a =", d_a, "  p_n_true =", p_n_true,
      "  a_d_fixed =", a_d_fixed, "  nreps =", nreps, "\n")
  cat("mean(a_t) =", mean(a_t), "   sd(a_t) =", sd(a_t), "\n")
  invisible(a_t)
}


# Problem 3b: Layer 2 -- a_d and p_n both noisy
section3_problem_3b_fish <- function(a_d_true, r_sample, d_a, p_n_true, nreps) {
  cat("\n----------------------------------\n")
  cat("Problem 3b: Layer 2 -- a_d and p_n both noisy\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  w_sim   <- rbinom(nreps, size = a_d_true, prob = r_sample)
  d_n_sim <- rbinom(nreps, size = d_a, prob = p_n_true)
  a_d_sim <- w_sim / r_sample
  p_n_hat <- d_n_sim / d_a
  a_t     <- a_d_sim / (1 - p_n_hat)
  cat("a_d_true =", a_d_true, "  r_sample =", r_sample,
      "  d_a =", d_a, "  p_n_true =", p_n_true,
      "  nreps =", nreps, "\n")
  cat("mean(a_t) =", mean(a_t), "   sd(a_t) =", sd(a_t), "\n")
  invisible(a_t)
}


# Problem 3c: Delta-method variance approximation
section3_problem_3c_fish <- function(a_d_true, r_sample, p_n_true, d_a) {
  cat("\n----------------------------------\n")
  cat("Problem 3c: Delta-method variance for a_t = a_d / (1 - p_n)\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  var_a_d <- a_d_true * (1 - r_sample) / r_sample^2
  var_p_n <- p_n_true * (1 - p_n_true) / d_a
  g_dX <- 1 / (1 - p_n_true)
  g_dY <- a_d_true / (1 - p_n_true)^2
  var_a_t_delta <- g_dX^2 * var_a_d + g_dY^2 * var_p_n
  cat("Var(a_d)             =", var_a_d, "\n")
  cat("Var(p_n)             =", var_p_n, "\n")
  cat("Var(a_t) [delta]     =", var_a_t_delta, "\n")
  cat("sd(a_t)  [delta]     =", sqrt(var_a_t_delta), "\n")
  invisible(list(var_a_d = var_a_d, var_p_n = var_p_n,
                 var_a_t_delta = var_a_t_delta))
}


# Problem 3d: escapeLGD fallback joint MLE
section3_problem_3d_fish <- function(pf_true, pre_f_true, dt) {
  cat("\n----------------------------------\n")
  cat("Problem 3d: escapeLGD fallback joint MLE\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  # Simulate data
  df_count  <- rbinom(1, dt, pf_true)
  dfr_count <- rbinom(1, df_count, pre_f_true)
  dr_count  <- rbinom(1, dt, pf_true * pre_f_true)

  # escapeLGD fallback_reascend_likelihood.R lines 13-15 verbatim:
  fallback_log_likelihood <- function(pf, pre_f, dfr, df, dr, dt) {
    dbinom(dfr, df, pre_f,    log = TRUE) +
    dbinom(dr,  dt, pf*pre_f, log = TRUE)
  }
  optimllh <- function(par, dfr, df, dr, dt) {
    fallback_log_likelihood(par[1], par[2], dfr, df, dr, dt)
  }

  opts <- optim(par = c(0.1, 0.9), fn = optimllh,
                dfr = dfr_count, df = df_count,
                dr  = dr_count,  dt = dt,
                control = list(fnscale = -1, maxit = 1000),
                method = "L-BFGS-B",
                upper = 1 - 1e-7, lower = 1e-7)
  cat("Simulated: dt =", dt, "  df =", df_count,
      "  dfr =", dfr_count, "  dr =", dr_count, "\n")
  cat("Truth: pf =", pf_true, "  pre_f =", pre_f_true, "\n")
  cat("MLE:   pf =", round(opts$par[1], 3),
      "  pre_f =", round(opts$par[2], 3), "\n")
  invisible(opts)
}


# Problem 3e: escapeLGD nighttime-rate bootstrap
section3_problem_3e_fish <- function(d_a, p_n_true, boots) {
  cat("\n----------------------------------\n")
  cat("Problem 3e: escapeLGD nighttime-rate bootstrap\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  # escapeLGD night_fall_reascend_wc_binom.R line 110 verbatim:
  p_night_boot <- rbinom(boots, d_a, p_n_true) / d_a
  cat("d_a =", d_a, "  p_n_true =", p_n_true, "  boots =", boots, "\n")
  cat("mean(p_night_boot) =", mean(p_night_boot),
      "   sd =", sd(p_night_boot), "\n")
  invisible(p_night_boot)
}


# Narrative walkthrough ----
# The wrapper functions above ARE the correct answers.
# What follows is plain runnable code that explains the two-layer variance
# structure and the delta-method approximation, then shows where each breaks.

set.seed(2026)

# ---- shared parameters (match stub argument names) ----
d_a       <- 200L    # PIT-tagged adults available to estimate the nighttime rate
p_n_true  <- 0.25    # true fraction of total daily escapement that passes at night
a_d_fixed <- 500L    # daytime escapement (held constant in layer 1)
a_d_true  <- 500L    # true daytime escapement (free in layer 2)
r_sample  <- 5/6     # window sampling fraction
nreps     <- 10000L  # Monte Carlo replicates
pf_true   <- 0.05    # true P(fallback) — fraction of fish that drop back
pre_f_true <- 0.60   # true P(reascend | fallback)
dt        <- 2000L   # total ladder ascensions observed
boots     <- 2000L   # bootstrap draws for the nighttime-rate bootstrap

# ---- Problem 3a: layer 1 — only p_n is noisy ----
# The total escapement is a_t = a_d / (1 - p_n).
# In layer 1 we pretend a_d is known exactly; only p_n carries uncertainty.
# p_n is estimated from PIT-tagged adults: d_n ~ Binomial(d_a, p_n_true),
# then p_n_hat = d_n / d_a.

d_n_sim <- rbinom(nreps, size = d_a, prob = p_n_true)
# d_n_sim[i] = number of tagged fish detected moving at night on replicate i.
p_n_hat <- d_n_sim / d_a
# p_n_hat = MLE for p_n; same w/r pattern as a_d_hat but for a proportion.
a_t     <- a_d_fixed / (1 - p_n_hat)
# Dividing by (1 - p_n_hat) grosses up the daytime count to a total.
# When p_n_hat is near 1, (1 - p_n_hat) → 0 and a_t explodes — see extension.

cat("Problem 3a (layer 1: only p_n noisy)\n")
cat("  mean(a_t) =", round(mean(a_t), 1),
    "  sd(a_t) =", round(sd(a_t), 1), "\n")
# mean ≈ a_d_fixed/(1-p_n_true) ≈ 667.  The estimator is slightly biased
# upward because E[1/(1-X)] > 1/(1-E[X]) (Jensen's inequality).

# ---- Problem 3b: layer 2 — a_d and p_n both noisy ----
# Now add the window-count noise on top of the p_n noise.
# Both w and d_n are random, so a_t = (w/r) / (1 - d_n/d_a) has two sources
# of variance that compound each other.

w_sim   <- rbinom(nreps, size = a_d_true, prob = r_sample)   # window counts
d_n_sim <- rbinom(nreps, size = d_a,      prob = p_n_true)   # night detections
a_d_sim <- w_sim / r_sample      # daytime MLE, one per replicate
p_n_hat <- d_n_sim / d_a         # nighttime-rate MLE, one per replicate
a_t     <- a_d_sim / (1 - p_n_hat)   # total escapement estimate, one per replicate

cat("Problem 3b (layer 2: both noisy)\n")
cat("  mean(a_t) =", round(mean(a_t), 1),
    "  sd(a_t) =", round(sd(a_t), 1), "\n")
# sd is larger than layer 1 because window-count noise adds on top.
# The ratio sd(layer2) / sd(layer1) tells you how much the window contributes
# relative to the nighttime-rate sensor.

# ---- Problem 3c: delta-method variance ----
# The delta method approximates Var(g(X, Y)) when g is a smooth function of
# two independent random variables X and Y:
#   Var(g(X,Y)) ≈ (∂g/∂X)² Var(X) + (∂g/∂Y)² Var(Y)
#
# Here g(X, Y) = X / (1 - Y)  where X = a_d, Y = p_n.
# Partial derivatives:
#   ∂g/∂X = 1 / (1 - p_n)          = g_dX
#   ∂g/∂Y = X / (1 - p_n)²  = a_d / (1 - p_n)²  = g_dY

var_a_d <- a_d_true * (1 - r_sample) / r_sample^2
# Binomial variance for w is a_d * r * (1-r), so Var(w/r) = a_d*(1-r)/r^2.

var_p_n <- p_n_true * (1 - p_n_true) / d_a
# Binomial variance for d_n/d_a: p*(1-p)/n.

g_dX <- 1 / (1 - p_n_true)
g_dY <- a_d_true / (1 - p_n_true)^2

var_a_t_delta <- g_dX^2 * var_a_d + g_dY^2 * var_p_n
# Each term = (how sensitive a_t is to one source) × (how variable that source is).
# The larger term tells you which sensor to upgrade first.

cat("Problem 3c (delta method)\n")
cat("  Var(a_d) =", round(var_a_d, 1), "\n")
cat("  Var(p_n) =", round(var_p_n, 5), "\n")
cat("  Var(a_t) delta =", round(var_a_t_delta, 1), "\n")
cat("  sd(a_t)  delta =", round(sqrt(var_a_t_delta), 1),
    "  vs simulated =", round(sd(a_t), 1), "\n")
# Delta-method sd should be close to the simulated sd from 3b.
# The approximation degrades when p_n is large (denominator near zero).

# ---- Problem 3d: escapeLGD fallback joint MLE ----
# Some adults ascend the ladder, drop back (fallback), then reascend.
# escapeLGD models this with two Binomials that share parameters:
#   df  ~ Binomial(dt, pf)          — fish that fall back at all
#   dfr ~ Binomial(df, pre_f)       — of those, the fraction that reascend
#   dr  ~ Binomial(dt, pf * pre_f)  — directly observed reascensions
# The joint log-likelihood from fallback_reascend_likelihood.R lines 13-15
# is the sum of two independent Binomial log-likelihoods.

df_count  <- rbinom(1, dt, pf_true)
# Simulate: out of dt total ladder crossings, df_count fell back.
dfr_count <- rbinom(1, df_count, pre_f_true)
# Of those fallbacks, dfr_count reascended.
dr_count  <- rbinom(1, dt, pf_true * pre_f_true)
# dr is the same event observed from the downstream detector.

fallback_log_likelihood <- function(pf, pre_f, dfr, df, dr, dt) {
  dbinom(dfr, df, pre_f,    log = TRUE) +
  dbinom(dr,  dt, pf*pre_f, log = TRUE)
}
optimllh <- function(par, dfr, df, dr, dt) {
  fallback_log_likelihood(par[1], par[2], dfr, df, dr, dt)
}
# par[1] = pf, par[2] = pre_f.  We maximise so fnscale = -1.
opts <- optim(par = c(0.1, 0.9), fn = optimllh,
              dfr = dfr_count, df = df_count,
              dr  = dr_count,  dt = dt,
              control = list(fnscale = -1, maxit = 1000),
              method = "L-BFGS-B",
              upper = 1 - 1e-7, lower = 1e-7)
# L-BFGS-B handles box constraints (parameters must stay in (0,1)).
# lower/upper bound away from 0 and 1 to avoid log(0) in dbinom.

cat("Problem 3d (fallback joint MLE)\n")
cat("  Truth: pf =", pf_true, "  pre_f =", pre_f_true, "\n")
cat("  MLE:   pf =", round(opts$par[1], 3),
    "  pre_f =", round(opts$par[2], 3), "\n")

# ---- Problem 3e: escapeLGD nighttime-rate bootstrap ----
# escapeLGD propagates uncertainty in p_n by parametric bootstrap:
# draw bootstrap p_n values from Binomial(d_a, p_n_true) / d_a,
# then pass each draw through the a_t formula.
# This is vectorised in escapeLGD line 110.

p_night_boot <- rbinom(boots, d_a, p_n_true) / d_a
# Each element is one bootstrap estimate of the nighttime rate.
# The spread of this vector directly determines the width of the CI on a_t.

cat("Problem 3e (nighttime-rate bootstrap)\n")
cat("  mean(p_night_boot) =", round(mean(p_night_boot), 4),
    "  sd =", round(sd(p_night_boot), 4), "\n")
# Compare sd to sqrt(p_n_true*(1-p_n_true)/d_a) — the theoretical SE of p_n_hat.
cat("  theoretical SE(p_n_hat) =",
    round(sqrt(p_n_true * (1 - p_n_true) / d_a), 4), "\n")

# ---- Extension: delta method breaks near p_n = 1 ----
# The delta method assumes the function g(X,Y) = X/(1-Y) is well approximated
# by its first-order Taylor expansion near the true parameter values.
# When p_n is large the denominator (1 - p_n) is small and the function is
# highly curved — the linear approximation fails badly.

p_n_vals <- c(0.10, 0.25, 0.50, 0.70, 0.85)
for (pn in p_n_vals) {
  vd <- (1 / (1 - pn))^2 * (a_d_true * (1 - r_sample) / r_sample^2) +
        (a_d_true / (1 - pn)^2)^2 * (pn * (1 - pn) / d_a)
  # Run a quick simulation to get the true sd
  d_n_s <- rbinom(5000L, d_a, pn)
  w_s   <- rbinom(5000L, a_d_true, r_sample)
  a_t_s <- (w_s / r_sample) / (1 - d_n_s / d_a)
  a_t_s <- a_t_s[is.finite(a_t_s)]
  cat("  p_n =", pn, ": delta-sd =", round(sqrt(vd), 0),
      "  sim-sd =", round(sd(a_t_s), 0), "\n")
}
# At p_n = 0.85 the delta sd underestimates the simulated sd substantially.
# escapeLGD avoids this by using the parametric bootstrap (3e) rather than the
# delta method for final inference.

# ---- Forward pointer ----
# Section 4 builds the full parametric bootstrap for the smolt-trap passage
# count, reproducing SCRAPI's bootsmolt daily-count loop (lines 139-145) and
# the vectorised escapeLGD equivalent.
