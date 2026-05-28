#---------------------------------------------------------
# File:   section02_mle_window_count-solutions.R
# Part I, Section 2 solutions
#---------------------------------------------------------
# Section 2 solutions ----

library(ggplot2)
library(tibble)


#--------------------------------------
# Problem 2a: Replicated MLE centering
section2_problem_2a_fish <- function(a_d_true, r_sample, nreps) {
  cat("\n----------------------------------\n")
  cat("Problem 2a: Replicated MLE centering check\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  estimates <- vapply(seq_len(nreps),
                      function(i) rbinom(1, a_d_true, r_sample) / r_sample,
                      numeric(1))
  cat("a_d_true =", a_d_true, "  r_sample =", r_sample,
      "  nreps =", nreps, "\n")
  cat("mean(a_d_hat) =", mean(estimates),
      "   sd(a_d_hat) =", sd(estimates), "\n")
  invisible(estimates)
}


# Problem 2b: Grid evaluation of the log-likelihood
section2_problem_2b_fish <- function(w_obs, r_sample, a_d_max) {
  cat("\n----------------------------------\n")
  cat("Problem 2b: Grid evaluation of the binomial log-likelihood\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  a_d_grid <- seq(max(w_obs, 1L), a_d_max, by = 1L)
  loglik   <- dbinom(w_obs, size = a_d_grid, prob = r_sample, log = TRUE)
  mle_grid <- a_d_grid[which.max(loglik)]
  cat("w_obs =", w_obs, "  r_sample =", r_sample, "\n")
  cat("Grid MLE         =", mle_grid, "\n")
  cat("Closed-form w/r  =", w_obs / r_sample, "\n")

  p <- ggplot(tibble(a_d = a_d_grid, loglik = loglik),
              aes(a_d, loglik)) +
    geom_line(colour = "steelblue", linewidth = 1) +
    geom_vline(xintercept = w_obs / r_sample,
               colour = "firebrick", linewidth = 1, linetype = "dashed") +
    labs(title = "Section 2 - Binomial log-likelihood for a_d",
         x = expression(a[d]), y = expression(log~L(a[d])))
  invisible(list(a_d_grid = a_d_grid, loglik = loglik, mle = mle_grid))
}


# Problem 2c: optim() on a continuous relaxation
section2_problem_2c_fish <- function(w_obs, r_sample) {
  cat("\n----------------------------------\n")
  cat("Problem 2c: optim() continuous-relaxation MLE\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  neg_loglik_relaxed <- function(a_d, w, r) {
    if (a_d < w) return(Inf)
    -(lgamma(a_d + 1) - lgamma(w + 1) - lgamma(a_d - w + 1) +
      w * log(r) + (a_d - w) * log(1 - r))
  }
  fit <- optim(par = w_obs * 1.1, fn = neg_loglik_relaxed,
               w = w_obs, r = r_sample,
               method = "Brent", lower = w_obs, upper = 5000)
  cat("w_obs =", w_obs, "  r_sample =", r_sample, "\n")
  cat("optim MLE        =", round(fit$par, 2), "\n")
  cat("Closed-form w/r  =", w_obs / r_sample, "\n")
  invisible(fit)
}


# Problem 2d: Inverse-SR weighting in SCRAPI's thetahat
section2_problem_2d_fish <- function(stocks, strats, SR, n_fish) {
  cat("\n----------------------------------\n")
  cat("Problem 2d: Inverse-SR weighting (SCRAPI thetahat pattern, SMOLT trap)\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  AllPrime <- tibble(
    Strat = sample(strats, size = n_fish, replace = TRUE),
    PGrp  = sample(stocks, size = n_fish, replace = TRUE),
    SR    = SR
  )
  # SCRAPI line 94 reproduced:
  Primarystrata <- tapply(1 / AllPrime$SR,
                          list(factor(AllPrime$Strat, levels = strats),
                               factor(AllPrime$PGrp,  levels = stocks)),
                          sum)
  Primarystrata[is.na(Primarystrata)] <- 0
  Primaryproportions <- prop.table(Primarystrata, margin = 1)

  cat("AllPrime rows =", n_fish, "  SR =", round(SR, 3), "\n")
  cat("Primarystrata (inverse-SR weighted counts):\n")
  print(round(Primarystrata, 1))
  cat("Primaryproportions (row-normalized within stratum):\n")
  print(round(Primaryproportions, 3))
  invisible(list(AllPrime = AllPrime,
                 Primarystrata = Primarystrata,
                 Primaryproportions = Primaryproportions))
}


# Narrative walkthrough ----
# The wrapper functions above ARE the correct answers.
# What follows is plain runnable code that explains why each answer works
# and shows where the MLE machinery breaks in real data.

set.seed(2026)

# ---- shared parameters (match stub argument names) ----
a_d_true <- 500L              # true daytime adult escapement
r_sample <- 5/6               # fraction of ladder width covered by the window
nreps    <- 1000L             # Monte Carlo repetitions for the centering check
a_d_max  <- 700L              # upper bound for the log-likelihood grid
stocks   <- c("LOSALM", "CHMBLN", "IMNAHA")   # stock identifiers (PGrp values)
strats   <- c("S1", "S2", "S3")               # temporal strata
SR       <- 5/6 * 0.45        # smolt-trap combined detection probability
n_fish   <- 90L               # fish sampled at the bypass trap

# ---- Problem 2a: the MLE is w / r — prove it by simulation ----
# Each day a_d_true fish pass through the ladder.  The window counts only the
# fraction r_sample of the ladder width, so the window count w is:
#   w ~ Binomial(n = a_d_true, p = r_sample)
#
# The log-likelihood for a_d given observed w is:
#   log L(a_d | w) = log C(a_d, w) + w log(r) + (a_d - w) log(1 - r)
# Differentiating with respect to a_d and setting to zero gives a_d_hat = w / r.
# Here we verify that formula centres on the truth across many replicates.

estimates <- vapply(seq_len(nreps),
                    function(i) rbinom(1, a_d_true, r_sample) / r_sample,
                    numeric(1))
# rbinom(1, a_d_true, r_sample) draws one window count w.
# Dividing by r_sample is the MLE for that draw.
# vapply repeats this nreps times and collects results into a numeric vector.

cat("Problem 2a\n")
cat("  mean(a_d_hat) =", round(mean(estimates), 1),
    "  sd =", round(sd(estimates), 1), "\n")
# mean ≈ 500 confirms the estimator is unbiased.
# sd quantifies the shot-to-shot noise at this sampling rate.

# ---- Problem 2b: see the log-likelihood surface directly ----
# Rather than trusting the algebra, we evaluate log L at every integer a_d
# from w_obs up to a_d_max.  The peak should land exactly at w / r.

w_obs <- rbinom(1, a_d_true, r_sample)   # one window count for this example

# The grid must start at w_obs: if a_d < w_obs the window count would exceed
# the number of fish that passed, which has probability zero.
a_d_grid <- seq(max(w_obs, 1L), a_d_max, by = 1L)

# dbinom evaluates the Binomial probability mass function at each grid point.
# log = TRUE gives log-probabilities, which are numerically stable for small p.
loglik   <- dbinom(w_obs, size = a_d_grid, prob = r_sample, log = TRUE)

mle_grid <- a_d_grid[which.max(loglik)]   # integer with highest log-likelihood

cat("Problem 2b\n")
cat("  w_obs =", w_obs, "\n")
cat("  Grid MLE    =", mle_grid,         "\n")   # nearest integer to w/r
cat("  w / r       =", w_obs / r_sample, "\n")   # analytic MLE
# They agree to within rounding — the surface peaks exactly where the
# algebra predicts.

# ---- Problem 2c: optim() when closed forms do not exist ----
# dbinom requires integer size, so we cannot differentiate through it with
# standard optimisers.  We relax by replacing the factorial ratio with lgamma,
# which extends log(n!) smoothly to all real n >= 0.
#   lgamma(n + 1) == log(n!)  for integer n
#   lgamma(x + 1) is smooth and differentiable for real x > 0

neg_loglik_relaxed <- function(a_d, w, r) {
  if (a_d < w) return(Inf)   # constraint: can't observe more fish than passed
  # Relaxed log C(a_d, w) = lgamma(a_d+1) - lgamma(w+1) - lgamma(a_d-w+1)
  -(lgamma(a_d + 1) - lgamma(w + 1) - lgamma(a_d - w + 1) +
    w * log(r) + (a_d - w) * log(1 - r))
}
# The negative sign turns maximisation into the minimisation that optim() does.

# method = "Brent" is for 1-D problems: efficient and exact on an interval.
fit <- optim(par = w_obs * 1.1, fn = neg_loglik_relaxed,
             w = w_obs, r = r_sample,
             method = "Brent", lower = w_obs, upper = 5000)

cat("Problem 2c\n")
cat("  optim MLE   =", round(fit$par, 2), "\n")
cat("  w / r       =", w_obs / r_sample,  "\n")
# The continuous optimiser recovers the analytic answer.
# This pattern — lgamma relaxation + Brent — is what escapeLGD uses when
# the likelihood has no closed-form derivative.

# ---- Problem 2d: inverse-SR weighting (SCRAPI thetahat pattern) ----
# At the smolt bypass, SR = SampleRate * GuidanceEfficiency is the probability
# that a passing fish is both guided into the trap and actually counted.
# A fish caught with SR = 0.20 represents 1/0.20 = 5 fish in the population;
# a fish caught with SR = 0.40 represents only 2.5 fish.
# Summing 1/SR within each Strat x PGrp cell is the Horvitz-Thompson estimator
# for the number of fish in that cell.

AllPrime <- tibble(
  Strat = sample(strats, size = n_fish, replace = TRUE),   # temporal stratum
  PGrp  = sample(stocks, size = n_fish, replace = TRUE),   # stock assignment
  SR    = SR                                               # constant here
)
# One row per sampled fish — mirrors SCRAPI's AllPrime data frame exactly.

# SCRAPI line 94: tapply sums 1/SR within every Strat x PGrp combination.
Primarystrata <- tapply(1 / AllPrime$SR,
                        list(factor(AllPrime$Strat, levels = strats),
                             factor(AllPrime$PGrp,  levels = stocks)),
                        sum)
Primarystrata[is.na(Primarystrata)] <- 0   # empty cells get weight zero

# prop.table with margin = 1 divides each row by its row total,
# converting weighted counts into within-stratum stock proportions.
Primaryproportions <- prop.table(Primarystrata, margin = 1)

cat("Problem 2d\n")
cat("  Primaryproportions (each row sums to 1.0):\n")
print(round(Primaryproportions, 3))
# When SR is constant across strata, the 1/SR factors cancel in the ratio
# and proportions equal the raw fish-count proportions.
# The weighting only changes results when SR varies across strata or time.

# ---- Extension: what breaks when SR varies across strata ----
# In real trap data, SR often changes across the season as operators calibrate
# the system or flows shift.  Ignoring that variation biases stock composition.

SR_by_strat <- c(S1 = 0.15, S2 = 0.35, S3 = 0.50)   # early season has low SR

set.seed(42)
AllPrime_var <- tibble(
  Strat = sample(strats, size = n_fish, replace = TRUE),
  PGrp  = sample(stocks, size = n_fish, replace = TRUE),
  SR    = SR_by_strat[sample(strats, size = n_fish, replace = TRUE)]
)

# Naïve: ignore SR entirely, just count fish per cell.
raw_props <- prop.table(table(AllPrime_var$Strat, AllPrime_var$PGrp), margin = 1)

# Correct: weight each fish by 1/SR before summing.
wt_counts <- tapply(1 / AllPrime_var$SR,
                    list(factor(AllPrime_var$Strat, levels = strats),
                         factor(AllPrime_var$PGrp,  levels = stocks)), sum)
wt_counts[is.na(wt_counts)] <- 0
wt_props  <- prop.table(wt_counts, margin = 1)

cat("\nExtension: variable SR across strata\n")
cat("  Naïve (raw count) proportions for S1:\n")
print(round(raw_props["S1", ], 3))
cat("  Inverse-SR weighted proportions for S1:\n")
print(round(wt_props["S1", ], 3))
# S1 has the lowest SR so fish there are hardest to catch.
# Naïve counts underrepresent S1 fish; inverse-SR weighting corrects this.
# The bias is proportional to how much SR varies — in real LGD data the
# difference can shift stock composition estimates by several percent.

# ---- Forward pointer ----
# Section 3 combines the a_d estimator (Section 1-2) with the nighttime-
# proportion estimator p_n_hat to form the total-escapement estimator
#   a_t = a_d / (1 - p_n)
# and derives its variance via the delta method.
