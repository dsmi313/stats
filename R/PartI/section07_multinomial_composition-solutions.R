#---------------------------------------------------------
# File:   section07_multinomial_composition-solutions.R
# Part I, Section 7 solutions
#---------------------------------------------------------
library(dplyr)

softmax <- function(x) { e <- exp(x - max(x)); e / sum(e) }

# Problem 7a
section7_problem_7a_fish <- function(n = 200, props = c(0.6, 0.3, 0.1)) {
  cat("\n----------------------------------\n")
  cat("Problem 7a: rmultinom — MLE is observed proportions\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  obs <- rmultinom(1, size = n, prob = props)[, 1]
  cat("True props:  ", round(props, 3), "\n")
  cat("MLE (sample):", round(obs / sum(obs), 3), "\n")

  reps <- replicate(1000, {
    x <- rmultinom(1, size = n, prob = props)[, 1]
    x[1] / sum(x)
  })
  hist(reps, breaks = 30, main = "Stock 1 proportion MLE (1,000 sims)",
       xlab = "Estimated proportion", col = "steelblue")
  abline(v = props[1], col = "red", lwd = 2)
  invisible(obs / sum(obs))
}


# Problem 7b
section7_problem_7b_fish <- function(counts = c(120L, 60L, 20L)) {
  cat("\n----------------------------------\n")
  cat("Problem 7b: multinomial log-likelihood; maximize with optim()\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  ll <- function(par, counts) sum(counts * log(softmax(par)))
  K   <- length(counts)
  fit <- optim(rep(0, K), ll, counts = counts,
               control = list(fnscale = -1), method = "BFGS")
  mle_opt  <- softmax(fit$par)
  mle_true <- counts / sum(counts)
  cat("optim MLE:    ", round(mle_opt, 4), "\n")
  cat("counts/n MLE: ", round(mle_true, 4), "\n")
  invisible(mle_opt)
}


# Problem 7c
section7_problem_7c_fish <- function(
    tagged_counts   = c(S1 = 42L, S2 = 35L, S3 = 18L),
    untagged_counts = c(S1 = 80L, S2 = 60L, S3 = 30L, Unassigned = 32L)) {

  cat("\n----------------------------------\n")
  cat("Problem 7c: product of two multinomials — skeleton of Delomas & Hess\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  K <- length(tagged_counts)
  joint_ll <- function(par, tc, uc) {
    p_tag  <- softmax(par[1:K])
    p_full <- softmax(c(par[1:K], par[K + 1]))
    sum(tc * log(p_tag)) + sum(uc * log(p_full))
  }
  fit <- optim(rep(0, K + 1), joint_ll,
               tc = tagged_counts, uc = untagged_counts,
               control = list(fnscale = -1), method = "BFGS")
  p_stocks <- softmax(fit$par[1:K])
  p_full   <- softmax(c(fit$par[1:K], fit$par[K + 1]))
  cat("Stock proportions:\n")
  print(round(setNames(p_stocks, names(tagged_counts)), 4))
  cat("Full composition (incl. Unassigned):\n")
  print(round(setNames(p_full, c(names(tagged_counts), "Unassigned")), 4))
  invisible(list(p_stocks = p_stocks, p_full = p_full))
}


# Problem 7d
section7_problem_7d_fish <- function(counts = c(S1 = 120L, S2 = 60L, S3 = 20L),
                                      boots  = 10000L) {
  cat("\n----------------------------------\n")
  cat("Problem 7d: nonparametric bootstrap CI for stock proportions\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  fish   <- rep(names(counts), counts)
  stocks <- names(counts)
  boot_mat <- matrix(NA_real_, nrow = boots, ncol = length(counts))
  colnames(boot_mat) <- stocks
  for (b in seq_len(boots)) {
    smp <- sample(fish, length(fish), replace = TRUE)
    tbl <- table(factor(smp, levels = stocks))
    boot_mat[b, ] <- as.numeric(tbl) / sum(tbl)
  }
  cis <- apply(boot_mat, 2, quantile, c(0.025, 0.975))
  cat("Bootstrap 95% CIs:\n"); print(round(cis, 3))
  cat("Observed proportions:\n"); print(round(counts / sum(counts), 3))
  invisible(cis)
}


# Narrative walkthrough ----
# The wrapper functions above ARE the correct answers.
# What follows is plain runnable code that builds the multinomial from first
# principles, shows why the MLE is always the observed proportion, traces the
# softmax reparametrisation trick, and derives the joint-likelihood skeleton
# used by the Delomas & Hess PBT estimator.

set.seed(2026)

# ---- shared parameters (match stub argument names and test values) ----
n              <- 200L          # fish in one sampling draw
props          <- c(0.6, 0.3, 0.1)   # true stock proportions
counts         <- c(120L, 60L, 20L)  # observed stock counts (sum = 200)
tagged_counts  <- c(S1 = 42L, S2 = 35L, S3 = 18L)   # PBT-resolved fish
untagged_counts <- c(S1 = 80L, S2 = 60L, S3 = 30L, Unassigned = 32L)
boots          <- 5000L         # nonparametric bootstrap iterations

# ---- Problem 7a: rmultinom — MLE is observed proportions ----
# The multinomial generalises the binomial from 2 to K outcomes.
# If n fish are independently drawn from a population where each stock has
# probability props[k], then the count vector c follows:
#   c ~ Multinomial(n, props)
# The log-likelihood is:
#   log L(props | c) = const + sum_k c[k] * log(props[k])
# Maximising subject to sum(props) = 1 gives props_hat = c / n.

obs <- rmultinom(1, size = n, prob = props)[, 1]
# rmultinom returns a K × 1 matrix; [,1] extracts the vector.
cat("Problem 7a\n")
cat("  True props:   ", round(props, 3), "\n")
cat("  MLE (obs/n):  ", round(obs / sum(obs), 3), "\n")

# Replicate 1000 times to see the sampling distribution of the stock-1 estimate.
reps <- replicate(1000, {
  x <- rmultinom(1, size = n, prob = props)[, 1]
  x[1] / sum(x)   # stock-1 MLE on this replicate
})
cat("  stock-1 MLE mean =", round(mean(reps), 3),
    "  sd =", round(sd(reps), 4), "\n")
# mean ≈ props[1] = 0.60 (unbiased); sd decreases as n grows.
# The theoretical SE is sqrt(p*(1-p)/n) = sqrt(0.6*0.4/200) ≈ 0.035.
cat("  theoretical SE =", round(sqrt(props[1]*(1-props[1])/n), 4), "\n")

# ---- Problem 7b: multinomial log-likelihood and optim() ----
# Writing the log-likelihood directly and maximising numerically recovers
# the closed-form answer — useful when the model gets more complex (7c).
#
# Constraint: props must be non-negative and sum to 1.
# Trick: reparametrise with softmax so optim() can work unconstrained.
#   softmax(u)[k] = exp(u[k]) / sum(exp(u))
# Any real vector u maps to a valid probability vector, so optim() is free
# to search all of R^K without hitting boundary constraints.

ll <- function(par, counts) sum(counts * log(softmax(par)))
# par is an unconstrained vector; softmax(par) gives the props.
# We maximise, so fnscale = -1 (optim() minimises by default).

K   <- length(counts)
fit <- optim(rep(0, K), ll, counts = counts,
             control = list(fnscale = -1), method = "BFGS")
mle_opt  <- softmax(fit$par)     # convert back from unconstrained to props
mle_true <- counts / sum(counts) # closed-form MLE
cat("Problem 7b\n")
cat("  optim MLE:    ", round(mle_opt, 4), "\n")
cat("  counts/n MLE: ", round(mle_true, 4), "\n")
# They agree — optim recovers the closed form through numerical search.

# ---- Problem 7c: joint log-likelihood for tagged and untagged samples ----
# In PBT, two independent samples share the same underlying stock proportions:
#   Tagged fish:   resolved to a stock → Multinomial(n_tag, p_tag)
#   Untagged fish: may be unassigned   → Multinomial(n_unt, p_full)
#                  where p_full = (p_tag, p_unassigned)
#
# Because the samples are independent, the joint log-likelihood is the sum.
# This is the structural skeleton of composition_estimation_utils.R lines 62-66.

K_tag  <- length(tagged_counts)   # number of fully resolved stocks (3)
# K_tag parameters for tagged proportions + 1 for the unassigned category = K_tag+1.

joint_ll <- function(par, tc, uc) {
  p_tag  <- softmax(par[1:K_tag])               # proportions among tagged fish
  p_full <- softmax(c(par[1:K_tag], par[K_tag + 1]))  # adds unassigned bucket
  sum(tc * log(p_tag)) + sum(uc * log(p_full))  # joint = sum of two log-likelihoods
}
fit_jt <- optim(rep(0, K_tag + 1), joint_ll,
                tc = tagged_counts, uc = untagged_counts,
                control = list(fnscale = -1), method = "BFGS")
p_stocks <- softmax(fit_jt$par[1:K_tag])
p_full   <- softmax(c(fit_jt$par[1:K_tag], fit_jt$par[K_tag + 1]))
cat("Problem 7c\n")
cat("  Stock proportions:\n")
print(round(setNames(p_stocks, names(tagged_counts)), 4))
cat("  Full composition (incl. Unassigned):\n")
print(round(setNames(p_full, c(names(tagged_counts), "Unassigned")), 4))
# The joint model borrows information from BOTH samples to estimate stock props.
# The unassigned fish inform the total but are not assigned to a specific stock.

# ---- Problem 7d: nonparametric bootstrap CI for stock proportions ----
# Instead of deriving a variance formula, we resample the individual fish
# with replacement and compute the MLE on each resample.
# This automatically captures the discrete, bounded nature of proportions
# without relying on normal approximations.

named_counts <- c(S1 = 120L, S2 = 60L, S3 = 20L)   # use named version from test
fish    <- rep(names(named_counts), named_counts)    # expand to individual fish
stock_names <- names(named_counts)
boot_mat <- matrix(NA_real_, nrow = boots, ncol = length(named_counts))
colnames(boot_mat) <- stock_names
for (b in seq_len(boots)) {
  smp <- sample(fish, length(fish), replace = TRUE)
  # Resample n fish with replacement — the bootstrap principle.
  tbl <- table(factor(smp, levels = stock_names))
  # factor() ensures all stocks appear even if none are drawn (shouldn't happen here).
  boot_mat[b, ] <- as.numeric(tbl) / sum(tbl)   # MLE on this bootstrap sample
}
cis <- apply(boot_mat, 2, quantile, c(0.025, 0.975))
cat("Problem 7d (boots =", boots, ")\n")
cat("  Bootstrap 95% CIs:\n"); print(round(cis, 3))
cat("  Observed proportions:\n")
print(round(named_counts / sum(named_counts), 3))

# ---- Extension: what breaks in the joint likelihood when sample sizes are small ----
# The joint MLE in 7c works well with large samples but can fail in two ways
# when one sample is small.

# Failure 1: when the tagged sample is tiny, optim() can get stuck at a
# saddle point because the tag-only likelihood surface is nearly flat.
tiny_tc  <- c(S1 = 2L, S2 = 1L, S3 = 1L)   # only 4 PBT-resolved fish
big_uc   <- c(S1 = 80L, S2 = 60L, S3 = 30L, Unassigned = 32L)
fit_tiny <- optim(rep(0, K_tag + 1), joint_ll,
                  tc = tiny_tc, uc = big_uc,
                  control = list(fnscale = -1, maxit = 5000), method = "BFGS")
cat("\nExtension: tiny tagged sample (n_tag = 4)\n")
cat("  MLE:", round(softmax(fit_tiny$par[1:K_tag]), 3),
    "  convergence:", fit_tiny$convergence, "\n")
# When convergence != 0, optim did not find the maximum — the estimate is
# driven almost entirely by the untagged data and the tag information is useless.

# Failure 2: the unassigned fraction can exceed 1 if the tag rates imply more
# unassigned fish than actually exist.  Happens when the tagged sample skews
# heavily toward one stock that is rare in the trap.
skewed_tc <- c(S1 = 90L, S2 = 5L, S3 = 5L)   # almost all tagged fish are S1
skewed_uc <- c(S1 = 20L, S2 = 40L, S3 = 30L, Unassigned = 32L)
fit_skew  <- optim(rep(0, K_tag + 1), joint_ll,
                   tc = skewed_tc, uc = skewed_uc,
                   control = list(fnscale = -1), method = "BFGS")
p_full_skew <- softmax(c(fit_skew$par[1:K_tag], fit_skew$par[K_tag + 1]))
cat("  Skewed tag sample: unassigned fraction =",
    round(p_full_skew[length(p_full_skew)], 3), "\n")
# A very large unassigned fraction signals that the tagged and untagged
# composition estimates are inconsistent — investigate the tagging protocol.

# ---- Forward pointer ----
# Part II applies these models to multi-year LGD data and introduces the
# hierarchical structure needed when stock proportions vary across years.
