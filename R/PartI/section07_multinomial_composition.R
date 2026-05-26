# Section 7 - Multinomial composition
# ---------------------------------------------------------------
# Generalize the Binomial (2 outcomes) to the Multinomial (many).
# This is the distribution under every composition estimator in EASE
# and SCRAPI, and the skeleton of the Delomas & Hess (2021) likelihood
# you will fill in during Sections 8-10.
#
# Reading: PLAN.md Section 7, MIT Class 4b (Expected Value),
# Class 6a (Continuous RVs), Class 7a (Joint Distributions), R Studio 8.

library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)

set.seed(2026)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# --- 1. Simulate a multinomial trap sample --------------------------------
true_p <- c(LOSALM = 0.6, CHMBLN = 0.3, IMNAHA = 0.1)
n      <- 200L
counts <- as.numeric(rmultinom(1, size = n, prob = true_p))
names(counts) <- names(true_p)

cat("Observed counts:\n");                       print(counts)
cat("Observed proportions (closed-form MLE):\n"); print(round(counts / n, 3))
cat("Truth:\n");                                  print(true_p); cat("\n")

# --- 2. Multinomial log-likelihood by hand + optim() ---------------------
multinom_loglik <- function(p, counts) {
  if (any(p <= 0) || abs(sum(p) - 1) > 1e-6) return(-Inf)
  sum(counts * log(p))
}

# Reparametrize to K-1 free parameters via softmax so optim() can search
# the unconstrained space and we still get a valid simplex back.
neg_loglik_optim <- function(par, counts) {
  z <- c(par, 0)
  p <- exp(z) / sum(exp(z))
  -multinom_loglik(p, counts)
}

fit <- optim(par = c(0, 0),
             fn  = neg_loglik_optim,
             counts = counts,
             method = "BFGS")
p_mle <- {
  z <- c(fit$par, 0); exp(z) / sum(exp(z))
}
names(p_mle) <- names(true_p)
cat("optim()-recovered MLE proportions:\n"); print(round(p_mle, 3)); cat("\n")
cat("Agreement with closed-form MLE: max abs diff =",
    max(abs(p_mle - counts / n)), "\n\n")

# --- 3. Product of two multinomials: tagged + untagged --------------------
# Tagged fish: known hatchery group via PBT. Two hatcheries (H1, H2).
# Untagged pool: a single multinomial across stocks with cell probabilities
# y_j that combine wild and untagged hatchery contributions (Section 9 builds
# y_j from w_j, t_i, D_ij; here we use a hand-picked vector that sums to 1).
true_hatch_p <- c(H1 = 0.55, H2 = 0.45)    # tagged hatchery proportions
y_j          <- c(LOSALM = 0.55, CHMBLN = 0.30, IMNAHA = 0.15)
stopifnot(abs(sum(true_hatch_p) - 1) < 1e-9,
          abs(sum(y_j) - 1) < 1e-9)

n_tagged   <- 80L
n_untagged <- 120L
tagged_counts   <- as.numeric(rmultinom(1, size = n_tagged,   prob = true_hatch_p))
untagged_counts <- as.numeric(rmultinom(1, size = n_untagged, prob = y_j))
names(tagged_counts)   <- names(true_hatch_p)
names(untagged_counts) <- names(y_j)

cat("Tagged counts:\n");   print(tagged_counts)
cat("Untagged counts:\n"); print(untagged_counts); cat("\n")

joint_loglik <- function(p_hatch, p_y, c_tag, c_unt) {
  multinom_loglik(p_hatch, c_tag) + multinom_loglik(p_y, c_unt)
}

cat("Joint log-likelihood at the truth =",
    joint_loglik(true_hatch_p, y_j, tagged_counts, untagged_counts), "\n\n")

# --- 4. Joint MLE via optim() --------------------------------------------
neg_joint_loglik <- function(par, c_tag, c_unt) {
  # par = c(logit hatchery free, two softmax free parameters for y_j)
  p_h <- {
    z <- c(par[1], 0); exp(z) / sum(exp(z))
  }
  p_y <- {
    z <- c(par[2:3], 0); exp(z) / sum(exp(z))
  }
  -(multinom_loglik(p_h, c_tag) + multinom_loglik(p_y, c_unt))
}

joint_fit <- optim(par = c(0, 0, 0),
                   fn  = neg_joint_loglik,
                   c_tag = tagged_counts,
                   c_unt = untagged_counts,
                   method = "BFGS")
p_hatch_mle <- {
  z <- c(joint_fit$par[1], 0); exp(z) / sum(exp(z))
}
p_y_mle <- {
  z <- c(joint_fit$par[2:3], 0); exp(z) / sum(exp(z))
}
names(p_hatch_mle) <- names(true_hatch_p)
names(p_y_mle)     <- names(y_j)

cat("Joint MLE hatchery proportions:\n"); print(round(p_hatch_mle, 3))
cat("Joint MLE untagged y_j:\n");           print(round(p_y_mle, 3)); cat("\n")

# --- 5. Nonparametric bootstrap CIs for multinomial proportions ----------
B <- 5000L
boot_props <- replicate(B, {
  smpl <- sample(rep(names(true_p), counts), replace = TRUE)
  vapply(names(true_p), function(s) mean(smpl == s), numeric(1))
})
cis <- apply(boot_props, 1, quantile, c(0.025, 0.975))
cat("Bootstrap 95% CIs for stock proportions:\n")
print(round(rbind(MLE = counts / n, cis), 3)); cat("\n")

# --- 6. Visualize the bootstrap distribution -----------------------------
draws_tbl <- as_tibble(t(boot_props)) |>
  pivot_longer(everything(), names_to = "stock", values_to = "p")

p_mn <- ggplot(draws_tbl, aes(p, fill = stock)) +
  geom_histogram(bins = 40, alpha = 0.8) +
  facet_wrap(~ stock, ncol = 1, scales = "free_y") +
  labs(title = "Section 7 - Bootstrap stock proportions",
       x = "proportion", y = "Bootstrap draws") +
  theme(legend.position = "none")
ggsave(file.path(plots_dir, "section07_bootstrap_props.png"), p_mn,
       width = 6, height = 7, dpi = 150)

# Section 7 payoff ---------------------------------------------------------
# Tagged + untagged as two multinomials with shared parameters is the
# Delomas & Hess (2021) likelihood skeleton. In Section 9 you will replace
# the placeholder y_j with y_j = w_j + sum_i(p_i * (1 - t_i) * D_ij) and
# recover the wild proportions w_j from a joint optim().
