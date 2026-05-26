# Studio 2 - Binomial distributions (fish version)
# ---------------------------------------------------------------
# Source: https://ocw.mit.edu/courses/18-05-introduction-to-probability-and-statistics-spring-2022/mit18_05_s22_studio2-instructions.pdf
#
# Wraps MIT Studio 2 problems in fish form. Each coin toss is one fish at
# the ladder during the sampled portion of the hour. Studio 2 is also the
# tutorial for Section 1 of Part I, so this file is the canonical port.

library(ggplot2)

set.seed(2026)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# Problem 1a: simulate P(W = k) and P(W <= k) for window count W
studio2_problem_1a_fish <- function(n_fish, r_sample, k, ntrials) {
  ws <- rbinom(ntrials, size = n_fish, prob = r_sample)
  p_eq <- mean(ws == k)
  p_le <- mean(ws <= k)
  cat("Simulated P(W =", k, ")    =", round(p_eq, 4), "\n")
  cat("Simulated P(W <=", k, ")   =", round(p_le, 4), "\n")
  invisible(list(p_eq = p_eq, p_le = p_le))
}

# Problem 1b: exact P(W = k) via the pmf formula (no dbinom)
studio2_problem_1b_fish <- function(n_fish, r_sample, k) {
  p_eq <- choose(n_fish, k) * r_sample^k * (1 - r_sample)^(n_fish - k)
  p_le <- sum(choose(n_fish, 0:k) *
              r_sample^(0:k) *
              (1 - r_sample)^(n_fish - 0:k))
  cat("Exact     P(W =", k, ")    =", round(p_eq, 4), "\n")
  cat("Exact     P(W <=", k, ")   =", round(p_le, 4), "\n")
  invisible(list(p_eq = p_eq, p_le = p_le))
}

cat("Studio 2 Problem 1, window count W ~ Binomial(60, 5/6), k = 50:\n")
studio2_problem_1a_fish(n_fish = 60, r_sample = 5/6, k = 50, ntrials = 1e5)
studio2_problem_1b_fish(n_fish = 60, r_sample = 5/6, k = 50)
cat("\n")

# Problem 2a: plot the marked-fish payoff k^2 - 7k for k = 0..10
fish_payoff <- function(k) k^2 - 7 * k

studio2_problem_2a_fish <- function() {
  ntosses <- 10
  ks <- 0:ntosses
  payoffs <- fish_payoff(ks)
  png(file.path(plots_dir, "studio02_payoff.png"),
      width = 600, height = 400)
  plot(ks, payoffs, type = "h", lwd = 4, col = "steelblue",
       xlab = "k = marked fish in daily 10-fish draw",
       ylab = "payoff ($)",
       main = "Studio 2 - Marked-fish daily payoff")
  abline(h = 0)
  dev.off()
  invisible(payoffs)
}
studio2_problem_2a_fish()

# Problem 2b: exact expected payoff under p = 0.6 mark rate
studio2_problem_2b_fish <- function() {
  ntosses <- 10
  phead   <- 0.6
  ks <- 0:ntosses
  exact_value <- sum(fish_payoff(ks) * dbinom(ks, size = ntosses, prob = phead))
  cat("Exact E[payoff | p = 0.6] =", round(exact_value, 3), "\n")
  cat("Decision: ",
      ifelse(exact_value > 0, "GOOD bet for the supervisor",
                              "BAD bet for the supervisor"), "\n\n")
  invisible(exact_value)
}
studio2_problem_2b_fish()

# Problem 2c: simulate the average payoff
studio2_problem_2c_fish <- function(ntrials) {
  ntosses <- 10
  phead   <- 0.6
  payoffs <- fish_payoff(rbinom(ntrials, ntosses, phead))
  avg <- mean(payoffs)
  cat("Simulated average payoff over", ntrials, "days =", round(avg, 3), "\n\n")
  invisible(avg)
}
studio2_problem_2c_fish(ntrials = 1e5)

# Problem 3 (optional): derangements of PIT tag assignments
# If n tagging crews each have a target fish to retag, a derangement is a
# reassignment where no crew gets its own original target back.
studio2_problem_3_fish <- function(n, ntrials) {
  hits <- replicate(ntrials, {
    perm <- sample(seq_len(n))
    !any(perm == seq_len(n))
  })
  est <- mean(hits)
  cat("Estimated P(derangement) for n =", n, ":", round(est, 4),
      "  (1/e =", round(1/exp(1), 4), ")\n")
  invisible(est)
}
studio2_problem_3_fish(n = 10, ntrials = 1e5)
