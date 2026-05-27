#---------------------------------------------------------
# File:   section01_binomial_window_count.R
# Part I, Section 1 - Binomial window count (adult ladder passage)
# Workshop flow: simulate -> inspect -> estimate -> diagnose -> explain -> EASE
#---------------------------------------------------------

# What this block teaches:
# We begin at fish level. Each fish is either detected (1) or not (0) during
# a window sample. That is a Bernoulli trial with probability p.
# Summing fish-level 0/1 outcomes gives the daily window count W.

set.seed(2026)
N <- 400
p <- 5 / 6
nsims <- 5000

det_vec <- rbinom(N, 1, p)
W <- sum(det_vec)

# Side-by-side diagnostic requested for workshop intuition.
# Panel 1: fish-level Bernoulli outcomes for one day.
# Panel 2: repeated Binomial window counts across many days.
# Panel 3: expanded abundance estimator N_hat = W / p.
W_sim <- replicate(nsims, sum(rbinom(N, 1, p)))
N_hat <- W_sim / p

op <- par(no.readonly = TRUE)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
plot(seq_len(N), det_vec,
     pch = 16, cex = 0.6,
     xlab = "Fish index", ylab = "Detected (0/1)",
     main = "One day: fish-level Bernoulli")
hist(W_sim, breaks = 35, col = "lightsteelblue", border = "white",
     xlab = "W (window count)", main = "Repeated W")
abline(v = N * p, col = "firebrick", lwd = 2)
hist(N_hat, breaks = 35, col = "honeydew3", border = "white",
     xlab = "N_hat = W / p", main = "Repeated abundance estimate")
abline(v = N, col = "firebrick", lwd = 2)
par(op)

cat("Section 01 workshop setup\n")
cat("N =", N, " p =", p, " one-day W =", W, "\n")
cat("Mean(W_sim) =", round(mean(W_sim), 2), " Theoretical N*p =", N * p, "\n")
cat("Var(W_sim)  =", round(var(W_sim), 2), " Theoretical N*p*(1-p) =", N * p * (1 - p), "\n")
cat("Mean(N_hat) =", round(mean(N_hat), 2), " True N =", N, "\n\n")

# Exact vs simulation:
# The Binomial model gives exact probabilities with dbinom/pbinom.
# Simulation should closely approximate those exact probabilities.
k_vals <- 0:N
sim_freq <- tabulate(W_sim + 1, nbins = N + 1) / nsims
exact_prob <- dbinom(k_vals, size = N, prob = p)

plot(k_vals, exact_prob, type = "l", lwd = 2, col = "black",
     xlab = "k", ylab = "Probability",
     main = "Exact Binomial vs simulation")
points(k_vals, sim_freq, pch = 16, cex = 0.45, col = "dodgerblue3")
legend("topright", legend = c("Exact dbinom", "Simulated frequency"),
       col = c("black", "dodgerblue3"), lty = c(1, NA), pch = c(NA, 16), bty = "n")

k_check <- floor(N * p)
sim_cdf <- mean(W_sim <= k_check)
exact_cdf <- pbinom(k_check, size = N, prob = p)
cat("P(W <=", k_check, ") simulated =", round(sim_cdf, 4),
    " exact =", round(exact_cdf, 4), "\n\n")

# Failure mode 1: unequal p among fish.
# Half fish are easy to detect, half are harder.
# This breaks the identical-p assumption of Binomial(N, p).
p_high <- 0.95
p_low <- 0.65
N_half <- N / 2
W_unequal <- replicate(nsims,
                      sum(rbinom(N_half, 1, p_high)) +
                      sum(rbinom(N_half, 1, p_low)))
p_bar <- (p_high + p_low) / 2
W_naive <- rbinom(nsims, size = N, prob = p_bar)

tail_k <- floor(N * p_bar + 15)
cat("Unequal-p lab\n")
cat("mean(W_unequal) =", round(mean(W_unequal), 2),
    " mean(W_naive) =", round(mean(W_naive), 2), "\n")
cat("var(W_unequal)  =", round(var(W_unequal), 2),
    " var(W_naive)  =", round(var(W_naive), 2), "\n")
cat("P(W >=", tail_k, ") unequal =", round(mean(W_unequal >= tail_k), 4),
    " naive =", round(mean(W_naive >= tail_k), 4), "\n\n")

op <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
hist(W_unequal, breaks = 35, col = "khaki2", border = "white",
     main = "Unequal p across fish", xlab = "W")
hist(W_naive, breaks = 35, col = "gray75", border = "white",
     main = "Naive Binomial(N, p_bar)", xlab = "W")
par(op)

# Failure mode 2: time-varying p across passage blocks.
# Connects to spill/flow/time conditions that change observation probability.
p_blocks <- c(0.95, 0.75, 0.40)
N_blocks <- c(160, 140, 100) # sums to 400
stopifnot(sum(N_blocks) == N)

W_timevarying <- replicate(nsims,
  sum(rbinom(1, N_blocks[1], p_blocks[1]),
      rbinom(1, N_blocks[2], p_blocks[2]),
      rbinom(1, N_blocks[3], p_blocks[3])))
p_bar_time <- sum(N_blocks * p_blocks) / N
W_naive_time <- rbinom(nsims, size = N, prob = p_bar_time)

tail_k2 <- floor(N * p_bar_time + 20)
cat("Time-varying p lab\n")
cat("mean(W_timevarying) =", round(mean(W_timevarying), 2),
    " mean(W_naive_time) =", round(mean(W_naive_time), 2), "\n")
cat("var(W_timevarying)  =", round(var(W_timevarying), 2),
    " var(W_naive_time)  =", round(var(W_naive_time), 2), "\n")
cat("P(W >=", tail_k2, ") time-varying =", round(mean(W_timevarying >= tail_k2), 4),
    " naive =", round(mean(W_naive_time >= tail_k2), 4), "\n")

op <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
hist(W_timevarying, breaks = 35, col = "lightcoral", border = "white",
     main = "Time-varying p blocks", xlab = "W")
hist(W_naive_time, breaks = 35, col = "gray75", border = "white",
     main = "Naive Binomial(N, p_bar)", xlab = "W")
par(op)

# EASE connection:
# In EASE-style expansion, N_hat = W/p assumes p is known, constant, and fish
# detections behave like independent Bernoulli draws. Failure labs show why
# misspecified p can give misleading uncertainty and tail behavior.

# Section summary:
# 1. Biological process simulated: fish detections at a ladder window.
# 2. Stochastic model: independent Bernoulli fish detections with count W.
# 3. Why Binomial: W is a sum of Bernoulli( p ) trials.
# 4. Why W/p estimates abundance: E[W]=Np so W/p centers near N.
# 5. Broken assumptions: unequal p across fish and time-varying p by blocks.
# 6. Why this matters for EASE: fixed-p expansion can miss variance/tail risk.
