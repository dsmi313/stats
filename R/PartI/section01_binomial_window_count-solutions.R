#---------------------------------------------------------
# File:   section01_binomial_window_count-solutions.R
# Part I, Section 1 concise instructor key
#---------------------------------------------------------
set.seed(2026)

N <- 400
p <- 5 / 6
nsims <- 5000

# Fish-level Bernoulli day and window count
Det <- rbinom(N, 1, p)
W <- sum(Det)

# Repeated Binomial days and abundance expansion
W_sim <- replicate(nsims, sum(rbinom(N, 1, p)))
N_hat <- W_sim / p

cat("One-day W:", W, "\n")
cat("Mean(W_sim) vs Np:", round(mean(W_sim), 2), "vs", N * p, "\n")
cat("Var(W_sim) vs Np(1-p):", round(var(W_sim), 2), "vs", N * p * (1 - p), "\n")
cat("Mean(N_hat) vs N:", round(mean(N_hat), 2), "vs", N, "\n\n")

par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
plot(seq_len(N), Det, pch = 16, cex = 0.6,
     xlab = "Fish index", ylab = "Detected (0/1)", main = "One day")
hist(W_sim, breaks = 35, col = "lightsteelblue", border = "white", main = "W_sim", xlab = "W")
abline(v = N * p, col = 2, lwd = 2)
hist(N_hat, breaks = 35, col = "honeydew3", border = "white", main = "N_hat", xlab = "W/p")
abline(v = N, col = 2, lwd = 2)
par(mfrow = c(1, 1))

# Exact-vs-simulation check
k <- 0:N
sim_prob <- tabulate(W_sim + 1, nbins = N + 1) / nsims
exact_prob <- dbinom(k, N, p)
plot(k, exact_prob, type = "l", lwd = 2, main = "dbinom vs simulation", xlab = "k", ylab = "Probability")
points(k, sim_prob, pch = 16, cex = 0.4, col = "dodgerblue3")

k0 <- floor(N * p)
cat("P(W <=", k0, ") sim:", round(mean(W_sim <= k0), 4),
    " exact:", round(pbinom(k0, N, p), 4), "\n\n")

# Failure lab 1: unequal p
p_high <- 0.95; p_low <- 0.65
W_unequal <- replicate(nsims, sum(rbinom(N/2, 1, p_high)) + sum(rbinom(N/2, 1, p_low)))
p_bar <- (p_high + p_low) / 2
W_naive <- rbinom(nsims, N, p_bar)
cat("Unequal p var vs naive:", round(var(W_unequal), 2), "vs", round(var(W_naive), 2), "\n")

# Failure lab 2: time-varying p
N_blocks <- c(160, 140, 100)
p_blocks <- c(0.95, 0.75, 0.40)
W_timevarying <- replicate(nsims,
                           rbinom(1, N_blocks[1], p_blocks[1]) +
                           rbinom(1, N_blocks[2], p_blocks[2]) +
                           rbinom(1, N_blocks[3], p_blocks[3]))
p_bar_time <- sum(N_blocks * p_blocks) / N
W_naive_time <- rbinom(nsims, N, p_bar_time)
cat("Time-varying p var vs naive:", round(var(W_timevarying), 2), "vs", round(var(W_naive_time), 2), "\n")
