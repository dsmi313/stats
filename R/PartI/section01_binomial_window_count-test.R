# Section 01 lightweight sanity checks
set.seed(2026)

N <- 400
p <- 5 / 6
nsims <- 2000

det_vec <- rbinom(N, 1, p)
W <- sum(det_vec)
stopifnot(length(det_vec) == N)
stopifnot(all(det_vec %in% c(0, 1)))
stopifnot(W == sum(det_vec))
stopifnot(W >= 0, W <= N)

W_sim <- replicate(nsims, sum(rbinom(N, 1, p)))
N_hat <- W_sim / p
stopifnot(length(W_sim) == nsims)
stopifnot(length(N_hat) == nsims)
stopifnot(abs(mean(N_hat) - N) < 20)

k_vals <- 0:N
db <- dbinom(k_vals, N, p)
stopifnot(all(db >= 0), abs(sum(db) - 1) < 1e-8)

k <- floor(N * p)
sim_cdf <- mean(W_sim <= k)
exact_cdf <- pbinom(k, N, p)
stopifnot(abs(sim_cdf - exact_cdf) < 0.05)

p_high <- 0.95
p_low <- 0.65
W_unequal <- replicate(nsims, sum(rbinom(N / 2, 1, p_high)) + sum(rbinom(N / 2, 1, p_low)))
stopifnot(length(W_unequal) == nsims, all(W_unequal >= 0 & W_unequal <= N))

N_blocks <- c(160, 140, 100)
p_blocks <- c(0.95, 0.75, 0.40)
W_timevarying <- replicate(nsims,
                           rbinom(1, N_blocks[1], p_blocks[1]) +
                           rbinom(1, N_blocks[2], p_blocks[2]) +
                           rbinom(1, N_blocks[3], p_blocks[3]))
stopifnot(length(W_timevarying) == nsims, all(W_timevarying >= 0 & W_timevarying <= N))
