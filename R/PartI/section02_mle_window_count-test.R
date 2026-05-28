#---------------------------------------------------------
# File:   section02_mle_window_count-test.R
# Test driver for Part I Section 2
# Source section02_mle_window_count.R or its solutions before running.
#---------------------------------------------------------

set.seed(2026)

# Problem 2a: replicated MLE centering check
section2_problem_2a_fish(a_d_true = 500L, r_sample = 5/6, nreps = 1000L)

# Problem 2b: grid evaluation of the log-likelihood for one observed w
w_obs <- rbinom(1, 500L, 5/6)
section2_problem_2b_fish(w_obs = w_obs, r_sample = 5/6, a_d_max = 700L)

# Problem 2c: optim continuous-relaxation MLE
section2_problem_2c_fish(w_obs = w_obs, r_sample = 5/6)
