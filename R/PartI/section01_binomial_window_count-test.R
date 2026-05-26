#---------------------------------------------------------
# File:   section01_binomial_window_count-test.R
# Test driver for Part I Section 1
# Before running: clean your environment and source either
#   section01_binomial_window_count.R           (your answers)
# or
#   section01_binomial_window_count-solutions.R (reference answers)
#---------------------------------------------------------

set.seed(2026)

# Problem 1a: one day at the ladder window
section1_problem_1a_fish(a_d_true = 500L, r_sample = 5/6)

# Problem 1b: build a 12-week wc tibble
wc <- section1_problem_1b_fish(n_weeks = 12L, wc_prop = 5/6,
                                mean_truth = 3500)

# Problem 1c: bootstrap the season total from the wc table built in 1b
section1_problem_1c_fish(wc = wc, wc_prop = 5/6, boots = 5000L)

# Problem 1d: replicated unbiasedness check
section1_problem_1d_fish(a_d_true = 500L, r_sample = 5/6, nreps = 1000L)
