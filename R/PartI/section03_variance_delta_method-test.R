#---------------------------------------------------------
# File:   section03_variance_delta_method-test.R
# Test driver for Part I Section 3
#---------------------------------------------------------

set.seed(2026)

section3_problem_3a_fish(d_a = 200L, p_n_true = 0.25,
                         a_d_fixed = 500L, nreps = 10000L)
section3_problem_3b_fish(a_d_true = 500L, r_sample = 5/6,
                         d_a = 200L, p_n_true = 0.25, nreps = 10000L)
section3_problem_3c_fish(a_d_true = 500L, r_sample = 5/6,
                         p_n_true = 0.25, d_a = 200L)
section3_problem_3d_fish(pf_true = 0.05, pre_f_true = 0.60, dt = 2000L)
section3_problem_3e_fish(d_a = 200L, p_n_true = 0.25, boots = 2000L)
