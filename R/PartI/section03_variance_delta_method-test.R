# Part I Section 3 sanity checks
set.seed(2026)

l1 <- section3_problem_3a_fish(d_a = 40L, p_n_true = 0.25, a_d_fixed = 500, nreps = 2000L)
stopifnot(length(l1) == 2000L)

l2 <- section3_problem_3b_fish(a_d_true = 500L, r_sample = 5/6, d_a = 40L, p_n_true = 0.25, nreps = 2000L)
stopifnot(length(l2) == 2000L)

dm <- section3_problem_3c_fish(a_d_true = 500L, r_sample = 5/6, p_n_true = 0.25, d_a = 40L)
stopifnot(dm$var_a_t_delta > 0)

fb <- section3_problem_3d_fish(pf_true = 0.35, pre_f_true = 0.6, dt = 1000L)
stopifnot(all(fb$par >= 0 & fb$par <= 1))

bt <- section3_problem_3e_fish(a_d_true = 500L, p_n_true = 0.25, d_a = 40L, r_sample = 5/6, boots = 1000L)
stopifnot(length(bt$ci) == 2L)
