# Part I Section 2 sanity checks
set.seed(2026)

est <- section2_problem_2a_fish(a_d_true = 500L, r_sample = 5/6, nreps = 3000L)
stopifnot(length(est) == 3000L, abs(mean(est) - 500) < 20)

w_obs <- rbinom(1, 500L, 5/6)
grid <- section2_problem_2b_fish(w_obs = w_obs, r_sample = 5/6, a_d_max = 800L)
stopifnot(length(grid$a_d_grid) == length(grid$loglik), grid$mle >= w_obs)

fit <- section2_problem_2c_fish(w_obs = w_obs, r_sample = 5/6)
stopifnot(is.list(fit), fit$convergence == 0)

wgt <- section2_problem_2d_fish(stocks = c("LOSALM","CHMBLN","IMNAHA"), strats = c("S1","S2","S3"), SR = 5/6*0.45, n_fish = 90L)
stopifnot(all(wgt$Primaryproportions >= 0), all(wgt$Primaryproportions <= 1))
