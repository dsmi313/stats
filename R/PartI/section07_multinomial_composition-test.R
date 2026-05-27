library(tibble)
set.seed(2026)
true_w <- 0.45
true_hatch <- c(H1 = 0.30, H2 = 0.25)
tagRates_tib <- tibble(group = names(true_hatch), tagRate = c(0.70, 0.80))
n_total <- 400L
true_obs_probs <- c(true_hatch * tagRates_tib$tagRate, true_w + sum(true_hatch * (1 - tagRates_tib$tagRate)))
names(true_obs_probs) <- c(names(true_hatch), "Unassigned")
obs <- rmultinom(1, n_total, true_obs_probs)[,1]
nGroups <- obs[c("H1","H2")]
nUntag <- obs["Unassigned"]

pr <- section7_problem_7a_fish(c(0,0,0))
stopifnot(abs(sum(pr)-1) < 1e-8, all(pr>=0))
ll <- section7_problem_7b_fish(true_hatch, true_w, nGroups, nUntag, tagRates_tib$tagRate)
stopifnot(is.numeric(ll), length(ll)==1)
mle <- section7_problem_7d_fish(nGroups, nUntag, tagRates_tib)
stopifnot(abs(sum(mle$prop)-1) < 1e-6)
tex <- section7_problem_7e_fish(nGroups, nUntag, tagRates_tib)
stopifnot(abs(sum(tex$prop)-1) < 1e-6)
out <- section7_problem_7f_fish(nGroups, nUntag, tagRates_tib, boots = 200L)
stopifnot(nrow(out)==200L)
