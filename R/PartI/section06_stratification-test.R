library(tibble)
library(dplyr)
set.seed(2026)
season <- section6_problem_6a_fish(n_weeks = 8L, days_per_wk = 7L, p_n_early = 0.1, p_n_late = 0.4)
res <- section6_problem_6b_fish(season)
stopifnot(res$a_t_pooled > 0, res$a_t_stratified > 0)

nFish <- 200L
fish_weeks <- sample(season$Week, size = nFish, replace = TRUE)
All <- tibble(CollectionDate = format(as.Date("2024-04-01") + (fish_weeks - 1L) * 7L, "%m/%d/%Y"), GenStock = sample(c("LOSALM","CHMBLN","IMNAHA"), nFish, TRUE))
pass <- season |> mutate(SampleEndDate = format(as.Date("2024-04-01") + (Week - 1L) * 7L, "%m/%d/%Y")) |> group_by(Week, Collaps, SampleEndDate) |> summarise(SampleCount = sum(a_d), .groups = "drop")

mapped <- section6_problem_6c_fish(pass, All)
stopifnot(nrow(mapped) == nFish)
section6_problem_6d_fish(pass, All)
ci <- section6_problem_6e_fish(season, B = 500L)
stopifnot(length(ci$ci_pool) == 2L, length(ci$ci_strat) == 2L)
