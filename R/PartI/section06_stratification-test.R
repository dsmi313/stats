#---------------------------------------------------------
# File:   section06_stratification-test.R
# Test driver for Part I Section 6
#---------------------------------------------------------

library(tibble)
library(dplyr)
set.seed(2026)

# Problem 6a: simulate a season with shifting p_n
season <- section6_problem_6a_fish(n_weeks = 8L, days_per_wk = 7L,
                                   p_n_early = 0.10, p_n_late = 0.40)

# Problem 6b: Cpattern + pooled vs stratified
section6_problem_6b_fish(season = season)

# Build a toy trap data set (All) and the SCRAPI-style pass tibble so we
# can exercise the Collaps assignment loop.
nFish <- 200L
fish_weeks <- sample(season$Week, size = nFish, replace = TRUE)
All <- tibble(
  CollectionDate = format(as.Date("2024-04-01") + (fish_weeks - 1L) * 7L,
                          "%m/%d/%Y"),
  GenStock       = sample(c("LOSALM", "CHMBLN", "IMNAHA"),
                          size = nFish, replace = TRUE)
)
pass <- season |>
  mutate(SampleEndDate = format(as.Date("2024-04-01") + (Week - 1L) * 7L,
                                "%m/%d/%Y")) |>
  group_by(Week, Collaps, SampleEndDate) |>
  summarise(SampleCount = sum(a_d), .groups = "drop")

# Problem 6c: SCRAPI Collaps assignment
section6_problem_6c_fish(pass = pass, All = All)

# Problem 6d: trigger Error 3 with duplicated FPC date
section6_problem_6d_fish(pass = pass, All = All)

# Problem 6e: parametric-bootstrap CI comparison
section6_problem_6e_fish(season = season, B = 1500L)
