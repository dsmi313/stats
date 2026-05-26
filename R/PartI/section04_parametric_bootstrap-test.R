#---------------------------------------------------------
# File:   section04_parametric_bootstrap-test.R
# Test driver for Part I Section 4
#---------------------------------------------------------

set.seed(2026)

# A 60-day smolt-trap toy season (SCRAPI passdata layout).
LGDdaily <- data.frame(
  Stratum = rep(seq_len(6), each = 10L),
  Tally   = rbinom(60, size = 200L, prob = (5/6) * 0.45),
  Ptrue   = 5/6 * 0.45    # SampleRate * GuidanceEfficiency (smolt trap)
)

# Problem 4a: SCRAPI bootsmolt daily-count loop (slow, line-by-line port)
out_a <- section4_problem_4a_fish(LGDdaily = LGDdaily, B = 500L)

# Problem 4b: vectorized escapeLGD equivalent
out_b <- section4_problem_4b_fish(LGDdaily = LGDdaily, B = 5000L)

# Problem 4c: coverage of 95% CI across simulated seasons
section4_problem_4c_fish(ndays = 20L, true_pass_per_day = 200L,
                         Ptrue = 5/6 * 0.45,
                         B = 1000L, nseasons = 100L)

# Problem 4d: plot bootstrap distribution of season total
section4_problem_4d_fish(LGDdaily = LGDdaily, B = 5000L)
