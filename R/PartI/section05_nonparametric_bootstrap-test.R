#---------------------------------------------------------
# File:   section05_nonparametric_bootstrap-test.R
# Test driver for Part I Section 5
#---------------------------------------------------------

set.seed(2026)

stocks <- c("LOSALM", "CHMBLN", "IMNAHA")
strats <- c("S1", "S2", "S3")

# Problem 5a: inverse-SR weighting (Horvitz-Thompson step)
section5_problem_5a_fish(stocks = stocks, strats = strats,
                         SR     = 5/6 * 0.45, n_fish = 90L)

# Problem 5b: build toy SCRAPI-shaped data frames
toy <- section5_problem_5b_fish(stocks = stocks, strats = strats, nfish = 90L)

# Problems 5c/5d: weighted resamples (one iteration to spot-check)
section5_problem_5c_fish(FishWH  = toy$RearData, strats = strats)
section5_problem_5d_fish(FishDat = toy$AllPrime, strats = strats)

# Problem 5e: thetahat_toy on the real data
section5_problem_5e_fish(passage = toy$passdata,
                         RearDat = toy$RearData,
                         Fish    = toy$AllPrime,
                         Pgrps   = stocks,
                         strats  = strats)

# Problem 5f: full bootsmolt loop
section5_problem_5f_fish(passdata = toy$passdata,
                         RearData = toy$RearData,
                         AllPrime = toy$AllPrime,
                         stocks   = stocks,
                         strats   = strats,
                         B        = 1000L)

# Problem 5g: trigger Error 1
section5_problem_5g_fish()
