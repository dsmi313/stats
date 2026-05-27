#---------------------------------------------------------
# File:   section07_multinomial_composition-test.R
# Test driver for Part I Section 7
#---------------------------------------------------------
set.seed(2026)

# Problem 7a: simulate three stocks; MLE = observed proportions
section7_problem_7a_fish(n = 200, props = c(0.6, 0.3, 0.1))

# Problem 7b: multinomial log-likelihood and optim() MLE
section7_problem_7b_fish(counts = c(120L, 60L, 20L))

# Problem 7c: joint log-likelihood for two multinomial samples
section7_problem_7c_fish()

# Problem 7d: nonparametric bootstrap CI
section7_problem_7d_fish(counts = c(S1 = 120L, S2 = 60L, S3 = 20L), boots = 5000L)
