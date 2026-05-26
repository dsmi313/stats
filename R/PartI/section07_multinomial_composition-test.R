#---------------------------------------------------------
# File:   section07_multinomial_composition-test.R
# Test driver for Part I Section 7
#---------------------------------------------------------

library(tibble)
set.seed(2026)

# Truth
true_w       <- 0.45
true_hatch   <- c(H1 = 0.30, H2 = 0.25)
tagRates_tib <- tibble(group = names(true_hatch),
                       tagRate = c(0.70, 0.80))

# Simulate observed counts
n_total <- 400L
true_obs_probs <- c(true_hatch * tagRates_tib$tagRate,
                    true_w + sum(true_hatch * (1 - tagRates_tib$tagRate)))
names(true_obs_probs) <- c(names(true_hatch), "Unassigned")
obs     <- rmultinom(1, size = n_total, prob = true_obs_probs)[, 1]
nGroups <- obs[c("H1", "H2")]
nUntag  <- obs["Unassigned"]
cat("Observed counts:\n"); print(obs); cat("\n")

# Problem 7a: softMax
print(section7_problem_7a_fish(c(0, 0, 0)))

# Problem 7b: log-likelihood at the truth
ll_truth <- section7_problem_7b_fish(
  pGroups = true_hatch, pW = true_w,
  nGroups = nGroups, nUntag = nUntag,
  tagRates = tagRates_tib$tagRate)
cat("Log-likelihood at truth =", ll_truth, "\n")

# Problem 7d: PBT_expand_calc_MLE
section7_problem_7d_fish(nGroups, nUntag, tagRates_tib)

# Problem 7e: PBT_expand_calc (TotEx)
section7_problem_7e_fish(nGroups, nUntag, tagRates_tib)

# Problem 7f: PBT_breakdown nonparametric bootstrap
section7_problem_7f_fish(nGroups, nUntag, tagRates_tib, boots = 500L)
