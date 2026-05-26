# Section 7 - Multinomial composition + escapeLGD's PBT_log_likelihood
# ---------------------------------------------------------------
# Goal: port escapeLGD's MLE composition likelihood verbatim and compare
# it to the "Accounting" (TotEx) approach used in PBT_expand_calc().
# Every function below has the exact signature you will see in
# escapeLGD/R/composition_estimation_utils.R.
#
# Repo pointers (escapeLGD/R/composition_estimation_utils.R):
#   line  15-18:  softMax(x)
#   lines 28-51:  PBT_expand_calc_MLE(values, tagRates)
#   lines 62-66:  PBT_log_likelihood(pGroups, pW, nGroups, nUntag, tagRates)
#                   y <- pW + sum((1 - tagRates) * pGroups)
#                   pGroups_tag <- pGroups * tagRates
#                   dmultinom(x = c(nGroups, nUntag),
#                             prob = c(pGroups_tag, y), log = TRUE)
#   lines 77-84:  PBT_optimllh(par, nGroups, nUntag, tagRates) -- softmax wrapper
#   lines 440-455: PBT_expand_calc(values, tagRates) -- accounting/TotEx method

library(ggplot2)
library(dplyr)
library(purrr)
library(tibble)

set.seed(2026)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# --- 1. softMax: numerically stable softmax (escapeLGD lines 15-18) --------
softMax <- function(x) {
  pr <- exp(x - max(x))
  pr / sum(pr)
}

# --- 2. PBT_log_likelihood (escapeLGD lines 62-66 verbatim) ----------------
PBT_log_likelihood <- function(pGroups, pW, nGroups, nUntag, tagRates) {
  y           <- pW + sum((1 - tagRates) * pGroups)
  pGroups_tag <- pGroups * tagRates
  dmultinom(x = c(nGroups, nUntag),
            prob = c(pGroups_tag, y), log = TRUE)
}

# --- 3. PBT_optimllh (escapeLGD lines 77-84 verbatim) ----------------------
PBT_optimllh <- function(par, nGroups, nUntag, tagRates) {
  par_prob <- softMax(par)
  pGroups  <- par_prob[1:(length(par_prob) - 1)]
  pW       <- par_prob[length(par_prob)]
  PBT_log_likelihood(pGroups = pGroups, pW = pW,
                     nGroups = nGroups, nUntag = nUntag,
                     tagRates = tagRates)
}

# --- 4. PBT_expand_calc_MLE (escapeLGD lines 28-51 verbatim) ---------------
PBT_expand_calc_MLE <- function(nGroups, nUntag, tagRates) {
  # nGroups: named vector of counts per PBT group
  # nUntag:  count of Unassigned (untagged) fish
  # tagRates: tibble with `group` and `tagRate` columns
  K <- length(nGroups)
  propsMLE <- optim(par = rep(1, K + 1),
                    fn  = PBT_optimllh,
                    nGroups = nGroups, nUntag = nUntag,
                    tagRates = tagRates$tagRate,
                    control = list(fnscale = -1, maxit = 1000),
                    method = "BFGS")
  if (propsMLE$convergence != 0) warning("Convergence error in PBT_expand_calc_MLE")
  par_prob <- softMax(propsMLE$par)
  tibble(group = c(tagRates$group, "Unassigned"), prop = par_prob)
}

# --- 5. PBT_expand_calc (TotEx accounting; escapeLGD lines 440-455) --------
PBT_expand_calc <- function(nGroups, nUntag, tagRates) {
  pr <- tibble(group   = tagRates$group,
               count   = as.numeric(nGroups),
               tagRate = tagRates$tagRate) |>
    mutate(expand = count / tagRate,
           diff   = expand - count)
  pr <- pr |> bind_rows(tibble(group = "Unassigned",
                               count = nUntag, tagRate = 1,
                               expand = 0, diff = 0))
  if (pr$count[pr$group == "Unassigned"] < sum(pr$diff)) {
    pr$diff   <- pr$count[pr$group == "Unassigned"] *
                 (pr$diff / sum(pr$diff))
    pr$expand <- pr$count + pr$diff
    pr$expand[pr$group == "Unassigned"] <- 0
  } else {
    pr$expand[pr$group == "Unassigned"] <-
      pr$expand[pr$group == "Unassigned"] - sum(pr$diff)
  }
  pr |> mutate(prop = expand / sum(expand)) |> select(group, prop)
}

# --- 6. Simulate a known mixture and recover with both methods ------------
# Two PBT hatchery groups + a wild pool. The wild pool also contains
# untagged hatchery fish (the bias D&H corrects for).
true_w        <- 0.45                 # truly wild fraction in run
true_hatch    <- c(H1 = 0.30, H2 = 0.25)  # truly hatchery in run
tagRates_tib  <- tibble(group   = names(true_hatch),
                        tagRate = c(0.70, 0.80))

n_total <- 400L
# Observed:
#   nGroups[i]: fish that PBT-assigned to hatchery group i (= hatch_i * t_i)
#   nUntag:   fish in the untagged pool (= sum(hatch_i*(1-t_i)) + w)
true_obs_probs <- c(true_hatch * tagRates_tib$tagRate,
                    true_w + sum(true_hatch * (1 - tagRates_tib$tagRate)))
names(true_obs_probs) <- c(names(true_hatch), "Unassigned")
obs <- rmultinom(1, size = n_total, prob = true_obs_probs)[, 1]
nGroups <- obs[c("H1", "H2")]
nUntag  <- obs["Unassigned"]

cat("Observed counts:  H1 =", nGroups["H1"],
    "  H2 =", nGroups["H2"], "  Unassigned =", nUntag, "\n")
cat("Truth:  wild prop =", true_w, "  H1 prop =", true_hatch["H1"],
    "  H2 prop =", true_hatch["H2"], "\n\n")

mle_fit <- PBT_expand_calc_MLE(nGroups, nUntag, tagRates_tib)
acc_fit <- PBT_expand_calc   (nGroups, nUntag, tagRates_tib)
cat("MLE method (D&H):\n");  print(mle_fit)
cat("\nAccounting (TotEx):\n"); print(acc_fit)
cat("\n")

# --- 7. Nonparametric bootstrap CIs (escapeLGD PBT_breakdown style) -------
# escapeLGD/R/composition_estimation_utils.R lines 466-486 reproduced.
PBT_breakdown <- function(nGroups, nUntag, tagRates, boots) {
  values <- c(rep(names(nGroups), nGroups), rep("Unassigned", nUntag))
  out <- matrix(NA_real_, nrow = boots,
                ncol = length(nGroups) + 1L,
                dimnames = list(NULL, c(names(nGroups), "Unassigned")))
  for (b in seq_len(boots)) {
    smp <- sample(values, length(values), replace = TRUE)
    counts <- table(factor(smp, levels = c(names(nGroups), "Unassigned")))
    fit <- PBT_expand_calc(counts[names(nGroups)],
                           counts["Unassigned"], tagRates)
    out[b, fit$group] <- fit$prop
  }
  out
}
boots <- PBT_breakdown(nGroups, nUntag, tagRates_tib, boots = 1500L)
cis <- apply(boots, 2, quantile, c(0.025, 0.975), na.rm = TRUE)
cat("Bootstrap 95% CIs from PBT_breakdown:\n")
print(round(cis, 3)); cat("\n")

# --- 8. Plot bootstrap distributions --------------------------------------
boots_long <- as_tibble(boots) |>
  tidyr::pivot_longer(everything(), names_to = "group", values_to = "p")
p_mn <- ggplot(boots_long, aes(p, fill = group)) +
  geom_histogram(bins = 50, alpha = 0.8) +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  labs(title = "Section 7 - Bootstrap of PBT_expand_calc",
       x = "proportion", y = "Bootstrap draws") +
  theme(legend.position = "none")
ggsave(file.path(plots_dir, "section07_pbt_bootstrap.png"), p_mn,
       width = 6, height = 7, dpi = 150)

# --- 9. End-of-section pointers --------------------------------------------
# You can now read every line of:
#   escapeLGD/R/composition_estimation_utils.R lines 1-110 (single-variable)
#   escapeLGD/R/composition_estimation_utils.R lines 114-223 (var2 nested)
# and recognize each function. Section 9 of PLAN.md picks up here and
# extends to the joint D&H likelihood with full (PBT, GenStock) structure.
