#---------------------------------------------------------
# File:   section04_parametric_bootstrap-solutions.R
# Part I, Section 4 solutions
#---------------------------------------------------------
set.seed(2026)

# Section 4 solutions ----

library(ggplot2)
library(tibble)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

#--------------------------------------
# Problem 4a: bootsmolt's daily-count step (SCRAPI verbatim)
section4_problem_4a_fish <- function(LGDdaily, B) {
  cat("\n----------------------------------\n")
  cat("Problem 4a: bootsmolt daily-count step (SCRAPI lines 139-145)\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  dailypass <- round(LGDdaily$Tally / LGDdaily$Ptrue)
  ndays     <- nrow(LGDdaily)
  out       <- matrix(NA_real_, nrow = B, ncol = ndays)
  for (b in seq_len(B)) {
    cntstar <- numeric(ndays)
    for (i in seq_len(ndays)) {
      if (dailypass[i] != 0) {
        cntstar[i] <- rbinom(1, dailypass[i], LGDdaily$Ptrue[i])
      }
    }
    dailyStar <- data.frame(Stratum = LGDdaily$Stratum,
                            Tally   = cntstar,
                            Ptrue   = LGDdaily$Ptrue)
    out[b, ] <- dailyStar$Tally / dailyStar$Ptrue
  }
  cat("B =", B, "  ndays =", ndays, "  point total =",
      round(sum(LGDdaily$Tally / LGDdaily$Ptrue)), "\n")
  invisible(out)
}


# Problem 4b: Vectorized escapeLGD-style daily-count bootstrap
section4_problem_4b_fish <- function(LGDdaily, B) {
  cat("\n----------------------------------\n")
  cat("Problem 4b: Vectorized escapeLGD daily-count bootstrap\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  dailypass <- round(LGDdaily$Tally / LGDdaily$Ptrue)
  out <- sapply(seq_len(nrow(LGDdaily)), function(i) {
    rbinom(B, dailypass[i], LGDdaily$Ptrue[i]) / LGDdaily$Ptrue[i]
  })
  cat("B =", B, "  ndays =", nrow(LGDdaily), "\n")
  invisible(out)
}


# Problem 4c: Coverage check across simulated seasons
section4_problem_4c_fish <- function(ndays, true_pass_per_day, Ptrue,
                                     B, nseasons) {
  cat("\n----------------------------------\n")
  cat("Problem 4c: 95% bootstrap-CI coverage across simulated seasons\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  hits <- vapply(seq_len(nseasons), function(s) {
    set.seed(s)
    Tally <- rbinom(ndays, true_pass_per_day, Ptrue)
    LD <- data.frame(Stratum = 1L, Tally = Tally, Ptrue = Ptrue)
    draws  <- section4_problem_4b_fish(LD, B = B)
    totals <- rowSums(draws)
    ci     <- quantile(totals, c(0.025, 0.975))
    truth  <- ndays * true_pass_per_day
    (truth >= ci[1]) && (truth <= ci[2])
  }, logical(1))
  cat("ndays =", ndays, "  true_pass_per_day =", true_pass_per_day,
      "  Ptrue =", Ptrue, "  B =", B,
      "  nseasons =", nseasons, "\n")
  cat("Coverage of 95% CI =", mean(hits), "\n")
  invisible(hits)
}


# Problem 4d: Plot the bootstrap distribution of the season total
section4_problem_4d_fish <- function(LGDdaily, B) {
  cat("\n----------------------------------\n")
  cat("Problem 4d: Plot bootstrap distribution of season total\n")

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

  draws  <- section4_problem_4b_fish(LGDdaily, B = B)
  totals <- rowSums(draws)
  ci     <- quantile(totals, c(0.025, 0.975))

  p <- ggplot(tibble(total = totals), aes(total)) +
    geom_histogram(bins = 60, fill = "steelblue", colour = "white") +
    geom_vline(xintercept = ci, colour = "firebrick",
               linewidth = 1, linetype = "dashed") +
    labs(title = "Section 4 - Bootstrap season totals (bootsmolt inner loop)",
         subtitle = "Dashed lines = 95% CI",
         x = "season total smolt passage", y = "Bootstrap draws")
  ggsave(file.path(plots_dir, "section04_boot_season.png"), p,
         width = 6, height = 4, dpi = 150)
  cat("Saved plot to docs/figures/PartI/section04_boot_season.png\n")
  invisible(totals)
}
