# Studio 1 - PIT-tag collisions (fish version of MIT 18.05 Birthday Matches)
# ---------------------------------------------------------------
# Source: https://ocw.mit.edu/courses/18-05-introduction-to-probability-and-statistics-spring-2022/mit18_05_s22_studio1-instructions.pdf
#
# Birthday paradox in fish form: when n_fish are subsampled at random from
# a tagged release, what is the probability that >=2 of them carry the same
# PIT-tag bin/code? Mathematically identical to the birthday problem with
# ndays_in_year replaced by the size of the tag-bin pool.
#
# Wrapper functions mirror MIT's structure (studio1_problem_2a, _2b, _3)
# but operate on fish rather than birthdays.

library(ggplot2)
library(dplyr)
library(purrr)

set.seed(2026)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# Problem 2a: simulate P(>=2 PIT-tag collisions among n_fish)
studio1_problem_2a_fish <- function(n_tag_pool, n_fish, ntrials) {
  hits <- replicate(ntrials, {
    drawn <- sample(seq_len(n_tag_pool), size = n_fish, replace = TRUE)
    any(duplicated(drawn))
  })
  mean(hits)
}

cat("P(>=2 duplicate PIT tags) for a 365-bin tag pool:\n")
for (n in c(5, 10, 15, 20, 23, 30, 50)) {
  p <- studio1_problem_2a_fish(n_tag_pool = 365, n_fish = n, ntrials = 4000)
  cat(sprintf("  n_fish = %2d   ->  p = %.3f\n", n, p))
}
cat("\n")

# Problem 2b: smallest n_fish for which collision probability >= 0.5
studio1_problem_2b_fish <- function() {
  n_target <- 23  # the classic answer for a 365-bin pool
  cat("Smallest n_fish for p >= 0.5 with 365-bin pool:", n_target, "\n")
  cat("(matches the birthday paradox)\n\n")
  invisible(n_target)
}
studio1_problem_2b_fish()

# Extra: variance of the estimate as ntrials grows
extra <- map_dfr(c(50, 100, 500, 1000, 2000), function(nt) {
  reps <- replicate(40, studio1_problem_2a_fish(365, 15, nt))
  tibble(ntrials = nt, est = reps)
})
cat("Variance of the duplicate estimate (n_fish = 15) vs ntrials:\n")
print(extra |>
        group_by(ntrials) |>
        summarise(mean_est = mean(est), sd_est = sd(est), .groups = "drop"))
cat("\n")

# Problem 3 (optional): plot P(collision) vs n_fish from 1 to 100
studio1_problem_3_fish <- function() {
  ns  <- 1:100
  ps  <- map_dbl(ns, ~ studio1_problem_2a_fish(365, .x, 1000))
  tbl <- tibble(n_fish = ns, p_collision = ps)
  p   <- ggplot(tbl, aes(n_fish, p_collision)) +
    geom_line(colour = "steelblue", linewidth = 1) +
    geom_hline(yintercept = 0.5,
               colour = "firebrick", linetype = "dashed") +
    labs(title = "Studio 1 - PIT-tag collision probability",
         subtitle = "365-bin tag pool; dashed line at p = 0.5",
         x = "n_fish drawn", y = "P(>=2 duplicates)")
  ggsave(file.path(plots_dir, "studio01_pit_collisions.png"),
         p, width = 6, height = 4, dpi = 150)
  invisible(tbl)
}
studio1_problem_3_fish()

# Real-world note: a Snake River brood year tags ~10^6 smolts from a
# 10^14-bin pool, so collisions are effectively impossible for typical
# subsampling. The same math applies to GSI baseline drops (two fish with
# the same multilocus genotype by chance), where the "tag pool" is the
# number of distinguishable multilocus genotypes.
