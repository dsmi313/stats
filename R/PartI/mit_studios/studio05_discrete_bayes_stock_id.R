# Studio 5 - Discrete Bayesian updating (stock identification from markers)
# ---------------------------------------------------------------
# Source: https://ocw.mit.edu/courses/18-05-introduction-to-probability-and-statistics-spring-2022/mit18_05_s22_studio5-instructions.pdf
#
# MIT Studio 5 uses five Platonic dice (4, 6, 8, 12, 20 sided) and updates a
# prior over the type of die from each roll. Our fish version:
#
# Five candidate Snake River stocks. Each fish carries a genetic marker
# whose distribution differs by stock (different "die"). For each marker
# observed we apply Bayes' rule to update the posterior over the stock of
# origin. This is exactly the iterative replacement scheme used for GSI
# uncertainty in Section 10 of the learning plan, with a uniform prior.

set.seed(2026)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# Five candidate stocks ----------------------------------------------------
STOCKS <- c("LOSALM", "CHMBLN", "IMNAHA", "TUCANO", "GRROND")
# Per-stock allele frequencies at a 4-allele marker.  Rows = stocks, cols =
# alleles. Each row sums to 1. Substitute real GSI baseline frequencies here.
LIK <- matrix(c(
  0.50, 0.25, 0.15, 0.10,
  0.40, 0.40, 0.10, 0.10,
  0.20, 0.30, 0.30, 0.20,
  0.10, 0.20, 0.40, 0.30,
  0.05, 0.10, 0.25, 0.60
), nrow = 5, byrow = TRUE,
   dimnames = list(STOCKS, paste0("A", 1:4)))

# Problem 0a/0b/0c: enumerate hypotheses, outcomes, likelihoods ------------
studio5_problem_0_fish <- function() {
  cat("Hypotheses (candidate stocks):\n");      print(STOCKS)
  cat("Possible outcomes (alleles):\n");        print(colnames(LIK))
  cat("Likelihood table P(allele | stock):\n"); print(round(LIK, 3))
  cat("\n")
}
studio5_problem_0_fish()

# Problem 1b: Bayesian updating from nrolls allele observations -----------
studio5_problem_1b_fish <- function(prior, n_alleles_observed,
                                    plot_individual = FALSE,
                                    true_stock = "IMNAHA") {
  true_freq <- LIK[true_stock, ]
  alleles <- sample(colnames(LIK), size = n_alleles_observed,
                    replace = TRUE, prob = true_freq)
  posteriors <- matrix(0, nrow = n_alleles_observed + 1,
                       ncol = length(STOCKS),
                       dimnames = list(NULL, STOCKS))
  posteriors[1, ] <- prior
  for (i in seq_len(n_alleles_observed)) {
    a <- alleles[i]
    lik <- LIK[, a]
    new <- posteriors[i, ] * lik
    posteriors[i + 1, ] <- new / sum(new)
  }
  png(file.path(plots_dir, "studio05_stacked_posteriors.png"),
      width = 800, height = 400)
  barplot(t(posteriors), beside = FALSE,
          col = c("steelblue", "tan", "darkseagreen",
                  "orchid", "firebrick"),
          legend.text = STOCKS,
          args.legend = list(x = "topright", bty = "n"),
          xlab = "marker number", ylab = "posterior P(stock | data)",
          main = paste0("Studio 5 - Posterior over stock (true = ",
                        true_stock, ")"))
  dev.off()
  invisible(posteriors)
}

# Problem 1c: compare uniform vs strongly informative prior ----------------
studio5_problem_1c_fish <- function() {
  cat("Run with uniform prior:\n")
  uniform <- rep(1/5, 5); names(uniform) <- STOCKS
  post_unif <- studio5_problem_1b_fish(prior = uniform,
                                       n_alleles_observed = 20,
                                       true_stock = "IMNAHA")
  print(round(post_unif[nrow(post_unif), ], 3))

  cat("\nRun with prior concentrated on LOSALM (0.001 elsewhere, 0.996 LOSALM):\n")
  biased <- c(LOSALM = 0.996, CHMBLN = 0.001, IMNAHA = 0.001,
              TUCANO = 0.001, GRROND = 0.001)
  post_biased <- studio5_problem_1b_fish(prior = biased,
                                         n_alleles_observed = 20,
                                         true_stock = "IMNAHA")
  print(round(post_biased[nrow(post_biased), ], 3))

  cat("\nA strongly informative wrong prior takes many more markers to\n")
  cat("be overwhelmed by data. This is exactly why partial-pooling priors\n")
  cat("regularize sparse-stock estimates more gracefully than bootstrap\n")
  cat("resampling does (Plan, Sections 6 and 19).\n\n")
}
studio5_problem_1c_fish()

# Problem 1d: pathological prior - zero on the true stock -----------------
studio5_problem_1d_fish <- function() {
  cat("Prior that excludes the true stock (IMNAHA gets zero):\n")
  no_imnaha <- c(LOSALM = 0.25, CHMBLN = 0.25, IMNAHA = 0.0,
                 TUCANO = 0.25, GRROND = 0.25)
  post <- studio5_problem_1b_fish(prior = no_imnaha,
                                  n_alleles_observed = 30,
                                  true_stock = "IMNAHA")
  print(round(post[nrow(post), ], 3))
  cat("Posterior on IMNAHA stays exactly zero no matter how much evidence\n")
  cat("arrives. A prior of zero is a permanent veto - Cromwell's rule.\n\n")
}
studio5_problem_1d_fish()
