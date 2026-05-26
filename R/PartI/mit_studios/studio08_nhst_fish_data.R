# Studio 8 - NHST: F-stat, z-test, chi-square, ANOVA on fish data
# ---------------------------------------------------------------
# Source: https://ocw.mit.edu/courses/18-05-introduction-to-probability-and-statistics-spring-2022/mit18_05_s22_studio8-instructions.pdf
#
# Four fish-themed NHST problems following MIT Studio 8:
#   1. F-statistic distribution across stock groups of fish lengths
#   2. z-test on whether mean smolt length equals a nominal target
#   3. chi-square test of independence: stock origin vs trap recapture
#   4. One-way ANOVA across three hatchery treatments

set.seed(2026)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# Problem 1: simulated f-statistic under H0 ------------------------------
studio8_problem_1_fish <- function(n, m, mu, sigma, n_trials) {
  fstats <- numeric(n_trials)
  for (t in seq_len(n_trials)) {
    groups <- matrix(rnorm(n * m, mean = mu, sd = sigma), nrow = m)
    grp_means <- colMeans(groups)
    grand_mean <- mean(groups)
    SSB <- m * sum((grp_means - grand_mean)^2)
    SSW <- sum(sweep(groups, 2, grp_means)^2)
    MSB <- SSB / (n - 1)
    MSW <- SSW / (n * (m - 1))
    fstats[t] <- MSB / MSW
  }
  png(file.path(plots_dir, "studio08_fstat.png"),
      width = 600, height = 400)
  hist(fstats, breaks = 40, freq = FALSE,
       col = "lightyellow", border = "steelblue",
       main = "Studio 8 - F-statistic under H0 (equal stock means)",
       xlab = "F statistic")
  xs <- seq(0, max(fstats), length.out = 300)
  lines(xs, df(xs, df1 = n - 1, df2 = n * (m - 1)),
        col = "firebrick", lwd = 2)
  dev.off()
  invisible(fstats)
}
studio8_problem_1_fish(n = 5, m = 8, mu = 100, sigma = 12, n_trials = 5000)

# Problem 2: z-test on mean smolt length ---------------------------------
studio8_problem_2_fish <- function(data, mu0, known_sigma, alpha) {
  xbar <- mean(data)
  z    <- (xbar - mu0) / (known_sigma / sqrt(length(data)))
  p    <- 2 * (1 - pnorm(abs(z)))
  cat("Studio 8 Problem 2: z-test on smolt length\n")
  cat("  n =", length(data),
      "  xbar =", round(xbar, 2),
      "  mu0 =", mu0,
      "  sigma =", known_sigma, "\n")
  cat("  z =", round(z, 3),
      "   p-value =", round(p, 4),
      "   alpha =", alpha, "\n")
  cat("  decision: ",
      ifelse(p < alpha,
             "REJECT H0 (lengths differ from target)",
             "FAIL TO REJECT H0"), "\n\n")
  invisible(list(z = z, p = p))
}
fake_lengths <- rnorm(40, mean = 98, sd = 10)
studio8_problem_2_fish(fake_lengths, mu0 = 100,
                       known_sigma = 10, alpha = 0.05)

# Problem 3: chi-square independence of stock origin vs trap recap -----
studio8_problem_3_fish <- function(contingency, alpha) {
  cat("Studio 8 Problem 3: chi-square independence test\n")
  cat("  Contingency table (rows = stock, cols = recap status):\n")
  print(contingency)
  fit <- chisq.test(contingency)
  cat("  chi-squared =", round(fit$statistic, 3),
      "   df =", fit$parameter,
      "   p-value =", round(fit$p.value, 4), "\n")
  cat("  decision: ",
      ifelse(fit$p.value < alpha,
             "REJECT independence (stock and recap are associated)",
             "FAIL TO REJECT independence"), "\n\n")
  invisible(fit)
}
sample_contingency <- matrix(c(40, 12,
                               25, 23,
                               18, 22),
                             nrow = 3, byrow = TRUE,
                             dimnames = list(
                               stock = c("LOSALM", "CHMBLN", "IMNAHA"),
                               recap = c("recovered", "not_recovered")))
studio8_problem_3_fish(sample_contingency, alpha = 0.05)

# Problem 4: ANOVA across three hatchery treatments ---------------------
studio8_problem_4_fish <- function(T1, T2, T3, alpha) {
  df <- data.frame(
    treatment = factor(rep(c("Hatch_A", "Hatch_B", "Hatch_C"),
                           c(length(T1), length(T2), length(T3)))),
    length    = c(T1, T2, T3)
  )
  fit  <- aov(length ~ treatment, data = df)
  smry <- summary(fit)
  p    <- smry[[1]][["Pr(>F)"]][1]
  cat("Studio 8 Problem 4: one-way ANOVA on three hatchery groups\n")
  print(smry)
  cat("  p-value =", round(p, 4), "\n")
  cat("  decision: ",
      ifelse(p < alpha,
             "REJECT H0 (hatchery treatment matters)",
             "FAIL TO REJECT H0"), "\n\n")
  invisible(list(fit = fit, p = p))
}
hatch_A <- rnorm(20, mean = 100, sd = 8)
hatch_B <- rnorm(20, mean = 104, sd = 8)
hatch_C <- rnorm(20, mean = 108, sd = 8)
studio8_problem_4_fish(hatch_A, hatch_B, hatch_C, alpha = 0.05)
