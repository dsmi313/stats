# Studio 3 - Histograms (fish ladder inter-arrival times)
# ---------------------------------------------------------------
# Source: https://ocw.mit.edu/courses/18-05-introduction-to-probability-and-statistics-spring-2022/mit18_05_s22_studio3-instructions.pdf
#
# MIT Studio 3 plots histograms of exponential samples and their averages,
# illustrating the Central Limit Theorem. We use times between successive
# fish passing the LWG ladder window (modeled Poisson, so inter-arrival
# times are Exponential).
#
# Wrapper functions mirror MIT's structure.

set.seed(2026)

plots_dir <- file.path("docs", "figures", "PartI")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# Problem 1a: frequency histogram of inter-arrival times -------------------
studio3_problem_1a_fish <- function(rate, nsamples) {
  data <- rexp(nsamples, rate = rate)
  bin_width <- 0.2
  bins <- seq(0, max(data) + bin_width, by = bin_width)
  png(file.path(plots_dir, "studio03_freq_exp.png"),
      width = 600, height = 400)
  hist(data, breaks = bins, freq = TRUE,
       col = "lightyellow", border = "steelblue",
       main = "Studio 3 - Frequency histogram of fish inter-arrival times",
       xlab = "minutes between successive fish at LWG ladder")
  dev.off()
  invisible(data)
}
studio3_problem_1a_fish(rate = 1.0, nsamples = 1000)

# Problem 1b: density histogram + true pdf overlay -------------------------
studio3_problem_1b_fish <- function(rate, nsamples) {
  data <- rexp(nsamples, rate = rate)
  bin_width <- 0.2
  bins <- seq(0, max(data) + bin_width, by = bin_width)
  png(file.path(plots_dir, "studio03_density_exp.png"),
      width = 600, height = 400)
  hist(data, breaks = bins, freq = FALSE,
       col = "lightyellow", border = "steelblue",
       main = "Studio 3 - Density histogram + Exponential pdf",
       xlab = "minutes between successive fish at LWG ladder")
  xs <- seq(0, max(data), length.out = 200)
  lines(xs, dexp(xs, rate = rate), col = "firebrick", lwd = 2)
  dev.off()
  invisible(data)
}
studio3_problem_1b_fish(rate = 1.0, nsamples = 2000)

# Problem 2a: density histogram of the average of 2 inter-arrival samples
studio3_problem_2a_fish <- function(rate, nsamples) {
  x1 <- rexp(nsamples, rate = rate)
  x2 <- rexp(nsamples, rate = rate)
  y  <- (x1 + x2) / 2
  bin_width <- 0.15
  bins <- seq(0, max(y) + bin_width, by = bin_width)
  png(file.path(plots_dir, "studio03_avg_n2.png"),
      width = 600, height = 400)
  hist(y, breaks = bins, freq = FALSE,
       col = "lightyellow", border = "steelblue",
       main = "Studio 3 - Average of 2 inter-arrival times",
       xlab = "(x1 + x2) / 2")
  dev.off()
  invisible(y)
}
studio3_problem_2a_fish(rate = 1.0, nsamples = 2000)

# Problem 2b: CLT - average of n_to_average exponentials + Normal overlay --
studio3_problem_2b_fish <- function(rate, nsamples,
                                    n_to_average, bin_width) {
  mat <- matrix(rexp(nsamples * n_to_average, rate = rate),
                ncol = n_to_average)
  y   <- rowMeans(mat)
  mu  <- 1 / rate
  sd_y <- (1 / rate) / sqrt(n_to_average)
  bins <- seq(max(0, mu - 4 * sd_y),
              mu + 4 * sd_y + bin_width, by = bin_width)
  bins <- c(bins, max(y) + bin_width)
  png(file.path(plots_dir,
                paste0("studio03_clt_n", n_to_average, ".png")),
      width = 600, height = 400)
  hist(y, breaks = sort(unique(bins)), freq = FALSE,
       col = "lightyellow", border = "steelblue",
       main = paste0("Studio 3 - CLT: average of ",
                     n_to_average, " inter-arrival times"),
       xlab = "average waiting time (min)")
  xs <- seq(mu - 4 * sd_y, mu + 4 * sd_y, length.out = 300)
  lines(xs, dnorm(xs, mean = mu, sd = sd_y),
        col = "firebrick", lwd = 2)
  dev.off()
  invisible(y)
}

# Run the CLT sweep
for (k in c(2, 9, 36, 100)) {
  studio3_problem_2b_fish(rate = 1.0, nsamples = 5000,
                          n_to_average = k, bin_width = 0.05)
}

cat("Studio 3 outputs saved to", plots_dir, "\n")
cat("As n_to_average grows the density of averaged inter-arrival times\n")
cat("converges to the Normal pdf with sd shrinking like 1/sqrt(n).\n")
