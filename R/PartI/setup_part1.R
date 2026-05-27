# Part I setup helper
# Usage:
#   source("R/PartI/setup_part1.R")
#
# Installs missing packages required by Part I and reports versions.

required_packages <- c("ggplot2", "dplyr", "tibble")

cat("\nPart I setup starting...\n")
cat(sprintf("R version: %s\n", R.version.string))

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("Installing missing packages:\n")
  cat(sprintf("  - %s\n", paste(missing_packages, collapse = ", ")))
  install.packages(missing_packages)
} else {
  cat("All required packages are already installed.\n")
}

cat("\nLoaded package versions:\n")
for (pkg in required_packages) {
  version <- as.character(utils::packageVersion(pkg))
  cat(sprintf("  - %s: %s\n", pkg, version))
}

cat("\nPart I setup complete.\n")
