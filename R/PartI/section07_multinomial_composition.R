#---------------------------------------------------------
# File:   section07_multinomial_composition.R
# Part I, Section 7 — Multinomial composition
#
# Generalises the binomial (two outcomes) to multinomial (many).
# Problem 7c builds the two-sample joint-likelihood that is the
# structural core of the Delomas & Hess estimator.
#
# Repo pointer (escapeLGD/R/composition_estimation_utils.R):
#   lines 62-66:   PBT_log_likelihood uses dmultinom() on the same
#                  (tagged, untagged) two-sample structure you build in 7c
#---------------------------------------------------------
# Section 7 ----

# softmax() helper -- needed by problems 7b and 7c.
# The solutions file defines it at the top level.  Define it here before
# sourcing, or uncomment the line below and paste it inside each function.
#
# softmax <- function(x) { e <- exp(x - max(x)); e / sum(e) }

# Problem 7a: simulate three stocks with rmultinom(); verify MLE = observed proportions
section7_problem_7a_fish <- function(n = 200, props = c(0.6, 0.3, 0.1)) {
  cat("\n----------------------------------\n")
  cat("Problem 7a: rmultinom — MLE is observed proportions\n")

  # Simulate one draw of n fish from three stocks using rmultinom().
  # Compute sample proportions (the MLE) and compare to props.
  # Replicate 1,000 times with replicate(). Plot the distribution of
  # the stock-1 estimate with hist() and mark the truth with abline().

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 7b: multinomial log-likelihood and MLE via optim()
section7_problem_7b_fish <- function(counts = c(120L, 60L, 20L)) {
  cat("\n----------------------------------\n")
  cat("Problem 7b: multinomial log-likelihood; maximize with optim()\n")

  # Write ll(props | counts) = sum(counts * log(props)).
  # Use a softmax reparametrization so optim() works over unconstrained space.
  # Verify the MLE equals counts / sum(counts).

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 7c: joint log-likelihood for tagged and untagged samples
section7_problem_7c_fish <- function(
    tagged_counts   = c(S1 = 42L, S2 = 35L, S3 = 18L),
    untagged_counts = c(S1 = 80L, S2 = 60L, S3 = 30L, Unassigned = 32L)) {

  cat("\n----------------------------------\n")
  cat("Problem 7c: product of two multinomials — skeleton of Delomas & Hess\n")

  # Tagged and untagged fish are separate multinomial samples that share
  # the same underlying stock proportions.
  # Joint log-likelihood = ll(tagged) + ll(untagged).
  # Maximize over the shared stock proportions with optim().

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}


# Problem 7d: bootstrapped CI for stock proportions
section7_problem_7d_fish <- function(counts = c(S1 = 120L, S2 = 60L, S3 = 20L),
                                      boots  = 10000L) {
  cat("\n----------------------------------\n")
  cat("Problem 7d: nonparametric bootstrap CI for stock proportions\n")

  # Build the individual-fish vector from counts.
  # sample(..., replace = TRUE) boots times; compute proportions each draw.
  # Return quantile(., c(0.025, 0.975)) for each stock.

  # Do not change the above code.
  # ********* YOUR CODE HERE ***********

}
