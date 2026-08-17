# Session 4: profile versus marginalize a nuisance parameter
# Objective: explain the difference between profiling a nuisance parameter out
# and integrating it out, and show they give nearly the same interval here.

# This is escapeLGD's fallback likelihood exactly. Parameter of interest is
# P(fallback); the nuisance is P(reascend | fallback). df, dfr are spillway fish
# that fell back and of those later reascended; dt, dr are ladder ascensions and
# of those the ones that were reascensions.
truth <- list(p_fallback = 0.08, p_reascend = 0.85, n_spill = 300L, n_ladder = 5000L)
set.seed(4)

dfr <- rbinom(1, truth$n_spill, truth$p_reascend)
dr  <- rbinom(1, truth$n_ladder, truth$p_fallback * truth$p_reascend)
df <- truth$n_spill; dt <- truth$n_ladder

# Build the two-parameter joint log-likelihood by hand (the production form).
loglik <- function(pf, pre) dbinom(dfr, df, pre, log = TRUE) +
                            dbinom(dr, dt, pf * pre, log = TRUE)

pf_grid  <- seq(0.02, 0.20, length.out = 300)
pre_grid <- seq(0.50, 0.999, length.out = 300)
surface  <- outer(pf_grid, pre_grid, Vectorize(loglik))

# Profile: at each P(fallback), keep the best nuisance value.
profile <- apply(surface, 1, max)
# Marginalize: at each P(fallback), integrate the nuisance out (flat prior).
dpre     <- diff(pre_grid)[1]
marginal <- log(rowSums(exp(surface - max(surface))) * dpre) + max(surface)

# Recover. Read a likelihood-ratio interval off each curve (drop of 1.92 = half
# the chi-squared(1) cutoff) and show they nearly coincide.
peak_pf <- pf_grid[which.max(profile)]
prof0 <- profile - max(profile); marg0 <- marginal - max(marginal)
cut <- -qchisq(0.95, 1) / 2
ci_prof <- range(pf_grid[prof0 >= cut]); ci_marg <- range(pf_grid[marg0 >= cut])
cat(sprintf("truth %.3f | MLE %.3f\n", truth$p_fallback, peak_pf))
cat(sprintf("profile   interval [%.3f, %.3f]\n", ci_prof[1], ci_prof[2]))
cat(sprintf("marginal  interval [%.3f, %.3f]\n", ci_marg[1], ci_marg[2]))

png("figs/session04_profile_vs_marginalize.png", width = 900, height = 600)
plot(pf_grid, prof0, type = "l", lwd = 3, col = "black", ylim = c(cut - 2, 0.3),
     xlab = "P(fallback)", ylab = "relative log-likelihood (peak at 0)",
     main = "Session 4: profiling vs marginalizing the nuisance parameter")
lines(pf_grid, marg0, lwd = 3, col = "darkorange", lty = 2)
abline(h = cut, col = "steelblue", lwd = 2, lty = 3)
abline(v = truth$p_fallback, col = "firebrick", lwd = 3)
legend("topright", c("profile", "marginal", "chi-sq cutoff"), bty = "n",
       lwd = 3, lty = c(1, 2, 3), col = c("black", "darkorange", "steelblue"))
dev.off()

# Exercise. Shrink n_spill to 30 so the nuisance is poorly pinned down. The
# marginal curve widens more than the profile; say out loud why integrating over
# a badly-known nuisance costs more width than profiling it away.

# Locate.
# escapeLGD (adults): fallback_log_likelihood() and its analytic
#   gradient_fallback_log_likelihood() in R/fallback_reascend_likelihood.R are this
#   exact two-parameter likelihood; nightFall() optimises it with optim and then
#   bootstraps, which is Monte Carlo marginalization rather than the grid integral.
# smoltEASE (smolts): has no two-parameter reascension likelihood at all, because
#   smolts do not fall back and reascend. That absence is itself a talk point.

endgap <- max(abs(ci_prof - ci_marg))
same <- if (endgap < 0.005) "almost the same interval" else "visibly different intervals"
writeLines(c(
"How I would explain session 4 in three minutes",
"",
"To estimate how many adults fall back over the dam, we also have to reckon with",
"a second unknown we do not really care about: of those that fall back, how many",
"climb the ladder again. That second number is a nuisance, and there are two",
"honest ways to deal with it.",
"",
"Profiling asks, for each fallback rate, the best-case nuisance value, and reads",
"the answer off that best-case ridge. Marginalizing instead averages over every",
sprintf("plausible nuisance value. Here profiling gives [%.3f, %.3f] and marginalizing",
        ci_prof[1], ci_prof[2]),
sprintf("[%.3f, %.3f], their endpoints within %.3f of each other, so the two give %s",
        ci_marg[1], ci_marg[2], endgap, same),
"and it does not matter which we quote. They only pull apart when the nuisance",
"itself is barely known, and then averaging over our ignorance costs extra width."
), "docs/session04_explain.md")
