# Session 8: composition two ways, accounting versus likelihood
# Objective: explain why splitting an expanded total into origin groups can be
# done by bookkeeping or by a likelihood, why they agree on the point estimate,
# and why only the likelihood tells you honestly how much an imperfect PBT tag
# rate costs your certainty. This is the sharpest adult/smolt divergence.

# Truth: one stratum. Known origin composition (two hatchery groups A and B plus
# wild) and known PBT tag rates. Tag rates are imperfect: some hatchery fish carry
# no PBT tag and land, undetected, in the Unassigned pile with the wild fish.
truth <- list(pA = 0.45, pB = 0.35, pW = 0.20, N = 400L, trA = 0.45, trB = 0.30)
set.seed(8)

# Simulate observed PBT calls. A hatchery fish of group g is PBT-assigned with
# probability tagRate_g; otherwise it is Unassigned. Every wild fish is Unassigned.
simulate <- function(tr) {
  o <- sample(c("A", "B", "W"), truth$N, replace = TRUE, prob = c(truth$pA, truth$pB, truth$pW))
  a <- ifelse(o == "A" & runif(truth$N) < tr["A"], "A",
       ifelse(o == "B" & runif(truth$N) < tr["B"], "B", "U"))
  c(A = sum(a == "A"), B = sum(a == "B"), U = sum(a == "U"))
}
tr  <- c(A = truth$trA, B = truth$trB)
obs <- simulate(tr)

# Build 1: accounting (escapeLGD PBT_expand_calc; smoltEASE thetahat is this style).
# Expand each tagged group by 1/tagRate, take the extra fish out of the Unassigned
# pile, and call whatever is left of that pile wild. Truncate if the pile runs out.
account <- function(nG, nU, tr) {
  nU <- unname(nU); expand <- nG / tr; d <- expand - nG
  if (nU < sum(d)) { d <- nU * d / sum(d); expand <- nG + d; wild <- 0 }
  else wild <- nU - sum(d)
  v <- c(expand, W = unname(wild)); v / sum(v)
}

# Build 2: likelihood (escapeLGD PBT_expand_calc_MLE). A multinomial over softmax
# proportions of (A, B, wild): an assigned-g fish has probability p_g * tagRate_g,
# and the Unassigned bucket has probability pW + sum((1 - tagRate) * p_g).
softmax <- function(x) { e <- exp(x - max(x)); e / sum(e) }
mle_ll  <- function(p, nG, nU, tr) {
  pg <- p[-length(p)]; pw <- p[length(p)]
  sum(nG * log(pg * tr)) + nU * log(pw + sum((1 - tr) * pg))
}
mle <- function(nG, nU, tr) {
  o <- optim(rep(1, length(nG) + 1), function(par) mle_ll(softmax(par), nG, nU, tr),
             method = "BFGS", control = list(fnscale = -1, maxit = 1000))
  setNames(softmax(o$par), c(names(nG), "W"))
}

ac <- account(obs[c("A", "B")], obs["U"], tr)
ml <- mle(obs[c("A", "B")], obs["U"], tr)
cat(sprintf("truth       A %.3f  B %.3f  W %.3f\n", truth$pA, truth$pB, truth$pW))
cat(sprintf("accounting  A %.3f  B %.3f  W %.3f\n", ac["A"], ac["B"], ac["W"]))
cat(sprintf("likelihood  A %.3f  B %.3f  W %.3f\n", ml["A"], ml["B"], ml["W"]))
cat(sprintf("largest gap between the two point estimates: %.4f\n", max(abs(ac - ml))))

# The point estimates coincide because accounting IS the interior MLE. The real
# divergence is that only the likelihood path carries an interval. Profile the
# multinomial over the wild fraction and read a chi-squared interval off it; run it
# at the real (poor) tag rates and again at near-perfect tags on the same truth, to
# show the interval widens exactly because the tag rate is imperfect.
profile_wild <- function(nG, nU, tr) {
  grid <- seq(0.001, 0.6, length.out = 500)
  ll <- sapply(grid, function(pw) {
    o <- optim(rep(0, length(nG)), function(par) {
      pg <- softmax(par) * (1 - pw); mle_ll(c(pg, pw), nG, nU, tr) },
      method = "BFGS", control = list(fnscale = -1, maxit = 500))
    o$value })
  ci <- range(grid[ll >= max(ll) - qchisq(0.95, 1) / 2])
  list(grid = grid, ll = ll, ci = ci, mle = grid[which.max(ll)])
}
set.seed(8); obs_good <- simulate(c(A = 0.95, B = 0.95))
pp <- profile_wild(obs[c("A", "B")], obs["U"], tr)               # poor tags (real)
pg <- profile_wild(obs_good[c("A", "B")], obs_good["U"], c(A = 0.95, B = 0.95))  # good tags
cat(sprintf("wild-fraction 95%% interval, poor tags [%.3f, %.3f] width %.3f\n",
            pp$ci[1], pp$ci[2], diff(pp$ci)))
cat(sprintf("wild-fraction 95%% interval, good tags [%.3f, %.3f] width %.3f\n",
            pg$ci[1], pg$ci[2], diff(pg$ci)))

png("figs/session08_composition.png", width = 900, height = 600)
plot(pg$grid, pg$ll - max(pg$ll), type = "l", lwd = 3, col = "black", ylim = c(-6, 0.3),
     xlim = c(0, 0.5), xlab = "wild fraction", ylab = "relative profile log-likelihood",
     main = "Session 8: the likelihood path prices imperfect tag rates")
lines(pp$grid, pp$ll - max(pp$ll), lwd = 3, col = "darkorange")
abline(h = -qchisq(0.95, 1) / 2, col = "steelblue", lwd = 2, lty = 3)
abline(v = truth$pW, col = "firebrick", lwd = 2)
abline(v = ac["W"], col = "grey50", lwd = 2, lty = 2)  # accounting: a point, no interval
legend("topright", c("good tags", "poor tags (real)", "chi-sq cutoff", "truth", "accounting point"),
       lwd = c(3, 3, 2, 2, 2), lty = c(1, 1, 3, 1, 2),
       col = c("black", "darkorange", "steelblue", "firebrick", "grey50"), bty = "n")
dev.off()

# Exercise. Drop trA and trB to 0.15 and rerun. The point estimates still roughly
# track truth, but the poor-tags profile goes nearly flat and its interval runs
# most of the axis; say out loud why a single accounting number would hide that.

# Locate.
# escapeLGD (adults): HNC_expand_unkGSI() in R/wrappers_HNC_expand.R takes
#   method = c("Account", "MLE") and branches at lines 85-95 to HNC_expand_one_strat()
#   or HNC_expand_one_strat_MLE(). The MLE path runs PBT_expand_calc_MLE()
#   (R/composition_estimation_utils.R:28): this exact multinomial over softmax
#   proportions with an analytic gradient (PBT_grad), optim BFGS. The accounting
#   path is PBT_expand_calc() at line 440: expand by 1/tagRate, borrow from
#   Unassigned, truncate, normalize, precisely the account() here.
# smoltEASE (smolts): has no likelihood path. thetahat() in R/SCRAPI2.R:173 does
#   inverse-sample-rate weighting and prop.table only, the accounting side of this
#   divergence with no MLE counterpart. That absence is the sharpest place the two
#   production tools part ways.

writeLines(c(
"How I would explain session 8 in three minutes",
"",
"Once we know the total run, we split it into origin groups using the PBT tags fish",
"carry. The catch is that tags are imperfect: a hatchery fish whose tag failed looks",
"exactly like a wild fish. There are two ways to undo that. The bookkeeping way",
"expands each tagged group by one over its tag rate and calls whatever is left over",
"wild. The likelihood way writes down the chance of every observed call and finds",
"the composition that makes the data most likely.",
"",
sprintf("On the point estimate they agree here to within %.3f, because the bookkeeping is", max(abs(ac - ml))),
"the likelihood answer whenever nothing hits a boundary. The difference is honesty",
sprintf("about certainty. With the real, poor tag rates the likelihood's wild-fraction"),
sprintf("interval is %.2f wide; give the same fish near-perfect tags and it shrinks to %.2f.",
        diff(pp$ci), diff(pg$ci)),
"The single accounting number cannot tell you that a low tag rate has made the split",
"much less certain. Adults have this likelihood path; smolts only keep the books."
), "docs/session08_explain.md")
