# Session 7: the real run is the same arithmetic
# Objective: show that on real MY2025 inputs the smolt machine is the session-1
# move and nothing more. Expand a real day by the plug-in formula, sum its week,
# and get exactly the weekly number SCRAPI2 prints, because it is the identical
# calculation surfaced, not an independent check of it.

# No simulation. This session reads the real MY2025 steelhead files in data/ and
# runs production smoltEASE::SCRAPI2() on them. The passage file dates are
# %m/%d/%Y and the trap file dates are %Y-%m-%d; SCRAPI2's tryFormats parses both,
# so we do not touch them. Guidance efficiency is the fixed column in the passage
# file (no geDraws), so the expansion is deterministic.
if (!requireNamespace("smoltEASE", quietly = TRUE))
  stop("install smoltEASE (github.com/dsmi313/smoltEASE) to run this session")

pass  <- read.csv("data/MY2025_STHD_passage.csv")
ntrap <- nrow(read.csv("data/MY2025_STHD_trapData.csv"))
cat(sprintf("real inputs: %d passage days, %d trapped fish\n", nrow(pass), ntrap))

# Build. The daily expansion is exactly session 1: count over the fraction seen,
# where the fraction seen is SampleRate times GuidanceEfficiency. This is line for
# line what SCRAPI2 does at SCRAPI2.R:306, pass$estimated <- SampleCount / true.
hand_day  <- pass$SampleCount / (pass$SampleRate * pass$GuidanceEfficiency)
sel       <- which(pass$SampleEndDate == "05/12/2025")  # a mid-run day near the peak
wsel      <- pass$Week[sel]
wk_rows   <- which(pass$Week == wsel)
hand_one  <- hand_day[sel]
hand_week <- sum(hand_day[wk_rows])
hand_tot  <- sum(hand_day)
cat(sprintf("chosen day %s (Week %d): %d counted / (%.2f * %.3f) = %.0f fish\n",
            pass$SampleEndDate[sel], wsel, pass$SampleCount[sel], pass$SampleRate[sel],
            pass$GuidanceEfficiency[sel], hand_one))
cat(sprintf("hand sum over Week %d (%d days): %.0f\n", wsel, length(wk_rows), hand_week))

# Recover. Run SCRAPI2 on the real files. B is dropped far below the 5000 default
# purely for runtime; a smaller B only widens the bootstrap CI, it does not move
# the deterministic point expansion checked here.
# SCRAPI2 returns list(CI, bootstrap) and neither element carries the per-week
# totals, so the only route to its weekly number is parsing the console output.
# This parsing is therefore coupled to SCRAPI2's print format: the by-week header
# and its two-line print at SCRAPI2.R:325, and the grand total at line 327.
options(width = 10000)
invisible(capture.output(
  res <- smoltEASE::SCRAPI2(smoltData   = "data/MY2025_STHD_trapData.csv",
                            passageData = "data/MY2025_STHD_passage.csv",
                            Primary = "GenStock", RTYPE = "W", alph = 0.1,
                            B = 200, Run = file.path(tempdir(), "scrapi_"), seed = 1)
) -> out)
gt <- as.numeric(sub("^Total smolts:\\s*", "", grep("^Total smolts:", out, value = TRUE)))
h  <- grep("Total smolts by week:", out)
wk <- setNames(scan(text = out[h + 2], quiet = TRUE), scan(text = out[h + 1], quiet = TRUE))
cat(sprintf("SCRAPI2 Week %d total: %d   hand sum: %.0f\n", wsel, wk[as.character(wsel)], hand_week))
cat(sprintf("SCRAPI2 total smolts: %d   sum of hand daily expansions: %.0f\n", gt, hand_tot))
stopifnot(wk[as.character(wsel)] == round(hand_week), abs(gt - hand_tot) <= 1)

png("figs/session07_real_data_scrapi.png", width = 900, height = 600)
cols <- rep("grey75", nrow(pass)); cols[wk_rows] <- "steelblue"
plot(hand_day, type = "h", lwd = 4, col = cols, xlab = "passage day (MY2025)",
     ylab = "hand-expanded fish",
     main = "Session 7: real daily expansions; one week summed to the SCRAPI2 total")
points(sel, hand_one, pch = 19, col = "firebrick")
legend("topright", bty = "n", col = c("steelblue", "firebrick"), pch = c(15, 19),
       legend = c(sprintf("Week %d days: sum %.0f = SCRAPI2's %d", wsel, hand_week, wk[as.character(wsel)]),
                  sprintf("chosen day %s = %.0f", pass$SampleEndDate[sel], hand_one)))
dev.off()

# Exercise. Swap RTYPE = "W" for RTYPE = "HNC" and rerun; the same daily expansion
# feeds a different composition split. Say out loud why the expansion step is
# shared while only the composition step downstream of it changes.

# Locate.
# smoltEASE (smolts): SCRAPI2() in R/SCRAPI2.R forms pass$estimated <- SampleCount /
#   (SampleRate * GuidanceEfficiency) at line 306, the exact hand_day formula, and
#   prints round(tapply(pass$estimated, Week, sum)) at line 325. Matching the hand
#   sum to that print is re-running the same arithmetic, not a second estimator.
# escapeLGD (adults): the analogue is expand_wc_binom_night(), round(wc / wc_prop)
#   summed over statistical weeks; same divide-by-the-fraction-seen move on the
#   adult window count, with nighttime passage standing in for guidance efficiency.

writeLines(c(
"How I would explain session 7 in three minutes",
"",
"Every number in the talk so far was simulated. This one is not: it is the real",
sprintf("MY2025 steelhead run, %d days of trap counts and %d fish. Be clear about what this", nrow(pass), ntrap),
"shows and what it does not: we did not independently verify SCRAPI2. We cannot, by",
"hand. What we can show is that the weekly number SCRAPI2 prints is nothing more",
"than the plug-in expansion, summed.",
"",
sprintf("Take %s: %d fish counted at a %.0f%% sample rate and %.0f%% guidance is that count",
        pass$SampleEndDate[sel], pass$SampleCount[sel], 100 * pass$SampleRate[sel],
        100 * pass$GuidanceEfficiency[sel]),
sprintf("over the fraction seen, about %.0f fish. Its week has %d days; add their expansions",
        hand_one, length(wk_rows)),
sprintf("and you get %.0f. SCRAPI2 prints exactly %d for Week %d, because line for line it is",
        hand_week, wk[as.character(wsel)], wsel),
"the same division summed the same way: count over sample rate times guidance,",
"summed by week and rounded (SCRAPI2.R line 306, printed at line 325). There is no",
"separate model doing the counting and no magic in the machine; the whole run is",
"this one division, repeated day after day and then grouped."
), "docs/session07_explain.md")
