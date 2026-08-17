# Session 7: the real run, checked by hand
# Objective: explain that the whole smolt machine is, at bottom, the session-1
# move applied to real MY2025 inputs: expand one day by hand and watch production
# reproduce it, then watch those daily expansions sum to the number SCRAPI2 prints.

# No simulation. This session reads the real MY2025 steelhead files in data/ and
# runs production smoltEASE::SCRAPI2() on them. The passage file dates are
# %m/%d/%Y and the trap file dates are %Y-%m-%d; SCRAPI2's tryFormats parses both,
# so we do not touch them. Guidance efficiency is the fixed column in the passage
# file (no geDraws), so the point expansion is deterministic and hand-checkable.
if (!requireNamespace("smoltEASE", quietly = TRUE))
  stop("install smoltEASE (github.com/dsmi313/smoltEASE) to run this session")

pass <- read.csv("data/MY2025_STHD_passage.csv")
cat(sprintf("real inputs: %d passage days, %d trapped fish\n",
            nrow(pass), nrow(read.csv("data/MY2025_STHD_trapData.csv"))))

# Build. The daily expansion is exactly session 1: count over the fraction seen,
# where the fraction seen is SampleRate times GuidanceEfficiency.
hand_day <- pass$SampleCount / (pass$SampleRate * pass$GuidanceEfficiency)
sel      <- which(pass$Week == 26)   # Week 26 is a single day, so it is its own week
hand_one <- hand_day[sel]
hand_tot <- sum(hand_day)
cat(sprintf("chosen day %s: %d counted / (%.2f * %.3f) = %.2f fish\n",
            pass$SampleEndDate[sel], pass$SampleCount[sel], pass$SampleRate[sel],
            pass$GuidanceEfficiency[sel], hand_one))

# Recover. Run SCRAPI2 on the real files. B is dropped far below the 5000 default
# purely for runtime (this keeps the session under a minute); a smaller B only
# widens and roughens the bootstrap CI, it does not move the point expansion we
# check here, which is deterministic given fixed GE.
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
cat(sprintf("SCRAPI2 for that day (its Week-26 total): %d   hand: %.0f\n", wk["26"], hand_one))
cat(sprintf("SCRAPI2 total smolts: %d   sum of hand daily expansions: %.0f\n", gt, hand_tot))
stopifnot(wk["26"] == round(hand_one), abs(gt - hand_tot) <= 1)

png("figs/session07_real_data_scrapi.png", width = 900, height = 600)
plot(hand_day, type = "h", lwd = 3, col = "grey70", xlab = "passage day (MY2025)",
     ylab = "hand-expanded fish", main = "Session 7: real daily expansions, summing to the SCRAPI2 total")
points(sel, hand_one, pch = 19, col = "firebrick")
legend("topright", sprintf("total = %.0f = SCRAPI2's %d", hand_tot, gt), bty = "n")
dev.off()

# Exercise. Swap Primary = "GenStock" for RTYPE = "HNC" and rerun; the same daily
# expansion feeds a different composition split. Say out loud why the expansion
# step is shared while only the composition step downstream of it changes.

# Locate.
# smoltEASE (smolts): SCRAPI2() in R/SCRAPI2.R forms pass$estimated <- SampleCount /
#   (SampleRate * GuidanceEfficiency) and sums it by week and stratum; that is the
#   exact hand_day formula above, run over all 61 days. Everything else SCRAPI2 does
#   (composition, wild split, bootstrap CI) sits on top of this one expansion.
# escapeLGD (adults): the analogue is expand_wc_binom_night(), round(wc / wc_prop)
#   summed over statistical weeks; same divide-by-the-fraction-seen move on the
#   adult window count, with nighttime passage standing in for guidance efficiency.

writeLines(c(
"How I would explain session 7 in three minutes",
"",
sprintf("Every number in the talk so far was simulated. This one is not: it is the real"),
sprintf("MY2025 steelhead run, %d days of trap counts and %d fish. Production SCRAPI2 is a",
        nrow(pass), nrow(read.csv("data/MY2025_STHD_trapData.csv"))),
"large program, but on any single day it does exactly what session 1 did: take the",
sprintf("count and divide by the fraction seen. On %s that is %d fish counted at a %.0f%%",
        pass$SampleEndDate[sel], pass$SampleCount[sel], 100 * pass$SampleRate[sel]),
sprintf("sample rate and %.0f%% guidance, giving %.0f fish, and SCRAPI2 reports the same %d.",
        100 * pass$GuidanceEfficiency[sel], hand_one, wk["26"]),
"",
sprintf("Do that for every day and the daily expansions sum to %.0f, which is exactly the", hand_tot),
sprintf("total smolts SCRAPI2 prints, %d. The whole production run is that one move applied", gt),
"day after day and then split into groups. Anyone in the room can check the day I",
"picked with a calculator, and that is the point: the machine has no magic in it."
), "docs/session07_explain.md")
