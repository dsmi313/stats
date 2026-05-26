# Part I -> SCOBI / escapeLGD line map

After working through each Part I section, the listed source lines in the
production repos should be readable without surprise. Every Part I file
includes a "Repo pointer" comment header citing the exact lines it
mirrors. This table is the consolidated index.

## SCOBI (https://github.com/mackerman44/SCOBI)

| Part I section | What you build | SCOBI source (file:lines) |
| --- | --- | --- |
| 2. MLE for a_d | smolt-trap `AllPrime` table with Strat, PGrp, SR (SCRAPI thetahat input) | `R/SCRAPI.r:284-313` |
| 2. MLE for a_d | `dailypass = Tally / Ptrue` (per-day MLE; SMOLT trap, Ptrue includes GE) | `R/SCRAPI.r:75` |
| 4. Parametric bootstrap | SCRAPI smolt `passdata = data.frame(Stratum, Tally, Ptrue)` | `R/SCRAPI.r:232` |
| 4. Parametric bootstrap | SCRAPI `pass$true = SampleRate * GuidanceEfficiency` (smolt only) | `R/SCRAPI.r:227` |
| 2. MLE for a_d | Inverse-SR weighting `mApply(1/SR, ..., sum)` | `R/SCRAPI.r:78, 94, 103` |
| 2. MLE for a_d | `prop.table(Freqs, margin = ...)` | `R/SCRAPI.r:96, 105` |
| 4. Parametric bootstrap | `cntstar <- rbinom(1, dailypass[i], Ptrue[i])` daily loop | `R/SCRAPI.r:141-145` |
| 4. Parametric bootstrap | `dailyStar <- data.frame(Stratum, Tally, Ptrue)` | `R/SCRAPI.r:145` |
| 5. Nonparametric bootstrap | `theta.b <- matrix(numeric(p*B), ncol = p)` setup | `R/SCRAPI.r:131` |
| 5. Nonparametric bootstrap | Weighted FishWH resample per stratum | `R/SCRAPI.r:148-157` |
| 5. Nonparametric bootstrap | Weighted FishDat resample per stratum | `R/SCRAPI.r:159-169` |
| 5. Nonparametric bootstrap | `theta.b[b, ] <- c(...)` row-store | `R/SCRAPI.r:175-176` |
| 5. Nonparametric bootstrap | The Error 1 (theta.b dimension mismatch) trigger | `R/SCRAPI.r:175` |
| 5. Nonparametric bootstrap | `thetahat_toy()` mirrors the full thetahat function | `R/SCRAPI.r:74-126` |
| 6. Stratification | `Cpattern <- unique(cbind(pass$Week, pass$Collaps))` | `R/SCRAPI.r:217` |
| 6. Stratification | Collaps assignment by date match (Error 3 trigger) | `R/SCRAPI.r:254-259` |
| 6. Stratification | True-rate assignment by date match | `R/SCRAPI.r:262-269` |

## escapeLGD (https://github.com/delomast/escapeLGD)

| Part I section | What you build | escapeLGD source (file:lines) |
| --- | --- | --- |
| 1. Binomial window count (adult) | `wc` tibble with sWeek, wc + expansion by wc_prop (no GE -- adults have one filter) | `R/night_fall_reascend_wc_binom.R:130-150` |
| 3. Variance + delta method | `fallback_log_likelihood(pf, pre_f, dfr, df, dr, dt)` | `R/fallback_reascend_likelihood.R:13-15` |
| 3. Variance + delta method | `gradient_fallback_log_likelihood()` analytical gradient | `R/fallback_reascend_likelihood.R:26-50` |
| 3. Variance + delta method | `optimllh()` softmax wrapper | `R/fallback_reascend_likelihood.R:61-70` |
| 3. Variance + delta method | Joint optim call (`L-BFGS-B`, `fnscale = -1`) | `R/night_fall_reascend_wc_binom.R:71-75` |
| 3. Variance + delta method | Nighttime bootstrap `rbinom(boots, totalPass, p_night)` | `R/night_fall_reascend_wc_binom.R:110` |
| 4. Parametric bootstrap | Vectorized daily-count bootstrap `rbinom(boots, wc, wc_prop)` | `R/night_fall_reascend_wc_binom.R:150` |
| 7. Multinomial composition | `softMax(x)` numerically stable softmax | `R/composition_estimation_utils.R:15-18` |
| 7. Multinomial composition | `PBT_log_likelihood(pGroups, pW, nGroups, nUntag, tagRates)` | `R/composition_estimation_utils.R:62-66` |
| 7. Multinomial composition | `PBT_optimllh()` softmax wrapper | `R/composition_estimation_utils.R:77-84` |
| 7. Multinomial composition | `PBT_expand_calc_MLE()` optim call | `R/composition_estimation_utils.R:38-41` |
| 7. Multinomial composition | `PBT_expand_calc()` accounting/TotEx method | `R/composition_estimation_utils.R:440-455` |
| 7. Multinomial composition | `PBT_breakdown()` nonparametric bootstrap | `R/composition_estimation_utils.R:466-486` |

## How to use this map

1. Clone the production repos somewhere local:

   ```bash
   git clone https://github.com/mackerman44/SCOBI ../SCOBI
   git clone https://github.com/delomast/escapeLGD ../escapeLGD
   ```

2. After finishing each Part I section, open the cited source file in
   parallel with your section script and walk through the cited lines.
   Every variable name, function signature, and loop structure in the
   Part I scripts is chosen to match the production code.

3. Sections that go further than Part I:
   - Joint D&H likelihood (PBT_var2, PBT_var3) -> Part II, PLAN.md
     Sections 8-10
   - Compound bootstrap with composition layer -> Part III, PLAN.md
     Section 12
   - lgr2SCRAPI data pipeline -> Part III, PLAN.md Section 13

## Naming key

| Symbol used in Part I | SCOBI/escapeLGD name | Meaning |
| --- | --- | --- |
| `a_d` | n/a (implicit) | truth: daytime adult escapement |
| `r` | `wc_prop` (escapeLGD); `SampleRate` (SCRAPI smolts) | sampling fraction; adults have only r, smolts have r and GE |
| `w` | `wc` (escapeLGD adults); `SampleCount` (SCRAPI smolts) | observed count |
| `e_sd` | `GuidanceEfficiency` (SCRAPI smolts only) | bypass guidance efficiency -- NOT used for adult window counts |
| `p_sd` | `Ptrue`, `SR` (SCRAPI smolts) | smolt-trap combined detection = SampleRate * GuidanceEfficiency |
| `p_n` | `p_night`, `nightPassage_rates` | nighttime passage proportion |
| `p_f` | `p_fa`, `pf` | fallback proportion |
| `c_sd` | `Tally` (in passdata) | daily trap count |
| `Strat` | `Strat`, `Stratum`, `Collaps`, `stratum` | collapsed week label |
| `PGrp` | `PGrp`, `GenStock` | primary stock group |
| `SGrp` | `SGrp`, `fwAge`, `GenSex` | secondary group |
| `w_j` | `pW` in `PBT_log_likelihood` | true wild proportion |
| `t_i` | `tagRate`, `tagRates` | PBT tag rate per hatchery |
| `theta` | `par` in optim() calls | parameter vector being estimated |
