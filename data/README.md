# Data

Real MY2025 steelhead inputs from Lower Granite Dam. The simulate-and-recover
sessions are tuned to look like these shapes but generate their own data from a
known truth at the top of each script, so they do not read these files. Session 7
is the exception: it runs production `smoltEASE::SCRAPI2()` on the two MY2025 CSVs
below and checks the result by hand, so those two files are load-bearing for it.

## Files

- `MY2025_STHD_passage.csv` — daily smolt passage at the trap: `SampleEndDate`,
  `GuidanceEfficiency`, `SampleRate`, `SampleCount`, `CollectionCount`,
  `Collapse` (stratum), `Week`. This is the `passageData` argument of
  `smoltEASE::SCRAPI2()`.
- `MY2025_STHD_trapData.csv` — one row per trapped smolt: `MasterID`,
  `CollectionDate`, `GenStock`, `GenSex`, `fwAge`, `Rear` (W or HNC),
  `LGDMarkAD` (AD or AI), and related fields. This is the `smoltData` argument
  of `smoltEASE::SCRAPI2()` and the composition sample.
- `LGR_spill_winter.csv`, `LGR_spill_spring.csv`, `LGR_spill_summer.csv` —
  hourly Lower Granite spill by season (`Site`, `DateTime`, `HourlySpill`,
  `HourlyFlow`, gas readings). These are the spill covariate for the guidance
  efficiency model (`smoltEASE::prep_ge_data()` / `fit_ge_model()`).

## Not included

The PIT tag detection extract (spillway GRS and Little Goose GOJ detections)
that feeds `prep_ge_data()` as `dat_up` was shown but not provided as a complete
file. Session 5, the guidance-efficiency fit, is simulate-only and does not need
it; add the full extract here if a real GE fit is wanted later.
