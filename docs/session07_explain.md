How I would explain session 7 in three minutes

Every number in the talk so far was simulated. This one is not: it is the real
MY2025 steelhead run, 61 days of trap counts and 1074 fish. Be clear about what this
shows and what it does not: we did not independently verify SCRAPI2. We cannot, by
hand. What we can show is that the weekly number SCRAPI2 prints is nothing more
than the plug-in expansion, summed.

Take 05/12/2025: 88 fish counted at a 5% sample rate and 18% guidance is that count
over the fraction seen, about 9888 fish. Its week has 5 days; add their expansions
and you get 47444. SCRAPI2 prints exactly 47444 for Week 20, because line for line it is
the same division summed the same way: count over sample rate times guidance,
summed by week and rounded (SCRAPI2.R line 306, printed at line 325). There is no
separate model doing the counting and no magic in the machine; the whole run is
this one division, repeated day after day and then grouped.
