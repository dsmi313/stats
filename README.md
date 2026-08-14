# stats

Working repository for a 45-minute talk to IDFG colleagues: *"uncertainty in
escapement estimation at Lower Granite: adults and smolts, one framework."*

## What this is

Nine self-contained R sessions that reverse engineer the two production
escapement tools at Lower Granite Dam — `escapeLGD` (EASE, adults) and
`smoltEASE` (SCRAPI2, smolts) — treating them as one estimator structure applied
at two life stages. Each session writes a toy version by hand, then points at
the exact production file and function that does the same job.

## Where to start

1. Read `PLAN.md`. It states the one idea the talk is built on, the five-part
   session format, and the nine sessions in order.
2. Read `docs/session02_shared_skeleton.md` for the component-by-component
   comparison of the adult and smolt estimators. That table is the spine.
3. Work the sessions in order in `R/`. Each writes a plain-English explanation
   to `docs/sessionNN_explain.md`.

## What is not here

Everything cut from the earlier construction-oriented plan is in `BACKLOG.md`,
each item with one line saying why it was cut. Nothing was deleted.

## Ground rules

No script over 80 lines or over a minute to run. No banner comments. No
placeholder headings. No dates — those get added by hand.
</content>
