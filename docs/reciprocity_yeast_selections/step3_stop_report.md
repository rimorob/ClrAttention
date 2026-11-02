# Yeast reciprocity preregistration: STOPPED at protocol step 3 (structural dye check)

*1 Nov 2026. No truth has been downloaded or opened. Preregistration fixed in 02cbcfd.*

## What ran

| Step | Result |
|---|---|
| 2 download (job 0980) | 920 MB; all files complete (see `download_log.md`) |
| 3 matrix (job 0982) | 6,169 ORFs × **652** mutants. Of the 700 responsive columns, 3 were dropped as unresolved (sus1, rsc30, tlc1) and **45 failed the own-M verification** (own M > −1 in all 45; median −0.68, range −0.98 to +1.95). 0 ORFs dropped for missing values. Logs: `dropped_columns.csv`, `masking_log.csv`; hash `deleteome_matrix.sha256` |
| 3 check (jobs 0982, 0983) | **STOPPED by the structural dye check** (preregistration, gate F item 1: "each dye-swap pair must contain one hybridization of each orientation. Any violation stops the stage and is logged.") |
| 5 R check (job 0981, synthetic) | passed (15/15 checks; synthetic metrics reproduced to about 2e-15). Kept out of git (`results/reciprocity_yeast/R_check_preliminary_20261101.csv`) until steps 3–4 are complete, to keep the protocol order |

## The failure

- **Rule part 1, exactly one reference channel per hybridization: holds.** 1,387/1,387 mutant hybridizations (GSE42527) and 755/755 wild-type hybridizations (GSE45115).
- **Rule part 2, one hybridization of each orientation per mutant: fails for 15 of 701 mutants** (`dye_check_failures.txt`):
  - 13 have a single hybridization: atg24, bcs1, bik1, pus9, rad28, sae2, scp160, set5, tif2, yil001w, ymr031c, yor338w, ypl199c;
  - cho2 has two hybridizations with the test strain in channel 1;
  - dbp3 has two with the reference in channel 1.
- **Wild types (not yet reached by the code; checked by hand):** GSE45115 titles are irregular (for example `wt-sv017_2-del`, `wt-ea024-1-ab`, single `-a` or `-b` cultures). Grouped by title stem there are 407 groups: 345 balanced pairs, 60 single hybridizations and 1 triple. So part 2 would also fail for wild types.
- **Code bug found and fixed (no effect on the stop):** the grouping key did not strip GEO title prefixes such as `[hs1991] `, and it paired mutant cultures by title stem instead of by deleted gene. Before the fix, 221 "pairs" failed. With the fix, the 15 genuine violations above remain. The fix also changes which mutants are "the first 10 in C-locale order" for the implementation check. The download fetched sample tables for the old, prefixed selection, so the fixed selection needs about 20 more sample tables.

## What a decision would need (PI; logged deviation before any truth)

The rule's purpose, orientation verification, is met by part 1. Part 2 fails because some mutants and wild-type cultures were hybridized only once or twice in the same orientation. Options:

- **(a)** Keep part 1 as the check. For the orientation-rule reproduction and the WT matrix, use only balanced groups. That drops 15 mutants from the check pool (the published-M matrix is unaffected) and about 61 WT groups from gate F's WT matrix.
- **(b)** Keep part 1 and average whatever hybridizations exist, after orientation.
- **(c)** Treat the stop as final: gate F is not computable, which gives verdict 1 for the whole study (yeast not consumed).

None of these is applied. Work stops here.
