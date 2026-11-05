# DREAM5 truth stage: items the lead must state in the truth-stage log

From the independent step-6 review of the DREAM5 code (branch accuracy-review). The design file (`docs/plans/dream5_design.md`) is FIXED and is not edited; these are the points the truth-stage log (dated change-log entry before any truth is opened) has to settle or acknowledge. Review item numbers in brackets.

## Settle before the first scoring run (switches have no defaults; `analysis/dream5_truth_run.R` refuses to start without them)

1. **Tie switches knowing GRNBoost2 (N4).** `--ties_at_cap` (keep_tied | drop_tied | list_order) and `--ties_in_list` (end_of_block | expected) interact with GRNBoost2 and every sparse list: GRNBoost2 writes only the edges it ranks and fills the unlisted edges with 0, while the paper says that pairs not in a submitted list are randomly ordered. Decide the tie switches knowing that (a) a zero-filled block is one huge tie, (b) `end_of_block` credits every positive in a tied block with the block-end precision (the `ap_fast` convention; higher than ranking the block's positives last, lower than the expected-tie AP in the cases tested), and (c) `keep_tied` can make the list much longer than the cap when the block that straddles rank 100,000 is the zero block. State what the official scorer does and which switch values reproduce it.
2. **TF-TF pairs and `collapse_max` order (N5).** `collapse_max` is applied to the score matrix before the truth jackknife mask is applied (the jackknife keeps/drops TFs and genes after one orientation per TF pair has been chosen). If the official scorer collapses after restriction, the two orders can differ; state which one is meant.
3. **Pair universe, cap order, jackknife cap.** `--universe`, `--cap_order` (only with `--universe gs`; with `--universe all` supplying it, or the universe files, is now an error) and `--jackknife_cap` (always required, because the truth jackknife always runs).

## Acknowledge in the log (known limits of the code; no change made)

4. **Replicate index range (N8).** Both seed streams run over replicate index b = 1..100: SE_B group-level refits use seed 20263000 + b, the supplementary series-level refits use seed 20264000 + b. The truth jackknife uses its own fixed seeds (`D5$seeds`). The driver validates `--refit_n` and `--series_n` in 0..100.
5. **PIDC (N9).** PIDC scores are used with NaN set to 0 and assuming non-negative weights; the pipeline does not check non-negativity. State this when PIDC results are reported.
6. **Job-script settings are not versioned (N10).** The refit and score jobs ran from `jobs_pc/` scripts (thread counts, worker counts, memory) that are not under version control. `results/dream5/PROVENANCE.txt` records the code commits, the environment and the thread settings, and the S3 check shows the refit hashes are reproduced at a different worker count (16 vs 8).
7. **Scorer gate and official numbers.** The official scorer gate (agreement of this scorer with the official scoring scripts) and the official p-value and overall score are NOT implemented. If the gate fails at 1e-3, a logged code edit is needed (including the AUPR integration convention: step interpolation vs trapezoid, the point at recall 0), before any verdict is read.
8. **Not implemented.** The replicate-averaged scoring option, non-contrast pair tables and exploratory variants. Only the pre-registered contrast (attention vs CLR) on Y-E and Y-M is produced, as `verdict.csv` plus supporting tables.
9. **All-genes tree runs are not scored.** The all-genes GENIE3 / GRNBoost2 runs (`*_all`) exist as label-free score files and are hashed in the manifest, but `d5_truth_pipeline` scores only the methods it is given (TF-regulator tree runs by default); the all-genes runs are not scored by the pipeline.

## Behaviour changed by the review fixes (for the log)

- Truth edge table: exactly two columns, or three with a 0/1 label (label-0 rows dropped and counted, never positives); self-edges dropped and counted; unknown IDs and repeats counted. The driver prints all counts.
- All options and switches are validated before the guards, before the net or truth file is read, and before any structure count is written (also with `--structure_only 1`, which now needs the switches too).
- verdict.csv and the primary outputs are written first; the supplementary series-level jackknife runs afterwards inside `tryCatch` and a failure is written to `SUPPLEMENTARY_SERIES_FAILED.txt` without touching the primary files. With `--series_n 0` the guard does not demand the series manifest.
- The guard also requires `analysis/dream5_lib.R`, `dream5_score_lib.R`, `dream5_truth_run.R`, `regulondb.R` and `regulons.R` to be tracked, unchanged and pushed.
- A net4 row whose d or SE is not finite reads `Inconclusive` (new column `verdict_note` states why), not the 'descriptive' label.
- The dream5 expected-tie AP is now `d5_ap_expected` (the `ap_expected` of `analysis/regulons.R` is no longer shadowed).
