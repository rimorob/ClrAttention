# Next steps (kept current by the lead; each item says what triggers it)

## After the collaborator approves the general direction
1. **Max-statistic bootstrap for the K562 selection.** tcaCLR was selected among about 20 read-outs on K562. Resample genes, recompute all read-outs on each resample, and take the maximum of the difference statistic to get simultaneous intervals and adjusted p-values for the selected read-out (handles the correlation among read-outs; the winner's-curse bias, about sqrt(2 ln m) standard errors for m candidates, is reported alongside). Keep a ledger of every method variant scored on each truth and correct over that ledger from now on. Does not apply to the HepG2 test (different truth, run once).
2. **Pre-specified test of tcaCLR2 on fresh data.** Freeze the rule (docs/plans/tcaclr2_definition.md), the dataset, the truth, the metric and the comparators in git before any score exists; read as non-inferiority against tcaCLR. Lead candidate: DepMap gene dependency scores with Complex Portal or CORUM complexes (check overlap with the data already used first). Report as a pre-specified test, not a preregistered one, unless it is registered externally.

## Pending decisions or small work
0. **scGPT extension (docs/plans/scgpt_extension_plan.md).** Plan written, nothing run. Needs the PI's decisions: download RPE1 and Jurkat single-cell files (8.7 GB each; n = 6 datasets makes every pre-specified contrast attainable, n = 4 allows no decision), permission to delete each raw h5ad after the control cells are extracted (disk budget), and whether Complex Portal complexes stay as a descriptive stratum. Expected outcome stated in the plan: CLR-NS vs the best scGPT readout is expected to fail (token-embedding cosine beat it on hESC, hHep and K562).
3. Register edits in the Evidence and Losses text of the short pre-read (colloquialisms such as "answer in hand", "used up", "after many looks"); numbers and claims unchanged. Waiting for the PI's go.
4. Verify the bibliographic details in the short pre-read's reference list (Schmidhuber is cited as 1992, the journal paper; the main write-up says 1991, the technical report).
5. tcaCLR2 minus aCLR intervals (experiment-refit jackknife) are not computed; only point estimates are shown.
6. Long document: figures `jim_diffusion.png` and `depth_log_summary.png` print under 7 pt; dates remain in its text; the "same reach" control rows were removed with the stacked variants.
7. The 2007 CLR statement (paper Euclidean-only, software both) remains in older repo files not yet updated (abstract_draft.md, figures/README.md, R/clr_attention.R, R/calibrate.R, docs/plans/fig1_review_round2.md, STATUS.md lines about CLR.zip).

## Unscheduled hypotheses and probes
8. **Realistic planted simulation** (hubs, overlapping modules, resampled real complex structure) before any simulation-based depth rule is trusted; the clean-clique simulation disagreed with the real truths about depth.
9. **aCLR retune.** Its best of 819 cells beat the registered cell by 14 to 69% relative (wide keff, unclipped first score); a hypothesis that needs a fresh truth.
10. **Exploratory read-out table**: needs tcaCLR_allentries added to the registered list first, logged before any score output.
11. **Translational-symmetry probe** (docs/plans/fable_translational_symmetry_review.md): planted ring through the pipeline first; stop if not recovered.
12. **Island/conductance diagnostic, StARS instability curve** from the adaptive-methods review: not built.

## Housekeeping
13. Disk cleanup per docs/disk_budget.md when the paper is ready.
