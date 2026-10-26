# Code review: `analysis/clr_logpt_test.R` against `docs/preregistration_clr_logpt.md`

*Independent review required by protocol step 5, 26 Oct 2026, following `docs/clr_logpt_review_brief.md`. Reviewer: Claude (Opus), who has seen no E. coli, RPE1 or K562 truth file and no score or metric produced on them. Read: the script, the preregistration, the brief, the helpers it calls or copies (`regulons.R`, `regulondb.R`, `human_truths.R`, `baselines.R`, `keff_test.R`, `chips_bootstrap.R`, `build_ecoli_complex_truth.R`, `corum_exploratory/corum_robust.R`, `R/clr_attention.R`), the staged jobs and the K562 repro check output. Reviewed commit: c64d2a5 (line numbers below refer to it). Fixes: bdcc05b, 825b801, 93e5bc8, 1f45cfc.*

## Verdict

**The script is ready for the truth phase.** One protocol item has to be done first (P1 below: the preregistration's change log does not yet record the complex source). Adding that entry changes the preregistration hash, so the K562 reproduction check has to be rerun afterwards. Order before queuing 0910:

1. Let job 0902 (RPE1 t\*) finish and commit `docs/clr_logpt_selections/rpe1_*`. Do this **before** amending the preregistration: the selections fingerprint includes the preregistration hash, so a resume after the amendment would refuse.
2. Add the change-log entry for P1, plus the one-line clarification of interpretation 1, and commit.
3. Rerun `--mode repro` (K562 only, in-sample). `gate_repro` requires the current hash, and the check must pass again.
4. Add `after: 0914_clrlogpt_truth_k562` to `0915_clrlogpt_verdict.sh` (R2), then queue 0910 onwards.

## Findings fixed in the script

| # | Severity | Line (c64d2a5) | Rule | Finding | Fix |
|---|---|---|---|---|---|
| F1 | High (crash) | 145, 537–550 | Truths, human; Endpoints S5 | `build_truths_human()` calls `corum_sets()`, `collectri_sets()`, `dorothea_sets()` and `response_sets()`. These are defined in `analysis/human_truths.R`, which the script never sourced. The synthetic smoke test bypasses them, so the RPE1 (S5) and K562 runs would have stopped with "could not find function". | bdcc05b: source `human_truths.R`. |
| F2 | Medium | 863 | Statistics: "Any claim about a single primary endpoint outside the adoption rule uses Holm-adjusted p over the four primary tests"; Decision rules, post-adoption bullet 3 | The note "small but significant regulon loss within the margin" tested the **raw** `p_loss < 0.025`. The code checked `loss_holm` only for NA, so the intent was clearly Holm. The note was also skipped while S5 was pending, although it does not depend on the secondary family. | 825b801: use `loss_holm < 0.025` on P3/P4, which matches the "Regulon loss" outcome row, and evaluate the note whenever adoption holds. |
| F3 | Low | 839 | Decision rules ("Adopt only if all of these hold") | `adopt` was NA (verdict "PENDING") whenever the bootstrap summary was missing, even when rule 1, 2 or 4 had already failed and the outcome was already determined. | 825b801: three-valued `r1 && r2 && r4 && r3`. The verdict is PENDING only if rules 1, 2 and 4 pass and the bootstrap is missing. |
| F4 | Low | 874–876 | Protocol steps 1 and 6 | The verdict read `endpoints_rpe1.csv` even if the RPE1 run crashed between writing endpoints and `DONE`. It did not check the preregistration hash of `boot_summary.csv`. A non-synthetic verdict would also have accepted outputs of a `--synthetic` run written to the same `--out`. | 93e5bc8: use endpoints only from finished runs (with `DONE`), check the bootstrap's hash, and refuse synthetic inputs (tested: refused). |
| F5 | Low (resumability) | 580, 620, 635, 641, 790 | Protocol step 5 (a rerun only after a crash) | `fit.rds`, `scores.rds` and `jack_*.rds` were written in place, so a crash mid-write left a truncated cache that broke every resume. A jackknife cache was reused on column names alone (not on replicate count, halves or pairs). Per-draw bootstrap files were reused with no fingerprint. | 1f45cfc: atomic writes; jackknife caches keyed (MD5) to halves, pairs and preregistration hash; bootstrap fingerprint (genes, chips, pairs, hash, operator settings). |

Tests: `decide()` was unit-tested on constructed endpoint tables covering adoption, a missing bootstrap with rules passing or failing, bootstrap failure, repro failure, regulon loss, the small-loss note, and the one-sided directions. The old `decide()` fails these tests and the new one passes them. Synthetic smoke test after the fixes: `--mode truth` for m3d, ecomac, rpe1 and k562ess (`--synthetic 1 --jack 10`), a second run refused by `DONE`, a resume after deleting `DONE` (jackknife caches reused, endpoints byte-identical), `--mode boot` (6 draws, then resumed with 0 to run) and `--mode verdict`. No truth was read.

## Protocol findings (not fixable in the script; for the PI)

**P1 (must be done before scoring).** The preregistration says the complex source "is chosen on availability alone … and logged **here** with its release and SHA-256" (Truths), and protocol step 4 says to "record source, release, SHA-256 of the raw file and of the built member list, and the evaluability counts (no scores) **in the change log**, and commit". The change log has only "2026-10-22: fixed." The information is in `docs/clr_logpt_ecoli_truth_log.md`, and the header of `build_ecoli_complex_truth.R` points to a change-log entry that does not exist. Add the entry (it can cite the build log) before the first truth-phase run. See the Verdict for the consequences of changing the hash: rerun repro, and wait for 0902.

## Recommendations (low; not fixed)

- **R1, run-once is per `--out`.** `DONE`, `gate_repro` and `gate_selections` all look under `opt$out`. A truth run with another `--out`, after copying `repro_check.csv` and `selections.csv` there, is not refused. Locking `--out` to `results/clr_logpt_test` outside `--synthetic` would make "the truth phase, once" mechanical. Until then, it rests on the staged jobs, which use the default.
- **R2, verdict job order.** `0915` waits for 0912 and 0913 but not 0914. If K562 is not finished, the post-adoption note "K562 merged-CORUM difference vs attention ≤ 0" cannot be evaluated.
- **R3, clean-tree gate.** `gate_script_clean` covers only the script itself, not the helpers it sources or copies (`regulons.R`, `regulondb.R`, `human_truths.R`, the `clr` package in `R/`). `git_head` is stamped, but a dirty tree is not detected. Consider running `git diff --quiet HEAD -- analysis R` in the gate.
- **R4, unlocked options.** `fix_check` does not lock `--hubstrat` (descriptive) or `--expr` (the K562 input). Neither can touch P1–P4 or S1–S5.
- **R5, synthetic `DONE`.** `gate_order` accepts any `DONE`, including one from a synthetic run written to the default `--out`. The real run of that dataset would then be refused loudly, and since F4 the verdict refuses synthetic inputs.

## Checklist of the brief (all checked; consistent unless stated)

1. **Gate order.** Every gate comes before any truth read: `fix_check`, `gate_script_clean`, `gate_repro`, `gate_selections`, `gate_ecoli_truth`, `gate_order`. `load_expr()` reads expression and h5 metadata only (the half split uses perturbation p-values, not a truth). `build_truths()` runs after `gate_all()` in truth and boot modes.
2. **Hash in every output.** `wout()` stamps every CSV. `DONE` and the log carry the hash. The caches carry it in their fingerprints.
3. **Refusals.** A second run is refused by `DONE` (tested). `fix_check` locks the seed, 100 replicates, 200 draws, t = 9, α, keff, cap and `--quick`.
4. **Verdict.** It is written by the script to `verdict/verdict.csv`.
5. **Operator.** The same call chain as repro. `f$operator` is the row-normalised softmax with self-loops for empty rows (`.finish_operator`). Pᵗ comes from repeated products with no floor. The readout functions (l. 108–126) are byte-identical to e891f99 and 5881335 and match `corum_robust.R`'s `rowz_na`. The only difference: non-finite rather than NA values are set to 0, which is what preregistration step 3 asks for.
6. **Depths.** The primary is t = 9 for a3r and attention on the same Pᵗ. t\* is used only for `a3r_tstar` and `attention_tstar`, read from the committed selections.
7. **E. coli complexes.** Members come from the SHA-gated file. The script stops if any member gene is outside the universe. The primary is operon-excluded. Evaluability uses structure only and is the builder's rule. When the member check passes, its counts depend only on members and the operon map, so they necessarily equal the builder's. The sensitivities and the rp module are as preregistered. The script's `merge_sets` is identical to the builder's. `cutree(h = 0.5)` joins at height ≤ 0.5, that is Jaccard ≥ 0.5.
8. **E. coli regulons.** `regulons_SC` is the same as `baselines.R` `bench$SC`. TF regulons, sigmulons and operons follow the definitions in `keff_test.R` `load_truths`.
9. **Human truths.** CORUM merged (3–200, then merge), raw, and no-RP (≥ 30%, the preregistration's regex), plus CollecTRI, DoRothEA and KD response (half B). The RPE1 gene set is checked against the selections fingerprint.
10. **AUPR.** `ap_fast`; NA ranks last; tie blocks end. Per-complex AP uses `regulon_ap_table` (`macro_aupr.R` definitions, ≥ 10 positives).
11. **Endpoints table.** It matches the preregistration row by row (dataset, truth, comparator, hypothesis). d = a3r_t9 minus the comparator.
12. **Jackknife.** As `baselines.R`: all halves are drawn up front after `set.seed(20261022)`, in bench order (complexes, then regulons SC, then the rest), with the same halves for every method, se = SD of the paired differences, and d from the full data.
13. **Directions.** Superiority `1 − Φ(d/se)` (H0: d ≤ 0). NI `1 − Φ((d + 0.005)/se)` (H0: d ≤ −0.005, so the margin sign is correct). Loss `Φ(d/se)` (H0: d ≥ 0). All rejections are at p < 0.025.
14. **`decide()`.** Rule 1 requires p < 0.025 and d ≥ 0.005 on both P1 and P2. Rule 2 is NI on both P3 and P4. Rule 3 is the P1 lower bound > 0 and the P3 lower bound > −0.005. Rule 4 is the repro check. Holm over P1–P4 is computed for single-endpoint claims. Holm over S1–S5 is computed only if adoption holds and all five S endpoints are present. The outcome rows mirror the table (F2 and F3 fixed).
15. **Bootstrap.** `draw_cols` is verbatim from `chips_bootstrap.R`. The gene filter and the jitter for constant genes are the same. The seed is `20260924 + b`, set in the master before drawing and again in the worker. MI, CLR, operator and readouts are refitted per draw. Percentiles of a3r − attention for P1 and P3.

**Leakage.** None found. No choice depends on a truth: the depths come from committed label-free selections, the evaluability fallback from structure only, and the operator settings are locked. Truths that are empty or have no positives are dropped on structure alone, and a missing primary truth stops the run. The logs print AUPRs after each dataset. Nothing downstream can change in response, because the remaining runs are gated and fixed.

## The seven interpretations in the brief

1. **One `set.seed(20261022)` per dataset run, with the draws in table order within the dataset.** Consistent in substance. It is fixed before any truth, cannot be chosen by outcome, and puts the primary truths first. It differs from the literal reading, which would be one stream across datasets in table order (P1 M3D, P2 EcoMAC, …). A single stream is impractical with separate resumable jobs, and it only changes which random halves are drawn. **Log it as a one-line clarification in the change log, together with P1, before scoring.**
2. **Holm over the four primary loss p-values, applied to P3 and P4.** Consistent. It is the natural reading of "Holm-adjusted one-sided p (for d < 0)" with the preregistration's four-test family, and it is conservative. F2 makes the post-adoption note use the same adjustment.
3. **S3 and S4 use the 0.005 NI margin; S5 is superiority vs attention.** Consistent: this is stated in the endpoint table.
4. **K562 t\* is skipped when its selections are not committed.** Consistent, since K562 is descriptive. If K562 t\* is ever reported, say that it was selected on the full half-A genes and scored on the 5,671-gene subset.
5. **Descriptive truths keep only attention, a3r, CLR and |r| (plus their hub-stratified versions).** Consistent: no claim depends on them. Note that this applies to the E. coli descriptive truths too, not only to the human ones.
6. **Evaluability is recomputed and logged.** Consistent. Whenever the run proceeds, the recomputation necessarily equals the builder's decision (see item 7 above).
7. **The verdict reads whatever exists.** Consistent with the gatekeeping. After F3, PENDING appears only when rules 1, 2 and 4 pass and the bootstrap is missing. The S family stays unevaluated until RPE1 is DONE (F4). See R2 for the job order.
