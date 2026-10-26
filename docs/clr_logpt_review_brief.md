# Review brief: `analysis/clr_logpt_test.R` against `docs/preregistration_clr_logpt.md`

For a reviewer who has seen no E. coli or RPE1 result. Do not look at `results/clr_logpt_test/truth/`, `docs/clr_logpt_ecoli_truth_log.md` counts for scores, or any AUPR on E. coli or RPE1 (there are none yet: the truth phase has not been run). The script was smoke-tested only on generated data with a generated truth (`--synthetic 1`). Reviewed commit: the one that last touched the script (`git log -1 -- analysis/clr_logpt_test.R`). Line numbers below are for commit `c64d2a5`.

The job is to find lines that violate the preregistration, not to improve the method. Name the line and the rule for each finding. A bug fix after the first truth-phase run is allowed only if it is shown to violate the file (protocol step 5), so anything wrong should be found now.

## What the script does

Modes: `repro` (done, passed), `selections` (done, committed), and new: `truth` (one dataset per run), `boot` (M3D cluster bootstrap), `verdict` (rules applied to the outputs). The new code starts at the banner `## ---- truth phase, bootstrap and verdict` (about l. 330) and runs to the end of the verdict block (about l. 910); everything above is the repro/selections code that already ran.

## Checklist (preregistration section, then where to look)

**Gates and protocol (Protocol and gate)**
1. Gate order. `gate_all()` (about l. 358): prereg gate (at script start, l. 137) then script committed, `gate_repro`, `gate_selections` (m3d, ecomac, rpe1; k562ess only if its selections file is committed), `gate_ecoli_truth` (E. coli), `gate_order` (E. coli before RPE1 before K562). No truth file is read before `gate_all()` returns: check that `load_expr()` reads expression only and `build_truths()` comes after the gates. `gate_script_clean` is an addition (audit trail), not in the file.
2. Hash in every output. `wout()` adds `prereg_commit`, `git_head`, `script_commit`, `synthetic` to every CSV; logs and `DONE` carry the hash too.
3. A second run of a finished dataset is refused (`DONE`); settings other than the preregistered ones are refused (`fix_check`: seed 20261022, 100 jackknife replicates, 200 draws, t = 9, alpha 0.5, keff 10, cap 50, no `--quick`).
4. Verdict written by the script to `verdict/verdict.csv` (`decide()` and the verdict block).

**Operator and readouts (Operator)**
5. Truth-phase fit (`score_dataset`): `estimate_mi(bins = "hg", transform = "none")`, `calibrate(normal, stouffer)`, `build_operator(alpha, softmax_keff = 10, softmax_cap = 50)`, the same call chain as `repro` mode (which reproduced 0.1054351). Pt by repeated `Pt %*% P` with no floor. a3r, a2r, a1 and attention are the functions at l. 108-126, unchanged since the repro check; confirm by diff against commit e891f99/5881335 if wanted.
6. Depths: primary t = 9 for the new readout and for attention (same Pt). t* only for the secondary `a3r_tstar` and `attention_tstar`, read from `docs/clr_logpt_selections/<dataset>_selections.csv`.

**Truths (Truths)**
7. E. coli complexes: members read from the committed, SHA-gated `ecoli_complexes_members.tsv` (variants main, size2, nomerge, noribo, rp_module); the merge was done by the builder, not here. Check `build_truths_ecoli`: operon-excluded pairs for the primary `complexes`; the evaluability fallback (>= 200 positives and >= 10 complexes with positives, else operon-included) is computed from structure only; sensitivities are `complexes_operon_incl`, `complexes_size2`, `complexes_nomerge`, `complexes_noribo`; the rp module is scored by `regulon_coherence` (descriptive).
8. E. coli regulons: primary `regulons_SC` = `make_regulons(conf_keep = c("C","S"))`, operon pairs excluded, default sizes 5-500; descriptive `regulons_all`, `TF_regulons`, `sigmulons`, `operons` with the `keff_test.R` / `granularity.R` definitions.
9. Human truths: `CORUM_merged` (`corum_sets` 3-200, then `merge_sets`, a verbatim copy of the builder's rule), `CORUM_raw`, `CORUM_noRP` (merged complexes with >= 30% RP members removed, RP regex as in the file), CollecTRI, DoRothEA A-B and A-C, KD response (half B). K562 input is the committed truth-gene expression file (as `repro`), RPE1 is the 8,747-gene half-A matrix (checked against the selections' gene fingerprint).
10. Pooled AUPR is `ap_fast` of `S[cbind(i, j)]` over `comembership_pairs`; non-finite scores rank last (NA); ties counted at block end.

**Endpoints and statistics (Endpoints, Statistics)**
11. `ENDPOINTS` table (about l. 665): P1 M3D complexes vs attention (sup), P2 EcoMAC same, P3 M3D regulons SC vs attention (NI), P4 EcoMAC same, S1/S2 regulons SC vs CLR (sup), S3/S4 complexes vs CLR (NI), S5 RPE1 merged CORUM vs attention (sup). New readout is `a3r_t9`.
12. Jackknife (`jack_dataset`): `set.seed(20261022)`, then for each truth in table order 100 draws of `sample(genes_in, floor(n/2))`; AUPR over pairs with both genes kept; same halves for every method; se = SD of the 100 paired differences; d = full-data AUPR difference (as `baselines.R`).
13. p-values: superiority `1 - pnorm(d/se)`, NI `1 - pnorm((d + 0.005)/se)`; `p_loss = pnorm(d/se)`.
14. Decision rules in `decide()`: rule 1 (P1 and P2 p < 0.025 and d >= 0.005), rule 2 (P3 and P4 NI p < 0.025), rule 3 (bootstrap P1 lower > 0 and P3 lower > -0.005), rule 4 (repro passed); adopt = all four. Holm over P1-P4 for single-endpoint claims; Holm over S1-S5 only when adopted (a pending S5 keeps the secondary family unevaluated). The outcome table rows and the post-adoption notes mirror the file's table.
15. Bootstrap (`boot` mode): `draw_cols` is a verbatim copy of `chips_bootstrap.R` cluster mode; seed `20260924 + b`, b = 1..200, `set.seed` before drawing and again in the worker; MI, CLR, operator and readouts refit per draw; percentile 2.5-97.5 intervals of `a3r - attention` on P1 and P3.

**Descriptive outputs** (no claims): per-complex AP (`regulon_ap_table`, >= 10 positives), hub-stratified AUPR (10 strata of pi_i + pi_j, teleport 1e-3, as `baselines.R`), a2r, a1, t* results, truth sensitivities, attention vs CLR on complexes.

## Places where I had to interpret the preregistration (confirm or overrule)

1. "Halves are drawn per truth, in the order of the endpoint tables, after `set.seed(20261022)`": implemented as one `set.seed(20261022)` per dataset run, then the draws for that dataset's truths in table order (primary truths first, then sensitivities and descriptive truths), because datasets run as separate resumable jobs. The stream is therefore per dataset, not one stream across datasets.
2. Holm-adjusted p for the "regulon loss" outcome: Holm over the four primary loss p-values (`pnorm(d/se)` for P1-P4), applied to P3 and P4.
3. S3 and S4 (NI against CLR) use the same margin 0.005 and p formula as P3; S5 is superiority against attention.
4. If the K562 selections file is not committed, K562 t* readouts are skipped (K562 is descriptive); the selections that exist were computed on the full half-A genes, the K562 truth run uses the 5,671 truth-gene subset.
5. Truths that are huge on human data (descriptive CollecTRI, DoRothEA, KD) keep only attention, a3r, CLR and |r| at t = 9 to bound memory.
6. Evaluability for the E. coli complex endpoint is recomputed in the script (same rule as the builder) and logged to `evaluability.csv`; the build log says no fallback is needed.
7. Verdict reads whatever of P1-P4, bootstrap, S1-S5 exists; without the bootstrap it reports PENDING; without RPE1 the adoption decision is made but the Holm family over S1-S5 is left unevaluated.

## Running it (after review, scheduler only)

Staged scripts are in `jobs_pc/staged/` (0910 M3D, 0911 EcoMAC, 0912 bootstrap, 0913 RPE1, 0914 K562, 0915 verdict; each with memory, CPU and `after` headers). They are not in `jobs_pc/sched` and must not be moved there until (a) this review is done, (b) any fixes are committed, (c) the RPE1 selections file is committed (0913 needs it). Smoke test without truth: `Rscript analysis/clr_logpt_test.R --mode truth --dataset m3d --synthetic 1 --jack 10 --out <scratch>` (then ecomac, boot, verdict).
