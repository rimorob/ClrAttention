# DREAM5 results report (protocol step 8)

Draft by Claude (Sonnet 5.5) for the PI, 2026-11-11. Branch `accuracy-review`. Facts and tables only. Every number below is copied from a result file named next to it, except where a line says "arithmetic" and shows the sum. Nothing was re-run for this report; no score, code, design file or job was touched. The design is `docs/plans/dream5_design.md` (FIXED 2026-11-04); the full tables are in `docs/dream5_selections/truth_stage_log.md` sections 6 and 7 and in `results/dream5/truth_net1`, `truth_net3`, `truth_net4`, `results/dream5/scorer_gate`.

Reading rule from the design: any text that cites Y-M or Y-E must cite both. They are shown together in section 2 and are never separated in this report.

## 1. Scope and what was done

**Question (design).** Whether CLR attention beats CLR on an organism and a truth the project had never seen (*S. cerevisiae*, DREAM5 net4), for regulon co-membership (Y-M, the headline construct) and for the benchmark's own directed TF-to-target readout (Y-E). A second, descriptive question: where CLR, CLR-NS, CLR attention, GENIE3, GRNBoost2, PIDC and |r| stand on the three DREAM5 networks (1, 3, 4) in the organisers' scoring.

**Timeline and commits.**

| Date | Step | Commit or file |
|---|---|---|
| 2026-11-04 | Design FIXED (rev 1, amended 2026-11-04); the PI read it and approved it as written | `docs/plans/dream5_design.md`, change log 2026-11-04 |
| 2026-11-05 | Four label-free clarifications from implementation, accepted by the PI | design change log 2026-11-05 |
| 2026-11-06 | Official scorer extracted and read; truth-stage definitions adopted before any score of ours existed | design change log 2026-11-06; commits f84d5e2, e7cf498 |
| 2026-11-10 | Scorer gate run | `results/dream5/scorer_gate/gate_report.md`; commit 95f78eb |
| 2026-11-11 | Truth stage: our methods scored once per net | net1 and net3 commit dd052ea; net4 commit 3f54c8b; design change log entry commit 7a20040 |

**Scorer gate (`results/dream5/scorer_gate/gate_report.md`).** PASS. 36 of 36 published methods (29 challenge teams, 6 off-the-shelf methods, Community integration) reproduced by the R scorer, the Python port and the matrix path, to the pass rule of 1e-3 on the six per-network AUPR and AUROC values and on the three aggregate scores (required: 3 methods). Largest difference to the published table: 5.79e-4 for the per-network AUPR/AUROC values (Correlation3), 4.93e-4 for the aggregate scores (MI2). Between our implementations: maximum |R - port| 8.82e-12, maximum |R - matrix path| 3.89e-15. The gate checks the scorer against the organisers' published predictions; it does not check any of our methods.

**Truth stage.** One primary run per net, with the switches fixed on 2026-11-06 (`run_info.txt` in each `truth_net*` directory): universe = gs, cap_order = cap_then_universe, tftf = keep_both, ties_at_cap = drop_tied, ties_in_list = expected, jackknife_cap = restrict_after_cap, official trapezoid readout for Y-E, 200 truth-jackknife replicates, net4 additionally 100 group-level refits (SE_B) and 100 series-level refits (supplementary). Scheduler jobs 3000 (net1), 3020 (net3), 3010 (net4). The run guard passed (design, data log, manifests; 21 score files hash-verified; all 100 group-level and 100 series-level refits recomputed and matching their committed hashes). Truth edges kept: net1 4,012; net3 2,066; net4 3,940 (`structure_counts.csv`), with 0 rows dropped for unknown IDs, label 0, self-edges or repeats (log section 7). No code, switch, seed or endpoint was changed after any score was seen; the one look at the truth has been used. No failure and no re-run.

**Y-M evaluability gate (net4, `results/dream5/truth_net4/ym_gate.csv`).** 90 regulons with 5 to 500 targets (required at least 20) and 141,943 positive co-membership pairs (required at least 1,000): evaluable = TRUE.

**Y-E pair universe (net4, log section 2 and `point_estimates.csv`).** 114 TFs with at least one gold edge times 1,994 genes in any gold edge, minus self pairs: N = 227,202 pairs, 3,940 positives.

## 2. Net4 (S. cerevisiae): the two co-primary endpoints together

Source: `results/dream5/truth_net4/verdict.csv` (verdicts written by the script before interpretation). Contrast: CLR attention versus CLR; d = log2(AUPR_attention / AUPR_CLR); delta = log2(1.10), the smallest effect of interest (10% relative AUPR).

| Quantity (column in `verdict.csv`) | Y-M: co-membership AUPR (`ap_fast`) | Y-E: official directed edge AUPR (trapezoid) |
|---|---|---|
| AUPR, CLR attention (`aupr_a`) | 0.1111 | 0.02164 |
| AUPR, CLR (`aupr_b`) | 0.105 | 0.02076 |
| Ratio attention / CLR (`ratio`) | 1.0583 | 1.0424 |
| Difference in AUPR (`difference`) | 0.00612 | 0.000879 |
| d = log2 ratio (`d`) | 0.0817 | 0.0599 |
| delta = log2(1.10) (`delta`) | 0.1375 | 0.1375 |
| SE_J, truth jackknife, 200 replicates (`SE_J`) | 0.0225 | 0.0684 |
| SE_B, replicate-group jackknife, 100 refits (`SE_B`) | 0.0173 | 0.0373 |
| SE = sqrt(SE_J^2 + SE_B^2) (`SE`) | 0.0284 | 0.0779 |
| 95% interval for d, d +- 1.96 SE (`ci_lo`, `ci_hi`) | 0.0260 to 0.1374 | -0.0928 to 0.2125 |
| Interval from SE_J alone (`ci_J_lo`, `ci_J_hi`; reference, not the deciding interval) | 0.0375 to 0.1259 | -0.0742 to 0.1939 |
| p_sup (`p_sup`) | 0.00203 | 0.221 |
| p_harm (`p_harm`) | 0.998 | 0.779 |
| Equivalence test, d > -delta (`p_noninf_lo`) | 6.33e-15 | 0.00564 |
| Equivalence test, d < +delta (`p_noninf_hi`) | 0.0248 | 0.159 |
| **Verdict as written by the script** (`verdict`) | **Equivalent** | **Inconclusive** |

Design rule applied, from the design's "Endpoints and decision rule": Superior needs p_sup < 0.025 and d >= delta; Equivalent within 10% needs both equivalence tests below 0.025 (the 95% interval inside +-delta); Inconclusive is the remainder. For Y-M, d = 0.0817 is below delta, so Superior is excluded whatever p_sup is; both equivalence p-values are below 0.025 (the upper one is 0.0248), giving Equivalent. For Y-E, the interval -0.0928 to 0.2125 contains 0 and +delta (its upper end is above delta; its lower end is above -delta); the superiority and harm conditions are not met and only one of the two equivalence tests is below 0.025 (d > -delta, 0.0056; d < +delta, 0.159): Inconclusive.

Two points stated literally from the table. (a) The Y-M interval 0.0260 to 0.1374 excludes 0 and lies inside +-delta; its upper end is below delta by 8.8e-05 (arithmetic: 0.1375035 - 0.1374158). (b) On the AUPR-ratio scale the same limits are 1.018 to 1.100 for Y-M and 0.938 to 1.159 for Y-E (arithmetic: 2 to the power of each limit of the d interval).

**Pre-registered interpretation lines that apply (design, "Interpretation lines fixed in advance").**

- **Y-M is Equivalent within 10%, i.e. not Superior.** The fixed wording is: *the co-membership headline did not replicate in yeast*, with the interval. Estimate and interval: d = 0.0817 (AUPR ratio 1.058), 95% interval for d 0.0260 to 0.1374. The superiority claim is not confirmed in yeast net4. The lines "module gain with a direct-edge cost replicates in yeast" and "attention improved both" both require Y-M Superior and therefore do not apply.
- **Y-E is Inconclusive.** Interval reported, no claim: d = 0.0599 (AUPR ratio 1.042), 95% interval for d -0.0928 to 0.2125.

**The prediction recorded before the truth stage (design, hypotheses and PI decision 2 of 2026-11-03): Y-M positive (d > 0), Y-E not positive (d <= 0).** Compared plainly:

- Y-M: the point estimate is positive (0.0817), as predicted, and the interval excludes 0. Its size is below the 10% margin (1.058 against 1.10), and the pre-registered verdict is Equivalent, not Superior.
- Y-E: the point estimate is positive (0.0599), which is the opposite sign from the prediction (d <= 0). The interval -0.0928 to 0.2125 is wide and contains 0, so the data do not separate the two; the verdict is Inconclusive and no claim is made either way.
- The pattern "Y-M positive and Y-E not" is therefore not what the point estimates show (both are positive, 0.0817 and 0.0599), and the table does not establish a difference between the two.

## 3. Supplementary series-level jackknife (one line, no role in any verdict)

Source: `results/dream5/truth_net4/supplementary_series_jackknife.csv` (delete-half over the 59 `#Experiment` series of net4, 100 refits, seeds 20264000 + b; descriptive; not used in `verdict.csv`).

| Endpoint | SD_series (series-level) | SE_B_group (deciding, replicate-group level) | Ratio SD_series / SE_B_group (arithmetic) |
|---|---|---|---|
| Y-M | 0.0402 | 0.0173 | 2.32 |
| Y-E | 0.0508 | 0.0373 | 1.36 |

Limitation sentence (design, "What is exploratory"): because 33 of the 59 series contain several replicate groups, the group-level SE_B may be somewhat small. In both endpoints SD_series is larger than SE_B_group. The Y-M interval 0.0260 to 0.1374 that excludes zero relies on the group-level resampling (SE_B = 0.0173). No alternative interval from the series-level number is computed here, and it enters no decision (PI decision of 2026-11-04).

## 4. Net3 and net1: descriptive tables (no verdicts)

**Labelling.** Net3 (*E. coli*) is the same lineage as the headline: its data come from the M3D-derived compendium and its truth is RegulonDB, the source of the headline truths. It is never pooled with net4 and is not a replication. Net1 is a GeneNetWeaver simulation (E. coli RegulonDB structure plus 10% random edges), a synthetic mechanism check only. The design says net3 and net1 results do not count as support. No verdict is defined for them (no refits were run, so SE_B does not exist; the `verdict` column reads "descriptive (no verdict: needs net4 and SE_B)").

**Intervals are understated.** For net1 and net3 the only interval available is the truth jackknife alone (`ci_J_lo` and `ci_J_hi` of `verdict.csv`, which equal d +- 1.96 SE_J, checked by arithmetic; `SE_B`, `SE`, `ci_lo` and `ci_hi` are NA). They omit the experiment-resampling component and are narrower than the deciding interval used on net4.

### 4.1 Net3 (same lineage as the headline), attention versus CLR

Source: `results/dream5/truth_net3/verdict.csv`.

| Readout | AUPR attention | AUPR CLR | Ratio | d | SE_J | Jackknife-only interval for d (d +- 1.96 SE_J) | `verdict` column |
|---|---|---|---|---|---|---|---|
| Y-E (official trapezoid) | 0.06655 | 0.07513 | 0.8857 | -0.1751 | 0.1590 | -0.4867 to 0.1365 | descriptive (no verdict: needs net4 and SE_B) |
| Y-M (co-membership) | 0.291 | 0.2443 | 1.1912 | 0.2524 | 0.0320 | 0.1897 to 0.3151 | descriptive (no verdict: needs net4 and SE_B) |

All methods, `results/dream5/truth_net3/point_estimates.csv` (point estimates only):

| Method | Y-E AUPR (official trapezoid) | Y-E AUROC | Y-E pairs listed after cap and universe (`K_listed`) | Y-M AUPR (`ap_fast`) | Y-M AUPR, expected ties (`aupr_sens`) | Y-M AUPR, regulons with >= 10 targets (`aupr_min10`) |
|---|---|---|---|---|---|---|
| \|r\| | 0.06359 | 0.5721 | 8,101 | 0.2269 | 0.2269 | 0.2637 |
| CLR | 0.07513 | 0.5897 | 10,059 | 0.2443 | 0.2478 | 0.2806 |
| CLR-NS | 0.05701 | 0.5774 | 10,620 | 0.2423 | 0.2457 | 0.2805 |
| CLR attention | 0.06655 | 0.5773 | 9,768 | 0.291 | 0.291 | 0.3321 |
| GENIE3 (TF regulators) | 0.09558 | 0.6173 | 11,509 | not scored | not scored | not scored |
| GRNBoost2 (TF regulators) | 0.09598 | 0.6242 | 12,121 | not scored | not scored | not scored |
| PIDC | 0.0687 | 0.5841 | 11,194 | 0.2398 | 0.2398 | 0.2771 |

### 4.2 Net1 (synthetic), attention versus CLR

Source: `results/dream5/truth_net1/verdict.csv`.

| Readout | AUPR attention | AUPR CLR | Ratio | d | SE_J | Jackknife-only interval for d (d +- 1.96 SE_J) | `verdict` column |
|---|---|---|---|---|---|---|---|
| Y-E (official trapezoid) | 0.296 | 0.2572 | 1.1511 | 0.2030 | 0.0909 | 0.0248 to 0.3813 | descriptive (no verdict: needs net4 and SE_B) |
| Y-M (co-membership) | 0.3548 | 0.3626 | 0.9786 | -0.0312 | 0.0237 | -0.0777 to 0.0153 | descriptive (no verdict: needs net4 and SE_B) |

All methods, `results/dream5/truth_net1/point_estimates.csv` (point estimates only):

| Method | Y-E AUPR (official trapezoid) | Y-E AUROC | Y-E pairs listed after cap and universe (`K_listed`) | Y-M AUPR (`ap_fast`) | Y-M AUPR, expected ties (`aupr_sens`) | Y-M AUPR, regulons with >= 10 targets (`aupr_min10`) |
|---|---|---|---|---|---|---|
| \|r\| | 0.1887 | 0.7567 | 88,374 | 0.3475 | 0.3475 | 0.3697 |
| CLR | 0.2572 | 0.7835 | 87,798 | 0.3626 | 0.3713 | 0.386 |
| CLR-NS | 0.2417 | 0.7555 | 87,530 | 0.3249 | 0.3327 | 0.3475 |
| CLR attention | 0.296 | 0.8147 | 85,923 | 0.3548 | 0.3548 | 0.3834 |
| GENIE3 (TF regulators) | 0.303 | 0.8161 | 87,371 | not scored | not scored | not scored |
| GRNBoost2 (TF regulators) | 0.3363 | 0.8245 | 88,306 | not scored | not scored | not scored |
| PIDC | 0.2253 | 0.7662 | 87,558 | 0.3547 | 0.3547 | 0.3793 |

### 4.3 Net4, all methods (descriptive point estimates, same file type)

Source: `results/dream5/truth_net4/point_estimates.csv`. This table adds no inference to section 2. No per-method intervals are shown; the per-replicate values behind them are in `jack_edge.csv` and `jack_module.csv` and are not summarised here. GENIE3 and GRNBoost2 have no Y-M entry because the all-genes-regulator runs were not scored (section 6).

| Method | Y-E AUPR (official trapezoid) | Y-E AUROC | Y-E pairs listed after cap and universe (`K_listed`) | Y-M AUPR (`ap_fast`) | Y-M AUPR, expected ties (`aupr_sens`) | Y-M AUPR, regulons with >= 10 targets (`aupr_min10`) |
|---|---|---|---|---|---|---|
| \|r\| | 0.01886 | 0.5099 | 12,114 | 0.1016 | 0.1016 | 0.1048 |
| CLR | 0.02076 | 0.5175 | 12,769 | 0.105 | 0.1057 | 0.1083 |
| CLR-NS | 0.02157 | 0.5217 | 13,008 | 0.1052 | 0.1059 | 0.1084 |
| CLR attention | 0.02164 | 0.5232 | 12,685 | 0.1111 | 0.1111 | 0.1145 |
| GENIE3 (TF regulators) | 0.02083 | 0.5184 | 11,387 | not scored | not scored | not scored |
| GRNBoost2 (TF regulators) | 0.02084 | 0.5196 | 11,568 | not scored | not scored | not scored |
| PIDC | 0.02131 | 0.522 | 12,440 | 0.104 | 0.104 | 0.1074 |

Reading note, from the same file: on net4 the Y-E AUPR of all seven methods lies between 0.01886 and 0.02164, against a base rate P/N = 3,940 / 227,202 = 0.0173 (arithmetic), i.e. 1.09 to 1.25 times the base rate (arithmetic); the Y-E AUROC lies between 0.510 and 0.523. The Y-M AUPR lies between 0.1016 and 0.1111, against a base rate of 141,943 / 1,838,403 = 0.0772 (arithmetic; `n_pos` and `N` of `point_estimates.csv`), i.e. 1.32 to 1.44 times.

## 5. Official-score table (descriptive)

Source: `official_scores.csv` of each `truth_net*` directory: the organisers' AUPR and AUROC (official trapezoid definitions, the scorer validated by the gate in section 1) and the p-values from the organisers' precomputed null densities. Descriptive only: no p-value of this table is used in any verdict, and the net4 contrast of section 2 is the only inferential result. The values are our implementations of each method, not the challenge-era submissions, and comparison with the published challenge table is context only (not blind, not like-for-like).

| Net | Method | AUPR | AUROC | P_AUPR | P_AUROC | -log10 P_AUPR | -log10 P_AUROC |
|---|---|---|---|---|---|---|---|
| 1 | \|r\| | 0.1887 | 0.7567 | 3.63e-13 | 1.59e-32 | 12.44 | 31.80 |
| 1 | CLR | 0.2572 | 0.7835 | 1.18e-64 | 1.08e-60 | 63.93 | 59.97 |
| 1 | CLR-NS | 0.2417 | 0.7555 | 1.51e-49 | 1.71e-31 | 48.82 | 30.77 |
| 1 | CLR attention | 0.296 | 0.8147 | 2.79e-111 | 4.86e-106 | 110.55 | 105.31 |
| 1 | GENIE3 (TF regulators) | 0.303 | 0.8161 | 5.6e-121 | 1.99e-108 | 120.25 | 107.70 |
| 1 | GRNBoost2 (TF regulators) | 0.3363 | 0.8245 | 1.6e-172 | 4.31e-123 | 171.80 | 122.37 |
| 1 | PIDC | 0.2253 | 0.7662 | 1e-35 | 2.23e-41 | 35.00 | 40.65 |
| 3 | \|r\| | 0.06359 | 0.5721 | 3.53e-05 | 0.96 | 4.45 | 0.02 |
| 3 | CLR | 0.07513 | 0.5897 | 6.56e-10 | 0.0455 | 9.18 | 1.34 |
| 3 | CLR-NS | 0.05701 | 0.5774 | 0.00353 | 0.762 | 2.45 | 0.12 |
| 3 | CLR attention | 0.06655 | 0.5773 | 3.03e-06 | 0.766 | 5.52 | 0.12 |
| 3 | GENIE3 (TF regulators) | 0.09558 | 0.6173 | 1.04e-21 | 3.69e-11 | 20.98 | 10.43 |
| 3 | GRNBoost2 (TF regulators) | 0.09598 | 0.6242 | 5.57e-22 | 8.53e-15 | 21.25 | 14.07 |
| 3 | PIDC | 0.0687 | 0.5841 | 4.34e-07 | 0.262 | 6.36 | 0.58 |
| 4 | \|r\| | 0.01886 | 0.5099 | 0.989 | 0.907 | 0.00 | 0.04 |
| 4 | CLR | 0.02076 | 0.5175 | 0.129 | 0.0128 | 0.89 | 1.89 |
| 4 | CLR-NS | 0.02157 | 0.5217 | 0.0195 | 3.3e-05 | 1.71 | 4.48 |
| 4 | CLR attention | 0.02164 | 0.5232 | 0.0168 | 2.5e-06 | 1.78 | 5.60 |
| 4 | GENIE3 (TF regulators) | 0.02083 | 0.5184 | 0.111 | 0.00427 | 0.95 | 2.37 |
| 4 | GRNBoost2 (TF regulators) | 0.02084 | 0.5196 | 0.109 | 0.000934 | 0.96 | 3.03 |
| 4 | PIDC | 0.02131 | 0.522 | 0.0376 | 2.3e-05 | 1.42 | 4.64 |

Aggregate scores from this table (arithmetic, same definition as the official `go_all.m:66-68`, quoted in log section 2: AUPR score = mean over nets 1, 3, 4 of -log10 P_AUPR, AUROC score likewise, overall = mean of the two; these three-network aggregates are not written by the pipeline, they are computed here from the table above):

| Method | AUPR score = mean over nets 1, 3, 4 of -log10 P_AUPR | AUROC score = mean over nets 1, 3, 4 of -log10 P_AUROC | Overall = mean of the two |
|---|---|---|---|
| \|r\| | 5.63 | 10.62 | 8.13 |
| CLR | 24.67 | 21.07 | 22.87 |
| CLR-NS | 17.66 | 11.79 | 14.72 |
| CLR attention | 39.28 | 37.01 | 38.15 |
| GENIE3 (TF regulators) | 47.40 | 40.17 | 43.78 |
| GRNBoost2 (TF regulators) | 64.67 | 46.49 | 55.58 |
| PIDC | 14.26 | 15.29 | 14.78 |

**p = 0 and Inf.** The official code has no floor or cap on p or on -log10 p (`probability.m:2-5`, `go_all.m:66-67`; log section 6), and p is exactly 0 only for a score above the support top of the null density (AUPR 0.4134, 0.3033, 0.1217 and AUROC 0.9082, 0.8199, 0.6174 for nets 1, 3, 4; log section 2). None occurred in our table: all 42 p-values (7 methods, 3 nets, 2 metrics) are positive and finite, the smallest being 1.6e-172 (net1, GRNBoost2, P_AUPR), and every -log10 value is finite, so no aggregate score is Inf. Had one occurred it would have been reported as p = 0 and the score as Inf, never capped silently (design change log 2026-11-06, decision 8). None of the 36 published methods reached p = 0 either (log section 6).

## 6. Disclosure table: what is not blind, not like-for-like, or not controlled

Base: the design's "Limitations to state in any report" and "Independence" sections. Items added from the logs carry their source.

| # | Item | What is known | Source |
|---|---|---|---|
| 1 | Not blind | Published challenge results and ranking are known to the team in outline. We reimplemented every method; results are not like-for-like with the challenge-era ones. No claim is made about the challenge ranking beyond the descriptive table of section 5. | design |
| 2 | Net3 shares lineage with the headline | M3D-derived data and RegulonDB truth. Labelled "same lineage as the headline", never pooled with net4, not a replication. | design |
| 3 | Net1 is simulated | Descriptive mechanism check only. | design |
| 4 | Net4 truth is binding evidence | ChIP-chip binding plus conservation: incomplete and a weak proxy for regulation; yeast was the hardest network of the challenge, so differences may be undetectable. Inconclusive was a likely and legitimate outcome. | design |
| 5 | Absolute performance on net4 is near chance for every method | Y-E AUPR 0.01886 to 0.02164 against a base rate 0.0173 (1.09 to 1.25 times); Y-E AUROC 0.510 to 0.523 (0.5 is chance). Y-M AUPR 0.1016 to 0.1111 against a base rate 0.0772. | `truth_net4/point_estimates.csv`; arithmetic |
| 6 | No power computation | Power was not computed (no yeast truth structure before the truth stage). The realised SE is reported; nothing is claimed from intervals that include both +-delta. The Y-E interval is wide: it contains 0 and +delta. | design |
| 7 | Familywise error across Y-M and Y-E is not controlled | Each endpoint has its own rule, tested at 0.025 per side; there is no combined claim. | design |
| 8 | Replicate-level chips are not independent samples | The input is the replicate-level matrix (536 chips on net4, in 321 replicate groups). The experiment jackknife deletes whole replicate groups. | design |
| 9 | SE_B may be somewhat small | The deciding SE_B is at replicate-group level; 33 of the 59 series contain several groups. The series-level supplementary SD is larger (section 3). | design; `supplementary_series_jackknife.csv` |
| 10 | PIDC has a tied block at the cap | A 2-entry tied block straddles rank 100,000 on net1 and net4 (`n_tied_at_cap` = 2), none on net3; the design assumed none. All six other methods have 0 on every net, including attention and CLR. The pre-registered handling, `drop_tied`, removed the 2 entries from PIDC's list only (`K_listed` PIDC: 87,558 net1, 12,440 net4). This affects PIDC only; the attention-versus-CLR verdict is unaffected. Table below. | `cap_boundary.csv`, `point_estimates.csv`; log 7.1, 7.3 |
| 11 | PIDC scores assume non-negative weights | PIDC scores are used with NaN set to 0 and assuming non-negative weights; the pipeline does not check non-negativity. | `docs/dream5_truth_stage_checklist.md` item 5 |
| 12 | PIDC Gamma fallback | `results/dream5/net4/pidc.log` contains 23,790 lines of the form "Gamma distribution failed for G_ and G_; used normal instead" (arithmetic: 23,790 of 17,698,275 edges, 0.134%). The net1 and net3 logs contain no such line (each holds one summary line). PIDC only. | `results/dream5/net{1,3,4}/pidc.log` (local; `results/` is gitignored) |
| 13 | Job-script settings are not versioned | The refit and score jobs ran from `jobs_pc/` scripts (thread counts, worker counts, memory) that are not under version control. `results/dream5/PROVENANCE.txt` records the thread settings (refit jobs: OMP_NUM_THREADS 8, OPENBLAS_NUM_THREADS 4, 8 workers; truth driver default 16 workers) and the environment; the truth-stage job commands are recorded in log section 7. The group-level refit code commit is inferred from commit times (no commit id recorded by the jobs); the series-level refit working tree was not checked clean at launch; the installed `clr` build is not tied to a commit hash (017e02b inferred). A re-computation of group and series refit b = 7 with 16 workers reproduced all 8 of 8 SHA-256 values, and all 100 + 100 refits matched their committed hashes at the truth stage. | `PROVENANCE.txt`; checklist item 6; log 7.3 |
| 14 | Y-E definition set before any score of ours | The official scorer has no cap or tie rule in code; the Y-E definition (gold-derived universe, cap then universe, TF-TF pairs both directions, `drop_tied`, expected ties, restrict after cap, official trapezoid AUPR/(1 - 1/P) with expected-fraction tail) was adopted on 2026-11-06 and logged in the design before any score existed; no method, endpoint, margin or seed changed. Y-E AUPR is therefore the trapezoid measure and differs slightly from step-wise average precision (`aupr_stepwise` is kept in `point_estimates.csv`; net4 attention 0.02164 trapezoid, 0.02171 step-wise). | design change log 2026-11-06; `point_estimates.csv` |
| 15 | Scorer-gate sheet precision | The organisers' sheet stores 3 decimals (a few cells 4). Differences of up to 5.79e-4 to the published table are consistent with that storage (2 of 324 metric values exceed 5e-4, all within 1e-3); the pass rule was 1e-3, so the gate shows agreement to about the third decimal. | `gate_report.md`; log section 6 |
| 16 | Not implemented in the truth-stage pipeline | The replicate-averaged scoring option, the all-genes-regulator GENIE3 and GRNBoost2 module readouts (score files hashed but not scored), and the exploratory variants of the design. Only the pre-registered contrast on Y-E and Y-M is in `verdict.csv`; the other tables are point estimates. | checklist items 8, 9 |
| 17 | Methods use expression only; decoys kept | Perturbation annotations are not used; the official files are used as provided, including decoy genes (which cannot enter the scorer's pair universe). t = 9 and all settings were transferred or fixed, not tuned on DREAM5. | design |

Cap-boundary check (`cap_boundary.csv` of each net; the design expected no tied block at rank 100,000):

| Net | Candidate pairs | Methods with a tied block at the cap (`n_tied_at_cap`) |
|---|---|---|
| 1 | 320,190 | PIDC (2) |
| 3 | 1,506,340 | none (0 for all 7 methods) |
| 4 | 1,981,017 | PIDC (2) |

## 7. What this does and does not support

**Supported literally by the tables above, in the design's wording.**

- On the independent yeast network (net4), under the pre-registered rule, Y-M (regulon co-membership, CLR attention versus CLR) is Equivalent within 10%, hence not Superior. The fixed wording applies: *the co-membership headline did not replicate in yeast*, with the estimate and interval: d = 0.0817 (AUPR ratio 1.058), 95% interval 0.0260 to 0.1374.
- On the same network, Y-E (official directed edge AUPR) is Inconclusive: d = 0.0599 (AUPR ratio 1.042), 95% interval -0.0928 to 0.2125. No claim.
- Both endpoints are to be cited together.

**Not supported.**

- Any claim of a module gain with a direct-edge cost in yeast (it requires Y-M Superior), or of an improvement on both (likewise).
- Any statement that Y-E differs, or does not differ, between the two methods: the interval -0.0928 to 0.2125 contains 0 and +delta.
- Any claim of "no effect" for either endpoint; the wording is "not confirmed", with estimate and interval.
- Support or replication from net3 (same lineage) or net1 (simulated), from pooling networks, or from results within +-10% (design, "What does not count").
- Any combined claim across Y-M and Y-E, or any control of the familywise error.
- Any statement that net4 supports a ranking of methods: all seven methods are near chance in absolute terms on net4 (disclosure item 5), and the other-method values are descriptive point estimates.
- Any comparison with the challenge's published ranking beyond the descriptive table of section 5.

**PI decision needed: framing in the paper.** (Stub. This report does not propose wording for the paper or changes to the paper's claims; how the net4 result, the two-endpoint rule and the section 6 disclosures enter the manuscript is left to the PI.)
