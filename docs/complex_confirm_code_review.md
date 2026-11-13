# Code review: yeast complex estimation (complex_confirm), design rev 2

Reviewer: Claude (Opus), launched as the design's protocol step 4 reviewer, 13 Nov 2026. Branch `accuracy-review`, HEAD 5e999e5. The code under review is byte-identical to the code that ran: there is no diff in `analysis/complex_confirm/`, `analysis/dream5_*.R`, `analysis/regulon*.R` between eabb2ee (the 6020 run) and HEAD, and the three libraries are unchanged since d8fdb58.

## Read this first

1. **No bug found that could change a reported number.** Within this design, the primary row P, rows C1 and C2, their gene, series and group intervals, and the fixed sentences are computed as the design specifies. Every check I could run on synthetic data passed (list below).
2. **This review is not fully blind (process disclosure).** I ran `grep -n` on the design file to find the line numbers of the change-log entries, so that I could skip the truth-stage entry. The grep printed that whole line (design line 525, "Truth stage run"). It shows the row-P point estimate and its total and gene-only intervals. This happened early, before I read any code. After that I opened none of the blinded files (`estimate*.csv`, `tables.md`, `net4/truth/*`, `net4/expjack/*`, STATUS estimates). Every finding below comes from the code, the design, the allowed gate and label-free outputs, and synthetic checks. Still, a strict reader should treat this as a review by someone who has seen the headline estimate.
3. **This review comes after the run, not before it.** Protocol step 4 puts the code review before steps 5–6. The executing agent disclosed that no independent pre-run review was done (change log, implementation note iv). This review covers the code that ran, after the fact.

## Summary

- The code implements design rev 2 faithfully on every point that affects a reported number: AUPR definition, gene / series / group jackknives and seeds, pairing, total-interval combination, sentence thresholds, a3r readout, pair-universe exclusions, hub stratifier, and the write-once guard.
- The findings are robustness and disclosure items: two process items (late review, partial blinding) and some minor guard and gate-coverage gaps. None changes a reported number in this run.
- The things I could not verify need the blinded outputs. Examples: replicate counts actually written, the per-replicate hash results, the git head recorded in the truth-stage outputs.

## Findings

| ID | Severity | File:line | What | Could it change a reported number? |
|---|---|---|---|---|
| F1 | major (process) | design, Protocol step 4; change log impl. note (iv) | The independent code review was not run before the truth stage. This review is post hoc. | No. It is a protocol deviation and should be reported as one: "code review after the run". |
| F2 | major (process) | (this review) | The reviewer accidentally saw the row-P estimate and intervals in the change log (see "Read this first"). | No, but the review cannot be called blind. |
| F3 | minor | `analysis/complex_confirm/gate_a_k562.R`, `gate_a_m3d.R`; `cc_lib.R:221-234`; `cc_run.R:282-291` | Ga parts 2–3 reproduce stored K562 and M3D values through the **old** drivers' scoring code (`ideas47_lib.R::score_readout`, `stack_lib.R::sk_make_scorer`). They do not exercise the new wrappers that produced the yeast numbers: `cc_score_truth`, `cc_interval`, and the summarize arithmetic. These wrappers were covered only by synthetic smoke tests. They share `ap_fast` and `sk_a3r` with the gated path. | Not in this run: I checked `cc_score_truth` against an independent loop of the `sk_make_scorer` / `jack_dataset` kind, serial and forked, with max diff 0. The gap is in gate coverage, not in the result. A cheap gate for later runs: re-score the stored M3D truths through `cc_score_truth`. |
| F4 | minor | `cc_run.R:22-39` | Guard Ga5 (no earlier a3r score on net4). (a) It reads only the first line of text files in directories whose *dirname* matches net4/dream5. (b) Binary files (`.rds`, 42 found) are checked by name only, so an `.rds` holding a3r scores under a neutral name would pass. (c) A file whose *filename* contains net4/dream5 but sits elsewhere gets only the name check. The design asks for "header or column names" in dream5/net4 directories. | No. It concerns only the independence claim. The CSV states what was searched, and the design's own scope (`results/` only) is implemented. |
| F5 | minor | `cc_run.R:230-238, 272-277` | In `--collect`, `replicates_*.csv` is written before the missing-replicate and hash checks. Summarize checks that the series set is complete (`series_ok`), but it does not check that there are 100 group replicates or that all hashes matched. | No. Every scored replicate is hash-verified first (the lane stops at `cc_run.R:259` before scoring on a mismatch), and the scheduler's `# after` requires rc 0 (`tools/scheduler.sh:8`). A failed collect therefore cannot be followed by 6040. Only a manual 6040 after a failed 6030 could use partial group replicates. |
| F6 | minor | `cc_run.R:315` | Write-once uses `file.link`, which is atomic and refuses an existing target (verified). But `estimate.csv` is linked before `estimate_final.csv`. If the second link failed, `estimate.csv` would exist alone, and the guard would block a repair. | No. |
| F7 | minor | `cc_run.R:327` | `tables.md` is written with `writeLines`, with no exists-guard. Write-once covers only the two CSVs. | No. `tables.md` is derived from `est`. |
| F8 | note | `cc_run.R:118-121` | Gb5 rows: (a) the two hash rows are hard-coded `1 / TRUE`; the real check is `cc_truth_hash_guard()` one line earlier, which stops on a mismatch. (b) n_eff is read from the frozen counts file, not recomputed. (c) "merged complexes with positives" counts sets with ≥ 2 genes, not sets with ≥ 1 non-paralog positive pair. The structure rows (positives and universe for 6 variants and 2 strata) are recomputed and match. | No. |
| F9 | note | `cc_lib.R:239-241`, `235` | The design says the percentile column is "centred on d". The code computes d + quantile(rep − mean(rep), type 7). The skewness is the moment estimator g1 (biased, no small-sample correction). Disclosed in change-log note (iii). | Descriptive columns only. Another reading (raw percentiles of the replicates) would shift both bounds by mean(rep) − d. The direction depends on the data. |
| F10 | note | `cc_run.R:162`, `cc_lib.R:217-219` | Each size stratum has its own 500 halves, drawn from the primary truth's 1,602 genes in table order. They are not the same halves as `main`. This follows the design's table-order rule. Within a stratum the three readouts are paired. | No (Monte Carlo only). |
| F11 | note | `robust_z_lib.R:134-139`; `cc_run.R:169-171` | The hub strata are deciles of H_i + H_j over **all** 17.7 M upper-triangle pairs of the 5,950 genes, including decoys and pairs outside the universe. The within-stratum percentile is relative to those pairs. This is `rz_hub_stratify` as the design names it, and as the baselines used it. | No (matches the design). Secondary lines only. |
| F12 | note | `cc_lib.R:142` | Smoke-test synthetic truths cap complex size at 40, whereas the design says "the real complex-size distribution". | No (smoke tests only). |
| F13 | note | `jobs/6030_cc_net4_expjack.sh` | The 6030 resource declaration differs from the design table (5 lanes, 30 GB, 40 threads instead of 1 process, 16 GB, 8). It is logged as an implementation note. The refits are deterministic and hash-checked. | No. |
| F14 | note | `cc_lib.R:54`, `cc_run.R:14` | `--t` is validated but never used. The depth comes from `D5$t`, which `stopifnot(D5$t == CC$t)` (`cc_lib.R:26`) ties to 9. | No. |
| F15 | note | (attribution) | The commit trailer requested for this file names Sonnet. The reviewer was Opus, so the commit carries an accurate trailer. | No. |

## Verified OK

**Code reading (design ↔ code):**

- **Operator and readouts** (`cc_lib.R:113-127`):
  - CLR and attention follow `d5_m_clr_attention` exactly (hg bins, no transform, normal + Stouffer; keff 10, cap 50, α 0.5; P⁹ by repeated product).
  - a3r is `sk_a3r` on **the same `Pt` object**. Ga 4/4b (allowed output) shows max |diff| 0 against the DREAM5 caches.
  - `stack_lib.R`, `ideas47_lib.R` and `robust_z_lib.R` are git-blob-identical to d8fdb58.
- **Hub stratifier = MI row mean.** `(rowSums(M) − diag(M))/(G − 1)` on `f$mi`, the raw MI active binding (`R/clr_attention.R:515`), over all 5,950 genes. It is applied identically to all three readouts.
- **Pair universe.** All pairs among the genes of merged complexes, built by `comembership_pairs`, minus the 49 paralog pairs (scope `all`). The labels are "share a merged complex". The same `p$i, p$j, p$label` is used for all three readouts in the pooled score, every gene replicate, every complex replicate and every refit.
- **Exclusions.**
  - `main_xsub_excl` removes `all` + `main_xsub_excl` keys.
  - `main_paralog_incl` removes nothing.
  - Structure only, from the frozen files: the 5,334 cross-subunit pairs are all in the main universe and none is a positive; `main_paralog_incl` has exactly 49 more pairs than `main`.
- **Size strata** (`cc_stratum`). Universe = main universe minus the pairs that are co-members only of complexes outside the stratum; positives = co-members of a stratum complex.
  - Structure only: in all three strata the negatives are identical to the main negatives, and every main positive is a positive in at least one stratum.
  - Strata: small 198 complexes / 1,252 positives; mid 51 / 2,781; large 13 / 9,531. Small and large equal the frozen counts.
  - Boundaries: 6 → small, 7 and 19 → mid, 20 → large.
- **Gene jackknife.**
  - Seeding: `set.seed(20261114)`, then the halves for all 11 truths in design table order (`CC$order`, ending with `dream5_ym`).
  - Each half is `sample(genes_in, floor(n/2))`, and a replicate scores the pairs with both genes kept. This is identical to `clr_logpt_test.R::jack_dataset`.
  - There are 500 replicates, and SE = sd of the paired differences (factor 1).
- **Complex jackknife.** `set.seed(20261115)`, 500 halves of the 262 merged complexes. The universe is rebuilt from the kept complexes minus the paralog pairs; the scores are paired.
- **Experiment jackknives.**
  - Group halves: `set.seed(20263000 + b)`, `sort(sample.int(321, 160))`. Series halves: `set.seed(20264000 + b)`, `sort(sample.int(59, 29))`; all chips of each kept series are kept. Both match `d5_group_half` / `d5_series_half`, and I verified them directly.
  - Each refit's CLR and attention SHA-256 must equal the DREAM5 manifest row (`net == 4`, `rep == b`, `method`), or the lane stops before scoring.
  - Because a3r comes from the same fit, a hash match proves that a3r was refitted on the same chips.
  - The truth structure in 6030 is md5-checked against the truth stage.
- **Totals.**
  - Total intervals: series total = √(SE_J² + SE_series²) and group total = √(SE_J² + SE_group²), each as d ± 1.96·total.
  - The sentences use the series total for row P; if the series set is incomplete, they fall back to the group total and set the flag.
  - Side-by-side lines are given for gene-only and gene ⊕ group.
- **Pairing in the experiment jackknives.** The differences are paired by replicate id (`order(rep)`), which is robust to row order (tested).
- **No label leakage.**
  - Operators, readouts and the hub stratifier depend only on expression (and MI).
  - The truth stage zeroes `nd$X` before building truths.
  - Labels enter only in `ap_fast`, `regulon_ap_table` and `regulon_coherence`, and in the choice of jackknife gene sets, which the design specifies.
- **Guards.** The truth stage refuses unless all of these hold:
  - the design is tracked, clean, pushed and carries `Status: FIXED`;
  - the code files are clean and pushed, `dream5_score_lib.R` included;
  - the truth hashes and the 6010 cache manifest match;
  - the settings are the preregistered ones;
  - no truth-stage output or estimate file exists.
- **Forked-worker death.** If a worker is OOM-killed, `mclapply` returns NULL. The `length(x) != ncol` and `!= 5` checks then stop the job, so a killed worker cannot leave a silently shortened replicate set.

**Synthetic unit checks (scripts in `/tmp/claude-1000/-home-boris-Projects-ClrAttention/4451c478-b0bc-4b2f-84d1-11ecf74c2d08/scratchpad/t1.R` and `t2.R`; all passed):**

- **AUPR.** `ap_fast` equals `pr_summary()["aupr"]` (the validated scorer; max diff 6e-17 over 200 random tied problems).
  - It also equals a brute-force AP that credits each tie block with the precision at its end (1e-16).
  - Without ties it equals standard step AP, with no interpolation (diff 0).
  - Hand cases: constant score → base rate; NA ranks last.
- **a3r.** `sk_a3r` equals the `docs/plans/corum_analysis.md` reference code (diff 0) and a brute-force per-row loop over reachable entries only (mean, SD with n − 1; 4e-16). It is exactly symmetric, and it differs from a 1e-300-floor variant, which shows that unreachable entries really are excluded from the row statistics.
- **Gene jackknife.** `cc_score_truth` matches an independent loop. Results are identical with 1 and 2 forked cores.
- **Halves and pairing.** The first half equals `sample(genes_in, floor(n/2))` right after the seed is set. (a3r − att) = (a3r − CLR) − (att − CLR) holds replicate by replicate, so the three contrasts share draws and universe.
- **Interval and skewness.** `cc_interval`: SE = sd(n − 1), d ± 1.96 SE, percentile and g1 formulas.
- **Sentences** (`cc_sentences`), at and 1e-12 either side of 0, +0.005 and +0.027:
  - L = 0 or U = 0 → "includes no difference";
  - L = 0.005 or U = 0.005 → "includes the smallest effect";
  - L = 0.027 or U = 0.027 → "includes the E. coli estimate";
  - strict inequalities otherwise, as the design table states.
- **Paralog exclusion and strata** on toy overlapping complexes.
- **Hub stratifier.** The row mean excludes the diagonal; the stratified matrix is symmetric, with values in (0, 1].
- **Write-once.** `file.link` onto an existing `estimate.csv` returns FALSE and leaves the file unchanged.

**Allowed outputs read:**
- `gate_a/equivalence_gate.csv` (43 checks pass);
- `net4/labelfree/equivalence_net4.csv` (Ga4/4b max |diff| 0);
- `feasibility_gate.csv` (all Gb pass; frozen structure reproduced);
- `score_manifest.tsv`, `refit_timing.csv` (CLR and attention hashes match DREAM5 group rep 1);
- `guard_no_prior_a3r.csv`;
- the frozen truth files (both SHA-256 values recomputed and matching).

## Not verifiable without opening blinded outputs

- That 6020 actually wrote 500 gene and 500 complex replicates per truth, and that `truth_counts.csv` structure md5s match what 6030 rebuilt. The code enforces both; the counts are in `net4/truth/*`.
- That `refit_hash_check.csv` has 200 rows, all TRUE, and that `estimate.csv` has `n_series_reps = n_group_reps = 100` and no fallback flag. The code makes each a precondition for a scored replicate, for collect success and for the series basis.
- The git head recorded in the 6020_retry and 6040 outputs (expected: eabb2ee).
- Whether the replicate distributions are degenerate in some way (for example many identical replicate AUPRs) that would make the SD or percentile columns unreliable.
- Any check on the numbers themselves: plausibility of AUPR levels, sign consistency between the gene and experiment replicates, or how the a3r gain splits across strata.
