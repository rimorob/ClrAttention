# GRNBoost2 diagnostic: does arboreto's broken early stopping explain last place?

*27 October 2026. Diagnostic read-out and rerun plan (job 0894 plus rescoring with the main evaluation). It does not edit the paper or tables. Every number below comes from `results/grnboost2_check/` (gitignored; committed copies are in `docs/grnboost2_diagnostic/`).*

> **Order of events (disclose).** We reported GRNBoost2 as a comparator first: in `baselines.R` runs from 16 October, in Table 2 / `sota_table.md` (21 October), in Supplementary Tables 9 and 10, and on the Collins checkpoint. The early-stopping fault was found afterwards, on 22 October. This diagnostic was designed after that and scored after that. The paper already marks every GRNBoost2 result as "not interpreted until the wrapper check" [PENDING].

## 1. Bottom line

1. **Yes, mostly. The broken early stopping explains most of GRNBoost2's deficit, but not all of it.**
   - The fix raises GRNBoost2 by about 35 seed-SDs on M3D. On M3D SC it goes from 1.18 to 1.36 times the base rate (Table 2 metric). On M3D all evidence it goes from 1.13 to 1.26.
   - On BEELINE it gains +0.12 log2 (hESC_1000) and +0.27 log2 (mHSC-E_500) against CLR.
   - It is no longer last: it passes |r| and the hub score on M3D and |r| on BEELINE. On hESC_1000 it reaches CLR's level.
   - It still trails CLR, PIDC, GENIE3 and attention on M3D. It also trails GENIE3 and PIDC on both BEELINE cases.
2. **Our wrapper is not at fault.** arboreto's own `grnboost2()` (dask) reproduces our networks: its scores match the current variant to within 0.002 on every M3D score, seed by seed. Seed variability is small (SD of the ratio 0.002-0.013), so neither the seed nor the wrapper explains last place.
3. **Legacy early stopping and no early stopping give the same accuracy.** On M3D, no_es scores within ±0.02 of the legacy_es mean on every Table 2 cell; on BEELINE it is within 0.024 log2. legacy_es costs about a tenth of the CPU time. Over-fitting is therefore not the issue; too few trees was.
4. **No conclusion about attention changes on M3D.**
   - The best method and the bold or tie pattern of Table 2 are unchanged in all four M3D cells.
   - GRNBoost2's row position does change: from 6th or 7th of 7 to 5th of 7.
   - Under hub control it becomes indistinguishable from CLR and PIDC on SC regulons.
5. **Claims that change:**
   - "GRNBoost2 ranks last everywhere";
   - "below CLR in 7 of 7 BEELINE datasets (P = 0.016)", which is likely to change;
   - the per-regulon lead of attention over GRNBoost2 (Holm P 0.004), which is at risk because the mean gap halves;
   - the sparse-network footnote, since 54% of M3D pairs are now non-zero against 6%.
6. **Recommendation:** rerun GRNBoost2 everywhere with **legacy early stopping**, which is arboreto's algorithm as published. It goes in every table and figure, and the current (broken) runs are disclosed in the supplement. no_es stays a sensitivity analysis, already in hand on M3D and two BEELINE cases. Fix this rule in writing before running EcoMAC, K562 or the other 12 BEELINE cases (section 6). The compute is about 15-25 h of wall time at 24 threads (section 5).

## 2. The fault (recap)

- `tools/grnboost2.py` runs arboreto 0.1.6's own per-target routine (`infer_partial_network`, GBM, `SGBM_KWARGS`, early-stopping window 25), and it does so faithfully. The `dask` variant (arboreto's `grnboost2()` itself) gives the same networks (table in §3), so the wrapper is not at fault.
- arboreto's `EarlyStopMonitor` stops a target's boosting when the mean of `oob_improvement_` over the last 25 stages is below 0. Before scikit-learn 1.3, `oob_improvement_[i]` was the out-of-bag loss drop of stage i, measured on the same OOB sample before and after the stage. From 1.3 on (we run 1.9.1) it is `oob_scores_[i-1] − oob_scores_[i]`, measured on two different random 10% OOB subsets. The 25-stage window mean therefore telescopes to the difference of two noisy subset losses, and the monitor stops at an essentially random early stage.
- Effect on network size (trees per target, median and IQR, from the variant logs):
  - **hESC_1000:** legacy_es 212 (171-262), no_es 5,000.
  - **M3D:** legacy_es 494 (451-536), no_es 5,000.
  - **current:** mostly 25-36 trees (spot check in `docs/paper_plan.md`). The current networks keep only 3.5% of M3D pairs non-zero, against 35% for legacy_es.
- Variants (`tools/grnboost2_variants.py`; 3 seeds unless stated):
  - **current:** `tools/grnboost2.py` as used for every reported result.
  - **dask:** arboreto's `grnboost2()` itself.
  - **legacy_es:** the arboreto profile with the pre-1.3 (same-sample) OOB improvement. This is the algorithm as published.
  - **no_es:** no early stopping (5,000 stages). On M3D only seed 1 was run, because it costs about 230 CPU-hours per seed.

## 3. Scores with the project's own evaluation

No new truth or metric is used. `analysis/grnboost2_check_score.R` (new; nice 19, one thread, about 15 min, peak under 3 GB) works as follows:

- **M3D:** it evaluates `analysis/baselines.R`'s own definitions from its source. These cover the data, the RegulonDB SC and all-evidence co-membership truths, the delete-half gene jackknife splits (seed 20261016), `evaluate()`, `sym_max()`, `hub_stratify()` and the walk's stationary distribution.
  - The script stops unless the reported GRNBoost2 edges reproduce the reported AUPR, the hub-stratified AUPR and the per-regulon AP. **These gates passed:** raw AUPR SC 0.111273801 and all 0.136021394; hub-stratified 0.120854297 and 0.141579650; every per-regulon AP; the per-regulon mean ΔAP vs attention also reproduces Supplementary Table 9 (−0.0269 SC, −0.0240 all).
  - Each variant takes GRNBoost2's place in its Table 2 row. Paired jackknife differences use the stored replicates of the other methods, and ties use `sota_table.R`'s Bonferroni rule.
- **BEELINE:** the diagnostic's `metrics.csv` already used `beeline_baselines.R`'s tie-aware `score_truth` with directed GRNBoost2 weights, as in the main run. Its `grnboost2_reported` rows equal `results/beeline_baselines/metrics.csv` (a gate in the script). Per case, the summary is Supplementary Note 13's statistic: the mean over BEELINE's primary truths of log2(AUPRC ratio / CLR's).

### 3a. E. coli M3D, Table 2 metric (pooled co-membership AUPR / base rate, sym-max)

Ratio is AUPR / base rate (base rate 0.094 for SC, 0.120 for all). Seeds 1-3 are shown as mean ± SD (range); no_es has seed 1 only. "Reported" is the Table 2 value (seed 20261016). Reference values come from the stored run.

| Scoring | Truth | reported | current | dask | **legacy_es** | no_es (s1) | Attention | GENIE3 | PIDC | CLR | \|r\| | Hub only |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| raw | SC | 1.18 | 1.181 ± 0.005 (1.177-1.186) | 1.181 ± 0.005 | **1.361 ± 0.002** (1.359-1.363) | 1.369 | 1.68 | 1.51 | 1.45 | 1.44 | 1.23 | 1.08 |
| raw | all | 1.13 | 1.129 ± 0.002 (1.127-1.131) | 1.129 ± 0.002 | **1.258 ± 0.001** (1.257-1.258) | 1.251 | 1.44 | 1.40 | 1.40 | 1.37 | 1.23 | 1.22 |
| hub-controlled | SC | 1.28 | 1.293 ± 0.007 | 1.293 ± 0.008 | **1.373 ± 0.003** | 1.367 | 1.56 | 1.47 | 1.40 | 1.40 | 1.22 | |
| hub-controlled | all | 1.18 | 1.202 ± 0.013 | 1.202 ± 0.013 | **1.276 ± 0.001** | 1.256 | 1.35 | 1.37 | 1.34 | 1.33 | 1.21 | |

The table below gives the place in the Table 2 row and the paired jackknife differences (variant − method, AUPR, ±1.96 SD, the same 100 delete-half splits). The legacy_es row is seed 1; seeds 2 and 3 differ by 0.0006 or less.

| Scoring, truth | Place in row: current → legacy_es | legacy_es − CLR | legacy_es − PIDC | legacy_es − GENIE3 | legacy_es − \|r\| | legacy_es − attention |
|---|---|---|---|---|---|---|
| raw SC | 6/7 → 5/7 | −0.0074 [−0.0114, −0.0033] | −0.0083 [−0.0125, −0.0041] | −0.0134 | +0.0127 [+0.0078, +0.0176] | −0.0294 |
| raw all | 7/7 → 5/7 | −0.0132 [−0.0170, −0.0095] | −0.0166 [−0.0204, −0.0128] | −0.0171 | +0.0037 [−0.0011, +0.0084] | −0.0218 |
| hub SC | 5/6 → 5/6 | −0.0027 [−0.0062, +0.0008] | −0.0023 [−0.0056, +0.0011] | −0.0088 | +0.0145 [+0.0098, +0.0192] | −0.0175 |
| hub all | 6/6 (seed 1: 5/6) → 5/6 | −0.0069 [−0.0104, −0.0034] | −0.0076 [−0.0107, −0.0045] | −0.0109 | +0.0083 [+0.0036, +0.0130] | −0.0089 |

- **Best method in the row:** attention (raw; hub-controlled SC) or GENIE3 (hub-controlled all), for every variant. GRNBoost2 is never tied with the best, so the bold and dagger pattern of the four M3D rows does not change.
- **Per regulon** (Supplementary Table 9 statistic: mean ΔAP of GRNBoost2 − attention over the regulons with ≥ 10 positive pairs; point estimates only, no interval or Holm here):
  - **SC (123 regulons):** reported −0.0269; current −0.0276 (−0.0296 to −0.0244); legacy_es **−0.0157** (−0.0163 to −0.0148); no_es −0.0149.
  - **All evidence (157 regulons):** reported −0.0240; current −0.0245 (−0.0271 to −0.0230); legacy_es **−0.0121** (−0.0138 to −0.0106); no_es −0.0116.
  - **Share of regulons where GRNBoost2 is ahead:** 22-32% (current) against 39-46% (legacy_es).
- **TF-node AUPR ratio** (`baselines.R`'s secondary TF → target benchmark, sym-max):
  - **SC:** 3.5-3.8 (current) → 5.1-5.4 (legacy_es); GENIE3 6.82, |r| 3.94.
  - **All evidence:** 2.6-2.8 → 3.7-3.9.
  - The job 0894 log reports the same on directed edges: AUPRC ratio 3.31 → 4.92, no_es 4.90.
- **Non-zero pairs after sym-max:** reported 6.0%; current 6.2-6.3%; legacy_es 54.4%; no_es 99.98%.

### 3b. BEELINE (two of the 14 cases), Supplementary Note 13 statistic

The statistic is the mean over primary truths (cell-type ChIP, non-specific ChIP, STRING; three truths in both cases) of log2(AUPRC ratio / CLR's). It is reported as mean ± SD over 3 seeds, with the range in brackets. The GRNBoost2 row is the reported run (seed 20261021).

| Case | current | dask | **legacy_es** | no_es | reported GRNBoost2 | GENIE3 | PIDC | CLR-NS | \|r\| |
|---|---|---|---|---|---|---|---|---|---|
| hESC_1000 | −0.128 ± 0.011 [−0.137, −0.116] | −0.128 ± 0.010 | **−0.009 ± 0.016** [−0.023, +0.008] | −0.013 ± 0.016 [−0.025, +0.005] | −0.127 | +0.036 | +0.102 | +0.009 | −0.139 |
| mHSC-E_500 | −0.559 ± 0.036 [−0.598, −0.527] | −0.559 ± 0.037 | **−0.293 ± 0.012** [−0.307, −0.285] | −0.317 ± 0.005 [−0.322, −0.312] | −0.571 | −0.261 | −0.192 | +0.086 | −0.411 |

Mean AUPRC ratio over the same truths:

| Case | current | legacy_es | no_es | GENIE3 | PIDC | CLR | \|r\| |
|---|---|---|---|---|---|---|---|
| hESC_1000 | 1.14 (1.14-1.15) | 1.28 (1.26-1.29) | 1.28 (1.26-1.29) | 1.35 | 1.41 | 1.27 | 1.13 |
| mHSC-E_500 | 2.25 (2.16-2.32) | 2.83 (2.78-2.86) | 2.75 (2.74-2.75) | 2.93 | 3.28 | 4.05 | 2.59 |

- **Early precision** (log2 EPR vs CLR):
  - **hESC_1000:** current −0.24 ± 0.10 → legacy_es −0.01 ± 0.04.
  - **mHSC-E_500:** −0.43 ± 0.05 → −0.31 ± 0.02.
- **Gain from the fix:** +0.12 log2 (hESC) and +0.27 log2 (mHSC-E). That is 8-25 times the seed SD.
- **Rank change:** GRNBoost2 moves from last (at or below |r|) to level with CLR on hESC_1000, and to just below GENIE3 on mHSC-E_500. It stays below PIDC and GENIE3 on hESC.
- **Consistency with BEELINE:** this matches BEELINE's placement of GRNBoost2 "near GENIE3".

## 4. Interpretation: rankings and claims affected

**Is the effect large?** Yes, compared with seed noise and with the gaps to the next methods. The fix closes about 40% of GRNBoost2's gap to attention on M3D SC (−0.047 → −0.029 AUPR) and half (mHSC-E_500) to nearly all (hESC_1000) of its gap to CLR on BEELINE. It is not large enough to change who wins: on no M3D cell does any GRNBoost2 variant come within the tie rule of the best method.

**Table 2 / `sota_table.md`.**
- **M3D, all four rows:**
  - The GRNBoost2 values change (1.18 → 1.36, 1.13 → 1.26, 1.28 → 1.37, 1.18 → 1.28).
  - GRNBoost2's position changes from last or second last to 5th of 7, ahead of |r| and the hub score.
  - The best method and the bold or dagger pattern do not change.
  - Under hub control (SC), GRNBoost2 is no longer distinguishable from CLR or PIDC.
- **EcoMAC (4 rows) and K562 (5 rows):** not tested here. Expect a similar rise, because the same code and fault apply.
  - On EcoMAC a rise of about 0.15 still leaves GRNBoost2 below attention (1.68 / 1.47).
  - On K562 the gaps are small. Several cells could change their bold pattern: a fixed GRNBoost2 could join the ties on DoRothEA, CORUM (now 2.18 against 6.3-7.0) or knockdown response (1.37 against 2.35-2.45). A K562 cell could even get a new best, which would change the other methods' daggers. This is unknown until the rerun.
- **Footnote:** "GRNBoost2 returns a sparse network (87-93% of pairs tied at zero) ... not like-for-like" becomes wrong. The sparsity was an artefact of the fault.

**Paper claims.**
- *"GRNBoost2 ranks last everywhere in our hands, bulk and single cell, whereas BEELINE placed it near GENIE3"* (`paper.md` Results, Supplementary Note 13). This is false after the fix. The discrepancy with BEELINE is explained by the early-stopping fault, so the claim should be replaced by the explanation, not by a new ranking claim.
- *BEELINE: "our GRNBoost2 runs score below [CLR] in 7 of 7 datasets (−0.35; P = 0.016)"* (paper and Supplementary Table 10). This is likely to change:
  - on hESC_1000 the fixed GRNBoost2 is at CLR's level (seed range −0.023 to +0.008);
  - mHSC-E gains +0.27.
  - The 0/7 count, P = 0.016 and the median AUPRC ratio of 1.26 need the full 14-case rerun.
  - GRNBoost2 will probably move above |r| (−0.23) and stay below GENIE3 and PIDC.
- *Per regulon, attention's mean AP exceeds GRNBoost2's in 4 of 4 cells (Holm P 0.004)* (paper Results, Supplementary Table 9 reading 2, `abstract_points.md`, the Collins checkpoint). This is at risk: on M3D the mean ΔAP halves (−0.027 → −0.016 SC, −0.024 → −0.012 all).
  - Whether it still passes both the regulon interval and the gene interval is unknown until `macro_aupr.R` is rerun.
  - **Side effect:** the E. coli Holm family has a fixed size of 40 and includes the GRNBoost2 tests. Holm is step-down, so changing those p-values can also move the adjusted p-values of other comparators (|r|, PIDC) in the same family.
- **Unaffected:**
  - the headline (attention vs CLR);
  - "attention has the highest pooled AUPR in all four E. coli cells" (holds on M3D; EcoMAC still needs its rerun);
  - per-regulon parity with CLR and GENIE3;
  - "CLR comparable to GENIE3 and PIDC on BEELINE";
  - every preregistered test (CLR(log Pᵗ), keff, scfm primary). None of them uses GRNBoost2 as a comparator.
- **Methods:** "GRNBoost2 (arboreto) ... an audit is running [PENDING]" must state the versions (arboreto 0.1.6, scikit-learn 1.9.1), the fault, and the fix.

## 5. Recommendation: how to rerun GRNBoost2 fairly

**Choice: legacy early stopping (legacy_es) as the GRNBoost2 of record, everywhere.** No-early-stopping (no_es) is a sensitivity analysis only; it is not rerun.
- **Fidelity.** legacy_es is the GRNBoost2 that Moerman et al. (2019) published and that BEELINE benchmarked: arboreto's profile and stopping rule, with the OOB improvement it was written for. no_es is a different method (5,000 stages at learning rate 0.01, no regularization by stopping).
- **Equal accuracy, a tenth of the cost.** On every cell scored, no_es is within seed noise or ±0.02 of legacy_es.
- **Do not choose between them by score.** Both are now seen on three cases. The choice above rests on fidelity and cost, and should be fixed before any new case is run.

**Steps (all through `jobs_pc/sched`):**
1. **Make it the default in `tools/grnboost2.py`:**
   - add `--early_stop legacy|sklearn|none`, default `legacy`, reusing `LegacyOOBMonitor` from `tools/grnboost2_variants.py`;
   - print the mode and the trees-per-target quantiles in the log line.
   - Keep `grnboost2_variants.py` as the diagnostic.
2. **Equivalence gate (cheap, minutes):** in a separate Python 3.11 venv with scikit-learn 1.2.2 and arboreto 0.1.6 (scikit-learn < 1.3 has no Python 3.12 wheels), run stock arboreto on one small case (hESC_500 or mDC_500, 3 seeds). Compare it with legacy_es: trees-per-target distribution and AUPRC ratio within the seed SD. This shows that our monitor reproduces the published behaviour. If it fails, use the pinned venv itself for the rerun.
3. **Move aside every cached GRNBoost2 output** (`grnboost2_edges.tsv`, `grnboost2_scores.rds`) to `*_sklearn13.*` so that no script reuses them. `baselines.R` reuses `<method>_scores.rds` when it exists.
4. **Bulk reruns:** `baselines.R --dataset m3d|ecomac|k562ess (--gene_subset truths for K562) --methods genie3,grnboost2,pidc`, at the production seed 20261016. GENIE3 and PIDC are reused from their `.rds` caches. This regenerates `metrics.csv`, `jack_replicates.csv` and `per_regulon*.csv` for all methods from the same splits.
   - Seed variability: M3D already has legacy_es seeds 1-3 from job 0894. Add a second seed on EcoMAC and K562 only if compute allows. The ratio SD on M3D is 0.002, so a production-seed table is adequate, and GENIE3 also uses a single seed.
5. **BEELINE:** `beeline_baselines.R` recomputes GENIE3 every time; it has no cache. Add a GENIE3 score cache, or score only GRNBoost2 and then rerun the summary block. Either way, check that the stored GENIE3, PIDC, CLR and |r| rows are reproduced exactly. All 14 cases at seed 20261021, plus seeds 1-3 for variability (cheap).
6. **Regenerate the tables:** `macro_aupr.R` (whole E. coli and K562 families), then `sota_table.R`, then Supplementary Tables 9 and 10, and build Fig. 7b with the new GRNBoost2. A fact-check pass (worker) then updates the text that changes (section 4). The PI edits the Collins checkpoint page, which is shared.

**Results, tables and figures that change:**
- `results/baselines_{m3d,ecomac,k562ess_truthgenes_genie3}/` (metrics, jackknife, per-regulon);
- `results/macro_aupr_v2` (or v3);
- `results/beeline_baselines/{metrics,summary,comembership}.csv`;
- Table 2 (paper) / `sota_table.md`: the GRNBoost2 column of all 9 rows in both panels, plus possibly the bold pattern of the K562 rows;
- Supplementary Table 9 (GRNBoost2 rows; Holm-adjusted P of the whole family);
- Supplementary Table 10 (GRNBoost2 row);
- Fig. 7b (not yet built);
- text in `paper.md` (Results "Against established methods", the BEELINE paragraph, Methods "Competitor methods", the Table 2 footnote), `supplement.md` (Notes 10 and 13), `abstract_points.md` (lines 53 and 62), `docs/checkpoint_collins.md` and its page;
- `preregistration_scfm.md` (descriptive baseline source; section 6).
- **No change:** the Shapley, MI vs ρ and lexA figures, which do not use GRNBoost2.

**Compute** (CPU-hours from job 0894's measured times; wall time at 24 threads):

| Item | legacy_es | (no_es, not recommended) |
|---|---|---|
| M3D (4,297 genes × 466) | 24 CPU-h, 1 h per seed | 230 CPU-h, 9.4 h per seed |
| EcoMAC (4,189 × 2,198) | about 150-250 CPU-h, 7-10 h per seed (current took 23 min; ×18 trees) | about 1,000+ CPU-h |
| K562 (5,671 × 707) | about 120-200 CPU-h, 5-8 h per seed (current took 17 min) | about 600+ CPU-h |
| BEELINE, 14 cases | about 15 CPU-h, under 1 h per seed (hESC_1000 took 3.7 min at 24 workers) | about 250 CPU-h per seed |
| Rescoring (`baselines.R` macro, `macro_aupr.R`, BEELINE summary) | about 1 h per bulk dataset at 8 workers, plus about 1 h | |

Total for legacy_es (production seed everywhere, plus BEELINE seeds 1-3): about 350-500 CPU-hours, or roughly 15-25 h of wall time at 24 threads. That fits next to the current queue within the 70-80% utilization rule. Memory per run is small: job 0894 declared 12 GB for 24 workers.

## 6. Preregistration and disclosure

- **No formal preregistration is required.** GRNBoost2 is a descriptive comparator in Table 2 and Supplementary Tables 9 and 10. No preregistered confirmatory test uses it.
- **But the rerun rule must be fixed before the new runs, because the variant results have been seen.** Variant results have been seen on M3D (TF-node metric, 26 Oct; Table 2 metric, 27 Oct) and on hESC_1000 and mHSC-E_500. Commit a one-paragraph rule (this memo, section 5: legacy_es, production seeds, scoring scripts unchanged, all GRNBoost2 outputs replaced, no_es as sensitivity only) before queueing the EcoMAC, K562 and remaining BEELINE runs. Whatever the new numbers are, the variant is not revisited.
- **`docs/preregistration_scfm.md`** (Analysis 4, "GENIE3, GRNBoost2 and PIDC from results/beeline_baselines, descriptive"): add a change-log entry before any scfm scoring. It should say that the GRNBoost2 baseline is the legacy-early-stopping rerun, with the reason and the date.
- **Macro Holm family:** the family membership and size (40) stay as fixed. Note in Supplementary Note 10 that the GRNBoost2 inputs were replaced after the first run, and that adjusted P values of other comparators can move as a result.
- **Disclosure text** (Methods plus a supplementary note), stating the order:
  1. GRNBoost2 results were first generated and reported (16-21 October) with arboreto 0.1.6 under scikit-learn 1.9.1.
  2. On 22 October we found that arboreto's early stopping is broken under scikit-learn ≥ 1.3 (`oob_improvement_` semantics), which stops most targets after 25-36 trees.
  3. A diagnostic designed after that discovery (3 cases; current, arboreto dask, legacy monitor, no early stopping; 3 seeds) showed that the fault, and not our wrapper or the seed, caused GRNBoost2's last place.
  4. All GRNBoost2 results were then rerun with the pre-1.3 stopping signal.
  5. The superseded numbers and the diagnostic are reported in the supplement.
- **Also say that this affects other users.** Any arboreto or pySCENIC GRNBoost2 run under scikit-learn ≥ 1.3 has the same fault. Consider reporting it upstream; check whether an arboreto issue already exists before claiming novelty.
