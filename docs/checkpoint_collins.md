# CLR as attention: checkpoint for Jim Collins

*Boris, 21 October 2026 (yeast test, item 8, added 11 November 2026; exploratory follow-ups A and B added 12 November 2026). Status of the project, what we have in hand, what did not work, and the decisions ahead. The key results in items 1 to 7 have been checked by an independent statistical review; numbers are from the reviewed outputs. Item 8 (DREAM5, run 11 Nov) and the exploratory follow-ups have had a separate statistical review (12 Nov, `docs/plans/stats_review_2026-11-12.md`), and its corrections are applied in this version (v15). **GRNBoost2 rows regenerated 29–31 Oct 2026 (legacy early stopping, after an early-stopping fault was found on 22 Oct); GRNBoost2 claims reviewed 31 Oct. The shared claude.ai page was republished on 12 Nov 2026 (version 14; it also carries the draft abstract, the graphical abstract and the method pseudocode, which are not in this file).***

**Update, 12 November.** New: exploratory follow-ups. (A) The log-scale depth re-test, preregistered before any compute: the pre-specified reading is (b) "not resolution-like", i.e. the fine-to-coarse ordering was not observed; with 3 and 5 truths and flat depth curves the test could neither establish nor exclude it. (B) The PI's idea of robust row-standardisation (median/IQR instead of mean/SD) in the CLR step: exploratory, no variant is better across the board, nothing goes to the paper. Also updated: the attention-width test (finished 27 October, headline unchanged) and the status lines for depth and robust standardisation, which were out of date.

**Update, 12 November (statistical review).** A separate statistical review of item 8 and the follow-ups found five points that needed correcting. All five are applied. (1) The Y-M verdict "Equivalent within 10%" sits on the boundary and is not a finding; the series-level interval is now given. (2) The yeast official-score sentence that singled out attention's AUROC significance is replaced by the full ranking, which supports no method ordering on yeast. (3) The statement that best depth scales with module size is no longer made: across the 7 truths the slope is not distinguishable from 0. (4) "Not resolution-like" is glossed with the limits of a test on 3 and 5 truths. (5) The summary of hub drift in follow-up B now matches its own table. Smaller points: a ledger of every preregistered test, including the CLR(log P^t) test of 26 October that failed to confirm, a note on the number of looks, and caveats on the intervals and the hub controls.

## In one paragraph

CLR (2007) standardizes each gene's mutual information against that gene's own background (the 2007 code estimates MI by B-spline smoothing, as we do, and combines the two directions by a Euclidean formula; the softmax and the walk are new). That is a context-normalized attention score, published a decade before transformers. It remains a standard baseline for de novo, unsupervised network inference. We complete it as an attention mechanism, with no learned parameters: a softmax over each gene's top CLR partners, then a short multi-step walk (depth t, chosen without labels on held-out experiments). On *E. coli* compendia this recovers strong-evidence transcriptional regulons better than CLR and better than GENIE3, PIDC and GRNBoost2 in pooled ranking. The gain requires context normalization, it survives controls for batch structure and, on strong-evidence regulons, for hub genes, and it concentrates in large regulatory programs. Per regulon, attention ties CLR and GENIE3. It loses on tight co-expression cliques such as ribosomal complexes. In the one independent, preregistered test (DREAM5 yeast, item 8), the co-membership gain was 5.8% (interval +1.8% to +10.0% with group-level resampling, including zero with series-level resampling), below the pre-specified 10% margin, and every method was near chance. A mechanism consistent with the *E. coli* results, though not isolated by them: the walk trades pair-level resolution for module-level recovery. That trade is not seen in the yeast or in silico networks (item 8).

## What we have in hand

**1. The headline (E. coli, M3D, 466 experiments; RegulonDB regulon co-membership)**
- Attention beats CLR by +0.022 AUPR [+0.019, +0.026]. Depth is re-selected without labels in each of 1,000 nested bootstrap draws, and attention is ahead in all 1,000.
- The result is consistent across 50 laboratory splits (+0.022 [+0.017, +0.027]) and across a replicate-aware bootstrap on the 907 chips (+0.020).
- The held-out depth is stable: median t = 9, with the criterion flat within 1% for t = 7–11.
- Scope: items 1 to 3 come from the development organism and lineage (M3D with RegulonDB, and the overlapping EcoMAC), after many looks (see the ledger below). The module-for-edge trade in item 4 is an *E. coli* pattern; it is not seen in the DREAM5 yeast or in silico networks (item 8).

**2. Why it works**
- Context normalization is essential. A softmax on uncalibrated |r| or MI gets about a third of the gain and none of the coherence gain.
- A Shapley decomposition over six components attributes the AUPR gain mainly to depth (+0.018) and to CLR normalization (+0.013).

**3. Controls**
- Removing 0–20 global principal components leaves the gain intact. The top components are biology, not noise: growth rate, motility (σ28) and the σS stationary-phase program.
- Hub control: the walk converges toward well-connected genes, so we score every method within strata of hubness. Attention stays first on strong-evidence truths (ahead of CLR in every hub decile on M3D); on all-evidence truths it ties GENIE3. Under a stratifier that does not use the walk (row mean of raw MI; follow-up B, exploratory) attention stays ahead of CLR on M3D strong evidence (+0.18 enrichment [+0.12, +0.24]). That covers attention against CLR only, not GENIE3 or PIDC.

**4. Against established methods**

Pooled AUPR over the base rate; bold is best in the full row (|r| is best on CORUM, 7.00) or statistically tied (Bonferroni; † tied only after correction).

| Data | Truth | Attention | CLR | GENIE3 | PIDC | GRNBoost2 |
|---|---|---|---|---|---|---|
| M3D | RegulonDB strong evidence | **1.68** | 1.44 | 1.51 | 1.45 | 1.36 |
| EcoMAC (2,198 arrays, mixed platforms) | strong evidence | **1.68** | 1.39 | 1.50 | 1.43 | 1.36 |
| M3D, hub-controlled | strong evidence | **1.56** | 1.40 | 1.47 | 1.40 | 1.37 |
| Human K562 Perturb-seq | DoRothEA A–B | **1.09** | **1.06**† | **1.08** | **1.07** | 1.05 |
| Human K562 Perturb-seq | CORUM complexes | 3.51 | **6.39** | 6.32 | **6.33** | 3.66 |

Full table, with all truths and footnotes: `docs/whitepapers/sota_table.md`.

EcoMAC is not an independent replication of M3D. It is a larger, partly overlapping compendium: 828 of its 2,198 arrays (37.7%) have a one-to-one M3D counterpart, covering 91% of M3D's chips (overlap audit, 1 Nov 2026). On the 1,370 EcoMAC arrays without an M3D counterpart, attention minus CLR is +0.0285 on strong-evidence regulons, the same as on all arrays (+0.0278); this is a post hoc sensitivity analysis that is not preregistered and reuses the same truth, so it is a robustness check and not a held-out test.

Intervals in this table resample genes only. In the DREAM5 yeast test the experiment-resampling component of the standard error was of the same size as the gene component or larger (77% of it at group level, 178% at series level, on Y-M), so these intervals understate the total uncertainty and some differences marked as wins may be ties.

- **Per regulon** (each regulon weighted equally): attention ties CLR and GENIE3. Its mean AP is higher than |r|'s in 3 of 4 comparisons after correction and than PIDC's and GRNBoost2's in 1 of 4 each; medians are about equal. Its win rate is highest in regulons of 25 or more genes.
- **Direct TF→target edges:** on *E. coli* (M3D; DREAM5 network 3, same lineage) attention trails CLR. The trade is not seen in DREAM5 yeast (both readouts point up) or in the in silico network (edge ratio 1.15, module ratio 0.98; descriptive).

**5. Human data**
- K562 Perturb-seq pseudobulk: every method is within 10% of random on TF truths. Attention's small DoRothEA lead over CLR survives the attention-derived hub control (it ties GENIE3 and PIDC); under the MI row-mean stratifier of follow-up B (exploratory) the lead is +0.016 enrichment [−0.010, +0.042]. It loses to CLR, GENIE3 and PIDC on protein complexes, and on knockdown-response modules per module and after hub control; 86% of the complex loss comes from ten ribosomal and pre-rRNA complexes. GRNBoost2, which does not aggregate, does about as poorly as attention on these complexes, so the loss is not specific to the walk.
- Single cells (BEELINE: human ESC and hepatocytes, plus five mouse datasets): CLR on normal scores (CLR-NS) has the highest per-dataset mean of the methods tested, including GENIE3, PIDC and GRNBoost2, but it was chosen on these truths, it is not significantly ahead of GENIE3 or PIDC, and GENIE3 has the higher median. Attention adds nothing there.
- **Against learned attention (preregistered, run once, 2 Nov):** on all seven BEELINE datasets, CLR-NS, a formula with no training, is non-inferior and superior to zero-shot attention from scPRINT, a single-cell foundation model trained for this task (mean log2 AUPRC-ratio difference +1.03, exact interval +0.34 to +1.90, positive in 7 of 7; holds when scPRINT gets its best readout). Read it with three qualifications: scPRINT's attention scores the same as the same model run on permuted expression, and below a data-free gene-embedding similarity, so the comparator is weak here; the advantage is absent on cell-type-specific ChIP-seq (scPRINT higher in 11 of 14 cases, exploratory); and scGPT (two datasets, descriptive) is behind CLR-NS on attention but ahead on its token embeddings. On K562 control cells the direction is the same but smaller (+0.27, 0.09 to 0.46, descriptive). Fine-tuning each model on its own cells (descriptive) did not help: scPRINT's R0 fell slightly (−0.05 log2) and CLR-NS stayed ahead of it (+1.11, 0.41 to 1.96). Details: `docs/scfm_scoring.md`.

**6. Perturbation targets (likely a second paper)**
- Attention-based influence ranks a perturbed TF's targets above differential expression in 12 of 15 perturbations (p = 0.015; 9 of 15, p = 0.095, under a laboratory-matched null). Suggestive, not established.
- It is complementary to SSEM-Lasso, re-implemented from Cosgrove et al. 2008: SSEM pinpoints the regulator itself, attention its regulon.
- Differential CLR is the natural method here and is next to port.

**7. Batch effects** (supplement)
- Multi-platform EcoMAC shows platform "Christmas trees" in gene-pair plots.
- Platform-level ComBat helps more than switching from correlation to MI.
- One clean platform (735 arrays) beats all 2,198 pooled.

**8. Yeast test: DREAM5 network 4 (fixed in advance, run once, 11 Nov)**
- DREAM5 yeast compendium (536 chips, 5,950 genes; the challenge's gold standard of 3,940 interactions). Design frozen 4 Nov, 10% relative-AUPR margin (0.1375 on d = log2 of the AUPR ratio), two co-primary readouts reported together, one look at truth. Our implementation of the official scorer reproduced the published scores of all 36 challenge methods (to 1e-3) before ours were scored.
- **Y-M (regulon co-membership):** attention 0.111 vs CLR 0.105, ratio 1.058, d = 0.082 [0.026, 0.137]. Preregistered verdict, as written by the script: Equivalent within 10%. The upper limit is 0.0001 below the margin (0.1374 against 0.1375), which is within the Monte Carlo error of the interval itself (about ±0.002), so the boundary between Equivalent and Inconclusive is not resolved. With the coarser series-level resampling the interval is −0.009 to 0.172 (Inconclusive, includes 0), and with a Bonferroni correction across the two co-primaries it is 0.018 to 0.145 (also Inconclusive). What holds under every version: the estimate is below the 10% margin, so the result is not Superior. Fixed wording: the co-membership headline did not replicate in yeast at the stated margin. The estimate is positive; the interval excludes zero only under group-level resampling.
- **Y-E (official directed-edge score):** attention 0.0216 vs CLR 0.0208, ratio 1.042, d = 0.060 [-0.093, 0.213]. Inconclusive. The design predicted no gain here (d ≤ 0, the M3D pattern); the estimate is a 4% gain (d = +0.060), but the interval runs from a 6% loss to a 16% gain, so the prediction is neither confirmed nor refuted.
- Every method is near chance on yeast edges (AUPR about 0.02 against a base rate of 0.017; AUROC 0.51 to 0.52). Series-level resampling gives a larger spread than the group-level SE_B (Y-M SD 0.040 vs 0.017); with it the Y-M interval for d is −0.009 to 0.172 and the Y-E interval is −0.107 to 0.227. The Y-M interval excluding zero rests on the group-level resampling.
- Networks 3 and 1, descriptive (jackknife-only intervals, understated): E. coli net3 Y-M ratio 1.19, Y-E 0.89 (the M3D pattern; same lineage as the headline, not independent); in silico net1 Y-M 0.98, Y-E 1.15. Official challenge score (descriptive; our reimplementations; mean over the three networks of the two -log10 p scores): GRNBoost2 and GENIE3 lead overall (55.6 and 43.8), attention is third (38.2) and CLR fourth (22.9). On yeast all seven methods lie within AUROC 0.510 to 0.523, so no ranking of methods is supported there.
- PIDC has a 2-entry tied block at the 100,000 cap on nets 1 and 4 (design expected none); dropped by the fixed rule, affects only PIDC.
- Review status: this item had a separate statistical review on 12 Nov, and its corrections are applied above. The design's exploratory list (replicate-averaged input, tree-method module read-outs, hub-stratified scores, other depths, per-TF AUROC) has not been run; anything run now is post-unblinding exploratory. Full report: `docs/dream5_selections/dream5_report.md`.

## Exploratory follow-ups

Two analyses run on 11 and 12 November. Neither enters the paper as it stands. Letters A and B are used so they are not confused with items 1 to 8.

### A. Log-scale depth re-test (preregistered before compute, 8861ae4, 11 Nov; not a blind test)

- **Question.** Is walk depth a resolution dial, with the best depth rising from fine truths (operons) to coarse ones (regulons, sigmulons)? The first test used a linear grid (t = 1 to 32) and two best depths sat at its edge. The re-test uses a geometric grid, t = 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, plus CLR as t = 0.
- **Status of the evidence.** The analysis and the three possible readings were committed before any new compute. It is still not blind or confirmatory: the truths and operators were already used in the linear-grid run and its curves had been seen. Power is small (3 truths in E. coli M3D, 5 in K562); the statistical review adds that reading (a) was nearly unattainable. The results document gives no interpretation beyond the pre-specified readings, and neither does this file.
- **Pre-specified reading: (b) "not resolution-like"** (Spearman rho between truth coarseness and log2 best depth at most 0.5 in at least one organism). E. coli rho = 0.500 (n = 3); K562 rho = 0.700 (n = 5). Reading (a) (rho at least 0.9 in both organisms, no edge censoring) was not met. 0 of 8 truths were censored at the grid edge t = 256; the within-1% interval of K562 CollecTRI reaches 256. The best depth moved from the edge of the linear grid to 96 (E. coli TF regulons) and 48 (K562 CollecTRI).
- **What "not resolution-like" can mean here.** It is the literal outcome of the pre-specified rule: the pre-specified ordering was not observed. With 3 and 5 truths and flat depth curves, the test could neither establish nor exclude it. In E. coli the rule fires unless all three truths are perfectly ordered (rho of at least 0.9 occurred in 10% of bootstrap draws, rho of at most 0.5 in 88%); the order is broken by the sigmulons (5 sets; bootstrap interval of best t 0 to 256). In K562, rho is 0.70 with CLR ranked lowest, and 0.20 over attention depths only (the preregistration records that the primary convention was chosen after the linear-grid curves were seen). The coarseness rank is nominal: the K562 KD-response modules (median 100 genes) are ranked finer than DoRothEA (median 13).
- **Equivalence gate, passed:** at t = 0, 1, 2, 4, 8, 16, 32 the new run reproduces the saved linear-grid curves (AUPR, base rate, coherence) to within 5.0e-16 (tolerance 1e-9; 21 E. coli rows, 35 K562 rows).

| Organism | Truth, fine to coarse | Median size (sets) | Best t | Within-1% interval of t | AUPR gain over CLR | Bootstrap 95% interval of best t |
|---|---|---|---|---|---|---|
| E. coli | 1 operons | 2 (786) | 4 | 4 to 4 | +0.208 | 4 to 4 |
| E. coli | 2 TF regulons | 11 (158) | 96 | 48 to 128 | +0.027 | 16 to 128 |
| E. coli | 3 sigmulons | 128 (5) | 12 | 6 to 24 | +0.043 | 0 to 256 |
| K562 | 1 CORUM complexes | 4 (1,186) | CLR (0) | 16 to 16* | 0.000 | 0 to 0 |
| K562 | 2 KD-response modules | 100 (585) | 24 | 24 to 32 | +0.017 | 24 to 32 |
| K562 | 3 DoRothEA A–B | 13 (105) | 4 | 3 to 12 | +0.008 | 3 to 24.2 |
| K562 | 4 DoRothEA A–C | 28 (269) | 8 | 3 to 16 | +0.005 | 3 to 32 |
| K562 | 5 CollecTRI | 16 (490) | 48 | 24 to 256 | +0.014 | 32 to 256 |

Best t is the highest-AUPR point over CLR and the grid. Within-1% interval: the run of grid depths around the best attention depth whose AUPR is at least 99% of the maximum (t of at least 1; widths are limited by the grid spacing of about half an octave). AUPR gain: absolute AUPR at the best depth minus CLR's. Bootstrap: 1,000 resamples of the modules in each truth, operators held fixed; a lower end of 0 means CLR was best in some draws. For comparison, the linear-grid run gave best depths of 4, 32 and 13 (E. coli) and CLR, 26, 4, 7 and 32 (K562). *K562 CORUM: CLR beats all 16 attention depths; the best attention depth is t = 16, with an AUPR 0.037 below CLR's (0.0894), so CLR is best at every depth.

- **Power.** Exact one-sided permutation P for rho: E. coli 0.500, K562 0.117. The smallest attainable one-sided P is 0.167 (n = 3) and 0.008 (n = 5). Descriptive; not part of the reading rule.
- **Best depth against module size** (secondary, truths with best t of at least 1; bootstrap over modules): pooled slope of log2 best t on log2 median module size 0.197 [0.078, 0.668] over 7 truths. E. coli 0.186 [0.054, 1.929] (3 truths; 952 of 1,000 draws usable), K562 0.380 [−0.203, 0.581] (4 truths). The interval resamples modules within fixed truths and operators; it does not reflect variation between truths. Least squares across the 7 truths gives 0.197 with a 95% interval of −0.77 to 1.16 (P = 0.62), and dropping the operon truth changes the slope to −0.22. The slope is descriptive.
- Figures: `results/granularity_log/depth_log_summary.png` (best depth by coarseness, best depth against module size, interval width) and `results/granularity_log/depth_log_curves.png` (AUPR gain over CLR at each depth, one panel per truth). Sources: `docs/depth_log_results.md`, `docs/preregistration_depth_log.md`.

### B. Robust row-standardisation of the CLR step (the PI's idea; EXPLORATORY, not for the paper)

**Exploratory.** Not preregistered. No variant was selected, no p-values are computed (none should be read into the intervals), K562 CORUM is not held out (it was used in the attention-width test and in earlier CORUM work), M3D RegulonDB is the project's tuning set for depth, and nothing here goes to the paper. Every variant, truth and read-out was fixed before the runs. Intervals are paired delete-half gene jackknife intervals (100 replicates) with the networks held fixed; they do not resample experiments, cells or seeds. The two `tables.md` files hold about 700 such intervals (711 counted); this page shows 56 of the paired contrasts. With that many, "the interval excludes 0" is not evidence of a real difference, and it is used here only to describe.

**Question.** CLR z-scores each gene's MI row with the row mean and SD. A gene in a module with many strong partners inflates its own row SD, which lowers the z-scores of its true partners; this could hurt large complexes such as the ribosome. Does z = (x − median) / (IQR/1.349) change the pairwise score or the attention score? Everything after the row z (truncation at 0, symmetrisation, softmax, lazy walk) is the default path; the package itself is unchanged.

| Variant | Row centre and scale |
|---|---|
| V0 | mean and SD of the MI row: default CLR (the zeroed diagonal is part of the row, as in the package) |
| V1 | median and IQR/1.349 of the raw MI row (the PI's proposal); diagonal excluded |
| V2 | median and IQR/1.349 of MI^(1/3) (cube-root scale); diagonal excluded |
| V3 | the existing rCLR: two-way median polish of MI^(1/3), per-gene MAD; as defined it is not truncated at 0 |
| V4 | median and 1.4826 × MAD of the raw MI row; diagnostic; diagonal excluded |

**Checks.** V0 as re-implemented equals the package (largest difference in CLR scores 7.1e-14 on M3D and 3.4e-14 on K562) and its pooled AUPR equals the project's baselines run to within 4e-16. No row had a zero scale under any variant.

**Pooled enrichment on the method-table truths** (pooled AUPR divided by the base rate, 1.00 = random; format of `docs/whitepapers/sota_table.md`):

| Truth | Base rate | CLR V0 | V1 | V2 | V3 | V4 | Att V0 | V1 | V2 | V3 | V4 | GENIE3 | PIDC | GRNBoost2 | \|r\| |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| M3D, RegulonDB strong evidence | 0.094 | 1.44 | 1.45 | 1.46 | 1.52 | 1.45 | 1.68 | 1.67 | 1.66 | 1.66 | 1.66 | 1.51 | 1.45 | 1.36 | 1.23 |
| M3D, all evidence | 0.120 | 1.37 | 1.42 | 1.41 | 1.47 | 1.43 | 1.44 | 1.44 | 1.42 | 1.43 | 1.45 | 1.40 | 1.40 | 1.26 | 1.23 |
| K562, CORUM all | 0.013 | 6.39 | 6.94 | 6.37 | 6.07 | 6.44 | 3.51 | 3.26 | 3.44 | 3.78 | 2.99 | 6.32 | 6.33 | 3.66 | 7.00 |
| K562, CORUM without ribosome | 0.012 | 2.18 | 2.28 | 2.34 | 2.57 | 2.17 | 2.05 | 1.96 | 2.09 | 2.27 | 1.86 | - | - | - | - |
| K562, CORUM large (20+ genes) | 0.126 | 2.29 | 2.31 | 2.30 | 2.39 | 2.27 | 1.93 | 1.80 | 1.89 | 2.08 | 1.75 | - | - | - | - |
| K562, KD response | 0.079 | 2.42 | 2.61 | 2.53 | 2.55 | 2.60 | 2.38 | 2.42 | 2.19 | 2.07 | 2.45 | 2.45 | 2.35 | 1.84 | 2.45 |
| K562, DoRothEA A–B | 0.183 | 1.06 | 1.06 | 1.06 | 1.08 | 1.05 | 1.09 | 1.08 | 1.08 | 1.08 | 1.07 | 1.08 | 1.07 | 1.05 | 1.05 |

Columns: pairwise score (CLR step) V0 to V4; attention at t = 9 V0 to V4; reference methods. M3D: 466 experiments, 4,297 genes. K562: essential Perturb-seq pseudobulk, network from knockdown half A, 5,671 truth genes; the KD-response truth comes from half B. "Without ribosome" drops the 8 complexes dominated by ribosomal proteins; "large" is complexes of 20 or more genes. The ribosome-only truth is saturated (base rate 0.75; pairwise enrichment 1.21 for every variant), so it is not shown. The V0 columns equal the method-table rows (Attention at t = 9, CLR). Reference columns are taken from the same runs and the same scoring (`results/baselines_m3d`, `results/baselines_k562ess_truthgenes_genie3`) and are directly comparable; "-" means the method was scored on the full truth only, not on this subset. Unlike the method table, nothing is bolded: the variants were not part of its Bonferroni-corrected tie family, so no best-or-tied marking is made. GRNBoost2 uses legacy early stopping (rerun 29–31 October). As in the method table, depth t = 9 was chosen on M3D and intervals resample genes only. The hub-only column of the full table is not repeated.

**Paired contrasts against V0** (difference in pooled enrichment, variant minus V0, with the 95% delete-half gene jackknife interval; descriptive, no p-values):

| Read-out | Data | Truth | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|---|
| Pairwise | M3D | strong evidence | +0.013 [−0.004, +0.029] | +0.022 [+0.003, +0.040] | +0.077 [+0.042, +0.112] | +0.009 [−0.008, +0.026] |
| Pairwise | M3D | all evidence | +0.047 [+0.035, +0.059] | +0.043 [+0.031, +0.055] | +0.102 [+0.072, +0.131] | +0.058 [+0.043, +0.073] |
| Pairwise | K562 | CORUM all | +0.547 [+0.345, +0.749] | −0.018 [−0.267, +0.232] | −0.318 [−0.784, +0.147] | +0.052 [−0.224, +0.327] |
| Pairwise | K562 | CORUM without ribosome | +0.100 [+0.019, +0.182] | +0.158 [+0.050, +0.266] | +0.381 [+0.182, +0.581] | −0.019 [−0.155, +0.117] |
| Pairwise | K562 | CORUM large | +0.023 [−0.000, +0.045] | +0.010 [−0.010, +0.030] | +0.097 [+0.046, +0.149] | −0.021 [−0.056, +0.013] |
| Pairwise | K562 | KD response | +0.189 [+0.158, +0.220] | +0.102 [+0.075, +0.128] | +0.127 [+0.057, +0.198] | +0.178 [+0.135, +0.221] |
| Pairwise | K562 | DoRothEA A–B | −0.002 [−0.007, +0.002] | +0.002 [−0.003, +0.006] | +0.025 [+0.013, +0.037] | −0.008 [−0.013, −0.003] |
| Attention t = 9 | M3D | strong evidence | −0.009 [−0.037, +0.020] | −0.014 [−0.047, +0.020] | −0.020 [−0.070, +0.030] | −0.011 [−0.038, +0.016] |
| Attention t = 9 | M3D | all evidence | +0.000 [−0.018, +0.019] | −0.023 [−0.046, +0.001] | −0.010 [−0.037, +0.018] | +0.011 [−0.009, +0.031] |
| Attention t = 9 | K562 | CORUM all | −0.250 [−0.539, +0.039] | −0.067 [−0.447, +0.313] | +0.272 [−0.351, +0.895] | −0.519 [−0.980, −0.059] |
| Attention t = 9 | K562 | CORUM without ribosome | −0.088 [−0.233, +0.057] | +0.040 [−0.112, +0.193] | +0.216 [+0.031, +0.401] | −0.196 [−0.390, −0.001] |
| Attention t = 9 | K562 | CORUM large | −0.127 [−0.171, −0.082] | −0.039 [−0.093, +0.015] | +0.158 [+0.051, +0.265] | −0.180 [−0.252, −0.108] |
| Attention t = 9 | K562 | KD response | +0.040 [+0.002, +0.078] | −0.181 [−0.231, −0.130] | −0.302 [−0.386, −0.218] | +0.071 [+0.019, +0.123] |
| Attention t = 9 | K562 | DoRothEA A–B | −0.017 [−0.027, −0.007] | −0.011 [−0.023, +0.001] | −0.013 [−0.040, +0.014] | −0.025 [−0.039, −0.012] |

Each variant's attention at t = 9 is compared with V0's attention at t = 9. At each variant's own label-free depth (6 to 7 for V1, V2 and V4; 5 to 6 for V3; V0 is 9 in both datasets) several cells are lower still, for example K562 CORUM all, V1: −0.398 [−0.704, −0.093]. The pooled CORUM-all numbers are dominated by the ribosomal complexes, whose own enrichment is saturated, so the without-ribosome and large rows are the split. In the per-complex (macro) tables, M3D mean per-regulon AP is lower than V0's for every robust variant in both truths, pairwise and attention (−0.005 to −0.014).

- **Attention against CLR, within each variant.** On K562 CORUM all and large complexes, attention is below V0 CLR for every variant, with intervals entirely below 0 (CORUM all, V0: −2.88 [−4.32, −1.44]; across variants −2.6 to −3.4 enrichment units). On M3D attention is above V0 CLR on both truths for every variant (strong evidence +0.21 to +0.23; all evidence +0.05 to +0.08).

**Hub-stratified contrasts: the two stratifiers disagree.** Hub control uses within-stratum percentiles in 10 strata of H_i + H_j. The primary stratifier is the stationary distribution pi of the V0 attention walk (attention-derived, as in the method table); the second is the row mean of the raw MI, which does not depend on any variant. It is not neutral with respect to the treatment, though: the variants differ precisely in how they remove the row location (mean, median, or two-way median polish), so neither stratifier is the hub-controlled answer. On K562 CORUM all, the V1 pairwise gain vanishes under pi (+0.014) and persists under the MI row mean (+0.483). The pi stratifier is derived from V0's own attention walk, and the two stratifiers disagree on this contrast.

| Variant | K562 CORUM all, pairwise: variant minus V0 | Unstratified | Stratified by pi (V0 attention) | Stratified by MI row mean |
|---|---|---|---|---|
| V1 | median / IQR | +0.547 [+0.345, +0.749] | +0.014 [−0.178, +0.206] | +0.483 [+0.234, +0.733] |
| V2 | cube-root median / IQR | −0.018 [−0.267, +0.232] | −0.210 [−0.478, +0.058] | +0.139 [−0.115, +0.394] |
| V3 | rCLR | −0.318 [−0.784, +0.147] | −0.220 [−0.624, +0.184] | +0.430 [+0.062, +0.799] |
| V4 | median / MAD | +0.052 [−0.224, +0.327] | −0.256 [−0.505, −0.006] | +0.195 [−0.033, +0.423] |

Under the MI row mean on M3D the pairwise gains stay close to the unstratified ones (all evidence: +0.044 to +0.073).

**Hub drift.** Share of the 100,000 top-scoring pairs that touch at least one of the 5% most connected genes (chance about 0.098), with connectedness taken from pi (V0 attention) and from the MI row mean:

| Variant | M3D pairwise, pi | M3D pairwise, MI | M3D attention, pi | M3D attention, MI | K562 pairwise, pi | K562 pairwise, MI | K562 attention, pi | K562 attention, MI |
|---|---|---|---|---|---|---|---|---|
| V0 (CLR) | 0.211 | 0.169 | 0.257 | 0.116 | 0.336 | 0.419 | 0.608 | 0.413 |
| V1 median / IQR | 0.277 | 0.175 | 0.348 | 0.129 | 0.373 | 0.507 | 0.711 | 0.433 |
| V2 cube-root | 0.279 | 0.130 | 0.328 | 0.112 | 0.368 | 0.419 | 0.636 | 0.331 |
| V3 rCLR | 0.270 | 0.034 | 0.294 | 0.077 | 0.328 | 0.259 | 0.543 | 0.171 |
| V4 median / MAD | 0.276 | 0.205 | 0.348 | 0.133 | 0.347 | 0.564 | 0.685 | 0.455 |

By pi, V1, V2, V3 and V4 touch the most connected genes more often than V0 in M3D; in K562 V1, V2 and V4 do and V3 does not. By the MI row mean, V3 is below V0 in both datasets and both read-outs, while V1 and V4 are above V0 in K562 (pairwise 0.507 and 0.564 against 0.419) and close to V0 in M3D (pairwise 0.175 and 0.205 against 0.169). The share of non-zero pairwise scores is higher for V1, V2 and V4 (0.62 to 0.64) than for V0 (0.53 to 0.54).

Caveats. The shares have no intervals, so small differences (for example 0.169 against 0.175) cannot be interpreted. pi is the stationary distribution of V0's own walk, so for V0's attention the pi share is partly built in, and for the other variants it measures overlap with V0's hubs. The attention columns are at the fixed t = 9, which mixes hub preference with mixing speed: V1, V2 and V4 have denser operators and shallower label-free depths (5 to 7). At each variant's own label-free depth, every variant's K562 attention share by pi is below V0's (V1 0.571, V2 0.502, V3 0.382, V4 0.559, against 0.608); in M3D V1, V2 and V4 stay above V0 by pi (0.305, 0.270, 0.283 against 0.257) and V3 falls below (0.244).

**What the row scales track** (Spearman correlation across genes; 4,297 in M3D, 5,671 in K562):

| Row scale | Strong partners, M3D | Strong partners, K562 | Own module size, M3D | Own module size, K562 |
|---|---|---|---|---|
| SD | 0.86 | 0.89 | −0.05 | +0.17 |
| IQR / 1.349 | 0.72 | 0.85 | −0.05 | +0.13 |
| MAD | 0.71 | 0.83 | −0.06 | +0.12 |

Strong partners: the number of partners whose MI is in the top 1% of all cells. Own module size: the largest RegulonDB strong-evidence regulon (M3D) or CORUM complex (K562) containing the gene (1,879 and 2,260 genes). The SD and both robust scales rise with the number of strong partners (Spearman 0.71 to 0.89), so the robust scales are not independent of it. They are only weakly related to a gene's own module size. SD relative to the IQR scale rises modestly from the lowest to the highest decile of strong partners (median 1.10 to 1.19 in M3D, 1.07 to 1.35 in K562); in the top K562 decile the SD is about 1.35 times the IQR scale and 1.6 times the MAD scale.

**Summary, as far as the tables go.** No variant is consistently better than V0 across datasets, truths, read-outs or variants. In the pairwise score, V1 is higher on K562 CORUM all (+0.547) and on KD response (+0.189), and every robust variant is at or above V0 on M3D in pooled enrichment but below V0 in mean per-regulon AP (pairwise −0.007 to −0.014), but the direction does not carry over to attention at t = 9, where V1, V2 and V4 are lower than V0 on large complexes and V3 is higher. By pi (from V0's own walk), the pairwise top pairs of V1, V2 and V4 touch V0's hub genes more often than V0's do; for attention this holds at t = 9 but not, in K562, at each variant's own label-free depth (V1 0.571, V2 0.502, V4 0.559, against V0 0.608). The two hub stratifiers disagree on at least one contrast (K562 CORUM all, V1: +0.014 under pi, +0.483 under the MI row mean). V3 (rCLR) lowers the hub share by the MI row mean in both datasets and by pi in K562; V2 lowers it by the MI row mean in M3D and in K562 attention. Both connectedness rankings are tied to what is being compared (pi to V0's walk, the MI row mean to row centring), so neither is a neutral hub measure, and the shares have no intervals. Nothing here is a held-out test, and the order of the variants is not a recommendation. Sources: `docs/robust_z_results.md`, `results/robust_z/{m3d,k562}/tables.md` and the CSVs next to it.

## Every preregistered test so far

Each outcome is taken from the verdict document of that test (`docs/clr_logpt_verdict/`, `docs/keff_verdict/`, `docs/scfm_scoring.md`, `docs/reciprocity_yeast_selections/`, `results/dream5/truth_net4/`, `docs/depth_log_results.md`), in the order they finished.

| Finished | Test (fixed design) | Outcome |
|---|---|---|
| 26 Oct | CLR applied to log P^t as a fix for protein complexes; confirmatory on *E. coli* M3D and EcoMAC | **Not adopted: failed to confirm.** Complex superiority over attention at t = 9 needed one-sided p < 0.025 on both datasets: M3D +0.026 AUPR (p = 0.079), EcoMAC +0.029 (p = 0.051). Regulon non-inferiority passed. Attention stays the method. |
| 27 Oct | Attention width: 8 variants, 32 primary tests, Holm-corrected | **No variant meets the decision rule;** the headline width (10 effective partners) is unchanged. K25 and O3A each won 2 of their 4 tests and O3B 1, none enough to meet the rule. |
| 2 Nov | Yeast reciprocity | **Not run.** The data-quality gate failed in all 10 draws, so by the preregistration no truth was opened. |
| 2 Nov | CLR-NS against zero-shot scPRINT attention, 7 BEELINE datasets | **Non-inferior and superior** (Hodges-Lehmann +1.03 log2 AUPRC ratio, exact interval +0.34 to +1.90, p = 0.0078, positive in 7 of 7). The comparator scores at the level of its own permutation control, so it is a weak comparator. |
| 11 Nov | DREAM5 yeast, network 4 (Y-M and Y-E) | **Not replicated at the stated margin.** Y-M: Equivalent within 10%, on the boundary with Inconclusive, not Superior. Y-E: Inconclusive. |
| 11 Nov | Log-scale depth re-test (not blind) | **Reading (b), "not resolution-like";** with 3 and 5 truths the test could neither establish nor exclude the ordering. |

**Many looks.** These six are the confirmatory tests with a fixed design, and one of them (scPRINT) was positive. The *E. coli* headline of items 1 to 4 is not one of them; it comes from development on M3D. Around the six sit dozens of exploratory analyses on the same M3D/RegulonDB and K562 truths (FDR operators, value transforms, multiscale algorithms, width rounds, rCLR, robust standardisation, macro averaging, hub controls, granularity). The nested bootstrap protects the depth choice. The other settings (10 effective partners, 50 candidates, lazy weight 0.5, the MI estimator, the truth construction with operon exclusion) were settled during development on the same data. Only the yeast test lies outside that lineage.

## What did not work (reported as negatives)
- Hard FDR-thresholded operators: they isolate 62–77% of genes.
- Yeast reciprocity test: not run; its preregistered data-quality gate failed in all 10 draws (median chance-corrected reproducibility 0.12 for mutants against 0.30 for wild types), so no truth was opened. Reciprocity stays a synthetic-only candidate (Supplementary Note 18).
- Value transforms and diffused-profile readouts.
- Wider gene contexts on single cells.
- Data-driven choice of attention width: all three options failed review before running, and were tested empirically anyway under a preregistration (finished 27 October). The headline width (10 effective partners) is unchanged and no variant meets the decision rule. A separate exploratory width study, on synthetic and label-free data with no truth opened, closed negatively on 31 October.
- A label-free multiscale algorithm based on spectral plateaus: there are no plateaus; the structure is a continuum.
- CLR applied to log P^t as a fix for complexes: preregistered, failed to confirm (26 October; complex superiority p = 0.079 on M3D and 0.051 on EcoMAC; regulon non-inferiority passed). Attention stays the method.
- Depth as a discrete "resolution dial" ordering operons → regulons → sigmulons: not supported on a linear depth grid, and the preregistered re-test on a geometric grid to t = 256 did not find the fine-to-coarse ordering either (pre-specified reading (b); with 3 and 5 truths it could not have established or excluded it, follow-up A). Across truths, the best depth is not clearly related to median module size (slope 0.20; module-bootstrap interval 0.08 to 0.67, which ignores variation between truths; across the 7 truths −0.77 to 1.16; the sign depends on the operon truth). Within regulons, an earlier per-set analysis on the linear grid found best depth rising with set size (t ∝ n^0.43) but not for random sets.

## In progress and recently completed
- Learned attention (done 2 Nov; see item 5 above).
- GENIE3, GRNBoost2 and PIDC on the BEELINE single-cell benchmark: done. CLR is comparable to GENIE3 and PIDC. GRNBoost2 (legacy early stopping) is below CLR in 6 of 7 datasets (mean log2 −0.15; Wilcoxon P = 0.031, sign P = 0.125; hESC is the exception), a descriptive count, not a significant difference.
- Preregistered attention-width test, including multi-head attention as the candidate fix for the complex loss: done (27 October; 8 variants, 32 primary tests, Holm-corrected). Headline unchanged at 10 effective partners; no variant meets the decision rule. Multi-head attention did not gain on K562 CORUM (−0.003 AUPR against the headline, Holm P = 0.54). The write-up and the list of individual tests that variants won are still to do (`docs/keff_verdict/`).
- Directed DREAM5 evaluation: done (item 8; ARACNE was not part of the design). The write-up is in progress.
- Log-scale depth analysis: done (follow-up A). Robust row-standardisation, the PI's median/IQR idea: done, exploratory (follow-up B).
- Still planned: perturbation-diversity learning curves, differential CLR, and a hub stratifier for the method table that does not depend on the walk (follow-up B's MI row mean is one, but it is tied to row centring and is not neutral).

## Decisions for you and me
1. **Venue:** Nature Machine Intelligence (150-word abstract) or Nature Communications.
2. **Scope:** one paper on CLR attention, with perturbation targets and differential CLR as a second paper.
3. **The central claim:** "better recovery of regulatory programs from expression alone", not "better network inference". Is that strong enough for the target venue? The statistical review suggests scoping the claim to *E. coli*, or presenting the work as a reframing of CLR with characterised trade-offs, since every positive result comes from the development organism and lineage.
4. **The human single-cell story:** CLR-NS rather than attention is the single-cell result. Does it belong in this paper?
5. **The yeast result:** Y-M has a positive estimate under the 10% margin (the verdict "Equivalent" is on the boundary with Inconclusive), Y-E is inconclusive, and the network 3 pattern (module gain, edge cost) is not independent of the headline. How should this enter the paper and the abstract, which is scoped to E. coli?
6. **AI disclosure:** the analysis code, runs and drafts were produced with Claude (Anthropic), with independent AI statistical reviews, under my direction. A contributions statement is drafted to match Nature Portfolio policy.
