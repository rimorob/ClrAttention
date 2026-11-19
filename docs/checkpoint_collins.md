# CLR as attention: checkpoint for Jim Collins

*Boris, 21 October 2026 (yeast test, item 8, added 11 November 2026; exploratory follow-ups A and B added 12 November 2026; follow-ups C and D added 12 and 13 November 2026; item 9 on tcaCLR and the confirmatory runs added 19 November 2026). Status of the project, what we have in hand, what did not work, and the decisions ahead. The key results in items 1 to 7 have been checked by an independent statistical review; numbers are from the reviewed outputs. Item 8 (DREAM5, run 11 Nov) and the exploratory follow-ups A and B have had a separate statistical review (12 Nov, `docs/plans/stats_review_2026-11-12.md`), and its corrections are applied in this version. Item 9 and follow-ups C and D are the scoring scripts' own outputs and have not had a separate statistical review. **GRNBoost2 rows regenerated 29–31 Oct 2026 (legacy early stopping, after an early-stopping fault was found on 22 Oct); GRNBoost2 claims reviewed 31 Oct. The shared claude.ai page was republished on 12 Nov 2026 (version 14; it also carries the draft abstract, the graphical abstract and the method pseudocode, which are not in this file).***

**Update, 12 November.** New: exploratory follow-ups. (A) The log-scale depth re-test, preregistered before any compute: the pre-specified reading is (b) "not resolution-like", i.e. the fine-to-coarse ordering was not observed; with 3 and 5 truths and flat depth curves the test could neither establish nor exclude it. (B) The PI's idea of robust row-standardisation (median/IQR instead of mean/SD) in the CLR step: exploratory, no variant is better across the board, nothing goes to the paper. Also updated: the attention-width test (finished 27 October, headline unchanged) and the status lines for depth and robust standardisation, which were out of date.

**Update, 12 November (statistical review).** A separate statistical review of item 8 and the follow-ups found five points that needed correcting. All five are applied. (1) The Y-M verdict "Equivalent within 10%" sits on the boundary and is not a finding; the series-level interval is now given. (2) The yeast official-score sentence that singled out attention's AUROC significance is replaced by the full ranking, which supports no method ordering on yeast. (3) The statement that best depth scales with module size is no longer made: across the 7 truths the slope is not distinguishable from 0. (4) "Not resolution-like" is glossed with the limits of a test on 3 and 5 truths. (5) The summary of hub drift in follow-up B now matches its own table. Smaller points: a ledger of every preregistered test, including the CLR(log P^t) test of 26 October that failed to confirm, a note on the number of looks, and caveats on the intervals and the hub controls.

**Update, 12 November (later).** New: follow-up C describes the log-then-CLR read-out (a3r, not adopted on 26 October) and a stacked version of it (a3r2, a3r3), run today. Stacking did not help: neither two nor three layers meets the pre-written rule on K562 or *E. coli*. Panel 5 of the graphical abstract now carries a3r, a3r2 and a3r3 as dashed bars.

**Update, 13 November.** New: follow-up D, an exploratory estimate of a3r (the log-then-CLR read-out) against attention on yeast protein complexes, run once on a dataset and truth a3r had never been scored on, under a plan fixed before compute. The total 95% interval for the pooled AUPR difference is +0.067 to +0.125 (estimate +0.096), about 3.5 times the *E. coli* estimate; the gain is carried by large complexes, and the interval that treats complexes as the sampling unit includes 0. It is an estimate, not a confirmatory test, and a3r stays not adopted. Also new: pages 3 and 3b of the graphical abstract, a scoreboard of every metric (wins and losses), and a3r rows in the page 2 table. BEELINE medians are over 44 case-truth cells (14 dataset-and-gene-set cases), not 14 cases, and the captions are corrected.

**Update, 19 November.** New: section 9, on tcaCLR (CLR applied to the log of the walk), and the confirmatory runs of the plan frozen on 17 November. The names aCLR and tcaCLR are those of the short pre-read; the older sections and the figures of this page still say "attention" and "a3r", and the mapping is given once, at the start of section 9. **Arm A** (HepG2 Perturb-seq pseudobulk, Complex Portal human complexes), run once: tcaCLR exceeds aCLR by +0.0268 AUPR, total 95% interval +0.0111 to +0.0424, one-sided p = 0.0004, so the preregistered test is met. The gain is not broad: three complexes carry 81% of it, the median complex is not improved, the interval that treats complexes as the sampling unit (−0.014 to +0.067) includes 0, and tcaCLR is not shown to be the best of the seven methods run. The point estimate is close to the two *E. coli* estimates that missed the threshold (+0.026, p = 0.079; +0.029, p = 0.051); what differs is the precision. **Arm B** (Kemmeren deletion compendium, CYC2008), estimate only: +0.0566, total interval +0.0392 to +0.0741. Also new: the planted-module simulation (the gain edge is 20 genes), a note on the CYC2008 licence, two ledger rows, and the decisions list (the question of a further evaluation of tcaCLR is answered; the headline choice and the wording of the complex result remain).

## In one paragraph

CLR (2007) standardizes each gene's mutual information against that gene's own background (the 2007 code estimates MI by B-spline smoothing, as we do, and combines the two directions by a Euclidean formula; the softmax and the walk are new). That is a context-normalized attention score, published a decade before transformers. It remains a standard baseline for de novo, unsupervised network inference. We complete it as an attention mechanism, with no learned parameters: a softmax over each gene's top CLR partners, then a short multi-step walk (depth t, chosen without labels on held-out experiments). On *E. coli* compendia this recovers strong-evidence transcriptional regulons better than CLR and better than GENIE3, PIDC and GRNBoost2 in pooled ranking. The gain requires context normalization, it survives controls for batch structure and, on strong-evidence regulons, for hub genes, and it concentrates in large regulatory programs. Per regulon, attention ties CLR and GENIE3. It loses on tight co-expression cliques such as ribosomal complexes. In the one independent, preregistered test (DREAM5 yeast, item 8), the co-membership gain was 5.8% (interval +1.8% to +10.0% with group-level resampling, including zero with series-level resampling), below the pre-specified 10% margin, and every method was near chance. A mechanism consistent with the *E. coli* results, though not isolated by them: the walk trades pair-level resolution for module-level recovery. That trade is not seen in the yeast or in silico networks (item 8). A second CLR applied to the log of the walk (tcaCLR; no fitted parameters, picked among about 20 read-outs on K562) missed the threshold of its preregistered *E. coli* test on protein complexes (+0.026 and +0.029 AUPR over aCLR, p = 0.079 and 0.051) and met the threshold of its preregistered HepG2 test (+0.027, one-sided p = 0.0004), where three complexes carry 81% of the gain and tcaCLR is not shown to be the best of seven methods (item 9).

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
- **Against learned attention (preregistered, run once, 2 Nov):** on all seven BEELINE datasets, CLR-NS, a formula with no training, is non-inferior and superior to zero-shot attention from scPRINT, a single-cell foundation model trained for this task (Hodges-Lehmann estimate of the log2 AUPRC-ratio difference +1.03, exact interval +0.34 to +1.90, positive in 7 of 7; holds when scPRINT gets its best readout). Read it with three qualifications: scPRINT's attention scores the same as the same model run on permuted expression, and below a data-free gene-embedding similarity, so the comparator is weak here; the advantage is absent on cell-type-specific ChIP-seq (scPRINT higher in 11 of 14 cases, exploratory); and scGPT (two datasets, descriptive) is behind CLR-NS on attention but ahead on its token embeddings. On K562 control cells the direction is the same but smaller (+0.27, 0.09 to 0.46, descriptive). Fine-tuning each model on its own cells (descriptive) did not help: scPRINT's R0 fell slightly (−0.05 log2) and CLR-NS stayed ahead of it (+1.11, 0.41 to 1.96). Details: `docs/scfm_scoring.md`.

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

**9. tcaCLR: preregistered tests, a simulation and the confirmatory runs (17 to 19 November)**

*E. coli* test, 26 Oct: threshold not met. HepG2 test, 19 Nov: met, not broad. Kemmeren, 19 Nov: estimate only.

**Names, and what tcaCLR is.** This section uses the names of the short pre-read: **aCLR** is CLR with attention (CLR, softmax over each gene's strongest partners, lazy walk, t = 9), which the older sections and the figures of this page call "attention"; **tcaCLR** ("twice cooked aCLR") is aCLR followed by CLR applied to the log of the walk, which follow-ups C and D and the figures call "a3r" (a3r2 and a3r3 are its stacked versions). The names are temporary. tcaCLR keeps steps 1 to 4 of the method, replaces the walk mass P⁹ by its logarithm, z-scores each row of the log matrix over the pairs the walk reached, and combines the two directions as CLR does: CLR, softmax, walk, log, CLR, with one normalization layer on each side of the attention operator. It has no fitted parameters and the depth is unchanged. It was picked among about 20 read-outs by its K562 scores (K562 CORUM 7.95 against 6.39 for CLR and 3.51 for aCLR), so its K562 numbers are a selection and not evidence for it; that is a selection effect, not overfitting.

**Every look at tcaCLR on a complex truth.** The difference is pooled co-membership AUPR, tcaCLR minus aCLR, at t = 9.

| Data and truth                                      | Status                                                | tcaCLR minus aCLR              | Interval and p-value                                              |
|-----------------------------------------------------|-------------------------------------------------------|--------------------------------|-------------------------------------------------------------------|
| K562 Perturb-seq, CORUM                             | Selection among about 20 read-outs (in-sample)        | 7.95 against 3.51 (enrichment) | Not an estimate of the difference                                 |
| RPE1 Perturb-seq, merged CORUM                      | Descriptive                                           | +0.036                         | [+0.013, +0.059]                                                |
| *E. coli* M3D, Complex Portal                       | Preregistered test, 26 Oct                            | +0.026                         | p = 0.079; threshold (one-sided p < 0.025) not met               |
| *E. coli* EcoMAC, Complex Portal                    | Preregistered test, 26 Oct                            | +0.029                         | p = 0.051; threshold not met                                      |
| DREAM5 net4 yeast, Complex Portal yeast             | Exploratory estimation, 13 Nov (follow-up D)          | +0.096                         | Total +0.067 to +0.125; complex-level −0.023 to +0.215            |
| **HepG2 Perturb-seq, Complex Portal human (Arm A)** | **Preregistered confirmatory test, run once, 19 Nov** | **+0.0268**                    | **Total +0.0111 to +0.0424; one-sided p = 0.0004; threshold met** |
| Kemmeren deletion compendium, CYC2008 (Arm B)       | Preregistered secondary estimate, no p-value, 19 Nov  | +0.0566                        | Total +0.0392 to +0.0741; complex-level −0.009 to +0.122          |

The *E. coli* figures are from the 26 October test (follow-up C and the ledger below); the net4 estimate is follow-up D. Total interval: gene and experiment resampling. Complex-level interval: delete-half jackknife with complexes as the sampling unit, a sensitivity that never decides. The two complex tests with a p-value gave nearly the same point estimate, +0.026 and +0.029 on *E. coli* (threshold missed) and +0.027 on HepG2 (threshold met); see the reading below.

#### The preregistered *E. coli* test (26 October)

Complex superiority of tcaCLR over aCLR at t = 9 had to reach one-sided p < 0.025 on both M3D and EcoMAC. It was +0.026 AUPR on M3D (p = 0.079) and +0.029 on EcoMAC (p = 0.051), so the threshold was not met. Regulon non-inferiority (margin 0.005 AUPR) passed on both (+0.004 and +0.003). The design could detect an effect of about 0.05; the observed ones were about half of that. Details are in follow-up C.

#### Planted-module simulation (17 November; synthetic data; rules written before running)

120 simulated compendia with planted co-expression modules of 5, 10, 20, 40, 70 and 100 genes (three compendium sizes, N = 150, 400 and 800; Gaussian and heavy-tailed noise; two marginal types; 10 seeds per cell). The statements below are for the primary arm (Gaussian noise, identity marginals); no method setting was changed (`docs/plans/sim_planted_results.md`; independent reading in `docs/plans/fable_sim_reading.md`).

- **The gain edge is 20 genes.** tcaCLR beat aCLR in 10 of 10 seeds at module sizes 20 and 40 at every N, and in 0 of 10 at sizes 5 and 10, where the mean AUPR difference was −0.26 to −0.53. By the pre-written rule (the smallest size at which tcaCLR is ahead in at least 9 of 10 seeds, at each N) the lower edge of the gain is 20 genes, the same bin edge as in the net4 yeast data. At 70 and 100 genes the result depends on the seed (ahead in 7 to 10 of 10 seeds at 70 genes and 3 to 9 of 10 at 100), and tcaCLR lost to CLR at those sizes.
- **No mechanism claim.** The account in which tcaCLR's score is nearly independent of module size was refuted (its within-module score fell by 1.2 to 2.0 z per nat of module size, against fewer than 0.3 predicted), so the registered plan makes no mechanism claim.
- **A confinement collapse.** When attention is closed around a strong module, a member row reaches only the module plus 2 to 3 genes (in-module share 0.96 to 0.97), the row z-score of the log is then taken over the module itself and has no background, and the member score falls to about 0.6 against 5 to 7 in open modules (module AP 0.35 against 0.98 or more). This was shown on one compendium (size 70 at module strength λ = 0.85, size 100 at λ = 0.70), so there is no seed variability.
- **None on real K562.** A structure-only check on the K562 matrix (no truth scored; the data were already spent) found no closed complex: median reachable 5,654 of 5,670 genes for ribosomal, 20-plus, 7 to 19 and 3 to 6 gene complexes alike. The same descriptor was reported with both new runs below, and no class was closed there either.
- **Limits.** The generator is linear and Gaussian (so \|r\| wins), and it has no ribosome-specific structure, no experiment-level variance and no K562 selection effect.

#### The confirmatory runs: design

The plan (`docs/plans/tcaclr_confirm_plan.md`) was frozen on 17 November, before any score existed on these data; production runs were scored once on 19 November through the standard gates, with no tuning. On 17 November the PI withdrew the earlier "last complex test" declaration for this plan, so this is the last complex test of tcaCLR in this project, and every report of its result states the *E. coli* miss in the same sentence.

- **Arm A (primary, confirmatory):** HepG2 Perturb-seq pseudobulk (9,622 genes by 739 perturbation profiles; Nadig et al. 2025) with EBI Complex Portal human complexes: 359 merged complexes, 17,417 positive pairs among 2,040 truth genes, n_eff 34.0 (registered sensitivity +0.03). Hypothesis H1: tcaCLR's pooled AUPR exceeds aCLR's. Decision rule: one-sided p < 0.025 and an estimate of at least +0.005. Independence is limited: a new cell type and a new truth, but the same laboratory, CRISPRi essential-gene library family and pipeline as K562, and about 60% of the positive pairs are also CORUM positives (counted on a proxy gene universe before the real one existed). The claim is therefore cross-cell-type and cross-truth transfer within one lineage, not cross-technology. A built-in data-truth link: 418 of the 2,040 truth genes are themselves perturbed in the panel.
- **Arm B (secondary, estimate only, no p-value and no claim rule):** the Kemmeren yeast deletion compendium (6,169 ORFs by 1,283 mutant profiles; Kemmeren et al. 2014) with CYC2008 complexes: 213 merged complexes, 10,308 positive pairs among 1,324 genes, n_eff 10.4. The data are independent of net4, but 89.5% of the positive pairs are also positives of the Complex Portal yeast truth already used in follow-up D, so it is new data with largely the same complexes. 340 of the 1,324 truth genes are deleted in the panel.
- **Common to both:** the operator and settings as in the K562 and net4 work (t = 9); the total standard error combines a delete-half gene jackknife (500 replicates) and a delete-half experiment jackknife (100 replicates, CLR, aCLR, tcaCLR and \|r\| refitted); seven methods in each arm (tcaCLR, aCLR, CLR, \|r\|, GENIE3, PIDC, GRNBoost2). For GENIE3, GRNBoost2 and PIDC only gene resampling was run (100 refits would cost about 8 to 700 hours per method), so the required note applies to every interval and sentence that uses them: "error bars for GENIE3, GRNBoost2 and PIDC comparisons include gene resampling only and are narrower than they should be."

#### The confirmatory runs: results

| tcaCLR minus aCLR, main truth                                      | Arm A: HepG2, Complex Portal (test)            | Arm B: Kemmeren, CYC2008 (estimate)                                                   |
|--------------------------------------------------------------------|------------------------------------------------|---------------------------------------------------------------------------------------|
| AUPR: tcaCLR against aCLR                                          | 0.0669 against 0.0401                          | 0.208 against 0.151                                                                   |
| Estimate d                                                         | +0.0268                                        | +0.0566                                                                               |
| Total 95% interval (gene + experiment)                             | +0.0111 to +0.0424                             | +0.0392 to +0.0741                                                                    |
| SE gene, SE experiment                                             | 0.0060, 0.0052                                 | 0.0075, 0.0048                                                                        |
| One-sided p                                                        | 0.0004                                         | none (by registration)                                                                |
| Against the smallest effect, +0.005                                | interval above it                              | interval above it                                                                     |
| Complex-level interval (sensitivity)                               | −0.014 to +0.067 (includes 0)                  | −0.009 to +0.122 (includes 0)                                                         |
| Share of d carried by the top three complexes                      | 81%                                            | 82%                                                                                   |
| The three complexes                                                | 60S and 40S ribosomal subunits, 26S proteasome | Cytoplasmic large and small ribosomal subunits, mitochondrial large ribosomal subunit |
| Share of d in complexes of 20 or more genes                        | 89%                                            | 95%                                                                                   |
| Complexes with an improved AP (of those with 10 or more positives) | 113 of 211 (sign test p = 0.17)                | 51 of 96 (sign test p = 0.31)                                                         |
| Without the ribosomal complexes                                    | +0.0041 [−0.0016, +0.0099]                   | +0.0351 [+0.0150, +0.0552]                                                          |
| Pairs where neither gene is perturbed (or deleted)                 | +0.0139 [+0.0033, +0.0244]                   | +0.0485 [+0.0304, +0.0667]                                                          |
| Pairs involving a perturbed (or deleted) gene                      | +0.0393 [+0.0150, +0.0637]                   | +0.0156 [−0.0031, +0.0343]                                                          |

Intervals in brackets are the gene-plus-experiment total intervals. Arm A: one-sided p from the normal approximation with the total SE (0.0080). "Without the ribosomal complexes" removes the pre-declared ribosome list; the other splits are registered secondaries with no p-values. In Arm A the registered claim rule applied; in Arm B no claim rule applies, so its interval is shown next to Arm A's and nothing is claimed from it.

**Arm A met the registered test, and the gain is not broad.** The total interval lies above the +0.005 smallest effect and above 0 (p = 0.0004). Three complexes (the 60S and 40S ribosomal subunits and the 26S proteasome) carry 81% of the gain; the median complex is not improved (113 of 211 improve; the mean per-complex AP change is +0.0023); and if the complex, not the gene or the experiment, is the unit of resampling, the interval is −0.014 to +0.067 and includes 0. Without the ribosomal complexes the difference is +0.004 [−0.002, +0.010]. 89% of the gain comes from complexes of 20 or more genes. The gain is also present in pairs where neither gene is perturbed in the panel (+0.014 [+0.003, +0.024]), so it is not confined to pairs involving perturbed genes (which would make it a property of the screen design). Under Arm B the same three-ribosome concentration appears (82%), but the difference stays positive without the ribosomal complexes (+0.035) and in pairs where neither gene is deleted (+0.049); with only 10.4 independent complexes the complex-level interval is wide.

**The *E. coli* miss and the HepG2 hit have similar point estimates.** tcaCLR minus aCLR on complexes was +0.026 and +0.029 on *E. coli* (threshold missed, p = 0.079 and 0.051) and +0.027 on HepG2 (threshold met, p = 0.0004). The estimates are the same size; what differs is the precision (HepG2 has about 34 independent complexes by the registered count, *E. coli* about 4.4; the gene-resampling SE is 0.0060 on HepG2 against 0.018 on *E. coli*). The HepG2 result is therefore best read as the same effect measured more precisely in a new cell type, not as a reversal, and the registered headline carries both results together.

**Size pattern (H2, descriptive; no p-value anywhere).** The prediction, written after the net4 look and informed by the simulation, was a positive gain in complexes of 20 or more genes and none below 20 genes.

| tcaCLR minus aCLR by complex size | Arm A: HepG2                 | Arm B: Kemmeren                                               |
|-----------------------------------|------------------------------|---------------------------------------------------------------|
| 3 to 6 genes                      | +0.0003 [−0.0007, +0.0014] | −0.0040 [−0.0132, +0.0051]                                  |
| 7 to 19 genes                     | +0.0014 [−0.0024, +0.0052] | +0.0114 [−0.0055, +0.0284]                                  |
| of which 7 to 10 genes            | +0.0001 [−0.0005, +0.0008] | −0.0019 [−0.0254, +0.0215]                                  |
| of which 11 to 19 genes           | +0.0019 [−0.0037, +0.0075] | +0.0133 [−0.0083, +0.0348]                                  |
| 20 or more genes                  | +0.0354 [+0.0116, +0.0593] | +0.0734 [+0.0478, +0.0990]                                  |
| Script verdict                    | The prediction held          | The prediction partly held (2 of 3 size classes as predicted) |

Gene-plus-experiment intervals. In Arm B the 7 to 19 class is above the registered "not positive" limit (d at most 0.01) while its interval includes 0, which is why the verdict is "partly held". Complexes of 20 or more genes: 21 in Arm A and 9 in Arm B, so these strata rest on few complexes.

**Comparators (H3).** The estimand is the pooled AUPR difference against each other method; "best of the 7 methods run" would require every one of the six comparisons to be won at one-sided 0.025 with no correction.

| Method    | Arm A: AUPR | Arm A: tcaCLR minus method     | Arm A: Holm one-sided p | Arm B: AUPR | Arm B: tcaCLR minus method     |
|-----------|-------------|--------------------------------|-------------------------|-------------|--------------------------------|
| tcaCLR    | 0.0669      |                                |                         | 0.208       |                                |
| PIDC      | 0.0660      | +0.0009 [−0.0071, +0.0089]\* | 0.41                    | 0.142       | +0.0654 [+0.0406, +0.0903]\* |
| GENIE3    | 0.0592      | +0.0077 [+0.0009, +0.0145]\* | 0.052                   | 0.109       | +0.0984 [+0.0706, +0.1262]\* |
| CLR       | 0.0582      | +0.0087 [−0.0032, +0.0206]   | 0.15                    | 0.121       | +0.0865 [+0.0575, +0.1155]   |
| \|r\|     | 0.0516      | +0.0153 [+0.0002, +0.0303]   | 0.07                    | 0.075       | +0.1322 [+0.0906, +0.1738]   |
| aCLR      | 0.0401      | +0.0268 [+0.0111, +0.0424]   | 0.002                   | 0.151       | +0.0566 [+0.0392, +0.0741]   |
| GRNBoost2 | 0.0392      | +0.0277 [+0.0149, +0.0405]\* | < 0.0001               | 0.078       | +0.1296 [+0.0972, +0.1620]\* |

AUPR of the main truth, raw scores. Intervals are total (gene + experiment) except those marked \*, which are gene resampling only: error bars for GENIE3, GRNBoost2 and PIDC comparisons include gene resampling only and are narrower than they should be. Arm A: Holm-adjusted one-sided p over the six comparisons. Arm B: estimates and intervals only, no p-value and no best-of statement by registration; tcaCLR has the highest AUPR of the seven there. The registered trap: tcaCLR minus CLR has a smaller standard error than tcaCLR minus aCLR on *E. coli* (0.016 to 0.027 against 0.037 to 0.041) and can exclude 0 where the aCLR contrast does not; that is never a basis for "tcaCLR beats CLR". The HepG2 AUPR order is tcaCLR, PIDC, GENIE3, CLR, \|r\|, aCLR, GRNBoost2, and the differences to PIDC, GENIE3, CLR and \|r\| are small against the uncertainty.

**The registered sentences (verbatim, as the scoring scripts wrote them; the plan fixes them, so they are not paraphrased).**

Arm A, headline claim (registered claim rule, accompanied by the *E. coli* failure):

- "tcaCLR exceeds aCLR in pooled co-membership AUPR on the Complex Portal human complexes in HepG2 Perturb-seq pseudobulk (same laboratory and screen family as the discovery data; cross-cell-type and cross-truth transfer, not cross-technology); the gain is carried by complexes of 20 or more genes; the E. coli preregistered test of the same readout did not meet the one-sided 0.025 threshold (M3D p = 0.079, EcoMAC p = 0.051)."
- "The total interval excludes no difference: tcaCLR exceeded aCLR on complexes (preregistered confirmatory test)."
- "The whole interval lies above the smallest effect of interest (+0.005)."
- "The interval includes the E. coli estimate (+0.027)."
- "The pooled gain is carried by 3 complexes; tcaCLR did not improve the median complex."
- "H2: the prediction held (large positive, mid near zero, small not positive)."
- "tcaCLR is not shown to be the best of the 7 methods run: it does not exceed CLR, PIDC at one-sided 0.025 (intersection-union)."
- "tcaCLR exceeds aCLR (Holm-adjusted one-sided p = 0.0020, 6 comparisons)."; "tcaCLR is not shown to exceed CLR (Holm-adjusted one-sided p = 0.1503, 6 comparisons)."; "tcaCLR is not shown to exceed \|r\| (Holm-adjusted one-sided p = 0.0699, 6 comparisons)."; "tcaCLR is not shown to exceed GENIE3 (Holm-adjusted one-sided p = 0.0516, 6 comparisons)."; "tcaCLR is not shown to exceed PIDC (Holm-adjusted one-sided p = 0.4143, 6 comparisons)."; "tcaCLR exceeds GRNBoost2 (Holm-adjusted one-sided p < 0.0001, 6 comparisons)."
- "Note: error bars for GENIE3, GRNBoost2 and PIDC comparisons include gene resampling only and are narrower than they should be."
- "Trap: tcaCLR minus CLR has a smaller standard error than tcaCLR minus aCLR (on E. coli, k 0.016 to 0.027 against 0.037 to 0.041) and may exclude 0 where tcaCLR minus aCLR does not; this is never a basis for 'tcaCLR beats CLR' (E. coli S3)."
- "Arm A (primary): the claim rule of plan section 2 applies; the interval sentences above are descriptive of the same interval."
- "Perturbed-gene split: the gain is present in pairs where neither gene is perturbed in the panel."
- "Power gate: n_eff = 34.0 (overlap-counted); the registered sensitivity is +0.03."
- "Built-in data-truth link: 418 of the 2040 truth genes are themselves perturbed in the panel."

Arm B (secondary estimate):

- "The total interval excludes no difference: tcaCLR exceeded aCLR on complexes (preregistered secondary estimate; not a confirmatory test)."
- "The whole interval lies above the smallest effect of interest (+0.005)."
- "The whole interval lies above the E. coli estimate (+0.027)."
- "Arm B (secondary estimate): no p-value and no claim rule apply to this arm (plan section 2)."
- "Trap: tcaCLR minus CLR has a smaller standard error than tcaCLR minus aCLR (on E. coli, k 0.016 to 0.027 against 0.037 to 0.041) and may exclude 0 where tcaCLR minus aCLR does not; this is never a basis for 'tcaCLR beats CLR' (E. coli S3)."
- "The pooled gain is carried by 3 complexes; tcaCLR did not improve the median complex."
- "H2: the prediction partly held (2 of 3 size classes as predicted: small, large)."
- "Note: error bars for GENIE3, GRNBoost2 and PIDC comparisons include gene resampling only and are narrower than they should be."
- "Perturbed-gene split: the gain is present in pairs where neither gene is perturbed in the panel."
- "Arm B: n_eff = 10.4 (overlap-counted); estimate only by registration."
- "Built-in data-truth link: 340 of the 1324 truth genes are themselves perturbed in the panel."

#### Reachability descriptor (structure only; no score or truth compared)

The simulation showed that tcaCLR can collapse when attention is closed around a module, so the plan required the reachable-gene count of member rows by complex size class. The scorer did not emit it, so a separate structure-only script produced it afterwards (`analysis/tcaclr_confirm/reach_descriptor.R`). No class is flagged and no complex is closed on either arm: the median member row reaches 9,542 of 9,621 other genes on Arm A and 6,150 of 6,168 on Arm B, and the in-module share of reachable entries stays at the module's share of the universe (for example 0.0036 for complexes of 20 or more genes on Arm A, against 0.96 to 0.97 in the simulated collapse). The confinement collapse of the simulation did not operate on either dataset.

#### Qualifications, deviations and what is not claimed

- **Not broad, not best.** The registered test is met, and it is a statement about pooled AUPR on one frozen truth. The pooled gain rests on a handful of large complexes; the median complex is not improved on either arm; tcaCLR is not shown to be the best of the seven methods on HepG2 (PIDC, GENIE3, CLR and \|r\| are not exceeded after Holm correction), and the comparisons with GENIE3, GRNBoost2 and PIDC have understated error bars.
- **Independence.** HepG2 shares its laboratory and screen family with K562, where tcaCLR was selected; Kemmeren shares about 90% of its positive pairs with a truth already used. Both truths have built-in links to the data (418 of 2,040 and 340 of 1,324 truth genes are perturbed or deleted in the panel), as for any perturbation compendium.
- **Not claimed, whatever the result:** that tcaCLR is a theoretically justified or optimal estimator; a replication of the *E. coli* test; cross-technology transfer; generalisation to complexes outside the truth; any statement about mechanism; any statement about direct regulatory edges (CLR and GENIE3 stay ahead of tcaCLR there on *E. coli*, above); any statement about a method on a dataset where it was not run.
- **Deviations from the registered procedure, none changing a registered statistic.** The collect line of both experiment-jackknife jobs lacked `--allow_missing 1` and stopped before writing, so standalone collect jobs ran afterwards on the complete 100 replicates; the reachability descriptor was added as a separate script; the smoke-subset score gate had no reference value for a new arm; and the registered table of exploratory read-outs is only partly built (hub-stratified read-outs only), so nothing in it is reported here. n_eff was quoted in the plan as 31.69 (Arm A) and 10.02 (Arm B); the gate's own definition gives 33.97 and 10.35, which corrects a statement of the plan and changes no rule or decision.
- **Review status.** These results are the scripts' own outputs. They have not yet had a separate statistical review. The plan had an independent review before it was frozen (`docs/plans/tcaclr_confirm_review.md`).

#### Data and truth sources, and the CYC2008 licence

- HepG2 Perturb-seq: Nadig et al. 2025, *Nature Genetics* (PMID 40259084), GEO GSE264667.
- Yeast deletion compendium: Kemmeren et al. 2014, *Cell* 157:740 (PMID 24766815), GEO GSE42528.
- Complexes: EBI Complex Portal, human (taxon 9606) and yeast (559292) files; CYC2008: Pu et al. 2009, "Up-to-date catalogues of yeast protein complexes", *Nucleic Acids Research* 37:825–831 (doi 10.1093/nar/gkn1005, PMC2647312). The authors' server did not answer, so the catalogue file was taken from an Internet Archive capture of 16 February 2020 of the authors' file.
- **Licence.** The source article is open access under Creative Commons Attribution Non-Commercial 2.0 UK, which permits non-commercial use, distribution and reproduction with citation; the PI judges this fit for purpose. The catalogue files themselves carry no separate licence statement that was found, so the article's terms are taken to apply. Any paper must cite Pu et al. and state the capture.

Sources: `docs/tcaclr_confirm_results/armA/` and `armB/` (estimate, h2, h3, breadth, secondaries, descriptors, reach_descriptor_by_class, sentences.txt), `docs/plans/tcaclr_confirm_plan.md` (frozen plan and change log), `docs/plans/tcaclr_prereg_candidates.md` (datasets and licence), `docs/plans/sim_planted_results.md`.

## Exploratory follow-ups

Four analyses run on 11 to 13 November. None enters the paper as it stands. Letters A to D are used so they are not confused with items 1 to 9.

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

### C. A second CLR after the walk (a3r), and stacking it

a3r: preregistered 26 Oct, not adoptedStacking: exploratory, plan committed before compute (98f3b0d)

In this follow-up "a3r" is tcaCLR and "attention" is aCLR (names explained at the start of section 9). The confirmatory runs of 17 to 19 November are in section 9.

**The problem it targets.** Attention conserves mass: each gene's row sums to 1. Inside a tight complex of n proteins that are copies of one latent factor, every member therefore receives about 1/n of the factor's affinity, however strongly they co-vary. CLR has no such budget, so a ribosomal pair keeps its score regardless of complex size. That is our account of why attention loses to CLR on complexes (follow-up and diagnosis in `docs/plans/corum_analysis.md`).

**a3r in one paragraph.** Keep steps 1 to 4 of the method unchanged. In step 5, replace the raw walk mass M = P^(t) by its logarithm, which is the log-probability of reaching j from i in t steps (a multi-hop attention logit). Then z-score each gene's row of the log matrix against its own background, using only pairs the walk reached, and combine the two directions as in CLR (z_(ij) + z_(ji))/√2, clipped at 0. The pipeline is CLR, softmax, walk, log, CLR: one normalization layer on each side of the attention operator. The log turns the 1/n dilution into an additive offset that the row z-score removes. No parameter is added and t = 9 is unchanged.

    Pt = P^9;  L = log(Pt)  # unreachable pairs (Pt = 0) are set to missing
    for each gene i:  Z[i,j] = max(0, (L[i,j] - mean_j L[i,:]) / sd_j L[i,:])   # over reachable j only
    score[i,j] = (Z[i,j] + Z[j,i]) / sqrt(2)

**Stacking (a3r2, a3r3, in general a3rN).** a3r's output is a symmetric z-scale matrix, the same kind of object CLR hands to the softmax in step 3. So it can be fed back in: layer k + 1 is the same pipeline with the score S_(k) in place of the CLR score. The question was whether a second or third layer pools evidence across overlapping modules. Primary setting: t = 9 at every layer. Controls in the same run: a single layer at deeper t = 18 and 27 (the same total reach), a plain repeat of the CLR step on a3r's output, a shorter second layer (t = 2, 4), and hub-stratified scores.

| Read-out                                   | K562 CORUM | *E. coli* regulons, strong | *E. coli* regulons, all | *E. coli* complexes |
|--------------------------------------------|------------|----------------------------|-------------------------|---------------------|
| Attention (t = 9)                          | 3.51       | 1.68                       | 1.44                    | 22.6                |
| CLR                                        | 6.39       | 1.44                       | 1.37                    | 21.8                |
| a3r, one layer                             | 7.95       | 1.72                       | 1.56                    | 26.1                |
| a3r2, two layers                           | 8.19       | 1.73                       | 1.51                    | 31.0                |
| a3r3, three layers                         | 8.15       | 1.66                       | 1.46                    | 32.5                |
| a3r at t = 18 (same reach as two layers)   | 8.00       |                            |                         |                     |
| a3r at t = 27 (same reach as three layers) | 7.86       |                            |                         |                     |

Pooled AUPR divided by the base rate (1.00 = random), raw scores, one run each. K562: essential Perturb-seq, all CORUM complexes. *E. coli*: M3D with RegulonDB regulons (strong evidence, all evidence) and the complex truth of the 26 October test (1,150 positive pairs). Matched-depth rows were scored on K562 CORUM only in this table.

#### Where a3r wins and loses on every metric we have (descriptive, page 3)

- **Ahead of attention and CLR:** K562 CORUM (7.95 against 3.51 and 6.39; \|r\| is 7.00), K562 DoRothEA A-B (1.11 against 1.09 and 1.06) and CollecTRI (1.13 against 1.10 and 1.03), *E. coli* regulons on M3D (strong 1.72, all 1.56 against attention 1.68 and 1.44) and EcoMAC (1.72 and 1.58 against 1.68 and 1.47), and *E. coli* complexes (26.1 against 22.6 and 21.8; M3D, tuning set). On BEELINE's primary truths, with the CLR-NS operator at t = 9, the median ratio is 1.73 against attention 1.54, CLR-NS 1.60 and GENIE3 1.69 (exploratory; not tuned).
- **Behind or level:** direct TF-to-target edges on M3D, where a3r (5.84) is below CLR (6.75) and GENIE3 (6.54) and above attention (5.17); per-regulon mean AP, where a3r (25.2) is level with attention (25.6) and above CLR (23.0); RPE1 merged CORUM, where CLR (6.77) is above a3r (6.14), which is above attention (3.33); and BEELINE cell-type ChIP truths, where every method is about 1.0.
- **Stacking trades away the other metrics.** a3r2 and a3r3 fall well below a3r on direct edges (3.95 and 3.30), and below it on K562 knockdown response, DoRothEA and CollecTRI, for a small gain on complexes. That is the reason a3r, not a3r2, is the candidate.
- Knockdown response (a3r 3.13 against attention 2.38) comes from a truth with near-duplicate modules and is not a clean win.

#### What the preregistered a3r test showed (26 October)

- Complex superiority over attention at t = 9 was required at one-sided p < 0.025 on both *E. coli* datasets. It was +0.026 AUPR on M3D (p = 0.079) and +0.029 on EcoMAC (p = 0.051), so it did not confirm. The jackknife standard error on complexes was 0.018, and the smallest effect the design could detect was about 0.05.
- Regulon non-inferiority (margin 0.005 AUPR) passed on both: +0.004 and +0.003.
- Verdict: not adopted. Attention stays the method. The K562 numbers in the table (7.95 against CLR 6.39) are selected, not fitted: a3r has no fitted parameters, so it cannot memorise the truth, but it was picked as the best of about 20 read-outs by its K562 scores, and the best of many noisy scores is biased upward. The 20 candidates were highly correlated and the gain is large, so the bias is probably a fraction of the gap, but it has not been measured. The effect points the same way on RPE1 and *E. coli*, neither confirmatory. On RPE1, another cell line and a descriptive look only, a3r is +0.036 [+0.013, +0.059] AUPR above attention on merged CORUM.

#### What stacking showed (exploratory)

- Neither two nor three layers meets the rule written before the run, on K562 or on *E. coli*, with either hub stratifier. The rule required a K562 CORUM gain of at least +0.005 over a3r and over the matched-depth control, with intervals excluding 0, and no loss over 0.005 on DoRothEA or CollecTRI.
- K562: a3r2 is +0.0031 [−0.0027, +0.0089] over a3r on CORUM and +0.0025 [−0.0036, +0.0087] over a single layer at t = 18. It loses on CollecTRI (−0.006), knockdown response (−0.044) and CORUM without ribosomal complexes (−0.003). Three layers lose more (CollecTRI −0.011, knockdown response −0.079).
- *E. coli*: a3r2 is non-inferior on strong-evidence regulons (+0.0006) but not on all-evidence regulons (−0.0055, lower bound −0.009). Three layers fail both (−0.006 and −0.012). Complex AUPR rises (+0.036 for two layers, +0.047 for three), but the interval for two layers includes 0, a single deeper layer gives only +0.001, and hub stratification roughly halves the gain. The truth is small.
- A shorter second layer is much worse (t = 2: 1.63 on K562 CORUM; t = 4: 3.10). A plain repeat of the CLR step on a3r's output is not idempotent and is worse than a3r (5.85 and 6.83 for the two variants of the repeat).
- What happens to the operator: nothing collapses, but the share of pairs the walk has not reached grows from 0.28% (one layer) to 5.7%, 12.3% and 17.2% (K562, layers 2 to 4), and the stationary distribution drifts away from layer 1's (Spearman 0.74, 0.52, 0.40).

**Summary.** A single layer already carries the K562 complex gain over attention. A second layer adds a gain that is small, inside the noise and no larger than the same reach from a single deeper layer, and it costs the regulon and knockdown truths. Read at the point estimates alone, a3r2 leads on K562 CORUM (8.19) and strong-evidence regulons (1.73) and the a3r family beats attention and CLR on both. The a3r2-over-a3r gaps (+0.003, +0.0006) are inside the noise, though, and a3r2 gives up all-evidence regulons, CollecTRI and knockdown response, so a3r, the simpler form, is the candidate. It has not been confirmed (K562 and M3D are both tuning sets for it), and stacking it is not a lever. The depths tried were t = 9 at every layer, with a second layer at t = 2, 4 and 9 only. Layers longer than 9 were not run. Full tables: `docs/explore_stack_results.md`.

### D. a3r on yeast protein complexes: an exploratory estimate

Plan fixed before compute (design rev 2, 099c26d)Exploratory estimation, not a confirmatory testRun once

**Question.** How large is the pooled AUPR difference between a3r and attention on protein complexes, in a dataset and truth where a3r had never been run? The 26 October test on *E. coli* could not confirm (+0.026 and +0.029, p = 0.079 and 0.051) and was underpowered. A new confirmatory test was designed first. Its projected precision (total standard error about 0.020) was too low to confirm a plausible effect, so on your decision of 13 November it was run as an estimation study instead: no p-value, no confirmed or not confirmed label, and a fixed set of sentences chosen by a script before the run.

**Design in brief.** Expression: the DREAM5 network 4 yeast compendium (536 arrays, 5,950 genes), which is the matrix of the section 8 test, so it is not new to the project; a3r had never been run on it (checked by a guard before the run). Truth: EBI Complex Portal yeast complexes, 262 merged complexes and 13,438 positive gene pairs, frozen with a checksum before any score, about 14.5 effectively independent complexes (*E. coli*: 54 complexes, about 4.4). The two cytoplasmic ribosome subunits are kept as separate complexes, as in the source, and carry 28.9% of the positive pairs. Estimand: pooled AUPR of a3r minus attention at t = 9, with a3r minus CLR and attention minus CLR alongside. Reported interval: gene jackknife plus experiment-series jackknife (59 series). Every implementation check passed, including exact reproduction of the K562 and *E. coli* numbers and of the stored yeast attention and CLR scores. Full plan: `docs/plans/complex_confirm_design.md`; independent review before the run: `docs/plans/complex_confirm_review.md`.

| Contrast (AUPR difference, main truth) | Estimate | Total interval (gene + series) | Gene-only interval | Complex-level interval (sensitivity) |
|----------------------------------------|----------|--------------------------------|--------------------|--------------------------------------|
| a3r minus attention                    | +0.096   | +0.067 to +0.125               | +0.075 to +0.118   | −0.023 to +0.215                     |
| a3r minus CLR                          | +0.039   | +0.010 to +0.067               | +0.028 to +0.049   | −0.026 to +0.103                     |
| attention minus CLR                    | −0.058   | −0.098 to −0.018               | −0.076 to −0.040   | −0.122 to +0.007                     |

Pooled enrichment over random: CLR 18.1, attention 12.6, a3r 21.7 (base rate 0.0105). Complex-level interval: delete-half jackknife with complexes as the sampling unit; sensitivity only. The total interval is the reported one. The a3r minus CLR contrast has a smaller standard error than a3r minus attention, so its gene-only interval excludes 0 more easily; that is not evidence that a3r beats CLR.

**The sentences the script wrote (verbatim, total interval):**

- "The total interval excludes no difference: a3r exceeded attention on yeast complexes (preregistered exploratory estimation; not a confirmatory test)."
- "The whole interval lies above the smallest effect of interest (+0.005)."
- "The whole interval lies above the *E. coli* estimate (+0.027)."
- "Exploratory estimation, not a confirmatory test, a replication or a validation. a3r is not adopted; this study cannot change that by itself (PI decisions 12 and 13 Nov)."

#### Qualifications, from the same outputs

- **Complexes as the sampling unit.** If the unit of generalisation is the complex and not the gene or the experiment, the interval is −0.023 to +0.215 and includes 0. With 14.5 effectively independent complexes that interval is wide, and it was a sensitivity line, not the reported one.
- **Where the gain is.** Complexes of 20 or more genes: +0.137 [+0.106, +0.168]. Complexes of 3 to 6 genes: −0.004 [−0.011, +0.004]; 7 to 19 genes: −0.003 [−0.009, +0.002] (gene-only intervals). Without the ribosomal complexes: +0.046 [+0.025, +0.066]. Per complex (125 complexes with at least 10 positive pairs), the mean AP change is +0.0018, the median is about 0, and 60 improve while 65 worsen. So the gain is about how large complexes compare with the rest of the genome, not about ranking within a complex; that is what the mass-conservation account predicts, since the 1/n dilution is largest in large complexes. It is consistent with that account, not a test of it.
- **The ribosome rule matters.** The ribosome subunits are separate complexes in the source. With cross-subunit pairs removed, the difference is +0.191; with the two subunits joined it is +0.241 (effective independent complexes fall to 4.9). Both were descriptive sensitivity truths, fixed before the run. The K562 CORUM truth, where a3r was selected, merged the subunits.
- **Direction elsewhere.** Attention is below CLR on yeast complexes (−0.058), in the same direction as K562 CORUM. On the touched DREAM5 regulon truth (descriptive) a3r minus attention is +0.011 [+0.007, +0.015]. Ribosomal-protein module coherence (AUROC): CLR 0.974, a3r 0.968, attention 0.953.
- **Realised precision.** The total interval's half-width is 0.029, against 0.040 projected.
- **What it uses up.** The Complex Portal yeast truth is now consumed for any later complex candidate. A code review of the code that ran was done after the run, not before it as the protocol required, and the reviewer saw the headline estimate by accident early on, so it is not fully blind. It found no bug that could change a reported number (two process findings, and minor robustness and gate-coverage gaps); `docs/complex_confirm_code_review.md`.

**Summary.** On yeast, a dataset and truth where a3r was never selected or tuned, a3r's pooled AUPR on complexes is above attention's by an estimated +0.096, with the interval for gene and experiment resampling well above both 0 and the *E. coli* estimate. The gain is concentrated in large complexes, is absent in small and mid-size ones, and does not survive treating complexes as the unit of resampling. It is an estimate in a second organism, not a confirmation, and a3r is not adopted.

## Every preregistered test so far

Each outcome is taken from the verdict document of that test (`docs/clr_logpt_verdict/`, `docs/keff_verdict/`, `docs/scfm_scoring.md`, `docs/reciprocity_yeast_selections/`, `results/dream5/truth_net4/`, `docs/depth_log_results.md`, `docs/tcaclr_confirm_results/`), in the order they finished.

| Finished | Test (fixed design) | Outcome |
|---|---|---|
| 26 Oct | CLR applied to log P^t as a fix for protein complexes; confirmatory on *E. coli* M3D and EcoMAC | **Not adopted: failed to confirm.** Complex superiority over attention at t = 9 needed one-sided p < 0.025 on both datasets: M3D +0.026 AUPR (p = 0.079), EcoMAC +0.029 (p = 0.051). Regulon non-inferiority passed. Attention stays the method. (The later tcaCLR runs are the rows of 13 and 19 Nov and item 9.) |
| 27 Oct | Attention width: 8 variants, 32 primary tests, Holm-corrected | **No variant meets the decision rule;** the headline width (10 effective partners) is unchanged. K25 and O3A each won 2 of their 4 tests and O3B 1, none enough to meet the rule. |
| 2 Nov | Yeast reciprocity | **Not run.** The data-quality gate failed in all 10 draws, so by the preregistration no truth was opened. |
| 2 Nov | CLR-NS against zero-shot scPRINT attention, 7 BEELINE datasets | **Non-inferior and superior** (Hodges-Lehmann +1.03 log2 AUPRC ratio, exact interval +0.34 to +1.90, p = 0.0078, positive in 7 of 7). The comparator scores at the level of its own permutation control, so it is a weak comparator. |
| 11 Nov | DREAM5 yeast, network 4 (Y-M and Y-E) | **Not replicated at the stated margin.** Y-M: Equivalent within 10%, on the boundary with Inconclusive, not Superior. Y-E: Inconclusive. |
| 11 Nov | Log-scale depth re-test (not blind) | **Reading (b), "not resolution-like";** with 3 and 5 truths the test could neither establish nor exclude the ordering. |
| 13 Nov | a3r against attention on yeast complexes (DREAM5 net4 expression, Complex Portal); plan fixed before compute, run once; **estimation, not a test** | **Estimate +0.096 AUPR**, total 95% interval +0.067 to +0.125 (gene + experiment series); the complex-level interval, −0.023 to +0.215, includes 0. Gain carried by complexes of 20 or more genes. a3r stays not adopted (follow-up D). |
| 19 Nov | tcaCLR against aCLR on human complexes: HepG2 Perturb-seq pseudobulk with Complex Portal human (Arm A); plan frozen 17 Nov, run once; **confirmatory** | **Preregistered test met:** +0.0268 AUPR, total 95% interval +0.0111 to +0.0424, one-sided p = 0.0004, above the +0.005 smallest effect. The gain is not broad: three complexes carry 81%, the median complex is not improved (113 of 211 improve), the complex-level interval, −0.014 to +0.067, includes 0, and tcaCLR is not shown to be the best of the 7 methods run (item 9). |
| 19 Nov | The same estimand on the Kemmeren deletion compendium with CYC2008 (Arm B); plan frozen 17 Nov, run once; **secondary estimate, not a test** | **Estimate +0.0566 AUPR**, total 95% interval +0.0392 to +0.0741; no p-value by registration. Ribosomal complexes carry 82% of the gain; the complex-level interval, −0.009 to +0.122, includes 0 (item 9). |

**Many looks.** Rows 1 to 6 and the HepG2 row of 19 November are tests with a fixed design (one of them was not run and one, the depth re-test, was not blind); two were positive (scPRINT and HepG2), and the 26 October complex test and the DREAM5 yeast test did not reach their thresholds. The 13 November yeast row and the Kemmeren row are estimation studies, not tests. The *E. coli* headline of items 1 to 4 is not one of them; it comes from development on M3D. Around them sit dozens of exploratory analyses on the same M3D/RegulonDB and K562 truths (FDR operators, value transforms, multiscale algorithms, width rounds, rCLR, robust standardisation, macro averaging, hub controls, granularity). The nested bootstrap protects the depth choice. The other settings (10 effective partners, 50 candidates, lazy weight 0.5, the MI estimator, the truth construction with operon exclusion) were settled during development on the same data. Only the yeast tests (section 8, the net4 estimate of follow-up D, which used the same expression matrix as section 8, and the Kemmeren estimate) and the HepG2 test lie outside that lineage. The HepG2 data are a new cell type with a new truth, but come from the same laboratory and screen family as K562, where tcaCLR was selected; the Kemmeren truth shares about 90% of its positive pairs with the Complex Portal yeast truth already used.

## What did not work (reported as negatives)

- Hard FDR-thresholded operators: they isolate 62–77% of genes.
- Yeast reciprocity test: not run; its preregistered data-quality gate failed in all 10 draws (median chance-corrected reproducibility 0.12 for mutants against 0.30 for wild types), so no truth was opened. Reciprocity stays a synthetic-only candidate (Supplementary Note 18).
- Value transforms and diffused-profile readouts.
- Wider gene contexts on single cells.
- Data-driven choice of attention width: all three options failed review before running, and were tested empirically anyway under a preregistration (finished 27 October). The headline width (10 effective partners) is unchanged and no variant meets the decision rule. A separate exploratory width study, on synthetic and label-free data with no truth opened, closed negatively on 31 October.
- A label-free multiscale algorithm based on spectral plateaus: there are no plateaus; the structure is a continuum.
- CLR applied to log P^t as a fix for complexes: preregistered, failed to confirm (26 October; complex superiority p = 0.079 on M3D and 0.051 on EcoMAC; regulon non-inferiority passed). Attention stays the method.
- Depth as a discrete "resolution dial" ordering operons → regulons → sigmulons: not supported on a linear depth grid, and the preregistered re-test on a geometric grid to t = 256 did not find the fine-to-coarse ordering either (pre-specified reading (b); with 3 and 5 truths it could not have established or excluded it, follow-up A). Across truths, the best depth is not clearly related to median module size (slope 0.20; module-bootstrap interval 0.08 to 0.67, which ignores variation between truths; across the 7 truths −0.77 to 1.16; the sign depends on the operon truth). Within regulons, an earlier per-set analysis on the linear grid found best depth rising with set size (t ∝ n^0.43) but not for random sets.
- tcaCLR on HepG2 (19 November): the preregistered test is met, but the gain is not broad. Three complexes carry 81% of it, the median complex is not improved, the complex-level interval includes 0, without the ribosomal complexes the difference is +0.004 [−0.002, +0.010], and tcaCLR is not shown to be the best of the 7 methods (PIDC, GENIE3, CLR and \|r\| are not exceeded after Holm correction; item 9).
- tcaCLR on small planted modules (simulation, 17 November): it was not ahead of aCLR in any of 10 seeds at 5 and 10 genes (mean AUPR difference −0.26 to −0.53), collapsed when attention was closed around a strong module, and its size-invariance account was refuted; it was seed-dependent and below CLR at 70 and 100 genes. Real K562, HepG2 and Kemmeren data show no closed complexes (item 9).
- Stacking a3r (a second and third attention layer on the log-then-CLR score): exploratory, 12 November. Neither two nor three layers meets the rule written beforehand on K562 or *E. coli* (follow-up C); the K562 complex gain over one layer is +0.003 [−0.003, +0.009] and a single deeper layer matches it.

## In progress and recently completed
- Learned attention (done 2 Nov; see item 5 above).
- GENIE3, GRNBoost2 and PIDC on the BEELINE single-cell benchmark: done. CLR is comparable to GENIE3 and PIDC. GRNBoost2 (legacy early stopping) is below CLR in 6 of 7 datasets (mean log2 −0.15; Wilcoxon P = 0.031, sign P = 0.125; hESC is the exception), a descriptive count, not a significant difference.
- Preregistered attention-width test, including multi-head attention as the candidate fix for the complex loss: done (27 October; 8 variants, 32 primary tests, Holm-corrected). Headline unchanged at 10 effective partners; no variant meets the decision rule. Multi-head attention did not gain on K562 CORUM (−0.003 AUPR against the headline, Holm P = 0.54). The write-up and the list of individual tests that variants won are still to do (`docs/keff_verdict/`).
- Directed DREAM5 evaluation: done (item 8; ARACNE was not part of the design). The write-up is in progress.
- tcaCLR confirmatory runs (HepG2 Arm A, Kemmeren Arm B): done, scored once on 19 November (item 9), with the planted-module simulation and the real-K562 reachability check of 17 November. Not done: a separate statistical review of the item 9 results; the registered table of exploratory read-outs, which is only partly built; and the decision on what to write.
- a3r on yeast complexes: done, one run (follow-up D), with a post hoc code review (no bug found; not fully blind). Stacked a3r layers: done (follow-up C). Log-scale depth analysis: done (follow-up A). Robust row-standardisation, the PI's median/IQR idea: done, exploratory (follow-up B).
- Still planned: perturbation-diversity learning curves, differential CLR, and a hub stratifier for the method table that does not depend on the walk (follow-up B's MI row mean is one, but it is tied to row centring and is not neutral).

## Decisions for you and me
1. **Venue:** Nature Machine Intelligence (150-word abstract) or Nature Communications.
2. **Scope:** one paper on CLR attention, with perturbation targets and differential CLR as a second paper.
3. **The central claim:** "better recovery of regulatory programs from expression alone", not "better network inference". Is that strong enough for the target venue? The statistical review suggests scoping the claim to *E. coli*, or presenting the work as a reframing of CLR with characterised trade-offs, since every positive result comes from the development organism and lineage.
4. **The human single-cell story:** CLR-NS rather than attention is the single-cell result. Does it belong in this paper?
5. **The yeast result:** Y-M has a positive estimate under the 10% margin (the verdict "Equivalent" is on the boundary with Inconclusive), Y-E is inconclusive, and the network 3 pattern (module gain, edge cost) is not independent of the headline. How should this enter the paper and the abstract, which is currently scoped to *E. coli*?
6. **What to do with tcaCLR (the headline method).** The evidence now: selected on K562 (in-sample); the preregistered *E. coli* test missed its threshold (+0.026, p = 0.079; +0.029, p = 0.051); an exploratory yeast estimate of +0.096, carried by large complexes, with a complex-level interval that includes 0; the preregistered HepG2 test met its threshold (+0.027, p = 0.0004; Kemmeren estimate +0.057), with three complexes carrying 81% of the gain, the median complex not improved and tcaCLR not shown to be the best of seven methods; and behind CLR and GENIE3 on direct edges in *E. coli*. Is the headline method aCLR, tcaCLR, or a presentation of both with their trade-offs? The earlier decision (not adopted, 12 November) stands until you change it; the short pre-read calls tcaCLR the lead candidate. The names are temporary.
7. **How to word the tcaCLR complex result.** The registered HepG2 sentence states the result together with the *E. coli* miss and the 20-gene carriers. Still to settle: whether the paper presents HepG2 as a confirmation of the *E. coli* effect, whose estimate it matches, or as a separate fixed-truth result; how much weight the breadth qualifications carry (three complexes, median complex not improved, complex-level interval includes 0, not shown best of seven); and how the Kemmeren estimate, which has no p-value, is shown next to it.
8. **AI disclosure:** the analysis code, runs and drafts were produced with Claude (Anthropic), with independent AI statistical reviews, under my direction. A contributions statement is drafted to match Nature Portfolio policy.

Answered since the last version: whether the yeast estimate was the last evaluation of tcaCLR on a new complex truth. On 17 November the PI withdrew the 1 November "last test" declaration for the HepG2 plan, which is therefore the last complex test of tcaCLR in this project, confirmatory or exploratory.
