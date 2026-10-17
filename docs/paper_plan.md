# Paper plan and results checklist

*Living document, started 2026-10-16. Organized by theme; every analysis we try gets a line, including the ones that did not work. Numbers are headline values. Details are in `docs/whitepapers/paper.md`, `supplement.md` and `paper_draft_notes.md`, and in `results/`.*

**Legend:** ✅ worked, ❌ did not work (kept as a reported negative), ⚠️ mixed or weaker than hoped, ⏳ running or queued, ☐ not started. "QC" marks results checked by an independent reviewer agent.

## Two papers

| | Paper 1: completed CLR attention | Paper 2: differential CLR for perturbation targets |
|---|---|---|
| Question | Does completing CLR as attention recover regulatory programs better, why, and at what resolution? | Which genes does a perturbation rewire? Targets are genes whose CLR context changes between conditions. |
| Data | M3D, EcoMAC; human Perturb-seq pseudobulk compendia; BEELINE (secondary) | M3D perturbation groups; Perturb-seq single cells (knockdown cells against controls) |
| Target venue | Nature Machine Intelligence if the multiscale variación and competitive baselines land, otherwise Nature Communications or Genome Biology | Genome Biology, Nature Methods or Cell Systems |

---

## Paper 1

### 1. Completing CLR as attention (M3D, E. coli)

- [x] ✅ Softmax attention plus depth beats CLR on regulon co-membership: +0.021 AUPR (one lab split), +0.022 [+0.017, +0.027] over 50 lab splits. QC
- [x] ✅ Nested bootstrap, depth re-selected in each of 1,000 draws: +0.022 [+0.019, +0.026], ahead in 1,000/1,000 draws; coherence +0.026 [+0.008, +0.045]. QC
- [x] ✅ Replicate-aware bootstrap on 907 chips (cluster 1,000 and replicate 200 draws): +0.020 [+0.018, +0.023].
- [x] ✅ Held-out depth is stable: median t* = 9 over 50 splits, criterion flat within 1% for t = 7-11. QC
- [x] ✅ Softmax candidate cap does not drive the result (caps 50/100/200 equal; no cap gains less AUPR). QC
- [x] ⚠️ All-evidence coherence gain is not significant in the nested bootstrap (+0.011 [-0.004, +0.025]). QC
- [x] ❌ Hard-thresholded FDR operator alone: 62-77% of genes isolated, unusable for propagation.
- [x] ❌ FDR graph plus top-10 neighbours (hard selection): +0.013 AUPR, no coherence gain.
- [ ] ⏳ FDR operators inside the chips bootstrap, per-draw null B = 100 (job 100).
- [ ] ☐ Repeat the core result on EcoMAC with held-out depth and bootstrap intervals (point estimate: +0.032 AUPR).

### 2. Mechanism: which components matter

- [x] ✅ Context normalization is essential: softmax on uncalibrated |r| or MI gets about a third of the AUPR gain and no coherence gain (kernel ablation at held-out depths). QC
- [x] ✅ Shapley decomposition over six components: depth +0.018 and CLR +0.013 (AUPR); CLR +0.025 (coherence); removing the softmax costs most coherence (-0.043).
- [x] ✅ Stouffer and Euclidean combination equivalent (Spearman 0.993).
- [x] ✅ Common low-bias bins beat per-gene Freedman-Diaconis bins (small).
- [x] ✅ CLR-score permutation null stays sparse under global confounders where an MI null does not.
- [x] ✅ Permutation null of B = 100 is enough (B = 1,000 gives the same edge sets).
- [x] ❌ Value transforms (signed value, conditional-expectation W_V): work on synthetic data, no gain on real data.
- [x] ❌ Diffused-profile readout |cor(P^t Z)|: little gain (best +0.007).
- [x] ❌ Rank (copula) transform before MI on bulk data: slightly worse.

### 3. Robustness and confounders

- [x] ✅ Attention gain survives removal of 0-20 global principal components (M3D +0.017 to +0.022; EcoMAC +0.020 to +0.029).
- [x] ✅ The top components are biology, not noise: growth rate (ribosomal proteins, AUROC 0.98), motility (σ28, FlhDC) and σS stationary phase.
- [x] ❌ `sva::num.sv` to choose how many components to remove: over-removes (k = 31 and 68).
- [x] ❌ Held-out reproducibility as the criterion for batch correction or k: always prefers keeping reproducible global programs.
- [ ] ☐ Blocked (per-dataset) tests wherever several comparisons share one dataset.

### 4. Modularity and resolution (the final variación)

- [x] ✅ Hints: operons gain sharply at t = 2; coherence peaks at t = 4-8 while pooled AUPR rises to about t = 27; attention helps TF regulons and hurts complexes in K562.
- [x] ⚠️ Granularity sweep (jobs 0873, 0870b; QC): only operons show a sharp interior optimum (AUPR 0.30 → 0.50 at t = 4-5, falling beyond t ≈ 8). Coarser truths are flat (AUPR changes ≤ 0.02 over t = 8-32) with edge-censored optima; no fine-to-coarse ordering (E. coli argmax 4/32/13; K562 0/26/4/7/32, p = 0.12). Dropping σS moves the sigmulon argmax 24 → 13, so single sets drive it. CORUM: no depth helps. Do not claim "depth is a resolution parameter" yet.
- [ ] ☐ Granularity follow-ups before any resolution claim: hub-baseline-adjusted curves; t to about 128 on a log grid; set-level bootstrap (and leave-one-set-out) CIs on band and gain-weighted centroid; size- and base-rate-matched random-set controls; tie-randomized AP or nonzero-support scoring (S_1 is mostly exact zeros); preregistered ordering test on held-out truths.
- [ ] ☐ Hub-bias check: score the degree baseline π_i + π_j against the truths; the pooled AUPR rise to t ≈ 27 could partly be convergence to hubs.
- [ ] ☐ Pilot (no bootstrap): trace(P^t) and its log-slope on a log-spaced t grid, for keff ∈ {3, 5, 10}. Plateaus = regimes, plateau height = module count. Stop here if there are no plateaus.
- [ ] ☐ Multiscale algorithm (label-free; design review 2026-10-16): per-gene module size k_i(t) = G·π_i / P^t_ii (no global k); membership by lift P^t_ij / π_j with a lowest-conductance sweep cut; regimes by kernel change-point detection on log-slopes over 50 lab half-splits; overlapping modules = consensus sweep-cut sets reproduced in ≥ 70% of splits; collapse when trace(P^t) − 1 < 1 or reproducibility reaches a rewired-operator null. Track births and deaths with `ModuleTester`.
- [x] ❌ Smarter keff, reviewed before running (2026-10-16): one global temperature on CLR z gives a bimodal width split (29% of rows < 3 partners, 16% near-uniform) tracking the top score, not regulon breadth, and loses on held-out data; per-gene local-FDR keff fails its assumptions (point mass at 0, heavy tail, dependent rows); held-out selection of keff is biased toward wide attention (keff 25 > 10 > 3 on one split), echoing the batch-correction negative.
- [ ] ☐ Try all three smarter-keff options empirically anyway (user, 2026-10-16), next to the pilot's keff 3/5/10 sweep: (1) one global temperature on CLR z, (2) per-gene keff from local FDR, (3) keff chosen on held-out data. Compare against keff = 10 and the multi-head operator; the statistics reviewer weighs in on the design before running and on the results after.
- [ ] ☐ Width as a resolution axis: keff ∈ {3, 5, 10, 25} × t surface, preregistered; truths only at the end, stratified by granularity (prediction: small keff wins on operons/CORUM, large on sigmulons/CollecTRI). keff = 10 stays the headline.
- [ ] ☐ Multi-head attention (the NLP advance): average the walks built at keff = 3, 10 and 30 (multiscale SNE / heads at different temperatures); no width selection. Preregistered against keff = 10 with paired jackknife, Holm-corrected, on M3D, EcoMAC and K562.
- [ ] ☐ Comparison with Markov stability (PyGenStability) on the same operators.

### 5. Performance against competitors

- [x] ⚠️ GENIE3 and GRNBoost2 on M3D and EcoMAC (jobs 0874, 0875; QC): attention has the highest pooled regulon co-membership AUPR (+0.010 to +0.017 over GENIE3, +0.008 to +0.028 over CLR, jackknife intervals exclude 0), but the gain is concentrated in global regulons: CRP, σS, Fur, IHF and FNR carry 62% of positive pair incidences. GENIE3 wins median per-regulon coherence and per-TF AUROC (EcoMAC 0.690 vs 0.596); attention is last of the top three on direct TF-gene pairs (TF-node AUPR ratio M3D: CLR 6.75, GENIE3 6.54, attention 5.17). Claim: better module recovery, not better GRN inference. GRNBoost2 is sparse (87-93% tied zeros), not like-for-like.
- [ ] ⏳ GENIE3 + GRNBoost2 on K562 rerun with duplicated symbols dropped, all rows regenerated (jobs 0878, 0879).
- [ ] ☐ Baseline follow-ups (QC): score PIDC on M3D; directed DREAM5-style TF→target evaluation (regulators = RegulonDB TFs, GENIE3 unsymmetrized); ARACNE; paired per-regulon Wilcoxon on coherence stratified by size; pooled AUPR without the top-5 regulons; intervals for coherence and TF metrics; regulon-cluster resample; 3 seeds for tree methods; Holm over the prestated family; t re-selected per dataset (t = 9 reported as transfer); check M3D/EcoMAC array overlap before calling EcoMAC a replication.
- [ ] ⏳ PIDC on the same three datasets (second job lane, single-threaded: about 7 h each for M3D and EcoMAC, about 40 h for K562).
- [ ] ☐ ARACNE (the classic MI method).
- [ ] ☐ Learned attention from scGPT or Geneformer on a shared task.

### 6. Bacteria to human

- [x] ✅ K562 essential perturbation compendium: attention beats CLR on DoRothEA and CollecTRI (+0.008/+0.009 AUPR, significant).
- [x] ❌ Attention on CORUM complexes: -0.041 AUPR (fine structure; see theme 4).
- [ ] ⏳ K562 genome-wide and RPE1 compendia (rerun after a missing-value fix; jobs 0871, 0872).
- [x] ✅ BEELINE: CLR is the strongest of the methods tested (median AUPRC ratio 1.58).
- [x] ✅ Normal-score correlation (CLR-NS) beats MI-based CLR on single-cell data in every gene context (+0.07 to +0.12). QC
- [x] ❌ Attention adds nothing over CLR-NS on BEELINE. QC
- [x] ❌ Wider gene contexts on BEELINE (5,000 genes, all expressed, driven genes): no gain; one dataset drives a nominal loss. QC
- [x] ✅ Driven-gene filter (multimodal or PC-driven) beats variance ranking, which drops well-expressed genes.
- [x] ❌ Random tie-breaking for normal scores: too noisy on small truths; average ranks preferred.
- [x] ❌ Pearson residuals: no better than log-normalization for single-cell co-expression.

### 7. Perturbation diversity

- [x] ⚠️ Single-condition K562 control cells scored about as well as the essential compendium on literature truths: cuts against the simple version of the claim.
- [ ] ☐ Learning curves: AUPR against the number and the diversity (labs, condition types, perturbation classes) of experiments at fixed sample count (M3D, EcoMAC, K562 genome-wide).

### 8. Dependence measure and batch effects (supplement)

- [x] ✅ EcoMAC platform mixtures ("Christmas trees"): 8.4% of gene pairs carry far more MI than their correlation implies, against 0.2-0.3% in single-platform data.
- [x] ⚠️ MI against correlation on uncorrected multi-platform data: MI wins on coherence only, not AUPR.
- [x] ✅ ComBat or centring by platform helps both measures more than switching measure; MI keeps a coherence edge after correction.
- [x] ✅ One clean platform (Genome 2.0, 735 arrays) beats all 2,198 arrays pooled: quality over quantity.
- [x] ✅ MI and correlation interchangeable on single-platform M3D.

### 9. Scale and engineering

- [x] ✅ Packed, band-parallel MI → CLR (from the original `miSubMarix`): 5x lower peak memory at 8,563 genes, bit-identical in double precision.
- [x] ✅ Attention walks restricted to scored genes; memory-lean BEELINE pipeline; BLAS threads per worker.
- [ ] ☐ Use the packed path for genome-wide runs and report memory and time.

### 10. Framing

- [x] ✅ CLR as the first attention operator widely used in computational biology; standard de novo, unsupervised baseline (DREAM5 context; prior-informed and regression methods distinguished).
- [ ] ☐ Review by Jim Collins.

---

## Paper 2: differential CLR for perturbation targets

- [x] ✅ Variance confound identified and controlled (variance-stratified AUROC, variance computed without the perturbation's experiments).
- [x] ⚠️ Attention-based influence against differential expression: 0.572, 12/15, p = 0.015 under the random null; 0.550, 9/15, p = 0.095 under the laboratory-matched null (conservative: it removes related knockouts). QC
- [ ] ⏳ Laboratory-matched null excluding perturbation experiments (job 089).
- [ ] ⏳ SSEM-Lasso re-implementation (job 095).
- [ ] ☐ Port the differential-CLR family into one module: delete-k influence, z-variance (Hayete 2006, `rankClrCluster.m`), Blom-calibrated differential correlation and Klotz scale tests (DiffClr 2022).
- [ ] ☐ Perturb-seq single-cell contrasts (knockdown cells against controls), thousands of knockdowns; TF knockdowns scored against DoRothEA/CollecTRI targets.
- [ ] ☐ Baselines: differential expression, SSEM-Lasso, foundation-model attention.

## Decisions for the authors

- [ ] Confirm the two-paper split and venues.
- [ ] Foundation-model attention comparisons: worth the engineering, and for which paper?
- [ ] BEELINE single-cell results: supplementary section of Paper 1 (recommended) or Paper 2?
