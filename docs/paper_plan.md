# Paper plan

*Living document, started 2026-10-16. Axes of advance, the evidence for each, gaps, and the analyses that close them. Status marks: ✅ established and QC-checked, 🟡 partial or mixed, ⏳ running or queued, ❌ not started. Detailed numbers are in `docs/whitepapers/` (manuscript, supplement, working notes) and `results/`.*

## 1. Two papers

The work splits naturally into two papers with different questions, data and baselines.

| | Paper 1: completed CLR attention | Paper 2: differential CLR for perturbation targets |
|---|---|---|
| Question | Does completing CLR as attention (softmax, depth, held-out depth) recover regulatory programs better, why, and at what resolution? | Which genes does a perturbation rewire? Targets as genes whose CLR context changes between conditions. |
| Core method | Calibrated CLR scores -> softmax attention -> iterated attention; final variación: a multiscale depth sweep that recovers overlapping modules from fine to coarse until a label-free quality collapses. | Differential CLR: delete-k influence (current), change in CLR z-variance (Hayete 2006, `rankClrCluster.m`), Blom-calibrated differential correlation and Klotz row/column scale tests (DiffClr, 2022). |
| Data | M3D and EcoMAC (E. coli); human Perturb-seq pseudobulk compendia (K562 essential, K562 genome-wide, RPE1); BEELINE as a secondary single-cell test. | M3D perturbation groups (15); Perturb-seq K562 genome-wide and RPE1 at single-cell level (knockdown cells against control cells), thousands of knockdowns. |
| Truth | RegulonDB; DoRothEA, CollecTRI, CORUM, knockdown-response modules. | Regulons of the perturbed regulator; knockdown-response sets held out from network building. |
| Baselines | CLR, correlation, MI, GENIE3/GRNBoost2, PIDC, Markov stability; learned attention (scGPT/Geneformer) if feasible. | Differential expression, gene variance, SSEM-Lasso (Cosgrove et al. 2008), foundation-model attention (Kendiukhov 2026 setting). |
| Target venue | Nature Machine Intelligence if the multiscale variación and competitive baselines land; otherwise Nature Communications or Genome Biology. | Genome Biology, Nature Methods or Cell Systems. |

Moving perturbation targets to Paper 2 also removes Paper 1's weakest claim (it does not survive the laboratory-matched null cleanly) and gives that question the data it needs (single-cell knockdowns with direct contrasts).

## 2. Paper 1: axes of advance

| Axis | Claim | Evidence so far | Status | Gap / next step |
|---|---|---|---|---|
| A. Improvement over CLR | Completing CLR as attention improves regulon recovery. | M3D: co-membership AUPR +0.022 [+0.019, +0.026], attention ahead in 1,000 of 1,000 nested-bootstrap draws; coherence +0.026 [+0.008, +0.045]. EcoMAC: +0.032 AUPR (+27%). | ✅ (bacteria) | Human compendia (axis F). |
| B. Mechanism | The gain needs CLR's per-gene context normalization, not a particular dependence measure. | Kernel ablation at held-out depths: softmax on uncalibrated scores gets about a third of the AUPR gain and no coherence gain. Shapley over six components: depth +0.018 and CLR +0.013 (AUPR); CLR +0.025 (coherence). | ✅ | — |
| C. Robustness | The gain is not an artefact. | 50 lab splits (t* = 9, flat optimum); caps 50/100/200; nested bootstrap; survives removal of global components (Supplementary Note 7); top components are growth, motility and σS programs (Note 8). | ✅ | Blocked tests where combinations share datasets. |
| D. Modularity and resolution (final variación) | Depth is a resolution parameter: shallow attention resolves complexes and specific regulons, deep attention union regulons and pathways; a depth sweep recovers overlapping modules at every scale. | Hints: operons gain sharply at t = 2; coherence peaks at t = 4-8 while pooled AUPR keeps rising to about t = 27; attention helps TF regulons and hurts complexes in K562. | ⏳ | Granularity sweep (0873, today). If confirmed: build the multiscale algorithm (modules tracked across t, overlaps, birth/death scales, label-free collapse criterion) on `find_modules()` / `ModuleTester`; compare with Markov stability. |
| E. Performance against competitors | Training-free calibrated attention is competitive with the best de novo methods. | None yet. | ❌ | GENIE3/GRNBoost2, PIDC (and ARACNE as the MI classic) on M3D, EcoMAC and the K562 compendium, same truths. Foundation-model attention on a shared task if feasible. |
| F. Bacteria and human | The result transfers to human regulation. | K562 essential compendium: attention beats CLR on DoRothEA and CollecTRI (+0.008/+0.009 AUPR, significant) but loses on CORUM complexes. K562 genome-wide and RPE1 rerunning (MI methods failed on missing values). | 🟡 ⏳ | 0871/0872 today; QC. Frame complexes against regulons through axis D. |
| G. Perturbation diversity | Inference power comes from the diversity of network states sampled, not from sample count. | Single-condition K562 control cells (10,691 cells) scored about as well as the essential-gene compendium on literature truths, which cuts against the simple version of the claim. BEELINE: wider gene contexts add nothing on single-cell data. | 🟡 | Proper test: learning curves of AUPR against the number and diversity (labs, condition types, perturbation classes) of experiments, at fixed sample count, on M3D, EcoMAC and K562 genome-wide. |
| H. Dependence measure | MI matters where data are heterogeneous; rank-based correlation where counts are zero-inflated. | EcoMAC: MI registers platform mixtures ("Christmas trees", 8.4% of pairs above the Gaussian MI curve) and keeps a coherence edge after correction. BEELINE: CLR-NS beats MI-CLR in every gene context. | ✅ (supplement) | — |
| I. Scalability | The method scales to genome-wide human compendia. | Packed, band-parallel MI -> CLR in the C++ core (from the original `miSubMarix`): 5x lower peak memory at 8,563 genes, bit-identical in double precision. | ✅ | Use it for the genome-wide runs; report memory and time. |
| J. History and framing | CLR was the first attention operator widely used in computational biology, and the standard de novo baseline. | Manuscript introduction (DREAM5 context, de novo vs prior-informed and regression-based methods). | ✅ | Jim Collins's review. |

### Paper 1 figure plan

1. CLR as attention: the components and their transformer counterparts.
2. Main result and depth (M3D, EcoMAC), with held-out depth distribution.
3. Shapley decomposition (done: `fig_shapley_SC`).
4. Resolution: best depth rises from fine to coarse truths; the multiscale module atlas (axis D).
5. Competitors on bacterial and human compendia (axis E).
6. Human transfer and perturbation diversity (axes F, G).

## 3. Paper 2: differential CLR for perturbation targets

**What exists.**

- Delete-k influence scores on M3D with random and laboratory-matched nulls. Attention-based influence beats differential expression after variance stratification: 0.572, 12/15, p = 0.015 under the random null. Under the laboratory-matched null this drops to 0.550, p = 0.095, and that null is conservative. The null excluding perturbation experiments (089) is running, and SSEM-Lasso is queued.
- Your 2006 z-variance method (`rankClrCluster.m`, `zMIvar.m` in `data/CLR versions/CLR.zip`): ranks genes by the change in their CLR z-score variance with and without the test experiments.
- DiffClr (2022, `~/iCloud/Documents/History/Prototypes/DiffClr`):
  - `dCLR.R`: Blom-calibrated row and column z-scores of the correlation difference between conditions, Stouffer-combined, with BH p-values and module discovery (Horn, clustering);
  - `diffClr.R`: per-pair Klotz scale tests on row and column correlation profiles, Stouffer-combined;
  - `glmBased.R`: interaction-model alternative.

**What it needs.**

1. **Port and unify** influence, z-variance, dCLR and Klotz diffCLR as one R module with a shared interface (condition A cells or experiments against condition B), using the packed CLR path for scale.
2. **Data with direct contrasts:** Perturb-seq K562 genome-wide and RPE1 single cells, knockdown cells against controls. For TF knockdowns, the truth is the TF's DoRothEA or CollecTRI targets. For all knockdowns, it's held-out response sets.
3. **Baselines:** differential expression (the Kendiukhov 2026 trivial baseline), SSEM-Lasso, and foundation-model attention.
4. **Evaluation:** variance-stratified AUROC (as now), top-100 sensitivity, laboratory- or batch-matched nulls, and scaling to thousands of knockdowns.

## 4. Priorities (next two weeks)

1. Finish today's queue: lab null without perturbation experiments (089), K562 genome-wide and RPE1 compendia, granularity sweep. QC each.
2. Competitor baselines (axis E): GENIE3/GRNBoost2, PIDC, ARACNE on M3D, EcoMAC, K562.
3. If the granularity sweep supports it: the multiscale algorithm (axis D), the main bet for Nature Machine Intelligence.
4. Perturbation-diversity learning curves (axis G).
5. Paper 2 groundwork: port the differential-CLR family and set up the single-cell knockdown contrasts.

## 5. Decisions for the authors

- Confirm the two-paper split and the venue strategy.
- Whether foundation-model attention comparisons (scGPT/Geneformer) are worth the engineering for Paper 1, Paper 2, or both.
- Whether the single-cell BEELINE results stay as a supplementary section of Paper 1 (current recommendation) or move to Paper 2.
