# CLR was an attention mechanism: completing it with softmax, aggregation and depth improves regulon recovery

*Working notes behind the manuscript (`paper_nmi.md`, Nature Machine Intelligence format): full findings tables, the running-analysis queue and the journal assessment. Quick draft, 2026-10-14.  It assembles the findings logs (`methods.md`, `ecoli_m3d.md`, `beeline.md`) and the runs finished on crayon up to 2026-10-14 08:30.  Every number is provisional until the QC pass; numbers marked ⏳ will change when the queued full-scale analyses finish (section 8).*

---

## Abstract

Context likelihood of relatedness (CLR; Faith, Hayete et al. 2007) infers gene regulatory networks by scoring each gene's mutual information with every other gene against that gene's own background.  Read as a machine-learning operator, those scores are attention scores:  each gene (the query) gets a calibrated, data-dependent weight over all other genes (the keys), normalized within its own context.  CLR stopped at the scores.  It never used the rest of the attention machinery that neural-network research later made standard:  weight normalization with a temperature, aggregation, and depth.

We complete the operator.  Calibrated CLR scores become softmax attention with a per-gene temperature, and attention is applied repeatedly, to a depth chosen on held-out laboratories without reference to any gold standard.

On the *E. coli* M3D compendium (466 experiments) scored against RegulonDB, the completed operator gives these gains over CLR on regulon co-membership:

- **Pooled AUPR:** 0.157 against 0.136, +15% relative, +0.021 [+0.015, +0.027].
- **Median per-regulon coherence:** 0.632 against 0.597, +0.031 [+0.018, +0.046].
- **Replicate-aware cluster bootstrap:** across 1,000 draws that resample whole experiments and recompute everything from scratch, the gain holds (+0.020 [+0.018, +0.023]).

Three controls show where the gain comes from:

- The gain requires CLR's context normalization, not MI as such.  The same propagation without CLR calibration loses most of the gain (section 3.6), while CLR on correlation works as well as CLR on MI in bulk data.
- A hard-thresholded operator gains less, so the softmax matters.
- Coherence peaks at intermediate depth, so depth trades sharpness against reach.

The same operator ranks the targets of perturbed regulators better than differential expression once gene variance is controlled for:  median AUROC 0.572 against 0.516, better in 12 of 15 perturbations (p = 0.015).  On single-cell data (BEELINE), CLR is the strongest of the methods we tested, but attention does not improve on it in BEELINE's narrow gene sets ⏳.

CLR's estimated, calibrated attention needs no training corpus or embeddings, and it carries a permutation null.  Its 2007 form is an early, independently derived instance of the attention-score pattern, and completing it is worthwhile.

---

## 1. Introduction

**The bridge.**  Attention computes, for each query, normalized compatibilities with a set of keys, then aggregates values with those weights.  CLR computes, for each gene, a dependence (MI) with every other gene, and normalizes it per query against that gene's own MI distribution.  That row-wise z-score plays the role softmax plays in attention:  it makes one query's weights comparable within its own context.  The compatibility is an *estimated* statistical dependence, not a learned dot product.  It is sign-blind, it handles non-monotone relationships, and it needs no embeddings, which matters when gene "meanings" depend on tissue.

**Priority, stated carefully.**  We do not claim CLR was the first attention mechanism.  Kernel smoothing (Nadaraya 1964; Watson 1964) and fast-weight programmers (Schmidhuber 1991; Schlag et al. 2021) are earlier instances of the pattern.  The defensible claim is narrower:  CLR (2007) independently computed attention-style scores over continuous molecular measurements, seven years before neural attention (Bahdanau et al. 2014), and stopped at the scores.

**The gap.**  CLR scores were only ever used to rank edges.  Three components of attention were missing:

- normalization with a temperature, which sets how many genes each gene attends to;
- aggregation, which reads out multi-hop attention mass or attention-weighted profiles;
- depth, meaning stacked or iterated attention.

**Relation to network propagation.**  Iterated attention through a lazy operator is a random walk on a graph.  It is mathematically related to diffusion kernels (Kondor & Lafferty 2002) and to label propagation (Zhou & Schölkopf 2004), and network propagation is a known "universal amplifier" in genomics (Cowen et al. 2017).  We accept this.  The contribution is not propagation itself.  It is:

1. the kernel:  context-normalized (CLR) dependence scores.  The same softmax propagation on uncalibrated scores loses most of the gain.  Whether the dependence is MI or correlation matters little in bulk data (section 3.6);
2. the attention design choices:  the softmax temperature, set by an effective neighbour count;
3. a label-free, held-out procedure for choosing depth.

**Contrast with learned attention.**  Transformer-based GRN methods and single-cell foundation models (scGPT, Geneformer) learn attention from large corpora.  A recent systematic evaluation found that their attention-derived edge scores add no value over trivial gene-level features for predicting perturbation targets (Kendiukhov 2026).  We test the same question with estimated attention, on the compendium where the perturbation-target task was first benchmarked (Cosgrove et al. 2008).

---

## 2. Method

1. **MI.**  B-spline MI (Daub et al. 2004), spline order 3, with one common bin count chosen for the joint histogram:  the Hacine-Gharbi et al. (2012) low-bias rule at ρ = 0, which gives 7 bins at N = 466.  Per-gene adaptive bins make heavy-tailed genes false hubs.  MI is computed on raw values.
2. **Calibration (CLR).**  Each row of the MI matrix is z-scored against its own background (negative z clipped to 0), and the two directions are combined by Stouffer's rule, (z_ij + z_ji)/√2.  The 2007 Euclidean rule gives an almost identical ranking (Spearman 0.993).
3. **Permutation null.**  Samples are shuffled independently per gene, and MI and CLR are recomputed.  Edges are selected by Benjamini–Hochberg on CLR scores.  This asks whether an edge is exceptional relative to each gene's own background, and it stays sparse under global confounders where an MI null does not.
4. **Attention operator (`soft10`).**  Each row is a softmax over the gene's top 50 CLR scores, with the temperature set per row by bisection so that the effective number of attended genes, exp(entropy), is 10.
5. **Depth.**  P = (1 − α)I + αA, with α = 0.5.  The primary readout is attention mass, (Pᵗ + Pᵗᵀ)/2.
6. **Held-out depth selection.**  Experiments are split by laboratory, and the operator is built on one half.  t* maximizes how well attention mass predicts the other half's top 1% of CLR pairs, so RegulonDB is never used.  For single-cell data, cells are split into random halves.

**Controls, fixed in advance.**  These separate the kernel, the attention design and the propagation:

- `pearson_top10`:  the identical propagation on |Pearson| top-10 neighbours.
- `fdr05_top10`:  a hard BH-selected CLR graph plus each gene's own top-10 neighbours.
- |Pearson|, raw MI, and the parity-2007 CLR.

---

## 3. Regulon recovery in *E. coli*

**Data.**  The M3D compendium (4,297 genes; 907 Affymetrix chips from 466 experiments) is the compendium CLR was originally validated on.  The gold standard is RegulonDB 14.5, with every regulator class expanded to its genes (TFs, sRNAs, sigma factors, ppGpp, other proteins).  Regulons of 5–500 genes are kept.

**Benchmark.**  Regulator-agnostic regulon co-membership.  A gene pair is positive if both genes share a regulon, and same-operon pairs are excluded.  Regulators are not used as nodes, because a regulator's activity is generally latent in its mRNA.  At strong/confirmed evidence there are 155 regulons, 1,879 genes and 1.76 M pairs, with a base rate of 9.4%.

### 3.1 Main result (466 experiments, strong/confirmed evidence)

| Method | Co-membership AUPR | Δ vs CLR [95% CI] | Median coherence | Δ vs CLR [95% CI] |
|---|---|---|---|---|
| \|Pearson\| | 0.116 | −0.020 [−0.024, −0.016] | | |
| Raw MI | 0.115 | −0.021 [−0.026, −0.016] | | |
| CLR, parity 2007 | 0.134 | −0.002 [−0.003, −0.001] | | |
| **CLR (HG, Stouffer)** | **0.136** | 0 | 0.597 | 0 |
| `pearson_top10`, t* = 6 | 0.135 | −0.001 [−0.005, +0.003] | | −0.014 [−0.025, −0.003] |
| `fdr05_top10`, t* = 7 | 0.148 | +0.013 | | ≈ 0 |
| **`soft10`, t* = 8** | **0.157** | **+0.021 [+0.015, +0.027]** | **0.632** | **+0.031 [+0.018, +0.046]** |

Intervals come from a paired delete-half jackknife over genes (AUPR) and a paired bootstrap over regulons (coherence).  For `soft10`, 69% of 144 regulons improve (Wilcoxon p = 4×10⁻⁴).  At all evidence levels the gains are +0.007 AUPR and +0.026 coherence, both significant.

**Held-out depth.**  The two lab-split directions chose t = 10 and 7, a geometric mean of 8.  Pooled AUPR keeps rising until about t = 32, but coherence peaks at t = 4–8 and collapses by t = 128.  The label-free criterion lands where both metrics improve.  ⏳ With 50 random lab splits, the full criterion curve at t = 1…64 will show how stable t* is and how flat the criterion is around it (`lab_splits.R`).

### 3.2 Robustness:  replicate-aware bootstrap on 907 chips

Every draw recomputes MI, CLR and attention from scratch.  The cluster draws (n = 1,000) resample experiments with replacement and use one replicate chip per occurrence.  The replicate draws (n = 200) take one chip per experiment.

| Method | ΔAUPR vs CLR, cluster | ΔCoherence vs CLR, cluster | ΔAUPR vs CLR, replicate |
|---|---|---|---|
| `soft10`, t = 4 | +0.014 [+0.012, +0.016] | +0.032 [+0.015, +0.050] | +0.014 [+0.014, +0.015] |
| `soft10`, t = 8 | **+0.020 [+0.018, +0.023]** | **+0.028 [+0.010, +0.046]** | +0.021 [+0.020, +0.022] |
| `soft10`, t = 16 | +0.027 [+0.025, +0.030] | +0.019 [−0.000, +0.038] | +0.028 [+0.027, +0.029] |
| `pearson_top10`, t = 6 | −0.001 [−0.004, +0.001] | −0.026 [−0.046, −0.005] | −0.002 [−0.003, −0.001] |
| Raw MI | −0.019 [−0.022, −0.017] | −0.009 [−0.022, +0.004] | −0.020 |
| \|Pearson\| | −0.019 [−0.021, −0.017] | −0.022 [−0.037, −0.007] | −0.020 |

(Strong/confirmed evidence; percentile intervals over 2.5–97.5% of draws.  At all evidence levels, `soft10` at t = 8 gains +0.007 [+0.005, +0.009] AUPR.)

**Reading.**

- The gain is not an artefact of pseudo-replication or of a few influential experiments.
- It needs CLR calibration:  hard top-10 propagation on raw correlation is flat on AUPR and costs coherence (section 3.6 separates calibration from softmax).
- ⏳ These intervals hold t fixed.  The nested bootstrap, which re-selects t inside each of 1,000 draws, will give intervals for "attention at the depth the procedure would choose".

### 3.3 Sensitivity analyses that did not change the conclusion

- **Permutation null size.**  Moving from B = 100 to B = 1,000 leaves the 5% FDR edge set unchanged (4,800 edges) and moves the 20% set by 0.2% (8,590 to 8,608 edges).  Every downstream metric is identical.
- **Binning.**  Coarse common bins (6–10) are indistinguishable from one another, and all beat per-gene Freedman–Diaconis bins.  CLR AUPR falls monotonically from 8 to 16 bins.
- **Combination rule.**  Stouffer and Euclidean combination are equivalent.
- **Transform.**  A rank (copula) transform helps on synthetic monotone distortions but is slightly worse on real data (0.127 against 0.134).
- ⏳ **Softmax cap.**  The cap of 50 candidates per row will be tested against caps of 100 and 200, and against no cap.

### 3.4 Negative results (reported deliberately)

- **Value transforms.**  A signed value and a conditional-expectation W_V, read off the MI estimator's own B-spline basis, restore repressed targets on synthetic data (+0.78 → −0.97).  They do not help on the real benchmark, because 99% of the selected edges are positively correlated.
- **Diffused-profile readout.**  The readout |cor(PᵗZ)| adds little; the best case is +0.007.
- **FDR operator alone.**  The permutation-FDR operator leaves 62–77% of genes isolated, which makes it unusable for propagation.
- **Regulator classes.**  sRNA regulons are nearly unrecoverable from mRNA co-expression under any method.  The diffuse ppGpp regulon *loses* under propagation.  Most of the gain comes from TF regulons:  AUPR 0.114 → 0.145, coherence 0.59 → 0.63.

### 3.5 Direct TF → target edges:  an edge detector against a target-set ranker

Under the 2007 convention, candidate edges are (TF gene, gene) pairs scored by the TF's own mRNA.  At strong/confirmed evidence there are 199 TFs and 3,971 positives, a base rate of 0.48%.

| Method | Pooled AUPR (× random) | Median per-TF AUROC | TFs beating CLR |
|---|---|---|---|
| CLR | **6.75** | 0.595 | — |
| `soft10`, t = 4 | 4.65 | **0.642** | 91/149 |
| `soft10`, t = 8 | 5.15 | 0.623 | 86/149 |
| `pearson_top10`, t = 6 | 4.07 | 0.535 | 61/149 |

Single-step CLR is the better *edge detector* at the top of the ranking.  Attention is the better *target-set ranker* across each TF's whole row.  That matches its gains on co-membership and coherence.  ⏳ The nested bootstrap adds intervals to both metrics.

### 3.6 What the gain needs:  context normalization, not MI (preliminary)

**Ablation on the 466 experiments.**  This is a point estimate at t = 8 on strong/confirmed evidence.  ⏳ Job 065 repeats it with each kernel's own held-out t* over 20 lab splits, with jackknife intervals.

| Operator | Co-membership AUPR | Median coherence |
|---|---|---|
| CLR on B-spline MI (no attention) | 0.136 | 0.597 |
| CLR on \|r\| (no attention) | 0.138 | 0.590 |
| CLR on \|r\| of normal scores (no attention) | 0.133 | 0.591 |
| Softmax attention on raw \|r\|, **no CLR** | 0.143 | 0.586 |
| Softmax attention on raw MI, **no CLR** | 0.143 | 0.583 |
| Softmax attention on CLR(\|r\|) | 0.155 | 0.622 |
| Softmax attention on CLR(normal-score \|r\|) | 0.157 | 0.626 |
| **Softmax attention on CLR(MI)** (`soft10`) | **0.157** | **0.632** |

**Reading.**

- **Propagation alone** buys about a third of the AUPR gain, and it *lowers* coherence.
- **The rest of the AUPR gain, and all of the coherence gain,** appears only when the softmax acts on context-normalized scores.
- **MI against correlation:**  the underlying dependence measure barely matters in this compendium, where 99% of selected edges are positively correlated.

This sharpens the paper's thesis:  CLR's per-gene context normalization is the attention-score step, and attention needs it.  It also weakens any claim that the MI estimator itself is essential.

---

## 4. Perturbation targets beyond differential expression

**Task.**  The task follows Cosgrove et al. (2008), on the same compendium:  rank the genes of a perturbed regulator's regulon.  There are 15 perturbation groups, mapped to RegulonDB regulators by a table fixed in advance.  The perturbed genes themselves are excluded from scoring.

**Score.**  The delete-k jackknife influence of the perturbation's k experiments on each network:  how much each gene's network row changes when those experiments are removed.  It is calibrated per gene against B = 500 random delete-k subsets.  No Gaussian assumption is made.

**The confound, and how we control it.**  Gene variance alone, which knows nothing about the perturbation, reaches a median AUROC of 0.646 and beats differential expression in 14 of 15 perturbations, simply because regulon members are more variable than other genes.  The primary metric therefore compares targets only with non-targets in the same variance decile, with variance computed without the perturbation's own experiments.

| Score | Raw median AUROC | Variance-stratified median AUROC | Wins vs DE (stratified) | p vs DE |
|---|---|---|---|---|
| **CLR-attention influence** | 0.583 | **0.572** | **12/15** | **0.015** |
| DE + CLR-attention (Stouffer) | 0.565 | 0.548 | 12/15 | 0.030 |
| CLR influence | 0.564 | 0.543 | 11/15 | 0.15 |
| \|Pearson\| influence | 0.534 | 0.530 | 8/15 | 0.80 |
| Gene variance | 0.646 | 0.521 | 8/15 | 0.76 |
| Differential expression | 0.519 | 0.516 | — | — |

**Reading.**  Once variance is held fixed, the variance prior's advantage disappears.  CLR-attention influence keeps a modest, significant advantage over differential expression.

- **Where it helps:**  the advantage is concentrated in HU, AppY, LexA/SOS and RyhB.
- **Where it fails:**  the diffuse global regulators ppGpp, CRP and FNR.
- **Pearson influence:**  adds nothing, so again the signal needs the CLR kernel.
- **Top-100 sensitivity:**  low for every method (at most 0.075).

⏳ Two analyses are pending:

- **A lab-matched null:** the random subsets are drawn from the perturbation's own labs, which asks whether *this* perturbation matters rather than any k experiments from that lab.
- **A like-for-like comparison with Cosgrove et al. (2008):** a re-implementation of their SSEM-Lasso network filter.

---

## 5. Single-cell data (BEELINE)

**Protocol.**  We follow BEELINE (Pratapa et al. 2020):

- **Datasets:**  the seven experimental scRNA-seq datasets.
- **Gene sets:**  the "TFs + 500" and "TFs + 1000" selections.
- **Metrics:**  AUPRC ratio and early precision ratio.
- **Ground truths:**  cell-type and non-specific ChIP-seq, STRING, and mESC LOF/GOF, plus two refinements we built:
  - evidence-filtered STRING (`regulatory`:  v11 transcriptional actions; `curated`:  experimental or database channels);
  - graded ChIP-Atlas binding (per-TF Spearman correlation with MACS2 peak scores).

**Depth.**  The held-out depth, chosen on random halves of cells, is t* = 5–11 for CLR attention in all 14 cases.  That brackets the *E. coli* t* = 8, although depth was never tuned on single-cell data.

**Median AUPRC ratio (7 datasets × 2 gene sets; ground truths as listed):**

| Truth (cases) | CLR | CLR attention | \|Pearson\| | Raw MI | Partial correlation (ridge) |
|---|---|---|---|---|---|
| Cell-type ChIP (14) | 1.00 | 0.93 | 1.02 | 1.01 | 1.01 |
| Non-specific ChIP (14) | 1.54 | **1.59** | 1.32 | 1.24 | 1.08 |
| STRING (14) | **2.11** | 1.76 | 1.90 | 1.77 | 1.19 |
| STRING, curated (14) | **1.40** | 1.37 | 1.35 | 1.21 | 1.18 |
| STRING, regulatory (14) | 1.91 | **2.19** | 1.64 | 1.56 | 1.02 |
| mESC LOF/GOF (2) | 1.45 | **1.54** | 1.34 | 1.41 | 0.93 |
| **All 72** | **1.58** | 1.54 | 1.32 | 1.20 | 1.08 |

**Reading.**

1. **CLR is the strongest method here,** about 20–45% above correlation, MI and the linear-model baseline, so the context-normalization idea carries over to single-cell data.
2. **Attention neither clearly helps nor clearly hurts.**
   - It is better on non-specific ChIP, on STRING-regulatory and on LOF/GOF.
   - It is worse on full STRING, whose co-expression channel rewards direct correlation.
   - It wins against CLR in 28% of cases.
3. **The cell-type ChIP truth and graded binding are at chance for every method.**  The median per-TF ρ with binding strength is between −0.04 and +0.02.  Binding is not regulation, and these truths cannot separate the methods.

**Zero inflation and the dependence measure (small test, `analysis/beeline_blom_test.R`).**  Single-cell matrices are zero-inflated; in BEELINE's selected genes 11–38% of entries are zero, and far more in wider gene sets.  We gave each gene normal scores (Blom, (rank − 3/8)/(n + 1/4)) and recalibrated.  Results over 72 case × truth combinations:

| Variant | ΔAUPRC ratio vs current CLR | Wins | p |
|---|---|---|---|
| CLR on \|r\| of normal scores, average ranks for ties ("CLR-NS") | **+0.115** | 57/72 | 1×10⁻⁶ |
| Same, ties broken at random (mean of 10 draws) | +0.071 | 50/72 | 0.002 |
| Same, ties broken at random (one draw) | +0.053 | 49/72 | 0.005 |
| CLR on B-spline MI of normal scores | +0.021 | 46/72 | 0.42 |
| CLR on Gaussian MI, −½ log(1 − r²), of normal scores | +0.013 | 37/72 | 0.40 |

- **On single-cell data, CLR does not need MI.**  CLR on correlation of normal scores is the best CLR we have tested.
- **Random tie-breaking is noisy** for truths with few edges, so average ranks are preferred.
- **Attention on CLR-NS** helps on regulatory STRING (2.17 against 1.74) and ties on non-specific ChIP.  It loses on full STRING, whose co-expression channel rewards direct correlation.
- **Status:**  these are point estimates, at depths chosen for raw CLR.  The queued runs add CLR-NS with its own held-out depth, in wider gene contexts and with bootstrap intervals.

**Caveat this draft must resolve ⏳.**  BEELINE builds every network on 500–1,500 selected genes, 30–60% of them TFs.  CLR's calibration is background-dependent:  each gene's MI is judged against its distribution over *all other genes*.  A thin, TF-heavy background handicaps CLR, and with it the attention built on CLR.  The queued runs rebuild CLR on the 5,000 most variable genes and on all expressed genes (up to about 18,000), still scoring only BEELINE's gene sets.  A nested cell bootstrap then gives intervals.  The single-cell conclusion should wait for these runs.

---

## 6. Discussion

**What attention adds to CLR.**  Attention adds reach without losing calibration.  Multi-hop attention mass recovers co-regulated genes that share few direct strong dependences.  The softmax temperature keeps each gene's attention focused on about 10 effective neighbours, which hard thresholds cannot do without isolating most genes.  The gains are modest in absolute AUPR (+0.02, about 15% relative), but they are consistent across:

- resampling schemes:  the gene jackknife, the regulon bootstrap, and the experiment and replicate bootstraps;
- evidence levels:  strong/confirmed and all evidence;
- tasks:  regulon co-membership, coherence, per-TF target ranking, and perturbation targets.

**What it does not do.**  Attention does not detect the strongest direct edges better than CLR.  It dilutes them.  It also fails on diffuse global programs (ppGpp, CRP).  The honest summary is that single-step CLR is the edge detector, and attention is the module and target-set ranker.

**Why estimated attention matters now.**  Learned attention from foundation models does not beat trivial features on perturbation targets (Kendiukhov 2026).  Estimated attention, with an explicit null and a held-out depth, does, though modestly.  This argues for treating attention as an *estimator* with statistical semantics, not only as a learned layer.

**Limitations.**

- RegulonDB evidence partly comes from expression studies.  The strong/confirmed tier limits this, and it applies equally to every method compared.
- The depth criterion predicts CLR structure in held-out labs, which is a proxy for biological accuracy, not the same thing.
- Experiments within a lab are correlated:  the lab split and the cluster bootstrap handle this, but the gene jackknife does not.
- One bacterial compendium carries the main result.  A second organism with a regulon-level gold standard would strengthen it.  Candidates are *B. subtilis* (SubtiWiki) and yeast (YEASTRACT with a large bulk compendium).

---

## 7. Figures (plan)

1. **The bridge:**  CLR as attention scores, and the three missing components.
2. **Main result:**  AUPR and coherence against depth for `soft10`, `pearson_top10`, `fdr05_top10` and CLR, with the held-out t* distribution overlaid (from the 50 lab splits).
3. **Bootstrap forest plot:**  Δ vs CLR with intervals across the jackknife, regulon, cluster, replicate and nested bootstraps.
4. **Edge detector against target-set ranker:**  pooled TF-node AUPR and per-TF AUROC.
5. **Perturbation targets:**  raw against variance-stratified AUROC per perturbation, with the random and lab-matched nulls and SSEM-Lasso.
6. **BEELINE:**  AUPRC ratio by truth for each gene context (selected genes, 5,000 genes, all genes).

---

## 8. Analyses running (2026-10-14, crayon) and what each settles

| # | Analysis | Settles |
|---|---|---|
| 060 | 50 random lab splits; criterion curve at t = 1…64; softmax caps 50, 100, 200 and none | Is t* stable and the optimum flat?  Is the cap of 50 driving the result? |
| 070 | Nested chips bootstrap (1,000 draws), with t re-selected on each draw's own lab split, plus TF-node metrics | Intervals that include depth-selection uncertainty; intervals on the edge-detector/ranker trade-off |
| 080 | Lab-matched perturbation null (500 draws per perturbation) | Is the perturbation signal perturbation-specific, or a lab effect? |
| 085, 087 | BEELINE with 5,000-gene and all-gene CLR context | Does a representative background restore CLR, and attention on top of it, in single-cell data? |
| 090 | BEELINE nested cell bootstrap (200 draws) and 20 cell halvings, 5,000-gene context | Intervals for the single-cell comparison; stability of the single-cell t* |
| 095 | SSEM-Lasso re-implementation on the same perturbations | A like-for-like comparison with the 2008 state of the art |
| 100 | FDR operators inside the chips bootstrap (200 draws, per-draw null with B = 100) | Intervals for the hard-threshold control |

---

## 9. Target journal (proposal)

The table below gives impact factors as rough recent values; check them before submitting.

| Journal | IF (approx.) | Fit | What it would take |
|---|---|---|---|
| *Nature Machine Intelligence* | ~20 | The "estimated, calibrated attention" angle and the contrast with learned foundation-model attention suit it.  Its readers are ML-first. | A result that generalizes beyond one bacterial compendium.  That means the single-cell context runs coming out positive, a second organism, and ideally a head-to-head with learned attention (scGPT/Geneformer attention maps) on the same perturbation task.  On current evidence, a +15% relative gain on one compendium is likely to be judged incremental. |
| *Nature Communications* | ~15 | Broad scope.  It accepts methods papers with solid benchmarks and a conceptual hook. | The current results plus the queued robustness analyses, with a clear single-cell story either way. |
| *Genome Biology* | ~10–12 | A natural home for GRN inference methods and benchmarks; reviewers know CLR and BEELINE. | The current results, the queued analyses, and the published BEELINE baselines tabulated (GENIE3, GRNBoost2, PIDC, PPCOR). |
| *Cell Systems* | ~9 | Likes conceptual bridges and systems-level framing. | Similar to *Genome Biology*, with more emphasis on the biological reading (edge detector against module ranker, the regulator classes). |
| *Nature Computational Science* | ~12 | Methods with cross-disciplinary framing. | Similar to *Nature Communications*. |

**Recommendation:**

- **Aim for *Nature Communications* or *Genome Biology*.**
- **Keep *Nature Machine Intelligence* open,** conditional on two results:
  1. the gene-context BEELINE runs showing that attention on a representative CLR background helps on single-cell data;
  2. a direct comparison in which estimated CLR attention beats learned foundation-model attention on perturbation targets.

  If both land, the paper becomes a general statement about attention as a calibrated estimator, which is an NMI paper.  If not, it is a strong, honest methods paper with a historical hook, which suits *Genome Biology* or *Nature Communications*.
