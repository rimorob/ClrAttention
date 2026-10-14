# Context-normalized attention without learning recovers gene regulatory programs

*Manuscript draft in* Nature Machine Intelligence *Article format, 2026-10-14.  Bracketed items marked [pending] depend on analyses now running; all numbers await an independent QC pass.  Working notes, the analysis queue and the journal assessment are in `paper_draft_notes.md`.*

---

## Abstract

Attention weights each element of a set by its normalized compatibility with every other element.  The idea predates neural networks; kernel regression and fast-weight programmers are earlier instances.  In computational biology, context likelihood of relatedness (CLR), a 2007 method that remains a standard for de novo gene-network inference, was arguably the first attention operator in wide use.  CLR computes attention scores from data alone, z-scoring each gene's dependence on every other gene against that gene's own background, but never used them as attention.  Here we bring back to CLR the advances attention research made in language modelling:  a temperature-controlled softmax, iterated attention, and depth chosen on held-out laboratories.  On *Escherichia coli*, regulon recovery improves by 15% and per-regulon coherence rises, robustly across 1,000 resamplings of experiments.  The gain requires CLR's context normalization, not a particular dependence measure.  The estimated attention also identifies perturbation targets beyond differential expression.

---

## Introduction

Attention is now the dominant way that machine-learning models relate the elements of a set.  For each query, a compatibility with every key is computed, normalized across keys (usually by a softmax with a temperature), and used to aggregate values.  Stacking or iterating the operation lets information flow along multi-step paths.  In transformers, the compatibilities are dot products of learned embeddings, and the network learns what to attend to from large corpora.

The pattern is older than its neural form, and we make no claim that CLR invented it.  Kernel regression (Nadaraya 1964; Watson 1964) weights observations by their normalized similarity to a query and has been a standard statistical tool for decades.  Fast-weight programmers (Schmidhuber 1991) are formally equivalent to linear attention (Schlag et al. 2021).  Our claim is narrower.  In computational biology, attention entered everyday practice through a different door, years before neural attention (Bahdanau et al. 2014; Vaswani et al. 2017).

Context likelihood of relatedness (CLR; Faith, Hayete et al. 2007) infers gene regulatory networks from expression compendia.  For every gene it estimates the mutual information (MI) with every other gene, then expresses each value as a z-score against the distribution of that gene's MI with all genes.  In attention terms:

- each gene is a query;
- every other gene is a key;
- the compatibility is an *estimated* statistical dependence, not a learned dot product;
- the normalization is per query, against the query's own context.

That last step is what made CLR work.  A gene that is weakly dependent on everything is not a hub, and a modest dependence that stands out from a gene's background is a candidate regulatory link.

CLR became a standard tool for *de novo* network inference, that is, inference from expression data alone without prior knowledge.  Methods that add informative priors, such as known pathways, binding motifs or curated interactions, can outperform it, but they answer a different question.  Among methods that start from expression alone, CLR has remained a standard reference for nearly two decades:

- it was one of the "well-established off-the-shelf" methods benchmarked in the DREAM5 network-inference challenge, and a component of several of that challenge's ensemble predictors (Marbach et al. 2012);
- on the single-cell benchmark we use here, it is the strongest of the methods we evaluated.

It is not the top scorer on every benchmark.  In DREAM5, regression-based methods such as tree ensembles scored higher on directed transcription-factor-to-target edges.  Those methods fit a predictive model for every gene from a supplied list of candidate regulators.  CLR is different in kind:  it is model-free, unsupervised estimation of pairwise dependence with per-gene calibration.  It needs neither a regulator list nor a predictive model, so it applies equally when regulators are unknown or latent, which is the setting of our main benchmark.  Its simplicity, speed and stability across data types have kept it the usual unsupervised baseline.  To our knowledge, CLR is the first attention-style operator to be widely adopted in computational biology.

CLR, however, stopped at the scores.  It ranked them as edges and never applied the rest of attention:  a temperature-controlled normalization, aggregation, or depth.  Over the same two decades, attention research in natural-language processing turned those components into a mature toolkit.  The aim of this work is to bring that toolkit back to the setting where CLR has proved durable, and to measure what each component adds.

Here we complete the operator and ask what each component contributes.  Three findings stand out:

- Completing CLR as attention improves recovery of bacterial regulatory programs on the compendium CLR was validated on, robustly to resampling of experiments.
- The improvement requires the per-query context normalization, not the specific dependence measure.  Softmax propagation on uncalibrated scores recovers little of it, while correlation and MI are interchangeable as inputs.
- Estimated attention carries information about perturbation targets beyond differential expression, a task on which learned attention maps from single-cell foundation models add no value over trivial features (Kendiukhov 2026).

Because the operator needs no training data, no embeddings and no labels, and comes with a permutation null, it is a calibrated, interpretable counterpart to learned attention for domains where the elements have no fixed meaning.

---

## Results

### CLR scores are attention scores, and completing them takes three components

Let D be a G × G matrix of pairwise dependences between genes, estimated from an expression compendium (Fig. 1).  CLR computes, for each row i, z_ij = (D_ij − μ_i)/σ_i, clips negative values to zero, and combines the two directions into a symmetric score S_ij = (z_ij + z_ji)/√2.  The row-wise standardization is a per-query normalization, the same role softmax plays in attention.

We add three components:

1. **Softmax with a per-query temperature.**  Each row of S becomes attention weights A_i· = softmax(S_i·/τ_i), over the gene's 50 highest scores.  The temperature τ_i is set so that the effective number of attended genes, exp(entropy), is 10 for every gene.  This makes attention comparably focused for genes with very different score distributions.
2. **Depth.**  Attention is iterated through a lazy operator, P = ½(I + A).  We read out the symmetric attention mass between genes after t steps, ½(Pᵗ + Pᵗᵀ).
3. **Label-free depth selection.**  t is chosen on held-out data.  Experiments are split by laboratory, the operator is built on one half, and t maximizes how well attention mass predicts the strongest CLR pairs of the other half.  No gold standard is used.

Iterated attention through a row-stochastic operator is a random walk.  The construction is therefore related to diffusion kernels (Kondor & Lafferty 2002) and to network propagation (Cowen et al. 2017).  We treat that relationship as a baseline to beat, not a novelty to claim:  every comparison below includes propagation on uncalibrated scores.

### Completed CLR attention improves recovery of *E. coli* regulons

We evaluated the operator on the M3D *E. coli* compendium (4,297 genes; 907 microarrays from 466 experiments; Faith et al. 2008) against RegulonDB 14.5.  The benchmark is regulator-agnostic:  a gene pair is positive if both genes belong to the same regulon, for any regulator class (transcription factors, small RNAs, sigma factors, ppGpp).  Same-operon pairs are excluded.  At strong or confirmed evidence there are 155 regulons, 1,879 genes and 1.76 million scored pairs, with a positive rate of 9.4%.  We report:

- pooled co-membership AUPR;
- per-regulon coherence, the AUROC of within-regulon pairs against pairs of a member and an unrelated gene.

On the 466 replicate-averaged experiments, held-out laboratory splits chose t* = 8 (the two directions gave 10 and 7).  At that depth, attention reached:

- **pooled AUPR:** 0.157 against 0.136 for CLR, +0.021 (95% CI +0.015 to +0.027; paired delete-half jackknife over genes);
- **median coherence:** 0.632 against 0.597, +0.031 (+0.018 to +0.046; paired bootstrap over regulons).

Coherence improved for 69% of 144 regulons (Wilcoxon P = 4 × 10⁻⁴).  Correlation and raw MI were well below CLR (AUPR 0.116 and 0.115).  Results at all evidence levels were directionally identical (+0.007 AUPR, +0.026 coherence; both significant).

Depth behaves as attention depth should (Fig. 2):

- pooled AUPR rises until about t = 32;
- coherence peaks at t = 4–8 and collapses by t = 128, as deep propagation blurs specific programs into diffuse ones.

The label-free criterion selects the range where both improve.  [pending: distribution of t* over 50 random laboratory splits and the full criterion curve, t = 1–64.]

### The gain survives resampling of experiments

Expression compendia contain technical replicates and clusters of related experiments, so gene-level resampling can understate uncertainty.  We therefore resampled at the level of experiments, recomputing dependences, CLR and attention from scratch in every draw.

- **Cluster bootstrap (1,000 draws):** the 466 experiments are drawn with replacement, with one replicate array per occurrence.  Attention at t = 8 improved AUPR by +0.020 (2.5–97.5% interval +0.018 to +0.023) and coherence by +0.028 (+0.010 to +0.046) (Fig. 3).
- **Replicate bootstrap (200 draws):** one array per experiment, removing pseudo-replication.  The result was +0.021 (+0.020 to +0.022).

[pending: nested bootstrap in which t is re-selected inside each of 1,000 draws, so the intervals include the uncertainty of choosing depth.]

### Context normalization, not the dependence measure, drives the gain

To isolate each component, we applied the same softmax and depth to five score matrices (Table 1).

**Table 1 | Component ablation (466 experiments, strong/confirmed evidence, t = 8).**

| Scores fed to the softmax | Context-normalized | Co-membership AUPR | Median coherence |
|---|---|---|---|
| None (CLR scores ranked directly) | yes | 0.136 | 0.597 |
| Raw \|Pearson r\| | no | 0.143 | 0.586 |
| Raw MI | no | 0.143 | 0.583 |
| CLR on \|r\| | yes | 0.155 | 0.622 |
| CLR on \|r\| of normal scores | yes | 0.157 | 0.626 |
| CLR on MI | yes | **0.157** | **0.632** |

[pending: each operator at its own held-out depth over 20 laboratory splits, with jackknife intervals.]

- **Propagation alone is not enough.**  Softmax propagation on uncalibrated scores recovers about a third of the AUPR gain and *reduces* coherence below that of CLR without attention.
- **Context normalization is necessary.**  The full gain on both metrics appears only when attention acts on context-normalized scores.
- **The dependence measure matters little in bulk data.**  CLR on correlation performs like CLR on MI.  In this compendium 99% of selected relationships are positive and monotone, which leaves little non-monotone dependence for MI to capture.

Two further controls agree:

- **Hard thresholding instead of a softmax:** a graph of permutation-FDR-selected CLR edges plus each gene's top 10 neighbours gains less (+0.013 AUPR) and does not improve coherence.  Hard selection alone leaves 62–77% of genes without neighbours, so they cannot propagate.
- **The same hard top-10 propagation on correlation:** no gain in AUPR, and a loss in coherence (−0.014; −0.025 to −0.003).

### Attention ranks targets, while CLR detects edges

Regulon co-membership asks which genes are regulated together.  The 2007 benchmark asks a different question:  is gene X a direct target of transcription factor Y, scored using Y's own mRNA as a proxy for its activity (199 TFs, 3,971 positives, positive rate 0.48%).  On that task the two operators split (Fig. 4):

- **Pooled AUPR, dominated by the strongest few hundred edges:** CLR wins (6.75× random against 5.15× for attention at t* = 8).
- **Each TF's full ranking of candidate targets:** attention wins (median per-TF AUROC 0.623 against 0.595; better for 86 of 149 TFs, and 0.642 at t = 4).

Propagation dilutes the sharpest direct dependences but lifts weaker targets out of the middle of each ranking.  Single-step CLR is the better edge detector, and attention is the better ranker of target sets.  The same pattern appears in single-cell data (below).  [pending: bootstrap intervals for both metrics.]

### Estimated attention identifies perturbation targets beyond differential expression

A stringent test of any network is whether it says something about interventions that the measurements of the intervention do not already say.  Learned attention from single-cell foundation models fails this test (Kendiukhov 2026).  We posed it on M3D, where perturbation-target identification was first benchmarked (Cosgrove et al. 2008).

There are 15 groups of experiments in which a regulator was deleted, overexpressed or mutated.  For each group, the targets are the regulator's RegulonDB regulon, mapped by a table fixed in advance.  Each gene is scored by the influence of the perturbation's experiments on its row of the network:  how much the row changes when those experiments are removed.  The score is calibrated per gene against 500 random removals of the same number of experiments.

A naive evaluation is confounded.  Gene variance alone, which ignores the perturbation entirely, reaches a median AUROC of 0.646 and beats differential expression in 14 of 15 perturbations, because regulon members are more variable than other genes.  We therefore compared targets only with non-targets of similar variance, computing variance without the perturbation's own experiments.  Under this control (Fig. 5):

- **Gene variance** falls to 0.521.
- **Differential expression** scores 0.516.
- **Attention-based influence** keeps a median AUROC of 0.572 and beats differential expression in 12 of 15 perturbations (paired Wilcoxon P = 0.015).
- **CLR influence** reaches 0.543 (P = 0.15).
- **Correlation influence** reaches 0.530 (P = 0.80).

The effect is concentrated in specific regulators (HU, AppY, LexA/SOS, RyhB).  It is absent for the global regulators ppGpp and CRP, whose regulons are diffuse.

[pending: a null drawn from the perturbation's own laboratories; a re-implementation of the SSEM-Lasso network filter of Cosgrove et al. (2008) on the same perturbations and truth.]

### On single-cell data, context normalization transfers and the dependence measure must change

We tested transfer on the BEELINE benchmark (Pratapa et al. 2020):

- seven single-cell RNA-seq datasets;
- its two standard gene sets ("TFs + 500" and "TFs + 1000");
- ground truths from ChIP-seq, STRING, and loss- and gain-of-function screens, plus evidence-filtered STRING subsets we derived.

Three observations emerge (Fig. 6):

- **CLR is the strongest of the methods we evaluated.**  Its median AUPRC ratio is 1.58 over 72 dataset × gene-set × truth combinations, against 1.32 for correlation, 1.20 for MI and 1.08 for regularized partial correlation.
- **The held-out depth transfers.**  It falls at t* = 5–11 in all 14 dataset × gene-set cases, bracketing the bacterial t* = 8, although depth was never tuned on single-cell data.
- **MI is a poor dependence measure for zero-inflated counts.**  The B-spline estimator spends its resolution on the point mass at zero.  Replacing MI by the absolute correlation of normal scores (rank-based Blom scores, with tied zeros sharing one score) improves CLR in 57 of 72 combinations (median +0.115 AUPRC ratio; P = 1 × 10⁻⁶).

Attention on these better-calibrated scores repeats the bacterial pattern:

- it improves on CLR for truths that reward co-regulated target sets, such as transcriptional-regulation STRING (2.17 against 1.74);
- it loses on truths dominated by direct co-expression evidence.

BEELINE also restricts each network to a few hundred to 1,500 genes, 30–60% of them TFs.  That gives CLR's per-gene background an unrepresentative context.  [pending: CLR and attention built on the 5,000 most variable genes and on all expressed genes, with nested cell-bootstrap intervals.]

---

## Discussion

CLR anticipated the central move of attention:  normalize each query's compatibilities within its own context, so that weights mean the same thing across queries with very different scales.  Our results suggest this move, not the choice of compatibility function, is what carries over.  Once the normalization is in place, the rest of the attention machinery adds measurable value:

- a temperature that fixes how many genes each gene attends to;
- aggregation of multi-step attention mass;
- a depth chosen without labels.

Without the normalization, the same machinery is ordinary network propagation, and it is markedly weaker.

Two properties distinguish estimated attention from learned attention:

1. **The weights carry statistical semantics.**  Each score is a standardized deviation from a gene's own background and has a permutation null, so attention can be thresholded with a controlled false-discovery rate and audited edge by edge.
2. **Nothing is learned.**  No corpus, no embeddings and no labels are needed.  Domains where elements have no fixed meaning, as genes do not across tissues and conditions, are exactly where learned embeddings are least transferable.

That learned attention from single-cell foundation models does not beat trivial features on perturbation targets, while estimated attention does here, suggests the two are complementary.  [pending, or future work: a direct comparison on a shared task.]

The gains are modest in absolute terms:  about 0.02 AUPR, 15% relative.  They are, however, consistent across resampling schemes, evidence levels and tasks.  The failures are informative too.  Attention dilutes the sharpest direct edges and blurs diffuse global programs (ppGpp, CRP), and it does not help when the truth rewards direct co-expression.  Practitioners should therefore use single-step CLR to nominate direct edges, and attention to rank modules and target sets.

**Limitations.**

- RegulonDB evidence partly derives from expression studies.  The strong/confirmed tier limits this, and the bias applies equally to every method compared.
- The depth criterion predicts held-out CLR structure, which is a proxy for biological accuracy, not the same thing.
- The main bacterial result rests on one compendium.  A second organism with regulon-level truth would strengthen it [pending].

---

## Methods

**Data.**

- *E. coli*:  M3D *E. coli* v4 build 6, 4,297 genes; 466 replicate-averaged experiments as the primary set, and 907 arrays for the replicate-aware bootstraps.  Ground truth is RegulonDB 14.5 (RISet, TUSet, OperonSet, NetworkSigmaGene), with 2008-era gene names mapped to b-numbers through NCBI synonyms.
- Single-cell:  BEELINE datasets and networks (Zenodo 3701939).  Graded binding comes from ChIP-Atlas.  Evidence-filtered STRING truths come from STRING v11.0 transcriptional actions and the v12.0 experimental and database channels.

**Dependence estimation.**

- Mutual information uses the B-spline estimator of Daub et al. (2004), spline order 3.  One common bin count is chosen for the joint histogram by the low-bias rule of Hacine-Gharbi et al. (2012); per-gene adaptive bins inflate MI for heavy-tailed genes.
- Normal scores use the Blom transform, Φ⁻¹((r − 3/8)/(n + 1/4)), with average ranks for ties.

**Context normalization (CLR).**  Each row of the dependence matrix is z-scored against that gene's own background (sample standard deviation, negative values clipped to zero).  The two directions are combined by Stouffer's rule, (z_ij + z_ji)/√2, which has a standard-normal null; the original Euclidean rule gives nearly identical rankings (Spearman 0.993).  For edge selection, samples are permuted independently per gene, CLR is recomputed, and Benjamini–Hochberg selection is applied at q = 0.05.  The edge sets are unchanged between 100 and 1,000 permutations.

**Attention operator.**  For each row:

- candidates are the 50 highest positive scores;
- weights are a softmax with temperature found by bisection so that exp(entropy) = 10;
- rows sum to one, and rows without candidates keep a self-loop;
- P = (1 − α)I + αA with α = 0.5, read out as ½(Pᵗ + Pᵗᵀ).

[pending: candidate caps of 100, 200 and none.]

**Held-out depth.**

- *Bulk data:*  experiments are split into two balanced halves by laboratory (experimenter), and each half is used in turn for training.  The criterion is the AUPR of the training-half attention mass for the top 1% of the held-out half's scores, computed from the same dependence measure.  t* is the geometric mean of the two directions' optima.
- *Single-cell data:*  cells are split into random halves.
- RegulonDB and BEELINE truths are never used to choose t.

**Uncertainty.**

- Paired delete-half jackknife over annotated genes (pooled AUPR).
- Paired bootstrap over regulons (coherence).
- Cluster and replicate bootstraps over experiments, with all quantities recomputed in each draw (percentile intervals).
- Nested cell bootstraps for BEELINE [pending].

**Perturbation influence.**  For a perturbation with k experiments, influence is the root-mean-square change in a gene's network row when those experiments are removed.  It is z-scored against 500 random removals of k experiments.  The variance-stratified AUROC pools Mann–Whitney statistics within deciles of gene variance, with variance computed without the perturbation's experiments.

**Software.**  This is the R package `clr`:  a C++17 core with OpenMP and a thin Rcpp layer.  Analyses are scripted in `analysis/` and run in parallel.  Code is available at [repository].

---

## Figure legends (planned)

**Fig. 1 | CLR as attention.**  Dependence estimation, per-query context normalization (CLR), softmax with a per-query temperature, and iterated attention, alongside the corresponding transformer components.

**Fig. 2 | Depth.**  Co-membership AUPR and regulon coherence against depth t for attention on CLR, attention on uncalibrated scores, and hard-threshold propagation.  The held-out t* distribution over laboratory splits is overlaid.

**Fig. 3 | Robustness.**  Differences from CLR with 95% intervals under the gene jackknife, regulon bootstrap, cluster and replicate bootstraps, and nested bootstrap.

**Fig. 4 | Edge detection against target ranking.**  Pooled TF→target AUPR and per-TF AUROC for CLR and attention across depths.

**Fig. 5 | Perturbation targets.**  Raw and variance-stratified AUROC per perturbation for differential expression, gene variance, and correlation-, CLR- and attention-based influence.

**Fig. 6 | Single-cell transfer.**  AUPRC ratio by ground truth for MI-based CLR, normal-score CLR and attention on each, across gene contexts.

---

## References (working list)

Bahdanau, Cho & Bengio 2014; Cosgrove et al. 2008, *Bioinformatics*; Cowen et al. 2017, *Nat. Rev. Genet.*; Daub et al. 2004, *BMC Bioinformatics*; Faith, Hayete et al. 2007, *PLoS Biol.*; Faith et al. 2008, *Nucleic Acids Res.* (M3D); Hacine-Gharbi et al. 2012; Kendiukhov 2026, *BMC Genomics*; Kondor & Lafferty 2002, *ICML*; Marbach et al. 2012, *Nat. Methods* (DREAM5); Nadaraya 1964; Pratapa et al. 2020, *Nat. Methods*; Schlag, Irie & Schmidhuber 2021, *ICML*; Schmidhuber 1991; Vaswani et al. 2017; Watson 1964; Zhou & Schölkopf 2004.
