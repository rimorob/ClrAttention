# Context-normalized attention without learning recovers gene regulatory programs

*Manuscript draft in* Nature Machine Intelligence *Article format, 2026-10-14.  Bracketed items marked [pending] depend on analyses now running; all numbers await an independent QC pass.  Working notes, the analysis queue and the journal assessment are in `paper_draft_notes.md`.*

---

## Abstract

Attention weights each element of a set by its normalized compatibility with every other element.  The idea predates neural networks; kernel regression and fast-weight programmers are earlier instances.  In computational biology, context likelihood of relatedness (CLR), a 2007 method that remains a standard for de novo gene-network inference, was arguably the first attention operator in wide use.  CLR computes attention scores from data alone, z-scoring each gene's dependence on every other gene against that gene's own background, but never used them as attention.  Here we bring back to CLR the advances attention research made in language modelling:  a temperature-controlled softmax, iterated attention, and depth chosen on held-out laboratories.  On *Escherichia coli*, regulon recovery improves by 15% and per-regulon coherence rises, robustly across 1,000 resamplings of experiments.  The gain requires CLR's context normalization, not a particular dependence measure.  The estimated attention also identifies perturbation targets beyond differential expression.

---

## Introduction

Attention is now the dominant way in which machine-learning models relate the elements of a set.  For each query, a compatibility with every key is computed, normalized across keys (usually by a softmax with a temperature) and used to aggregate values; stacking the operation lets information travel along multi-step paths.  In transformers the compatibilities are dot products of learned embeddings, and what to attend to is learned from large corpora (Bahdanau et al. 2014; Vaswani et al. 2017).

The pattern is older than its neural form, and we do not claim that it began in biology.  Kernel regression weights observations by their normalized similarity to a query and has been a standard statistical tool since the 1960s (Nadaraya 1964; Watson 1964), and fast-weight programmers are formally equivalent to linear attention (Schmidhuber 1991; Schlag et al. 2021).  Our claim is narrower:  in computational biology, attention entered everyday practice independently, and years before neural attention, through a method for inferring gene regulatory networks.

Context likelihood of relatedness (CLR; Faith, Hayete et al. 2007) estimates the mutual information between every pair of genes in an expression compendium and expresses each value as a z-score against the distribution of that gene's dependences on all other genes.  Read as an operator, each gene is a query, every other gene a key, and the compatibility an estimated statistical dependence rather than a learned dot product, normalized per query against the query's own context.  That normalization is what made CLR work:  a gene weakly dependent on everything is not a hub, whereas a modest dependence that stands out from a gene's own background is a candidate regulatory link.

CLR became a standard tool for de novo network inference, that is, inference from expression data alone.  Methods that bring informative priors, such as pathways, binding motifs or curated interactions, can outperform it, but they answer a different question.  Among expression-only methods, CLR has remained the usual unsupervised baseline for nearly two decades.  It was one of the 'well-established off-the-shelf' methods of the DREAM5 network-inference challenge and a component of several of its ensemble predictors (Marbach et al. 2012), and it is the strongest of the methods we evaluate on a current single-cell benchmark.  It is not the top scorer everywhere:  in DREAM5, regression-based tree ensembles, which fit a predictive model for each gene from a supplied list of candidate regulators, scored higher on directed edges from transcription factors to targets.  CLR is different in kind.  It is model-free, unsupervised estimation of pairwise dependence with per-gene calibration, and it needs neither a regulator list nor a predictive model, so it applies equally when regulators are unknown or latent.  To our knowledge it is the first attention-style operator to be widely adopted in computational biology.

CLR, however, stopped at the scores.  It ranked them as edges and never used the rest of attention:  a temperature-controlled normalization, aggregation and depth.  Over the same two decades, attention research in natural-language processing turned these components into a mature toolkit.  Here we bring that toolkit back to the setting in which CLR has proved durable, and measure what each component adds.  Completing CLR as attention improves recovery of bacterial regulatory programs on the compendium on which CLR was validated, robustly to resampling of experiments.  The improvement requires CLR's per-query context normalization:  softmax propagation on uncalibrated scores recovers little of it, whereas the dependence measure matters less than expected (Supplementary Note 1).  And the estimated attention carries information about the targets of perturbed regulators beyond differential expression, a task on which learned attention from single-cell foundation models adds nothing over trivial features (Kendiukhov 2026).  Because the operator needs no training data, embeddings or labels, and comes with a permutation null, it is a calibrated, interpretable counterpart to learned attention for domains in which the elements have no fixed meaning.

---

## Results

### CLR scores are attention scores, and completing them takes three components

Let D be a G × G matrix of pairwise dependences between genes estimated from an expression compendium (Fig. 1).  CLR computes for each row i the z-scores z_ij = (D_ij − μ_i)/σ_i, clips negative values to zero, and combines the two directions into a symmetric score S_ij = (z_ij + z_ji)/√2.  The row-wise standardization is a per-query normalization and plays the role that softmax plays in attention.

We add the three components that CLR lacked.  First, each row of S becomes attention weights through a softmax with a per-query temperature, A_i· = softmax(S_i·/τ_i), taken over the gene's 50 highest scores; τ_i is chosen so that the effective number of attended genes, exp(entropy), equals 10 for every gene, which makes attention comparably focused for genes whose score distributions differ widely.  Second, attention is iterated through the lazy operator P = ½(I + A), and the readout is the symmetric attention mass between genes after t steps, ½(Pᵗ + Pᵗᵀ).  Third, the depth t is selected without labels:  experiments are split by laboratory, the operator is built on one half, and t is chosen to maximize how well attention mass predicts the strongest CLR pairs of the other half.  No gold standard is used at any stage of model selection.

Iterated attention through a row-stochastic operator is a random walk, and the construction is therefore related to diffusion kernels (Kondor & Lafferty 2002) and to network propagation (Cowen et al. 2017).  We treat this relationship as a baseline to beat rather than a novelty to claim, and every comparison below includes propagation on uncalibrated scores.

### Completed CLR attention improves recovery of *E. coli* regulons

We evaluated the operator on the M3D *E. coli* compendium (4,297 genes; 907 microarrays from 466 experiments; Faith et al. 2008) against RegulonDB 14.5.  The benchmark is regulator-agnostic:  a gene pair is positive if both genes belong to the same regulon of any regulator class (transcription factors, small RNAs, sigma factors and ppGpp), and pairs within an operon are excluded.  At strong or confirmed evidence this gives 155 regulons, 1,879 genes and 1.76 million scored pairs with a positive rate of 9.4%.  We report pooled co-membership AUPR and per-regulon coherence, the AUROC of within-regulon pairs against pairs of a member and an unrelated gene.

On the 466 replicate-averaged experiments, held-out laboratory splits chose t* = 8 (the two directions gave 10 and 7).  At this depth attention reached a pooled AUPR of 0.157 against 0.136 for CLR, a difference of +0.021 (95% CI +0.015 to +0.027, paired delete-half jackknife over genes), and a median coherence of 0.632 against 0.597 (+0.031; +0.018 to +0.046, paired bootstrap over regulons).  Coherence improved for 69% of 144 regulons (Wilcoxon P = 4 × 10⁻⁴).  Correlation and raw mutual information fell well below CLR (AUPR 0.116 and 0.115), and results at all evidence levels were directionally identical (+0.007 AUPR and +0.026 coherence, both significant).

Depth behaves as attention depth should (Fig. 2).  Pooled AUPR rises until about t = 32, whereas coherence peaks at t = 4–8 and collapses by t = 128, as deep propagation blurs specific programs into diffuse ones.  The label-free criterion selects the range in which both improve.  Over 50 random laboratory splits it chose a median t* of 9 (interquartile range 9–10, range 8–11; individual split directions 6–14), and the criterion lay within 1% of its maximum for t = 7–11, so the choice is stable and the optimum flat (Supplementary Note 4).  At t* = 9 the gains were +0.022 AUPR (+0.017 to +0.027) and +0.031 coherence (+0.019 to +0.043).  The held-out depth is conservative relative to the benchmark, whose pooled AUPR peaks near t = 27, so attention is not tuned to the gold standard.

### The gain survives resampling of experiments

Expression compendia contain technical replicates and clusters of related experiments, so resampling genes can understate uncertainty.  We therefore resampled experiments and recomputed dependences, CLR scores and attention from scratch in every draw.  In a cluster bootstrap that draws the 466 experiments with replacement and one replicate array per occurrence (1,000 draws), attention at t = 8 improved AUPR by +0.020 (2.5–97.5% interval +0.018 to +0.023) and coherence by +0.028 (+0.010 to +0.046) (Fig. 3).  A replicate bootstrap that removes pseudo-replication by using one array per experiment (200 draws) gave +0.021 (+0.020 to +0.022) [pending: a nested bootstrap that re-selects t within each of 1,000 draws, so that the intervals include the uncertainty of choosing depth].

### Context normalization drives the gain

To isolate each component we applied the same softmax and depth to five score matrices (Table 1).  Softmax propagation on uncalibrated |r| or mutual information recovered about a third of the AUPR gain and reduced coherence below that of CLR without attention.  The full gain on both metrics appeared only when attention acted on context-normalized scores.  Two further controls agree.  Replacing the softmax by hard selection, a graph of permutation-FDR-selected CLR edges plus each gene's ten strongest neighbours, gained less (+0.013 AUPR) and did not improve coherence; hard selection alone leaves 62–77% of genes without neighbours, which cannot propagate.  The same hard top-10 propagation on correlation gave no gain in AUPR and lost coherence (−0.014; −0.025 to −0.003).

**Table 1 | Component ablation (466 experiments, strong/confirmed evidence, t = 8).**

| Scores fed to the softmax | Context-normalized | Co-membership AUPR | Median coherence |
|---|---|---|---|
| None (CLR scores ranked directly) | yes | 0.136 | 0.597 |
| Raw \|Pearson r\| | no | 0.143 | 0.586 |
| Raw mutual information | no | 0.143 | 0.583 |
| CLR on \|r\| | yes | 0.155 | 0.622 |
| CLR on \|r\| of normal scores | yes | 0.157 | 0.626 |
| CLR on mutual information | yes | **0.157** | **0.632** |

[pending: each operator at its own held-out depth over 20 laboratory splits, with jackknife intervals.]

On this single-platform compendium the dependence measure mattered little:  CLR on correlation performed like CLR on mutual information, consistent with 99% of selected relationships being positive and monotone.  Heterogeneous data behave differently.  On EcoMAC (Carrera et al. 2014), which pools 2,198 arrays from two Affymetrix platforms and other sources, platform offsets bend pairwise relationships into composite shapes that mutual information registers and correlation largely misses.  Without batch correction CLR on mutual information gives higher regulon coherence than CLR on correlation, and it keeps a smaller advantage after per-platform correction; attention helps in every configuration (Supplementary Note 1).

### Attention ranks targets, whereas CLR detects edges

Regulon co-membership asks which genes are regulated together.  The benchmark of the original CLR study asks instead whether gene X is a direct target of transcription factor Y, scored using Y's own mRNA as a proxy for its activity (199 transcription factors, 3,971 positives, positive rate 0.48%).  On this task the two operators divide the labour (Fig. 4).  Pooled AUPR, which is dominated by the strongest few hundred edges, favours CLR (6.75 times random against 5.15 for attention at t* = 8).  Each transcription factor's ranking of all candidate targets favours attention (median per-factor AUROC 0.623 against 0.595, better for 86 of 149 factors, and 0.642 at t = 4).  Propagation dilutes the sharpest direct dependences but lifts weaker targets out of the middle of each ranking:  single-step CLR is the better edge detector and attention the better ranker of target sets.  The same division appears in single-cell data [pending: bootstrap intervals for both metrics].

### Estimated attention identifies perturbation targets beyond differential expression

A stringent test of a network is whether it says something about an intervention that the measurements of the intervention do not already say.  Learned attention from single-cell foundation models fails this test (Kendiukhov 2026).  We posed it on M3D, where the identification of perturbation targets was first benchmarked (Cosgrove et al. 2008).  For 15 groups of experiments in which a regulator was deleted, overexpressed or mutated, the targets are the regulator's RegulonDB regulon, mapped by a table fixed in advance.  Each gene is scored by the influence of the perturbation's experiments on its row of the network, that is, by how much the row changes when those experiments are removed, calibrated per gene against 500 random removals of the same number of experiments.

A naive evaluation is confounded.  Gene variance, which ignores the perturbation entirely, reaches a median AUROC of 0.646 and beats differential expression in 14 of 15 perturbations, because regulon members are more variable than other genes.  We therefore compared targets only with non-targets of similar variance, computing variance without the perturbation's own experiments.  Under this control gene variance falls to 0.521 and differential expression scores 0.516, whereas attention-based influence keeps a median AUROC of 0.572 and beats differential expression in 12 of 15 perturbations (paired Wilcoxon P = 0.015); influence on CLR and on correlation reaches 0.543 (P = 0.15) and 0.530 (P = 0.80) (Fig. 5).  The effect is concentrated in specific regulators (HU, AppY, LexA/SOS and RyhB) and absent for the global regulators ppGpp and CRP, whose regulons are diffuse [pending: a null drawn from the perturbation's own laboratories, and a re-implementation of the SSEM-Lasso network filter of Cosgrove et al. (2008) on the same perturbations and truth].

### On single-cell data, context normalization transfers and the dependence measure must change

We tested transfer on the BEELINE benchmark (Pratapa et al. 2020), which comprises seven single-cell RNA-sequencing datasets, two standard gene sets and ground truths from ChIP-seq, STRING and loss- and gain-of-function screens, to which we added evidence-filtered STRING subsets.  CLR is the strongest of the methods we evaluated, with a median AUPRC ratio of 1.58 over 72 combinations of dataset, gene set and truth, against 1.32 for correlation, 1.20 for mutual information and 1.08 for regularized partial correlation (Fig. 6).  The held-out depth transfers, falling at t* = 5–11 in all 14 dataset and gene-set cases and bracketing the bacterial t* = 8, although depth was never tuned on single-cell data.

Mutual information is, however, a poor dependence measure for zero-inflated counts, because the B-spline estimator spends its resolution on the point mass at zero.  Replacing it by the absolute correlation of normal scores (Blom scores, with tied zeros sharing one score) improves CLR in 57 of 72 combinations (median +0.115 AUPRC ratio; P = 1 × 10⁻⁶).  Attention on these scores repeats the bacterial division of labour:  it improves on CLR for truths that reward co-regulated target sets, such as transcriptional-regulation STRING (2.17 against 1.74), and loses on truths dominated by direct co-expression evidence.  BEELINE also restricts each network to a few hundred to 1,500 genes, 30–60% of them transcription factors, which gives CLR's per-gene background an unrepresentative context [pending: CLR and attention built on the 5,000 most variable genes and on all expressed genes, with nested cell-bootstrap intervals].

---

## Discussion

CLR anticipated the central move of attention:  normalize each query's compatibilities within its own context, so that weights mean the same thing across queries on very different scales.  Our results suggest that this move, more than the choice of compatibility function, is what carries over.  Once the normalization is in place, the rest of the attention machinery adds measurable value, namely a temperature that fixes how many genes each gene attends to, aggregation of multi-step attention mass and a depth chosen without labels.  Without the normalization, the same machinery is ordinary network propagation and is markedly weaker.  The dependence measure becomes important where the data demand it:  mutual information where batch structure bends relationships, and rank-based correlation where counts are zero-inflated.

Estimated attention differs from learned attention in two respects.  Its weights carry statistical semantics:  each score is a standardized deviation from a gene's own background with a permutation null, so attention can be thresholded at a controlled false-discovery rate and audited edge by edge.  And nothing is learned:  no corpus, embeddings or labels are needed, which matters most where elements have no fixed meaning, as genes do not across tissues and conditions.  That learned attention from single-cell foundation models does not beat trivial features on perturbation targets, while estimated attention does here, suggests that the two are complementary [pending, or future work: a direct comparison on a shared task].

The gains are modest in absolute terms, about 0.02 AUPR or 15% relative on M3D, but consistent across resampling schemes, evidence levels and tasks, and larger on the heterogeneous EcoMAC compendium.  The failures are informative too.  Attention dilutes the sharpest direct edges, blurs diffuse global programs such as those of ppGpp and CRP, and does not help when the truth rewards direct co-expression.  In practice, single-step CLR should be used to nominate direct edges and attention to rank modules and target sets.

Several limitations apply.  RegulonDB evidence derives in part from expression studies; the strong/confirmed tier limits this, and any bias applies equally to every method compared.  The depth criterion predicts held-out CLR structure, a proxy for biological accuracy rather than the thing itself.  And the main bacterial results rest on two compendia from one organism; a second organism with regulon-level truth would strengthen them.

---

## Methods

**Data.**  The *E. coli* analyses use M3D *E. coli* v4 build 6 (4,297 genes), with the 466 replicate-averaged experiments as the primary set and the 907 arrays for the replicate-aware bootstraps, and EcoMAC (Carrera et al. 2014; 2,198 arrays, 4,189 genes), with each array's Affymetrix platform taken from its GEO record.  The ground truth is RegulonDB 14.5 (RISet, TUSet, OperonSet and NetworkSigmaGene), with 2008-era gene names mapped to b-numbers through NCBI synonyms.  Single-cell data and networks are from BEELINE (Zenodo 3701939); graded binding is from ChIP-Atlas, and the evidence-filtered STRING truths use STRING v11.0 transcriptional actions and the experimental and database channels of v12.0.

**Dependence estimation.**  Mutual information uses the B-spline estimator of Daub et al. (2004) with spline order 3 and one common bin count chosen for the joint histogram by the low-bias rule of Hacine-Gharbi et al. (2012); per-gene adaptive bins inflate the mutual information of heavy-tailed genes.  Normal scores use the Blom transform, Φ⁻¹((r − 3/8)/(n + 1/4)), with average ranks for ties.

**Context normalization.**  Each row of the dependence matrix is z-scored against the gene's own background (sample standard deviation, negative values clipped to zero), and the two directions are combined by Stouffer's rule, (z_ij + z_ji)/√2, which has a standard-normal null; the original Euclidean rule gives nearly identical rankings (Spearman 0.993).  For edge selection, samples are permuted independently for each gene, CLR is recomputed and Benjamini–Hochberg selection is applied at q = 0.05; the selected edge sets are unchanged between 100 and 1,000 permutations.

**Attention operator.**  For each row, the candidates are the 50 highest positive scores and the weights are a softmax whose temperature is found by bisection so that exp(entropy) = 10.  Rows sum to one, rows without candidates keep a self-loop, and P = (1 − α)I + αA with α = 0.5 is read out as ½(Pᵗ + Pᵗᵀ).  Caps of 100 and 200 candidates give the same depth and gain as 50; an uncapped softmax selects shallower depths and gains less AUPR (Supplementary Note 4).

**Held-out depth.**  For bulk data, experiments are split into two balanced halves by laboratory (experimenter) and each half is used in turn for training.  The criterion is the AUPR of the training-half attention mass for the top 1% of the held-out half's scores computed from the same dependence measure, and t* is the geometric mean of the two directions' optima.  For single-cell data, cells are split into random halves.  RegulonDB and BEELINE truths are never used to choose t.

**Uncertainty.**  Pooled AUPR is assessed with a paired delete-half jackknife over annotated genes and coherence with a paired bootstrap over regulons.  Cluster and replicate bootstraps over experiments recompute every quantity in every draw and report percentile intervals [pending: nested cell bootstraps for BEELINE].

**Perturbation influence.**  For a perturbation with k experiments, influence is the root-mean-square change in a gene's network row when those experiments are removed, z-scored against 500 random removals of k experiments.  The variance-stratified AUROC pools Mann–Whitney statistics within deciles of gene variance, computed without the perturbation's experiments.

**Software.**  The R package `clr` has a C++17 core with OpenMP and a thin Rcpp layer; analyses are scripted in `analysis/` and run in parallel.  Code is available at [repository].

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

Bahdanau, Cho & Bengio 2014; Cosgrove et al. 2008, *Bioinformatics*; Cowen et al. 2017, *Nat. Rev. Genet.*; Daub et al. 2004, *BMC Bioinformatics*; Faith, Hayete et al. 2007, *PLoS Biol.*; Carrera et al. 2014, *Mol. Syst. Biol.* (EcoMAC); Faith et al. 2008, *Nucleic Acids Res.* (M3D); Hacine-Gharbi et al. 2012; Kendiukhov 2026, *BMC Genomics*; Kondor & Lafferty 2002, *ICML*; Marbach et al. 2012, *Nat. Methods* (DREAM5); Nadaraya 1964; Pratapa et al. 2020, *Nat. Methods*; Schlag, Irie & Schmidhuber 2021, *ICML*; Schmidhuber 1991; Vaswani et al. 2017; Watson 1964; Zhou & Schölkopf 2004.
