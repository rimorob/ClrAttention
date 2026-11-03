# CLR as attention: checkpoint for Jim Collins

*Boris, 21 October 2026. Status of the project, what we have in hand, what did not work, and the decisions ahead. The key results below have been checked by an independent statistical review; numbers are from the reviewed outputs. **GRNBoost2 rows regenerated 29–31 Oct 2026 (legacy early stopping, after an early-stopping fault was found on 22 Oct); GRNBoost2 claims reviewed 31 Oct. The shared claude.ai page was republished with these values on 31 Oct 2026 (version 3; it also carries the draft abstract and graphical abstract, which are not in this file).***

## In one paragraph

CLR (2007) standardizes each gene's mutual information against that gene's own background (the 2007 code estimates MI by B-spline smoothing, as we do, and combines the two directions by a Euclidean formula; the softmax and the walk are new). That is a context-normalized attention score, published a decade before transformers. It remains a standard baseline for de novo, unsupervised network inference. We complete it as an attention mechanism, with no learned parameters: a softmax over each gene's top CLR partners, then a short multi-step walk (depth t, chosen without labels on held-out experiments). On *E. coli* compendia this recovers strong-evidence transcriptional regulons better than CLR and better than GENIE3, PIDC and GRNBoost2 in pooled ranking. The gain requires context normalization, it survives controls for batch structure and, on strong-evidence regulons, for hub genes, and it concentrates in large regulatory programs. Per regulon, attention ties CLR and GENIE3. It loses on tight co-expression cliques such as ribosomal complexes. A mechanism consistent with these results, though not isolated by them: the walk trades pair-level resolution for module-level recovery.

## What we have in hand

**1. The headline (E. coli, M3D, 466 experiments; RegulonDB regulon co-membership)**
- Attention beats CLR by +0.022 AUPR [+0.019, +0.026]. Depth is re-selected without labels in each of 1,000 nested bootstrap draws, and attention is ahead in all 1,000.
- The result is consistent across 50 laboratory splits (+0.022 [+0.017, +0.027]) and across a replicate-aware bootstrap on the 907 chips (+0.020).
- The held-out depth is stable: median t = 9, with the criterion flat within 1% for t = 7–11.

**2. Why it works**
- Context normalization is essential. A softmax on uncalibrated |r| or MI gets about a third of the gain and none of the coherence gain.
- A Shapley decomposition over six components attributes the AUPR gain mainly to depth (+0.018) and to CLR normalization (+0.013).

**3. Controls**
- Removing 0–20 global principal components leaves the gain intact. The top components are biology, not noise: growth rate, motility (σ28) and the σS stationary-phase program.
- Hub control: the walk converges toward well-connected genes, so we score every method within strata of hubness. Attention stays first on strong-evidence truths (ahead of CLR in every hub decile on M3D); on all-evidence truths it ties GENIE3.

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

- **Per regulon** (each regulon weighted equally): attention ties CLR and GENIE3. Its mean AP is higher than |r|'s in 3 of 4 comparisons after correction and than PIDC's and GRNBoost2's in 1 of 4 each; medians are about equal. Its win rate is highest in regulons of 25 or more genes.
- **Direct TF→target edges:** attention trails CLR. The walk trades edge precision for module recovery.

**5. Human data**
- K562 Perturb-seq pseudobulk: every method is within 10% of random on TF truths. Attention's small DoRothEA lead over CLR survives hub control (it ties GENIE3 and PIDC). It loses to CLR, GENIE3 and PIDC on protein complexes, and on knockdown-response modules per module and after hub control; 86% of the complex loss comes from ten ribosomal and pre-rRNA complexes. GRNBoost2, which does not aggregate, does about as poorly as attention on these complexes, so the loss is not specific to the walk.
- Single cells (BEELINE: human ESC and hepatocytes, plus five mouse datasets): CLR on normal scores (CLR-NS) has the highest per-dataset mean of the methods tested, including GENIE3, PIDC and GRNBoost2, but it was chosen on these truths, it is not significantly ahead of GENIE3 or PIDC, and GENIE3 has the higher median. Attention adds nothing there.
- **Against learned attention (preregistered, run once, 2 Nov):** on all seven BEELINE datasets, CLR-NS, a formula with no training, is non-inferior and superior to zero-shot attention from scPRINT, a single-cell foundation model trained for this task (mean log2 AUPRC-ratio difference +1.03, exact interval +0.34 to +1.90, positive in 7 of 7; holds when scPRINT gets its best readout). Read it with three qualifications: scPRINT's attention scores the same as the same model run on permuted expression, and below a data-free gene-embedding similarity, so the comparator is weak here; the advantage is absent on cell-type-specific ChIP-seq (scPRINT higher in 11 of 14 cases, exploratory); and scGPT (two datasets, descriptive) is behind CLR-NS on attention but ahead on its token embeddings. On K562 control cells the direction is the same but smaller (+0.27, 0.09 to 0.46, descriptive). Fine-tuning a foundation model was not tested. Details: `docs/scfm_scoring.md`.

**6. Perturbation targets (likely a second paper)**
- Attention-based influence ranks a perturbed TF's targets above differential expression in 12 of 15 perturbations (p = 0.015; 9 of 15, p = 0.095, under a laboratory-matched null). Suggestive, not established.
- It is complementary to SSEM-Lasso, re-implemented from Cosgrove et al. 2008: SSEM pinpoints the regulator itself, attention its regulon.
- Differential CLR is the natural method here and is next to port.

**7. Batch effects** (supplement)
- Multi-platform EcoMAC shows platform "Christmas trees" in gene-pair plots.
- Platform-level ComBat helps more than switching from correlation to MI.
- One clean platform (735 arrays) beats all 2,198 pooled.

## What did not work (reported as negatives)
- Hard FDR-thresholded operators: they isolate 62–77% of genes.
- Value transforms and diffused-profile readouts.
- Wider gene contexts on single cells.
- Data-driven choice of attention width: all three options failed review before running, and are being tested empirically anyway under a preregistration.
- A label-free multiscale algorithm based on spectral plateaus: there are no plateaus; the structure is a continuum.
- Depth as a discrete "resolution dial" ordering operons → regulons → sigmulons: not supported on a linear depth grid. The best depth does scale with module size, roughly as size^0.4, but that is generic to diffusion. A log-scale re-test is planned.

## In progress
- Learned attention (done 2 Nov; see Results below).
- GENIE3, GRNBoost2 and PIDC on the BEELINE single-cell benchmark: done. CLR is comparable to GENIE3 and PIDC. GRNBoost2 (legacy early stopping) is below CLR in 6 of 7 datasets (mean log2 −0.15; Wilcoxon P = 0.031, sign P = 0.125; hESC is the exception), a descriptive count, not a significant difference.
- Preregistered attention-width test, including multi-head attention as the candidate fix for the complex loss.
- Directed DREAM5-style evaluation with ARACNE; log-scale depth analysis; perturbation-diversity learning curves; differential CLR.

## Decisions for you and me
1. **Venue:** Nature Machine Intelligence (150-word abstract) or Nature Communications.
2. **Scope:** one paper on CLR attention, with perturbation targets and differential CLR as a second paper.
3. **The central claim:** "better recovery of regulatory programs from expression alone", not "better network inference". Is that strong enough for the target venue?
4. **The human single-cell story:** CLR-NS rather than attention is the single-cell result. Does it belong in this paper?
5. **AI disclosure:** the analysis code, runs and drafts were produced with Claude (Anthropic), with independent AI statistical reviews, under my direction. A contributions statement is drafted to match Nature Portfolio policy.
