# CLR as attention: checkpoint for Jim Collins

*Boris, 21 October 2026. Status of the project, what we have in hand, what did not work, and the decisions ahead. The key results below have been checked by an independent statistical review; numbers are from the reviewed outputs.*

## In one paragraph

CLR (2007) standardizes each gene's mutual information against that gene's own background. That is a context-normalized attention score, published a decade before transformers. It remains a standard baseline for de novo, unsupervised network inference. We complete it as an attention mechanism, with no learned parameters: a softmax over each gene's top CLR partners, then a short multi-step walk (depth t, chosen without labels on held-out experiments). On *E. coli* compendia this recovers strong-evidence transcriptional regulons better than CLR and better than GENIE3, PIDC and GRNBoost2 in pooled ranking. The gain requires context normalization, it survives controls for batch structure and, on strong-evidence regulons, for hub genes, and it concentrates in large regulatory programs. Per regulon, attention ties CLR and GENIE3. It loses on tight co-expression cliques such as ribosomal complexes. The working mechanism: the walk trades pair-level resolution for module-level recovery.

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
| M3D | RegulonDB strong evidence | **1.68** | 1.44 | 1.51 | 1.45 | 1.18 |
| EcoMAC (2,198 arrays, mixed platforms) | strong evidence | **1.68** | 1.39 | 1.50 | 1.43 | 1.21 |
| M3D, hub-controlled | strong evidence | **1.56** | 1.40 | 1.47 | 1.40 | 1.28 |
| Human K562 Perturb-seq | DoRothEA A–B | **1.09** | **1.06**† | **1.08** | **1.07** | 1.02 |
| Human K562 Perturb-seq | CORUM complexes | 3.51 | **6.39** | 6.32 | **6.33** | 2.18 |

Full table, with all truths and footnotes: `docs/whitepapers/sota_table.md`.

- **Per regulon** (each regulon weighted equally): attention ties CLR and GENIE3. Its mean AP is higher than GRNBoost2's (4 of 4 comparisons), |r|'s (3 of 4) and PIDC's (1 of 4 after correction); medians are about equal. Its win rate is highest in regulons of 100 or more genes, though the gaps there are tiny.
- **Direct TF→target edges:** attention trails CLR. The walk trades edge precision for module recovery.

**5. Human data**
- K562 Perturb-seq pseudobulk: every method is within 10% of random on TF truths. Attention's small DoRothEA lead over CLR survives hub control (it ties GENIE3 and PIDC). It loses on protein complexes, and on knockdown-response modules per module and after hub control; 86% of the complex loss comes from ten ribosomal and pre-rRNA complexes.
- Single cells (BEELINE: human ESC and hepatocytes, plus five mouse datasets): CLR on normal scores (CLR-NS) is the best of the correlation- and MI-based methods tested so far. Attention adds nothing there.

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
- Learned attention: scPRINT (built for network inference from attention) and scGPT, against CLR-NS and CLR attention, on human and mouse single cells. To be preregistered before scoring, with leakage audits. About two weeks.
- GENIE3, GRNBoost2 and PIDC on the BEELINE single-cell benchmark (queued).
- Preregistered attention-width test, including multi-head attention as the candidate fix for the complex loss.
- Directed DREAM5-style evaluation with ARACNE; log-scale depth analysis; perturbation-diversity learning curves; differential CLR.

## Decisions for you and me
1. **Venue:** Nature Machine Intelligence (150-word abstract) or Nature Communications.
2. **Scope:** one paper on CLR attention, with perturbation targets and differential CLR as a second paper.
3. **The central claim:** "better recovery of regulatory programs from expression alone", not "better network inference". Is that strong enough for the target venue?
4. **The human single-cell story:** CLR-NS rather than attention is the single-cell result. Does it belong in this paper?
5. **AI disclosure:** the analysis code, runs and drafts were produced with Claude (Anthropic), with independent AI statistical reviews, under my direction. A contributions statement is drafted to match Nature Portfolio policy.
