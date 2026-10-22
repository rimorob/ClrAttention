# Learned-attention comparator: plan (2026-10-21)

*Design memo by the statistics/scFM reviewer (Opus, high effort), at the PI's request. Citations were retrieved by the reviewer from the web and must be verified (DOIs, versions) before they enter the paper.*

## Choice of comparators
- **Primary: scPRINT v1** (Kalfon, Samaran, Peyré & Cantini 2025, Nat Commun 16:3607, doi:10.1038/s41467-025-58699-1): the scFM built for GRN inference from attention; human + mouse pretraining, so all 7 BEELINE datasets; runs on Pascal GPUs without flash attention. In its own benchmark it beat scGPT, Geneformer v2, GENIE3 and DeepSEM on McCalla hESC/mESC truths, but GENIE3 beat it on Replogle genome-wide.
- **Secondary: scGPT whole-human** (Cui et al. 2024, Nat Methods 21:1470): the model referees will name; scored by attention and by token-embedding cosine (the best scFM signal in the independent BEELINE benchmark of Ye et al. 2026, bioRxiv doi:10.64898/2026.05.17.725083). Token embeddings are dataset-independent: a learned prior, not inference from the data; say so.
- Optional: Geneformer V2-104M, scPRINT-2. Not worth running: TranscriptFormer (too large; prompting), GeneCompass (feeds prior GRNs: leakage), UCE/CellFM/Nicheformer (no GRN readout).
- Sceptical literature to cite: Kendiukhov 2026 (arXiv:2602.17532), VCBench (Weidener et al. 2026), Qi et al. ICML 2026 (arXiv:2605.08128), Ahlmann-Eltze et al. 2025, Boiarsky et al. 2024, Kedzierska et al. 2025.

## Design (to preregister before scoring)
1. Same cells and genes for every method: all BEELINE cells; 5 seeded subsamples of 1,024 K562 control cells. Evaluate on genes every model covers; report dropped genes. CLR, CLR-NS, |r| rerun on exactly those inputs.
2. Same truths and scoring: analysis/beeline_baselines.R (directed TF-node AUPRC ratio and tie-aware EPR; primary vs derived truths).
3. Orientation fixed in advance: primary score(TF→g) = A[g, TF] (target as query, TF as key); secondary the reverse and max(A, Aᵀ).
4. Arms: primary zero-shot attention (mean over all heads and layers as each model's default; special tokens removed; rows renormalized); secondary zero-shot embeddings (scGPT token/hidden-state cosine); tertiary self-supervised fine-tuning per dataset (no labels). No head selection trained on truths (leakage).
5. Controls: same architecture with random weights (fully random; pretrained gene embeddings kept); expression permuted across cells (if the network does not change, it is a prior); gene-embedding-only prior; existing baselines (GENIE3, GRNBoost2, PIDC, |r|); hybrid CLR applied to learned attention (thesis test).
6. Statistics: per-dataset blocking (n = 7 BEELINE datasets), exact Wilcoxon and sign test (minimum two-sided p 0.0156); Holm over {scPRINT, scGPT} × {CLR-NS, CLR attention}; preregistered non-inferiority margin. **Power problem to fix before preregistering (flagged 2026-10-21):** with 7 datasets the smallest two-sided p (0.0156) exceeds Holm's first threshold for 4 tests (0.0125), so no comparison can reach significance. Options: fewer primary tests (one primary pair, scPRINT vs CLR-NS, others secondary), one-sided tests with a stated direction, adding independent blocks (K562 subsamples, other human datasets), or a mixed model over datasets × truths with the dataset as a random effect.
7. Mouse: scPRINT on all 7; scGPT on the 2 human datasets (ortholog mapping exploratory).

## Risks
Pretraining leakage (audit GEO/SRA accessions of BEELINE and Replogle against cellxgene census and scBaseCount; report), truth leakage via head selection, vocabulary loss (Ensembl mapping; Geneformer drops zero-count genes), attention sinks and hub keys (remove special tokens; CLR hybrid; stratify by cluster), normalization (scFMs want raw UMI counts; BEELINE ships log TPM/FPKM, so fetch raw counts from GEO).

## Steps (about 2 weeks of work, under 2 GPU-days)
1. Separate conda env with pinned PyTorch (CUDA build supporting Pascal), weights, smoke test on hESC: 1-2 days.
2. Raw counts, gene-ID mapping, coverage table, leakage audit: 2 days.
3. Extraction script through the scheduler (GPU jobs declared): about 1 GPU-day with controls: 2 days.
4. Preregistration committed before scoring: 0.5 day.
5. Scoring, statistics, figures: 2 days.
6. Optional fine-tuning arm: 3 days.
