# Preregistration: learned-attention comparator (single-cell foundation models)

*Fixed 2026-10-21, before any foundation-model network was computed or scored against a truth. The design memo of 2026-10-21 was revised after a reviewer found that the earlier plan could not reach significance: 4 Holm-corrected tests on 7 datasets have a smallest exact two-sided p of 0.0156, above Holm's first threshold of 0.0125. This version replaces that plan with a single primary contrast, read from one confidence interval. Any change after this commit is logged below with its date and reason.*

**Disclosure.** CLR-NS was chosen as our single-cell variant on these BEELINE truths (results/beeline_blom_test), and its scores are already known. CLR, which was fixed in advance, is therefore a key secondary reference. Every foundation-model readout below is fixed without looking at any truth.

## Resolution of the power problem

| Option | Decision | Reason |
|---|---|---|
| One primary contrast | **Adopted**: CLR-NS vs scPRINT, one endpoint | No multiplicity in the primary family. Everything else is a gatekept secondary or descriptive. |
| Non-inferiority (NI) with a margin on log2 AUPRC ratio | **Adopted** as the primary hypothesis | The claim is that a training-free operator is not worse than learned attention, which is an NI claim. Superiority is tested next in a fixed sequence (no α spent), and every other conclusion is read off the same interval. |
| One-sided tests with a stated direction | **Adopted** through the interval | The NI hypothesis is directional. The interval is the intersection of two one-sided tests, each at exact level 5/128 = 0.039. |
| More independent blocks (K562 subsamples, RPE1, new human ChIP datasets) | **Rejected as primary**. K562 is a replication stratum. | Subsamples of K562 control cells repeat one cell line and one truth set, so they are pseudo-replicates, not independent blocks. RPE1 single cells are not on disk. Adding human datasets after the fact invites selection. |
| Mixed model over datasets × truths with dataset random effects | **Rejected as primary**. Sensitivity only. | A variance component estimated from 7 levels is unstable and its degrees of freedom are ambiguous. Truths within a dataset share TFs and genes, so treating them as replicates overstates the information. |

**Attainability.** With n = 7, the smallest exact one-sided p is 1/128 = 0.0078. The one-sided test at α = 0.05 rejects iff W⁺ ≥ 25, with exact size 5/128 = 0.039. Simulation of the decision rule (normal per-dataset differences) gives:

| Per-dataset SD | True difference | P(NI shown) |
|---|---|---|
| 0.15 | 0 | 0.69 |
| 0.15 | +0.1 | 0.97 |
| 0.10 | 0 | 0.94 |

The SD of 0.05–0.17 comes from existing CLR-family contrasts on these 14 cases.

## Data

- **BEELINE** (Pratapa et al. 2020).
  - Datasets, with all cells (no subsampling): hESC 758, hHep 425, mDC 383, mESC 421, mHSC-E 1,071, mHSC-GM 889, mHSC-L 847.
  - Gene sets: TFs+500 and TFs+1000, from `select_genes()` in analysis/beeline.R.
  - Gene filter: genes with zero range are dropped, as in `case_X()` of analysis/beeline_baselines.R.
  - All 14 cases are used.
- **Model input.**
  - Values are the full ExpressionData matrix (all genes), back-transformed to the linear scale with BEELINE's documented transform. The transform is checked at the smoke test and logged.
  - Each model then applies its own normalization.
  - The input token set is the case gene set, including zeros.
  - Raw GEO counts are an optional sensitivity analysis.
- **Coverage.** Case genes are mapped to each model's vocabulary (scPRINT: Ensembl IDs per organism; scGPT: HGNC symbols).
  - All methods in a contrast are fitted and scored on the intersection G∩ of covered genes, using the same cells. CLR, CLR-NS and |r| are recomputed on G∩.
  - Report the dropped genes and dropped truth edges per case.
  - A dataset is dropped only if the model cannot run on it. If that happens, the interval below uses the largest C with P(W ≤ C − 1) ≤ 0.05 for the new n.
- **K562 replication stratum** (Replogle et al. 2022, data/perturbseq/K562_essential_raw_singlecell_01.h5ad).
  - Cells: non-targeting control cells only, 5 seeded subsamples of 1,024 cells (seeds 1–5).
  - Genes: CollecTRI TFs with ≥ 10 targets that are detected in ≥ 5% of control cells, plus the 1,000 most variable detected genes.
  - Truths: DoRothEA A–B, CollecTRI, and knockdown-response modules (analysis/human_truths.R, analysis/perturbseq.R).

## Models and weights

- **Primary: scPRINT v1** (Kalfon et al. 2025, *Nat Commun*).
  - Checkpoint `medium-v1.5.ckpt` from HF `jkobject/scPRINT`.
  - `large-v1.ckpt` and `small-v1.ckpt` are secondary.
  - Run with `transformer="normal"` in fp32 on Pascal GPUs.
  - Organism NCBITaxon:9606 for hESC and hHep, NCBITaxon:10090 for the mouse datasets.
  - Record the HF commit, checkpoint SHA-256 and package version in the change log before scoring.
- **Secondary: scGPT whole-human** (Cui et al. 2024, *Nat Methods*).
  - Checkpoint `scGPT_human` (`best_model.pt`, `vocab.json`, `args.json`); record SHA-256.
  - Human cases only: hESC, hHep and K562. No ortholog mapping in the scored analyses.
- **Optional, descriptive:** Geneformer V2-104M and scPRINT-2 (`small-v2.ckpt`).

## Readouts

- **R0 (primary, scPRINT).**
  1. Per cell, compute the attention of every layer and head in fp32.
  2. Remove the rows and columns of non-gene (special, cell or label) tokens.
  3. Renormalize each row to sum 1.
  4. Average over heads, layers and cells to get A (genes × genes, row-stochastic, query × key).
- **Orientation (primary).** score(TF → g) = A[g, TF]: the target is the query and the TF is the key, as scPRINT does with its transpose. The scored matrix is therefore S = Aᵀ (rows are regulators).
- **Secondary readouts:**
  - R1: the scPRINT package's own all-head-average gene-network output, without head selection, with its settings logged.
  - R2: last layer only.
  - Orientations S = A and S = max(A, Aᵀ).
- **scGPT.**
  - Attention: last layer, mean over heads and cells, `<cls>` removed, rows renormalized, same orientation.
  - Token-embedding cosine (data-free).
  - Final-layer hidden-state cosine, averaged over cells.
- **Excluded:** head selection fitted on any truth (OmniPath, McCalla or gwps heads).

## Arms and controls

1. **Zero-shot (primary arm).**
2. **Fine-tuning (descriptive).** Continued self-supervised training on each dataset's own cells, with the package's default losses for 10 epochs and no truths.
3. **Cross-fitted head selection (descriptive).** TFs are split into 2 folds; heads are selected on one fold's edges with scPRINT's ridge procedure and scored on the other fold.
4. **Controls (scPRINT, R0):**
   - Same architecture with random weights, 3 seeds, in two versions: pretrained ESM2 gene embeddings kept, and gene embeddings random.
   - Expression permuted independently per gene across cells.
   - ESM2 gene-embedding cosine (data-free prior).
   - CLR applied to A (hybrid).
   - |r|, CLR and CLR-NS on G∩.
   - GENIE3, GRNBoost2 and PIDC from results/beeline_baselines, descriptive.
5. **CLR attention (CLR-NS soft10 at its held-out depth, analysis/beeline.R) vs scPRINT:** descriptive.

## Scoring

The definitions of analysis/beeline_baselines.R are evaluated, not copied:

- `score_truth`: directed TF-node AUPRC ratio and tie-aware EPR.
- Primary truths: celltype_ChIP, nonspecific_ChIP, STRING, and LOFGOF (mESC only).
- Derived truths (STRING_regulatory, STRING_curated): secondary.
- Only cells (dataset × gene set × truth) complete for both arms of a contrast are used.
- Self-pairs score 0; non-finite values become 0.

## Primary test

d_i = for dataset i, the mean over gene sets and primary truths of log2(AUPRC ratio of CLR-NS / AUPRC ratio of scPRINT R0). Higher values favour CLR-NS.

- **Margin.** δ = 0.15 log2 units, about 11%. This is half of CLR-NS's gain over |r| on these cases (0.30), following the convention that an NI margin preserves 50% of the active control's effect.
- **Interval.** The exact Wilcoxon (Hodges–Lehmann) interval [L, U], where L and U are the 4th and 25th of the 28 ordered Walsh averages of d_1, …, d_7. Its coverage is 1 − 2·5/128 = 92.2%.
- **Primary p.** The exact one-sided signed-rank p for H0: median ≤ −δ, computed on d_i + δ.

## Decision rules

| Result | Conclusion |
|---|---|
| L > −δ | CLR-NS is non-inferior to learned attention (CLR-NS ≥ learned attention within the margin). |
| L > 0 | CLR-NS is superior (tested next in the fixed sequence). |
| −δ < L and U < δ | Tie (equivalence). |
| U < 0 | Learned attention wins. |
| U < −δ | Learned attention wins materially; the paper says so. |
| Any other interval | Inconclusive, reported as such and never as a tie. |

- **Robust NI.** NI is called robust only if it also holds when scPRINT takes, per dataset, the best of {R0, R1, R2} × 3 orientations. That maximum is chosen on the truths, so it is optimistic for scPRINT. It is used only to stress our claim, never to claim a scPRINT win.

## Multiplicity

- **Primary:** one interval and one family.
- **Key secondary.** Interpreted confirmatorily only if primary NI holds; Holm at 0.05 over two tests, each n = 7. The smallest p is 0.0078, below the first Holm threshold of 0.025, so both are attainable.
  - S1: CLR (fixed in advance) vs scPRINT R0, NI.
  - S2: CLR-NS vs scPRINT R0 on EPR, NI.
- **Everything else is descriptive**, with intervals and no claims:
  - K562: mean over subsamples, with a 95% bootstrap over TFs (B = 1,000).
  - scGPT, controls and other readouts.
  - Mixed model `lmer(log2 AUPRC ratio ~ method + N + (1|dataset) + (1|dataset:truth))`.
  - Leaked datasets excluded.

## Leakage audit (done before scoring)

1. Record the evaluation accessions: hESC GSE75748, hHep GSE81252, mDC GSE48968, mESC GSE98664, mHSC GSE81682, and the Replogle K562 deposit.
2. Search each model's pretraining manifest by accession, DOI and title:
   - scPRINT v1: its cellxgene dataset list.
   - scGPT: the CELLxGENE census snapshot it used.
   - scPRINT-2 / Geneformer: their lists, including scBaseCount.
3. Classify each model × dataset pair as in, out, or undeterminable.
4. The primary analysis keeps n = 7, because leakage favours the scFM and so makes NI conservative. Add a sensitivity analysis that excludes leaked datasets, and flag any scFM win on a leaked dataset.
5. Truth leakage: no head selection on truths. CLR-NS's selection on these truths is disclosed above.

## Compute

- GPU jobs go only through the systemd queue, one GPU per job, never started by hand.
- Estimate: about 10–30 GPU-minutes per dataset per readout, and under 2 GPU-days in total including controls and K562.
- Environment: a separate venv with a PyTorch build that supports sm_61.

## Change log

- 2026-10-21: fixed (single primary NI contrast replaces 4 Holm tests).
- 2026-10-22: environment and checkpoints recorded (no truth or scoring touched; the preregistered design is unchanged).
  - Environment: `~/venvs/scfm`, torch 2.3.1+cu121 (sm_60 binary for Pascal, unverified until the smoke test), scprint 1.0.13, scgpt 0.2.4 (`--no-deps`); full freeze in `docs/scfm_env_freeze.txt`.
  - scPRINT `medium-v1.5.ckpt` from HF `jkobject/scPRINT`, commit 8e1ad1a1bc40c1fc2f4f1700fe27923eea9a2bed, SHA-256 96ff7c5d9e1ad5076e818bde3c38f498c95cd354cc42029d90691b5ba2964ff8, stored in `~/scfm_models/scprint/`.
  - scGPT `scGPT_human` (whole-human) from the authors' own release, the Google Drive folder linked in the scGPT repository README (https://github.com/bowang-lab/scGPT, "whole-human"): https://drive.google.com/drive/folders/1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y . The README lists no other source. Downloaded 22 Oct 2026 to `~/scfm_models/scgpt_human/`, SHA-256: `best_model.pt` 6cb5d451ab5c4b33eb673adbe4fddc61d2389df1b89b7651a9fe2e557572b922; `vocab.json` acca93d114ca62c3f0f50debbd23e8c87f0714f4737764454f6b2b13f2e8580f; `args.json` c18e075e018140cb8b2d9029387b9de26607a5ce6a8ccabd6ead70cd76b95d60.
