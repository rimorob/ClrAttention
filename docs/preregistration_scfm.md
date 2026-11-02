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
- 2026-10-31: GRNBoost2 baseline replaced (no truth, scoring rule, test or decision rule touched; no foundation-model network or scfm score existed at the time of this entry). Analysis 4 uses GENIE3, GRNBoost2 and PIDC from `results/beeline_baselines` as descriptive comparators only. The GRNBoost2 rows in that directory are now the legacy-early-stopping rerun (arboreto 0.1.6, `tools/grnboost2.py --early_stop legacy`, production seed 20261021; BEELINE rerun finished 28 Oct 2026, jobs 0922-0924), not the first run reported 16-21 Oct 2026.
  - Reason and order of events: arboreto's early-stopping rule is broken under scikit-learn 1.3 and later (found 22 Oct 2026; the first GRNBoost2 run used scikit-learn 1.9.1 and stopped most targets after 25-36 trees). A diagnostic designed after that discovery (job 0894) showed the wrapper and the seed were not at fault. The rerun rule (legacy early stopping everywhere; no early stopping as sensitivity only; no other comparator rerun) was committed on 27 Oct 2026 (050695d) before any rerun, with the variant results already seen. Legacy early stopping reproduced stock arboreto 0.1.6 under scikit-learn 1.2.2 on BEELINE mDC_500 (commit 45a183b).
  - Gate: in the BEELINE rerun (job 0923) the |r|, CLR, CLR-NS, GENIE3 and PIDC rows reproduced the earlier results exactly (observed difference 0, tolerance 1e-6). The superseded GRNBoost2 outputs are kept in `results/grnboost2_superseded/`; full account in `docs/plans/grnboost2_diagnostic.md`, paper Methods and Supplementary Note 13.
  - GRNBoost2 on the 14 BEELINE cases (primary truths, vs CLR): median AUPRC ratio 1.47, mean log2 -0.15, better than CLR in 1 of 7 datasets, Wilcoxon P = 0.031, sign P = 0.125 (`results/beeline_baselines/summary.csv`); statistical review of the regenerated GRNBoost2 results pending.
- 2026-11-01: step-2 decisions (PI confirmed; no truth file opened, no foundation-model network computed for scoring). The preparation steps 1 and 2 (`docs/scfm_smoke_test.md`, `docs/scfm_step2.md`, `results/scfm_step2/`) opened no truth file and scored no network. The eight decisions D1 to D8 of `docs/plans/scfm_step2_prereg_draft.md` (commit 338f80b) were confirmed by the PI on 2026-11-01, as proposed, with two explicit calls: (1) **no SRA re-quantification** of hHep, mDC and mESC before scoring (so arm B has n = 4 and the primary stays arm A); (5) **R0' (heads with special-token mass < 0.5, at least 4 heads) is included in the robust-NI maximum** (E11). All other items are as proposed. The entries E1 to E12 below are committed before any production extraction and are the freeze point for these choices; the primary test, the margin, the interval, the decision table and the key secondaries are unchanged. Editorial changes to the draft wording: the "proposed, PI to confirm" tags were removed from the headings, and the E6 sentence on SRA re-quantification now records the PI decision; nothing else.
  Order of work, step 3 (no step opens a truth file): production scripts; truth-free gates (shim parity <= 1e-6 on mHSC-E_500, determinism, secondary-checkpoint load tests, GNInfer setup); scPRINT arm A, arm B, R1; controls; scGPT; K562; secondary checkpoints; QC document and the "step-3 freeze" entry (E13), which is committed before the scoring code opens any truth.
  **E1. Checkpoints (hashes and revisions)**

  - scPRINT `medium-v1.5.ckpt`, HF `jkobject/scPRINT`, commit `8e1ad1a1bc40c1fc2f4f1700fe27923eea9a2bed` (2026-05-05), 221,742,980 bytes, SHA-256 `96ff7c5d9e1ad5076e818bde3c38f498c95cd354cc42029d90691b5ba2964ff8` (unchanged from the 22 Oct entry; re-verified 31 Oct and 1 Nov).
  - **New fact.** The same file name had different content before 2026-05-05.
    - History: uploaded on 2025-03-04 as `v2-medium.ckpt` (commit `f60c8715`), renamed `medium-v1.5.ckpt` on 2025-12-11 (commit `fed0f6bb`), and replaced in commit `8e1ad1a1` ("update medium").
    - Earlier content: revision `d661d779e3614c77305e2fbe0a19312eda2794b8`, 221,775,592 bytes, SHA-256 `a4cf0753270d4ff451a5dbddadd4c59a53e10e4c3e96a8af0db407ad893c36c5`.
    - Comparison tensor by tensor (`torch.load`, 147 state-dict tensors): **identical (max abs difference 0)**; epoch 19, global step 400,000 in both. Only hyperparameter metadata differ (`classes`: a list of names before, a name-to-count dict after; `prenorm`: True before, absent after).
    - So the weights analysed are the same under either hash. The earlier revision is kept in `~/scfm_models/raw/scprint_medium_v15_rev_d661d779/`.
  - Secondary scPRINT files at the same HF commit (not downloaded; hashes from the HF LFS metadata, to verify after download):

    | File | SHA-256 | Bytes |
    |---|---|---|
    | `large-v1.ckpt` | `c0c53989dfc0121cb7502652f643b491e43538e5bd616b9334c504ac6be03ef9` | 1,237,154,305 |
    | `small-v1.ckpt` | `269bc14be52fd0264631facecef400cd95bb1e03f36c23ff400ba02d7c9cd090` | 111,120,529 |
    | `small-v2.ckpt` (scPRINT-2) | `2b586c144cf9a1f638b4b3e803554ebf6ce389f81c5f34179288a970addfd822` | 889,706,878 |

  - scGPT whole-human: hashes as in the 22 Oct entry, re-verified at every smoke run.
    - `best_model.pt` 6cb5d451ab5c4b33eb673adbe4fddc61d2389df1b89b7651a9fe2e557572b922
    - `vocab.json` acca93d114ca62c3f0f50debbd23e8c87f0714f4737764454f6b2b13f2e8580f
    - `args.json` c18e075e018140cb8b2d9029387b9de26607a5ce6a8ccabd6ead70cd76b95d60

  **E2. scGPT checkpoint source pin**

  Source: the "whole-human" link in the README of https://github.com/bowang-lab/scGPT, i.e. Google Drive folder `1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y`. The README lists no other source for this model.

  Drive file IDs, listed from the folder on 2026-11-01 (no download):

  - `best_model.pt` `14AebJfGOUF047Eg40hk57HCtrb0fyDTm`
  - `vocab.json` `1H3E_MJ-Dl36AQV6jLbna2EdvgPaqvqcC`
  - `args.json` `1hh2zGKyWAx3DyovD30GStZ3QlzmSqdk1`

  Drive keeps no version history that we can pin, so the SHA-256 values in E1 are the pin.

  **E3. Package versions and the scPRINT extraction route**

  - Main environment `~/venvs/scfm` (freeze: `docs/scfm_env_freeze.txt`): torch 2.3.1+cu121 (sm_60 binary running on the Pascal cards, verified in the smoke test), scprint 1.0.13, scgpt 0.2.4 (`--no-deps`), lamindb 0.69.9, scdataloader 1.0.1, httpx 0.27.2 and setuptools 80.10.2 (the two step-1 fixes).
  - **Finding.** scprint 1.0.13 and 1.6.1 do not match the checkpoint.
    - They inject a depth token at the encoder input and have a 256x256 decoder layer. The checkpoint has `depth_atinput = False`, no `depth_encoder` weights and a 256x257 decoder layer.
    - The first PyPI release that loads the checkpoint is **scprint 2.2.0**, released 2025-03-06, two days after the checkpoint upload. `scPrint.load_from_checkpoint(ckpt, precpt_gene_emb=None, transformer="normal", strict=True)` loads with no missing or unexpected key; the package's own `on_load_checkpoint` renames `gene_encoder.embedding.weight` and drops the training-only discriminator keys.
    - Versions checked: 2.1.1 fails (hyperparameter parsing); 2.2.0, 2.2.2, 2.3.0 and 2.3.8 load strictly; 2.0.0 was not tested.
    - Environment: `~/venvs/scfm_alt`. scprint 2.2.0 and simpler-flash 1.0.7 are installed with `--no-deps`; torch and the rest are borrowed read-only from `~/venvs/scfm` through a `.pth` file (`docs/scfm_env_freeze_alt.txt`).
  - **Parity.** On the 64 hESC cells of the smoke test, the layer x head attention computed with scprint 2.2.0 equals the step-1 shim result exactly (max abs difference 0.0, relative L1 0.0, all 8 layers).
  - **Entry.**
    - R0, R0' and R2 are extracted with scprint 2.2.0, loaded through the package's own `load_from_checkpoint` (no shim), in `~/venvs/scfm_alt`. The depth token is not fed (`depth_atinput = False`), as in the checkpoint's training.
    - The 1.0.13 shim is retained as a cross-check only. Before production it must reproduce 2.2.0 on mHSC-E_500 (arm A, all cells) to max abs difference ≤ 1e-6; otherwise production stops.
    - `large-v1` and `small-v1` are loaded with the same route, and any deviation is logged.
  - **R1 settings (entry).** R1 = `scprint.tasks.GNInfer` from scprint 2.2.0, after `scdataloader.utils.populate_my_ontology()` in a separate lamin instance. Settings:

    | Group | Settings |
    |---|---|
    | Genes and cells | `how="given"`, `genes` = the covered case genes (Ensembl IDs), `num_genes = max_len = len(genes)`, `max_cells=0` (all cells), `drop_unexpressed=False` |
    | Aggregation | `layer=None` (all layers), `head_agg="mean"`, `preprocess="softmax"` (package default) |
    | Graph output | `filtration="none"` (the default `"thresh"` would replace the ranking by ties), `apc=False`, `symmetrize=False`, `known_grn=None` |
    | Precision | `precision="32"`, `dtype=torch.float32` (fp32 as preregistered) |
    | Other | `forward_mode="none"`, `doplot=False`, `comp_attn` default, a single constant cell-type label so that one network is produced per case |

    The full kwargs dict is logged with each output. Any setting that has to change to make GNInfer run is logged before scoring.

  **E4. scGPT binning seed**

  scGPT's `_digitize` breaks bin-edge ties at random (`side='both'`). Without a seed, two identical runs differed by 0.5% (relative L1). With `np.random.seed(seed)` immediately before the binning of each run, they are bitwise identical.

  Entry: "The numpy global seed is set to **20261021** before binning in every scGPT run. A seed-1 rerun of hESC_1000 is reported as a truth-free stability check. Binning is over the case genes, as in scGPT's own downstream `Preprocessor` (normalize, log1p, gene subset, then bin), whereas pretraining binned over each cell's full gene set."

  **E5. Gene-mapping rule, all 14 cases**

  Fixed rule: `analysis/scfm_step2_genemap.py` (full text in its docstring). Sources, downloaded 2026-11-01:

  | Source | SHA-256 |
  |---|---|
  | Ensembl release 110 GTF, human GRCh38 | `5d3c363835bda4cf8ba025f728b190f69d454e894514bb3c81123efdb15ef1b9` |
  | Ensembl release 110 GTF, mouse GRCm39 | `fd89f56d03871fac5b05a699fa2ea0ce46245d35b64248cf74d9f70b59b61618` |
  | HGNC complete set | `d93ce6b22d01dc1e1ded28857246951d605ead47ed17ef973ad1142b1d13bed0` |
  | MGI `MRK_List2.rpt` | `00c4dd2abd9e2fe70b7810ee2914a6af65015dfd620fd493de2ce7f811cc915f` |
  | MGI `MRK_ENSEMBL.rpt` | `1c8ee6e696bc30dbe3eea7d73e58a1fe9c6ec38674db8ddf2f762babcae096eb` |

  1. Primary assembly = chromosomes 1-22 (human) or 1-19 (mouse), X, Y, MT. Non-primary (patch, haplotype, scaffold) IDs are never chosen, even though the scPRINT vocabulary contains 3,141 human and 51 mouse such IDs.
  2. Exact symbol match among primary genes (human: case-sensitive; mouse: case-insensitive, because the mDC and mHSC files use upper-case mouse symbols). If several primary genes share the name, tie-break in this order:
     - keep those in the scPRINT vocabulary;
     - then the HGNC/MGI-curated Ensembl ID;
     - then protein_coding;
     - then the longest genomic span;
     - if still not unique, drop the gene (ambiguous).
  3. Otherwise, the curated HGNC/MGI ID of an approved symbol.
  4. Otherwise, a unique HGNC previous symbol, then a unique HGNC alias / MGI synonym, resolved to the current symbol and mapped by steps 2-3.
     - Ambiguous synonyms are dropped.
     - ERCC spike-ins are dropped.
     - Excel-date symbols are mapped (none occurred).
  5. The ID must be in the scPRINT checkpoint vocabulary of that organism (23,149 human and 21,607 mouse IDs; the checkpoint has one joint vocabulary). If two symbols resolve to one ID, the exact-match symbol is kept; otherwise both are dropped.
  6. scGPT (human cases): the case symbol if it is in the 60,697-symbol vocabulary, else the resolved current symbol. No ortholog mapping for mouse.

  **Result.** Outputs: `results/scfm_step2/coverage_by_case.csv`, per-gene tables in `results/scfm_step2/gene_map/`, and dropped genes with reasons in `dropped_genes.csv`.

  - scPRINT covers 85.5% to 99.6% of the case genes. The lowest are mESC_1000 (1,306 of 1,528) and mHSC-L_1000 (647 of 692); the loss is non-coding or predicted genes missing from the vocabulary.
  - TFs: 2 of 321 lost in mDC (TSC22D4 absent from Ensembl 110, MLL2 ambiguous); all others covered.
  - scGPT intersection: hESC 841 / 1,285 genes and hHep 905 / 1,360 genes (500 / 1000 sets).
  - Non-exact routes: 3.3% of mapped genes. The 40 distinct TFs among them were checked by name, and all are standard renamings.
  - RSEM expected counts (hESC, arm B) are used without rounding.
  - Dropped truth edges per case are not reported yet. That needs the truth files, which stay closed until scoring; the per-gene mapping tables allow the count to be computed then.

  **E6. Input-scale arms**

  The preregistered input ("BEELINE's documented transform, back-transformed") is `2^y - 1` of the BEELINE ExpressionData. The transform check (all 7 datasets, `results/scfm_step2/counts_alignment.json`) shows that this is not TPM in general. BEELINE's `y = log2(1 + s_c * x)`, with a per-cell scale `s_c` and `x` the linear value of the GEO deposit. This holds to 1e-3 relative precision in 99.98-100% of entries (mESC 98.8%). Step 1's description of the input as TPM is corrected accordingly.

  | Dataset | GEO linear scale x | s_c median (range) |
  |---|---|---|
  | hESC | RSEM expected counts (not UMI) | 0.86 (0.34 to 13.5) |
  | hHep | log2 values only; BEELINE equals them exactly | 1 (identity) |
  | mDC | RSEM TPM | 0.020 (0.0006 to 0.035) |
  | mESC | sailfish transcript TPM summed to gene | 0.097 (0.065 to 0.18) |
  | mHSC (E, GM, L) | HTSeq integer counts | 1.02 to 1.25 (0.07 to 7.2) |

  Entry:

  - **Arm A (primary, all models, 14 cases):** `2^y - 1` of the BEELINE matrix.
  - **Arm B (scPRINT only, descriptive sensitivity):** GEO counts, where the deposit holds them: hESC and mHSC-E/GM/L, 8 cases, n = 4 datasets.
    - Reported: per-dataset d_i under B, the paired Δ_i (scPRINT R0 under B minus under A, log2 AUPRC ratio), and the HL estimate with its range.
    - No p-value is quoted.
    - If the median Δ_i exceeds δ = 0.15, the paper states that the primary result depends on the input scale on the count datasets.
  - **No arm B for scGPT.** Its preprocessing is invariant to a per-cell factor, so arm B equals arm A on hESC up to bin-edge ties.
  - **Arm C (GEO linear) is not run.**
  - **No re-quantification** from SRA is done before scoring (PI decision, 2026-11-01). Should it ever be ordered later, its pipeline is committed before scoring and the primary stays arm A.
  - The organism-head accuracy under A and B is reported as a truth-free check that the input is in distribution.

  **E7. Sink-head reporting**

  The step-1 smoke test found heads that put 90-100% of a gene query's attention on the 9 special tokens. R0 remains the mean over all heads and layers, as preregistered.

  Entry:

  - For every case, record the cell-averaged special-token mass and the relative row entropy per layer x head next to the readout.
  - **R0'** (descriptive) = the mean over heads with special-token mass < 0.5. If fewer than 4 of the 32 heads pass, use the 4 heads with the lowest mass.
  - The threshold and the minimum are fixed here, from the smoke-test distribution, and use no truth.
  - Each case's head set is computed in step 3 from expression only and committed before any truth file is opened.
  - R0' enters the robust-NI maximum (E11).
  - scGPT has no R0'; its `<cls>` mass is < 0.02 in the last layer.

  **E8. Data notes (no effect on the design)**

  - All 14 case matrices equal the corresponding block of the BEELINE ExpressionData (all cells, case genes).
  - mHSC: 67 (E), 45 (GM) and 17 (L) BEELINE cell labels do not carry the expression of the GEO column with the same name. They were re-matched by content to a different GEO column (identical zero pattern and a constant ratio, 100% agreement). Counts are therefore matched for every cell of every dataset.
  - Leakage audit (metadata only): see E9.

  **E9. Leakage audit outcome, done before scoring**

  **Corpora searched.**

  - scPRINT v1: CELLxGENE census LTS 2023-12-15 (Kalfon et al. 2025, Data availability).
  - scGPT whole-human: "33 million normal human cells" from CELLxGENE. The snapshot is not stated in any reachable source, so all three LTS releases (2023-05-15, 2023-07-25, 2023-12-15) were checked.
  - scPRINT-2: scBaseCount.

  **Result.**

  - None of the six accessions (GSE75748, GSE81252, GSE48968, GSE98664, GSE81682, Replogle PRJNA831566) is cited by any CELLxGENE public collection. Searched by accession, BioProject, DOI, and first author and year.
  - None has an SRX in scBaseCount 2025-02-25.
  - Replogle K562 has 6,237 of 7,883 SRX in scBaseCount 2026-01-12.

  **Classification.**

  - **Out (moderate confidence)** for every scPRINT v1 and scGPT pair.
  - **Undeterminable** for K562 x scPRINT-2, because the scBaseCount snapshot that scPRINT-2 used is unknown. It is flagged in the optional, descriptive scPRINT-2 arm.
  - Residual limits: re-use of the same cells inside an atlas that does not cite the study would go undetected. 12 census datasets without a current public collection were identified by title, and none is related.
  - Table: `results/scfm_step2/leakage_audit.csv`.
  - Under this classification, the preregistered leaked-dataset sensitivity analysis is empty for scPRINT v1 and scGPT.

  **E10. Cells, batches and control seeds**

  **Cells and batches.**

  - All cells of every BEELINE dataset, as preregistered: 4,794 per gene set.
  - Batch size 8. All case genes are tokens, so there is no padding.
  - Truth-free split-half stability per case: cells split once with seed 20261021; relL1 and Spearman between the half-means are reported in the step-3 QC.

  **Controls.**

  - Random weights: seeds 1-3, two versions (ESM2 gene embeddings kept; gene embeddings random).
  - Expression permuted per gene across cells: seeds 1001-1003.
  - K562: 5 subsamples of 1,024 non-targeting control cells (seeds 1-5), as preregistered.
  - Fine-tuning arm (descriptive): package default losses, 10 epochs, seed 20261021, no truths.

  **E11. Robust-NI set**

  Robust NI is called only if L > -δ also holds when scPRINT takes, per dataset, the best of {R0, R0', R1, R2} x {S = Aᵀ, A, max(A, Aᵀ)} x {arm A, arm B where available}. That maximum is chosen on the truths, so it is optimistic for scPRINT. This extends the preregistered maximum ({R0, R1, R2} x 3 orientations), which can only make robust NI harder to reach. It is the pre-specified fair-to-scFM reading. It is used only to stress our claim, never to claim a scPRINT win.

  **E12. Descriptive arms and per-contrast gene sets**

  **Gene sets.** Each contrast is scored on its own G∩:

  - CLR-NS, CLR or |r| vs scPRINT: genes covered by scPRINT.
  - vs scGPT: genes covered by scGPT.
  - scPRINT vs scGPT: the three-way intersection.

  **Descriptive arms.**

  - **scGPT** (hESC and hHep, plus K562): per-dataset d_i, plus the preregistered K562 bootstrap. No p-value and no NI claim.
  - **Arm B:** as in E6.
  - **CLR-NS on the full gene set** is reported next to its G∩ value.

  The primary test, the margin, the interval, the decision table and the key secondaries are unchanged.
- 2026-11-02: step-3 freeze (E13). Production extraction and its truth-free QC are finished; **no truth file has been opened and nothing has been scored**. The extraction code is `analysis/scfm_extract_{scprint,scgpt,r1}.py` (commits 84586ef to 45157e3); the QC document is `docs/scfm_step3.md`, the QC tables and the SHA-256 manifest are in `results/scfm_step3/`. Entry committed after the QC and before the scoring code opens any truth.
  **Files to be scored (SHA-256).** `results/scfm_step3/sha256_manifest.csv` lists the size and SHA-256 of all 239 network and attention files (scPRINT medium-v1.5 arm A and B, R1 arm A and B, controls, secondary checkpoints, scGPT, ESM2 cosines) under `~/scfm_models/step3/`; the SHA-256 of that manifest is `e27a8422f5cf868459c694407a0cd0e900a18ce8b3aaf0c33ee179260c2e3f15`. The manifest is the frozen list; a file whose hash differs from it is not to be scored.
  **R0' head sets (from expression only; layer L, head H, 1-based; rule: cell-averaged special-token mass < 0.5, 4 lowest-mass heads if fewer than 4 pass, which happened in no case).**

    | Case | Heads | Head set |
    |---|---|---|
    | hESC_500 | 13 | L1H1;L1H2;L1H3;L2H3;L3H2;L3H3;L4H1;L4H2;L5H1;L7H1;L7H3;L8H1;L8H2 |
    | hESC_1000 | 13 | L1H1;L1H2;L1H3;L2H3;L3H2;L3H3;L4H1;L4H2;L5H1;L7H1;L7H3;L8H1;L8H2 |
    | hHep_500 | 12 | L1H1;L1H2;L1H3;L2H3;L3H2;L4H1;L4H2;L5H1;L7H1;L7H3;L8H1;L8H2 |
    | hHep_1000 | 13 | L1H1;L1H2;L1H3;L2H3;L3H2;L4H1;L4H2;L5H1;L6H3;L7H1;L7H3;L8H1;L8H2 |
    | mDC_500 | 14 | L1H1;L1H2;L1H3;L2H3;L3H2;L3H4;L4H1;L4H2;L5H1;L6H3;L7H1;L7H3;L8H1;L8H2 |
    | mDC_1000 | 13 | L1H1;L1H2;L1H3;L2H3;L3H2;L4H1;L4H2;L5H1;L6H3;L7H1;L7H3;L8H1;L8H2 |
    | mESC_500 | 14 | L1H1;L1H2;L1H3;L2H3;L2H4;L3H2;L4H1;L4H2;L5H1;L5H2;L7H1;L7H3;L8H1;L8H2 |
    | mESC_1000 | 15 | L1H1;L1H2;L1H3;L2H1;L2H3;L2H4;L3H2;L4H1;L4H2;L5H1;L5H2;L7H1;L7H3;L8H1;L8H2 |
    | mHSC-E_500 | 14 | L1H1;L1H2;L1H3;L2H3;L2H4;L3H2;L3H3;L4H1;L4H2;L5H1;L7H1;L7H3;L8H1;L8H2 |
    | mHSC-E_1000 | 14 | L1H1;L1H2;L1H3;L2H3;L2H4;L3H2;L3H3;L4H1;L4H2;L5H1;L7H1;L7H3;L8H1;L8H2 |
    | mHSC-GM_500 | 13 | L1H1;L1H2;L1H3;L2H3;L2H4;L3H2;L4H1;L4H2;L5H1;L7H1;L7H3;L8H1;L8H2 |
    | mHSC-GM_1000 | 13 | L1H1;L1H2;L1H3;L2H3;L2H4;L3H2;L4H1;L4H2;L5H1;L7H1;L7H3;L8H1;L8H2 |
    | mHSC-L_500 | 14 | L1H1;L1H2;L1H3;L2H3;L3H2;L3H3;L4H1;L4H2;L5H1;L7H1;L7H3;L8H1;L8H2;L8H4 |
    | mHSC-L_1000 | 14 | L1H1;L1H2;L1H3;L2H3;L3H2;L3H3;L4H1;L4H2;L5H1;L7H1;L7H3;L8H1;L8H2;L8H4 |

  The 12 heads L1H1, L1H2, L1H3, L2H3, L3H2, L4H1, L4H2, L5H1, L7H1, L7H3, L8H1 and L8H2 are in every set. Special-token mass and entropy of every layer x head and case: `results/scfm_step3/sink_mass_A.csv`. R0 stays the all-head, all-layer mean (primary); R0' and R2 are derived from the stored per-head tensors.
  **R1 kwargs actually used** (scprint 2.2.0 `GNInfer`, identical for all cases; `genes` = the covered case genes as Ensembl ids): how="given", num_genes = max_len = len(genes), max_cells=0, drop_unexpressed=False, layer=None, head_agg="mean", preprocess="softmax", filtration="none", apc=False, symmetrize=False, known_grn=None, precision="32", dtype=torch.float32, forward_mode="none", doplot=False, comp_attn default, cell_type_col="cell_type" with one constant label, batch_size=8, num_workers=0 (the last two were not in E3). Settings that had to be arranged to make GNInfer run, not in E3: `populate_my_ontology()` raised a TypeError at its mouse developmental-stage step in scdataloader 1.0.1; the steps before it completed, the mouse-stage step was skipped, and the Disease and Gene steps were run by `analysis/scfm_r1_setup.py` in the separate lamin instance (Gene registry filled in chunks of 400 with the 44,756 ids of the scPRINT vocabulary); GNInfer sends batches to the bare device "cuda", so the R1 runs used `CUDA_VISIBLE_DEVICES=1`. R1 row sums are 0.25 to 0.62 (the package strips the special tokens without renormalising); this is the package rule and is not altered.
  **Other choices fixed in step 3 that E1-E12 did not cover.** (a) Arm B uses the case genes that have deposited counts: one gene fewer than arm A in mHSC-E_1000 and mHSC-GM_1000. (b) large-v1 and small-v1 have `depth_atinput` True; the depth token is fed as the per-cell sum of the model input over the case genes (the Collator rule); they are descriptive, mean readout only. (c) Random-weight controls use the package's native initialisation (`scPrint(**checkpoint hyper-parameters)`, `torch.manual_seed(seed)`); in the "random gene embeddings" version the embeddings keep `nn.Embedding`'s N(0,1) scale (std 1.0 against 0.069 pretrained). (d) The npz files hold string object arrays (`genes`, `ensembl`) and need `allow_pickle=True`.
  **QC summary.** Shim (scprint 1.0.13) vs scprint 2.2.0 on mHSC-E_500, arm A, all 1,071 cells: max abs difference 0.0 (bitwise), gate 1e-6 passed. Unit check: the 64-cell hESC_1000 step-1 arrays are reproduced bitwise by both extractors. Determinism: scPRINT and scGPT hESC_500 twice, bitwise identical. No NaN or Inf in any output; attention row sums within 2.4e-7 of 1. Split-half relative L1 of R0: 0.08 to 0.17 (14 cases). Permutation control (seeds 1001-1003) moves R0 by relative L1 0.07 to 0.33 (Spearman 0.91 to 0.99), 0.8 to 2.1 times the split-half distance; random-weight networks are near uniform. Arm A vs arm B R0: relative L1 0.04 to 0.09, Spearman 0.988 to 0.996. scGPT `<cls>` mass at most 0.0012 in the last layer; seed-1 rerun of hESC_1000 differs from seed 20261021 by relative L1 0.0023 (last-layer head mean). Organism-head accuracy 0.99 to 1.00 in all cases except mHSC-E_500 (0.67 under arm A, 0.69 under arm B), a flagged observation with no gate attached; the secondary checkpoints classify every case correctly.
  **Not in this freeze.** (i) The K562 stratum (step 3.8) has not been extracted: its preregistered gene set needs CollecTRI regulon sizes, which cannot be read before scoring under the no-truth rule; it will be frozen in a separate entry, with its gene list and file hashes, before any K562 truth is opened. (ii) The fine-tuning arm (step 3.10) has not been run; its settings stay as in E10 and it will be logged separately. (iii) No SRA re-quantification (PI decision, 2026-11-01).
- 2026-11-02: K562 gene-set definition (E14). The PI cleared scoring (GO, 2 Nov 2026) and allowed CollecTRI to be read **only to count targets per TF** for the K562 gene set. This entry records exactly what was read and the gene set that results; it is committed **before the K562 extraction**. **No K562 truth has been opened or scored**: no DoRothEA file, no knockdown-response (bulk) file, no CollecTRI target name, no network.
  **What was read (script `analysis/scfm_k562_prep.py`, scheduler job 1012, 2026-11-02 11:07).**
  - `data/human_truth/collectri.tsv` (SHA-256 `98ef1163146d2e28ee413e7c909174998cd9d742ea93d97ccfe93f3a14c160c1`, 64,516 rows): only the columns `source_genesymbol` and `target_genesymbol` were loaded, de-duplicated and aggregated to the number of distinct targets per source symbol (a regulon size). 1,201 source symbols, of which 587 have at least 10 targets (source symbols that are not a single detected deposit gene drop out in the next step). No target name was inspected, printed or stored. (A `head` of the file's first three lines had been displayed once, by accident, while checking the column layout, before the script was written; it carried no information used here.)
  - `data/perturbseq/K562_essential_raw_singlecell_01.h5ad` (10,661,879,995 bytes, SHA-256 `3e5a63a9e892b21029bb55fca4e12517a49aad7af6c14133ca63d12cf68c6cee`): obs `gene` (to find the non-targeting cells), obs `cell_barcode`, var `gene_name`/`gene_id`, and the expression rows of the **10,691 non-targeting control cells only**; no perturbed cell's expression was read.
  **Rules fixed here (the preregistration leaves these open; none uses a truth other than the regulon sizes).**
  - *Detected*: nonzero UMI count in at least 5% of all 10,691 control cells (pooled, not per subsample). 8,561 of the 8,563 matrix genes pass.
  - *TF part*: CollecTRI source symbols with at least 10 distinct targets (counted over the whole CollecTRI table, not restricted to the gene set), whose symbol equals a deposit `gene_name` exactly, that are detected, and whose symbol is unique in the matrix (4 matrix rows carry a duplicated symbol and are excluded). Result: **278 TFs** (regulon size 10 to 1,333, median 40).
  - *Variable part*: the **1,000** detected genes (unique symbols) with the largest variance of log1p(count / cell total x median cell total) over all control cells (cell total = sum over the 8,563 matrix genes, as in `analysis/perturbseq.R`; median total 13,714 UMI).
  - *Gene set*: union of the two parts (overlap 34), **1,244 genes**, deposit order. File `results/scfm_k562/geneset.tsv` (SHA-256 `22caabb6e23868cf495f55f4c8e33c95a029c1632be74decf4d1e3ba1cc0a550`) with symbol, deposit Ensembl id, detection fraction, variance, regulon size and part flags.
  - *Subsamples*: k = 1..5, `numpy.random.default_rng(k).choice(<control cell indices in file order>, 1024, replace=False)`, sorted. File `results/scfm_k562/subsamples.tsv` (SHA-256 `b4d6047b62e95c7fbde8cdc50b20d7a7a5011605b543fa72edc51fe57fd4dd1a`) lists the barcodes and file rows.
  - *Gene mapping*: the E5 rule (human). scPRINT covers **1,218 of 1,244** (all 278 TFs); scGPT 1,240; the intersection 1,217. The 26 scPRINT drops are 21 genes not in the checkpoint vocabulary (mostly lncRNAs such as SNHG1, NEAT1, MALAT1, XIST) and 5 with no match. For 1,217 of the 1,218 mapped genes the E5 Ensembl id equals the deposit's own `gene_id`. Gene maps `results/scfm_step2/gene_map/hK562_s1..5.tsv` (identical; SHA-256 of s1 `0d30d852f440e5f566d95b2de861195ef7147bea11c89c1387e7bf2695717079`).
  - *Model input*: raw UMI counts of the subsample (each model applies its own normalization); extractor case names `hK562_s1..5` (the prefix `h` selects the human organism in the unchanged extractors); the comparators (|r|, CLR, CLR-NS) are computed in the scoring step on log1p(count / cell total x 13,714) of the same cells and the same G-intersection.
  **K562 extraction plan (GPU device 1, one job at a time, scheduler `# after:` chain).** scPRINT medium-v1.5 R0, R0', R2 (arm A) and R1 (E3 settings) on the five subsamples; permutation control (seeds 1001-1003) and random-weight controls (seeds 1-3, two versions), R0 only; large-v1 and small-v1, R0 only; scGPT (attention, token-embedding and hidden-state cosine, binning seed 20261021). No ESM2-cosine prior for K562. The extraction is frozen in a separate entry (SHA-256 manifest) before any K562 truth is opened.
  **K562 scoring definition (fixed here, before any K562 truth is opened).** Truths: DoRothEA A-B and CollecTRI regulons as TF -> target edge lists, and the knockdown-response modules of `response_sets()` in `analysis/human_truths.R` (defaults) as knocked-down gene -> module-gene edge lists, each restricted to the scored gene set; metric: directed TF-node AUPRC ratio and tie-aware EPR (`score_truth` of `analysis/beeline_baselines.R`), as for BEELINE. Per truth and subsample, d = log2(AUPRC ratio of CLR-NS / AUPRC ratio of the model readout) on the model's G-intersection; reported: the mean over the five subsamples and a 95% percentile bootstrap over regulators (B = 1,000, regulators resampled with replacement, the same resample for all five subsamples and both methods), per truth and as the mean over the three truths. Descriptive only: no p-value, no NI claim.
- 2026-11-02: K562 extraction freeze (E15). The K562 extraction (E14) and its truth-free QC are finished; **no K562 truth file has been opened and no K562 network has been scored**. Committed before the K562 scoring code opens any K562 truth.
  **Files to be scored (SHA-256).** `results/scfm_k562/sha256_manifest.csv` lists size and SHA-256 of all 75 network and attention files under `~/scfm_models/step3/k562/` (scPRINT medium-v1.5 R0/R0'/R2 and the per-head tensors, R1, scGPT, permutation seeds 1001-1003, random weights seeds 1-3 in two versions, large-v1, small-v1; five subsamples each; paths relative to that directory). The SHA-256 of the manifest is `787cbf8749994df2f21946667cd28acd450e333b8ecebf161f53f70d5976080b`. A file whose hash differs is not to be scored; the exporter (`analysis/scfm_score_export.py --k562 --manifest-sha ...`) stops on any mismatch.
  **Extraction.** Scheduler jobs 1013 to 1017 (chained with `# after:`, GPU device 1 only, one at a time), the unchanged extractors of step 3 (`analysis/scfm_extract_{scprint,r1,scgpt}.py`, scprint 2.2.0 in `~/venvs/scfm_alt`, scGPT 0.2.4 in `~/venvs/scfm`), case names `hK562_s1..5`, 1,024 cells and 1,218 (scPRINT) or 1,240 (scGPT) genes per case, input = raw UMI counts, binning seed 20261021. R1 kwargs as in E13. The ESM2-cosine prior was not extracted for K562 (E14).
  **QC (truth-free).** No NaN or Inf and attention row sums within 1e-4 of 1 in every file (the job gate; all 14 output directories x 5 cases passed; wording corrected after commit 8f85558). scPRINT organism-head accuracy 1.00 in all five subsamples and for large-v1 and small-v1. Split-half relative L1 of R0: 0.077 to 0.088 (medium), 0.050 to 0.054 (large), 0.050 to 0.057 (small). R0' head sets: 13 heads for medium-v1.5 (same rule as E7). scGPT `<cls>` mass in the last layer at most 0.0007. R1 row sums 0.28 to 0.57 (the package rule, as in E13). Two details: the K562 scPRINT runs took 58 s per case (hESC_1000 in step 3: 11 s), cause not investigated; disk of `~/scfm_models` is 13 GB of the 20 GB budget.
  **Not in this freeze.** The fine-tuning arm (step 3.10) was not run. The scoring of K562 follows the definition fixed in E14.
- 2026-11-02: scoring done (E16). Scoring code committed before any real truth was opened (ee9de39; smoke test on synthetic truth only). BEELINE scoring run once on the real truths (job 1020, 08fa10b) after re-verifying the E13 manifest (0 mismatches); K562 scoring run once (job 1023) after re-verifying the E15 manifest. Primary: HL +1.032 log2 units, exact interval [+0.343, +1.897], p = 0.00781; rows of the decision table that apply: NI and superior; robust NI holds (interval [+0.233, +1.656]); S1 and S2 both pass Holm. Full results, deviations and caveats: `docs/scfm_scoring.md`. One post-run correction: a gene-set filter missing from some descriptive rows of `analysis/scfm_score_stats.R` was added and the statistics stage rerun on the same metric table (changes at most 1.3e-4 in one descriptive row; no decision changed). Interpretations of undefined details (best configuration in the robust maximum, hybrid CLR on the symmetrised R0, `nlme` for the mixed model, baselines scored on the G-intersection without refitting) are listed in section 9 of that document.
