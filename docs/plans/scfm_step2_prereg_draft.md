# Draft change-log entries for `docs/preregistration_scfm.md` (step 2, 2026-11-01)

**Status: DRAFT, not part of the preregistration.** `docs/preregistration_scfm.md` has not been edited. The eight open decisions of `docs/scfm_step2.md` are settled below. Each one is **proposed, PI to confirm**. Once the PI confirms (or edits) them, the entries E1 to E13 are committed as a single dated change-log entry before any production network is computed. Everything here was prepared without opening any truth file and without scoring any network. The evidence is in `docs/scfm_step2.md`, `docs/scfm_smoke_test.md` and `results/scfm_step2/`.

## Plain-language summary (for a big-picture reader)

1. We test whether CLR attention, a 2007 formula with no training, finds gene-regulatory links as well as the attention inside single-cell foundation models trained on tens of millions of cells (scPRINT is the main comparator, scGPT the second).
2. The test was fixed in advance: 7 standard BEELINE datasets, one yardstick (accuracy against known regulator-target links), one confidence interval, and a margin of about 11%.
3. Every choice that could tilt the result (model code version, gene naming, input scale, which attention heads count) is written down before anyone looks at the answers.
4. If CLR attention comes out non-inferior, an untrained, interpretable operator matches learned attention on this task. That would suggest the regulatory signal these models carry in their attention is largely available from the correlation structure that CLR computes directly.
5. If the foundation models win clearly, the paper says so and gives the size of the gap. CLR attention then becomes the transparent baseline these models have to beat. That is still a clean, publishable result.
6. Either way the comparison is fair to the models. A robustness check gives them their best readout, raw counts where these exist, and removal of noise heads. So a CLR win cannot be blamed on handicapping them, and an unclear result is reported as unclear, never as a tie.

## Decisions (proposed, PI to confirm)

| # | Decision | PI judgment call? | Committed before scoring |
|---|---|---|---|
| 1 | Input scale. Primary stays **arm A** (BEELINE scale, n = 7). **Arm B** (GEO counts, n = 4) is a pre-specified descriptive sensitivity analysis. Arm B also enters the robust-NI maximum (the fair-to-scFM reading). Arm C is dropped. No SRA re-quantification before scoring. | **Yes**: whether to commission re-quantification of hHep, mDC and mESC from SRA (and, if so, whether that blocks scoring) | E6, E11: arm definitions, the n = 4 reporting rule, robust-NI set |
| 2 | All cells per case, as preregistered. No subsampling. A split-half stability check (truth-free) is reported. | no | E10: all cells, batch size, split seed, control seeds |
| 3 | The primary is unchanged: one exact Wilcoxon NI interval on n = 7, δ = 0.15. scGPT (n = 2 plus K562) and arm B (n = 4) are descriptive only: per-dataset values and Hodges-Lehmann (HL) estimates, with no p-value and no NI claim. | confirm only | E12: descriptive reporting rules, per-contrast G∩ |
| 4 | Production extraction uses **scprint 2.2.0 `load_from_checkpoint`** (no shim) in `~/venvs/scfm_alt`. The 1.0.13 shim is a cross-check, with a parity gate on one mouse case. Both checkpoint hashes are recorded. R1 runs through 2.2.0 `GNInfer` with the settings in E3. | no | E1, E3: route, versions, parity tolerance, R1 settings |
| 5 | Special-token mass and entropy are reported for every layer x head. Add the descriptive **R0'** (heads with mass < 0.5, at least 4 heads). R0' also enters the robust-NI maximum. R0 stays primary. | **Yes**: whether R0' enters the robust-NI maximum (recommended) | E7, E11: threshold, minimum head count. Per-case head sets go in the step-3 entry |
| 6 | The scGPT binning seed is **20261021** (the project production seed), set before binning in every run. A seed-1 rerun on hESC_1000 is a truth-free stability check. scGPT is run (descriptive). | no | E4 |
| 7 | Accept the E5 mapping rule as it stands: alias/previous routes, the mESC loss, and non-integer RSEM counts used without rounding. | no | E5 |
| 8 | Accept "out, moderate confidence" for every scPRINT v1 and scGPT pair. K562 x scPRINT-2 becomes **undeterminable**, not "in". The leaked-dataset sensitivity analysis is empty for the two scored models. | optional: email the scGPT authors for the census snapshot (not blocking) | E9 |

### D1. Input-scale arm

**Decision.** The primary contrast stays on arm A, the preregistered input: `2^y - 1` of the BEELINE matrix, n = 7. Arm A is the only arm that keeps n = 7, and n = 7 is what makes the exact NI interval attainable. Switching the primary input now would also depart from the preregistered text. Arm B (GEO counts; hESC and mHSC-E/GM/L, n = 4) is run for scPRINT as a pre-specified descriptive sensitivity analysis. It is also added to the robust-NI maximum (E11), so the fair-to-scFM reading gives scPRINT the better of A and B wherever counts exist. Arm C (GEO linear: TPM or log2-derived values) is dropped. For hHep it equals A, and for mDC and mESC it differs from A only by a per-cell scale, so it is neither counts nor the preregistered input.

**Why arm A is a mild handicap at most, and where it is not one at all.**

- **scGPT.** Its own preprocessing (`normalize_total`, `log1p`, then per-cell quantile binning; `scgpt/preprocess.py` steps 3 to 6) is invariant to any positive per-cell factor. On hESC, arm A is counts x s_c with a cell-constant s_c (const-ratio 1.0000), so arms A and B give the same tokens up to bin-edge tie-breaking. scGPT therefore needs no counts arm.
- **scPRINT.** The encoder sees `log2(1 + x)` and gets no depth token (`depth_atinput = False`). On the 4 count datasets, arm A differs from counts only by the per-cell factor s_c (median 0.86 to 1.25). On mDC and mESC, arm A is TPM divided by about 50 and 10.
- **Arm B is not the pretraining distribution either.** Its counts are plate-based read counts (RSEM expected counts, HTSeq), not UMIs. The scPRINT assay head classifies the hESC input as Smart-seq2 (EFO:0008931, smoke test), consistent with plate-based data being present in the census.

**What a fair-to-scFM reading needs, and what is committed.**

1. Robust NI (E11): per dataset, scPRINT takes the best of {R0, R0', R1, R2} x {S = Aᵀ, A, max(A, Aᵀ)} x {arm A, arm B where available}, chosen on the truths. NI is called robust only if L > -δ still holds.
2. A paired, descriptive input-sensitivity readout on n = 4: Δ_i = log2 AUPRC ratio of scPRINT R0 under B minus under A, plus d_i under B. Report each value, the HL estimate and the range.
3. A pre-stated reading rule. If the median Δ_i exceeds δ (counts help scPRINT by more than the NI margin), the paper states that the primary result is input-dependent on the count datasets. The robust-NI verdict then governs the wording.
4. A truth-free check that the input is in distribution: the organism-head accuracy per case under A and B (step-3 QC).
5. An explicit statement that hHep, mDC and mESC have no deposited counts.

**Rejected alternative.** A hybrid primary (counts where available, A elsewhere) mixes input definitions across the 7 blocks. On the count datasets it differs from A only by s_c, and its information is already in the robust-NI maximum.

**PI judgment call.** Re-quantify hHep, mDC and mESC from SRA FASTQ, so that arm B covers 7 of 7 and becomes a powered sensitivity analysis with the same interval? Compute is free, but this means hundreds of GB of FASTQ (mDC alone has 2,425 SRX, of which 383 cells are used) and our own alignment choices. Recommendation: do not block scoring on it. If the PI wants it, its pipeline must be committed before scoring, and the primary stays arm A whatever it shows.

### D2. Cells per case

**Decision.** All cells (4,794 per gene set; 9,588 cell-runs over the 14 cases), as preregistered. Step 1 showed the 64-cell mean has not converged (64 vs 256 cells: relL1 0.20), and all cells cost about 5 GPU-minutes in total for scPRINT. A seeded subsample would add a choice and buy nothing.

- Batch size 8. All case genes are tokens, so there is no padding.
- Truth-free diagnostic: a split-half stability check per case (cells split once with seed 20261021; relL1 and Spearman between the two half-means).
- Controls as preregistered: random weights with 3 seeds (1 to 3) in 2 versions. The permutation control uses 3 seeds (1001 to 1003), because the preregistration does not fix the number and step 1 recommended several.
- K562 stays at the preregistered 5 seeded subsamples of 1,024 control cells.

### D3. Power and primary test structure

**Decision.** Keep the preregistered primary unchanged: one exact Wilcoxon/HL interval on n = 7 (4th and 25th Walsh averages, 92.2% coverage), δ = 0.15, the fixed sequence NI then superiority, and the decision table as written. Power by the preregistered simulation is 0.69 to 0.97. All 7 datasets run, so the drop rule is not triggered. Changing the primary statistic after CLR-NS's scores are known would cost more credibility than it could gain in power.

- **Key secondaries S1 and S2:** Holm over two tests at 0.025 stays attainable (p_min 0.0078).
- **scGPT** (hESC and hHep, n = 2; plus K562): descriptive only. Report the per-dataset d_i and the K562 bootstrap interval. No p-value is quoted for n = 2 (p_min = 1/4) and no NI claim is made.
- **Arm B** (n = 4): descriptive only. Report the HL estimate and the range (p_min = 1/16 is not quoted as evidence).

**Per-contrast gene sets (E12).** Each contrast is scored on its own G∩:

- CLR-NS/CLR/|r| vs scPRINT: genes covered by scPRINT.
- vs scGPT: genes covered by scGPT.
- scPRINT vs scGPT: the three-way intersection (within 4 genes of the scPRINT set).

### D4. scPRINT extraction route

**Decision.** The production route is scprint 2.2.0, `scPrint.load_from_checkpoint(ckpt, precpt_gene_emb=None, transformer="normal", strict=True)`, in `~/venvs/scfm_alt`. It is the released code path, loads with no missing or unexpected key, and reproduced the step-1 shim bit for bit (max abs difference 0, all 8 layers). A referee can rerun it without our patches.

**Gate (step 3).** The shim is a cross-check. Before production it is run on one mouse case (mHSC-E_500, arm A, all cells); max abs difference ≤ 1e-6 is required, otherwise stop and report.

**Checkpoint hashes (E1).** Record both: the pinned `96ff7c5d...4ff8` (used) and the earlier revision `a4cf0753...36c5`, whose weights are identical.

**R1.** The robust-NI maximum needs R1, so R1 runs through 2.2.0 `GNInfer` after `populate_my_ontology()` in a separate lamin instance, with the settings in E3. If R1 cannot be made to run, it is reported as not run and dropped from the robust-NI maximum. That makes robust NI less stringent, and this is stated.

**Secondary checkpoints.** `large-v1` and `small-v1` get a load test under 2.2.0; any other version they need is logged.

### D5. Sink heads

**Decision.** R0 stays the preregistered all-head, all-layer mean. For every case, report the special-token mass and the relative row entropy per layer x head. Add **R0'** as a descriptive readout: the mean over heads whose cell-averaged special-token mass is < 0.5. If fewer than 4 of the 32 heads pass, use the 4 heads with the lowest mass.

The threshold and the minimum are fixed now from the step-1 distribution. Each case's head set is computed from expression data only in step 3 and committed (step-3 entry) before any truth is opened. R0' answers the obvious referee point that averaging heads with 90 to 100% sink mass dilutes scPRINT's signal.

**PI judgment call.** Should R0' enter the robust-NI maximum? Recommended: yes. Adding R0' can only make our robust-NI claim harder to reach. scGPT needs no R0', because its `<cls>` mass is < 0.02 in the last layer it is read from.

### D6. scGPT binning seed and whether scGPT runs

**Decision.**

- **Seed.** `np.random.seed(20261021)` immediately before binning in every scGPT run. This is the project production seed (also used for GRNBoost2). Any fixed value is equally valid as long as it is fixed before scoring.
- **Stability check.** A truth-free rerun of hESC_1000 with seed 1 is reported as relL1 (step 1: 0.5%).
- **Binning.** Binning is over the case genes. This matches scGPT's own downstream `Preprocessor` order (normalize, log1p, gene subset, then bin; `scgpt/preprocess.py`), although pretraining binned over each cell's full gene set.
- **Whether scGPT runs.** Yes. It is preregistered, costs about 5 GPU-minutes, and referees will ask for it. It stays descriptive (D3).

### D7. Gene-mapping rule

**Decision.** Accept E5 as it stands.

**Alias and previous-symbol routes.** Non-exact routes cover 3.3% of mapped genes. The 40 distinct TFs mapped by them were reviewed by name: for example SFPI1→Spi1, CCRN4L→Noct, MLL3→Kmt2c, TCEB1/2/3→Eloc/Elob/Eloa, PHF17→JADE1, ZNF238→ZBTB18, PTRF→CAVIN1, DDX58→Rigi, Btbd11→Abtb3. All are standard renamings. An exact-only rule would drop them, including 21 mDC TFs, and shrink truth coverage for every method. A wrong alias could only hurt scPRINT (wrong gene embedding), so the remaining risk runs against our claim, not for it. The `route` column is published.

**mESC 15% loss.** Accepted. The lost genes are non-coding or predicted genes absent from the vocabulary, and all methods are scored on the same G∩. For continuity, report CLR-NS on the full gene set next to its G∩ value. Report dropped truth edges at scoring.

**Non-integer RSEM counts (hESC, arm B).** Used as is, without rounding. Rounding would change the zero pattern for values < 0.5. scPRINT takes `log2(1 + x)` of floats, and scGPT is rank-based.

### D8. Leakage classification

**Decision.** Accept "out, moderate confidence" for all scPRINT v1 and scGPT pairs. Leakage would favour the foundation models, so a missed "in" makes NI conservative and cannot inflate our claim. It matters only for interpreting a foundation-model win, which the preregistration already flags.

K562 x scPRINT-2 is classified **undeterminable** (a preregistered category), not "in", because the scBaseCount snapshot scPRINT-2 used is unknown. It is flagged in the optional, descriptive scPRINT-2 arm.

**Optional, PI's call.** Email the scGPT authors for the census snapshot. This does not block scoring. If an answer arrives before scoring, it is logged.

---

## Entries (suggested wording, to be dated on the day they are committed)

### E1. Checkpoints (hashes and revisions)

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

### E2. scGPT checkpoint source pin

Source: the "whole-human" link in the README of https://github.com/bowang-lab/scGPT, i.e. Google Drive folder `1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y`. The README lists no other source for this model.

Drive file IDs, listed from the folder on 2026-11-01 (no download):

- `best_model.pt` `14AebJfGOUF047Eg40hk57HCtrb0fyDTm`
- `vocab.json` `1H3E_MJ-Dl36AQV6jLbna2EdvgPaqvqcC`
- `args.json` `1hh2zGKyWAx3DyovD30GStZ3QlzmSqdk1`

Drive keeps no version history that we can pin, so the SHA-256 values in E1 are the pin.

### E3. Package versions and the scPRINT extraction route (proposed, PI to confirm)

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

### E4. scGPT binning seed (proposed, PI to confirm)

scGPT's `_digitize` breaks bin-edge ties at random (`side='both'`). Without a seed, two identical runs differed by 0.5% (relative L1). With `np.random.seed(seed)` immediately before the binning of each run, they are bitwise identical.

Entry: "The numpy global seed is set to **20261021** before binning in every scGPT run. A seed-1 rerun of hESC_1000 is reported as a truth-free stability check. Binning is over the case genes, as in scGPT's own downstream `Preprocessor` (normalize, log1p, gene subset, then bin), whereas pretraining binned over each cell's full gene set."

### E5. Gene-mapping rule, all 14 cases (proposed, PI to confirm: accepted as it stands)

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

### E6. Input-scale arms (proposed, PI to confirm: arm A primary)

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
- **No re-quantification** from SRA is done before scoring, unless the PI orders it. If ordered, its pipeline is committed before scoring, and the primary stays arm A.
- The organism-head accuracy under A and B is reported as a truth-free check that the input is in distribution.

### E7. Sink-head reporting (proposed, PI to confirm: R0' adopted)

The step-1 smoke test found heads that put 90-100% of a gene query's attention on the 9 special tokens. R0 remains the mean over all heads and layers, as preregistered.

Entry:

- For every case, record the cell-averaged special-token mass and the relative row entropy per layer x head next to the readout.
- **R0'** (descriptive) = the mean over heads with special-token mass < 0.5. If fewer than 4 of the 32 heads pass, use the 4 heads with the lowest mass.
- The threshold and the minimum are fixed here, from the smoke-test distribution, and use no truth.
- Each case's head set is computed in step 3 from expression only and committed before any truth file is opened.
- R0' enters the robust-NI maximum (E11).
- scGPT has no R0'; its `<cls>` mass is < 0.02 in the last layer.

### E8. Data notes (no effect on the design)

- All 14 case matrices equal the corresponding block of the BEELINE ExpressionData (all cells, case genes).
- mHSC: 67 (E), 45 (GM) and 17 (L) BEELINE cell labels do not carry the expression of the GEO column with the same name. They were re-matched by content to a different GEO column (identical zero pattern and a constant ratio, 100% agreement). Counts are therefore matched for every cell of every dataset.
- Leakage audit (metadata only): see E9.

### E9. Leakage audit outcome, done before scoring (proposed, PI to confirm)

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

### E10. Cells, batches and control seeds (new; proposed, PI to confirm)

**Cells and batches.**

- All cells of every BEELINE dataset, as preregistered: 4,794 per gene set.
- Batch size 8. All case genes are tokens, so there is no padding.
- Truth-free split-half stability per case: cells split once with seed 20261021; relL1 and Spearman between the half-means are reported in the step-3 QC.

**Controls.**

- Random weights: seeds 1-3, two versions (ESM2 gene embeddings kept; gene embeddings random).
- Expression permuted per gene across cells: seeds 1001-1003.
- K562: 5 subsamples of 1,024 non-targeting control cells (seeds 1-5), as preregistered.
- Fine-tuning arm (descriptive): package default losses, 10 epochs, seed 20261021, no truths.

### E11. Robust-NI set (new; proposed, PI to confirm)

Robust NI is called only if L > -δ also holds when scPRINT takes, per dataset, the best of {R0, R0', R1, R2} x {S = Aᵀ, A, max(A, Aᵀ)} x {arm A, arm B where available}. That maximum is chosen on the truths, so it is optimistic for scPRINT. This extends the preregistered maximum ({R0, R1, R2} x 3 orientations), which can only make robust NI harder to reach. It is the pre-specified fair-to-scFM reading. It is used only to stress our claim, never to claim a scPRINT win.

### E12. Descriptive arms and per-contrast gene sets (new; proposed, PI to confirm)

**Gene sets.** Each contrast is scored on its own G∩:

- CLR-NS, CLR or |r| vs scPRINT: genes covered by scPRINT.
- vs scGPT: genes covered by scGPT.
- scPRINT vs scGPT: the three-way intersection.

**Descriptive arms.**

- **scGPT** (hESC and hHep, plus K562): per-dataset d_i, plus the preregistered K562 bootstrap. No p-value and no NI claim.
- **Arm B:** as in E6.
- **CLR-NS on the full gene set** is reported next to its G∩ value.

The primary test, the margin, the interval, the decision table and the key secondaries are unchanged.

### E13. Step-3 freeze (to be committed after the step-3 QC and before any truth is opened)

- Per-case R0' head sets.
- The R1 kwargs actually used.
- The SHA-256 of every network file to be scored.
- The step-3 QC summary (`docs/scfm_step3.md`): NaN and row-sum checks, determinism, shim parity, split-half stability, permutation and random-weight distances, sink-mass tables, the A-vs-B distance and organism-head accuracy.

---

## Change-log entries to be committed (in this order)

1. **Entry "2026-11-xx: step-2 decisions"**, committed after PI confirmation and before any production extraction. It contains E1 to E12 in the confirmed wording.
   - E1 checkpoints, both revisions.
   - E2 scGPT source pin.
   - E3 route, versions, parity gate, R1 settings.
   - E4 binning seed.
   - E5 gene-mapping rule.
   - E6 input arms.
   - E7 sink heads and R0'.
   - E8 data notes.
   - E9 leakage.
   - E10 cells and seeds.
   - E11 robust-NI set.
   - E12 descriptive arms and G∩.

   Add one sentence: "no truth file opened, no foundation-model network computed for scoring".
2. **Entry "2026-11-yy: step-3 freeze" (E13)**, committed after production extraction and QC, and before the scoring code opens any truth file.
3. *(Only if the PI orders it)* **SRA re-quantification pipeline entry**, committed before scoring.

## Step 3: production extraction, exact order and cost

GPU jobs run through `jobs_pc/sched/`, chained with `# after:` so that only one GPU job runs at a time.

- Device 1 only (`CUDA_VISIBLE_DEVICES=1`); never device 0.
- Each job declares `# mem: 16` and `# cpus: 2`.
- Outputs go to `~/scfm_models/step3/`.
- No step opens a truth file.

The cost estimates below scale the measured rates with G²: scPRINT 0.05 s per cell at G = 1,161, and scGPT 0.15 s per cell at G = 1,221.

| Step | What | GPU time | Notes |
|---|---|---|---|
| 3.0 | PI confirms D1-D8. Commit change-log entry 1 to `docs/preregistration_scfm.md` and push. | 0 | The freeze point for all choices above. |
| 3.1 | Build the production scripts from the smoke scripts. | minutes | Scripts: `analysis/scfm_extract_scprint.py` (2.2.0 route, arms A and B, stores the per-layer x head cell-mean attention, special-token mass, entropy, organism-head predictions and split-half means), `analysis/scfm_extract_scgpt.py`, `analysis/scfm_controls.py`. Unit check: hESC_1000 at 64 cells reproduces the step-1 arrays exactly. About half a day of engineering. |
| 3.2 | Truth-free gates. Stop and report on any failure. | about 15 min | (a) Shim vs 2.2.0 on mHSC-E_500, all cells, ≤ 1e-6. (b) Determinism: hESC_500 run twice, bitwise identical. (c) Load test of `large-v1` and `small-v1` under 2.2.0. (d) `populate_my_ontology()` in a separate lamin instance, then GNInfer on hESC_500 returns a G x G matrix. |
| 3.3 | **Primary:** scPRINT medium-v1.5, arm A, 14 cases, all cells (9,588 cell-runs). | about 5 min plus 4 min of loads | Per-head tensors: 1.7 GB. |
| 3.4 | scPRINT, arm B, 8 cases. | about 3.5 min | 0.8 GB. |
| 3.5 | R1 (GNInfer) on arm A (14 cases) and arm B (8 cases). | ≤ 1 h (estimate; GNInfer untimed) | Kwargs logged. |
| 3.6 | scPRINT controls on arm A, 14 cases, R0 only. | about 45 min | Random weights 3 seeds x 2 versions (about 30 min); permutation 3 seeds (about 15 min). ESM2 cosine on CPU (seconds). Hybrid CLR is applied at scoring. |
| 3.7 | scGPT on 4 human cases, arm A: attention, token-embedding cosine and hidden-state cosine in one pass. Plus the seed-1 rerun of hESC_1000. | about 7 min | |
| 3.8 | K562: 5 x 1,024 control cells, scPRINT and scGPT. | about 20 min | Reuse the gene list of `results/baselines_k562ess/expression.tsv` if its export script implements the preregistered rule (check the script, not the truth file); otherwise rebuild by the rule at the scoring stage. Map its genes by E5. |
| 3.9 | Secondary checkpoints `large-v1` and `small-v1`, arm A, 14 cases, mean only. | about 30-60 min | |
| 3.10 | Fine-tuning arm (descriptive): scPRINT on 7 datasets and scGPT on 2, 10 epochs, then R0. | about 2-4 h | About 1 day of engineering. Its settings are fixed in E10. Because it is descriptive, it may finish after primary scoring. |
| 3.11 | Step-3 QC document `docs/scfm_step3.md`, then commit change-log entry 2 (E13), then push. | 0 | Only after this does step 4 (scoring) open the truths. |

**Total.**

- GPU: about 2.5-3 h for all zero-shot work (3.2-3.9) plus 2-4 h for fine-tuning, so 5-7 GPU-hours. That is well inside the preregistered estimate of under 2 GPU-days.
- Disk: about 4-5 GB new (`~/scfm_models` is at 4.1 of 20 GB).
- Host: ≤ 16 GB RAM per job.
- Wall time: 1-2 days, dominated by engineering (3.1, the GNInfer setup in 3.5, and 3.10), not compute.
- The optional scPRINT-2 and Geneformer arms are not scheduled.
