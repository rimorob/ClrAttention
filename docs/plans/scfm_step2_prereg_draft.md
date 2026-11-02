# Draft change-log entries for `docs/preregistration_scfm.md` (step 2, 2026-11-01)

**Status: DRAFT, not part of the preregistration.** `docs/preregistration_scfm.md` has not been edited. The PI / Opus decide which entries to commit, in what words, and which of the options marked **[DECISION]** to take. Everything below was prepared without opening any truth file and without scoring any network; the evidence is in `docs/scfm_step2.md` and `results/scfm_step2/`.

Suggested wording of the entries, to be dated on the day they are committed:

## E1. Checkpoints (hashes and revisions)

- scPRINT `medium-v1.5.ckpt`, HF `jkobject/scPRINT`, commit `8e1ad1a1bc40c1fc2f4f1700fe27923eea9a2bed` (2026-05-05), 221,742,980 bytes, SHA-256 `96ff7c5d9e1ad5076e818bde3c38f498c95cd354cc42029d90691b5ba2964ff8` (unchanged from the 22 Oct entry; re-verified 31 Oct and 1 Nov).
- **New fact.** The same file name had different content before 2026-05-05. The file was uploaded on 2025-03-04 as `v2-medium.ckpt` (commit `f60c8715`), renamed `medium-v1.5.ckpt` on 2025-12-11 (commit `fed0f6bb`), and replaced in commit `8e1ad1a1` ("update medium"). The earlier content (revision `d661d779e3614c77305e2fbe0a19312eda2794b8`, 221,775,592 bytes) has SHA-256 `a4cf0753270d4ff451a5dbddadd4c59a53e10e4c3e96a8af0db407ad893c36c5`. Compared tensor by tensor (`torch.load`, 147 state-dict tensors): **identical (max abs difference 0)**; epoch 19, global step 400,000 in both. Only hyperparameter metadata differ (`classes`: list of names before, name-to-count dict after; `prenorm`: True before, absent after). The weights analysed are therefore the same under either hash. The earlier revision is kept in `~/scfm_models/raw/scprint_medium_v15_rev_d661d779/`.
- Secondary scPRINT files at the same HF commit (not downloaded): `large-v1.ckpt` SHA-256 `c0c53989dfc0121cb7502652f643b491e43538e5bd616b9334c504ac6be03ef9` (1,237,154,305 bytes); `small-v1.ckpt` `269bc14be52fd0264631facecef400cd95bb1e03f36c23ff400ba02d7c9cd090` (111,120,529); `small-v2.ckpt` (scPRINT-2) `2b586c144cf9a1f638b4b3e803554ebf6ce389f81c5f34179288a970addfd822` (889,706,878). (Hashes from the HF LFS metadata; verify after download.)
- scGPT whole-human: hashes as in the 22 Oct entry (`best_model.pt` 6cb5d451ab5c4b33eb673adbe4fddc61d2389df1b89b7651a9fe2e557572b922, `vocab.json` acca93d114ca62c3f0f50debbd23e8c87f0714f4737764454f6b2b13f2e8580f, `args.json` c18e075e018140cb8b2d9029387b9de26607a5ce6a8ccabd6ead70cd76b95d60), re-verified at every smoke run.

## E2. scGPT checkpoint source pin

Source: the "whole-human" link in the README of https://github.com/bowang-lab/scGPT, i.e. Google Drive folder `1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y`. Drive file IDs listed from the folder on 2026-11-01 (no download): `best_model.pt` `14AebJfGOUF047Eg40hk57HCtrb0fyDTm`, `vocab.json` `1H3E_MJ-Dl36AQV6jLbna2EdvgPaqvqcC`, `args.json` `1hh2zGKyWAx3DyovD30GStZ3QlzmSqdk1`. The README lists no other source for this model. Drive has no version history that we can pin, so the SHA-256 values in E1 are the pin.

## E3. Package versions and the scPRINT extraction route **[DECISION]**

- Main environment `~/venvs/scfm` (freeze: `docs/scfm_env_freeze.txt`): torch 2.3.1+cu121 (sm_60 binary running on the Pascal cards, verified in the smoke test), scprint 1.0.13, scgpt 0.2.4 (`--no-deps`), lamindb 0.69.9, scdataloader 1.0.1, httpx 0.27.2 and setuptools 80.10.2 (the two step-1 fixes).
- **Finding.** scprint 1.0.13 (and 1.6.1) do not match the checkpoint: they inject a depth token at the encoder input and have a 256x256 decoder layer, whereas the checkpoint has `depth_atinput = False`, no `depth_encoder` weights and a 256x257 decoder layer. The first PyPI release that loads the checkpoint is **scprint 2.2.0** (released 2025-03-06, two days after the checkpoint upload): `scPrint.load_from_checkpoint(ckpt, precpt_gene_emb=None, transformer="normal", strict=True)` loads with no missing or unexpected key (the package's own `on_load_checkpoint` renames `gene_encoder.embedding.weight` and drops the training-only discriminator keys). Checked: 2.1.1 fails (hyperparameter parsing); 2.2.0, 2.2.2, 2.3.0 and 2.3.8 load strictly; 2.0.0 not tested. Environment: `~/venvs/scfm_alt` (scprint 2.2.0 and simpler-flash 1.0.7 installed with `--no-deps`; torch and the rest borrowed read-only from `~/venvs/scfm` through a `.pth` file; `docs/scfm_env_freeze_alt.txt`).
- **Parity.** On the 64 hESC cells of the smoke test, the layer x head attention computed with scprint 2.2.0 equals the step-1 shim result exactly (max abs difference 0.0, relative L1 0.0, all 8 layers). The shim in `analysis/scfm_smoke_scprint.py` is therefore a faithful reproduction of the released code path.
- Proposed entry: "R0 is extracted with scprint 2.2.0 loaded through the package's own `load_from_checkpoint` (no shim), in `~/venvs/scfm_alt`; the 1.0.13 shim is retained as a cross-check only. The depth token is not fed (`depth_atinput = False`), as in the checkpoint's training." Alternative: keep the 1.0.13 shim as primary (identical numbers on the one case tested). R1 (`GNInfer`) is not run unless wanted; it needs `populate_my_ontology()`.

## E4. scGPT binning-seed fix

scGPT's `_digitize` breaks bin-edge ties at random (`side='both'`). Without a seed two identical runs differed by 0.5% (relative L1); with `np.random.seed(seed)` immediately before the binning of each run they are bitwise identical. Entry: "the numpy global seed is set to the run seed before binning; the production seed is [DECISION: 1, as in the smoke test, or a fixed seed drawn before scoring]; binning is over the case genes (the pretraining binned over each cell's full gene set)".

## E5. Gene-mapping rule (all 14 cases)

Fixed rule, `analysis/scfm_step2_genemap.py` (full text in its docstring). Sources, downloaded 2026-11-01: Ensembl release 110 GTFs (human GRCh38, SHA-256 `5d3c363835bda4cf8ba025f728b190f69d454e894514bb3c81123efdb15ef1b9`; mouse GRCm39, `fd89f56d03871fac5b05a699fa2ea0ce46245d35b64248cf74d9f70b59b61618`), HGNC complete set (`d93ce6b22d01dc1e1ded28857246951d605ead47ed17ef973ad1142b1d13bed0`), MGI `MRK_List2.rpt` (`00c4dd2abd9e2fe70b7810ee2914a6af65015dfd620fd493de2ce7f811cc915f`) and `MRK_ENSEMBL.rpt` (`1c8ee6e696bc30dbe3eea7d73e58a1fe9c6ec38674db8ddf2f762babcae096eb`).

1. Primary assembly = chromosomes 1-22/19, X, Y, MT. Non-primary (patch, haplotype, scaffold) IDs are never chosen, even though the scPRINT vocabulary contains 3,141 human and 51 mouse such IDs.
2. Exact symbol match among primary genes (human: case-sensitive; mouse: case-insensitive, because the mDC and mHSC files use upper-case mouse symbols). If several primary genes share the name: keep those in the scPRINT vocabulary, then the HGNC/MGI-curated Ensembl ID, then protein_coding, then the longest genomic span; if still not unique, drop (ambiguous).
3. Else the curated HGNC/MGI ID of an approved symbol.
4. Else a unique HGNC previous symbol, then a unique HGNC alias / MGI synonym, resolved to the current symbol and mapped by 2-3; ambiguous synonyms are dropped. ERCC spike-ins are dropped. Excel-date symbols are mapped (none occurred).
5. The ID must be in the scPRINT checkpoint vocabulary of that organism (23,149 human, 21,607 mouse IDs; the checkpoint is one joint vocabulary). Two symbols resolving to one ID: the exact-match symbol is kept, otherwise both are dropped.
6. scGPT (human cases): the case symbol if in the 60,697-symbol vocabulary, else the resolved current symbol. No ortholog mapping for mouse.

Result (`results/scfm_step2/coverage_by_case.csv`, per-gene tables in `results/scfm_step2/gene_map/`, dropped genes with reasons in `dropped_genes.csv`): scPRINT covers 85.5% to 99.6% of the case genes (lowest mESC_1000 1,306 of 1,528 and mHSC-L_1000 647 of 692; the loss is non-coding or predicted genes missing from the vocabulary). TFs covered: 2 of 321 lost (TSC22D4 absent from Ensembl 110, MLL2 ambiguous) in mDC, all others covered. scGPT intersection: hESC 841 / 1,285 genes and hHep 905 / 1,360 genes (500 / 1000 sets). Dropped truth edges per case are not reported: that needs the truth files, which stay closed until scoring (the per-gene mapping tables allow it to be computed then).

## E6. Input-scale arms **[DECISION: which is primary]**

The preregistered input ("BEELINE's documented transform, back-transformed") is `2^y - 1` of the BEELINE ExpressionData. The transform check (all 7 datasets, `results/scfm_step2/counts_alignment.json`) shows that this is not TPM in general: BEELINE `y = log2(1 + s_c * x)` with a per-cell scale `s_c` and `x` the linear value of the GEO deposit, to 1e-3 relative precision in 99.98-100% of entries (mESC 98.8%). Linear scale of the GEO source and typical `s_c` (median):

| Dataset | GEO linear scale x | s_c median (range) |
|---|---|---|
| hESC | RSEM expected counts (not UMI) | 0.86 (0.34 to 13.5) |
| hHep | log2 values only; BEELINE equals them exactly | 1 (identity) |
| mDC | RSEM TPM | 0.020 (0.0006 to 0.035) |
| mESC | sailfish transcript TPM summed to gene | 0.097 (0.065 to 0.18) |
| mHSC (E, GM, L) | HTSeq integer counts | 1.02 to 1.25 (0.07 to 7.2) |

Arms prepared (`~/scfm_models/step2/inputs/<case>.npz`, arrays `tpm_beeline`, `geo_linear`, `counts`, masks): **A, BEELINE-scale** (`2^y - 1`, all 14 cases); **B, GEO counts** (only where the deposit holds counts: hESC and the three mHSC datasets, 4 of 7; unavailable for hHep, mDC, mESC); optionally **C, GEO linear scale** (TPM or FPKM-like, where the deposit has it). Neither arm is chosen here. Facts for the choice: scPRINT and scGPT were pretrained on raw counts; arm B covers only 4 of 7 datasets, so a primary contrast on n = 7 needs arm A; arm B is natural as a sensitivity analysis on n = 4 (which cannot reach p < 0.1 two-sided: smallest exact one-sided p is 1/16).

## E7. Sink-head reporting

The step-1 smoke test found heads that put 90-100% of a gene query's attention on the 9 special tokens. R0 remains an average over all heads and layers as preregistered. Proposed addition (reporting only, no new test): for every case record the special-token mass per layer x head and the per-head row entropy next to the readout, and add as a descriptive sensitivity R0' the average over heads with special-token mass below 0.5 (the threshold is fixed here, before scoring, from the smoke-test distribution and uses no truth). **[DECISION: adopt R0' or report mass only.]**

## E8. Data notes (no effect on the design)

- All 14 case matrices equal the corresponding block of the BEELINE ExpressionData (all cells, case genes).
- mHSC: 67 (E), 45 (GM) and 17 (L) BEELINE cell labels do not carry the expression of the GEO column with the same name; they were re-matched by content (identical zero pattern and a constant ratio, 100% agreement) to a different GEO column. Counts are therefore matched for every cell of every dataset.
- Leakage audit (metadata only): see E9.

## E9. Leakage audit outcome (done before scoring)

Corpora: scPRINT v1 = CELLxGENE census LTS 2023-12-15 (Kalfon et al. 2025, Data availability); scGPT whole-human = "33 million normal human cells" from CELLxGENE, snapshot not stated in any reachable source (all three LTS releases 2023-05-15, 2023-07-25, 2023-12-15 checked); scPRINT-2 = scBaseCount. Result: none of the six accessions (GSE75748, GSE81252, GSE48968, GSE98664, GSE81682, Replogle PRJNA831566) is cited by any CELLxGENE public collection (accession, BioProject, DOI, first author and year); none has an SRX in scBaseCount 2025-02-25; Replogle K562 has 6,237 of 7,883 SRX in scBaseCount 2026-01-12. Classification: **out** for every model x dataset pair, with the stated limits; the K562 x scPRINT-2 pair is "in" if scPRINT-2 used the 2026-01-12 scBaseCount (snapshot unknown). Residual limits: undetected re-use of cells inside an atlas that does not cite the study; 12 census datasets without a current public collection were identified by title and none is related. Table: `results/scfm_step2/leakage_audit.csv`. The preregistered leaked-dataset sensitivity analysis is empty under this classification.
