# Learned-attention comparator, step 2: preparation (2026-11-01)

Scope: gene-ID mapping, raw-count inputs, leakage audit, scPRINT checkpoint/code check, draft change-log entries. **Nothing was scored, no network was computed for any new case, and no truth file (BEELINE networks, Replogle, DoRothEA, CollecTRI, STRING, ChIP) was opened.** Gene lists come from `results/beeline_baselines/<case>/expression.tsv` and the TF lists `data/beeline/{human,mouse}-tfs.csv` (gene-set definitions, not truths). The preregistration `docs/preregistration_scfm.md` was not edited; proposed entries are in `docs/plans/scfm_step2_prereg_draft.md`. Step 1 (smoke test): `docs/scfm_smoke_test.md`.

Disk: `~/scfm_models` 3.1 GB before, 4.1 GB after (budget 20 GB); `~/venvs/scfm_alt` 3 MB. New downloads (about 0.9 GB, all under `~/scfm_models/raw/`, one directory each): `ensembl110_human_gtf` (GTF resumed from a truncated partial file of the interrupted earlier attempt; gzip -t ok), `ensembl110_mouse_gtf`, `hgnc`, `mgi`, `geo_GSE75748`, `geo_GSE81252`, `geo_GSE48968`, `geo_GSE98664`, `geo_GSE98664_mergedGTF`, `geo_GSE81682`, `cxg_api`, `scbasecount_meta`, `sra_runinfo`, `scprint_wheels`, `scprint_medium_v15_rev_d661d779`. The yeast human files (HepG2/Jurkat) were not touched. Nothing was deleted.

## Results at a glance

| Item | Result |
|---|---|
| 1. Gene mapping | Done for 14 cases. scPRINT covers 85.5% to 99.6% of case genes; all 7 datasets can run. TFs covered: all except 2 in mDC (TSC22D4 absent from Ensembl 110, MLL2 ambiguous). The 60 "ambiguous" genes of step 1 (POU5F1, GATA4, SOX7, ...) now map to their primary-assembly ID. |
| 2. Raw counts | **Counts exist for 4 of 7 datasets** (hESC, mHSC-E/GM/L). **Not available** for hHep, mDC, mESC (the GEO deposits hold only log2 values or TPM). All 5 deposits aligned; every BEELINE cell matched. Transform check: BEELINE `y = log2(1 + s_c x)` per cell, with `x` counts or TPM; so "BEELINE-scale" is not TPM in general. |
| 3. Leakage | Out for all model x dataset pairs on public metadata, moderate confidence; limits stated. Exception: Replogle K562 is in scBaseCount 2026-01-12 (relevant only to the optional scPRINT-2). |
| 4. scPRINT code/checkpoint | **scprint 2.2.0 (2025-03-06) loads `medium-v1.5` strictly with the package's own loader and reproduces the step-1 shim attention exactly (max abs diff 0).** 2.1.1 fails; 2.2.0, 2.2.2, 2.3.0, 2.3.8 load. New finding: the `medium-v1.5.ckpt` file was replaced on HF on 2026-05-05; weights identical, hash differs. |
| 5. Prereg draft | `docs/plans/scfm_step2_prereg_draft.md` (E1 to E9). |

## 1. Gene-ID mapping

Script `analysis/scfm_step2_genemap.py`, outputs `results/scfm_step2/{coverage_by_case.csv,dropped_genes.csv,genemap_meta.json,gene_map/<case>.tsv}`. The rule is documented in the script docstring and in the draft entry E5 (primary chromosomes only; exact symbol; curated HGNC/MGI ID; unique previous symbol, then unique alias/synonym; vocabulary membership; tie-break for duplicated names: in vocabulary, then curated ID, then protein_coding, then longest span). Mouse matching is case-insensitive because the mDC and mHSC files carry upper-case mouse symbols (QK, PARK2, CLEC1B).

The scPRINT checkpoint has one joint vocabulary of 44,756 IDs (23,149 ENSG, 21,607 ENSMUSG), so NCBITaxon:10090 uses the same gene list and the same encoder; only the organism class token differs. 3,141 of the human and 51 of the mouse IDs are non-primary contigs; the rule never selects them.

| case | genes | scPRINT covered | frac | TFs | TFs covered | scGPT covered | intersection | TFs in intersection |
|---|---|---|---|---|---|---|---|---|
| hESC_500 | 873 | 842 | 0.965 | 410 | 410 | 841 | 841 | 410 |
| hESC_1000 | 1323 | 1288 | 0.974 | 410 | 410 | 1285 | 1285 | 410 |
| hHep_500 | 911 | 906 | 0.995 | 448 | 448 | 907 | 905 | 447 |
| hHep_1000 | 1367 | 1362 | 0.996 | 448 | 448 | 1362 | 1360 | 447 |
| mDC_500 | 807 | 798 | 0.989 | 321 | 319 | n/a | n/a | n/a |
| mDC_1000 | 1270 | 1247 | 0.982 | 321 | 319 | n/a | n/a | n/a |
| mESC_500 | 1071 | 941 | 0.879 | 620 | 620 | n/a | n/a | n/a |
| mESC_1000 | 1528 | 1306 | 0.855 | 620 | 620 | n/a | n/a | n/a |
| mHSC-E_500 | 666 | 660 | 0.991 | 204 | 204 | n/a | n/a | n/a |
| mHSC-E_1000 | 1135 | 1122 | 0.989 | 204 | 204 | n/a | n/a | n/a |
| mHSC-GM_500 | 593 | 584 | 0.985 | 132 | 132 | n/a | n/a | n/a |
| mHSC-GM_1000 | 1055 | 1030 | 0.976 | 132 | 132 | n/a | n/a | n/a |
| mHSC-L_500 | 514 | 489 | 0.951 | 60 | 60 | n/a | n/a | n/a |
| mHSC-L_1000 | 692 | 647 | 0.935 | 60 | 60 | n/a | n/a | n/a |

(scGPT is human-only; mouse intersections are not defined and there is no ortholog mapping, per the preregistration.) Across the 14 cases the genes mapped by route are: exact primary 12,770; previous symbol 164; alias/synonym 269; tie-break 19. Drop reasons per case are in `dropped_genes.csv`: spike-ins (hESC only, 28 and 32), no match (e.g. read-through names such as `APOC4.APOC2`, `Gm`/`AK` clones), ambiguous synonym (at most 2), and "not in the scPRINT vocabulary" (mESC 95 and 170 genes, mostly `Gm`, `Mir`, `Snord` and other non-coding genes; mHSC-L 23 and 42).

Caveats. (a) Alias/previous-symbol mapping is by unique entry in the 2026 HGNC/MGI tables; a unique but wrong alias would silently remap a gene (269 + 164 genes use these routes; the list is in the per-case tables, column `route`, for inspection). (b) The mESC case loses 15% of its genes and mHSC-L 6% to 7%; the preregistration requires all methods to be scored on G, the covered intersection, so CLR, CLR-NS and |r| shrink to the same set. (c) **Dropped truth edges per case are not reported here** (they need the truth files); the per-gene tables let scoring code report them when the truths are opened. (d) The intersection with scGPT for the human cases is within 4 genes of the scPRINT set.

## 2. Raw counts and the transform check

Script `analysis/scfm_step2_counts.py`; outputs `results/scfm_step2/counts_alignment.{csv,json}`; model inputs `~/scfm_models/step2/inputs/<case>.npz` (47 MB; not in the repo). Each file holds `genes`, `cells`, `tpm_beeline` (`2^y - 1`, the preregistered input), `geo_linear` (the deposit's linear values where it has any, else NaN rows), `counts`, `counts_ok` (row has counts), `cells_ok`.

| Dataset | Accession (BioProject) | Deposit used | Counts? | BEELINE cells matched (verified) | Genes in BEELINE / matched / row-verified |
|---|---|---|---|---|---|
| hESC | GSE75748 (PRJNA305280) | `sc_time_course_ec`: RSEM expected counts, non-integer | **yes** (expected counts, not UMI) | 758 / 758 | 17,735 / 17,735 / 17,735 |
| hHep | GSE81252 (PRJNA322308) | `data.cast.log2.lineage`: log2 values | **no** | 425 / 425 | 11,515 / 11,432 / 11,432 |
| mDC | GSE48968 (PRJNA212439) | `allgenesTPM_GSM1189042_GSM1190902`: RSEM TPM | **no** | 383 / 383 | 7,371 / 7,371 / 7,366 |
| mESC | GSE98664 (PRJNA385832) | `tpm_sailfish_mergedGTF_RamDA_mESC_differentiation_time_course` (transcript TPM, summed to gene via the deposit's mergedGTF) | **no** (the count files cover only 10 pg and C1 cells) | 421 / 421 | 18,385 / 17,648 / 17,326 |
| mHSC-E / GM / L | GSE81682 (PRJNA322317) | `HTSeq_counts` (integer) | **yes** | 1,071 / 889 / 847, all matched | 4,762 / 4,683, 4,676, 4,670 / 4,682, 4,675, 4,670 |

Datasets where counts are unavailable: **hHep, mDC, mESC** (a raw-read re-quantification from SRA FASTQ would be needed; not attempted). The mESC 14 GB `RAW.tar` lists bigWig tracks, not counts; not downloaded.

Transform check (BEELINE `y` vs the deposit's linear values `x`; "const-ratio" = fraction of nonzero entries with `(2^y - 1)/x` within 1e-3 of that cell's median ratio `s_c`):

| Dataset | BEELINE linear scale | s_c median (min, max) | const-ratio | zero-pattern agreement |
|---|---|---|---|---|
| hESC | counts x s_c | 0.857 (0.344, 13.5) | 1.0000 | 1.0000 |
| hHep | identity (BEELINE = GEO log2 values, max abs difference 0) | 1 | 1.0000 | 1.0000 |
| mDC | TPM x s_c | 0.0200 (0.0006, 0.035) | 0.99998 | 0.99999 |
| mESC | summed TPM x s_c | 0.097 (0.065, 0.18) | 0.988 | 0.9991 |
| mHSC-E / GM / L | counts x s_c | 1.02 / 1.08 / 1.25 (0.07, 7.2) | 0.9998 / 0.9998 / 1.0000 | 0.9999 |

So BEELINE is `log2(1 + s_c x)` with a cell-specific `s_c`, not log2(TPM + 1) as step 1 assumed from hESC (there the column sums of `2^y - 1` were 0.83e6 to 5.7e6 because `s_c` is not 1e6/total). In mDC and mESC the BEELINE linear values are TPM divided by about 50 and 10. In hESC and mHSC they are within a factor of about 1 (median) of the counts, with a cell-wise spread of 0.3 to 14 (hESC). 322 mESC genes do not verify (their summed transcript TPM does not reproduce the BEELINE row; cause not investigated, plausibly a different transcript-to-gene assignment); they are excluded from `geo_linear`.

Cell labels: in mHSC, 67 (E), 45 (GM) and 17 (L) BEELINE cell names do not carry the expression of the GEO column with the same name (the labels follow a shifted chain, e.g. BEELINE `LT_HSC_012` has the expression of GEO `Prog_012`). They were re-matched by content (identical zero pattern, constant ratio), so counts are available for all cells. This is a BEELINE labelling issue, not a problem for the expression-only models, but it matters if cell labels are ever used.

Arms prepared (no choice made): **A** BEELINE-scale (`tpm_beeline`, all 14 cases); **B** GEO counts (4 of 7 datasets: hESC and mHSC x3); **C** GEO linear scale (`geo_linear`: expected counts, log2-derived FPKM-like values, TPM). A and B coincide up to the per-cell factor `s_c` for hESC and mHSC.

## 3. Leakage audit

Script `analysis/scfm_step2_leakage.py`; outputs `results/scfm_step2/leakage_*.csv|json`. Metadata only.

Method. scPRINT v1: the paper (Kalfon et al. 2025, *Nat Commun*, Data availability) states the pretraining set is CELLxGENE census **LTS 2023-12-15**; the README adds that some datasets were dropped for quality and some cells held out. scGPT whole-human: "33 million normal human cells" from CELLxGENE (README); **the census snapshot is not stated in any source reachable here**, so the three LTS releases (2023-05-15, 2023-07-25, 2023-12-15) were all checked. For each release the S3 listing of `h5ads/<dataset_id>.h5ad` (804, 850, 1,113 datasets) was joined to the CELLxGENE Discover API (404 public collections, 2,269 datasets) and every collection was searched for each GEO accession, BioProject, DOI and first-author-and-year. Positive control of the membership logic: 88.9%, 91.1% and 94.5% of datasets from collections published before 2023 are in the 2023-05-15, 2023-07-25 and 2023-12-15 listings. scPRINT-2: scBaseCount (SRX-keyed; 10x data only) sample metadata, releases 2025-02-25 and 2026-01-12, joined to the SRA run tables of each BioProject.

| Study (accession) | scGPT whole-human | scPRINT v1 | scPRINT-2 (optional) |
|---|---|---|---|
| hESC (GSE75748) | no, moderate confidence | no, moderate | no (0 of 1,810 SRX in either release) |
| hHep (GSE81252) | no, moderate | no, moderate | no (0 of 770) |
| mDC (GSE48968) | no (human-only corpus) | no, moderate | no (0 of 2,425) |
| mESC (GSE98664) | no (human-only corpus) | no, moderate | no (0 of 1,507) |
| mHSC (GSE81682) | no (human-only corpus) | no, moderate | no (0 of 3,840) |
| Replogle K562 (PRJNA831566) | no (normal cells only; not in any release) | no, moderate | no for 2025-02-25; **yes in 2026-01-12** (6,237 of 7,883 SRX); snapshot used by scPRINT-2 unknown |

What this does and does not show. No CELLxGENE collection cites any of the six studies; the five BEELINE studies are plate-based (SMART-seq, C1, RamDA), which scBaseCount does not contain. Evidence is "absent from the public metadata", which cannot exclude re-analysis of the same cells inside a large atlas that does not cite the study. 12 census datasets (19 dataset-release entries) have no current public collection (withdrawn); they were identified by title from the h5ad `uns` (CV19 T/B/M cells, hippocampus dissection, mouse motor cortex ATAC, pancreatic alpha and beta cells, cortex ATAC, great apes) and none is related (`leakage_unresolved_census_datasets.json`). Cell-line level overlap (hESC lines H1/H9 or K562 in other deposits) is not leakage of these cells and was not assessed. Under this audit the preregistered leaked-dataset sensitivity analysis is empty for scPRINT v1 and scGPT.

## 4. scPRINT checkpoint/code mismatch

Script `analysis/scfm_step2_scprint_alt.py`, output `results/scfm_step2/scprint_alt_check.json`, environment `~/venvs/scfm_alt` (`docs/scfm_env_freeze_alt.txt`). The checkpoint hyperparameters say `depth_atinput = False`, the state dict has no `depth_encoder` and `expr_decoder.fc.0.weight` is 256 x 257. Releases 1.0.13 and 1.6.1 have no `depth_atinput` option (always a depth token at the input, 256 x 256 decoder); from 2.0.0 the option exists.

- scprint **2.2.0** (released 2025-03-06; the checkpoint was uploaded on HF on 2025-03-04): `scPrint.load_from_checkpoint(ckpt, precpt_gene_emb=None, transformer="normal", strict=True)` succeeds with no missing and no unexpected key; the decoder layer is 256 x 257; 9 prefix tokens (base cell token plus 8 class tokens) before the genes, no depth token. Versions tested: 2.1.1 fails (`label_counts.cell_type_ontology_term_id` rejected), 2.2.0, 2.2.2, 2.3.0 and 2.3.8 load strictly; 2.0.0 was not tested.
- Parity on the 64 hESC cells of the smoke test (same cells, genes and input): attention with 2.2.0 versus the step-1 shim: max abs difference 0.0, relative L1 0.0 in all 8 layers. The step-1 shim reproduces the released code path exactly for the encoder.
- The alt environment installs only scprint 2.2.0 (`--no-deps`), simpler-flash 1.0.7 (1.3.5 needs `torch.amp.custom_bwd`, torch >= 2.4) and `sparse`; torch 2.3.1 and the rest are borrowed read-only from `~/venvs/scfm` through a `.pth` file. `~/venvs/scfm` was not modified. `populate_my_ontology` was not run; R1 (`GNInfer`) is untested.
- Finding about the file: HF `medium-v1.5.ckpt` changed on 2026-05-05 (commit `8e1ad1a1`, "update medium"); the earlier content (uploaded as `v2-medium.ckpt` 2025-03-04) has SHA-256 `a4cf0753...36c5`. All 147 state-dict tensors are identical; only hyperparameter metadata differ (`classes`, `prenorm`). The pinned file (`96ff7c5d...4ff8`) is the newer one. This also means scprint 1.0.13's rejection of `label_counts` came from this metadata format.

## Reproduce (exact commands)

```
source ~/venvs/scfm/bin/activate; export OMP_NUM_THREADS=2
cd ~/scfm_models/lamin_instance; REPO=/home/boris/Projects/ClrAttention
R="systemd-run --user --scope -q -p MemoryMax=16G -p MemorySwapMax=0 -p CPUWeight=10 nice -n 19"
# downloads (each into its own directory under ~/scfm_models/raw; -C - only when resuming a partial file)
B=https://ftp.ncbi.nlm.nih.gov/geo/series
curl -fsSL -o ensembl110_human_gtf/Homo_sapiens.GRCh38.110.gtf.gz https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
curl -fsSL -o ensembl110_mouse_gtf/Mus_musculus.GRCm39.110.gtf.gz https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz
curl -fsSL -o hgnc/hgnc_complete_set.txt https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
curl -fsSL -o mgi/MRK_List2.rpt https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt   # and MRK_ENSEMBL.rpt
# GEO: GSE75748_sc_time_course_ec.csv.gz, GSE81252_data.cast.log2.lineage.csv.gz, GSE48968_allgenesTPM_GSM1189042_GSM1190902.txt.gz,
#      GSE98664_tpm_sailfish_mergedGTF_RamDA_mESC_differentiation_time_course.txt.gz and GSE98664_mergedGTF.gtf.gz, GSE81682_HTSeq_counts.txt.gz
#      from $B/GSE<nnn>nnn/GSE<id>/suppl/<file>
# census listing: https://cellxgene-census-public-us-west-2.s3.us-west-2.amazonaws.com/?list-type=2&prefix=cell-census/<release>/h5ads/
# CELLxGENE: https://api.cellxgene.cziscience.com/curation/v1/collections?visibility=PUBLIC
# scBaseCount: https://storage.googleapis.com/arc-institute-virtual-cell-atlas/scbasecount/<release>/metadata/GeneFull/<Species>/sample_metadata.parquet
# SRA run tables: esearch db=sra term=<PRJNA>[BioProject] with usehistory, then efetch rettype=runinfo
$R python -I $REPO/analysis/scfm_step2_genemap.py        # step 1; writes results/scfm_step2 and ~/scfm_models/step2/scprint_vocab_ids.txt
$R python -I $REPO/analysis/scfm_step2_counts.py         # step 2
python -I $REPO/analysis/scfm_step2_leakage.py           # step 3 (metadata only)
$R ~/venvs/scfm_alt/bin/python -I $REPO/analysis/scfm_step2_scprint_alt.py --device 1   # step 4 (GPU device 1 only)
```

The alt venv was built with `uv venv ~/venvs/scfm_alt --python ~/venvs/scfm/bin/python`, `uv pip install --no-deps scprint==2.2.0 simpler-flash==1.0.7 sparse`, and a one-line `.pth` file in its site-packages pointing to `~/venvs/scfm/lib/python3.10/site-packages`.

## Decisions needed from the PI / Opus (nothing below was decided here)

1. **Input-scale arm.** A (BEELINE-scale, n = 7, matches the preregistration text but is neither TPM nor counts in general: counts x s_c for hESC and mHSC, TPM/50 for mDC, TPM/10 for mESC), B (GEO counts, only 4 of 7 datasets: hESC, mHSC-E/GM/L), C (GEO linear). The primary contrast needs all 7 datasets, which only A (or C) provides; B can only be a sensitivity analysis on n = 4. If B is wanted for hHep, mDC, mESC the raw reads must be re-quantified from SRA. Neither arm was run.
2. **Cells per case.** The preregistration says all cells (758 hESC, 425 hHep, 383 mDC, 421 mESC, 1,071 / 889 / 847 mHSC; 4,794 cells x 2 gene sets). Cost at the measured 0.05 s per cell (scPRINT) and 0.15 s (scGPT) is minutes per case. Step 1 showed the mean attention at 64 cells has not converged (64 vs 256 cells differ by relative L1 0.20), so all cells are recommended; the alternative is a fixed seeded subsample, which would have to be fixed now.
3. **Power.** The preregistered primary is one NI interval on n = 7 (smallest exact one-sided p 1/128 = 0.0078, attainable). The reviewer's problem (smallest two-sided p 2/128 = 0.0156 above the first Holm threshold 0.05/4 = 0.0125) applies only to the superseded 4-test plan or if the secondary tests are ever made two-sided or enlarged; the key-secondary Holm over two tests (threshold 0.025) is fine. Facts from this step that bear on it: no dataset has to be dropped (all 7 run, coverage at least 85.5%), so n = 7 stands; the human-only scGPT arm has n = 2 datasets (hESC, hHep) plus K562, so it can only be descriptive (smallest one-sided p 1/4 at n = 2); a counts-arm sensitivity has n = 4 (smallest one-sided p 1/16). Whether to adopt a different primary statistic is a design decision left to Opus.
4. **scPRINT route.** Freeze scprint 2.2.0 via `load_from_checkpoint` (no shim; identical numbers on hESC) as the extraction, with the 1.0.13 shim as cross-check, or keep the shim. Also whether to record both checkpoint hashes (E1).
5. **Sink heads.** Report special-token mass per layer x head only, or add the descriptive R0' (heads with mass below 0.5, threshold fixed before scoring) (E7).
6. **scGPT binning seed** for production (E4) and whether scGPT runs at all given n = 2 (plus K562).
7. **Gene-mapping rule** (E5): accept as is, in particular the use of synonym/alias routes (433 genes) and the 15% loss in mESC; alternative is exact symbols only (then more TFs are lost).
8. **Leakage classification** (E9): accept "out, moderate confidence", or require stronger evidence for scGPT (the census snapshot is not stated anywhere reachable; asking the authors or the scGPT preprint's methods would settle it).
