# Freeze record for the tcaCLR confirmatory plan: label-free build stage (17 Nov 2026)

Plan: `docs/plans/tcaclr_confirm_plan.md` (draft 3). **Nothing here was scored.** No AUPR, enrichment or per-complex or per-pair statistic of CLR, aCLR, tcaCLR, |r|, GENIE3, PIDC or GRNBoost2 against any truth exists. Operator matrices were computed only for the two label-free gates of section 6 and were not saved. Structure counts of truths (complexes, pairs, n_eff, overlaps) are the only truth-derived numbers.

Code commit used for everything below: `1e9178c` (branch accuracy-review). Small files are in `docs/tcaclr_confirm_selections/`; raw data and matrices are under `data/` and `results/` (gitignored).

## 0. What is filled and what is OPEN

| Item | State |
|---|---|
| Arm B (Kemmeren plus CYC2008): matrix hash, truth files and hashes, counts, gates | **FILLED** |
| Arm A data download, SHA-256 of the raw file | **FILLED** (no publisher checksum exists) |
| Arm A pseudobulk (HepG2): built, hashed, exported, gated | **FILLED but PROVISIONAL: it rests on forced deviations D1 to D5 (section 2) that need the PI's decision before this is a freeze.** If any of D1 to D5 is changed, every Arm A hash, count and gate below must be redone. |
| Arm A truth files, counts, n_eff, N and M, power gate | FILLED on the provisional matrix (the truth depends only on the gene universe, 9,622 genes, which D1 to D5 do not change; the perturbed-gene counts do depend on the profile set) |
| Operator gates P1 and reachability, both arms | FILLED (both pass) |
| Seeds and replicate counts | FILLED; the per-replicate seed derivation for SE_exp is OPEN (proposal in section 6) |
| Scheduling estimates | FILLED from the real dimensions (section 8); two plan statements are contradicted (section 8) |
| CYC2008 licence | OPEN (plan open item 2) |
| Planted-module simulation and H2's lower bin edge | OPEN (separate job, not part of this stage) |
| Method equivalence gates (plan gate 4), score gates (gate 5), GENIE3, GRNBoost2 and PIDC runs | OPEN (after freezing, by design) |
| Statement of the final freeze (PI sign-off, commit of this file) | OPEN |

## 1. Arm A data (GEO GSE264667, HepG2)

| Item | Value |
|---|---|
| File | `data/hepg2_perturbseq/GSE264667_hepg2_raw_singlecell_01.h5ad`, 5,614,460,941 bytes (equals the size on the GEO supplementary listing and the HTTP Content-Length) |
| Source | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/GSE264667_hepg2_raw_singlecell_01.h5ad`; downloaded in eight byte-range parts and concatenated |
| SHA-256 (computed here) | `e1ad7c3c5a201c861a207a858aa7e59f5e6ac1955674c415f7de0d1dadadb52e` (the GEO page and FTP directory publish no checksum; `md5sum.txt` does not exist) |
| Content | 145,473 cells by 9,624 genes, raw integer counts, 56 gem groups; 4,976 non-targeting cells (130 labels, 3 gem groups with fewer than 30 of them); 2,549 targeted labels (median 43 cells, minimum 1); one duplicated gene symbol (HSPA14) |
| Pseudobulk or energy-test p-values | **none exist for HepG2** (the K562 and RPE1 recipe starts from the authors' normalized-bulk h5ad with `obs/energy_test_p_value`; GEO has only the single-cell h5ad). The Jurkat fallback is not needed: the file is usable. |

## 2. Arm A pseudobulk: forced deviations, decided on label-free criteria, PROPOSED (need the PI's decision)

The K562 recipe (`analysis/perturbseq_bulk.R`, lines 1 to 71: significant-phenotype filter p < 1e-3, seeded half split, missing-value and constant-gene filters) is run **byte-identically**: `analysis/tcaclr_confirm/hepg2_expression_export.R` is those 71 lines unchanged (checked by `diff`, empty) followed by the duplicate-symbol rule and the export; only the command-line `--h5` and `--out` differ. The full diff is `docs/tcaclr_confirm_selections/armA_recipe_diff.txt` (the original's truth and scoring lines are not run; an export block is appended). What cannot be identical is its input file, which does not exist for HepG2. It is rebuilt from the single-cell file by `analysis/tcaclr_confirm/build_pseudobulk_from_sc.R`:

| # | Deviation | Label-free criterion and check |
|---|---|---|
| D1 | The normalized-bulk object is rebuilt: counts / UMI_count, z-score of every gene against the non-targeting cells of the same gem group, profile = mean z over the cells of a label | Defined by reverse-engineering the shipped K562 file from the K562 single-cell file. On K562 the rebuild matches the shipped profiles with per-profile correlation median 0.99999 (minimum 0.9988, 1st percentile 0.99984), regression slope median 1.0004, all 2,273 shared labels (`pseudobulk_validation_k562_profile_agreement.csv`). Gem-group-pooled controls instead of per-group controls gave correlation 0.988 (minimum 0.875); dividing by the file's `core_scale_factor` gave 0.94, so neither is the shipped definition. |
| D2 | The energy-test p-value is recomputed (HepG2 has none): energy distance in the first 50 principal components of the z-normalised cells (PCA on 25,000 random cells) against 1,000 random non-targeting cells, 10,000 permutations, p = (1 + exceedances) / (10,001) | Settings fixed before the first K562 run, not tuned afterwards. On K562 the call p < 1e-3 agrees with the shipped call in 97.3% of the 2,176 targeted labels (1,397 both, 17 only shipped, 42 only rebuilt; 1,439 rebuilt against 1,414 shipped); the selected set K has Jaccard 0.959 with the shipped one; Spearman of log p 0.967. The rebuild is slightly more liberal. Agreement is 95.7 to 97.8% in every cell-count class. |
| D3 | Gem groups with fewer than 30 control cells use the pooled control mean and SD | K562 has none (minimum 116 controls); HepG2 has 3 of 56 groups (minimum 5). The threshold 30 was fixed after seeing only the distribution of control cells per gem group (minimum 5, median 91.5) and before any expression value was read or built; it was not tuned, and it is the only rule invented for HepG2. |
| D4 | A label carried by a single cell cannot be tested and gets p = 1 | One HepG2 label; never selected. No K562 value is affected (the K562 rerun gave identical counts). |
| D5 | The HepG2 file has no `core_scale_factor` or `core_adjusted_UMI_count`; size normalisation uses `UMI_count` | This is the definition that reproduces K562 (D1). |

The pseudobulk builder's own constants (seed 20261120, 25,000 PCA cells, 50 components, 1,000 control cells, 10,000 permutations, ddof 1) are part of the frozen code. Logs: `pseudobulk_build_K562.log`, `pseudobulk_build_HepG2.log`.

### Arm A matrix (gate 1 values)

| Item | Value |
|---|---|
| Rebuilt bulk file | `results/tcaclr_confirm/pseudobulk/HepG2_rebuilt_bulk.h5ad`, SHA-256 `da8db71373debea5e1e97e11b094374c83bc29b08b20d5f1bd08e0c10e356314` (2,679 labels, 130 non-targeting) |
| Significant targeted perturbations (p < 1e-3) | 1,479, split with `set.seed(20261015)` into network half A = **739** and half B = 740 |
| Genes in file / after missing (> 1%) and constant filters | 9,624 / 9,624 (nothing dropped) |
| Duplicated gene names | 2 rows, symbol **HSPA14**, both dropped by the rule of `baselines.R` (all rows of a duplicated name; as for RPE1) |
| **Gene universe, profiles** | **9,622 genes x 739 profiles** (K562 had 8,563 genes; RPE1 8,745 x 1,023) |
| Expression matrix | `results/tcaclr_confirm/armA/expression.tsv`, 133,673,436 bytes, SHA-256 `bffc0e930393ecad2a36bca0ce1ccf2270374196f3673de8dddba1af0fc01074` |
| Profile list, universe list | `armA_profiles.txt` (739), `armA_universe_genes.txt` (9,622) |

## 3. Arm A truth: Complex Portal 9606 on the real HepG2 universe (counts only)

Built by `analysis/tcaclr_confirm/build_human_truth.R`, an adapted copy of the net4 builder with no rule change (recursive flattening, paralog molecule sets expanded and their 8 in-universe pairs excluded, size 3 to 200, overlap merge, ribosome list at 30% or more r-protein genes by `^RP[LS][0-9]+[A-Z]?$`, subunit variants). The only changes are the identifier map (UniProt reviewed human, primary gene symbol) and the four ribosomal-subunit complexes, located by their Complex Portal ids (40S CPX-5223, 60S CPX-5183, 28S CPX-5225, 39S CPX-5226). The net4 rule on one accession mapping to several genes does not apply to human (one primary symbol per accession). Source file `data/complex_portal_human_raw/9606.tsv`, SHA-256 `f9e977e0cce4f518de066eb267c70cabc2f8c1e659cfc1c164fbee41d4380d53`.

| Variant | Complexes | Truth genes | Pairs in universe | Positives | n_eff | Largest share |
|---|---|---|---|---|---|---|
| **main** | **359** | **2,040** | 2,079,772 | **17,417** | **31.69** | 0.073 |
| main_xsub_excl | 359 | 2,040 | 2,076,576 | 17,417 | 31.69 | 0.073 |
| main_ribo_joined | 357 | 2,040 | 2,079,772 | 20,613 | 15.11 | 0.165 |
| noribo | 356 | 1,956 | 1,911,986 | 15,681 | 31.20 | 0.081 |
| small (3 to 6) | 237 | 807 | 325,218 | 1,496 | 177.2 | 0.010 |
| mid (7 to 19) | 101 | 819 | 334,967 | 4,681 | 59.6 | 0.036 |
| mid 7 to 10 | 67 | 419 | 87,568 | 1,648 | 56.7 | 0.027 |
| mid 11 to 19 | 34 | 427 | 90,950 | 3,082 | 29.1 | 0.055 |
| large (20 or more) | 21 | 553 | 152,624 | 11,592 | 14.4 | 0.110 |

- **Against the proxy of the plan** (K562 gene names: 329 complexes, n_eff 31): the real universe gives 359 complexes and n_eff 31.69. The universe is larger (9,622 against 8,561 genes) and the largest share stays 0.073.
- Ribosome list: CPX-5183 (60S, 49 genes) and CPX-5223 (40S, 34 genes) and one 3-gene complex (CPX-2476); together 10.0% of positives. The 40S and 60S are separate complexes, so their cross-subunit pairs are negatives in `main`, the opposite of the K562 CORUM truth.
- **Data-truth link (counted label-free).** N = 2,040 truth genes. M = **418** of them are themselves perturbed in the profiles of the network matrix (half A). Under the wider definitions: 815 perturbed in any of the 1,479 significant knockdowns, 1,101 perturbed anywhere in the screen (`armA_complex_portal_perturbed_in_panel.csv`, all variants).
- **Power gate (plan section 6, gate 2): n_eff 31.69 is at least 27, so the registered sensitivity of +0.03 stands.** Power at +0.03 with k = 0.0417: 0.98 gene-only, 0.86 with the net4 experiment ratio (rho 0.893), 0.51 at rho 1.78; at +0.05: 1.00, 1.00, 0.91 (`power_gate.csv`). The plan's sentence in section 2 ("power 0.80 gene-only and 0.56 with rho 0.89") is the n_eff 15.6 row, not n_eff 31; it should be corrected.
- Files: `armA_complex_portal_members.tsv`, `armA_complex_portal_excluded_pairs.tsv`, counts, sizes, source log. Format as the net4 files (members: variant, dataset, complex, orf; excluded pairs: scope, orf_a, orf_b, reason; the `orf` columns hold gene symbols in the spelling of the expression matrix). `armA.sha256` (members and excluded pairs, the format of the score stage's hash guard) and `armA_expression.sha256`.
- Cross-check: the scoring library (`tc_lib.R`) reads these files independently and reproduces every positive and pair count above.

## 4. Arm B: Kemmeren deleteome and CYC2008

| Item | Value |
|---|---|
| Matrix | `data/yeast/deleteome_all_1484.tsv`, 6,169 ORFs x 1,283 profiles, 152,096,830 bytes, SHA-256 `a2bcc6aa9f5d5e36e98716aed1919bbf26d7b669237ea12bcd7a65f2a6266e63` (equals `data/yeast/SHA256SUMS`); no missing values, no duplicated names, no constant rows |
| Truth | `data/cyc2008_raw/CYC2008_complex.tab` (Internet Archive capture of 16 Feb 2020, SHA-256 `cdd17db9c04f969d3e16259e573415d7f7433e6acb42ee8067864813b61ffbd6`): 408 complexes. Paralog sets and the paralog-pair rule from the Complex Portal yeast release exactly as net4 (57 sets); ribosome rule `^RP[LS][0-9]+[AB]?$` as net4; the four ribosomal subunits located by their CYC2008 names. Builder: `analysis/tcaclr_confirm/build_cyc2008_truth.R`. |

| Variant | Complexes | Truth genes | Positives | n_eff | Largest share |
|---|---|---|---|---|---|
| **main** | **213** | **1,324** | **10,308** | **10.02** | 0.253 |
| main_xsub_excl | 213 | 1,324 | 10,308 | 10.02 | 0.253 |
| main_ribo_joined | 211 | 1,324 | 15,731 | 3.34 | 0.514 |
| noribo | 211 | 1,196 | 6,238 | 23.29 | 0.152 |
| small (3 to 6) | 155 | 582 | 951 | 109.1 | 0.016 |
| mid (7 to 19) | 49 | 481 | 2,659 | 34.9 | 0.058 |
| mid 7 to 10 | 27 | 214 | 779 | 25.1 | 0.058 |
| mid 11 to 19 | 22 | 268 | 1,880 | 19.8 | 0.081 |
| large (20 or more) | 9 | 313 | 6,800 | 4.43 | 0.383 |

- **Reproduces the candidates memo exactly:** 213 complexes, 10,308 positives, n_eff 10.02, and against the consumed Complex Portal yeast truth on the same universe (260 complexes, 13,223 positives, n_eff 15.61): 9,229 shared positives, **89.5% of CYC2008 positives shared**, 1,079 novel positives with n_eff 7.29, 172 of 213 complexes with a consumed complex at Jaccard 0.5 or more (`armB_overlap_with_consumed_cp_yeast.csv`).
- **Deleted genes that are truth members.** All 1,283 mutant columns map to a deleted ORF (SGD standard name, then systematic name, then unique alias; `armB_panel_deleted_orfs.tsv`). 340 of the 1,324 CYC2008 truth genes are deleted in the panel (for the consumed Complex Portal truth, 401 of 1,600). The memo's earlier 321 and 376 used standard names only and left 136 columns unmapped.
- Power (n_eff 10.02): at +0.03, 0.62 (gene-only), 0.40 (rho 0.893), 0.20 (rho 1.78); at +0.05, 0.97, 0.81, 0.46. No gate and no p-value in this arm.
- Files: `armB_cyc2008_members.tsv`, `armB_cyc2008_excluded_pairs.tsv`, counts, sizes; `armB.sha256` and `armB_expression.sha256`.

## 5. Hash table (all small files are in `docs/tcaclr_confirm_selections/`)

| File | SHA-256 |
|---|---|
| `armA_complex_portal_members.tsv` | `4a88e9461b99533a3b54b6a540c0df0c1cbccdb016a1124e77b23b050bef1b0d` |
| `armA_complex_portal_excluded_pairs.tsv` | `a5f79fde7a7fa13ab24b835fd298b527919f1b1140190b0d76205353a594d8fe` |
| `armA` expression.tsv | `bffc0e930393ecad2a36bca0ce1ccf2270374196f3673de8dddba1af0fc01074` |
| `armB_cyc2008_members.tsv` | `70af5c2063a6e57f26cd5a6bada9167cfcb4953c21267b62c1053684d1649d1e` |
| `armB_cyc2008_excluded_pairs.tsv` | `111e599096aee35ef978923a7e40b8a6417c3172a5300024826088fefabf0f0c` |
| `armB` deleteome_all_1484.tsv | `a2bcc6aa9f5d5e36e98716aed1919bbf26d7b669237ea12bcd7a65f2a6266e63` |

## 6. Seeds, replicates and operator settings

| Item | Value |
|---|---|
| seed_gene (delete-half over truth genes; halves drawn for all truths in table order) | 20261117, 500 replicates |
| seed_exp (delete-half over perturbations or mutants, CLR, aCLR and tcaCLR refitted) | 20261118, 100 replicates. **OPEN:** the per-replicate derivation; proposal: replicate b uses seed 20261118 + b - 1 (as net4 used 20264000 + b) |
| seed_complex (complex-level delete-half, sensitivity only) | 20261119, 500 replicates |
| Recipe split (Arm A) | `set.seed(20261015)`, p_max 1e-3 |
| GENIE3, GRNBoost2 | `set.seed(20261016)` (GENIE3, `doRNG`), seed 20261016 (GRNBoost2, legacy early stopping) |
| Pseudobulk builder (Arm A) | seed 20261120 and the constants of section 2 |

Operator (package code, `R/clr_attention.R`, as in `analysis/explore_stack/stack_lib.R::sk_a3r`): `ClrAttention$estimate_mi(bins = "hg", spline_order = 3, transform = "none")`, `calibrate(method = "normal", combine = "stouffer")` (row z clipped at 0, Stouffer), `build_operator(alpha = 0.5, softmax_keff = 10, softmax_cap = 50)`, lazy walk P = 0.5 I + 0.5 A, t = 9. aCLR = symmetrised P^9; tcaCLR = `sk_a3r(P^9)`: log P^9 with unreachable entries NA, row mean and SD over reachable off-diagonal entries, second CLR clipped at 0, Stouffer. The implementation gate confirms the package operator equals `rz_operator(CLR)` exactly (max difference 0 on both arms).

## 7. Label-free operator gates (no truth; score matrices not saved)

| Gate | Arm A (9,622 x 739) | Arm B (6,169 x 1,283) | Threshold |
|---|---|---|---|
| P1: tcaCLR at t = 1 against Stouffer of within-candidate row-z of CLR, max abs difference | 2.7e-15 | 3.6e-15 | 1e-10 (pass) |
| Rows with 10 or fewer candidates | 0 | 0 | |
| Unreachable share of P^9 pairs | 0.82% | 0.29% | at most 5% (pass) |
| Median row SD of log P^9 over reachable entries | 2.96 nats (IQR 2.48 to 3.43, CV 0.20) | 2.54 nats (IQR 2.40 to 2.71, CV 0.09) | at least 1 nat (pass) |
| Median row skewness, log P^9 against linear P^9 | -0.42 against 31.0 | -0.22 against 26.6 | descriptive (P3) |

Both arms are evaluable. Files: `armA_label_free_gates.csv`, `armB_label_free_gates.csv` with logs. Measured cost of the fit (MI, CLR, operator, P^9, both checks, 8 threads): Arm A 282 s, 7.8 GB; Arm B 114 s, 3.4 GB.

## 8. Scheduling estimates from the real dimensions (`schedule_estimates.csv`)

Anchors (32 workers): RPE1 8,745 x 1,023: GENIE3 435.6 min, peak 34.6 GB; GRNBoost2 221.7 min, 4.3 GB; PIDC 1,160 s, 10.5 GB. Second anchor net4 5,950 x 536 at 16 threads (GENIE3 12,435 s, GRNBoost2 11,764 s, halved for 32 threads). Scaling laws as in `tcaclr_prereg_candidates.md`. Declared memory = 1.4 x the scaled peak; the scheduler retries an out-of-memory kill once.

| Method | Arm A (9,622 x 739) | Arm B (6,169 x 1,283) |
|---|---|---|
| CLR, aCLR, tcaCLR, \|r\| (measured, 8 threads) | 5 to 7 min, 8 GB (declare 11) | 2 to 3 min, 3.4 GB (declare 5) |
| PIDC (`tools/pidc_chunked.jl`, 32 threads) | 20 to 27 min, 12.7 GB (declare 18) | 9 to 12 min, 5.2 GB (declare 8) |
| GENIE3 (32 workers) | **5.1 to 5.8 h**, peak about 28 to 31 GB; declare **48 GB** (1.4 x the RPE1 measurement; the size scaling alone gives 39 to 43) | **5.0 to 5.6 h**, same memory rule |
| GRNBoost2 (32 workers) | **3.2 to 5.9 h** (central 4.4), 4.7 GB (declare 7) | **2.3 to 4.2 h** (central 3.1), 3.0 GB (declare 5) |

Per-arm figures differ from the plan's (Arm A: GENIE3 about 7 h) because HepG2 has 739 profiles against RPE1's 1,023.

**Two plan statements are contradicted by these numbers.** (i) Plan section 5 calls PIDC refits for SE_exp "seconds to minutes each". One PIDC fit on the half-profile matrix takes about 11 to 18 min (Arm A) and 5 to 8 min (Arm B), so 100 replicates are about **19 to 30 h (Arm A) and 8 to 13 h (Arm B)** of 32-thread wall time. The cheap block (CLR, aCLR, tcaCLR, |r|) refits take 2 to 3 min (Arm A) and 1 min (Arm B). (ii) The scheduler budget (50 threads, 90 GB) allows the cheap methods, PIDC and GRNBoost2 to overlap, but two 32-worker runs cannot.

## 9. New disk use (not entered in `docs/disk_budget.md`)

| Location | Size | Note |
|---|---|---|
| `data/hepg2_perturbseq/` | 5.6 GB | raw HepG2 h5ad plus its SHA-256 file; needed to reproduce the pseudobulk |
| `results/tcaclr_confirm/pseudobulk/` | 317 MB | HepG2 rebuilt bulk 181 MB; K562 rebuilt bulk 136 MB (validation only, deletable) |
| `results/tcaclr_confirm/armA/` | 128 MB | the expression matrix (134 MB TSV) and lists |
| `results/tcaclr_confirm/armB/` | 12 KB | |
| `docs/tcaclr_confirm_selections/` | about 2.5 MB | committed |
