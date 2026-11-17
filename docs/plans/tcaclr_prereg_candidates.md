# Candidates for a new preregistered complex test of tcaCLR (internal a3r): data, truth, power and cost

Fact-finding only. Retrieved 2026-11-16 (web checks, local counts). **Nothing was scored**: no AUPR, enrichment or per-complex statistic of any method was computed against any truth in this work. Counted structure only (complexes, sizes, genes, pairs, n_eff, overlap). Part C re-reads numbers that already exist in `results/complex_confirm/` and `docs/plans/complex_confirm_design.md`; it computes no new score.

Tags: **[verified]** = checked in this session against a primary source or a file on disk; **[not verified]** = taken from memory, a secondary source, or an inference; **[assumed]** = a stated modelling assumption.

Code written for this memo (originals untouched; no commit):
- `analysis/prereg_candidates/count_candidate_truths.R`: adapted copy of `analysis/complex_confirm/build_yeast_complex_truth.R` (same flatten, paralog-set, size 3-200, overlap-merge and `structure_of` code). It reproduces the frozen consumed truth exactly on the net4 universe (262 complexes, 1,602 genes, 13,438 positives, n_eff 14.48) **[verified]**, then applies the same code to the Kemmeren universe.
- `analysis/prereg_candidates/count_human_complex_portal.R`: same code for human Complex Portal.
- Outputs (counts only): `results/prereg_candidates/*.csv`.

Convention (as the design doc): complexes of 3-200 in-universe genes, overlap merge (drop duplicates and proper subsets, complete linkage on 1 - Jaccard cut at 0.5, union), paralog pairs of the Complex Portal paralog sets (57 pairs) removed from the pair universe, ribosome subunits kept as in the source. **n_eff** = inverse Herfindahl index of positive-pair shares across merged complexes (fractional assignment of shared pairs). Size classes are those of the design (3-6, 7-19, 20+ genes of the merged complex); n_eff per class is computed on the class's own positives.

## 0. Summary

| Question | Answer |
|---|---|
| Recommended primary pair | **Kemmeren deleteome (on disk, 6,169 ORFs x 1,283 mutant profiles) + CYC2008 complexes** (408 source complexes, 213 merged on the universe, 10,308 positive pairs, n_eff 10.0). Data are genuinely independent of net4 (different lab, platform and design). **The truth is not: 89.5% of its positive pairs are also positives of the consumed Complex Portal truth.** |
| Recommended backup | Human Complex Portal (9606) with HepG2 or Jurkat Perturb-seq (GSE264667), but it needs a 5.6 GB or 9.4 GB single-cell download and a new pseudobulk build, and the data come from the same lab and pipeline as consumed K562/RPE1. Not run here. |
| Power | Adequate for +0.05 pooled AUPR (0.97 gene-only, 0.81 with an experiment component like net4's) on CYC2008. **Not adequate for +0.03** (0.63 / 0.40) and **not adequate for the 20-plus stratum** (n_eff 4.4; about 16 needed for +0.03 on the gene-only SE). |
| Blocked / needs the PI | GENIE3 and GRNBoost2 on Kemmeren are each estimated above 2 hours (section D). CYC2008 had to be taken from an Internet Archive capture because the authors' server does not answer. |

## A. Yeast

### A1. Kemmeren et al. 2014 deleteome (Cell 157:740, PMID 24766815)

| Item | Finding | Tag |
|---|---|---|
| Accession | **GEO SuperSeries GSE42528** ("Compendium of deletion mutant gene expression profiles"), made of **GSE42527** (700 responsive mutants, sub-series label hs1991) and **GSE42526** (784 non-responsive mutants, hs1990). The recalled **GSE42215 is not the deleteome**: it is "yeast wt pool" (200 hybridisations). | verified (GEO text pages, FTP listing, local series-matrix headers) |
| Platform | GPL11232 (*S. cerevisiae*, Holstege lab two-channel array); probe IDs `SCAB######`; annotation file on disk `data/yeast/GPL11232_full.txt` | verified |
| Samples | SuperSeries lists 2,633 GSM (hybridisations): 1,387 in GSE42527 + 1,246 in GSE42526 | verified |
| Design | Two independent cultures per mutant, hybridised on two dye-swapped two-channel arrays against one common wild-type reference RNA; each gene spotted twice; processed value per mutant = limma average log2(mutant/reference) "M" and an adjusted p-value | verified (GEO overall design, protocols file) |
| Mutants | 1,484 deletion strains (700 + 784) per the series summary | verified |
| Processed matrices (GEO) | `GSE42527_final_GeneExpressionMatrix.txt.gz` 61 MiB (63,790,953 B; 1,401 columns = 700 mutants x (M, p) + ID) and `GSE42526_final_GeneExpressionMatrix.txt.gz` 67,507,834 B; probe-level, mapped to ORFs through GPL11232. GSE42527_RAW.tar 917 MB, GSE42528_RAW.tar 1.7 GB (not needed). | verified |
| Wild-type / control arrays | **Yes.** GSE45115: 376 wild-type cultures grown in parallel with the mutants (755 hybridisations); GSE42215 / GSE42217: wild-type pools (200 + 200); GSE42240 / GSE42241: MATa wild types. The mutant series themselves contain no wild-type arrays. | GSE45115 and GSE42215 verified today; 42217/42240/41 as in `docs/preregistration_reciprocity_yeast.md`, not re-checked |
| Already on disk | **Yes, no download needed.** `data/yeast/deleteome_all_1484.tsv` (152,096,830 B; 6,169 ORF rows x 1,283 mutant columns), `deleteome_responsive_700.tsv` (76,842,569 B; 6,169 x 652 columns), `wt_controls_GSE45115.tsv`, SHA-256s in `data/yeast/SHA256SUMS`. Built by `analysis/yeast_prep.R` (neighbour masking within 500 genes of the deleted gene, columns kept only if the deleted gene's own M is at most -1 and below the column's 1% quantile, genes with more than 5% missing dropped). The file name "1484" is historical: 1,283 columns survive (652 + 631; 48 responsive mutants dropped in `docs/reciprocity_yeast_selections/dropped_columns.csv`, the dropped non-responsive ones are not listed there). The 6,169 row names are all ORF-format. | verified |
| Gene universe | 6,169 ORFs (net4 has 5,667 non-decoy ORFs) | verified |

**Risks a reviewer will raise (all from this repository's own record, verified):**
1. **Gate F of the 1 Nov reciprocity preregistration failed on this matrix**: median chance-corrected reproducibility 0.122 for mutants against 0.296 for wild types (`docs/reciprocity_yeast_selections/step4_gateF_fail_report.md`). That gate belonged to another hypothesis, but it says the deleteome's mutant network is less reproducible than its wild-type network. The same report records a hybridisation-orientation problem in sub-series hs1991, which is the 700-responsive series: a fitted orientation rule reproduced the published M poorly there (held-out median r 0.42, against 0.93 in the other sub-series). Whether the published M (as on disk) or a rebuilt M is "right" for hs1991 is **not verified**. The 1,283-column matrix includes both sub-series; the 652-column matrix is hs1991 only.
2. **The Kemmeren expression has been opened before, label-free only** (gate F reproducibility, structure counts of its row names). No tcaCLR/a3r, aCLR or CLR score has ever been computed on it, and no complex truth was scored on it. This is a counting statement from the design documents, not a re-audit of every results directory beyond what was grepped (`docs/`, `analysis/`): **[not verified beyond the grep]**.
3. **Built-in truth-data link.** The panel deletes genes of "transcription, mRNA turnover, signalling, nucleus". Of the 1,147 mutants whose symbol maps to an ORF (via SGD standard name; the other 136 columns did not map), **376 delete a gene that is itself a member of a consumed-Complex-Portal complex** (23.5% of the 1,600 truth genes), 321 for CYC2008 (24.2% of 1,324), 43 and 28 in the 20-plus complexes (`results/prereg_candidates/candidate_truth_deleted_gene_membership.csv`). Deleting one subunit perturbs the others, so co-membership can be partly built into the design. It is not expression-derived truth, but it is a data-truth dependence the DREAM5 net4 compendium did not have to the same degree.
4. Replicate structure for experiment-level resampling: the unit is the mutant (independent culture pair); there is no "series" structure like net4's 59 series, so net4's ratio of experiment SE to gene SE (below) cannot be transferred and is **not estimated for Kemmeren**.

### A2. Candidate yeast complex truths (independent of Complex Portal?)

| Candidate | Obtainable? | One-line verdict | Tag |
|---|---|---|---|
| **CYC2008** (Pu et al. 2009, *Nucleic Acids Res* 37:825, PMC2647312): 408 manually curated heteromeric complexes from small-scale literature experiments | The authors' site `http://wodaklab.org/cyc2008/` (file `resources/CYC2008_complex.tab`) gave **no HTTP response within 25 s on 16 Nov** (also two timeouts on 13 Nov per the design doc). **Obtained from the Internet Archive**, capture 2020-02-16 15:57:07 of that exact URL: `https://web.archive.org/web/20200216155707/http://wodaklab.org/cyc2008/resources/CYC2008_complex.tab`. 108,324 B, SHA-256 `cdd17db9c04f969d3e16259e573415d7f7433e6acb42ee8067864813b61ffbd6`, 1,920 ORF-complex rows, **408 complexes, 1,627 ORFs** (paper: 408 complexes, 1,630 genes). Licence: **not stated** on the paper page (the article is CC BY-NC); the wayback copy is a third-party archive of the authors' file. | **Usable; the best available independent-source yeast truth, but largely redundant with Complex Portal in content (section A3).** Uses literature small-scale evidence (affinity capture-Western most common), no expression data. | availability and counts verified; licence not verified |
| YHTP2008 (second catalogue of the same paper: 400 complexes from high-throughput purifications) | The guessed file name `YHTP2008_complex.tab` is not in the archive (404). | Not obtained; high-throughput-derived (Gavin/Krogan data), not a clean truth. | not verified |
| Pu 2009 beyond the two catalogues | n/a | CYC2008 *is* the Pu et al. 2009 curated set; there is no further independent file. | verified (paper) |
| Gavin 2006 / Krogan 2006 consolidated sets | Not looked for as ORF lists. | **Not recommended:** raw TAP-MS complexes, noisy, overlapping, and not an independent truth by design. | not verified |
| Benschop et al. 2010 "consensus of core protein complexes" (Mol Cell) | Not obtained. | **Not recommended:** from the Holstege lab that produced the deleteome (shared-lab dependence; authorship [not verified]). | not verified |
| GO cellular-component complexes (SGD GAF) | Not obtained or counted. | Possible but derived from the same curation literature as Complex Portal and CYC2008; adds no independence and its evidence codes include electronic annotation. Not counted. | not verified |
| hu.MAP (human, for completeness) | n/a for yeast | see B | |

### A3. Structure counts on the Kemmeren universe (6,169 ORFs); counts only

Source counts before merge: Complex Portal 649 complexes in file, 378 with at least 3 universe genes; CYC2008 408 complexes, 235 with at least 3 universe genes.

| Truth (Kemmeren universe) | Class | Merged complexes | Genes | Positive pairs | n_eff | Largest share | Ribosome share |
|---|---|---|---|---|---|---|---|
| **Consumed reference: Complex Portal 559292** | all | 260 | 1,600 | 13,223 | 15.61 | 0.171 | 0.273 |
| | 3-6 | 196 | 712 | 1,226 | 141.6 | | |
| | 7-19 | 52 | 527 | 3,011 | 36.6 | | |
| | 20+ | 12 | 415 | 9,100 | 7.54 | 0.248 | 0.397 |
| **CYC2008** | all | **213** | 1,324 | **10,308** | **10.02** | 0.253 | 0.395 |
| | 3-6 | 155 | 582 | 951 | 109.1 | | |
| | 7-19 | 49 | 481 | 2,659 | 34.9 | | |
| | 20+ | 9 | 313 | 6,800 | 4.43 | 0.383 | 0.599 |
| CP + CYC2008 union (the 1 Nov fallback rule) | all | 278 | 1,687 | 14,372 | 15.25 | 0.181 | 0.283 |
| | 3-6 | 204 | 721 | 1,235 | 143.5 | | |
| | 7-19 | 61 | 586 | 3,171 | 42.1 | | |
| | 20+ | 13 | 455 | 10,114 | 7.66 | 0.258 | 0.402 |

[verified; counts only]. The consumed reference on Kemmeren (260 complexes, 13,223 positives, n_eff 15.6) matches the rev 1 design-doc count for this universe within rounding of the paralog rule (260, 13,225, 15.6). The same code on the net4 universe returns the frozen 262 / 13,438 / 14.48. Per class, positives are those of complexes in that class only (the design's secondary item 3).

**How independent is each truth of the consumed Complex Portal truth** (same universe, `results/prereg_candidates/candidate_truth_overlap_with_consumed.csv`):

| Measure | CYC2008 | CP + CYC2008 union |
|---|---|---|
| Gene Jaccard with consumed truth genes | 0.733 (1,237 of 1,324 genes also in it) | 0.948 |
| Positive pairs shared with consumed truth | 9,229 | 13,223 (all of the consumed ones, by construction) |
| **Share of this truth's positives that are consumed positives** | **0.895** | 0.920 |
| Share of consumed positives present in this truth | 0.698 | 1.000 |
| Pair-set Jaccard | 0.645 | 0.920 |
| Merged complexes with a consumed complex at Jaccard 0.5 or more / 0.8 or more | 172 of 213 / 142 | 255 of 278 / 244 |
| Merged complexes with no consumed match at Jaccard 0.5 | **41 (19%)** | 23 (8%) |
| Positives **not** in the consumed truth | **1,079 (10.5%)** | 1,149 (8.0%) |
| n_eff of those novel positives alone (complex-wise) | 7.29 | 7.77 |

**Plain statement.** CYC2008 is an independent *source* (different curation, small-scale literature evidence, built 2008-2009), but in *content* about nine tenths of its positive pairs and four fifths of its complexes are already in the consumed truth. A tcaCLR-versus-aCLR result on CYC2008 is therefore a test on new expression data with largely the same complexes, not on new complexes. The independence is on the data side, not the truth side. Using the union is not independent at all (92% of its positives are consumed). A positives-not-in-Complex-Portal truth (1,079 pairs, n_eff about 7.3) is the only truly new part and is too small to power a +0.03 or +0.05 test on its own (SE by the k rule about 0.0154 gene-only, MDE about 0.043; section C arithmetic).

## B. Human backup

### B1. Perturb-seq in a cell type other than K562/RPE1

| Dataset | Accession, size, format | Cells / perturbations | Pseudobulk provided? | Verdict | Tag |
|---|---|---|---|---|---|
| **Nadig et al. 2025, HepG2 and Jurkat essential-gene Perturb-seq** (Nat Genet, PMID 40259084; GEO series title "Transcriptome-wide characterization of genetic perturbations") | **GSE264667**. Supplementary files: `GSE264667_hepg2_raw_singlecell_01.h5ad` **5,614,460,941 B (5.6 GB)**, `GSE264667_jurkat_raw_singlecell_01.h5ad` **9,366,490,264 B (9.4 GB)**, `GSE264667_RAW.tar` 16,321,884,160 B (16.3 GB). Raw reads on SRA PRJNA1100571. URL: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/` (directory listing is forbidden; file URLs answer HEAD). **Not downloaded (over 1 GB).** | CRISPRi, DepMap common-essential library, 10x v3, day 7; HepG2 and Jurkat run in parallel. Counts of cells, perturbed genes and measured genes **not verified** (a search snippet gave unlabelled numbers that I did not rely on; design doc says about 8,000 genes [verify]). | **No.** The paper's data statement: "Aligned sequencing data and processed single-cell populations are available on GEO at GSE264667". No pseudobulk or normalised-bulk h5ad was found (the K562 and RPE1 `*_normalized_bulk_01.h5ad` files in `data/perturbseq` came from the Replogle figshare, which has no HepG2/Jurkat counterpart that I could find). A pseudobulk must be built from the single-cell h5ad (the project's `analysis/perturbseq_bulk.R` exists for that). | Feasible but heavy (15 GB download, a build at tens of GB RAM). Weaker on independence: same lab, essential-gene panel family, pipeline and ribosome-dominated truth as consumed K562/RPE1. | accession, files, sizes, absence of pseudobulk on GEO verified; counts not verified |
| X-Atlas/Orion (Xaira; HCT116 and HEK293T, genome-wide FiCS Perturb-seq, about 8 M cells) | Figshare+ DOI 10.25452/figshare.plus.29190726 (also Hugging Face `Xaira-Therapeutics/X-Atlas-Orion`), CC BY-NC-SA 4.0. Files: `HCT116_filtered_dual_guide_cells.h5ad` **209.4 GB**, `HEK293T_filtered_dual_guide_cells.h5ad` **350.2 GB**, guide tables. | not verified | Not in the file list | **Not feasible** (over 1 GB by two orders; single cells only). | file list verified via figshare API |
| Zhu et al., genome-scale CRISPRi Perturb-seq in primary human CD4+ T cells (bioRxiv Dec 2025; GEO GSE314342) | GSE314342; supplementary `RAW.tar` **159 GB** plus small bulk-RNA-seq validation CSVs. | 22 M cells, four donors, three stimulation conditions (from the abstract, [not verified in the data]) | Only small validation bulk tables seen | **Not feasible.** | listing verified; contents not verified |

I did not search exhaustively (scPerturb collection, Frangieh melanoma, Tian iPSC-neuron and others were not checked).

### B2. Human complex truths: is anything untouched?

Grep of `docs/`, `analysis/`, `tools/` (16 Nov): **hu.MAP never appears as a project truth.** The only mentions are as an annotation source column inside `data/human_truth/complexes.tsv`. Complex Portal human appears only in `complex_confirm_design.md` ("available; downloaded to check and then deleted"), never as a truth or scored.

**Caveat found:** the consumed human "CORUM" truth (`analysis/human_truths.R::corum_sets`) is `data/human_truth/complexes.tsv` filtered to rows whose `sources` contain "CORUM". All 2,734 rows do. The same file also tags 318 of them as ComplexPortal, 104 hu.MAP2, 51 hu.MAP3 and 46 hu.MAP **[verified]**. So Complex Portal human and hu.MAP are **not untouched content-wise**: part of what they contain is already inside the consumed CORUM truth.

| Resource | Obtainable | Verdict |
|---|---|---|
| **EBI Complex Portal human (taxon 9606)** | `https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv`, 5,362,905 B, Last-Modified 2026-10-17, 2,533 complexes; **downloaded** to `data/complex_portal_human_raw/9606.tsv` (SHA-256 `f9e977e0cce4f518de066eb267c70cabc2f8c1e659cfc1c164fbee41d4380d53`) with a reviewed UniProt accession-to-gene map (UniProt release 2026_04, `uniprot_9606_reviewed.tsv`, 585,523 B). Never scored. | Usable. Partly overlaps CORUM (below). |
| hu.MAP 3.0 | `https://humap3.proteincomplexes.org/download`: complexes file `hu.MAP3.0_complexes_wConfidenceScores_total15326_wGenenames_20240922.csv` 1.27 MB; 15,326 predicted complexes. **Not downloaded.** | **Not recommended**: machine-learning clusters of mass-spectrometry/co-fractionation networks, and its training sets are CORUM and (reduced) Complex Portal (the download page lists `ComplexPortal_reduced_20230309.train/test`), so it is not independent of either label set [from the download page; the CORUM part from the hu.MAP2 literature, not verified here]. |

**Counts (Complex Portal human, merged), on a PROXY universe**: the 8,561 gene symbols of the consumed K562-essential pseudobulk (`data/perturbseq/K562_essential_normalized_bulk_01.h5ad`, names only). This is a proxy; the HepG2/Jurkat universe needs the matrix (what is needed: the gene list of the dataset's pseudobulk, i.e. the 5.6 or 9.4 GB h5ad `var`, then rerun `count_human_complex_portal.R`). The consumed CORUM merged truth is rebuilt on the same universe by the project's own `corum_sets` (636 complexes, close to the 634 of the design doc).

| Truth (proxy universe, 8,561 genes) | Class | Merged complexes | Genes | Positive pairs | n_eff | Largest share |
|---|---|---|---|---|---|---|
| Complex Portal 9606 | all | **329** | 1,949 | **17,050** | **31.4** | 0.075 |
| | 3-6 | 216 | 742 | 1,408 | 158.4 | |
| | 7-19 | 92 | 792 | 4,622 | 58.5 | |
| | 20+ | 21 | 548 | 11,388 | 14.4 | 0.112 |
| CORUM merged (consumed, same universe) | all | 636 | 2,262 | 34,974 | 9.26 | 0.249 |
| | 20+ | 21 | 652 | 25,828 | 5.02 | 0.338 |

Overlap of human Complex Portal with the consumed CORUM truth: gene Jaccard 0.509; positive pairs shared 10,225 = **60.0% of Complex Portal's positives** and 29.2% of CORUM's; pair Jaccard 0.245; 172 of 329 Complex Portal complexes have a CORUM complex at Jaccard 0.5 or more, **157 do not**. So human Complex Portal is considerably more independent of the consumed truth than CYC2008 is of Complex Portal yeast (40% novel positives against 10.5%), and has about three times the n_eff of CYC2008-on-Kemmeren. Its data side is the weak part (same lab as consumed K562/RPE1; a build of 15 GB; a hepatoma and a T-cell line).

## C. Variance inputs and power (existing files only)

### C1. Yeast exploratory estimation (net4, Complex Portal main truth), from `results/complex_confirm/`

Complexes by size class on the net4 truth (`truth_counts.csv`; n_eff from the same builder, recomputed and equal to the frozen 14.48 for all): 262 merged complexes in total.

| Class | Merged complexes | Positive pairs | n_eff | Share of positives |
|---|---|---|---|---|
| 3-6 | 198 | 1,252 | 144.1 | 9.3% |
| 7-19 | 51 | 2,781 | 36.3 | 20.7% |
| 20+ | 13 | 9,531 | 7.39 | 70.9% |
| all | 262 | 13,438 | 14.48 | |

Estimates and SEs of a3r minus attention (`estimate.csv`, `estimate_final.csv`, `strata.csv`; exploratory, already published in the project):

| Row | d | SE (gene delete-half, 500 reps) | n_eff | k = SE x sqrt(n_eff) |
|---|---|---|---|---|
| main (all) | +0.0962 | 0.01095 | 14.48 | 0.0417 |
| 3-6 stratum | -0.0037 | 0.00386 | 144.1 | 0.0463 |
| 7-19 stratum | -0.0034 | 0.00285 | 36.3 | 0.0172 |
| **20+ stratum** | **+0.1374** | 0.01578 | 7.39 | 0.0429 |

- **Total interval resampling (as fixed in the design):** SE_total = sqrt(SE_gene^2 + SE_series^2), interval d +/- 1.96 SE_total. SE_gene = SD of 500 paired delete-half replicates over the primary truth's genes (seed 20261114), SE_series = SD of 100 paired delete-half replicates over the 59 `#Experiment` series with all methods refitted (seeds 20264000 + b). The group-level (321 groups) version is reported beside it. Observed on main: SE_gene 0.01095, SE_series 0.00978, SE_group 0.00638, so SE_total (series) = 0.0147 (interval [+0.0675, +0.1250]); **rho = SE_series / SE_gene = 0.893** (group 0.582). The design had projected 1.78 from DREAM5 Y-M; the observed ratio on this study is half of that.
- **Complex-level SD of per-complex differences** (`per_complex_ap.csv`, a3r AP minus attention AP, complexes with 10 or more positive pairs; the only complexes with an AP): 125 complexes, mean +0.0018, SD **0.0306** (3-6: n 61, mean -0.0035, SD 0.0373; 7-19: n 51, mean +0.0017, SD 0.0105; **20+: n 13, mean +0.0273, SD 0.0376**). Median per-complex difference is about 0 in every class; the 20+ mean comes from five complexes (CPX-1601, 1599, 1602, 2262, 1604: +0.072, +0.057, +0.077, +0.106, +0.028; the first, second and fourth-listed ribosome subunits are CPX-1601 and 1599 cytoplasmic and 1602 mitochondrial per the design doc, the nature of CPX-2262 was not checked); CPX-1603 is +0.017 and the other seven 20+ complexes are within 0.001 of 0.
- **Contribution of the 20+ complexes:** they carry 71% of positives. In the positive-weighted mean of per-complex differences (+0.0330), the 20+ complexes contribute +0.0327, 3-6 contribute -0.0003 and 7-19 +0.0005. In the pooled estimate, the 20+ stratum gives +0.137 against -0.004 and -0.003 for the others (the three strata estimates are for strata universes defined in the design; they do not add to the pooled d). **The whole pooled effect is the large-complex (ribosome-dominated) stratum.** On the noribo variant the pooled d is +0.0455 (SE 0.0106).
- **Complex-level delete-half jackknife** (labelled sensitivity in the design, 500 reps): SE_C = **0.0608**, interval [-0.023, +0.215]. This is 5.6 times the gene-jackknife SE, because the effect rides on a handful of large complexes. The preregistered intervals use the gene jackknife (genes are the right unit for overlapping complexes) and add the experiment unit; SE_C shows how the answer would change if complexes were the unit. A delete-half SD approximates the full-sample SE (design: "jackknife factor 1").

### C2. What is detectable (one-sided 0.025, 80% power)

Required SE for a true difference d: SE <= d / (z_0.975 + z_0.80) = d / (1.9600 + 0.8416) = d / 2.8016.
- d = +0.03: SE <= 0.01071.
- d = +0.05: SE <= 0.01785.

**[assumed] k rule** (design doc): SE_gene = k / sqrt(n_eff), with k observed on net4 = 0.0417 (all) and 0.0429 (20+ stratum) [verified from the files above], inside the 0.037-0.041 of the four earlier anchors. The experiment component is added as SE_total = SE_gene x sqrt(1 + rho^2). rho = 0 (gene only), 0.893 (net4 series, observed) and 1.78 (design projection). **rho for Kemmeren is unknown** (no series structure; see A1.4), so all three are shown. Required n_eff = (k sqrt(1 + rho^2) / SE_max)^2.

| Scope | rho | d = +0.03: n_eff needed | d = +0.05: n_eff needed |
|---|---|---|---|
| Pooled (all complexes), k = 0.0417 | 0 | (0.0417 / 0.01071)^2 = **15.1** | (0.0417 / 0.01785)^2 = **5.5** |
| | 0.893 | 15.1 x 1.797 = **27.2** | 5.5 x 1.797 = **9.8** |
| | 1.78 | 15.1 x 4.17 = **63.3** | 5.5 x 4.17 = **22.8** |
| 20+ stratum, k = 0.0429 | 0 | (0.0429 / 0.01071)^2 = **16.0** | (0.0429 / 0.01785)^2 = **5.8** |
| | 0.893 | **28.8** | **10.4** |
| | 1.78 | **67.1** | **24.1** |

Converting n_eff to complex counts at the observed ratio **[assumed constant]**: all complexes 262 / 14.48 = 18.1 complexes per n_eff (so 274 / 99 complexes for gene-only at +0.03 / +0.05, 492 / 177 at rho 0.893); 20+ stratum 13 / 7.39 = 1.76 per n_eff (28 / 10 complexes gene-only, 51 / 18 at rho 0.893). By the complex-level SE (SE_C 0.0608 at 262 complexes, SE scaling with 1 / sqrt(complexes)): 262 x (0.0608 / 0.01071)^2 = **8,450 complexes** for +0.03 and 262 x (0.0608 / 0.01785)^2 = **3,040** for +0.05, i.e. impossible; the pooled test only exists on the gene-jackknife convention.

**What each candidate then gives** (power = Phi(d / SE - 1.96), SE = k sqrt(1 + rho^2) / sqrt(n_eff)):

| Candidate (n_eff) | rho | SE | MDE at 80% | Power at +0.03 | Power at +0.05 |
|---|---|---|---|---|---|
| Kemmeren + CP (15.6, consumed truth) | 0 / 0.893 / 1.78 | 0.0105 / 0.0141 / 0.0216 | 0.030 / 0.040 / 0.060 | 0.81 / 0.56 / 0.28 | 1.00 / 0.94 / 0.64 |
| **Kemmeren + CYC2008 (10.0)** | 0 / 0.893 / 1.78 | 0.0132 / 0.0176 / 0.0269 | 0.037 / 0.049 / 0.075 | **0.63 / 0.40 / 0.20** | **0.97 / 0.81 / 0.46** |
| Kemmeren + union (15.3) | 0 / 0.893 / 1.78 | 0.0107 / 0.0143 / 0.0218 | 0.030 / 0.040 / 0.061 | 0.80 / 0.55 / 0.28 | 1.00 / 0.94 / 0.63 |
| CYC2008 20+ stratum (4.43) | 0 / 0.893 / 1.78 | 0.0204 / 0.0273 / 0.0417 | 0.057 / 0.077 / 0.117 | 0.31 / 0.19 / 0.11 | 0.69 / 0.45 / 0.22 |
| CP 20+ stratum on Kemmeren (7.54) | 0 / 0.893 / 1.78 | 0.0156 / 0.0209 / 0.0319 | 0.044 / 0.059 / 0.090 | 0.48 / 0.30 / 0.15 | 0.89 / 0.67 / 0.35 |
| Human CP, proxy (31.4); gene SE only, pseudobulk has no experiment unit | 0 | 0.0074 | 0.021 | 0.98 | 1.00 |

Reading: the pooled test on CYC2008 is adequate for +0.05 and for the +0.096 seen on net4 (power about 1 for all rows), marginal for +0.03, and the 20+ stratum cannot be tested on its own at either size. The union (n_eff 15.3) would power +0.03 only in the gene-only case, but it is 92% consumed. The human figure uses the yeast k [assumed to transfer; the K562/RPE1 anchors had k 0.038-0.041] and the proxy truth/universe; treat it as indicative. Because the net4 effect lives entirely in the 20+ stratum (and the small/mid strata showed -0.004, -0.003), a Kemmeren pooled result would mostly re-measure the ribosome-dominated stratum again, in a truth where the 20+ stratum has 9 merged complexes and n_eff 4.4.

## D. Cost estimate on Kemmeren (32 threads, `nice 19`, through `clr-sched`)

Matrix A = 6,169 genes x 1,283 profiles (all, on disk); matrix B = 6,169 x 652 (responsive only). Basis timings are measured values from the project (net4 = 5,950 x 536 at 16 threads, `results/dream5/SCORE_MANIFEST.tsv`; RPE1 = 8,745 x 1,023, `docs/plans/scoreboard_competitors_plan.md` and `results/scoreboard/complex_competitors/`). Scaling factors are stated; they are estimates, not measurements on Kemmeren.

| Method | Basis (measured) | Scaling to Kemmeren | Wall time at 32 threads | Peak memory | Over 2 h? |
|---|---|---|---|---|---|
| \|r\| | net4: 18 s (16 threads) | cor of 6,169 x 1,283, a few seconds plus load | **under 1 min** | under 1 GB | no |
| CLR (MI + Stouffer CLR) | net4: 21.7 s at 16 threads; M3D rule "3.6 core-min per refit at 4,297 x 466", cost proportional to genes^2 x samples | genes^2 x samples ratio A/net4 = 1.075 x 2.394 = 2.57 (B: 1.075 x 1.216 = 1.31) | **about 0.5-1 min** (A) | 3-5 GB (net4 refit peak 3.33 GB; the G x G matrices scale by 1.075) | no |
| aCLR (adds operator and walk) | net4: 51.8 s at 16 threads | MI part as above, walk about genes^3 (1.11) | **about 1-2 min** (A) | 4-6 GB | no |
| tcaCLR (adds log, second CLR on the same operator) | net4 a3r adds about 1.5 GB over the operator (design Gc table) and seconds | as aCLR | **about 1.5-3 min for CLR + aCLR + tcaCLR together** (A) | about 6-8 GB (declare 10 GB) | no |
| PIDC (`tools/pidc_chunked.jl`, bit-identical to the package) | RPE1 real run: 8,745 genes, 32 threads, **1,160 s total**, peak **10.48 GB** (read 178 s, node pairs 311 s cumulative, puc stage 845 s) | puc stage cubic in genes: (6,169 / 8,745)^3 = 0.351 x 845 s x about 1.1 for more Bayesian-block bins at larger N, about 330 s; read and node pairs scale about with genes x samples and genes^2 | **about 10-15 min** (A and B similar) | **about 5-6 GB** (declare 9 GB; specific-information store scales with genes^2 x bins) | no |
| GENIE3 (R package, 1,000 trees, K = sqrt, all genes regulators) | net4 (all genes): **12,435 s at 16 threads = 55 core-h**; RPE1: 16 targets in 180.5 s on 8 workers, extrapolated peak 14 GB at 32 workers (declared 20 GB) | per target cost about samples x log(samples) x sqrt(regulators): A 2.78, targets 1.037 so x2.88 = **about 159 core-h**; B x1.32 = about 73 core-h. RPE1 basis gives 167 core-h for A | **A: about 5.0 h (4.9-5.2). B: about 2.3 h** | about 12-16 GB (mostly 32 worker copies; G x G output 0.3 GB) | **YES for A and B (PI must be asked)** |
| GRNBoost2 (`tools/grnboost2.py`, legacy early stopping, seed 20261016) | net4 (all genes): **11,764 s at 16 threads = 52 core-h**; RPE1 estimate 130 core-h (range 100-200), 52 core-s per target, trees per target vary 330-840 | per target about samples x regulators x trees: from net4 x2.6 = **135 core-h**; from RPE1 x0.62 = 81 core-h; B: 68 core-h (net4 basis) or 41 core-h (RPE1 basis); tree-count variation (about 0.7-2.5) not included | **A: about 2.5-4.2 h (range 3-10 h with tree-count noise). B: about 1.3-2.1 h** | about 1-2 GB for the edge list at this size (RPE1 measured 0.4-0.7 GB for subsets; the full edge-list peak was scaled in the RPE1 plan section 8b) | **YES for A (PI must be asked); B borderline (1.3-2.1 h)** |

**Flags.** GENIE3 and GRNBoost2 on the full matrix A are each above 2 hours (about 5 h and about 2.5-4 h), and GENIE3 on B is about 2.3 h; none should start without the PI's approval. Everything else is minutes. Cheaper alternatives, **different definitions, not run, for the PI to consider**: (a) matrix B (half the profiles); (b) all methods on the truth genes only (about 1,324 genes for CYC2008; the K562 precedent `baselines.R --gene_subset truths`, which would need CLR/aCLR/tcaCLR/\|r\| rescored on the same subset to stay like-for-like; GENIE3 on 1,324 targets would be about 17 core-h, about 0.5 h at 32 threads, by the same scaling).

Experiment-level resampling (refits): net4 refit with all of CLR, aCLR and tcaCLR took 57 s at 8 threads and 3.3 GB per refit; on Kemmeren about 1-2 min per refit at 8 threads, about 5 lanes in parallel as in job 6030. 100-200 refits are a few hours of wall time, not over 2 h per job if chunked, but the unit (mutants versus sub-series hs1990/hs1991) is a design decision [not made here].

## E. New files and disk (for `docs/disk_budget.md`; that file was not edited)

| Path | Size | Content |
|---|---|---|
| `data/cyc2008_raw/CYC2008_complex.tab` | 108,324 B (dir 112 KB) | Internet Archive capture 2020-02-16 of the authors' file; SHA-256 above |
| `data/complex_portal_human_raw/` | 5.7 MB | `9606.tsv` (5.36 MB), UniProt human reviewed map (0.59 MB), two small header files |
| `results/prereg_candidates/` | under 100 KB | six counts-only CSVs plus one deleted-gene membership CSV |
| `data/yeast` (already counted) | 2.1 GB | Kemmeren matrices already on disk; nothing new |

Not downloaded because over 1 GB (size and URL only): HepG2 h5ad 5.6 GB, Jurkat h5ad 9.4 GB, GSE264667_RAW.tar 16.3 GB (`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/`); X-Atlas/Orion HCT116 209 GB and HEK293T 350 GB (figshare.plus 29190726); GSE314342_RAW.tar 159 GB.

## F. Sources (retrieved 2026-11-16)

- GEO GSE42528, GSE42527, GSE42526, GSE45115, GSE42215 (series text pages and FTP listings at `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/`; local series matrices in `data/yeast`).
- Kemmeren et al. 2014, Cell 157:740, PMID 24766815 (via GEO summaries).
- Pu et al. 2009, "Up-to-date catalogues of yeast protein complexes", NAR 37:825, `https://pmc.ncbi.nlm.nih.gov/articles/PMC2647312` (counts: 408 complexes, 1,630 genes; YHTP2008 400 complexes; site `http://wodaklab.org/cyc2008/`).
- Internet Archive capture of `http://wodaklab.org/cyc2008/resources/CYC2008_complex.tab`, 2020-02-16.
- EBI Complex Portal complextab `https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv` and `.../559292.tsv`.
- UniProt REST (release 2026_04, 2026-11-04), reviewed human proteins.
- hu.MAP 3.0 download page `https://humap3.proteincomplexes.org/download`.
- GEO GSE264667 (Nadig et al. 2025, Nat Genet, PMID 40259084; Europe PMC full text of preprint PMC11244993 for the data statement).
- X-Atlas/Orion, figshare API article 29190726; GEO GSE314342.
- Project files: `docs/plans/complex_confirm_design.md`, `docs/plans/reciprocity_heldout_design.md`, `docs/plans/scoreboard_competitors_plan.md`, `docs/preregistration_reciprocity_yeast.md`, `docs/reciprocity_yeast_selections/step4_gateF_fail_report.md`, `results/complex_confirm/*`, `results/dream5/SCORE_MANIFEST.tsv`, `results/scoreboard/complex_competitors/*`.
