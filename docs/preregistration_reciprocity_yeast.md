# Preregistration: reciprocal attention and headline replication on a held-out yeast deletion compendium

> **STATUS: DRAFT rev 3, under review (1 Nov 2026). Not fixed.** This file becomes a preregistration only when this block is replaced by "*Fixed YYYY-MM-DD*" in a commit. Every truth-touching gate in `analysis/yeast_prep.R` and in the test script refuses to run while this file contains the string `STATUS: DRAFT`.

*Drafted by Claude (Opus, high effort, statistics lead).*

- *Rev 1 (18babdb) came from the design memo `docs/plans/reciprocity_heldout_design.md` and the PI decisions of 1 Nov 2026.*
- *Rev 2 answers the independent review `docs/plans/yeast_prereg_review.md` (fd9c827): three blockers (B1–B3), should-fix items S1–S9 and notes N1–N10. Where the review offered options, the PI chose them on 1 Nov; see "Changes from rev 1" at the end.*
- *Rev 3 records the PI decisions of 1 Nov 2026 on the four flags raised after rev 2 (see the disclosure rows "Primary depth rule" and "Flags resolved").*
- *No yeast or HepG2/Jurkat truth has been downloaded, opened or scored.*
- *No reciprocity score has been computed on any real dataset with any truth.*
- *After this file is fixed, any change is logged in the change log with its date and reason. A change made after the first truth-phase run is a deviation and is reported as one.*

## Question

1. **Reciprocity.** Does reciprocal attention improve recovery of small protein complexes over attention (K10) on data and truths never used in this project, without a regulon cost? Reciprocal attention keeps attention from i to j only when each gene is among the other's softmax candidates, and leaves unreturned mass on i.
2. **Headline replication (H1).** Does attention beat CLR on regulon co-membership in a second organism?

## Disclosure: what has been touched, and what was decided after seeing data

| Item | Status |
|---|---|
| **Reciprocity evidence (synthetic only)** | Round 2 (`analysis/synthetic_width/`, design memo section 8) used the round-2 construction: all genes in the universe and all non-complex pairs as negatives, at t = 9. At partner-less fractions of 5–50% and batch ×0, ×1, ×2, R vs K10 gave: small complexes ×1.86–6.24 (4/4 seeds, all 12 standard cells); regulons ×1.05–1.09; regulon-heavy ×1.02–1.04; TF–target edges ×1.01–1.02; ribosome ×0.80–0.83; all complexes 0.94–1.04. **Round 1, clean scenario** (Gaussian noise, no global programs; not covered by round 2): reciprocity cost regulons **×0.95 (0/4 seeds)**. The variant was added after one exploratory seed. |
| **Synthetic rescore with this file's construction** (review S1; run 1 Nov after rev 1, synthetic only) | The 56 round-2 runs were regenerated from their seeds and rescored with the P1 and S1 constructions below, at t = 9 and at each operator's criterion-A t\* on lab halves (`analysis/synthetic_width/sim_p1_rescore.R`, `out/p1_rescore/`). Results in the "Synthetic rescore" section. **That section is the basis of the planning effect**, not the round-2 ratios. |
| **CLR vs K10 in the synthetic data** | In the round-2 construction CLR beats K10 heavily on small complexes (e.g. `std_f10_b1_s1`: CLR `small_x` 1,877, K10 1,092, R 2,089). Much of R's synthetic gain over K10 restores CLR-level small-complex performance. The P1-construction comparison of R vs CLR is in the rescore section. R vs CLR on complexes is **descriptive only** here. This file does not claim that R beats CLR on complexes. |
| **Rule change (all complexes)** | Reciprocity **failed** the pre-set carry-forward rule on condition C2 (all complexes > 1.00 with ≥ 3/4 seeds), because of the large complexes. Testing it means treating C2 as non-binding, given the PI's choice of small complexes as the target. That decision was made **after** the round-2 tables were seen. Here all complexes are tested first among the key secondary tests (S1, non-inferiority), not as part of adoption. |
| **Estimand change** | The round-2 metric and P1 differ (universe, negatives, depth). P1's construction was written in rev 1 before the rescore was run. The rescore was run after rev 1 and before this revision. **No endpoint definition was changed after the rescore was seen.** Only the planning effect and the disclosure were updated from it. |
| **Gate variants** | The reproducibility gate, gate plus reciprocity and the empirical-null gate failed round 2. They are not tested. The round-1 gate gains came from a generator with 51% partner-less genes. |
| **Serial testing** | This is the **third** complex-oriented operator test, after keff (`docs/keff_verdict/`: headline unchanged) and CLR(log Pᵗ) (`docs/clr_logpt_verdict/`: not adopted). **It is declared the last.** Whatever the verdict, no further budget- or readout-level operator change aimed at complexes is tested confirmatorily in this project. |
| **Timing (PI decision, 1 Nov 2026)** | The design memo made "runs after the Collins readout" a condition. **The PI relaxed it on 1 Nov 2026:** the Collins message goes out on Tuesday at the earliest, and this test may run before it. |
| **K562 in-sample** | About 20 operators and readouts were scored on K562 CORUM, including `d_geomean_t9` (√(PᵢⱼPⱼᵢ), a mutual-attention *readout*, not this operator; `analysis/corum_exploratory/metrics_exploratory.csv`). The drafter saw that row. K562 shaped the PI's focus on complexes. K562 is descriptive here. |
| **No R score on a real truth before the yeast phase** | No reciprocity score is computed against any real truth (K562, *E. coli*, human or yeast) before the yeast truth phase (protocol step 10). Label-free R quantities (mass removed, criterion-A curves) are allowed. |
| **E. coli** | M3D and EcoMAC truths are consumed (headline, keff, CLR(log Pᵗ)); descriptive only (regulon no-harm). **EcoMAC overlaps M3D in 828 of its 2,198 arrays (37.7%)**, so all *E. coli* evidence (M3D and EcoMAC) counts as **one compendium**. A pass of H1 on yeast would therefore be the **first organism-independent replication** of the headline. |
| **Primary depth rule (PI decision, 1 Nov 2026)** | The synthetic rescore showed that much of R's own-t\* advantage is a depth effect, and that R's deeper t\* costs regulons under strong batch. **The primary depth rule (each operator at its own label-free t\*, for P1, P2 and H1) was reviewed after the rescore and kept unchanged by PI decision on 1 Nov 2026 (option A).** A same-depth contrast (both operators at t\*_K10) was added as a fixed, pre-specified **descriptive** secondary with no inferential role (see "Depth"). |
| **Flags resolved (PI decision, 1 Nov 2026)** | Kept as written: (i) the S3 "stop" branch (S2–S6 are tested only if S1 is rejected, so under verdicts 3a and 3b the human arms are descriptive); (ii) the power assumption ρ = SE_B/SE_J ≈ 0.44 is disclosed as is, with the ρ = 1 and ρ = 2 sensitivity rows (power at +0.010: 0.50 and 0.24); (iii) the gate F thresholds (median R′_mut > median R′_wt; Q ≤ 0.5) stay as judgment calls, as the reviewer noted. |
| **Yeast is single-use** | The deleteome expression and every yeast truth below are untouched in this project. Once the truth phase runs, yeast is consumed for good. The a3r follow-up (interpretation memo section 7.3) and the C3 asymmetry edge test (`docs/plans/asymmetry_expert.md`) are **not** included and therefore can never use these yeast truths (PI decision). A failure of gate F or of the power gate does **not** consume yeast, because no truth is opened before them. |
| **Truth–data independence (review B1)** | YEASTRACT's expression evidence is mostly wild-type vs TF-mutant microarrays and may include the deleteome itself or its constituent Holstege-lab datasets. The primary regulon truth therefore requires expression evidence from a reference outside a fixed exclusion list, or falls back to DNA-binding evidence alone (see "Truths"). A contamination audit is committed before scoring. |
| **What the drafter has seen of the yeast data** | The column headers and the first 1–3 rows (probes SCAB000001–3; M and p for about 10 mutants) of `GSE42527_final_GeneExpressionMatrix.txt.gz` and `GSE42526_final_GeneExpressionMatrix.txt.gz`; 35 rows of one wild-type sample table (GSM630359); GEO series metadata. That was to verify formats. No truth was seen. |
| **Human (secondary)** | HepG2 and Jurkat expression (GSE264667) has never been downloaded. CORUM and CollecTRI were used on K562 (in-sample) and RPE1 (CORUM descriptive in CLR(log Pᵗ), S5). |
| **Fixed without labels on M3D** | keff 10, cap 50, α = 0.5, B-spline MI (`bins = "hg"`), Stouffer CLR. Unchanged. |

## Data

### Expression (primary, confirmatory): Kemmeren et al. 2014 deleteome (*Cell* 157:740; PMID 24766815)

- **Design** (GEO overall design, verified 1 Nov 2026):
  - BY4742 deletion strains in synthetic complete medium, mid-log;
  - two independent cultures per mutant, on two two-channel arrays against one common wild-type reference RNA, with a dye swap;
  - each gene spotted twice, so four measurements per gene per mutant;
  - wild-type cultures grown in parallel with the mutants every day.
- **Processed mutant data** (GEO SuperSeries GSE42528):
  - `GSE42527_final_GeneExpressionMatrix.txt.gz` (63.8 MB): **700 "responsive" mutants**, defined as ≥ 4 significant changes against wild type;
  - `GSE42526_final_GeneExpressionMatrix.txt.gz` (67.5 MB): 784 non-responsive mutants, 6,182 probes × (M, p) pairs;
  - per mutant: M (limma average log2 mutant/reference) and a p-value (limma 2.12.0, Benjamini–Hochberg-adjusted; `GSE42527_protocols.txt.gz`);
  - platform GPL11232 (Holstege lab); probe IDs `SCAB######`.
- **Wild-type replicate controls (verified):**

  | Series | What it is | Hybridizations |
  |---|---|---|
  | GSE45115 | 376 wild-type cultures grown in parallel with the mutants | 755 |
  | GSE42215 | wild-type pool, Erlenmeyer protocol | 200 |
  | GSE42217 | wild-type pool, plate-reader protocol | 200 |
  | GSE42240 / GSE42241 | MATa wild types | 8 / 20 |

  - Per-hybridization VALUE = normalized log2 ratio test/ref, in the GEO sample tables (checked for GSM630359) and the `_RAW.tar` archives.
  - The GSE45115 series-matrix table is empty, so the wild-type matrix is assembled from the sample tables.
- **What the controls do here (stated plainly, review S2).**
  - (i) They are the reference against which the authors called mutants "responsive", so the primary sample selection rests on calibrated significance, not on effect size alone.
  - (ii) They supply a wild-type-only network for a **sanity check** (gate F) that the mutant network is not dominated by technical or day-to-day structure.
  - (iii) They give descriptive per-gene wild-type variation.
  - They do **not** calibrate MI, CLR or the operators, and F is not an inferential estimability test.
- **Primary sample set: the 700 responsive mutants (GSE42527).**
  - Rationale, fixed now and label-free: profiles with three or fewer significant changes are wild-type-like as columns.
  - The K562 compendium likewise used significant knockdowns only, and N ≈ 700 matches K562 network half A (707).
  - All 1,484 mutants: descriptive sensitivity.
- **Preprocessing (`analysis/yeast_prep.R`, stage `matrix`; fixed, label-free; review S6):**
  1. **Probes to ORFs.** Map probes to systematic ORF names through the GPL11232 annotation; average duplicate probes of one ORF.
  2. **Deleted gene.** Parse it from each column label (`<gene>-del-…`) and resolve it to an ORF via `SGD_features.tab`: standard name, then systematic name, then aliases (column 6). **Columns whose deleted gene cannot be resolved to an ORF in the matrix are dropped**, and their labels are logged.
  3. **Masking.** In each mutant's column, mask the deleted ORF and every ORF whose SGD coordinates overlap it or lie within **500 bp** of it on the same chromosome (cassette and neighbouring-gene effects). Replace masked entries with that gene's median over the other mutants.
  4. **Missing values.** Drop ORFs with missing values in more than 5% of profiles (after masking, masked entries do not count as missing); impute the rest with the gene median.
  5. Drop constant rows.
  6. **No batch, slow-growth or SVD correction.** Rationale: the headline (M3D) and K562 compendia were analysed without such correction, and any choice here would be a new fork. The authors' growth-corrected matrix (`deleteome_all_mutants_svd_transformed.txt`) is a descriptive sensitivity only.
- **Gene universe:** the rows that remain (expected about 6,000).
- **Mutant clusters (label-free; used by the experiment jackknife; committed with the selections; review B3).**
  - Complete-linkage clustering of the mutant profiles (columns) on 1 − Pearson r over all universe genes, cut at height 0.5, so every pair within a cluster has r ≥ 0.5.
  - The number of clusters and the cluster-size distribution are committed as a label-free diagnostic.

### Gate F (label-free sanity check with the wild-type controls; computed and committed before any truth is downloaded; review S2)

1. **Matrices.**
   - WT: one profile per GSE45115 culture. Each profile is the mean of its dye-swap pair, oriented to log2(wt/ref), over duplicate spots, with the same probe mapping and gene filter.
   - The orientation rule is chosen by the **implementation check**: the candidate rules "as deposited" and "negate the hybridization whose channel 1 is the test strain", applied to the sample tables of the first 10 GSE42527 mutants (alphabetical by deleted gene), must reproduce the published M with Pearson r ≥ 0.99 and median |difference| ≤ 0.02.
2. **Failed implementation check (fallback, implemented).** If no rule passes (the published M is a limma estimate with gene-specific dye correction, so a simple mean can fail for benign reasons):
   - download the sample tables of a fixed random 376 mutants (seed 20261105);
   - rebuild those mutants with the best rule (highest median r);
   - run F on rebuilt mutants against WT rebuilt the same way, so both arms are on one scale.

   The confirmatory analysis still uses the published M of all 700 mutants.
3. **Statistic.**
   - Take N_w = the number of WT profiles (about 376) columns from the mutants (10 draws, seed 20261199; in the fallback, the 376 rebuilt mutants).
   - Split each draw into random halves and fit Stouffer CLR on each half (`build_fit` as `analysis/lab_splits.R`).
   - R = |top-1% pairs of half 1 ∩ top-1% pairs of half 2| / |top-1% pairs|. The chance-corrected version is R′ = (R − 0.01)/0.99.
   - Do the same on 10 random halvings of the WT profiles.
   - **Overlap share** Q = the share of the full mutant network's top-1% pairs that are also top-1% pairs in the full WT network.
4. **Gate.** F **passes** if median R′_mut > median R′_wt **and** Q ≤ 0.5.
   - Justification of 0.5: if more than half of the mutant network's strongest pairs are also strongest in a network of wild-type replicates alone, the network mostly reflects technical or day-to-day structure, and an operator comparison would mostly measure that.
   - Also reported, not gating: R′_mut, R′_wt and their spread; the number of genes per mutant beyond the gene's WT 0.5th–99.5th percentile range.
5. **If F fails:** verdict 1. No truth is opened, so yeast is not consumed.

### Truths (downloaded only after this file is fixed; structure counted, audited and committed before any scoring)

**Exclusion list E (fixed now; review B1).** E is the PubMed IDs of Holstege/Kemmeren-lab expression studies and the deleteome's constituent datasets.

- **(a) Every PubMed ID linked by GEO to a series on platform GPL11232** (query of 1 Nov 2026; 34 series, 20 PMIDs): 21145464 (van Wageningen 2010, kinases/phosphatases), 21596317 (Lenstra 2011, chromatin), 21976730, 21998594, 22697265, 23599000, 23785440, 23932120, 24009531, 24037263, 24324601, **24766815 (Kemmeren 2014)**, 24952590 (O'Duibhir 2014), 25157590, 26700642 (Sameith 2015), 27979920, 29979971, 31506296, 34520050, 36066422.
- **(b) PubMed IDs of *S. cerevisiae* GEO series with Holstege F or Kemmeren P as author, on other platforms:** 19713935, 21971086, 28049706, 28637236, 30203048.
- **(c) Reanalysis screen** (mechanical, applied by script to every PubMed ID cited as expression evidence): any reference whose PubMed title or abstract contains "deleteome", "GSE42526", "GSE42527", "GSE42528" or "GSE45115" is added to E.

The list is committed as `docs/reciprocity_yeast_selections/excluded_pmids.txt` with this file. E restricts only the **expression** leg. DNA-binding evidence from any source counts.

**Yeast regulons (P2, H1, S6).**

- **Source.** YEASTRACT+, documented regulations, all environmental conditions, all TFs by all genes.
- **Timing of the export.** It is performed **after** this file is fixed. Its timestamp, release, form and parameters are recorded in `docs/reciprocity_yeast_selections/yeastract_export_log.md` and committed before the export is parsed.
- **Mode, decided on availability alone, before parsing, and logged.**
  - **Mode "refs" (primary):** used if YEASTRACT returns, for each association, its evidence type (DNA binding or expression) and supporting PubMed IDs.
    - **Primary truth:** TF–target associations with DNA-binding evidence (any reference) **and** expression evidence from **at least one PubMed ID not in E**.
    - "DNA binding and any expression evidence" is a descriptive comparator, labelled "possibly circular".
  - **Mode "binding_only" (fallback):** used if per-association references cannot be obtained. The primary truth is DNA-binding-evidence associations only (ChIP, ChIP-chip, ChIP-seq, ChIP-exo, EMSA, footprinting, as classified by YEASTRACT), and the paper says so. "DNA binding and expression" is then descriptive only, labelled "possibly circular".
- **Encoding.** Any non-zero matrix entry, including −1 for negative regulation, counts as an association. The value table of the export is logged.
- **Co-membership.** Regulator-agnostic, as `make_regulons()` / `comembership_pairs()` adapted to yeast: regulator = TF ORF, target = gene ORF, every association counted as one evidence class. Regulons of 5–500 in-universe targets.
- **TF–target edges (S6).**
  - Universe: pairs with exactly one regulator gene, where a regulator gene is a TF with a 5–500-target regulon in the primary truth.
  - Positives: primary associations.
  - S6 is non-inferiority only. With the deleted gene masked, the deleteome's causal evidence for X → Y (column X) is invisible to the row-based edge (X, Y), so the edge signal comes only from co-variation of X's mRNA across other columns. Many yeast TFs are regulated post-translationally, and S is symmetric, so direction is not identifiable.
- **Contamination audit (pre-scoring, truth structure only; committed with the counts).** In mode "refs", report:
  - (a) the share of "DNA binding and any expression" associations whose expression evidence comes **only** from E (removed from the primary truth);
  - (b) the share of all expression-evidence associations whose only references are in E;
  - (c) the same shares in co-membership positive pairs.

  In mode "binding_only", report the share of binding associations that also carry expression evidence (without references this cannot be split by source).

**Yeast complexes (P1, S1).**

- **Source.** EBI Complex Portal `complextab/559292.tsv`, the release current on the download date (the 1 Nov listing shows 2026-10-17, 913 KB), with UniProt 559292 accession → ordered-locus (ORF) mapping.
- **Build.** Flattening, non-protein removal and isoform handling exactly as `analysis/build_ecoli_complex_truth.R`; size filter 3–200 in-universe genes; overlap merge (absorb proper subsets, then complete linkage at Jaccard ≥ 0.5). There is no operon rule in yeast.
- **Evidence filter (B1).** Drop any complex whose cited PubMed IDs (complextab cross-reference column) are all in E. Log the count, and record the ECO evidence codes of every complex.
- **P1 truth:**
  - universe: all pairs of genes in some merged complex;
  - positives: co-members of a merged complex with 3–6 in-universe genes (condition C1's size class);
  - excluded: pairs that are co-members only of larger merged complexes;
  - negatives: the rest.
- **S1 truth:** all merged complexes; positives = co-members of any merged complex.
- **Fallback, fixed now.** If n_eff(P1) < 100 with Complex Portal, use the union of Complex Portal and CYC2008 (Pu et al. 2009; ORF lists), merged by the same rule, decided from counts alone before scoring. If CYC2008 cannot be obtained, or the union still has n_eff < 100, the result is verdict 1. The unused source is never scored.
- **Ribosome (descriptive).** Merged complexes with ≥ 30% cytoplasmic r-protein ORFs: the ORFs of `SGD_features.tab` standard names matching `^RP[LS][0-9]+[AB]?$`.

**Human (S2–S5).** The existing CORUM file (`data/human_truth/complexes.tsv`, as `corum_sets()`, then the overlap merge, then 3–6-gene merged complexes as P1) and CollecTRI (`data/human_truth/collectri.tsv`, as in K562).

**Evaluability (computed in `truth_structure`, before scoring).**

| Truth | Requirement |
|---|---|
| P1, S1 | ≥ 200 positive pairs and ≥ 30 merged complexes with positives |
| P2 / H1 | ≥ 200 positive pairs |
| S6 | ≥ 200 positive edges |
| Human truths | ≥ 200 positive pairs |

A primary truth that fails gives verdict 1. A secondary truth that fails is "not evaluable", counts in its family, and is reported.

**Descriptive truth sensitivities:**

- regulons with divergent (shared-intergenic-region) gene pairs excluded (SGD coordinates; adjacent genes on opposite strands with their 5′ ends facing; fixed, review S5);
- pairs on the same chromosome excluded (review N5);
- complex pairs whose two ORFs come from a single UniProt accession removed (paralog cross-hybridization, review S6);
- all-evidence regulons;
- Complex Portal without the evidence filter.

## Operators (fixed)

**Upstream, unchanged.** B-spline MI (`estimate_mi(bins = "hg", transform = "none")`) and Stouffer CLR (`calibrate(method = "normal", combine = "stouffer")`) as `build_fit` in `analysis/lab_splits.R`.

**K10 (comparator).** A = per-row softmax over the top 50 positive CLR scores at perplexity 10 (`clr:::.softmax_rows(C, keff = 10, cap = 50)`). Rows with no positive candidate get a self-loop **before** reciprocity is applied.

**R, reciprocal attention (candidate; no free parameter; review S8).**

- Let Kᵢ be row i's **candidate set**: the up to 50 genes with the largest positive CLR scores, the same set the softmax runs over.
- Write the indicator matrix as ψᵢⱼ = 1[j ∈ Kᵢ].
- Then:
  1. B = A ⊙ (ψ ⊙ ψᵀ), with the diagonal of B set to 0.
  2. Rows of B with no remaining mass get a self-loop.
  3. B_ii ← B_ii + max(0, 1 − Σⱼ B_ij).
- Support is defined by candidate membership, not by numeric positivity of A. The two differ only if exp() underflows inside the cap.
- **Implementation checks** (synthetic only, before any real data are fitted):
  - (i) a unit test of B on a hand-built 6-gene example, including an underflow case;
  - (ii) on the synthetic seed `std_f10_b1_s1`, the RECIP rows of `analysis/synthetic_width/out/round2/runs/std_f10_b1_s1_metrics.csv` are reproduced to a relative 1e-8. The synthetic `recip()` uses positivity, which is identical there because no underflow occurs; that is checked and logged.
  - The memo's 1e-10 tolerance is superseded by 1e-8.

**Walk and readout, both operators.** P = ½(I + B) (α = 0.5). S = (Pᵗ + Pᵗᵀ)/2 with diagonal 0. **Attention readout only.**

**CLR (H1 comparator).** The Stouffer CLR matrix that both operators are built on.

## Depth (label-free, fixed)

- **Splits.** 10 random halves of the mutants that survive preprocessing (seed 20261101; written by `yeast_prep.R`, stage `splits`, and committed). Each split is used in both directions.
- **Criterion A** (as `analysis/keff_test.R`): the AP of S on the training half for the other half's top-1% CLR pairs.
  - Grid {1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32}.
  - For each split-direction, the argmax.
  - For each split, the geometric mean of the two directions.
  - The median over splits, snapped to the grid.
- **Each operator gets its own t\*:** t\*_K10 and t\*_R.
- **Primary comparison:** R at t\*_R vs K10 at t\*_K10.
- **H1:** K10 at t\*_K10 vs CLR. **This is the headline replication** (review N6). K10 at the transferred t = 9 vs CLR is descriptive.
- **Same-depth contrast (fixed, pre-specified, descriptive; PI decision 1 Nov 2026).** R and K10 both at t\*_K10. It is reported for P1, S1 and P2 with the deciding SE's intervals (√(SE_J² + SE_B²), same jackknife halves and cluster replicates), and has **no inferential role**: no p-value, no verdict, and it cannot rescue or overturn a primary result.
- **Other secondary depths (descriptive):** both operators at t = 9.
- **House practice, stated (N2, N3).**
  - t\* is selected on halves of about 350 profiles and applied at the full N.
  - t\* is fixed inside the experiment jackknife, so its selection uncertainty is not propagated.
- **Flags** (reported, not changing anything): t\* at a grid edge; t\*_R > 2·t\*_K10.
- **No truth curve over t** is computed for yeast.
- **HepG2 and Jurkat.** Build the significant-knockdown pseudobulk as `analysis/perturbseq_bulk.R` from `GSE264667_{hepg2,jurkat}_raw_singlecell_01.h5ad`, splitting knockdowns into halves A and B (seed 20261015) as for K562 and RPE1. The network is fitted on half A, and t\* is selected on 10 random halves of half A's knockdowns, as for RPE1 in CLR(log Pᵗ).

## Endpoints

All endpoints are pooled co-membership AUPR, R minus K10 unless stated, at the primary depths.

**Primary (adoption), an intersection–union test (Berger 1982), each test at one-sided α = 0.025:**

| Code | Truth | Hypothesis |
|---|---|---|
| P1 | yeast small complexes (3–6, merged) | superiority: p < 0.025 **and** d ≥ +0.005 |
| P2 | yeast regulons (primary truth, mode as logged) | non-inferiority, margin 0.005 |

**Headline replication, a separate family H at one-sided α = 0.025, tested whatever P1 and P2 show.** The paper never combines family H with the P/S families (N7).

| Code | Truth | Contrast | Hypothesis |
|---|---|---|---|
| H1 | yeast regulons (P2's truth) | K10 at t\*_K10 minus CLR | superiority: p < 0.025 and d ≥ +0.005 |

**Key secondary: a gatekept fixed sequence (review S3), entered only if P1 and P2 both pass.**

1. **S1 first, at the full one-sided α = 0.025:** yeast, all merged complexes, non-inferiority with margin 0.005.
2. **Only if S1 is rejected:** Holm over S2–S6 at one-sided familywise α = 0.025. A test that cannot be run still counts in the family.
3. **If S1 is not rejected (or not evaluable):** S2–S6 are **not tested** and are reported as descriptive intervals without p-values. This is the "stop" branch of the review's S3 option, chosen so that no α is spent twice.

| Code | Data | Truth | Hypothesis |
|---|---|---|---|
| S2 | HepG2 | CORUM small (3–6, merged) | superiority: p < 0.025 (Holm) **and** d ≥ +0.005 |
| S3 | Jurkat | CORUM small (3–6, merged) | superiority: p < 0.025 (Holm) **and** d ≥ +0.005 |
| S4 | HepG2 | CollecTRI | non-inferiority, margin 0.005 |
| S5 | Jurkat | CollecTRI | non-inferiority, margin 0.005 |
| S6 | yeast | TF–target edges | non-inferiority, margin 0.005 |

**Descriptive (intervals only, never p-values):**

- per-complex AP (complexes with ≥ 3 positive pairs);
- relative changes (ratio of AUPRs) for every endpoint (N8);
- the ribosome alone; small complexes of 3–10 genes;
- the same-depth contrast (both at t\*_K10) for P1, S1 and P2, with intervals;
- R vs CLR on P1 and S1; K10 vs CLR on P1;
- the secondary depths;
- all 1,484 mutants; the growth-corrected matrix;
- the truth sensitivities listed above;
- K562 essential half A (CORUM small, CORUM all, CollecTRI; a fresh operator on a used truth);
- *E. coli* M3D and EcoMAC SC regulons (no-harm);
- the label-free diagnostics below.

**Label-free diagnostics** (reported with the selections, before any truth), per dataset including M3D and K562:

- the share of K10 off-diagonal mass that R removes, and its row distribution;
- the share of rows left as pure self-loops;
- the reproducible-partner count from the 5-fold k-means cluster split (design memo, section 8.3);
- the Spearman correlation between the mass removed and that count;
- for yeast, the mutant-cluster count and size distribution.

## Statistics

- **AUPR:** `ap_fast()` / `pr_summary()` in `analysis/regulondb.R`. Tie blocks are counted at their end, as everywhere in the paper; the expected-tie AP is reported descriptively (Supplementary Note 16).
- **Gene jackknife (SE_J).**
  - The paired delete-half gene jackknife, 100 replicates: in each, a random half of the truth's genes is kept, and AUPR is computed over the truth's pairs with both genes kept. The same halves serve every operator.
  - "The truth's genes" means all genes in that truth's universe. For P1 that is all genes in any merged complex, not only small-complex genes.
  - **H1 uses exactly P2's halves.**
  - Seed 20261102, set once per dataset run; truths are drawn in endpoint-table order (P1, P2, S1, S6); never rerun with another seed.
  - SE_J = the SD of the 100 paired differences.
- **Experiment jackknife (SE_B, yeast only; review B3).**
  - The delete-half jackknife over the label-free **mutant clusters**: 200 replicates. Replicate b keeps a random half of the clusters (rounded down; seed 20261103 + b) and all their mutants.
  - MI, CLR, K10 and R are refitted in every replicate at the committed t\* values; truths and gene universe are fixed.
  - With d = n/2 clusters deleted, the delete-d variance factor (n − d)/d is 1, so SE_B = the SD of the 200 paired differences.
  - No duplicate columns arise, and the clustering of near-identical deletion profiles is respected.
  - It is conservative, because half-N networks are noisier than the full-N network.
- **Deciding statistic (yeast P1, P2, H1, S1, S6).** SE = √(SE_J² + SE_B²).
  - superiority: p = 1 − Φ(d/SE);
  - non-inferiority: p = 1 − Φ((d + 0.005)/SE);
  - harm, for reporting a significant loss: p = 1 − Φ(−d/SE) for d < 0, together with d ≤ −0.005.
  - In a crossed genes × experiments model, the sum counts the interaction component twice (Owen 2007), so it is conservative and does not need independence.
- **Human tests (S2–S5).** SE = SE_J × √(1 + 0.44²), a pre-fixed inflation for the experiment component that is not estimated there (the ratio observed in M3D for CLR(log Pᵗ)). They are labelled accordingly.
- **Individual-endpoint claims outside the IUT** (for example "P1 superiority shown" under verdict 5, or "a significant regulon loss" under verdict 4) use Holm over {P1-superiority, P2-NI, P2-harm} at one-sided 0.025 (review S9).
- **No pooling** across datasets or across families.
- **Smallest effect of interest:** 0.005 AUPR (absolute). Relative changes are reported descriptively.

## Synthetic rescore and planning effect (review S1)

*Run 1 Nov 2026 after rev 1, synthetic only (scheduler jobs 0970–0971; `analysis/synthetic_width/sim_p1_rescore.R`, `summarize_p1.R`; outputs `out/p1_rescore/summary_p1.csv`, `tstar_by_cell.csv`).*

**Method.**

- All 56 round-2 runs were regenerated from their seeds. The regenerated K10 and R `small_x` match the stored round-2 metrics to a relative 4e-15.
- The runs were rescored with the P1 and S1 constructions of this file, and with round-2 regulon co-membership (REG).
- Depths: t = 9 and each operator's own criterion-A t\* on the lab halves (labs 1–5 vs 6–10).

**Results.** Geometric-mean ratios over 4 paired seeds; d is the mean absolute AUPR difference. Twelve standard cells (f = 5–50%, batch ×0–×2) and two regulon-heavy cells.

| Contrast | P1 (small complexes) | S1 (all complexes) | REG (regulons) |
|---|---|---|---|
| R vs K10, t = 9 | ×1.03–1.64 (median ×1.05), d +0.020 to +0.103; 4/4 seeds in all 14 cells | ×0.99–1.01, d −0.011 to +0.004 | ×1.02–1.09, 4/4 in all 14 cells |
| R vs K10, own t\* (**the primary contrast's form**) | ×1.09–2.41 (median ×1.36), d +0.050 to +0.171; 4/4 in all 14 cells | ×1.00–1.09 | **×0.86–1.16**; at batch ×2, **×0.86–0.96 (0/4 seeds) in all four cells** (d −0.001 to −0.005, at regulon AUPR levels of 0.02–0.04) |
| R (own t\*) vs CLR | ×0.88–1.21 | ×1.01–1.46 | ×0.79–1.69 |
| K10 (own t\*) vs CLR | **×0.39–0.96**, 0/4 to 1/4 seeds | ×1.01–1.39 | ×0.91–1.45 |
| K10 (t = 9) vs CLR, f = 10% | ×0.84–1.10 | ×0.98–1.29 | ×1.00–1.71 |

**What this means, disclosed.**

1. **In the P1 construction (complex genes only) the round-2 picture "CLR ≫ K10 on small complexes" largely disappears at t = 9.** K10 and CLR are close (×0.84–1.10 at f = 10%), and R adds ×1.03–1.64.
2. **Criterion A selects deep depths in these synthetic data:** K10 t\* 10–24, R t\* 16–32, so R's t\* is typically 4–8 steps deeper. At its own t\*, K10 loses heavily to CLR on small complexes (×0.39–0.96).
   - Much of R's own-t\* advantage over K10 (×1.09–2.41) is therefore a **depth effect**: R's unreturned mass stays home, so its walk degrades small complexes more slowly with t.
   - The common-depth contrast (both at t = 9) is the smaller ×1.03–1.64.
   - This is a known property of the primary contrast as written in rev 1. **The contrast was not changed after the rescore;** the PI reviewed it and kept it on 1 Nov 2026. The same-depth contrast (both at t\*_K10) is a fixed descriptive secondary, and both at t = 9 is descriptive.
3. **Regulon risk.** Under strong batch effects (×2) R's deeper t\* costs regulons ×0.86–0.96 (0/4 seeds in all four cells). On a real regulon AUPR of about 0.1, a 14% relative loss would exceed the 0.005 margin. This is the main synthetic risk for P2.
4. **R vs CLR on small complexes** is ×0.88–1.21. R is not reliably better than CLR on P1. This file does not test or claim it.

**Planning effect: +0.010 AUPR (P1, R vs K10 at own t\*).**

- Synthetic absolute differences (+0.05 to +0.17) are **not transferable**: synthetic P1 AUPR levels are 0.14–0.70, while real small-complex AUPRs are expected near 0.05–0.15 (K562 CORUM attention: 0.047, all complexes).
- The planning value is the relative effect at the bottom of the synthetic own-t\* range (×1.09) applied to a K10 level of about 0.10, or ×1.2 at 0.05.
- The test is not powered for the minimum effect of 0.005.

## Power (structure-only gates; computed and committed before scoring)

- **Small-complex arm.**
  - n_eff(P1) = (Σₖ pₖ)² / Σₖ pₖ², where pₖ is the number of positive pairs of merged complex k after merging and restriction to the universe.
  - **SE projection.** SE_J,proj = 0.036/√n_eff. The anchor is the *E. coli* complex arm of CLR(log Pᵗ): SE 0.018 at n_eff ≈ 4 (n_eff there was derived from assumed ribosome block shares, not computed). SE_proj = SE_J,proj × √(1 + ρ²), where ρ = SE_B/SE_J.
  - **Gate:** n_eff ≥ 100, so SE_J,proj ≤ 0.0036.

  | ρ = SE_B/SE_J | SE at n_eff = 100 | MDE (80%) | Power at the planning effect +0.010 | Power at +0.015 |
  |---|---|---|---|---|
  | 0.44 (M3D, CLR(log Pᵗ)) | 0.0039 | 0.011 | 0.72 | 0.97 |
  | 1 | 0.0051 | 0.014 | 0.50 | 0.84 |
  | 2 | 0.0080 | 0.023 | 0.24 | 0.46 |

  - The 0.72 figure assumes the M3D ratio. In a perturbation compendium a small complex's co-expression can hinge on a few deletion columns, so ρ may be 1 or more, and the cluster jackknife is conservative by construction. **The test is therefore adequately powered only under favourable assumptions. This is disclosed, and the realised SE_J and SE_B are reported against these projections.** The MDE is reported, never observed power.
  - **Label-free SE_B proxy** (selection stage, reported, not gating): the SD of the criterion-A AP difference R vs K10 over the first 50 cluster-jackknife replicates.
- **Regulon arm (B2).** n_eff(P2) is computed the same way (Σ over regulons of within-regulon positive pairs) and reported in `truth_structure`, not gating. *E. coli* mild K10 tweaks had paired regulon SEs of 0.0005–0.0022; at 0.0015 and ρ = 0.44, NI power at d = 0 is 0.86. **A failure to show non-inferiority under low power is not evidence of a regulon cost** (verdict 5a vs 5b).
- **If the power gate fails** (after the CYC2008 fallback): verdict 1, and no score is computed.

## Decision rules and verdicts

The verdict is computed by the script and written to `verdict.csv` before anyone interprets the numbers.

| # | Condition | Verdict (one line) | Reporting |
|---|---|---|---|
| — | Synthetic implementation check of R fails (protocol step 5) | **No verdict.** Fix the implementation and log the cause before any real data are fitted | — |
| 1 | Gate F fails, or the power gate fails after the fallback, or a primary truth fails evaluability | **Not run:** no adequately powered held-out test was possible | Reciprocity stays a synthetic-only, explained candidate; the gate numbers are reported; yeast is not consumed if no truth was opened |
| 2 | P1 and P2 pass, and S1 is rejected | **Confirmed, no large-complex cost:** reciprocal attention becomes the method's operator from here on | The *E. coli* headline keeps K10, and the order of events is stated; S2–S6 per Holm |
| 3a | P1 and P2 pass, and S1 is evaluable but not rejected | **Confirmed as a complex-mode variant:** it improves small complexes; non-inferiority on all complexes was not shown | K10 stays the default; if S1 shows a significant loss (harm test, Holm as above), "large-complex cost" is said; S2–S6 descriptive |
| 3b | P1 and P2 pass, and S1 is not evaluable | **Confirmed as a complex-mode variant; the large-complex question could not be evaluated** | as 3a, without any large-complex statement |
| 4 | P1 fails, whatever P2 shows | **Not confirmed** | Wording as in CLR(log Pᵗ): "not confirmed", with estimate and interval, never "no effect"; a significant regulon loss (Holm over {P1, P2-NI, P2-harm}) is added if present |
| 5a | P1 passes, P2 non-inferiority not shown, no significant regulon loss | **Not adopted: small-complex gain with an unresolved regulon question** (non-inferiority not shown) | Not a trade-off claim; the realised P2 SE is reported against the projection |
| 5b | P1 passes and there is a significant regulon loss (d ≤ −0.005, harm p < 0.025, Holm over {P1, P2-NI, P2-harm}) | **Regulon trade-off:** small complexes improve at a regulon cost | Not adopted |

- **Family H, reported alongside:**
  - H1 holds: "the headline replicates in yeast".
  - Otherwise: "attention did not beat CLR on yeast regulons", with the interval. A significant CLR advantage (d ≤ −0.005, harm p < 0.025) is said plainly.
  - H1 never changes verdicts 1–5. If gate F passes but the regulon truth fails evaluability, H1 is "not evaluable".
  - In mode "binding_only", every regulon claim says "DNA-binding regulons".
- **Human lines (only under verdict 2):** S2/S3 give "replicated in HepG2 / Jurkat" or "not replicated"; S4/S5 give regulon no-harm in human. Under any other verdict, human results are descriptive intervals only.
- **What does not count:**
  - a win at any depth other than the primary ones, including the same-depth contrast;
  - a win on all 1,484 mutants, on the growth-corrected matrix, on the 3–10 size class, on per-complex AP, on all-evidence or "possibly circular" regulons, or on any sensitivity truth when the primary fails;
  - a win after swapping the complex source (except the fixed fallback), the regulon mode or evidence filter, the exclusion list, the size class or the merge rule;
  - a win by a gate variant;
  - K562, RPE1 or *E. coli* results as support;
  - pooling yeast with human data;
  - results within ±0.005;
  - any p-value for an untested gated contrast.

## Protocol and gates (in order)

1. **Review, then fix this file and commit `excluded_pmids.txt`.** Replace the DRAFT block and commit. All later truth-touching steps check that the file is tracked, identical to HEAD (`git diff --quiet HEAD -- docs/preregistration_reciprocity_yeast.md`) and free of `STATUS: DRAFT`. They write the hash of the last commit touching the file into every output.
2. **Download** (`yeast_prep.R download`): expression, wild-type and annotation data only, with SHA-256 in `data/yeast/SHA256SUMS`. Commit the hashes.
3. **Matrices** (`matrix`, `check`, then `rebuild` only if the check fails, then `wt`): the processed mutant and WT matrices with their SHA-256, the dropped-column and masking logs, and the implementation-check table. Commit.
4. **Splits, clusters and gate F** (`splits`, `feasibility`): split assignments, mutant clusters and the F statistics with the pass/fail line, in `docs/reciprocity_yeast_selections/`. Commit. **Stop if F fails** (verdict 1).
5. **The synthetic implementation checks of R** pass (test script `--check_R`). Commit `R_check.csv`. If they fail: no verdict; fix and log.
6. **Label-free selections:** t\*_K10 and t\*_R per dataset (yeast, HepG2, Jurkat, K562), the label-free diagnostics, the SE_B proxy and the jackknife seeds, in `docs/reciprocity_yeast_selections/*_selections.csv`. Commit. The truth phase refuses to run unless the selections are byte-identical to the committed copies.
7. **The complete test script**, smoke-tested on synthetic and decoy truths (randomly drawn gene sets of the same sizes, never a real truth), is committed. Its git blob hash is recorded. Any change after step 8 is a logged deviation.
8. **Truth build, audit and power gate** (`yeast_prep.R truth_structure`). It refuses to run before steps 1 and 4–7 are committed and records the test-script hash. It:
   - performs and logs the YEASTRACT export, then commits the log;
   - downloads Complex Portal 559292 and the UniProt mapping;
   - builds the member lists with the exclusion and evidence filters;
   - computes evaluability counts, n_eff for P1, S1 and P2, and the contamination audit;
   - applies the CYC2008 fallback rule;
   - records SHA-256 of raw and built files;
   - reads membership only and never loads a network.

   Commit. **Stop if a primary truth fails evaluability or the power gate fails** (verdict 1).
9. **Code review** of the committed test script against this file, by a reviewer who has seen no yeast or HepG2 score. A required change is a logged deviation and is committed before step 10.
10. **Truth phase, once:** yeast (P1, P2, H1, S1, S6; gene and experiment jackknives), then HepG2 and Jurkat, then the descriptive sets. The first complete run is the confirmatory one. A rerun is allowed only after a crash, or for a bug shown to violate a named line of this file; both runs and the reason are reported.

## Compute (for `jobs_pc/sched`)

| Block | Core-hours | Memory |
|---|---|---|
| Yeast: matrix and WT build | under 1 | under 4 GB |
| Yeast: gate F (40 half-fits at N ≈ 188) | 2–4 | under 4 GB |
| Yeast: t\* (20 half-fits plus walks for two operators) | 15–25 | under 4 GB |
| Yeast: SE_B proxy (50 half-fits) | 8–12 | under 4 GB |
| Yeast: gene jackknife | 2–4 | under 4 GB |
| Yeast: experiment jackknife (200 refits at N ≈ 350, plus walks) | 40–60 | 4 GB per worker |
| HepG2 + Jurkat: pseudobulk build from the single-cell h5ad files (HepG2 5.6 GB; Jurkat size not retrieved) | 2–4 | about 60 GB, one job |
| HepG2 + Jurkat: t\*, walks, jackknife | 70–90 | 10–12 GB per worker |
| K562 descriptive | about 10 | 6 GB |

- **Total:** about 150–210 core-hours (supersedes the memo's 180–240).
- **Disk:** about 30 GB.
- All jobs run with `# mem` and `# cpus` headers, at about 70–80% utilization.

## Changes from rev 1 (review fd9c827)

| Finding | Change |
|---|---|
| B1 | Exclusion list E (GPL11232-linked PMIDs, Holstege/Kemmeren series PMIDs, reanalysis screen); primary regulon truth = DNA binding plus non-E expression evidence; fallback DNA-binding-only, decided on availability; complex evidence filter; contamination audit before scoring |
| B2 | Verdict 5 split into 5a (unresolved) and 5b (significant loss); structure-only n_eff(P2); explicit "low power is not a cost" |
| B3 | Experiment bootstrap replaced by the delete-half jackknife over label-free mutant clusters (r ≥ 0.5 complete linkage); √(SE_J² + SE_B²) kept |
| S1 | Synthetic rescore with the P1 and S1 constructions at t = 9 and own t\*; planning effect from it; CLR > K10 synthetic pattern and the round-1 clean ×0.95 regulon cost disclosed |
| S2 | F reworded as a sanity check; chance-corrected R′; overlap share Q ≤ 0.5 gates; failed-check fallback implemented; failed F does not consume yeast |
| S3 | S1 first at full α, then Holm over S2–S6 only if S1 is rejected; verdict 3b for S1 not evaluable |
| S4 | Power sensitivity to ρ = 1 and 2; realised SEs reported; SE_B proxy |
| S5 | Divergent-pair exclusion as a fixed descriptive sensitivity |
| S6 | Unresolvable columns dropped (SGD aliases used); overlapping and within-500-bp ORFs masked; paralog sensitivity; slow-growth rationale; growth-corrected matrix as a descriptive sensitivity |
| S7 | Gates for selections, the R check and the test script before the truth build; evaluability counts in `truth_structure`; −1 entries kept; export after fixing, logged; no R score on any real truth before the yeast phase |
| S8 | Support = cap-50 candidate set, with a unit test of B; self-loop order stated; jackknife gene set defined; H1 reuses P2's halves; S6 regulator defined; `make_regulons` mapping; S2/S3 need d ≥ 0.005; human SE inflated; human t\* on half A; ribosome rule via SGD names |
| S9 | No verdict for a failed R check; Holm over {P1, P2-NI, P2-harm} for individual claims; Collins timing relaxed by the PI (1 Nov) |
| N1–N10 | F seed 20261199; selection uncertainty and half-N t\* stated; S6 identifiability; same-chromosome sensitivity; H1 depth labelled; H not combined with P/S; relative changes reported; compute and tolerance reconciled; implementation-check fallback |

## Change log

- 2026-11-01: DRAFT rev 1 (18babdb).
- 2026-11-01: DRAFT rev 2: review fd9c827 addressed (see "Changes from rev 1"). Not fixed.
- 2026-11-01: DRAFT rev 3: PI decisions recorded. Primary depth rule kept after review of the synthetic rescore (option A); same-depth contrast (both at t\*_K10) added as a fixed descriptive secondary; S3 stop branch, ρ ≈ 0.44 disclosure and gate F thresholds kept as written; E. coli/EcoMAC counted as one compendium, so H1 on yeast is the first organism-independent replication. No endpoint, threshold or verdict rule changed. Not fixed.
