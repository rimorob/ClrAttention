# Preregistration: reciprocal attention and headline replication on a held-out yeast deletion compendium

> *Fixed 2026-11-01*

*Drafted by Claude (Opus, high effort, statistics lead).*

- *Rev 1 (18babdb) came from the design memo `docs/plans/reciprocity_heldout_design.md` and the PI decisions of 1 Nov 2026.*
- *Rev 2 answers the independent review `docs/plans/yeast_prereg_review.md` (fd9c827): three blockers (B1–B3), should-fix items S1–S9 and notes N1–N10. Where the review offered options, the PI chose them on 1 Nov; see "Changes from rev 1" at the end.*
- *Rev 3 records the PI decisions of 1 Nov 2026 on the four flags raised after rev 2 (see the disclosure rows "Primary depth rule" and "Flags resolved").*
- *Rev 5 answers the third review `docs/plans/yeast_prereg_review3.md` (858c08d): Holm on p-values only, with verdict precedence; a self-test that works on a fixed file; consistent SF3 wording; stricter commit-order code; the converter label fix. SF3 remains **provisional**.*
- *Rev 4 answers the second independent review `docs/plans/yeast_prereg_review2.md` (226e53e): blockers BL1–BL2, should-fix items SF1–SF8 and notes N1–N11. One item, SF3, is **provisional pending PI confirmation** (see "Decisions" below).*
- *No yeast or HepG2/Jurkat truth has been downloaded, opened or scored.*
- *No reciprocity score has been computed on any real dataset with any truth.*
- *After this file is fixed, any change is logged in the change log with its date and reason. A change made after the first truth-phase run is a deviation and is reported as one.*

## Question

1. **Reciprocity.** Does reciprocal attention improve recovery of small protein complexes over attention (K10) on data and truths never used in this project, without a regulon cost? Reciprocal attention keeps attention from i to j only when each gene is among the other's softmax candidates, and leaves unreturned mass on i.
2. **Headline replication (H1).** Does attention beat CLR on regulon co-membership in a second organism?

## Decisions

| Decision | Status |
|---|---|
| Primary depth rule kept after the synthetic rescore (option A); same-depth contrast descriptive | PI, 1 Nov 2026 |
| S3 stop branch; ρ ≈ 0.44 disclosed; gate F thresholds as judgment calls | PI, 1 Nov 2026 |
| Collins timing condition relaxed | PI, 1 Nov 2026 |
| "Consumed" counted per truth (not per compendium), with the disclosure condition for the unscored complex truth | PI-confirmed, 1 Nov 2026 |
| **SF3: the complex power gate is computed from complex sizes before any regulon membership is built. If it fails, the reciprocity arms are "not run" (verdict 1): no R score is computed anywhere, and the human and descriptive arms do not run. H1 alone is still tested, which consumes the yeast regulon truth. Consumption is counted per truth, so the unscored complex truth is not consumed** | **PI-confirmed, 1 Nov 2026** (drafted provisionally in rev 4–5) |

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
| **Yeast is single-use** | The deleteome expression and every yeast truth below are untouched in this project. The a3r follow-up (interpretation memo section 7.3) and the C3 asymmetry edge test (`docs/plans/asymmetry_expert.md`) are **not** included and therefore can never use these yeast truths (PI decision). **"Consumed" is counted per truth:** a yeast truth is consumed once any operator has been scored against it.
- A full truth phase consumes the regulon truth, the complex truth and the TF-edge truth.
- A failure of gate F consumes nothing, because no truth is downloaded before it.
- A failure of the complex power gate (or of P1 evaluability) builds complex membership and its size structure but scores nothing against it. Under the SF3 rule, H1 is then still tested, which consumes the yeast regulon truth (and the TF-edge truth is not scored).
- The unscored complex truth stays **unconsumed**. It may be used later only under a new preregistration that discloses three things: its size structure (n_eff) was seen here; the deleteome expression was already used for H1; and this file's reciprocity arms were "not run". |
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
  2. **Deleted gene.** Parse it from each column label (`<gene>-del-…`) and resolve it to an ORF via `SGD_features.tab`: standard name, then systematic name, then unique aliases (column 6; ambiguous aliases are logged and not used).
     - **Columns whose deleted gene cannot be resolved to an ORF in the matrix are dropped**, and their labels are logged.
     - **Verification (review2 SF6).** Before masking, the published M of the resolved ORF in its own column must be ≤ −1 and among that column's lowest 1% of values. **Columns that fail are dropped**, because a reassigned name may point to the wrong ORF. Each column's own-ORF M value is logged either way.
  3. **Masking.** In each mutant's column, mask the deleted ORF and every ORF whose SGD coordinates overlap it or lie within **500 bp** of it on the same chromosome (cassette and neighbouring-gene effects). Replace masked entries with that gene's median over the other mutants.
  4. **Missing values.** Drop ORFs with missing values in more than 5% of profiles (after masking, masked entries do not count as missing); impute the rest with the gene median.
  5. Drop constant rows.
  6. **No batch, slow-growth or SVD correction.** Rationale: the headline (M3D) and K562 compendia were analysed without such correction, and any choice here would be a new fork. The authors' growth-corrected matrix (`deleteome_all_mutants_svd_transformed.txt`) is a descriptive sensitivity only.
- **Gene universe:** the rows that remain (expected about 6,000).
- **Mutant clusters (label-free; used by the experiment jackknife; committed with the selections; review B3).**
  - Complete-linkage clustering of the mutant profiles (columns) on 1 − Pearson r over all universe genes, cut at height 0.5, so every pair within a cluster has r ≥ 0.5.
  - The number of clusters K, the cluster-size distribution and the range of replicate N are committed as a label-free diagnostic.
  - **Fallback, fixed now (review2 SF5).** If K < 20, or the largest cluster holds more than 25% of the mutants, the experiment jackknife instead keeps, in each cluster of size m, a random ⌊m/2⌋ or ⌈m/2⌉ members (each with probability ½): a delete-half over mutants stratified by cluster.
  - **Caveat, disclosed.** Complete linkage is the least inclusive linkage, so mutants sharing the slow-growth program can fall into many separate clusters, and SE_B may understate the variance that program contributes.
- **Reproducibility of random draws.** `RNGkind("Mersenne-Twister", "Inversion", "Rejection")`; R version and RNG kind are recorded in every output; collation is pinned with `Sys.setlocale("LC_COLLATE", "C")`, so every "alphabetical" order below is C-locale order; draws use `sample.int()` over objects sorted by ID.

### Gate F (label-free sanity check with the wild-type controls; computed and committed before any truth is downloaded; review S2)

1. **Matrices.**
   - WT: one profile per GSE45115 culture. Each profile is the mean of its dye-swap pair, oriented to log2(wt/ref), over duplicate spots, with the same probe mapping and gene filter.
   - The orientation rule is chosen by the **implementation check** among three candidates (review2 SF4):
     - "as deposited";
     - "negate the hybridization whose channel 1 is the test strain";
     - "negate the hybridization whose channel 1 is the reference".

     Applied to the sample tables of the first 10 GSE42527 mutants (C-locale order of deleted gene), a rule **passes** if it reproduces the published M with Pearson r ≥ 0.99 and median |difference| ≤ 0.02.
   - **Floor.** If no rule passes, the best rule (highest median r) is used only if its median r ≥ 0.8. Otherwise gate F is not computable, which gives verdict 1.
   - **Structural dye check (mutants and WT).** In every group, each hybridization must have exactly one channel whose source matches the reference pattern (`ref`, case-insensitive), and each dye-swap pair must contain one hybridization of each orientation. Any violation stops the stage and is logged. WT orientation is thereby verified with the same source fields and the same rule as the mutants.
2. **Failed implementation check (fallback, implemented).** If no rule passes (the published M is a limma estimate with gene-specific dye correction, so a simple mean can fail for benign reasons):
   - download the sample tables of a fixed random 376 mutants (seed 20261105);
   - rebuild those mutants with the best rule (highest median r, subject to the floor above);
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
   - Also reported, not gating: R′_mut, R′_wt and their spread; the number of genes per mutant beyond the gene's WT 0.5th–99.5th percentile range (computed in stage `feasibility`).
5. **If F fails:** verdict 1. No truth is opened, so yeast is not consumed.

### Truths (downloaded only after this file is fixed; structure counted, audited and committed before any scoring)

**Exclusion list E (fixed now; review B1).** E is the PubMed IDs of Holstege/Kemmeren-lab expression studies and the deleteome's constituent datasets.

- **(a) Every PubMed ID linked by GEO to a series on platform GPL11232** (query of 1 Nov 2026; 34 series, 20 PMIDs): 21145464 (van Wageningen 2010, kinases/phosphatases), 21596317 (Lenstra 2011, chromatin), 21976730, 21998594, 22697265, 23599000, 23785440, 23932120, 24009531, 24037263, 24324601, **24766815 (Kemmeren 2014)**, 24952590 (O'Duibhir 2014), 25157590, 26700642 (Sameith 2015), 27979920, 29979971, 31506296, 34520050, 36066422.
- **(b) PubMed IDs of *S. cerevisiae* GEO series with Holstege F or Kemmeren P as author, on other platforms:** 19713935, 21971086, 28049706, 28637236, 30203048.
- **(c) Reanalysis screen** (mechanical, applied by script to every PubMed ID cited as expression evidence): any reference whose PubMed XML record contains, anywhere (title, abstract, data-bank links or reference list), "deleteome", "GSE42526", "GSE42527", "GSE42528" or "GSE45115" is added to E.
- **Recall of E is unknown (review2 N3; stated in Limitations).** Lists (a) and (b) catch studies deposited on GPL11232 or in Holstege/Kemmeren GEO series, and the screen catches records that name the deleteome or its accessions. A study that reanalysed the deleteome without a deposit and without naming it in its PubMed record would stay in the primary truth. The audit measures what E removes, not what remains. The query strings for (a) and (b) are recorded in `excluded_pmids.txt`.

The list is committed as `docs/reciprocity_yeast_selections/excluded_pmids.txt` with this file. E restricts only the **expression** leg. DNA-binding evidence from any source counts.

**Yeast regulons (P2, H1, S6).**

- **Source.** YEASTRACT+, documented regulations, all environmental conditions, all TFs by all genes.
- **Timing of the export (review2 SF2, review3).** The export is protocol step 8b. It is performed only after steps 1–7 and 8a are committed **and pushed** to the upstream branch.
  - Its timestamp, release, form fields and parameters are recorded in `docs/reciprocity_yeast_selections/yeastract_export_log.md`, together with one line `raw_export: <sha256> <kind> <path>` per raw export file.
  - `regulon_structure` (step 8c) enforces the order in code. The commit that **adds** the export log (the earliest commit that adds it, `git log --diff-filter=A`) must be **strictly** later than (not equal to) the last commits of gate F, the selections, `R_check.csv`, the test script, `analysis/yeast_prep.R` and the step-8a counts. The log must not change after the commit that added it. Every one of these commits must be reachable from the upstream branch (run `git fetch origin` first). The stage refuses to run in a shallow clone.
  - The stage records the commit hashes it checked, and the check time, in its output. The remote's push records give the external timestamps.
  - **What code cannot prove:** git ancestry proves the order of commits, not the time at which the YEASTRACT web export was performed. That rests on the export log's timestamp and on the push order.
- **Form fields, fixed now (review2 SF1 iv):**
  - documented regulations;
  - evidence: DNA-binding evidence **or** expression evidence (both requested);
  - regulation type: activator, inhibitor and undefined (all);
  - environmental conditions: all groups;
  - every other field (including any high-/low-throughput or method filters) left at its default.

  The value of every form field is recorded in the log.
- **Converter (review2 SF1 iii).** `analysis/yeastract_convert.R` turns the export into the normalized table (tf, target, evidence ∈ {binding, expression}, pmids attributed to that evidence row, pmids_unattributed, value).
  - Evidence labels are mapped by a fixed table in the script.
  - Labels are split into tokens at `; , / |`, and each token must match exactly one **anchored** entry, as a whole token (for example "Direct" → binding, "Indirect" → expression, "ChIP-chip" → binding). A row's classes are the union of its tokens' classes, so "Expression; ChIP-chip" counts as both, with its references treated as a combined, unattributed list. **Any unrecognised token stops the converter**; nothing is silently dropped.
  - It is committed **before** the export (step 7). Step 8c **runs the committed converter itself** on the raw export files listed, with SHA-256, in the export log. A normalized table supplied by hand is never accepted. Its hash goes into `yeast_regulons.sha256`.
  - **Converter deviation route.** If the committed converter cannot parse the real export, the fix is a logged deviation. The fixed converter and a `converter_deviation.md` (what changed and why) are committed **after** the export log and before any association row is used. The 8c gate accepts a converter changed after the export only through this route. The converter blob recorded in 8a stays the reference.
- **Mode (review2 SF1 i–ii), computed by `regulon_structure` from the normalized table, never chosen by hand.**
  - **"refs"** applies if at least **95%** of the associations that carry expression evidence have at least one PMID that YEASTRACT attributes to that expression evidence.
  - Otherwise, including when PMIDs come only as one combined list per association (`pmids_unattributed`), the mode is **"binding_only"**.
  - The command-line flag is a cross-check that must agree.
  - **Mode "refs" (primary):**
    - **Primary truth:** TF–target associations with DNA-binding evidence (any reference) **and** expression evidence from **at least one PubMed ID not in E**.
    - "DNA binding and any expression evidence" is a descriptive comparator, labelled "possibly circular".
  - **Mode "binding_only" (fallback):** used if the 95% rule is not met. The primary truth is DNA-binding-evidence associations only (ChIP, ChIP-chip, ChIP-seq, ChIP-exo, EMSA, footprinting, as classified by YEASTRACT), and the paper says so. "DNA binding and expression" is then descriptive only, labelled "possibly circular".
- **Encoding.** Any non-zero matrix entry, including −1 for negative regulation, counts as an association. The value table of the export is logged.
- **Co-membership.** Regulator-agnostic, as `make_regulons()` / `comembership_pairs()` adapted to yeast: regulator = TF ORF, target = gene ORF, every association counted as one evidence class. Regulons of 5–500 in-universe targets.
- **TF–target edges (S6).**
  - Universe: pairs with exactly one regulator gene, where a regulator gene is a TF with a 5–500-target regulon in the primary truth.
  - Positives: primary associations whose target is not itself a regulator gene. TF–TF pairs are outside the universe (review2 SF8d).
  - S6 is non-inferiority only. With the deleted gene masked, the deleteome's causal evidence for X → Y (column X) is invisible to the row-based edge (X, Y), so the edge signal comes only from co-variation of X's mRNA across other columns. Many yeast TFs are regulated post-translationally, and S is symmetric, so direction is not identifiable.
- **Contamination audit (pre-scoring, truth structure only; committed with the counts).** In mode "refs", report:
  - (a) the share of "DNA binding and any expression" associations whose expression evidence comes **only** from E (removed from the primary truth);
  - (b) the share of all expression-evidence associations whose only references are in E;
  - (c) the same shares in co-membership positive pairs.

  In mode "binding_only", report only the share of binding associations that also carry expression evidence; the removal shares are NA, because without references they cannot be split by source (review2 N4).

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
- **Fallback, fixed now.** If n_eff(P1) < 100 with Complex Portal, use the union of Complex Portal and CYC2008 (Pu et al. 2009; ORF lists), merged by the same rule, decided from counts alone before scoring.
  - **CYC2008 procedure (two runs of 8a).** The first 8a run is committed with the fallback status "CYC2008 required but not available". Then the CYC2008 file is downloaded, and its log `cyc2008_download_log.md` (source URL, date, SHA-256) is committed **after** that first run under the same strict ancestry rule. Then 8a is re-run with `--cyc2008`. The gate checks the file's SHA-256 against the log. Both runs' counts stay in the git history.
  - If CYC2008 cannot be obtained, or the union still has n_eff < 100, the reciprocity arms get verdict 1 (H1 per the SF3 rule).
  - The unused source is never scored.
- **Ribosome (descriptive).** Merged complexes with ≥ 30% cytoplasmic r-protein ORFs: the ORFs of `SGD_features.tab` standard names matching `^RP[LS][0-9]+[AB]?$`.

**Human (S2–S5).** The existing CORUM file (`data/human_truth/complexes.tsv`, as `corum_sets()`, then the overlap merge, then 3–6-gene merged complexes as P1) and CollecTRI (`data/human_truth/collectri.tsv`, as in K562).

**Evaluability (computed in `complex_structure` and `regulon_structure`, before scoring).**

| Truth | Requirement |
|---|---|
| P1, S1 | ≥ 200 positive pairs and ≥ 30 merged complexes with positives |
| P2 / H1 | ≥ 200 positive pairs |
| S6 | ≥ 200 positive edges |
| Human truths | ≥ 200 positive pairs |

- A failing P1 truth gives verdict 1 for the reciprocity arms (H1 per the SF3 rule).
- A failing regulon truth gives verdict 1 for the reciprocity arms and makes H1 "not evaluable".
- A secondary truth that fails is "not evaluable", counts in its family, and is reported.
- Verdict 3b (S1 not evaluable) is practically reachable only after a crash, because S1 contains P1's complexes (review2 N6).

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
| P2 | yeast regulons (primary truth, mode as computed) | non-inferiority, margin 0.005 |

**Headline replication, a separate family H at one-sided α = 0.025.**
- H1 is tested whatever P1 and P2 show. Under the SF3 rule it is also tested when the reciprocity arms are not run because the complex power gate or P1 evaluability failed, provided gate F passed and the regulon truth is evaluable. In that case H1 is the only yeast score computed.
- The paper never combines family H with the P/S families (N7).
- H1 superiority and the "significant CLR advantage" statement have nulls that are both true at d = 0, so family H is in effect a two-sided 0.05 test with directional conclusions (review2 N7).

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
  - **Seeds per truth (review3):** P1 20261102, P2 and H1 20261112, S1 20261122, S6 20261132. The halves of a truth therefore do not depend on whether any other truth is evaluated or scored. Human truths: 20261102, set once per human dataset run, drawn in table order. Never rerun with another seed.
  - SE_J = the SD of the 100 paired differences.
- **Experiment jackknife (SE_B, yeast only; review B3).**
  - The delete-half jackknife over the label-free **mutant clusters**: 200 replicates. Replicate b keeps `sample.int(K, floor(K/2))` of the clusters sorted by ID, with seed 20262000 + b (b = 1, …, 200; outside every other seed used here, review2 N1), and all their mutants. Under the fallback (K < 20 or a cluster above 25%), the stratified rule of "Data" applies, with the same seeds.
  - MI, CLR, K10 and R are refitted in every replicate at the committed t\* values; truths and gene universe are fixed.
  - With d = n/2 clusters deleted, the delete-d variance factor (n − d)/d is 1, so SE_B = the SD of the 200 paired differences.
  - No duplicate columns arise, and the clustering of near-identical deletion profiles is respected.
  - It is conservative, because half-N networks are noisier than the full-N network.
- **Deciding statistic (yeast P1, P2, H1, S1, S6).** SE = √(SE_J² + SE_B²).
  - superiority: p = 1 − Φ(d/SE);
  - non-inferiority: p = 1 − Φ((d + 0.005)/SE);
  - harm (a significant loss): p_harm = 1 − Φ(−d/SE), with null d ≥ 0, and a loss is stated only if also d ≤ −0.005.
  - In a crossed genes × experiments model, the sum counts the interaction component twice (Owen 2007), so it is conservative and does not need independence.
- **Human tests (S2–S5).** SE = SE_J × √(1 + 0.44²), a pre-fixed inflation for the experiment component that is not estimated there (the ratio observed in M3D for CLR(log Pᵗ)). They are labelled accordingly.
- **Non-inferiority and harm have disjoint nulls** (d ≤ −0.005 and d ≥ 0), so at most one of them is true. For the same endpoint, each can therefore be tested at the full level assigned to that endpoint without inflating the familywise error (review2 BL2).
- **Claims outside the IUT (verdicts 4 and 5; review2 BL2, review3).** These use **Holm over two members** at one-sided 0.025:
  - H_P1: P1 superiority;
  - H_P2: the P2 direction, either P2-NI or P2-harm. These are disjoint nulls, each tested at the level Holm assigns to H_P2, and the member's p-value for ordering is min(p_NI, p_harm).
- **The step-down uses p-values only, with "p < level" throughout.** The effect-size thresholds never enter the step-down. They decide only the wording of a claim:
  - d ≥ +0.005 for "small-complex gain";
  - d ≤ −0.005 for "loss".
- **Loss definition.** A "significant regulon loss" means P2-harm is rejected under this Holm rule **and** d₂ ≤ −0.005. If P2-harm is rejected with −0.005 < d₂ < 0, the result is reported as "a regulon decrease smaller than the 0.005 margin". That is not a loss claim and not a verdict change.
- **FWER of the combined rule (disclosed).** The procedure is the IUT if both components pass, else Holm over {H_P1, H_P2}. It controls the familywise error at 0.025 under independence or positive dependence of the two statistics. Under adversarial negative dependence the worst case is 1.5 × 0.025.
- **S1 harm.** Inside the fixed sequence, S1 occupies one slot at 0.025. S1-harm is tested at 0.025 only when S1-NI is not rejected (disjoint nulls). "Large-complex cost" means S1-harm rejected with d ≤ −0.005; S1-harm rejected with −0.005 < d < 0 is "a decrease smaller than the margin".
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
| R vs K10, own t\* (**the primary contrast's form**) | ×1.09–2.41 (median ×1.27 over the 14 cells; ×1.36 over the 12 standard cells), d +0.050 to +0.171; 4/4 in all 14 cells | ×1.00–1.09 | **×0.86–1.16**; at batch ×2, **×0.86–0.96 (0/4 seeds) in all four cells** (d −0.001 to −0.005, at regulon AUPR levels of 0.02–0.04) |
| R (own t\*) vs CLR | ×0.88–1.21 | ×1.01–1.46 | ×0.79–1.69 |
| K10 (own t\*) vs CLR | **×0.39–0.96**, 0/4 to 2/4 seeds | ×1.01–1.39 | ×0.91–1.45 |
| K10 (t = 9) vs CLR, f = 10% | ×0.84–1.10 | ×0.98–1.29 | ×1.00–1.71 |

**What this means, disclosed.**

1. **In the P1 construction (complex genes only) the round-2 picture "CLR ≫ K10 on small complexes" largely disappears at t = 9.** K10 and CLR are close (×0.84–1.10 at f = 10%), and R adds ×1.03–1.64.
2. **Criterion A selects deep depths in these synthetic data:** K10 t\* 10–24, R t\* 16–32. R's t\* is deeper by 0–20 in t units per seed (cell medians +2 to +16). At its own t\*, K10 loses heavily to CLR on small complexes (×0.39–0.96).
   - Much of R's own-t\* advantage over K10 (×1.09–2.41) is therefore a **depth effect**: R's unreturned mass stays home, so its walk degrades small complexes more slowly with t.
   - The common-depth contrast (both at t = 9) is the smaller ×1.03–1.64.
   - This is a known property of the primary contrast as written in rev 1. **The contrast was not changed after the rescore;** the PI reviewed it and kept it on 1 Nov 2026. The same-depth contrast (both at t\*_K10) is a fixed descriptive secondary, and both at t = 9 is descriptive.
3. **Regulon risk.** Under strong batch effects (×2) R's deeper t\* costs regulons ×0.86–0.96 (0/4 seeds in all four cells). On a real regulon AUPR of about 0.1, a 14% relative loss would exceed the 0.005 margin. This is the main synthetic risk for P2.
4. **R vs CLR on small complexes** is ×0.88–1.21. R is not reliably better than CLR on P1. This file does not test or claim it.

**Planning effect: +0.010 AUPR (P1, R vs K10 at own t\*).**

- Synthetic P1 AUPR levels (cell means) are CLR 0.19–0.69, K10 at t = 9 0.17–0.73, K10 at own t\* 0.075–0.62, R at t = 9 0.25–0.75, and R at own t\* 0.17–0.71. Real small-complex AUPRs are expected near 0.05–0.15 (K562 CORUM attention: 0.047, all complexes).
- The synthetic absolute differences (+0.05 to +0.17) are **not a prediction** for real data. But they are not uniformly larger only because of higher levels: in the batch ×2 cells, K10's own-t\* level (0.075–0.14) lies inside the expected real range, and there R gained +0.09 to +0.11 (×1.8–2.4).
- **+0.010 is a deliberately conservative planning value, not an expected effect.** It is the bottom of the synthetic own-t\* ratio range (×1.09, `regheavy f10 b1`; the standard-cell minimum is ×1.12) applied to a K10 level of about 0.10, or ×1.2 at 0.05.
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
  - **Label-free SE_B proxy** (selection stage, reported, not gating): the SD of the criterion-A AP difference R vs K10 over the first 50 cluster-jackknife replicates. In replicate b the networks are fitted on the kept clusters at the committed t\* values, and the target is the top-1% CLR pairs of the deleted clusters' mutants (review2 N9).
  - **Order (SF3).** n_eff(P1) is computed in step 8a from Complex Portal membership and complex sizes only (and CYC2008 under the fallback), before any regulon membership is built. No expression value or network is used.
- **Regulon arm (B2).** n_eff(P2) is computed the same way (Σ over regulons of within-regulon positive pairs) and reported by `regulon_structure`, not gating. *E. coli* mild K10 tweaks had paired regulon SEs of 0.0005–0.0022; at 0.0015 and ρ = 0.44, NI power at d = 0 is 0.86. **A failure to show non-inferiority under low power is not evidence of a regulon cost** (verdict 5a vs 5b).
- **If the power gate fails** (after the CYC2008 fallback): verdict 1 for the reciprocity arms. Under the SF3 rule (review3; PI-confirmed 1 Nov 2026):
  - **no R score is computed against any truth**, not even descriptively;
  - the HepG2, Jurkat, K562 and *E. coli* arms and every descriptive set do not run;
  - H1 (family H, K10 vs CLR on yeast regulons) is the only score computed, if gate F passed and the regulon truth is evaluable. It consumes the yeast regulon truth.

## Decision rules and verdicts

The verdict is computed by the script and written to `verdict.csv` before anyone interprets the numbers.

| # | Condition | Verdict (one line) | Reporting |
|---|---|---|---|
| — | Synthetic implementation check of R fails (protocol step 5) | **No verdict.** Fix the implementation and log the cause before any real data are fitted | — |
| 1 | Gate F fails or is not computable; or the power gate fails after the fallback; or a primary truth fails evaluability | **Not run:** no adequately powered held-out test of reciprocity was possible | Reciprocity stays a synthetic-only, explained candidate; the gate numbers are reported. No R score is computed; no human or descriptive arm runs. If F failed, nothing else runs and no truth is consumed. Otherwise H1 alone is tested under the SF3 rule if the regulon truth is evaluable, consuming the regulon truth only |
| 2 | P1 and P2 pass, and S1 is rejected | **Confirmed, no large-complex cost:** reciprocal attention becomes the method's operator from here on | The *E. coli* headline keeps K10, and the order of events is stated; S2–S6 per Holm |
| 3a | P1 and P2 pass, and S1 is evaluable but S1-NI is not rejected | **Confirmed as a complex-mode variant:** it improves small complexes; non-inferiority on all complexes was not shown | K10 stays the default; if S1-harm is rejected at 0.025 with d ≤ −0.005 (disjoint nulls), "large-complex cost" is said; S2–S6 descriptive |
| 3b | P1 and P2 pass, and S1 is not evaluable | **Confirmed as a complex-mode variant; the large-complex question could not be evaluated** | as 3a, without any large-complex statement |
| 4 | P1 fails at the IUT level (p ≥ 0.025 or d < +0.005), whatever P2 shows | **Not confirmed** | Wording as in CLR(log Pᵗ): "not confirmed", with estimate and interval, never "no effect". **Regulon line, always evaluated:** Holm over {H_P1, H_P2} on p-values; if P2-harm is rejected, add "significant regulon loss" (d₂ ≤ −0.005) or "a regulon decrease smaller than the margin" (−0.005 < d₂ < 0) |
| 5a | P1 passes at the IUT level, P2-NI is not shown, H_P1 is rejected under Holm over {H_P1, H_P2} (p-values only), and there is **no significant regulon loss (as defined at the loss-definition line in "Statistics")** | **Not adopted: small-complex gain with an unresolved regulon question** (non-inferiority not shown) | Not a trade-off claim. If P2-harm was rejected with −0.005 < d₂ < 0, add "a regulon decrease smaller than the margin". The realised P2 SE is reported against the projection |
| 5b | P1 passes at the IUT level, P2-NI is not shown, H_P1 is rejected under Holm, and there is a significant regulon loss (P2-harm rejected under Holm and d₂ ≤ −0.005) | **Regulon trade-off:** small complexes improve at a regulon cost | Not adopted |
| 5c | P1 passes at the IUT level, P2-NI is not shown, and H_P1 is not rejected under Holm (then P2-harm cannot be rejected either) | **Not confirmed as a separate small-complex claim; regulon non-inferiority not shown** | P1 estimate and interval only; no gain and no cost claimed |

- **Family H, reported alongside:**
  - H1 holds: "the headline replicates in yeast".
  - Otherwise: "attention did not beat CLR on yeast regulons", with the interval. A significant CLR advantage (d ≤ −0.005, harm p < 0.025) is said plainly.
  - H1 never changes verdicts 1–5. If gate F passes but the regulon truth fails evaluability, H1 is "not evaluable". If gate F fails, H1 is not run.
  - In mode "binding_only", every regulon claim says "DNA-binding regulons".
- **Precedence.** Rows are applied in this order, and the first that matches is the verdict: no verdict (failed R check), 1, 4, 2, 3a, 3b, 5c, 5b, 5a. All conditions use "p < level". The rows are mutually exclusive and exhaustive.
- Verdicts 5a, 5b and 5c partition the case "P1 passes at the IUT level and P2-NI is not shown". Given p₁ < 0.025, if H_P1 is not rejected under Holm, H_P2 cannot be either.
- **Worked example (review3).** p₁ = 0.020, d₁ = +0.008; P2-NI not shown; p_harm = 0.005, d₂ = −0.004.
  - Holm on p-values rejects H_P2 at 0.0125, then H_P1 at 0.025.
  - d₂ > −0.005, so there is no significant regulon loss.
  - The verdict is **5a**, with the added sentence "a regulon decrease smaller than the margin".
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

1. **Review, then fix this file and commit `excluded_pmids.txt`.** Replace the status block at the top with the single fixed-date line and commit.
   - All later truth-touching steps check that the file is tracked, identical to HEAD (`git diff --quiet HEAD -- docs/preregistration_reciprocity_yeast.md`), and carries the positive fixed-date marker described in the header.
   - They write the hash of the last commit touching the file into every output.
   - Before committing the fixed file, run `Rscript analysis/yeast_prep.R selftest_gate` on it; it must print "fixed" (it prints "draft" on this DRAFT).
   - `excluded_pmids.txt` must be added at or before the fixing commit; the gate checks this by ancestry.
   - **Every gate commit (steps 1, 4, 6, 7, 8a and the export log) is pushed to the upstream branch before the next step begins.** The truth stages check that each one is reachable from the upstream branch.
2. **Download** (`yeast_prep.R download`): expression, wild-type and annotation data only, with SHA-256 in `data/yeast/SHA256SUMS`. Commit the hashes.
3. **Matrices** (`matrix`, `check`, then `rebuild` only if the check fails, then `wt`): the processed mutant and WT matrices with their SHA-256, the dropped-column and masking logs, and the implementation-check table. Commit.
4. **Splits, clusters and gate F** (`splits`, `feasibility`): split assignments, mutant clusters and the F statistics with the pass/fail line, in `docs/reciprocity_yeast_selections/`. Commit. **Stop if F fails** (verdict 1).
5. **The synthetic implementation checks of R** pass (test script `--check_R`). Commit `R_check.csv`. If they fail: no verdict; fix and log.
6. **Label-free selections:** t\*_K10 and t\*_R per dataset, the label-free diagnostics, the SE_B proxy, the jackknife seeds, and the R version and RNG kind, in `docs/reciprocity_yeast_selections/{yeast,hepg2,jurkat,k562}_selections.csv`. Commit. The truth phase refuses to run unless every one of these files is committed and byte-identical to the working copy.
7. **The complete test script and the YEASTRACT converter**, smoke-tested on synthetic and decoy truths (randomly drawn gene sets of the same sizes, never a real truth), are committed. `analysis/yeast_prep.R` itself is gated too: it must be committed and clean, and its blob is recorded in both count files. Their git blob hashes are recorded. Any change after step 8a is a logged deviation; for the converter, see the deviation route.
8. **Truth structure, in three sub-steps; each reads membership only and never loads a network or expression values.**
   - **8a. Complexes and the power gate** (`yeast_prep.R complex_structure`; refuses to run before steps 1 and 4–7 are committed):
     - downloads Complex Portal 559292 and the UniProt mapping;
     - applies the evidence filter, the merge and the CYC2008 fallback (its download is logged and committed under the ancestry rule);
     - computes P1/S1 evaluability counts and n_eff(P1);
     - records SHA-256 of raw and built files and the test-script and converter hashes.

     Commit and push. If the power gate or P1 evaluability fails, the reciprocity arms get verdict 1 (SF3: only H1 continues, through 8b, 8c and step 10).
   - **8b. YEASTRACT export:** performed only after 8a is committed and pushed, with the form fields above. Write the export log, including the `raw_export:` hash lines, then commit and push it in a commit of its own. The code enforces strict ancestry.
   - **8c. Regulons** (`yeast_prep.R regulon_structure`): checks the raw export files against the log's SHA-256, runs the committed converter on them, computes the mode, applies E and the reanalysis screen, and builds regulons and TF edges. It also computes the contamination audit, the P2/S6 evaluability counts and n_eff(P2), and records SHA-256. It verifies the expression matrix against `deleteome_matrix.sha256` and gates `masking_log.csv` and `dropped_columns.csv`. Commit.
9. **Code review** of the committed test script against this file, by a reviewer who has seen no yeast or HepG2 score. A required change is a logged deviation and is committed before step 10.
10. **Truth phase, once:** yeast (P1, P2, H1, S1, S6; gene and experiment jackknives), then HepG2 and Jurkat, then the descriptive sets. **Under verdict 1 from step 8a, only H1 runs** (K10 and CLR on yeast regulons, both jackknives). No R score, no human arm and no descriptive set run. The first complete run is the confirmatory one. A rerun is allowed only after a crash, or for a bug shown to violate a named line of this file; both runs and the reason are reported.

## Limitations to state in any report

- The recall of the exclusion list E is unknown (see "Truths").
- One lab, one platform, one perturbation design, with the responsive-only selection.
- Power is adequate only under the ρ ≈ 0.44 assumption.
- R's own-t\* advantage may be partly a depth effect (rescore section).
- Complete-linkage clusters may split shared programs, so SE_B may understate their variance.

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

## Changes from rev 4 (review3 858c08d)

| Finding | Change |
|---|---|
| BL2 residual | Holm step-down on p-values only; d-thresholds qualify wording only; loss definition with a "decrease smaller than the margin" wording; 5a reads "no significant regulon loss (as defined)"; verdict 4 regulon line always evaluated; precedence order stated; worked example (p₁ 0.020, p_harm 0.005, d₂ −0.004) → 5a; 1.5α worst case disclosed |
| Self-test | `selftest_gate` reports "fixed"/"draft" for any file (`--prereg`), validates a real date and the absence of a DRAFT line, and tests derived copies; passes on this draft and on a fixed copy |
| SF3 consistency | Under verdict 1 from 8a: no R score anywhere, no human or descriptive arm, H1 only; consumption counted per truth; the unscored complex truth stays unconsumed (with a disclosure condition); per-truth jackknife seeds; still provisional in rev 5 (PI-confirmed when fixed, 1 Nov 2026) |
| SF2 code | Adding commit of each log, strictly later than the earlier commits; log must not change after it was added; shallow clones refused; every gate commit must be on the upstream branch; `yeast_prep.R` gated and its blob recorded; 8c runs the committed converter on the hashed raw export; converter deviation route; CYC2008 two-run procedure; ancestry check for `excluded_pmids.txt`; `selftest_git` tests these rules in a throw-away repository |
| Converter | Anchored whole-token label table ("Indirect" → expression; "Expression; ChIP-chip" → both, unattributed); stop on any unrecognised label; value column never an "evidence" column |
| Notes | Stage name in `excluded_pmids.txt` corrected; what code cannot prove about export timing stated |

## Changes from rev 3 (review2 226e53e)

| Finding | Change |
|---|---|
| BL1 | Gate uses a positive fixed-date marker (regex) instead of the absence of a string; step 1 reworded; `selftest_gate` stage tests the check on this file (expects "draft") and on a fixed copy (expects "fixed") |
| BL2 | NI and harm declared disjoint-null tests at full level; claims outside the IUT use Holm over {H_P1, H_P2}; 5a/5b conditions wired to it; new 5c; S1-harm defined at 0.025 when S1-NI is not rejected |
| SF1 | YEASTRACT mode computed by code (≥ 95% of expression-evidence associations with attributed PMIDs); combined-only PMIDs give binding_only; converter committed before the export with a fixed label table; form fields fixed |
| SF2 | Step 8 split into 8a (complexes, power gate), 8b (export), 8c (regulons); ancestry enforced with `git merge-base --is-ancestor`; CYC2008 under the same rule |
| SF3 | Provisional in rev 4 (PI-confirmed when fixed, 1 Nov 2026): complex power gate before any regulon membership; on failure the reciprocity arms are verdict 1 but H1 is still tested and consumes the regulon truth; "consumed" defined |
| SF4 | Third dye rule; median r floor 0.8 (else F not computable, verdict 1); structural dye-balance check for mutants and WT |
| SF5 | Fallback for K < 20 or a cluster > 25% (cluster-stratified delete-half); RNG and R version pinned; complete-linkage caveat disclosed |
| SF6 | Resolved deleted ORF must have M ≤ −1 and lie in its column's lowest 1%, else the column is dropped |
| SF7 | Median ×1.27 over 14 cells; 0/4–2/4 seeds; AUPR levels corrected; "not transferable" restated; t\* gaps in t units |
| SF8 | Screen text = whole PubMed XML record; all four selection files gated; WT percentile diagnostic computed; TF–TF pairs excluded from S6; matrix hash and logs gated; prep header rev 4 and universe read from row names only |
| N1–N11 | Experiment-jackknife seeds moved to 20262001–20262200; C-locale collation; recall of E disclosed and query strings recorded; binding_only audit shares NA; 3b reachability noted; family H effectively two-sided; ambiguous aliases logged; SE_B proxy target defined; dubious ORFs skipped in divergent-pair adjacency; hashes recorded after `rebuild` and `wt` |

## Change log

- 2026-11-01: DRAFT rev 1 (18babdb).
- 2026-11-01: DRAFT rev 2: review fd9c827 addressed (see "Changes from rev 1"). Not fixed.
- 2026-11-01: DRAFT rev 3: PI decisions recorded. Primary depth rule kept after review of the synthetic rescore (option A); same-depth contrast (both at t\*_K10) added as a fixed descriptive secondary; S3 stop branch, ρ ≈ 0.44 disclosure and gate F thresholds kept as written; E. coli/EcoMAC counted as one compendium, so H1 on yeast is the first organism-independent replication. No endpoint, threshold or verdict rule changed. Not fixed.
- 2026-11-01: **FIXED.** The PI confirmed SF3 (complex power gate from Complex Portal sizes before any regulon truth is built; on failure the reciprocity arms get verdict 1 but H1 is still scored, consuming the yeast regulon truth) and the per-truth definition of "consumed" with its disclosure condition. Provisional flags removed; status block replaced by the fixed-date marker. No other change from DRAFT rev 5 (ba73a50). `excluded_pmids.txt` was committed earlier (18babdb onward).
- 2026-11-01: DRAFT rev 5: review3 858c08d addressed (see "Changes from rev 4"); no endpoint, margin or threshold changed; verdict 4's regulon line and the Holm rule made explicit; SF3 still provisional. Not fixed.
- 2026-11-01: DRAFT rev 4: review2 226e53e addressed (see "Changes from rev 3"). Verdict 5c added and the claim rule wired into verdicts 4–5; the order of the truth-structure steps changed (8a–8c); SF3 provisional pending PI confirmation. Not fixed.
