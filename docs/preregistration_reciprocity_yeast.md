# Preregistration: reciprocal attention and headline replication on a held-out yeast deletion compendium

> **STATUS: DRAFT, under review (1 Nov 2026). Not fixed.** This file becomes a preregistration only when this block is replaced by "*Fixed YYYY-MM-DD*" in a commit. Every gate in `analysis/yeast_prep.R` (truth stage) and in the future test script refuses to run while this file contains the string `STATUS: DRAFT`.

*Drafted by Claude (Opus, high effort, statistics lead) from the design memo `docs/plans/reciprocity_heldout_design.md` and the PI decisions of 1 Nov 2026:*

- *GO on yeast;*
- *the same file carries H1;*
- *reciprocity alone, with no gate;*
- *a structure-only power gate n_eff ≥ 100, with a CYC2008 fallback fixed in advance;*
- *this is the last complex-oriented operator test;*
- *HepG2 and Jurkat are secondary;*
- *a deletion compendium with calibrating wild-type controls;*
- *the combined jackknife + experiment-bootstrap SE is the deciding statistic.*

*No yeast or HepG2/Jurkat truth has been downloaded, opened or scored. No reciprocity score exists on any real dataset. Any change after this file is fixed is logged in the change log with its date and reason; a change made after the first truth-phase run is a deviation and is reported as one.*

## Question

1. **Reciprocity.** Does reciprocal attention improve recovery of small protein complexes over attention (K10) on data and truths never used in this project, without a regulon cost? Reciprocal attention keeps attention from i to j only when j also attends to i, and leaves unreturned mass on i.
2. **Headline replication (H1).** Does attention beat CLR on regulon co-membership in a second organism?

## Disclosure: what has been touched, and what was decided after seeing data

| Item | Status |
|---|---|
| **Reciprocity evidence** | Synthetic only (`analysis/synthetic_width/`, design memo section 8). Round 2, at partner-less fractions of 5–50% and batch strengths ×0, ×1, ×2: small complexes ×1.86–6.24 vs K10 (4/4 seeds, all 12 standard cells); regulons ×1.05–1.09; regulon-heavy generator ×1.02–1.04; TF–target edges ×1.01–1.02; ribosome ×0.80–0.83; all complexes 0.94–1.04. The variant was added after one exploratory seed. |
| **Rule change (all complexes)** | Reciprocity **failed** the pre-set carry-forward rule on condition C2 (all complexes > 1.00 with ≥ 3/4 seeds), because of the large complexes. Testing it means treating C2 as non-binding, given the PI's choice of small complexes as the target. That decision was made **after** the round-2 tables were seen. Here all complexes are a key secondary non-inferiority test (S1), not part of adoption. |
| **Gate variants** | The reproducibility gate, gate plus reciprocity and the empirical-null gate failed round 2. They are not tested here. The round-1 gate gains came from a generator with 51% partner-less genes. |
| **Serial testing** | This is the **third** complex-oriented operator test, after keff (`docs/keff_verdict/`: headline unchanged) and CLR(log Pᵗ) (`docs/clr_logpt_verdict/`: not adopted). **It is declared the last.** Whatever the verdict, no further budget- or readout-level operator change aimed at complexes is tested confirmatorily in this project. |
| **K562 in-sample** | About 20 operators and readouts were scored on K562 CORUM, including `d_geomean_t9` (√(PᵢⱼPⱼᵢ), a mutual-attention *readout*, not this operator; `analysis/corum_exploratory/metrics_exploratory.csv`). The drafter saw that row. K562 shaped the PI's focus on complexes. K562 is descriptive here. |
| **E. coli** | M3D and EcoMAC truths are consumed (headline, keff, CLR(log Pᵗ)). They are descriptive only (regulon no-harm). The M3D–EcoMAC array overlap does not matter here, because neither is confirmatory. |
| **Yeast is single-use** | The deleteome expression and every yeast truth below are untouched in this project. Once the truth phase runs, yeast is consumed for good. The a3r follow-up (interpretation memo section 7.3) and the C3 asymmetry edge test (`docs/plans/asymmetry_expert.md`) are **not** included and therefore can never use these yeast truths (PI decision). |
| **What the drafter has seen of the yeast data** | The column headers and the first 1–3 rows (probes SCAB000001–3, M and p-values of about 10 mutants) of `GSE42527_final_GeneExpressionMatrix.txt.gz` and `GSE42526_final_GeneExpressionMatrix.txt.gz`, and the first 35 rows of one wild-type sample table (GSM630359). That was to verify the format. No truth was seen. |
| **Human (secondary)** | HepG2 and Jurkat expression (GSE264667) has never been downloaded. CORUM and CollecTRI were used on K562 (in-sample) and RPE1 (CORUM descriptive in CLR(log Pᵗ), S5). |
| **Fixed without labels on M3D** | keff 10, cap 50, α = 0.5, B-spline MI (`bins = "hg"`), Stouffer CLR. Unchanged. |

## Data

### Expression (primary, confirmatory): Kemmeren et al. 2014 deleteome (*Cell* 157:740; PMID 24766815)

- **Design** (GEO overall design, verified 1 Nov 2026):
  - deletion strains of BY4742 in synthetic complete medium, mid-log;
  - two independent cultures per mutant, on two two-channel arrays against one common wild-type reference RNA, with a dye swap;
  - each gene spotted twice, so four measurements per gene per mutant;
  - wild-type cultures grown in parallel with the mutants every day.
- **Processed mutant data**, in the GEO supplementary files of SuperSeries GSE42528:
  - `GSE42527_final_GeneExpressionMatrix.txt.gz` (63.8 MB): **700 "responsive" mutants**, defined as ≥ 4 significant changes;
  - `GSE42526_final_GeneExpressionMatrix.txt.gz` (67.5 MB): 784 "non-responsive" mutants, 6,182 probes × (M, p-value) pairs, 1,569 columns;
  - per mutant: M (limma average log2 mutant/reference) and a p-value (limma 2.12.0, Benjamini–Hochberg-adjusted, per `GSE42527_protocols.txt.gz`);
  - platform GPL11232 (A-UMCU-Y16k-1.3, Holstege lab); probe IDs `SCAB######`.
- **Calibrating wild-type controls (verified):**

  | Series | What it is | Hybridizations |
  |---|---|---|
  | GSE45115 | 376 wild-type cultures grown in parallel with the 1,484 mutants, to measure day-to-day variance; dye-swap pairs `wt-SV###-a/b` | 755 |
  | GSE42215 | wild-type pool, Erlenmeyer protocol (100 wt vs refpool plus 100 dye-swapped) | 200 |
  | GSE42217 | wild-type pool, plate-reader protocol | 200 |
  | GSE42240 / GSE42241 | MATa wild types | 8 / 20 |

  - Per-hybridization values (VALUE = normalized log2 ratio test/ref, with the same print-tip LOESS and gene-specific dye correction) are in the GEO **sample tables**, checked for GSM630359, and in the `_RAW.tar` archives (GSE45115: 522 MB; GSE42215: 138 MB).
  - The GSE45115 series-matrix table is empty, so the wild-type matrix must be assembled from the sample tables.
  - The authors' site (deleteome.holstegelab.nl/downloads.php) offers `deleteome_all_mutants_controls.txt` (M, A and p for all mutants plus three control experiments), but not the wild-type replicate set.
- **Primary sample set: the 700 responsive mutants (GSE42527).**
  - Rationale, fixed now and label-free: profiles with three or fewer significant changes are wild-type-like noise as columns of a gene × sample matrix.
  - The K562 compendium likewise used significant knockdowns only, and N ≈ 700 matches K562 network half A (707).
  - All 1,484 mutants: descriptive sensitivity only.
- **Preprocessing (`analysis/yeast_prep.R`, stage `matrix`):**
  1. Map probes to systematic ORF names through the GPL11232 annotation.
  2. Average duplicate probes of one ORF.
  3. Mask each deleted gene's own entry in its mutant's profile and replace it with that gene's median over the other mutants. The deleted ORF is parsed from the column label; labels that cannot be resolved are logged.
  4. Drop ORFs with missing values in more than 5% of profiles; impute the rest with the gene median.
  5. Drop constant rows.
  6. No batch, slow-growth or SVD correction. The authors' growth-corrected file is not used.
- **Gene universe:** the rows that remain (expected about 6,000).

### Feasibility gate F (label-free, uses the wild-type controls; computed and committed before any truth is downloaded)

This answers the PI's concern about compendia without calibrating controls, such as Rosetta (Hughes 2000).

1. **Matrices.**
   - WT matrix: one profile per wild-type culture in GSE45115. Each profile is the mean of its dye-swap pair, oriented to log2(wt/ref), over duplicate spots, after the same probe-to-ORF mapping and gene filter.
   - **Implementation check:** the same rule, applied to the sample tables of the first 10 GSE42527 mutants (alphabetical), must reproduce the published M with Pearson r ≥ 0.99 and median |difference| ≤ 0.02.
   - If the check fails, the gate uses mutant profiles rebuilt by the rule for every mutant, so both arms are on one scale. The confirmatory matrix stays the published M.
2. **Statistic.** Cross-half reproducibility of the network.
   - Draw N_w = number of WT profiles (about 376) from the 700 mutants (10 draws, seed 20261104).
   - Split each draw at random into halves.
   - Fit Stouffer CLR on each half (`build_fit` as `analysis/lab_splits.R`).
   - R = |top-1% pairs of half 1 ∩ top-1% pairs of half 2| / |top-1% pairs|.
   - Do the same on 10 random halvings of the WT profiles.
3. **Gate.**
   - **Passes** if median R_mut ≥ 2 × median R_wt and min R_mut > max R_wt.
   - Also reported, not gating: R_wt itself (technical and day-to-day structure); the share of full-data mutant top-1% CLR pairs that are also top-1% in the full WT network; the median number of genes per mutant beyond the gene's WT 0.5th–99.5th percentile range.
   - If F fails, the result is **verdict 1** and nothing further is done.

### Truths (downloaded only after this file is fixed; structure counted before any scoring)

- **Yeast complexes (P1, S1).**
  - Source: EBI Complex Portal `complextab/559292.tsv`, the release current on the download date (the 1 Nov listing shows 2026-10-17, 913 KB), with UniProt 559292 accession → ordered-locus (ORF) mapping.
  - Build: flattening, non-protein removal and isoform handling exactly as `analysis/build_ecoli_complex_truth.R`; size filter 3–200 in-universe genes; overlap merge (absorb proper subsets, then complete linkage at Jaccard ≥ 0.5). There is no operon rule in yeast and no adjacency exclusion.
  - **P1 truth:**
    - universe: all pairs of genes in some merged complex;
    - positives: co-members of a merged complex with 3–6 in-universe genes (condition C1's size class, unchanged);
    - excluded: pairs that are co-members only of larger merged complexes;
    - negatives: the rest.
  - **S1 truth:** all merged complexes, positives = co-members of any merged complex.
  - **Fallback, fixed now.** If n_eff(P1) < 100 with Complex Portal, use the union of Complex Portal and CYC2008 (Pu et al. 2009; ORF lists), merged by the same rule, decided from counts alone before scoring. If CYC2008 cannot be obtained, or the union still has n_eff < 100, the result is **verdict 1**. The unused source is never scored.
- **Yeast regulons (P2, H1).**
  - Source: YEASTRACT+, documented regulations, all environmental conditions, all TFs by all genes, exported on the download date. Record the date, the form parameters and the release.
  - **Primary:** TF–gene pairs supported by both DNA-binding ("direct") **and** expression ("indirect") evidence. This is the analogue of RegulonDB strong/confirmed. If the form cannot filter conjunctively, export the two evidence classes separately and intersect them.
  - Co-membership as `make_regulons()` / `comembership_pairs()`: regulator-agnostic, regulons of 5–500 in-universe targets.
  - All-evidence regulons: descriptive.
- **Yeast TF–target edges (S6).** Universe: pairs with exactly one regulator gene (as in `docs/plans/asymmetry_expert.md`). Positives: the primary YEASTRACT edges.
- **Human (S2–S5).** The existing CORUM file (`data/human_truth/complexes.tsv`, as `corum_sets()`, then the overlap merge, then 3–6-gene merged complexes as P1) and CollecTRI (`data/human_truth/collectri.tsv`, as in K562).
- **Evaluability.** Each confirmatory truth needs ≥ 200 positive pairs (P1 and S1 also need ≥ 30 merged complexes with positives). If a truth fails, its test is not run, still counts in its Holm family, and is reported.

## Operators (fixed)

**Upstream, unchanged.** B-spline MI (`estimate_mi(bins = "hg", transform = "none")`) and Stouffer CLR (`calibrate(method = "normal", combine = "stouffer")`) as `build_fit` in `analysis/lab_splits.R`.

**K10 (comparator).** A = per-row softmax over the top 50 positive CLR scores at perplexity 10 (`clr:::.softmax_rows(C, keff = 10, cap = 50)`), with empty rows given a self-loop.

**R, reciprocal attention (candidate; no free parameter).**

1. B = A ⊙ 1[Aᵀ > 0], with the diagonal of B set to 0.
2. Rows of B with no remaining mass get a self-loop.
3. B_ii ← B_ii + max(0, 1 − Σⱼ B_ij).

This is `recip()` in `analysis/synthetic_width/sim_width.R` (line 259). The new script must reproduce the RECIP rows of `analysis/synthetic_width/out/round2/runs/std_f10_b1_s1_metrics.csv` to 1e-8 before any real data are fitted (implementation check; synthetic only).

**Walk and readout, both operators.** P = ½(I + B) (α = 0.5). S = (Pᵗ + Pᵗᵀ)/2 with diagonal 0. **Attention readout only.**

**CLR (H1 comparator).** The Stouffer CLR matrix that both operators are built on.

## Depth (label-free, fixed)

- **Splits.** 10 random halves of the 700 mutants (seed 20261101; written by `yeast_prep.R`, stage `splits`, and committed). Each split is used in both directions.
- **Criterion A** (as `analysis/keff_test.R`): the AP of S on the training half for the other half's top-1% CLR pairs.
  - Grid {1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32}.
  - For each split-direction, the argmax.
  - For each split, the geometric mean of the two directions.
  - The median over splits, snapped to the grid.
- **Each operator gets its own t\*:** t\*_K10 and t\*_R.
- **Primary comparison:** R at t\*_R vs K10 at t\*_K10. **H1:** K10 at t\*_K10 vs CLR.
- **Secondary (descriptive):** both at t\*_K10; both at t = 9.
- **Flags** (reported, not changing anything): t\* at a grid edge; t\*_R > 2·t\*_K10.
- **No truth curve over t** is computed for yeast.

For HepG2 and Jurkat, the same procedure applies on 10 random halves of the significant-knockdown pseudobulk (built as `analysis/perturbseq_bulk.R` from `GSE264667_{hepg2,jurkat}_raw_singlecell_01.h5ad`; the half-A/half-B convention and seed 20261015 as for K562 and RPE1).

## Endpoints

All endpoints are pooled co-membership AUPR, R minus K10 unless stated, at the primary depths.

**Primary (adoption), an intersection–union test (Berger 1982), each test at one-sided α = 0.025:**

| Code | Truth | Hypothesis |
|---|---|---|
| P1 | yeast small complexes (3–6, merged) | superiority: p < 0.025 **and** d ≥ +0.005 |
| P2 | yeast regulons (YEASTRACT, DNA-binding and expression evidence) | non-inferiority, margin 0.005 |

**Headline replication, a separate family H at one-sided α = 0.025, tested whatever P1 and P2 show:**

| Code | Truth | Contrast | Hypothesis |
|---|---|---|---|
| H1 | yeast regulons (P2's truth) | K10 minus CLR | superiority: p < 0.025 and d ≥ +0.005 |

**Key secondary, tested only if P1 and P2 both pass. Holm over S1–S6, one-sided familywise α = 0.025. A test that cannot be run still counts in the family:**

| Code | Data | Truth | Hypothesis |
|---|---|---|---|
| S1 | yeast | all merged complexes | non-inferiority, margin 0.005 |
| S2 | HepG2 | CORUM small (3–6, merged) | superiority |
| S3 | Jurkat | CORUM small (3–6, merged) | superiority |
| S4 | HepG2 | CollecTRI | non-inferiority, margin 0.005 |
| S5 | Jurkat | CollecTRI | non-inferiority, margin 0.005 |
| S6 | yeast | TF–target edges | non-inferiority, margin 0.005 |

**Descriptive (intervals only, never p-values):**

- per-complex AP (complexes with ≥ 3 positive pairs);
- the ribosome alone (merged complexes with ≥ 30% cytoplasmic r-protein genes, RPL/RPS ORFs);
- small complexes of 3–10 genes;
- all-evidence regulons;
- the secondary depths;
- all 1,484 mutants;
- R vs CLR on complexes;
- K562 essential half A (CORUM small, CORUM all, CollecTRI; a fresh operator on a used truth);
- *E. coli* M3D and EcoMAC SC regulons (no-harm);
- the label-free diagnostics below.

**Label-free diagnostics** (reported with the selections, before any truth), per dataset including M3D and K562:

- the share of K10 off-diagonal mass that R removes, and its row distribution;
- the share of rows left as pure self-loops;
- the reproducible-partner count from the 5-fold k-means cluster split (design memo, section 8.3);
- the Spearman correlation between the mass removed and that count.

## Statistics

- **AUPR:** `ap_fast()` / `pr_summary()` in `analysis/regulondb.R`. Tie blocks are counted at their end, as everywhere in the paper; Supplementary Note 16 covers the sensitivity to tie conventions, and the expected-tie AP is reported descriptively.
- **Gene jackknife.**
  - The paired delete-half gene jackknife, 100 replicates: in each, a random half of the truth's genes is kept, and the same halves serve every operator.
  - Seed 20261102, set once per dataset run; truths are drawn in endpoint-table order; never rerun with another seed.
  - SE_J = the SD of the 100 paired differences.
- **Experiment bootstrap (yeast only).**
  - Paired cluster bootstrap over the 700 mutant profiles, drawn with replacement, 200 draws, seeds 20261103 + b (b = 1, …, 200).
  - MI, CLR, K10 and R are refitted in every draw at the committed t\* values; truths are fixed.
  - SE_B = the SD of the 200 paired differences.
- **Deciding statistic (yeast P1, P2, H1, S1, S6).** SE = √(SE_J² + SE_B²): genes in the truth and experiments in the compendium are independent sources of variance that add (interpretation memo, section 3.4).
  - superiority: p = 1 − Φ(d/SE);
  - non-inferiority: p = 1 − Φ((d + 0.005)/SE).
  - The bootstrap percentile interval is reported descriptively. No bootstrap-only p-value is computed.
- **Human tests (S2–S5).** Jackknife SE only. The cost of a bootstrap is not justified for secondary tests.
- **No pooling** across datasets or across families.
- **Smallest effect of interest:** 0.005 AUPR.

## Power (structure-only gate; computed and committed before scoring)

- **Effective number of complexes.** n_eff = (Σₖ pₖ)² / Σₖ pₖ², where pₖ is the number of positive pairs of merged complex k in the P1 truth (after merging and restriction to the universe).
- **Projected SE.** SE_proj = 0.036/√n_eff × √(1 + 0.44²).
  - The anchor is the *E. coli* complex arm of CLR(log Pᵗ): SE 0.018 at n_eff ≈ 4.
  - The factor √(1 + 0.44²) is the bootstrap/jackknife SE ratio observed there.
  - This is conservative: mild K10 tweaks had about a third of a3r's paired SE on K562 CORUM.
- **Gate:** n_eff ≥ 100, so SE_proj ≤ 0.0039 and the MDE at 80% power ≤ 0.011.
  - Power at a true d = +0.010 is 0.72; at +0.015 it is 0.97.
  - The regulon NI arm, at a paired SE of 0.0015 (keff tweaks: 0.0005–0.0022), has power 0.86 at d = 0.
- **If the gate fails** (after the CYC2008 fallback): verdict 1, and no score is computed.
- **Planning effect: +0.010 AUPR.** This is about ×1.1–1.2 at a K10 small-complex AUPR of 0.05–0.10, a heavy discount on the smallest synthetic gain (×1.86). The test is not powered for the minimum effect of 0.005 (power about 0.25). A failure to confirm an effect between 0.005 and 0.010 is therefore expected and is not evidence of no effect.

## Decision rules and verdicts

The verdict is computed by the script and written to `verdict.csv` before anyone interprets the numbers. **Exactly one of verdicts 1–5 applies.**

| # | Condition | Verdict (one line) | Reporting |
|---|---|---|---|
| 1 | Feasibility gate F fails, or the power gate fails after the fallback, or a primary truth fails evaluability | **Not run:** no adequately powered held-out test was possible | Reciprocity stays a synthetic-only, explained candidate; the gate numbers are reported |
| 2 | P1 and P2 pass, and S1 non-inferiority is shown (Holm) | **Confirmed, no large-complex cost:** reciprocal attention becomes the method's operator from here on | The *E. coli* headline keeps K10, and the order of events is stated |
| 3 | P1 and P2 pass, and S1 is not shown | **Confirmed as a complex-mode variant:** it improves small complexes, with the large-complex limitation documented | K10 stays the default; R is reported as the confirmed variant for complexes |
| 4 | P1 fails, whatever P2 shows | **Not confirmed** | Wording as in CLR(log Pᵗ): "not confirmed", with the estimate and interval, never "no effect"; a significant regulon loss (d ≤ −0.005, p < 0.025 for d < 0) is added if present |
| 5 | P1 passes and P2 fails | **Regulon trade-off:** small complexes improve at a regulon cost | Not adopted; reported as a trade-off |

- **Family H, reported alongside:**
  - H1 holds: "the headline replicates in yeast".
  - Otherwise: "attention did not beat CLR on yeast regulons", with the interval. A significant CLR advantage (d ≤ −0.005, p < 0.025) is said plainly. H1 never changes verdicts 1–5.
  - If F passes but a regulon truth fails evaluability, H1 cannot be tested; this is reported.
- **Human lines (only under verdict 2 or 3):** S2/S3 give "replicated in HepG2 / Jurkat" or "not replicated"; S4/S5 give regulon no-harm in human. Under any other verdict, human results are descriptive intervals only.
- **What does not count:**
  - a win at any depth other than the primary ones;
  - a win on all 1,484 mutants, on the 3–10 size class, on per-complex AP or on all-evidence regulons when the primary fails;
  - a win after swapping the complex source (except the fixed fallback), the regulon evidence filter, the size class or the merge rule;
  - a win by a gate variant;
  - K562, RPE1 or *E. coli* results as support;
  - pooling yeast with human data;
  - results within ±0.005;
  - any p-value for an untested gated contrast (S1–S6 when P1 and P2 fail).

## Protocol and gates (in order)

1. **Review, then fix this file.** Replace the DRAFT block and commit. All later scripts check that the file is tracked, identical to HEAD (`git diff --quiet HEAD -- docs/preregistration_reciprocity_yeast.md`) and free of `STATUS: DRAFT`. They write the hash of the last commit touching the file into every output.
2. **Download** (`analysis/yeast_prep.R download`): expression and wild-type data only, with SHA-256 for every file in `data/yeast/SHA256SUMS`. Commit the hashes.
3. **Matrices** (`matrix`, `wt`, `check`): the processed mutant and WT matrices, their SHA-256, the unresolved-label log and the implementation-check table. Commit.
4. **Splits and feasibility** (`splits`, `feasibility`): split assignments in `docs/reciprocity_yeast_selections/` and the F statistics with a pass/fail line. Commit. **Stop if F fails.**
5. **The synthetic implementation check of R** passes (test script). Commit.
6. **Label-free selections**: t\*_K10 and t\*_R per dataset (yeast, HepG2, Jurkat, K562), the label-free diagnostics, and the bootstrap seeds, in `docs/reciprocity_yeast_selections/*_selections.csv`. Commit. The truth phase refuses to run unless the selections are byte-identical to the committed copies.
7. **Truth build and power gate** (`yeast_prep.R truth_structure`; refuses to run before step 1):
   - download Complex Portal 559292 and the UniProt mapping; take the YEASTRACT export (manual web form, saved with a parameter log);
   - build the member lists;
   - count positives, complexes and n_eff for P1, S1 and P2;
   - apply the fallback rule;
   - record SHA-256 of raw and built files;
   - the script reads membership only and never loads a network.

   Commit the counts and the power-gate line. **Stop if it fails** (verdict 1).
8. **Code review** of the test script against this file, by a reviewer who has seen no yeast or HepG2 score.
9. **Truth phase, once:** yeast (P1, P2, H1, S1, S6; jackknife and bootstrap), then HepG2 and Jurkat, then the descriptive sets. The first complete run is the confirmatory one. A rerun is allowed only after a crash, or for a bug shown to violate a named line of this file; both runs and the reason are reported.

## Compute (for `jobs_pc/sched`)

| Block | Core-hours | Memory |
|---|---|---|
| Yeast: matrix and WT build | under 1 | under 4 GB |
| Yeast: feasibility, 40 half-fits at N ≈ 188 | 2–4 | under 4 GB |
| Yeast: t\* (20 half-fits plus walks for two operators) | 15–25 | under 4 GB |
| Yeast: jackknife | 2–4 | under 4 GB |
| Yeast: 200-draw bootstrap | 60–80 | 4–6 GB per worker |
| HepG2 + Jurkat: pseudobulk build from the single-cell h5ad files (HepG2 5.6 GB; Jurkat size not retrieved) | 2–4 | about 60 GB, one job |
| HepG2 + Jurkat: t\*, walks, jackknife | 70–90 | 10–12 GB per worker |
| K562 descriptive | about 10 | 6 GB |

- **Total:** about 170–220 core-hours.
- **Disk:** about 30 GB (raw single-cell files plus the wild-type archives).
- All jobs run with `# mem` and `# cpus` headers, at about 70–80% utilization.

## Change log

- 2026-11-01: DRAFT written for review. Not fixed.
