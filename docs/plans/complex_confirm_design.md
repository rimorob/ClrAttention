# Preregistration (DRAFT): does CLR(log Pᵗ) (a3r) beat CLR attention on protein complexes in an untouched yeast compendium?

Status: **DRAFT rev 1, 13 Nov 2026. NOT FIXED.** Three authorisation items (gate G0) must be recorded in the change log before the truth stage. The label-free stages (gates Ga and Gb) may run before that. Any change after FIXED is a dated change-log entry. A change after the first truth-stage run is a deviation and is reported as one.

*Drafted by Claude (Opus, design stage) for the PI, who is away until 14 Nov and delegated "draft the design and run the compute". An independent reviewer reads this file before any compute. A third agent runs it.*

- **No readout has been scored against any candidate truth.** No AUPR, enrichment or per-complex statistic of a3r, attention or CLR exists on yeast complexes. Only truth *structure* was computed: gene counts, complex counts, pair counts, size and overlap. For the power analysis, existing outputs of spent runs were re-read. Nothing was re-scored.
- **What the drafter has seen.**
  - The DREAM5 net4 gene-ID table (names only).
  - The first column (ORF names) of `data/yeast/deleteome_responsive_700.tsv`.
  - The Complex Portal yeast file and its structure counts (below).
  - The committed outputs of the 26 Oct test (`results/clr_logpt_test/truth/*/jack_replicates.csv` and `per_complex_ap.csv`, all spent truths).
  - The DREAM5 net4 truth-stage results for attention vs CLR on regulons (`docs/dream5_selections/truth_stage_log.md`, as reported in `docs/STATUS.md`).
  - No expression value of net4 or Kemmeren was read.
- **Verification markers.** Every claim that was not checked in this session is marked **[verify]**. The section "What was verified, and how" lists the checks that were made.

## Question

a3r is a fixed label-free transfer. Its operator is MI → CLR (Stouffer) → per-row softmax (keff 10, cap 50) → lazy walk α = 0.5, t = 9. Its readout is bilateral CLR on log P⁹, with row statistics over reachable pairs only, clipped at 0, combined by Stouffer.

Does a3r beat plain CLR attention on protein-complex co-membership in an expression compendium and a complex truth never used for any tuning, selection or descriptive look? Attention uses the same operator with the readout (P⁹ + P⁹ᵀ)/2. CLR is a secondary comparator: it is reported but plays no part in the decision.

**History that this test must carry in every write-up.**

1. a3r was chosen from about 20 readouts on K562 truths (`docs/plans/corum_analysis.md`). K562 is in-sample.
2. Its preregistered E. coli test (fixed 22 Oct, run 26 Oct) was **not confirmed** (`docs/clr_logpt_verdict/`):
   - P1 M3D +0.0257, jackknife SE 0.0182, one-sided p 0.079;
   - P2 EcoMAC +0.0285, SE 0.0174, p 0.051;
   - regulon non-inferiority passed.
3. The PI decided on 12 Nov that a3r is **not adopted** and that "any revival needs a new preregistration on an untouched, larger complex truth" (`docs/STATUS.md`).
4. Since then a3r was used in exploratory work on K562 and M3D only (stacked layers, `docs/explore_stack_results.md`).

This file is that new preregistration. **It is proposed as the last confirmatory test of a3r** (interpretation memo section 7.3, item 7). If it does not confirm, a3r stays an explanatory ablation for good.

## G0. Prior commitments that conflict with this test (must be resolved by the PI before the truth stage)

These two items are binding text in a FIXED preregistration. The design cannot lift them by itself.

1. **"Yeast is single-use" (`docs/preregistration_reciprocity_yeast.md`, disclosure table, fixed 1 Nov).** The text reads: "The a3r follow-up (interpretation memo section 7.3) and the C3 asymmetry edge test are **not** included and therefore can never use these yeast truths (PI decision)." "These yeast truths" were three:
   - the EBI Complex Portal 559292 complexes, with CYC2008 as fallback;
   - YEASTRACT regulons;
   - the TF-edge truth.

   The primary truth below is that same Complex Portal source, used on a different expression matrix. Two readings exist:
   - *(a) Premise reading.* The exclusion followed from single use. The reciprocity test was supposed to consume the truths. Gate F failed on 2 Nov, so no yeast truth was downloaded or consumed (verdict 1). The same cell of that file then says that an unconsumed complex truth "may be used later only under a new preregistration that discloses three things". Under this reading the a3r test is allowed with that disclosure. Two of the three things are moot: the complex size structure was never computed there, and H1 was never run. The third is given here: the reciprocity arms were "not run".
   - *(b) Literal reading.* "Can never use" is absolute.

   **Only the PI can choose between these readings.** The drafter's reading is (a). It is not assumed.
2. **Serial testing.** The reciprocity preregistration also states: "This is the third complex-oriented operator test … It is declared the last. Whatever the verdict, no further budget- or readout-level operator change aimed at complexes is tested confirmatorily in this project."
   - This test would be the **fourth** complex-oriented confirmatory test, after keff, a3r on E. coli and reciprocity (not run). It would also be the **second** confirmatory test of a3r itself.
   - The PI's 12 Nov sentence ("any revival needs a new preregistration …") reopens the question, but it does not explicitly withdraw the 1 Nov declaration.
   - The PI must withdraw that declaration for this one test, in writing, and every write-up must disclose the serial testing.
3. **Power below the memo's own threshold** (see "Power"). The interpretation memo (section 7.3, item 2) says the follow-up should be confirmatory only if the projected SE is at most 0.009. No projection reaches that. The best-calibrated projection is 0.0094 to 0.0106. The PI must accept this explicitly, or downgrade the test to a preregistered estimation study. The verdict categories are the same either way, but "CONFIRMED" would then be reworded as "estimate excludes 0 (preregistered estimation)".

**G0 passes** only when a dated change-log entry records the PI's decision on items 1–3, in the PI's own words or as a quotation of the PI's delegation message. Gates Ga and Gb, which are label-free, do not need G0. **Recommendation to the lead:** run Ga and Gb tonight if wanted. Hold the truth stage until the PI confirms G0 on 14 Nov. The cost of waiting is one day. The cost of a test that later has to be disowned is the whole test.

## Candidates (feasibility, independence, power), structure only

For all candidates the filters are those of the spent truths:

- recursive flattening of sub-complexes;
- non-protein participants dropped;
- size 3–200 in-universe genes;
- overlap merge: exact duplicates and proper subsets dropped, then complete linkage on 1 − Jaccard cut at 0.5, union;
- the code is `analysis/build_ecoli_complex_truth.R`, copied into `analysis/complex_confirm/build_yeast_complex_truth.R`.

Definitions used in the table:

- **n_eff** is the inverse Herfindahl index of positive-pair shares across merged complexes. It is the effective number of independent complexes for a pooled pair metric; the K562 loss, for example, was "one module counted six times".
- **Components** are connected components of the gene-sharing graph of merged complexes, a stricter count of independent complexes.
- **SE** is the projected paired jackknife SE of a3r − attention (see "Power").

| | Genes (universe) | Merged complexes with positives | Truth genes | Pairs | Positives | Components | n_eff | Largest share | Ribosome share | Projected SE | Feasible tonight? | Independence |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **A. DREAM5 net4 + Complex Portal 559292** | 5,667 ORFs (+283 decoys kept in the operator); 536 chips, 321 replicate groups | **262** | 1,602 | 1,282,353 | **13,439** | 180 | **14.5** (15.1 anchor-comparable) | 0.188 (60S) | 0.289 | **0.0094–0.0106** (central 0.0099); surrogate range 0.0105–0.0226 | yes: data on disk (DREAM5 Synapse files, checksummed), truth 0.9 MB | Expression never used for any a3r tuning or look. Same matrix as the DREAM5 Y-M/Y-E test (attention vs CLR on binding regulons and edges, 11 Nov); a3r never computed on it. Truth never scored. G0 item 1 applies. |
| B. Kemmeren responsive-700 deleteome + the same truth | 6,169 ORFs; 700 mutants | 260 | 1,600 | 1,279,157 | 13,225 | 181 | 15.6 | 0.171 | 0.273 | 0.0091–0.0102 | yes: on disk (`data/yeast`, matrix SHA-256 `8cace0e6…`) | Never scored. Gate F of the reciprocity test failed: mutant chance-corrected reproducibility 0.122 vs wild type 0.296, a label-free warning about this matrix. G0 item 1 applies with full force (it is the matrix that preregistration named). The wild types alone (GSE45115, 376 cultures) carry mostly technical and day effects and are not a useful co-expression compendium. |
| B′. Kemmeren wild types only | about 6,170 | as B | | | | | | | | not projected | on disk | As B, plus little biological variation. Rejected. |
| C. Human, not K562/RPE1: HepG2 and Jurkat Perturb-seq (Nadig et al., GSE264667) with CORUM or human Complex Portal | about 8,000 [verify] | | | | | | | | | not projected (no universe without the data) | **no**: NCBI FTP lists only raw single-cell h5ad files, HepG2 5.2 GB and Jurkat 8.7 GB, plus a 15 GB RAW.tar (verified 13 Nov). A pseudobulk pipeline would have to be written and preregistered first. A figshare pseudobulk release [verify]. | Same lab, design and essential-gene panel family as K562. CORUM content is the ribosome again. The interpretation memo counts this as "barely" independent. Reactome and co-essentiality are pathway and module truths, not complex co-membership. |
| D. *B. subtilis* | | | | | | | | | | | **no** | Complex Portal has no taxon 224308 file (404 on 13 Nov; absent from the `complextab/` listing). SubtiWiki complexes [verify] would be small, likely E. coli-sized and underpowered. Compendia (Nicolas et al. 2012; Arrieta-Ortiz et al. 2015) [verify]. |
| CYC2008 as truth (for A or B) | | | | | | | | | | | **no** | `wodaklab.org/cyc2008/resources/CYC2008_complex.tab`: connection timed out twice, over http and https, 13 Nov. Not used. No third-party mirror is used (no provenance). |

**Choice.**

- **Primary: A** (DREAM5 net4 replicate-level expression, Complex Portal 559292 merged complexes).
- **Pre-declared fallback: B** (Kemmeren responsive-700 matrix, the same truth built on its universe). It is used only if A fails the label-free feasibility gate Gb. B is never run in addition to A, and never after any A score exists.
- **Why A:**
  - the data are already on disk and were never touched by a3r;
  - the truth has 262 merged complexes and an n_eff of about 15, against 54 complexes and an n_eff of 4.4 in E. coli, so the projected SE is about half the E. coli SE;
  - ribosome dominance is lower (29% of positives, against about 60% in E. coli);
  - the operator outputs for CLR and attention already exist with SHA-256 sidecars, so the implementation gate can be exact.
- **Why not C or D:** infeasible without new downloads and new pipelines, with weak independence (C) or a small truth (D).

## Data (candidate A)

- **Expression.** `data/dream5/synapse/net4/net4_expression_data.tsv`: 536 chips × 5,950 genes, replicate level. This is the primary input of the DREAM5 design (PI decision 3, 3 Nov). Checksums are in `data/dream5/synapse/SHA256SUMS`, verified there on 4 Nov.
  - No normalisation, log transform, batch correction or filtering is applied. There are no zero-variance genes.
  - **Decoys are kept** in the operator: 283 genes named `decoyN`, as in the DREAM5 design. They cannot be truth genes.
  - The replicate-averaged matrix (321 groups) is a descriptive sensitivity only.
- **Provenance and credentials.** The files came from the PI's earlier staged Synapse download (`analysis/dream5_download.py`). **No Synapse credential is used, requested or recorded by this test.** Only files already on disk are read.
- **Gene universe for the truth:** the 5,667 non-decoy ORF names of `net4_gene_ids.tsv`.

## Truth (fixed now; built, frozen and committed before any score)

- **Source.** EBI Complex Portal, *S. cerevisiae* S288C curated set.
  - File: `https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/559292.tsv`, downloaded 13 Nov 2026.
  - Server Last-Modified: 17 Oct 2026. 649 complexes.
  - SHA-256 `c104ffe324f69c3e9820e9b7ab9b648cd7177b390e402529ed7446ae59ab7ec5`.
  - The `559292_predicted.tsv` file is not used.
- **Accession map.** UniProt accession → ordered locus (ORF), from `https://rest.uniprot.org/uniprotkb/stream?query=organism_id:559292&fields=accession,gene_oln,gene_names,gene_primary,reviewed`.
  - Release 2026_04 (4 Nov 2026), 6,733 entries.
  - SHA-256 `f8e4cc8f6fb793922dea27df658baa6cc5178af8f059d748d364c260302a4aee`.
  - 1 of 2,936 flattened participants has no ORF.
- **Raw files** are in `data/complex_confirm_raw/` (gitignored, 1.3 MB; `SHA256SUMS`).
- **Build.** `analysis/complex_confirm/build_yeast_complex_truth.R` (committed with this file). The rules of the 26 Oct E. coli truth apply, plus two yeast-specific rules fixed before any score:
  1. **Paralog molecule sets.** The two cytoplasmic ribosome subunits (CPX-1599 40S, CPX-1601 60S) list interchangeable paralogous subunits as sets, `[ACC1,ACC2](n)`. No other complex in this release does, and the E. coli file has none. Each set is **expanded to all its members**.
     - Without this rule the cytoplasmic ribosome silently drops out: in a first count it was 0.8% of positives instead of 29%. This was found while counting structure, before any score existed, and it is disclosed.
  2. **Same-set paralog pairs are excluded** from the pair universe (48 pairs on net4). Near-identical paralogs are trivially co-measured on arrays (cross-hybridisation). This is the yeast analogue of the same-operon exclusion of the E. coli truth.
     - Implementation in the run: pass the paralog-set map as the `operon_of` argument of `comembership_pairs()`, so that same-set pairs are removed exactly as same-operon pairs are.
- **Frozen member list.** `docs/complex_confirm_selections/yeast_complexes_members.tsv`, SHA-256 `146030f76b57310049440c265a3d190f2a1f892558dc5ec2faab2f94433e6e3e` (`yeast_complexes.sha256`). The truth stage refuses to run on any mismatch. Counts are in `yeast_complexes_counts.csv`, per-complex sizes in `yeast_complexes_sizes.csv`, and coverage in `yeast_complexes_source_log.csv`.
- **Primary truth ("main"):**
  - universe: all pairs of genes in at least one merged complex, minus same-set paralog pairs;
  - positives: pairs that share a merged complex;
  - net4 numbers: 262 merged complexes, 1,602 genes, 1,282,353 pairs, 13,439 positives, base rate 0.0105;
  - 98.7% of the ORFs in complexes with at least 3 ORFs are in the net4 universe.
- **Pre-declared ribosome list ("ribosome-dominated").** Merged complexes with at least 30% cytoplasmic r-protein ORFs. The r-protein ORFs are those whose SGD standard name matches `^RP[LS][0-9]+[AB]?$` (`data/yeast/SGD_features.tab`; the rule of the reciprocity preregistration; 124 such ORFs in the net4 universe).
  - On net4 the list is exactly two complexes, CPX-1601 (72 genes) and CPX-1599 (53 genes), carrying 3,886 positives (28.9%).
  - Mitochondrial ribosome subunits (CPX-1602, CPX-1603) and the small-subunit processome (CPX-1604) do not match the r-protein rule and stay in the noRP truth. This matches the K562 noRP rule, which did not remove MRPL/MRPS genes.
- **Variants (sensitivity truths, built in the same file):**
  - noribo: 260 complexes, 9,553 positives;
  - nomerge: 378 complexes, 13,298 positives;
  - size2 (minimum size 2): 480 complexes, 13,657 positives;
  - small 3–6 genes: 198 complexes, 1,253 positives;
  - large ≥ 20 genes: 13 complexes, 9,531 positives;
  - the r-protein module (124 genes, descriptive coherence);
  - paralog pairs included (sensitivity, built at run time).

## Methods (fixed; nothing is tuned on net4)

All three readouts are computed from one operator, built once on the full 536 × 5,950 matrix.

| Readout | Definition | Code |
|---|---|---|
| CLR | `estimate_mi(bins = "hg", transform = "none")`, `calibrate("normal", "stouffer")`, truncated at 0 | as the DREAM5 run (`analysis/dream5_run.R`, method `clr`) |
| **attention (comparator)** | `build_operator(alpha = 0.5, softmax_keff = 10, softmax_cap = 50)` on CLR, lazy walk, **t = 9**, S = (P⁹ + P⁹ᵀ)/2, zero diagonal | as DREAM5 `attention`; `analysis/clr_logpt_test.R::attention` |
| **a3r (tested)** | L = log P⁹ (zeros → NA, diagonal NA); each row standardised over its reachable off-diagonal entries (mean, SD with n − 1); negative, missing and non-finite z → 0; S = (Z + Zᵀ)/√2, zero diagonal | `analysis/explore_stack/stack_lib.R::sk_a3r`, unchanged (byte-identical to commit d8fdb58; checked by Ga) |

- t = 9 is the transfer value. **No depth rule is applied.** net4 has no lab structure and no established label-free depth rule for this readout. The a3r t\* of the 26 Oct test was secondary and is not used.
- **No new hyperparameter.**
- Hub stratifier H_i (MI row mean) = (rowSums(M) − M_ii)/(G − 1) on the raw MI matrix M of all 5,950 genes. It is used as in `docs/robust_z_results.md`: 10 strata of H_i + H_j; within-stratum percentile (`rz_hub_stratify`); both methods under the same adjustment.

## Power (structure plus existing outputs only)

### Anchors

The four spent complex truths of the 26 Oct test have known jackknife SEs for a3r − attention (`jack_replicates.csv`, 100 replicates) and known per-complex structure (`per_complex_ap.csv`: sizes and positive pairs). The script `analysis/complex_confirm/power_projection.R` produces `docs/complex_confirm_selections/power_anchors.csv`, `power_projection.csv` and `power_summary.csv`.

| Anchor | Truth | Merged complexes | n_eff | SE a3r − att | k = SE·√n_eff | SE a3r − CLR | k_CLR |
|---|---|---|---|---|---|---|---|
| M3D | Complex Portal 83333, merged | 129 | 4.44 | 0.0182 | 0.0383 | 0.0130 | 0.0274 |
| EcoMAC | same | 124 | 4.42 | 0.0174 | 0.0367 | 0.0127 | 0.0267 |
| K562 essential | CORUM merged | 634 | 10.35 | 0.0128 | 0.0410 | 0.0050 | 0.0161 |
| RPE1 | CORUM merged | 678 | 10.44 | 0.0116 | 0.0375 | 0.0054 | 0.0176 |

The jackknife SE of the a3r − attention difference follows **SE ≈ k / √n_eff** with k = 0.037 to 0.041, a spread of about ±5%, across two organisms, three truths and n_eff from 4.4 to 10.4. The number of complexes alone does not predict the SE (K562 has five times as many complexes as E. coli but only 30% less SE). Concentration does: the jackknife moves whole blocks of a few large complexes in and out. That concentration is what made the E. coli test weak.

### Projection for the primary truth

net4 main has n_eff = 15.1, using the anchor definition (positive pairs per complex counted with overlap, as in `per_complex_ap.csv`), or 14.5 with fractional assignment.

| Quantity | Best-calibrated projection (k rule) | Surrogate models (less well calibrated) |
|---|---|---|
| SE of a3r − attention | **0.0099** (0.0094–0.0106) | 0.0105–0.0226 (median of 8 projections 0.0145) |
| Minimum detectable effect, 80% power, one-sided α = 0.025 (2.80 × SE) | **0.028** (0.026–0.030) | 0.029–0.063 |
| Power at +0.005 / +0.010 / +0.015 | 0.07 / 0.17 / 0.33 | |
| **Power at +0.025** (about the E. coli estimate) | **0.72** (0.66–0.75) | **0.20–0.66; 0.41 at the median SE** |
| Power at +0.035 / +0.050 | 0.94 / 1.00 | 0.34–0.92 / 0.60–1.00 |
| SE of a3r − CLR (secondary) | 0.0041–0.0071 | 0.0040–0.0062 |
| With experiment variance added (SE_B ≈ 0.41 SE_J, E. coli ratio): combined SE | about 0.0107; power at +0.025 about 0.65 | |

The surrogate models simulate delete-half gene draws on the truth's complex sizes, with per-complex effects that are either iid or resampled from the anchors' per-complex AP differences within size class. They calibrate less consistently across the anchors (the scale factor varies by a factor of 1.6 to 2.1 between E. coli and human), so they are reported as the pessimistic envelope, not as the projection. They are in `power_projection.csv`.

### Plain statement

- **The test can reliably detect effects of about +0.03 AUPR or more.** It detects an effect of the E. coli size (+0.025 to +0.027) with probability about 0.7 under the best-calibrated projection, and with probability between 0.2 and 0.66 under the pessimistic surrogates.
- **Under the median of all projections, power at +0.025 is 0.41, below 60%.** The memo's target, SE ≤ 0.009 (80% power at +0.027), is **not reached** by any projection (best 0.0094). This is G0 item 3.
- The test **cannot** detect effects below about +0.015: power is at most 0.33 at +0.015, and 0.07 at the +0.005 minimum effect. A NOT CONFIRMED verdict therefore cannot say that a3r "does not help". It can only exclude effects above the interval's upper end.
- **The expected yeast effect may be smaller than the E. coli one.** The mechanism a3r repairs (dilution in large cliques) acts mostly on large complexes, and the ribosome carries 29% of yeast positives against about 60% in E. coli. Large complexes (≥ 20 genes) carry 71% of yeast positives: the cytoplasmic ribosome and, outside it, the mitochondrial ribosome, the small-subunit processome, the proteasome, the nuclear pore and spliceosomal complexes.
- **Options considered to raise power, and why they were rejected:**
  - *Pooling A and B.* Both share one truth, and the SE is dominated by truth sampling (the k rule, and the interpretation memo section 7.2), so pooling would barely shrink the SE. It would also put the restricted deleteome into play and re-use a matrix that failed gate F.
  - *Making noribo the primary.* Its SE is lower (n_eff 22, SE about 0.008), but it removes most of the mechanism's expected effect. On K562, noRP a3r − attention was only +0.006.
  - *Equivalence framing.* A ±0.005 margin needs an SE of about 0.0015. Impossible.
  - *Adopted framing:* superiority with a three-way verdict that separates "excluded the planning effect" from "inconclusive", plus estimation reporting.

## Estimands and decision rule (fixed before any score)

- **Primary endpoint (P).** d = AUPR(a3r) − AUPR(attention):
  - pooled co-membership AUPR on the primary truth (net4 main);
  - `ap_fast` with the repo tie convention (ties at the end of their block), raw scores, t = 9.
- **Test.** One-sided superiority at α = 0.025, with SE = SE_J from the paired delete-half gene jackknife:
  - p_sup = 1 − Φ(d/SE_J);
  - interval d ± 1.96 SE_J.
- **Multiplicity.** One primary test and no family. Every other quantity is an estimate with an interval. No p-value is reported for it, and no claim is made from it.

| Verdict | Condition | Licenses | Does not license |
|---|---|---|---|
| **CONFIRMED** | p_sup < 0.025 **and** d ≥ +0.005 | "In a preregistered test on an independent yeast compendium and complex truth, a3r exceeded attention on complex co-membership by d (95% CI …), after failing to confirm on E. coli (P1 p 0.079, P2 p 0.051)." Makes a3r *eligible* for adoption. Adoption is the PI's decision, with the E. coli result, the serial testing and the absence of an untouched yeast regulon no-harm arm disclosed. | "a3r improves regulons"; "no regulon cost" in yeast (no untouched regulon truth; see below); "a3r beats CLR on complexes" (secondary, no test); any pooling with E. coli as confirmation |
| **NOT CONFIRMED** | not CONFIRMED, **and** either d ≤ 0 or the 95% upper bound < +0.025 | "Not confirmed in yeast: d (95% CI …)". If the upper bound is below +0.025: "effects the size of the E. coli estimate are excluded". If p_loss = Φ(d/SE_J) < 0.025: say plainly that a3r was **worse** than attention. a3r stays an explanatory ablation for good (stopping rule). | "a3r has no effect"; "shown not to work" |
| **INCONCLUSIVE** | not CONFIRMED, d > 0, and the 95% upper bound ≥ +0.025 | "Not confirmed; the interval (…) includes both no effect and the E. coli-size effect." a3r stays an ablation (stopping rule). | "trend", "marginal", "nearly significant", "underpowered but likely real" |

- The verdict is written by the script to `verdict.csv` before anyone interprets it.
- **What does not count:** a win by a3r at any other t; a win on a sensitivity truth (noribo, nomerge, size2, paralog pairs included); on per-complex AP; under the hub stratifier; under the complex-level or experiment resampling; on the fallback matrix when A was evaluable; after swapping the complex source, the size, merge or paralog rule; or a pooled analysis with E. coli.

### SE: one decisive procedure, pre-declared

- **Decisive: SE_J**, the paired delete-half gene jackknife of the 26 Oct test (`clr_logpt_test.R::jack_dataset`).
  - Each replicate keeps a random half (floor) of the primary truth's genes (`pairs$genes_in`, 1,602 genes) and scores the pairs with both genes kept.
  - The same halves are used for all methods.
  - SE_J = SD of the paired replicate differences (delete-half makes the jackknife factor 1).
  - **500 replicates** (the E. coli test used 100), with seed `set.seed(20261114)`. The halves are drawn for the truths in the fixed table order below, right after setting the seed.
  - The larger count only cuts the Monte Carlo error of SE_J from about 7% to about 3%. The DREAM5 Y-M verdict was a boundary label within Monte Carlo error (stats review M1, 12 Nov). This does not change the estimand.
- **Why the gene jackknife, and not the experiment bootstrap** (the 0.018 vs 0.0075 gap of 26 Oct):
  - the two measure different variance. The gene jackknife measures which genes and complexes make up the benchmark. The experiment bootstrap measures which chips make up the compendium, with the truth fixed;
  - the claim is about complexes in general, and truth sampling was the larger component (about 2.4×);
  - the gene jackknife is the paper's standard and the procedure of the test being followed up;
  - the two components add up. They are not alternatives. Choosing the smaller one would be choosing the answer.
  - Known caveat (interpretation memo 3.2): for a pair statistic, the delete-half gene jackknife may overstate the purely pair-level component (heuristically about 3×), so it may be somewhat conservative. This is accepted and not corrected.
- **Sensitivity 1, complex-level resampling.** A delete-half jackknife over **merged complexes**:
  - each replicate keeps a random half (floor) of the 262 merged complexes;
  - the truth is restricted to the kept complexes: universe = pairs among their genes minus paralog pairs; positives = pairs that share a kept complex;
  - 500 replicates, `set.seed(20261115)`, paired; SE_C = SD of the paired differences.
  - The p-value and verdict under SE_C are written to `verdict.csv` as a labelled sensitivity line. They never change the verdict.
- **Sensitivity 2, experiment jackknife.** A delete-half jackknife over the 321 **replicate groups**, built from `Repeat` as in the DREAM5 design:
  - CLR, attention and a3r are refitted in every replicate; 100 replicates; seeds 20263000 + b, the same halves as the DREAM5 group-level refits;
  - SE_B = SD of the paired differences. The combined SE is √(SE_J² + SE_B²);
  - both are reported as sensitivity lines with the verdict recomputed under the combined SE, labelled, never deciding;
  - the CLR and attention refits are checked against `results/dream5/net4/refit/REFIT_MANIFEST.tsv` hashes (an implementation check).

### Secondary and descriptive quantities (intervals from SE_J halves; no p-values, no claims)

1. a3r − CLR, raw (the comparator the lead named; S1).
2. Hub-stratified a3r − attention and a3r − CLR (MI row-mean stratifier, as above).
3. Size strata: positives restricted to merged complexes of 3–6, 7–19 and ≥ 20 genes. The universe is the primary universe, and pairs that are co-members only of complexes outside the stratum are excluded from that stratum's scoring.
4. No-ribosome: the noribo truth (the pre-declared two-complex list removed).
5. Sensitivity truths: nomerge, size2, paralog pairs included.
6. Attention − CLR on the primary truth: is there a loss to rescue in yeast? This was the motivating premise. E. coli showed none.
7. Per-complex AP (macro), complexes with ≥ 10 positive pairs, as `per_complex_summary`.
8. **Regulon non-inferiority: none.** No untouched regulon truth exists for net4.
   - The DREAM5 net4 gold standard (binding plus conservation) was scored for attention vs CLR on 11 Nov.
   - YEASTRACT falls under G0 item 1 and would need its own contamination audit.
   - Descriptively only: a3r − attention on the DREAM5 Y-M regulons (TFs with 5–500 targets, as the DREAM5 design), labelled "touched truth; a3r never scored on it; no-harm description, not a test".
9. Replicate-averaged input (321 × 5,950): the primary contrast, descriptive.
10. Descriptive combination with E. coli: inverse-variance fixed-effect estimate of yeast d and the E. coli mean of P1 and P2 (d = 0.0271, SE 0.0178, the two treated as fully correlated). Labelled "includes a known, previously unconfirmed estimate; not a test". Never the verdict.

Truth table order (this fixes the jackknife draw order): main, noribo, nomerge, size2, main_paralog_incl, small_3to6, mid_7to19, large_ge20, then the DREAM5 Y-M regulons.

## Gates (in order; each writes a CSV; any failure stops the job before any truth is read)

- **Ga. Implementation equivalence (spent data only).** Tolerance 1e-9 unless stated, as in `docs/explore_stack_results.md` section 7.
  1. `stack_lib.R`, `ideas47_lib.R` and `robust_z_lib.R` are byte-identical to their versions at commit d8fdb58 (SHA-256 recorded at the run).
  2. **K562:** the new driver's code path reproduces, against `results/corum_expert2/k562/t9_jackknife_EXPLORATORY.csv`:
     - a3r CORUM_all AUPR 0.1054351, printing as 7.95×;
     - CLR 0.0847 (6.39×) and attention 0.0465 (3.51×).

     It also reproduces a3r − K10 point differences and jackknife intervals for DoRothEA A–B, CollecTRI, CORUM and KD response against `analysis/corum_exploratory/metrics_robust.csv`, with the CORUM line printing as +0.059 [+0.033, +0.084]. For CORUM noRP only the point difference is compared (explore_stack deviation 1).
  3. **M3D:** a3r complexes AUPR 0.19209 and SC 0.16222, and the P1 d and SE (0.025659, 0.018169), against `docs/clr_logpt_verdict/endpoints_all.csv`.
  4. **net4, label-free:** the CLR and attention score matrices from the new run equal the DREAM5 caches `results/dream5/net4/clr__replicate.scores` and `attention__replicate.scores` (upper-triangle doubles, SHA-256 in `SCORE_MANIFEST.tsv`). Equality means identical SHA-256 or max |diff| ≤ 1e-12.

  Ga may reuse the explore_stack K562 driver unchanged for parts 2–3, but the a3r function must be the one the net4 run calls.
- **Gb. Label-free feasibility on net4 (no truth enrichment, no score against any truth).** All must hold:
  1. the operator builds without error, and no row has an empty candidate list (no self-loop rows); count reported;
  2. unreachable share of off-diagonal P⁹ pairs ≤ 5% (K562 0.28%, M3D 0.07%; the explore_stack layers degraded from 5.7%);
  3. 0 rows with fewer than 2 reachable entries;
  4. 0 all-zero rows in the a3r matrix S;
  5. truth structure from the frozen member file: SHA-256 matches; coverage of complex ORFs ≥ 0.90 (0.987); merged complexes with positives ≥ 100 (262); n_eff ≥ 10 (14.5); positives ≥ 5,000 (13,439). These pass now.

  Reported, not gating: π max/median, the share of in-weight held by the top 1% of genes, the median candidate count and the σ_i quantiles (`operator_checks`). **If Gb fails, the fallback B is used** with identical definitions:
  - universe of 6,169 ORFs; frozen member rows `dataset == "kemmeren_700"`; paralog rule as above;
  - experiment jackknife over the committed mutant clusters in `docs/reciprocity_yeast_selections/mutant_clusters.csv`, or its stratified fallback rule.

  If B also fails Gb, the result is "not evaluable" and nothing is scored.
- **G0.** Recorded in the change log before the truth stage (see above).
- **Gc. Memory declarations** for `clr-sched`: 1.4 × best estimate; headers `# mem`, `# cpus`, `# after`; `gc()` in the parent before every fork (`docs/STATUS.md`, resource facts of 11 Nov). Past OOM kills happened at 16.4 GB with 14 GB declared, at 16.7 GB (robust-z M3D) and at forked bootstraps. Machine: 125 GB RAM; scheduler budget 90 GB.

| Job | What | Best estimate | Declared (`# mem`) | `# cpus` | `# after` |
|---|---|---|---|---|---|
| 6000_cc_gate_a | Ga parts 1–3 (K562 rebuild and scoring of spent truths, M3D) | 18 GB (explore_stack K562 peak 25.7 GB with 11 readouts; 3 readouts here) | 26 | 8 | – |
| 6010_cc_net4_labelfree | operator, P⁹, CLR/attention/a3r caches with SHA-256, Ga part 4, Gb | 6.5 GB (DREAM5 net4 truth stage peak 6.4 GB) | 10 | 8 | 6000 |
| 6020_cc_net4_truth | **after G0 only:** truth SHA gate, primary, 500 gene-jackknife and 500 complex-jackknife replicates, secondaries, verdict | 11 GB with 8 forked workers | 16 | 8 | 6010 |
| 6030_cc_net4_expjack | 100 replicate-group refits (3 readouts), SE_B, combined-SE sensitivity line | 8 GB (refit 1.5–1.9 min per replicate measured for DREAM5; about 3–5 h) | 12 | 8 | 6020 |
| 6040_cc_summarize | tables, `verdict.csv` final with the sensitivity lines | 2 GB | 4 | 2 | 6030 |

## Protocol (in order)

1. Commit this file and the helper scripts. The reviewer reads them; reviewer changes are dated change-log entries.
2. FIXED: the PI (or the lead, quoting the PI's delegation) replaces the status line and records G0. Commit and push.
3. The runner writes `analysis/complex_confirm/cc_run.R` (one driver; modes `gate_a`, `labelfree`, `truth`, `expjack`, `summarize`). The truth mode refuses to run unless all of these hold:
   - this file is tracked, identical to HEAD, carries `Status: FIXED`, and is pushed;
   - the member-file SHA-256 matches;
   - the 6010 cache SHA-256s are committed;
   - `--jack 500`, `--t 9` and the seeds are the preregistered ones;
   - no `verdict.csv` exists yet. **One run; no re-run after any score exists**, except after a crash or a bug shown to violate a named line of this file, in which case both runs are reported.
   - Smoke tests use **synthetic truths only**: random partitions of the net4 genes with the real complex-size distribution, never the real truth.
4. Code review by someone who has seen no score. Required changes are logged.
5. Run 6000 and 6010 (label-free). Commit the gate CSVs and cache checksums.
6. Run 6020 → 6030 → 6040 once, after G0. Commit `verdict.csv` and all tables.
7. Report: the verdict row with its fixed wording, all sensitivity lines, the disclosure table, and the E. coli history.

## Outputs and file layout

`results/complex_confirm/` (gitignored; small CSVs and markdown force-added):

- `gate_a/equivalence_gate.csv`, `run.log`
- `net4/labelfree/`: `operator_checks.csv`, `feasibility_gate.csv`, `score_manifest.tsv` (SHA-256 of the CLR, attention and a3r caches), `cache/` (score caches, not committed; a3r about 0.14 GB, as the DREAM5 caches)
- `net4/truth/`:
  - `truth_counts.csv`, `metrics.csv` (pooled AUPR and enrichment per method × truth × scope)
  - `jack_gene_replicates.csv`, `jack_complex_replicates.csv`
  - `endpoints.csv` (P and secondaries with d, SE, interval; p only for P)
  - `per_complex_ap.csv`, `strata.csv`, `verdict.csv`, `run.log`
- `net4/expjack/`: `replicates.csv`, `se_b.csv`, `refit_hash_check.csv`
- `tables.md` (all tables). A results document `docs/complex_confirm_results.md` is written after the verdict.

Design-stage files (committed with this file):

- `analysis/complex_confirm/build_yeast_complex_truth.R` (build and count only)
- `analysis/complex_confirm/power_projection.R` (existing outputs only)
- `docs/complex_confirm_selections/`: `yeast_complexes_members.tsv`, `yeast_complexes_counts.csv`, `yeast_complexes_sizes.csv`, `yeast_complexes_source_log.csv`, `yeast_complexes.sha256`, `power_anchors.csv`, `power_projection.csv`, `power_summary.csv`

**Disk** (row for `docs/disk_budget.md`):

| Item | Size | Status |
|---|---|---|
| `data/complex_confirm_raw` | 1.3 MB | counted, re-downloadable |
| `docs/complex_confirm_selections` | 0.8 MB | committed |
| `results/complex_confirm` | expected under 1 GB: three net4 score caches of 0.14 GB each plus small tables; refits keep hashes only | counted; caches deletable after the numbers are archived |

The counted total stays at about 22–27 GB of 50 GB. Nothing larger than 5 GB is downloaded. Candidate C would have needed about 14 GB and was rejected partly for that reason.

## Independence and disclosure

| Item | Status |
|---|---|
| net4 expression | Used once before, for the DREAM5 test: attention vs CLR on binding regulons (Y-M Equivalent within 10%, at the boundary) and on official edges (Y-E Inconclusive). a3r was never computed on it. Nothing was tuned on it: t = 9, keff 10, cap 50 and α = 0.5 were fixed on M3D. The attention and CLR networks on net4 are therefore known objects, and their regulon AUPRs are known. Their complex AUPRs are not. |
| Complex Portal 559292 | Never downloaded before 13 Nov in this project (the reciprocity truths were never fetched). Never scored. Its structure was counted here. It is the source named in the reciprocity preregistration (G0 item 1). |
| Truth–data dependence | Complex Portal evidence is mostly physical-interaction evidence (purification, structure). Complexes are not inferred from these expression chips [verify for a sample of entries' ECO codes; the build log records ECO codes per complex if wanted]. |
| a3r history | Chosen in-sample on K562 among about 20 readouts. Not confirmed on E. coli (26 Oct). Explored further on K562 and M3D (stacking, 12 Nov). The stacked variants are not tested here. |
| Serial testing | Fourth complex-oriented confirmatory test, second for a3r (G0 item 2). Proposed as the last for a3r. |
| Drafter's exposure | Listed at the top. No a3r, attention or CLR score on any yeast complex truth has been seen by anyone. |

## What is NOT claimed, whatever the verdict

- No claim about regulons in yeast. No untouched regulon truth exists for this dataset; the DREAM5 regulon line is descriptive on a touched truth.
- No claim against CLR on complexes. a3r − CLR is a secondary estimate without a test.
- No claim about depths other than t = 9, about a2r, a1 or stacked layers, or about per-complex AP or any sensitivity truth.
- No claim of "consistency across datasets" that tallies K562 (in-sample), RPE1 (gated, same CORUM content), the two E. coli compendia (not independent of each other) and yeast.
- No re-opening of the 26 Oct E. coli verdict, which stands as "not confirmed" whatever happens here.
- A NOT CONFIRMED or INCONCLUSIVE verdict is not evidence that a3r has no effect. A CONFIRMED verdict does not by itself make a3r the method; adoption stays with the PI.

## Limitations to state in any report

- The power is modest and uncertain: about 0.7 at +0.025 under the best-calibrated projection, as low as 0.2 under pessimistic surrogates. Effects below +0.015 are essentially undetectable.
- The yeast microarray compendium (DREAM5 net4) was the hardest network of the DREAM5 challenge. All methods were near chance on its edges (Y-E AUPR 0.019–0.022, base rate 0.017). Complex signal may still be strong, because complexes are co-expressed more tightly than TF targets, but this is unknown until scored.
- Decoys are kept in the operator as provided. They can absorb attention mass, but this is equal for both readouts.
- The ribosome is still the largest single block, at 29% of positives.
- Replicate-level chips are not independent. The gene jackknife does not resample chips; the experiment jackknife does, as a sensitivity.

## What was verified, and how (13 Nov 2026)

- **Complex Portal yeast file:** downloaded by curl from the EBI FTP (HTTP 200, Last-Modified 17 Oct 2026, 935,090 bytes, 649 complexes, SHA-256 above).
  - The format was inspected: participant column, molecule-set notation in exactly 2 complexes (CPX-1599, CPX-1601).
  - The E. coli file `data/ecoli_complexes/83333.tsv` has 0 molecule-set entries, so the expansion rule does not change any spent truth.
- **Complex Portal taxa:** human (9606) is available (downloaded to check, then deleted). *B. subtilis* 224308 returned 404.
- **UniProt map:** downloaded by REST, release 2026_04 (response header), 6,733 rows.
- **CYC2008:** two attempts, http and https, both timed out after 60 s. Not available.
- **GSE264667:** the NCBI GEO web page returned a bot check. The NCBI FTP `suppl/` listing shows only raw single-cell h5ad files (HepG2 5.2 GB, Jurkat 8.7 GB) and RAW.tar (15 GB).
- **net4 gene IDs:** 5,667 ORF-format names plus 283 `decoyN` names, no others. The DREAM5 score caches and SHA-256 sidecars for CLR and attention exist (`SCORE_MANIFEST.tsv`), as do the group-level refit hashes (`REFIT_MANIFEST.tsv`).
- **Kemmeren matrix:** 6,169 ORF rows. Matrix SHA-256 committed in `docs/reciprocity_yeast_selections/deleteome_matrix.sha256`.
- **Anchor SEs:** recomputed from the stored replicate files. They match the committed verdict: M3D 0.01817 vs `endpoints_all.csv` 0.018169; EcoMAC 0.01744 vs 0.017445; RPE1 0.01160 vs S5 0.011595.
- **Not verified:** see each [verify] marker. These cover the HepG2/Jurkat gene count and any figshare pseudobulk, SubtiWiki complexes and the *B. subtilis* compendia, and the ECO-code profile of the yeast complexes.

## Change log

- 2026-11-13: DRAFT rev 1 (design stage; no score of any readout on any candidate truth exists). Truth built and frozen (member SHA-256 `146030f7…`); power projected from existing outputs; G0 items 1–3 raised for the PI. Two rules were fixed during structure counting, before any score: molecule-set expansion and same-set paralog exclusion.
- (skeleton) YYYY-MM-DD: reviewer findings and changes.
- (skeleton) YYYY-MM-DD: **G0** PI decisions on items 1 (yeast-truth exclusion: reading a or b), 2 (serial-testing declaration withdrawn for this test or not), 3 (power accepted, or test downgraded to estimation). FIXED.
- (skeleton) YYYY-MM-DD: Ga result (max |diff| per check); Gb result (unreachable share, rows, coverage); A or fallback B.
- (skeleton) YYYY-MM-DD: DEVIATION n (before / after the first truth-stage run): what, why, who decided.
- (skeleton) YYYY-MM-DD: truth stage run: commit, design SHA-256, member SHA-256, cache SHA-256s, job IDs, peaks; verdict as written by the script.
