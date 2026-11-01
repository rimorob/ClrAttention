# Reciprocal attention: design of a held-out confirmatory test

*1 November 2026. Design memo (statistics lead, Claude Opus, high effort) for the PI. **Not a preregistration. Design only:** nothing was scored, and no truth file was opened. Sources read: `docs/plans/background_width_exploration.md` (sections 4 and 8), the two preregistrations (keff, CLR(log Pᵗ)), `docs/plans/clr_logpt_interpretation.md`, `docs/paper_plan.md`, `docs/STATUS.md`, `analysis/synthetic_width/` (code and round-2 outputs), and the committed keff verdict SEs (`docs/keff_verdict/primary_holm.csv`). Data availability was checked on disk and, for public sources, on the web on 1 Nov 2026. Marks: **[ran]** comes from committed outputs, **[derived]** is arithmetic shown here, **[assumed]** is a stated assumption.*

## 0. Two corrections to the brief

1. **The quoted gate numbers are superseded.** "Regulons ×1.19 / ×1.40, complexes ×1.08 / ×1.04" are round-1 numbers for the reproducibility gate plus reciprocity, and the ×1.08 / ×1.04 is *all* complexes, not small ones. In round 2 (section 8.5, partner-less fraction f = 10%) gate plus reciprocity **failed** the regulon condition (×0.77–0.94). The round-1 gain came from a generator in which 51% of genes had no partner. **So the gate is out, and the candidate is reciprocity alone.**
2. **The 10.1% / 0.3% statistic is about the sink, not about reciprocity.** It counts genes with no reproducible partner, which is what the gate and sink act on. Reciprocity's synthetic small-complex gain held at every f, including 5% (×1.86–6.24, 4/4 seeds, all 12 standard cells) **[ran]**. Its expected effect therefore does not depend on the partner-less fraction.

## 1. Verdict and power

**Verdict: worth one preregistered test, but only under four conditions. There is no clean held-out test with the data now on disk.**

The four conditions:

- (a) The test runs after the Collins readout. This respects the PI's 31 Oct decision: no more confirmatory complex attempts before that readout.
- (b) The data are *S. cerevisiae*: a new organism, with a compendium and truths that have never been touched here.
- (c) A power gate computed from truth structure alone (§4, step 6) passes before any score is computed.
- (d) The same preregistration also contains the headline replication on yeast regulons (§2.6), because the yeast truths can be touched only once.

Disclosures the preregistration must carry:

- **Reciprocity failed its own carry-forward rule.** It failed C2 (all complexes) because of the ribosome. Treating C2 as non-binding, since the PI targets small complexes, is a rule change made after seeing the round-2 tables.
- **This is the third complex-oriented operator test**, after keff (no change) and CLR(log Pᵗ) (not adopted). It must be declared the last test of a budget-level operator.

**Power.** Small-complex arm (superiority, one-sided α = 0.025, smallest effect of interest +0.005 AUPR):

- **[assumed] A1.** The paired gene-jackknife SE scales as 1/√n_eff. Here n_eff = (Σₖ pₖ)² / Σₖ pₖ², the effective number of complexes, where pₖ is the number of positive pairs of merged complex k.
- **[assumed] A2.** The anchor is the *E. coli* a3r complex arm: SE 0.018 at n_eff ≈ 4 (about 60% of positives sit in the r-protein blocks; the 50S and 30S blocks are taken as about 46% and 18%) **[derived]**. That gives the conservative projection SE ≈ 0.036/√n_eff.
- An optimistic projection uses the observation that a mild tweak of K10 has a smaller paired SE than a3r: on K562 CORUM, K25 vs K10 had SE 0.0044, against about 0.013 for a3r vs attention. That gives SE ≈ 0.0126/√n_eff **[ran, derived]**.
- "Combined" adds the experiment bootstrap variance at the ratio observed for CLR(log Pᵗ) (bootstrap SE ≈ 0.44 × jackknife SE), so combined SE = √(SE_J² + SE_B²) (§2.5).

| n_eff (small complexes) | SE, conservative | SE, combined | SE, optimistic | MDE at 80% power (combined) | Power at d = +0.010 (combined) | Power at d = +0.015 (combined) |
|---|---|---|---|---|---|---|
| 25 | 0.0072 | 0.0079 | 0.0025 | 0.022 | 0.25 | 0.48 |
| 50 | 0.0051 | 0.0056 | 0.0018 | 0.016 | 0.44 | 0.77 |
| **100** | **0.0036** | **0.0039** | 0.0013 | **0.011** | **0.72** | **0.97** |
| 200 | 0.0025 | 0.0028 | 0.0009 | 0.008 | 0.95 | 1.00 |

**Planning effect.** +0.010 AUPR. This is a small-complex AUPR ratio of about ×1.1–1.2 at a plausible K10 level of 0.05–0.10, heavily discounted from the smallest synthetic gain (×1.86). It must not be the synthetic ratio itself.

**Gate.** The test runs as confirmatory only if n_eff ≥ 100. That is the conservative MDE of about 0.010 AUPR, twice the smallest effect of interest. *E. coli* had n_eff ≈ 4, which is why its SE was 0.018.

**Expected n_eff.** For Complex Portal yeast (taxon 559292, about 900 KB, so about 600 complexes by analogy with the 324 complexes in the 554 KB *E. coli* file), the 3–6-gene merged complexes should give roughly 100–160. That is borderline but plausible, and it is the reason for the gate.

**Regulon non-inferiority arm (margin 0.005).** Paired SEs of mild K10 tweaks on regulons were 0.0005–0.0022 (K25 0.0008, MH 0.0005, O3-A 0.0012) **[ran]**.

- At SE 0.0015 and a true d = 0, power is 0.92 (0.86 with the combined SE).
- At SE 0.0025 it falls to 0.52.
- Synthetic regulons favour reciprocity (×1.02–1.09), except round 1's clean scenario (×0.95).
- **Joint power** of the two primary arms at d = +0.010 is about 0.6–0.7. That is adequate, not generous.

## 2. Design

### 2.1 Operator (fixed; no free parameter)

- **K10, unchanged.** Stouffer CLR, a softmax over the top 50 positive candidates at keff = 10, and a self-loop for empty rows.
- **Reciprocal attention R.** B = A ⊙ 1[Aᵀ > 0]. Set the diagonal of B to 0. Rows of B that are now empty get a self-loop. Add each row's removed mass to its diagonal, B_ii += max(0, 1 − Σⱼ B_ij).
  - In words: i keeps attention to j only if j is in i's cap-50 support and i is in j's. Unreturned mass stays on i.
  - This is `recip()` in `analysis/synthetic_width/sim_width.R` (line 259).
- **Walk and readout.** P = ½(I + B), readout S = (Pᵗ + Pᵗᵀ)/2. The attention readout only: no a3r and no cosine. a3r was not adopted, and it is fragile on top of budget changes.
- **Not tested confirmatorily.** The reproducibility gate (GR), gate plus reciprocity, and the empirical-null gate. They may appear only as descriptive results.

### 2.2 Depth

- **Primary.** Each operator at its own label-free t\*: criterion A (AP against the other half's top-1% CLR pairs) on 10 random halves of the profiles, grid {1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32}, median over splits, aggregated as in `analysis/keff_test.R`. This compares R as deployed with K10 at K10's label-free t\*, the brief's comparison.
- **Secondary.** Both operators at t\*_K10 (isolates the operator) and both at t = 9.
- **Flag.** If t\*_R > 2·t\*_K10, flag it. Unreturned mass makes the walk lazier.

### 2.3 Endpoints (pooled co-membership AUPR, R minus K10)

**Primary: adoption by an intersection–union test (Berger), each test at one-sided α = 0.025.**

| Code | Data | Truth | Hypothesis |
|---|---|---|---|
| P1 | Yeast compendium | **Small complexes**: merged Complex Portal complexes with 3–6 in-universe genes (C1's size class, unchanged) | superiority: p < 0.025 and d ≥ +0.005 |
| P2 | Yeast compendium | YEASTRACT regulons, documented with DNA-binding **and** expression evidence | non-inferiority, margin 0.005 |

**P1 truth construction.**

- Merge first, as in the CLR(log Pᵗ) preregistration: absorption, then complete-linkage clustering at Jaccard ≥ 0.5. Classify complexes by size afterwards.
- Universe: pairs of genes that are both in some merged complex.
- Positives: co-members of a 3–6-gene merged complex.
- Excluded: pairs that are co-members only of larger complexes.
- Negatives: the rest.

**Key secondary: tested only if adoption holds. Holm over S1–S6, one-sided familywise α = 0.025.**

| Code | Contrast | Hypothesis |
|---|---|---|
| S1 | Yeast, all merged complexes (every size) | non-inferiority, margin 0.005. This decides "default operator" or "complex-mode variant" |
| S2, S3 | HepG2, Jurkat: CORUM small complexes (3–6, merged) | superiority |
| S4, S5 | HepG2, Jurkat: CollecTRI | non-inferiority, margin 0.005 |
| S6 | Yeast TF–target edges (pairs with exactly one regulator; condition C4) | non-inferiority, margin 0.005 |

**Descriptive only:**

- per-complex AP;
- the ribosome alone (expected about ×0.8, the documented limitation);
- small complexes at sizes 3–10;
- yeast all-evidence regulons;
- the depth sensitivities;
- K562 essential half A (CORUM small, CORUM all, CollecTRI): a fresh operator on a used truth;
- *E. coli* M3D and EcoMAC regulons: a no-harm check on consumed truths;
- gate variants.

### 2.4 Label-free diagnostics (reported before any truth is read)

Per dataset (yeast, HepG2, Jurkat, K562, M3D):

- the share of K10 off-diagonal mass that R removes, and its distribution by row;
- the share of rows that R leaves as pure self-loops;
- the reproducible-partner count from the 5-fold k-means cluster split (PI decision 3);
- the Spearman correlation between the mass removed and that count.

These describe the mechanism. They gate nothing.

### 2.5 Resampling and decision rule

- **Decisive statistic.** z = d / √(SE_J² + SE_B²), with p = 1 − Φ(z) for superiority and p = 1 − Φ((d + 0.005)/SE) for non-inferiority.
  - SE_J comes from the paired delete-half gene jackknife: 100 replicates, one fixed seed, never rerun with another.
  - SE_B comes from the paired cluster bootstrap over yeast profiles: 200 draws, MI, CLR, K10 and R refitted in every draw, at the committed t\* values.
- **Why one combined statistic.** The interpretation memo showed that the two variance sources (which genes are in the truth, which experiments are in the compendium) add rather than compete (§3.4 of that memo). One combined statistic replaces the earlier "jackknife decisive, bootstrap necessary" double gate.
- **Human arms (S2–S5).** Jackknife only. Their bootstrap would double the compute.
- **Inference rules.** Holm and gatekeeping as above. No pooling across datasets. The verdict is computed by script into `verdict.csv`.

**One-line verdicts the test can produce:**

1. **Not run, underpowered.** n_eff < 100 at the structure check. Reciprocity stays a synthetic-only explained candidate.
2. **Confirmed, no large-complex cost.** P1, P2 and S1 hold. Reciprocal attention becomes the method's default operator from this point on. The *E. coli* headline stays K10, with the order of events stated.
3. **Confirmed, complex-mode variant.** P1 and P2 hold and S1 is not shown. It improves small complexes at a documented large-complex cost.
4. **Not confirmed.** P1 fails, whatever P2 shows. Use the CLR(log Pᵗ) wording: "not confirmed", never "no effect".
5. **Regulon trade-off.** P1 holds, P2 fails, or there is a significant loss of at least 0.005.
6. **Human replication** (only after verdict 2 or 3): S2 and S3 are reported as "replicated in HepG2 / Jurkat" or "not replicated".

### 2.6 Headline replication in the same preregistration (separate family H, its own α = 0.025)

**H1.** Yeast, YEASTRACT regulons (P2's truth): K10 at t\*_K10 minus Stouffer CLR, superiority with d ≥ +0.005.

- **Why it belongs here.** CLR is scored anyway as a comparator, so leaving H1 out would compute the number without a test.
- **What it would add.** A second organism for the paper's main claim. The draft notes already flag that one bacterial compendium carries the result.
- **What it risks.** A yeast failure would weaken the generality of the headline. That risk is the price of a held-out test (§5).

## 3. Datasets and truths: availability check

**On disk (`data/`, `results/`): nothing can serve as confirmatory.**

| Path | What it is | Status |
|---|---|---|
| `data/E_coli_v4_Build_6/`, `data/ecomac/`, `data/RegulonDBExtract/`, `data/ecoli_complexes/` | M3D, EcoMAC, RegulonDB, Complex Portal 83333 | Consumed by the headline, keff and CLR(log Pᵗ); no-harm descriptive only |
| `data/perturbseq/K562_essential_*.h5ad`, `data/human_truth/{complexes,collectri,dorothea}.tsv` | K562 essential and its truths | In-sample: about 20 operators scored. `d_geomean_t9`, a mutual-attention readout √(PᵢⱼPⱼᵢ), was scored on K562 CORUM in-sample (`analysis/corum_exploratory/metrics_exploratory.csv`; the drafter saw this row) |
| `data/perturbseq/K562_gwps_normalized_bulk_01.h5ad` | K562 genome-wide | Same cells, same truths |
| `data/perturbseq/rpe1_normalized_bulk_01.h5ad` | RPE1 | CORUM already scored (CLR(log Pᵗ) S5); same study as K562 |
| `data/beeline/`, `data/chipatlas/`, `data/string/` | BEELINE single cell | CLR-NS was chosen on these truths; no complex truth; attention adds nothing there |
| `data/CLR versions/CLR.zip` → `CLR/InfoKit2/ClusterScript2/yeast.csv.gz` | A legacy 6,194 × 515 matrix | **Unusable:** no gene or sample IDs, unknown provenance |

Reciprocity (S2) has never been computed on any real dataset (only `analysis/synthetic_width/` contains it).

**Public sources (checked 1 Nov 2026), all to be downloaded:**

- **Yeast expression, primary choice: the Kemmeren et al. 2014 deleteome.**
  - GEO GSE42527 (700 responsive mutants) and GSE42526 (784 non-responsive mutants): 1,484 deletion profiles, two-channel arrays against a common wild-type reference, about 6,100 genes, one lab and one platform.
  - Why it fits: it is a perturbation compendium like K562, the setting where the complex loss was seen, and it has a strong slow-growth global program, like the "messy" generator.
  - Choice rule: the compendium is chosen on design grounds only, before any truth is read. A condition compendium (several labs, as M3D) is the alternative (open decision 3).
- **Yeast complexes.** EBI Complex Portal `complextab/559292.tsv` (913 KB, release 2026-10-17). It uses the same builder as *E. coli* (`analysis/build_ecoli_complex_truth.R`), with ORF mapping from UniProt 559292.
  - Fallback rule, fixed now: if Complex Portal gives n_eff < 100, take its union with CYC2008 (if obtainable) **before any scoring**. If n_eff is still < 100, the result is verdict 1.
- **Yeast regulons.** YEASTRACT+ "Generate Regulation Matrix" (documented regulations, filtered by evidence type). The bulk export format is unverified. Log the export date and parameters.
- **Human (key secondary).** Nadig et al. 2025 HepG2 and Jurkat essential Perturb-seq. *Correction (1 Nov, later the same day):* GSE220095, the accession a secondary source gave, is an unrelated prostate-cancer study. The correct series is **GSE264667** (PMID 40259084). It distributes only single-cell files, `GSE264667_hepg2_raw_singlecell_01.h5ad` (5.6 GB) and `GSE264667_jurkat_raw_singlecell_01.h5ad`, plus a 16 GB RAW tar, so the pseudobulk must be built. Truths: the existing CORUM and CollecTRI files, unchanged.

**What must be built:**

1. `analysis/yeast_deleteome.R`: download and parse the data, average replicates and dye swaps, and set each deleted gene's own entry to its median over profiles (a label-free rule, fixed now). Gene filter: complete in ≥ 95% of profiles. Systematic ORF IDs.
2. `analysis/build_yeast_truths.R`: complexes (flatten, 3–200 genes, merge, size class) and regulons (DNA-binding and expression evidence), plus the structure counts and n_eff.
3. A port of `recip()` into `analysis/reciprocity_test.R`, with a unit test that reproduces `sim_round2.R` RECIP metrics on one synthetic seed to 1e-10.
4. A HepG2 and Jurkat pseudobulk via `analysis/perturbseq_bulk.R`, if the data are verified.

## 4. What must be committed before any truth is read

1. **`docs/preregistration_reciprocity.md`** (this design, with α, margins, seeds and the outcome table), plus the git gate: the script refuses to read a truth unless the file is identical to HEAD, and it writes that commit's hash into every output.
2. **Data provenance**: GEO accessions, URLs, download dates and SHA-256 of the raw files; SHA-256 of the processed matrix; the preprocessing script's hash.
3. **The implementation check** (synthetic reproduction) passed and committed. It reads no real truth.
4. **The label-free phase**:
   - split seeds and memberships;
   - t\*_K10 and t\*_R for each dataset, in `docs/reciprocity_selections/`, byte-identical gate;
   - the §2.4 diagnostics;
   - bootstrap draw seeds.
5. **The truth builds**: raw truth files and their SHA-256, release dates, YEASTRACT export parameters; built member lists and their SHA-256.
6. **The power gate**: structure-only counts (complexes, positives, n_eff for P1 and S1, regulon counts), the projected SE from §1, and a pass or fail line, committed before scoring. If it fails, stop and report verdict 1.
7. **Code review** by a reviewer who has seen no yeast or HepG2 score. Then one truth phase: yeast, then human, then the descriptive sets. The first complete run is the confirmatory one.

## 5. What a positive and a negative result would mean

**Positive (verdict 2 or 3, and H1 passes):**

- It is the paper's first *confirmed* operator change. A mechanism (unreturned attention from junk rows competes with small-complex pairs) was found on synthetic data, fixed without labels, and confirmed in a new organism on truths that had never been touched.
- For a Nature-family reader, the arc becomes: a classic method recast as attention, its failure mode explained by mass conservation, a parameter-free fix predicted, and that fix confirmed out of sample. It also maps onto a known transformer idea (mutual-neighbour attention and hubness reduction), which strengthens the "CLR is attention" framing.
- H1 adds a second organism to the headline. **That is probably worth more for impact than reciprocity itself.**
- Caveats:
  - Verdict 3 limits the claim to small complexes; the ribosome limitation stays.
  - One organism plus gated human arms is "confirmed in yeast", not "general".

**Negative (verdict 4):**

- The third complex-oriented operator test is not adopted. The complex section becomes "the loss is explained, and the budget-level fixes tested (width, log readout, reciprocity) were not confirmed".
- That is honest and keeps the headline intact, but cumulatively it reads as an unsolved limitation, not a feature.
- The synthetic reciprocity result may still be reported, labelled as not confirmed.

**Headline failure (H1 fails):** attention does not beat CLR on yeast regulons. That would be a material hit to generality, and the paper would have to report it. Running yeast exposes the headline. The PI should choose this knowingly; it is also exactly what a Nature-family reviewer would ask for.

**Verdict 1 (not run):** nothing changes. The paper says that reciprocity is a synthetic-only candidate because no held-out truth with enough small complexes was available.

## 6. Compute estimate (for `jobs_pc/sched`)

Scaling from M3D, where one MI + CLR refit takes about 3.6 core-minutes (200 refits took 3 h at 4 workers); cost scales with genes² × samples **[derived]**.

| Block | Core-hours | Peak RAM | Disk |
|---|---|---|---|
| Yeast deleteome (about 6,100 × 1,484): full fit 0.4; t\* (20 half-fits plus walks on the t grid for 2 operators) 15–25; jackknife (100 replicates × 4 truths × 2 operators, plus CLR) 2–4; **cluster bootstrap, 200 refits plus walks, 80–100** | **100–130** | about 4 GB per worker (dense 6,100² is 0.3 GB) | raw GEO 1–3 GB |
| HepG2 + Jurkat (about 8,700 genes × about 2,400 profiles each): full fit about 1.3 each; t\* about 15 each; walks (8,700³ multiplications) about 20 each; jackknife about 3 each | **70–90** | 10–12 GB per worker (as RPE1) | 20–40 GB if built from single cells, plus a 1–2 h build at about 60 GB |
| K562 descriptive (R only; K10 cached) and label-free diagnostics on M3D and K562 | 10–15 | 6 GB | none |
| Downloads, truth builds, implementation check | 2–5 | under 4 GB | under 1 GB |
| **Total** | **about 180–240** | budget 8 workers × 6 GB (yeast) or 4 × 12 GB (human) | **about 45 GB at worst** (323 GB free) |

At 40 threads (70–80% of 56), the wall time is about 6–8 h of compute, after 1–2 days of building, provenance work and review.

## 7. Open decisions for the PI

1. **Go or no-go, and what the single yeast shot contains.** Running the test means accepting three things: C2 is relaxed after the fact (disclosed), the test is the last of its kind, and it runs after the Collins readout. The yeast truths can be touched once, so the same preregistration must also decide whether it includes H1 (recommended), a3r's last test (the interpretation memo's §7.3 asked for yeast), and the C3 asymmetry edge test. Whatever is left out cannot use yeast later.
2. **The human arms.** Keep HepG2 and Jurkat as gated key secondary (recommended: they reuse CORUM, and an extra primary arm costs about 20% of joint power), or promote HepG2 to co-primary for a stronger human readout. PI decision (1 Nov): secondary. The data are GSE264667 (see the correction in §3).
3. **The yeast compendium and the power gate.** The deleteome (perturbation design, single lab) or a condition compendium (closer to M3D). The n_eff ≥ 100 threshold (MDE about 0.010) and the CYC2008 fallback rule. And whether to adopt the combined jackknife + bootstrap SE as decisive (recommended) instead of the house "jackknife decisive, bootstrap necessary" rule.
