# Second independent review: DRAFT rev 3 of the yeast reciprocity / H1 preregistration

*1 Nov 2026. Second, independent reviewer (Claude). Reviewed: `docs/preregistration_reciprocity_yeast.md` (DRAFT rev 3, commit 1e23c71), `analysis/yeast_prep.R`, `docs/reciprocity_yeast_selections/excluded_pmids.txt`, `analysis/synthetic_width/sim_p1_rescore.R`, `summarize_p1.R` and `out/p1_rescore/` (56 run files, `summary_p1.csv`, `tstar_by_cell.csv`), and the first review `docs/plans/yeast_prereg_review.md` (fd9c827). Read-only: nothing downloaded, no truth opened, nothing scored, the preregistration not edited. The only computation was re-summarising the committed synthetic rescore files. "L" = line of the rev 3 preregistration; "prep L" = line of `yeast_prep.R`.*

**Overall: NO-GO to fix the file as is. GO once the two blockers below are fixed (both are text edits of a few lines, pre-truth) and the SHOULD-FIX (high) items are fixed or explicitly decided by the PI.** Rev 3 fixes every blocker and nearly every item of the first review in substance. The remaining problems are new, small and cheap. Most of them came from the rewrite (the verdict split, the cluster jackknife, the exclusion list) or are mechanics that the first review did not reach.

---

## Part 1. First-review findings, one by one

| Finding | Status | Where in rev 3 | Comment |
|---|---|---|---|
| **B1** truth–data circularity | **FIXED in substance; PARTIAL in mechanics** | L37, L112–118, L124–140, L146, L171–175; prep L318–365, L377 | Exclusion list E (20 GPL11232-linked PMIDs, 5 Holstege/Kemmeren PMIDs on other platforms, a mechanical reanalysis screen), committed. Primary regulon truth = binding AND expression evidence from at least one non-E PMID; binding-only fallback; complex evidence filter; contamination audit; descriptive "possibly circular" comparator. **Residual gaps:** the mode is passed as a CLI flag by whoever has seen the export, with no mechanical rule; the export → normalised TSV converter is unspecified and uncommitted; "per-association references" does not say whether the PMIDs must be attributable to the *expression* evidence. See **SF1**. The recall of E is unknown (**N3**). |
| **B2** verdict 5 reads failed NI as a cost | **FIXED**, but the split created a new gap | L350, L365–366 | 5a/5b split; "low power is not a cost" stated; structure-only n_eff(P2) computed (prep L414). **New:** the 5a/5b wording asserts a small-complex gain on the raw IUT-level P1, while L295 requires Holm for that claim. See **BL2**. |
| **B3** experiment bootstrap mis-specified | **FIXED** | L84–86, L283–288; prep L264–267 | Delete-half jackknife over label-free complete-linkage clusters (1 − r, h 0.5), 200 replicates, SE_B = SD of the replicate differences (factor (n−d)/d = 1), with no duplicates. Clusters are computed at step 4, before any truth, and committed. Residual specification gaps: **SF5**. |
| **S1** synthetic estimand ≠ P1 | **FIXED** | L23–27, L299–333 | Rescore done with the P1/S1 construction at t = 9 and own t\*. The planning effect now comes from it. The CLR > K10 round-2 pattern and the round-1 clean ×0.95 regulon cost are disclosed. A few quoted numbers are off: **SF7**. |
| **S2** gate F uncalibrated | **FIXED** (one new code risk) | L67–71, L88–108; prep L208–298 | Reworded as a sanity check; R′ = (R − 0.01)/0.99; Q ≤ 0.5 gates; median-vs-median replaces min > max, which mostly answers the asymmetric-dispersion point; the rebuild fallback is implemented; a failed F does not consume yeast. **New:** the orientation candidate set is incomplete and has no floor, so F can run on a dye-cancelled matrix. See **SF4**. |
| **S3** S1 inside a Holm family | **FIXED** | L238–243, L362–363 | S1 first at full α, then Holm over S2–S6 only if S1 is rejected (the "stop" branch); verdict 3b added. **New:** the S1 harm claim under 3a uses "Holm as above", which is undefined for S1 (**BL2**). 3b is practically unreachable (**N6**). |
| **S4** power depends on the imported ratio | **FIXED** | L337–350 | ρ = 1 and 2 rows; the 0.72 figure is labelled as assuming the M3D ratio; realised SEs are reported; MDE, not observed power; SE_B proxy added (its target is undefined, **N9**). I recomputed the power table: every cell matches (SE 0.0039/0.0051/0.0080; MDE 0.011/0.014/0.023; power 0.72/0.50/0.24 at +0.010 and 0.97/0.84/0.46 at +0.015; P2 NI power at d = 0 is 0.86). |
| **S5** divergent pairs | **FIXED** (descriptive option chosen) | L171; prep L91–98 | The strand logic is correct (C then W in coordinate order means the 5′ ends face each other). Dubious ORFs between two genes break adjacency (**N10**). |
| **S6** masking holes | **FIXED** | L78–82, L173; prep L78–90, L184–192 | Unresolvable columns dropped; SGD aliases (unique only); overlap/500 bp masking; paralog sensitivity; slow-growth rationale; growth-corrected matrix descriptive. **New:** no check that the resolved ORF is actually the deleted one (**SF6**). |
| **S7** enforcement gaps | **FIXED in text; one enforcement gap** | L32, L123, L129, L158–167, L386–404; prep L303–317, L328 | Selections, R check, test script and export log are all gated; the test-script blob hash is recorded; evaluability counts are computed; `value != 0` keeps −1; the "no R score on a real truth" line is added. **Gaps:** step 8's text contradicts the code on who performs the YEASTRACT export, and nothing enforces that the export comes after steps 1–7 (**SF2**). |
| **S8** under-specification | **FIXED** | L130–134, L181–195, L211, L219, L246–247, L279–280, L294, L154 | All nine sub-items answered. |
| **S9** verdict edge cases | **FIXED** (new gap from the fix) | L295, L359, L30 | No verdict for a failed R check; Holm over {P1, P2-NI, P2-harm} for individual claims; the Collins timing relaxation is logged as a PI decision. The Holm rule is not wired into the 5a/5b conditions (**BL2**). |
| **N1** seed collision | **NOT FIXED (moved)** | L100, L284; prep L234, L286 | The F seed is now 20261199 = 20261103 + 96, so it collides with experiment-jackknife replicate b = 96. The rebuild seed 20261105 = replicate b = 2. Harmless, but the change log calls it fixed. |
| N2 selection uncertainty | FIXED | L214–216 | |
| N3 t\* at half N | FIXED | L215 | |
| N4 S6 identifiability | FIXED | L134 | |
| N5 strain artefacts | FIXED (whole chromosome rather than arm) | L172 | |
| N6 H1 depth label | FIXED | L211 | |
| N7 family H separate | FIXED | L232 | See **N7** below on H1-harm. |
| N8 absolute margins | FIXED | L255 | |
| N9 housekeeping | FIXED | L195, L420 | The compute table sums to 150–210 core-hours, as stated. |
| N10 implementation-check threshold | FIXED (the fallback exists) | L93–98 | See **SF4** for its weakness. |

---

## Part 2. New problems

### BLOCKERS (must be fixed before the file is fixed; all are text-only)

**BL1. The DRAFT gate can never open on the fixed file.** L386 (protocol step 1) contains the literal string `STATUS: DRAFT` ("…and free of `STATUS: DRAFT`"). `yeast_prep.R` (prep L305) refuses to run if any line contains that string (`grepl(..., fixed = TRUE)`), and the test script is meant to do the same (L3). Step 1 says to replace only the DRAFT block (L3), so L386 would survive, and `truth_structure` would stop forever. Fixing it afterwards means editing the fixed file before step 8, a deviation on day one. **Fix:** reword L386, for example "and no longer carries the DRAFT status line". Grep the fixed file for the marker before committing it. A better long-term check is a positive marker: require a line matching `^> \*Fixed [0-9]{4}-[0-9]{2}-[0-9]{2}\*`.

**BL2. The verdict conditions and the claim rule disagree; two claims have no defined test.**

- L295 says individual-endpoint claims outside the IUT, explicitly "P1 superiority shown under verdict 5", use Holm over {P1-sup, P2-NI, P2-harm} at 0.025. But 5a/5b are *triggered* by "P1 passes", which L229 defines as raw p < 0.025 and d ≥ 0.005. Their verdict lines then assert "small-complex gain" (5a) and "small complexes improve" (5b).
- Reachable case: P1 p ∈ (0.0083, 0.025), P2 NI not shown, no harm. That gives verdict 5a, which states a gain that the file's own multiplicity rule does not license.
- In 5b, the harm leg may be Holm-significant while P1 is not, so "improve at a regulon cost" is half unlicensed.
- 5a's "no significant regulon loss" does not say whether it means raw or Holm. The 5a/5b partition is exhaustive only if both use the same definition.
- Verdict 3a's "large-complex cost" uses "harm test, Holm as above", but S1-harm is not in any defined family.
- **Fix (pick one, now):**
  - (a) Put P1, P2-NI and P2-harm under Holm in the 5a/5b conditions. Add **5c**: "P1 significant at the IUT level only; P2 NI not shown → not confirmed as a separate claim; estimate and interval only."
  - Or (b) declare 5a/5b statements per-test-level, as the CLR(log Pᵗ) preregistration allowed, and delete the Holm sentence for verdict 5.
- For S1-harm and P2-harm, note that the NI null (d ≤ −0.005) and the harm null (d ≥ 0) are disjoint, so at most one is true. NI and harm can then each be tested at full one-sided 0.025 without inflating FWER. Say so, and define S1-harm as "tested at 0.025 when S1-NI is not rejected (disjoint nulls)". With that observation the Holm family for verdict 5 could also shrink to {P1, P2-(NI or harm)}, but either rule is fine as long as it is written into the verdict conditions.

### SHOULD-FIX (high)

**SF1. The YEASTRACT mode decision and the export normalisation are human forks taken after the truth is visible (B1 mechanics).** L124–128; prep L26–28, L314.

- The mode is a CLI flag (`--yeastract_mode`), chosen by whoever has inspected the export. "Per-association references" is not defined: it does not say whether YEASTRACT must attribute the PMIDs to the *expression* evidence, or what happens if refs exist for some associations only.
- The export → normalised TSV step (columns `tf, target, evidence, pmids, value`) is unspecified and not committed. It classifies evidence types and assigns PMIDs to them, and those choices decide the primary truth.
- **Fix, in text now:**
  - (i) The mode is computed by `truth_structure` from the normalised file. "refs" applies iff at least 95% (fix the number) of associations that carry expression evidence have at least one PMID that YEASTRACT attributes to the expression evidence. Otherwise "binding_only". The flag becomes a cross-check that must agree.
  - (ii) If PMIDs come only as one combined list per association, not split by evidence type, the mode is "binding_only".
  - (iii) The converter is a committed script, written and committed before the export is opened, that maps YEASTRACT's evidence labels by a fixed table. Its hash goes into `yeast_truths.sha256`.
  - (iv) Record the exact YEASTRACT form fields now: documented regulations; DNA-binding OR expression evidence; activator, inhibitor and undefined; all environmental groups; any "high-/low-throughput" options left at their defaults (list them).

**SF2. Export timing contradicts the code and is not enforced.**

- L394 (step 8) says `truth_structure` "performs and logs the YEASTRACT export, then commits the log". The code (prep L307–308, L315) requires the log to be committed *before* it runs, and takes the normalised TSV as input. L123 agrees with the code.
- Nothing stops the YEASTRACT export (a truth download) from happening right after step 1, before gate F, the selections and the test script are committed.
- **Fix:** rewrite step 8 as "8a: after steps 1–7 are committed, perform the export, write and commit the log; 8b: run `truth_structure`". In code, check that the export-log commit descends from the commits of the selections, `R_check.csv` and the test script (`git merge-base --is-ancestor`). The CYC2008 file (prep L405) is also a truth obtained outside the code: give it the same rule.

**SF3. Is yeast consumed, and does H1 run, when the power gate fails?**

- L36 says a power-gate failure does not consume yeast "because no truth is opened before them". But the power gate is computed *inside* `truth_structure`, which downloads and builds every truth membership list (L393–402).
- L351 says a power-gate failure means "no score is computed", which includes H1. L232, however, says H1 is "tested whatever P1 and P2 show", and H1 does not depend on the complexes.
- If the complex arm is underpowered, the file does not say whether the regulon truth is spent or held for a later H1-only test. That is exactly the kind of choice a referee will suspect was made later.
- **Fix (PI decision):** define "consumed" as "any operator scored against a yeast truth". Then state either "H1 runs even if the P1 power gate or P1 evaluability fails (F passed and P2/H1 evaluable)" or "H1 does not run, and yeast regulons remain unused"; pick one. Correct L36 to "no truth is *scored* before them".

**SF4. The orientation check can choose a dye-cancelled rule, and the WT orientation is never verified.** L92–96; prep L142–150, L209–224.

- The candidates are "asis" and "flip_ch1" (negate the hybridisation whose channel 1 is the test). If VALUE were log2(ch1/ch2), the correct rule would be to negate the hybridisation whose channel 1 is the *reference*, which is not a candidate. Then "asis" averages a dye-swapped pair to about 0, and "flip_ch1" gives about −M (r ≈ −1).
- On failure, the code takes the rule with the highest median r with no floor (prep L220), so it would choose the dye-cancelled "asis". Gate F would then compare two near-noise matrices and could pass by chance.
- For WT, ch1/ch2 are classified by `grepl("ref", src1)`. That was validated only on mutant source names, and nothing compares the WT orientation with anything.
- The drafter's check of GSM630359 (VALUE = test/ref) makes "asis" likely correct, so this is a guard against a silent failure, not a known bug.
- **Fix:**
  - (i) Add the third candidate ("negate the hybridisation whose channel 1 is the reference").
  - (ii) Require the chosen rule to reach a median r ≥ 0.8 (fix the number), else gate F is not computable → verdict 1.
  - (iii) A structural check in `check` and `wt`: in every mutant and WT group, exactly one channel per hybridisation matches the reference pattern, and each dye-swap pair has one hybridisation of each orientation. Otherwise stop and log.

**SF5. Cluster jackknife: no rule for few or very unequal clusters, and the RNG is not pinned.** L84–86, L283–288; prep L264–267.

- With complete linkage at r ≥ 0.5, most responsive mutants will probably be singletons, which is fine. But nothing covers a small cluster count (with K ≤ 12 the delete-half SD from 200 draws is unstable), or one cluster holding a large share of mutants. In the second case, "half the clusters" varies N widely across replicates and SE_B becomes bimodal, driven by one cluster's in/out status.
- **Fix:** pre-state "if K < 20 or the largest cluster holds > 25% of mutants, then [fallback, e.g. delete-half over mutants stratified by cluster]". Commit K, the size distribution and the replicate-N range as already planned.
- Pin the RNG: record `RNGkind()` and R version; draw with `sample.int(K, floor(K/2))` over clusters sorted by ID.
- Note: complete linkage is the *least* conservative linkage for between-cluster correlation, so the shared slow-growth program is split across many "independent" clusters. That is acceptable, but disclose it in one sentence.

### SHOULD-FIX (normal)

**SF6. Deleted-gene resolution is never verified.** prep L186–190. SGD names have been reassigned since 2014, so a label can resolve to a *different* existing ORF. The wrong gene would then be masked, and the real deleted gene's value (about −5 log2) would stay in the column, which is the artefact S6 of the first review targeted. **Fix (label-free, cheap):** log the published M of the resolved ORF in its own column (before masking). Flag columns where it is not among the column's 1% lowest values, or above −1. Fix now whether flagged columns are dropped (recommended) or only reported.

**SF7. Synthetic rescore numbers: three small misstatements** (checked against `out/p1_rescore/`; details in Part 3).

- The R vs K10 own-t\* P1 median "×1.36" (L314) is the median of the 12 standard cells, while the range ×1.09–2.41 covers all 14 cells. The 14-cell median is ×1.27.
- K10 own-t\* vs CLR on P1: "0/4 to 1/4 seeds" (L316) should be 0/4 to 2/4 (`std f05 b0` has 2/4).
- "Synthetic P1 AUPR levels are 0.14–0.70" (L331) matches no file. Cell-mean levels are 0.17–0.75 for R and 0.075–0.73 for K10; K10 at own t\* is 0.075–0.62. In the batch ×2 cells the K10 own-t\* level (0.075–0.14) lies *inside* the expected real range (0.05–0.15), and there R gained +0.09 to +0.11 (×1.8–2.4). The "not transferable" argument should be restated honestly. The +0.010 planning value is still conservative.
- "R's t\* is typically 4–8 steps deeper" (L322): these are t units, not grid steps, and the per-seed range is 0–20.

**SF8. Inconsistencies between text and code to align before freezing.**

- (a) L116 describes the reanalysis screen as a title/abstract search. The code (prep L337–339) greps the whole PubMed XML record, including reference lists and data-bank links. That is broader, so harmless, but the text must match the code ("anywhere in the PubMed XML record") or a referee can call it a deviation.
- (b) L391 gates `*_selections.csv` per dataset; the code gates only `yeast_selections.csv` (prep L307). Add the HepG2, Jurkat and K562 selection files.
- (c) Gate F's "number of genes per mutant beyond the WT 0.5th–99.5th percentile" (L107) is not computed in `feasibility`.
- (d) The S6 universe (L132: pairs with exactly one regulator gene) excludes TF–TF pairs, but `edges` (prep L356) keeps TF → TF positives. Say which applies.
- (e) `truth_structure` reads the universe from `data/yeast/deleteome_responsive_700.tsv` without checking it against the committed `deleteome_matrix.sha256`, and it does not gate `masking_log.csv` or `dropped_columns.csv`. Add both.
- (f) The prep header says "DRAFT rev 2", and its comment says it never loads expression values, but `read.delim` loads the whole matrix (prep L319). Trivial.

---

## NOTES

- **N1.** The seed collisions (Part 1, N1): pick seeds outside 20261104–20261303.
- **N2. Locale-dependent order.** `sort(unique(key))` decides the "first 10 alphabetically" mutants (prep L176, L213) and the candidate order from which `sample()` draws the 376 rebuild mutants (prep L233–234). `sort` follows LC_COLLATE, so a different locale gives different mutants. Set `Sys.setlocale("LC_COLLATE", "C")` at the top of the script, and say "C-locale order" in L92.
- **N3. Recall of E is unknown.** (a) and (b) catch datasets deposited on GPL11232 or under Holstege/Kemmeren GEO series. The screen catches records that mention "deleteome" or the accessions. A paper that reanalysed the deleteome without deposit, and without naming it in its PubMed record, would stay in the primary truth. The audit measures only what E removes, not what is left. Say this in Limitations. Optionally widen the screen to "Kemmeren" in the record and to PubMed→GEO elinks hitting GSE42526/7/8/45115. The query string for list (b) is not recorded in `excluded_pmids.txt`; add it.
- **N4.** In binding_only mode, `share_comem_pos_removed` and `share_circular_removed` (prep L360, L364) compare a subset with its superset in the wrong direction and give meaningless or negative values. Label them NA in that mode.
- **N5. Depth effect (disclosed, PI decision).** The rescore shows that much of the own-t\* P1 advantage is K10 degrading at its deeper criterion-A t\*: K10 at own t\* vs CLR is ×0.39–0.96. Verdict 2 ("becomes the method's operator") could therefore rest on depth selection more than on reciprocity. This is disclosed (L34, L322–325), the same-depth contrast is pre-specified, and the PI kept the rule. No change requested. Expect it as referee point 1 if P1 passes.
- **N6. 3b is practically unreachable.** If P1 is evaluable, S1 (a superset of P1's complexes and pairs) is too. 3b can arise only from a crash. Harmless.
- **N7. H1-harm.** H1 superiority (null d ≤ 0) and the "significant CLR advantage" claim (null d ≥ 0) are both true at d = 0. Family H is therefore effectively a two-sided 0.05 test with directional conclusions. This is conventional; say it in one clause. (P2-NI/P2-harm and S1-NI/S1-harm have disjoint nulls and do not have this issue.)
- **N8. YEASTRACT vs SGD names.** prep L329 maps YEASTRACT TF and target names through the same alias resolver. Unmapped counts are logged; also log the ambiguous aliases that were dropped as duplicates (prep L82).
- **N9.** The SE_B proxy (L349) does not define the criterion-A target inside a cluster-jackknife replicate: the deleted clusters' top-1% CLR pairs? It is non-gating, but define it.
- **N10.** `divergent_pairs` (prep L91–98) and `neighbours` use all ORF-named SGD features, including dubious ORFs. A dubious ORF between two divergent genes hides the pair. This affects a descriptive sensitivity only.
- **N11.** Sample tables fetched later by `rebuild` (prep L235) are not added to `SHA256SUMS`. Call `record_hashes()` at the end of `rebuild` and `wt`.

---

## Part 3. Synthetic rescore spot-check (`analysis/synthetic_width/out/p1_rescore/`)

56 run files are present (12 standard cells × 4 seeds + 2 regulon-heavy cells × 4 seeds), and the sanity maximum relative difference is 4.4e-15 (L305 says 4e-15: ✓). From `summary_p1.csv` (14 cells):

| Draft statement | File | Verdict |
|---|---|---|
| R vs K10 t = 9, P1 ×1.03–1.64, median ×1.05, d +0.020 to +0.103, 4/4 in all 14 | 1.027–1.637, median 1.052, d 0.0198–0.1031, seeds 4–4 | ✓ |
| S1 at t = 9 ×0.99–1.01, d −0.011 to +0.004 | 0.987–1.008, d −0.0106 to 0.0039 | ✓ |
| REG at t = 9 ×1.02–1.09, 4/4 | 1.018–1.095, 4–4 | ✓ (1.0949 rounds to 1.09) |
| R vs K10 own t\*, P1 ×1.09–2.41, **median ×1.36**, d +0.050 to +0.171, 4/4 | 1.089–2.407; median **1.273 over 14 cells** (1.360 over the 12 standard cells only); d 0.0502–0.1711; 4–4 | **median mislabelled** (SF7) |
| own t\* S1 ×1.00–1.09 | 1.001–1.091 | ✓ |
| own t\* REG ×0.86–1.16; batch ×2 ×0.86–0.96, 0/4 in all four cells, d −0.001 to −0.005, levels 0.02–0.04 | 0.862–1.163; b2 cells 0.862, 0.864, 0.867, 0.959, all 0/4; d −0.0007 to −0.0049; levels 0.017–0.035 | ✓ |
| R (own t\*) vs CLR: P1 ×0.88–1.21, S1 ×1.01–1.46, REG ×0.79–1.69 | 0.885–1.208; 1.012–1.464; 0.786–1.690 | ✓ |
| K10 (own t\*) vs CLR: P1 ×0.39–0.96, **0/4 to 1/4 seeds**; S1 ×1.01–1.39; REG ×0.91–1.45 | 0.389–0.964, **seeds 0–2** (`std f05 b0`: 2/4); 1.009–1.387; 0.912–1.453 | **seed count wrong** (SF7) |
| K10 (t = 9) vs CLR, f = 10%: P1 ×0.84–1.10, S1 ×0.98–1.29, REG ×1.00–1.71 | 0.841–1.096; 0.979–1.290; 1.000–1.711 | ✓ |
| t\*: K10 10–24, R 16–32 | K10 10–24, R 16–32 (`tstar_by_cell.csv`) | ✓ |
| Synthetic P1 AUPR levels 0.14–0.70 | Cell means: CLR 0.19–0.69, K10 (t = 9) 0.17–0.73, K10 (own t\*) 0.075–0.62, R (t = 9) 0.25–0.75, R (own t\*) 0.17–0.71 | **not reproducible** (SF7) |
| Planning effect ×1.09 × 0.10 ≈ +0.010 | The bottom of the own-t\* range is `regheavy f10 b1` (×1.089); the standard-cell minimum is ×1.117 | ✓ (conservative) |

The rescore code matches its description: the P1 universe is complex-gene pairs, larger-only pairs are excluded, criterion A uses one lab split (disclosed at L307), and the sanity check is against the stored round-2 `small_x`. One difference from the preregistration's depth rule (10 random halves rather than 1 lab split) is disclosed. It matters only for planning.

---

## Ranked list

1. **BL1** (BLOCKER, trivial): the literal `STATUS: DRAFT` at L386 makes the gate refuse forever after fixing.
2. **BL2** (BLOCKER): the 5a/5b conditions use raw P1 while L295 requires Holm for the same claim; "significant loss" is undefined between raw and Holm; S1-harm has no family.
3. **SF1** (high): the YEASTRACT mode is a human flag; the converter is unspecified; PMID attribution is undefined.
4. **SF3** (high, PI decision): whether yeast is consumed, and whether H1 runs, when the P1 power gate fails; L36 is wrong.
5. **SF2** (high): step 8 contradicts the code on the export; nothing enforces that the export follows steps 1–7.
6. **SF4** (high): the orientation candidates are incomplete, there is no floor, and WT orientation is unverified, so F may run on a dye-cancelled matrix.
7. **SF5** (high): the cluster jackknife has no rule for small or dominated K; the RNG is not pinned.
8. **SF6**: verify the resolved deleted ORF from its own M value.
9. **SF7**: correct three rescore numbers and the "not transferable" sentence.
10. **SF8**: text/code alignment (screen scope, selection gating, the WT percentile diagnostic, the S6 universe, the matrix-hash check, the header).
11. **N1–N11**.

**Go/no-go:** **NO-GO to fix the file as is.** **GO** after BL1 and BL2 are fixed, SF1–SF3 are written into the text (SF3 needs a one-line PI decision), and SF4–SF5 are either specified in the text or logged as code items that must pass before step 4. All of this is pre-truth, needs no data, and amounts to about 20 lines of text plus small code edits. SF6–SF8 and the notes can go in the same pass or in the change log before step 2. Nothing found affects the integrity of what has been done so far: no truth has been opened, the exclusion list and the synthetic rescore are committed, and the quoted planning numbers stand apart from the three small misstatements.
