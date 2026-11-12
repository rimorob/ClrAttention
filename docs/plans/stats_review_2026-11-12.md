# Independent statistical review of the Collins pre-read (12 Nov 2026)

Reviewer: independent statistical review (Claude Opus 5.5), no prior involvement. Scope: `docs/checkpoint_collins.md` as of commit a2800ac, mainly item 8 (DREAM5), "Exploratory follow-ups" A and B, and the earlier sections where items 8, A and B change what they can say. Read-only. No file other than this one was written; no job, compute or credential was touched. Scratch scripts: `/tmp/claude-1000/-home-boris-Projects-ClrAttention/4451c478-b0bc-4b2f-84d1-11ecf74c2d08/scratchpad/{d5.R,slope.R}` (they read only the committed CSVs).

Sources checked: `docs/plans/dream5_design.md` (with its change log), `docs/dream5_selections/truth_stage_log.md`, `dream5_report.md`, `docs/dream5_truth_stage_checklist.md`, `results/dream5/truth_net{1,3,4}/*.csv`, `docs/preregistration_depth_log.md`, `docs/depth_log_results.md`, `results/granularity_log/*.csv`, `docs/robust_z_results.md`, `results/robust_z/{m3d,k562}/*`, `docs/whitepapers/sota_table.md`, `docs/paper_plan.md`, `docs/STATUS.md`, `docs/keff_verdict/*`, and the git history of `analysis/robust_z/`.

---

## Verdicts per question

### Q1. DREAM5 net4 (item 8)

**Interval construction: defensible as preregistered, with one anti-conservative part.** I checked the replicate files: SD of the 200 truth-jackknife d values is 0.02255 and SD of the 100 group-level refit d values is 0.01731. Both equal `SE_J` and `SE_B` in `verdict.csv`. The replicate distributions are close to symmetric (skew between −0.09 and +0.23), so the normal interval is reasonable. Delete-half resampling with variance factor 1 is correct for smooth statistics. Adding the two variances double-counts the gene × experiment interaction, which is conservative (Owen 2007). The weak point is the experiment unit. Chips within a `#Experiment` series (time courses, perturbation series) are correlated, and the series-level SD is 2.3 times SE_B for Y-M (0.0402 against 0.0173) and 1.4 times for Y-E. The design states this as a limitation and made the group-level SE_B the only one that decides, so this is not a deviation. A referee will still compute the series-level interval, and the pre-read should give it (finding M1). Half-data refits also show a smaller effect (mean refit d 0.067 against the full-data 0.082). That is expected and does not change the SE.

**"Equivalent within 10%, only just": the label is not robust and should not be read as a finding.** The upper limit 0.137416 is 8.8e-5 below delta = 0.137504. Each SD comes from a finite number of replicates (200 and 100), so the interval endpoint carries Monte Carlo error. Resampling the replicate values gives an MC SD of the upper limit of **0.0022**, about 25 times the gap. In **44%** of those resamples the upper limit reaches or passes delta. Three reasonable alternative analyses would each have produced "Inconclusive":
- series-level SE_B instead of group-level: d 95% interval **−0.009 to 0.172**, which includes 0 and passes +delta;
- Bonferroni across the two co-primaries (0.0125 per side, 97.5% interval): **0.018 to 0.145**, which passes +delta;
- the MC uncertainty alone, as above (a coin flip).

One sensitivity moves the other way. With the preregistered expected-tie AUPR, d = 0.072 and the interval is 0.017 to 0.128, inside the margin. CLR's tied zeros are counted at the end of their block in the primary `ap_fast`, which favours attention slightly. Under the conventional 90% TOST the interval would be 0.035 to 0.129. The design's 95% (0.025 per side) is stricter than usual, which is one reason the result sits on the boundary.

**What is robust:** the point estimate, 5.8% relative gain, is below the 10% margin. Superior was therefore out of reach whatever the SE, and the preregistered conclusion "the co-membership headline did not replicate in yeast" does not depend on the knife-edge. The superiority rule (point estimate at least delta, and p < 0.025 against zero) is lenient, and the gain still fell short. That strengthens the negative. Proposed wording in M1.

**Multiplicity.** Two co-primaries, each tested one-sided at 0.025 per direction, with no combined claim and an uncontrolled familywise error that the design states, is acceptable. The combined interpretation lines ("module gain with edge cost") are intersection-union statements and need no adjustment. Multiplicity changes no conclusion here, because neither endpoint produced a positive claim. It does matter for the "Equivalent" label, which disappears under Bonferroni (above). The pre-read should not let "Equivalent" stand as a finding.

**Prediction mismatch (Y-E).** The facts are reported, but the sentence "The estimate has the opposite sign from the prediction (not positive)" is garbled: it reads as if the estimate were not positive. Proposed wording in S1. Fairly stated, the prediction d ≤ 0 was neither confirmed nor refuted (interval −0.093 to 0.213). The point estimate is +0.060.

**Deviations.** I found none in the primary analysis. The guard, scorer gate (36 of 36 methods), step-6 code review, seeds and switches are all as logged. Three things the pre-read could state more precisely:
- "Design frozen 4 Nov" is true for methods, endpoints, margin and seeds. The definition of the Y-E scorer (trapezoid AUPR, `drop_tied`, the gold-universe box, lookup p-values) was fixed on 6 Nov, after the scoring archive, including the gold standards, was extracted, and before any score of ours. The design allowed this and it is logged. Say it once.
- The PIDC tie at the cap was handled by a pre-specified rule. It is not a deviation.
- The design's exploratory list (replicate-averaged input, tree-method module readouts, hub-stratified scores, t = 3 and 27, CLR-NS attention, per-TF AUROC) has not been run (`dream5_report.md` disclosure item 16). If any of it is run now, it is run after unblinding and must be labelled that way.

### Q2. Depth re-test (follow-up A)

**The reading rule is not statistically meaningful at n = 3 and n = 5. "Not resolution-like" overstates what the data show.**
- With n = 3, rho ≤ 0.5 holds for 5 of the 6 possible orderings, so reading (b) fires in E. coli unless the order is perfect. Under random orderings, reading (b) occurs with probability about 0.97 overall (1 − (1/6) × P(K562 rho ≥ 0.6)).
- In the module bootstrap, E. coli reached rho ≥ 0.9 in only 10% of draws (`depth_log_boot_rho.csv`). The design could essentially return only (b) or (c).
- The order in E. coli is broken by one truth: sigmulons, which has **5 sets** and a bootstrap interval of the best depth of **0 to 256**, i.e. undetermined.
- The AUPR-against-depth curves are flat. Within-1% intervals span up to 3.4 octaves, and K562 gains are 0.005 to 0.017 AUPR. An argmax over a flat curve is a noisy statistic.
- K562 rho is 0.70 with CLR-best ranked lowest. It is **0.20** over attention depths only (`depth_log_primary.csv`). The preregistration says that convention was chosen after seeing the linear grid. The pre-read does not mention the 0.20.
- The coarseness rank is nominal. K562 KD-response modules (median 100 genes) are ranked second, finer than DoRothEA (median 13).

The literal reading must be reported. "Not resolution-like" should be glossed as "the pre-specified ordering was not observed; with 3 and 5 truths the test could neither establish nor exclude it" (M4).

**Pooled slope with the module bootstrap: not valid as evidence that depth scales with size.**
- The module bootstrap holds the 7 truths and the operators fixed. It measures only within-truth module-sampling noise and ignores between-truth scatter, which is the variance that matters for a slope fitted across 7 truths.
- An ordinary least-squares fit on the same 7 points gives slope **0.197, 95% CI −0.77 to 1.16 (P = 0.62)**. With a separate intercept per organism it is 0.23 (−0.96 to 1.42).
- Leave-one-truth-out slopes run from **−0.22 (dropping operons) to 0.35**. The positive slope rests on the operon truth.
- Truths overlap: DoRothEA A–B is nested in A–C, CollecTRI overlaps both, and E. coli operons share genes with regulons. Modules are resampled independently across truths, which ignores this.
- 48 of 1,000 E. coli draws were dropped when sigmulons had best t = 0, so that interval is conditional.
- E. coli and K562 operators differ in spectrum and mixing time, so pooling them with one intercept is hard to justify.

Recommended: report the slope as descriptive, give the between-truth interval, and drop any statement that depth "does scale" with size (M3).

### Q3. Robust row-standardisation (follow-up B)

**"No variant is consistently better" is supported.** It is the right reading of the tables. Variants, truths and read-outs were fixed in code before any result: `analysis/robust_z/*` was committed in ed64bfa at 17:31 on 11 Nov, the results in 06e7ef9 at 21:37, and the code did not change in between. The forking-paths risk is in **presentation**, not in analysis:
- `tables.md` holds about **735 jackknife intervals** (170 for M3D, 565 for K562). The pre-read shows 56 contrasts and highlights the most favourable pairwise cell (V1 on CORUM all, +0.547).
- The text repeatedly uses "interval excluding 0" as a criterion. That is a 5% test in all but name, applied about 700 times, and it contradicts "none should be read into the intervals". State the count and drop "excluding 0" as a criterion (S8).
- In M3D the summary says every robust variant is "at or above V0" in the pairwise score. That is true of pooled enrichment, but mean per-regulon AP is lower for every robust variant (−0.008 to −0.014). The summary should carry both (S9).

**"Robust variants drift toward hubs by pi" is weaker than stated, and one summary sentence is factually wrong.**
- (i) pi is the stationary distribution of the **V0** walk. For V0's attention read-out the pi share is partly tautological, because the walk converges to pi. For the other variants it measures overlap with V0's hubs, not hubness as such.
- (ii) For attention at fixed t = 9, the comparison mixes hub preference with mixing speed. V1, V2 and V4 have more positive scores (0.62 to 0.64 against 0.53 to 0.54), so their operators are denser, and their label-free depths are 5 to 7. At t = 9 their walks are relatively further along. **At each variant's own t*, the K562 attention pi shares are 0.571 (V1), 0.502 (V2), 0.382 (V3) and 0.559 (V4), all below V0's 0.608** (`hub_drift.csv`, kind `attT`). In K562 the direction for attention reverses. In M3D at t*, V1, V2 and V4 stay above V0 and V3 falls below.
- (iii) The MI row mean ranks genes independently of any variant, but it is **not neutral** with respect to the treatment. The variants differ precisely in how they remove the row location (mean, median, or two-way median polish). Median centring penalises high-mean rows less than mean centring, and median polish removes row effects outright. That is why V1 and V4 go up and V3 goes down under this metric. The disagreement between stratifiers (+0.014 under pi, +0.483 under the MI row mean) is therefore expected from the construction of both, and neither is "the" hub-controlled answer.
- (iv) The drift shares have no intervals, so differences such as 0.169 against 0.175 cannot be interpreted.
- (v) **Mismatch:** the summary says "Only V3 (rCLR) moves the hub share down, and only by the MI row mean". The table above it shows V3 also lower by pi in K562 (pairwise 0.328 against 0.336; attention 0.543 against 0.608), and V2 lower by MI row mean in M3D pairwise (0.130 against 0.169), M3D attention (0.112 against 0.116) and K562 attention (0.331 against 0.413). See M5.

**Jackknife intervals on pooled enrichment: acceptable as labelled (descriptive, gene sampling only), with three caveats.**
- The networks are held fixed. Item 8 shows that experiment-level variance is of the same order as gene-level variance or larger (net4 Y-M SE_B 0.017 against SE_J 0.023; series-level 0.040).
- Pooled CORUM-all is dominated by 8 RP complexes. With 100 replicates and normal ±1.96 SD intervals, a heavy-tailed replicate distribution may be poorly covered (MC error of an SD from 100 replicates is about 7%).
- Delete-half over genes keeps about a quarter of the pairs. Using the variance factor 1 is fine for these near-smooth statistics.

**Does the text avoid implying unearned results?** Mostly yes. The labelling (exploratory, not for the paper, not held out) is exemplary. The exceptions are the hub-drift summary ((ii) and (v)) and the M3D "at or above" sentence.

### Q4. Spot checks and sentence-level overstatement

Thirteen numbers were checked (table below): all match the CSVs except one summary sentence (M5). Sentence-level fixes are in the findings table.

### Q5. The project as a whole

**Strongest claim the evidence supports.** On *E. coli* expression compendia (M3D and the partly overlapping EcoMAC), a parameter-free completion of CLR as per-gene softmax attention plus a short walk improves pooled recovery of strong-evidence regulon co-membership over CLR (+0.022 AUPR; enrichment 1.68 against 1.44, stable over nested, laboratory-split and replicate-aware resampling) and ranks first in pooled AUPR against GENIE3, PIDC and GRNBoost2. Per regulon it ties, and it loses on protein complexes and on direct E. coli edges. In the one independent preregistered test (DREAM5 yeast) the co-membership gain was small, 5.8% (95% CI 1.8% to 10.0% with group-level resampling; it includes zero with series-level resampling), below the 10% margin, with near-chance absolute performance. That is a well-characterised method contribution with an E. coli-specific gain, not "better recovery of regulatory programs from expression alone" in general. Decision 3's proposed central claim should be scoped to E. coli, or framed as a reframing of CLR with characterised trade-offs.

**Main vulnerability a referee would press.** Every positive result comes from the development organism and data lineage: M3D, EcoMAC (37.7% shared arrays), DREAM5 net3 (same lineage), all with the RegulonDB truth.
- The nested bootstrap protects the depth choice. It does not protect the other settings (keff 10, cap 50, alpha 0.5, Stouffer, the MI estimator, truth construction with operon exclusion), which were settled during development on the same data.
- The independent test did not replicate at the stated margin.
- In the in silico network, where the truth is exact, the E. coli pattern reverses (edge ratio 1.15, jackknife interval 0.025 to 0.381; module ratio 0.98).
- Human results are null on TF truths and negative on complexes and knockdown modules.
- The gain is pooled, concentrated in a few global regulons (62% of positive pairs; `paper_plan.md` section 5), and hub-related.

A referee will say the gain is a property of pooled AUPR on hub-centred E. coli regulons in one compendium lineage.

**Many looks.** Yes, this should be disclosed. Preregistered tests to date:
- keff test: 32 Holm-corrected tests; no variant met the decision rule; K25 and O3A each won 2 of their 4 tests;
- CLR(log P^t): **not adopted**, failure to confirm (26 Oct). This is **absent from the pre-read**;
- yeast reciprocity: not run, the gate failed;
- scPRINT: superior, against a comparator at permutation level;
- DREAM5: not replicated at the margin;
- depth re-test: reading (b).

On top of these come dozens of exploratory analyses on the same M3D/RegulonDB and K562 truths (FDR operators, value transforms, multiscale, width rounds, rCLR, robust z, macro, hub, granularity). The pre-read is commendably candid about negatives, but it should carry a one-table ledger of every confirmatory test and its outcome, and add the CLR(log P^t) result (S6).

---

## Findings

Severity: **must-fix** = a reader would be misled, or a number or statement is wrong; **should-fix** = a hostile referee would seize on it; **note** = worth knowing, optional.

| ID | Severity | Location (pre-read) | Issue | Proposed fix (replacement wording) |
|---|---|---|---|---|
| M1 | must-fix | Item 8, Y-M bullet and the series-level bullet | "Equivalent within 10%, and only just" invites reading Equivalent as a finding. The gap to the margin (8.8e-5) is 1/25 of the Monte Carlo SD of the endpoint (0.0022; crossing in 44% of replicate resamples). Series-level resampling gives −0.009 to 0.172 and Bonferroni across the co-primaries gives 0.018 to 0.145, both Inconclusive. The series-level interval is not given. | "**Y-M (regulon co-membership):** attention 0.111 vs CLR 0.105, a 5.8% relative gain (d = 0.082, 95% interval 0.026 to 0.137, i.e. +1.8% to +10.0%). Preregistered verdict, as written by the script: Equivalent within 10%. The upper limit is 0.0001 below the margin, which is within the Monte Carlo error of the interval (about ±0.002), so the boundary between Equivalent and Inconclusive is not resolved; with the coarser series-level resampling the interval is −0.009 to 0.172 (Inconclusive, includes 0). What is robust: the estimate is below the 10% margin, so the result is not Superior under any of these intervals. Fixed wording: the co-membership headline did not replicate in yeast." Delete "the direction is positive" or replace it with "the estimate is positive; the interval excludes zero only under group-level resampling". |
| M2 | must-fix | Item 8, networks 3 and 1 bullet, last sentence | "attention has the best yeast AUROC significance of our seven methods (−log10 p 5.6)" picks the one favourable cell of 42 (AUROC 0.523 against 0.522 for PIDC). It contradicts `dream5_report.md` section 7, which lists "any statement that net4 supports a ranking of methods" as not supported. | "Official challenge score (descriptive; our reimplementations): GRNBoost2 and GENIE3 lead overall (55.6 and 43.8), attention is third (38.2) and CLR fourth (22.9); on yeast all seven methods lie within AUROC 0.510 to 0.523, so no ranking is supported there." |
| M3 | must-fix | "What did not work", last bullet; follow-up A "Best depth against module size" bullet | "The best depth does scale with module size, roughly as size^0.4 … (pooled slope 0.197 [0.078, 0.668]) … generic to diffusion." The module-bootstrap interval ignores between-truth variance. Least squares over the same 7 truths gives 0.197 (−0.77 to 1.16), P = 0.62, and dropping operons flips the sign (−0.22). The 0.4 is a different statistic (per set, within regulons, linear grid; `paper.md` line 144). "Generic to diffusion" is not tested here. | Bullet: "… and the preregistered re-test on a geometric grid to t = 256 did not find the fine-to-coarse ordering either (pre-specified reading (b); with 3 and 5 truths it could not have established or excluded it, follow-up A). Across truths, the best depth is not clearly related to median module size (slope 0.20; module-bootstrap interval 0.08 to 0.67, which ignores between-truth variation; across the 7 truths −0.77 to 1.16; the sign depends on the operon truth). Within regulons, an earlier per-set analysis on the linear grid found best depth rising with set size (t ∝ n^0.43) but not for random sets." In A, add after the slope numbers: "The interval resamples modules within fixed truths and operators; it does not reflect variation between truths (least squares across the 7 truths: −0.77 to 1.16)." |
| M4 | must-fix | Update paragraph; A "Pre-specified reading"; "What did not work" | "Not resolution-like" reads as a demonstrated negative. Under random orderings reading (b) fires with probability about 0.97. E. coli needed a perfect order of three truths, which occurred in 10% of bootstrap draws, and the breaking truth (sigmulons, 5 sets) has a best-depth interval of 0 to 256. | Update paragraph: "(A) The log-scale depth re-test, preregistered before any compute: the pre-specified reading is (b) 'not resolution-like', i.e. the fine-to-coarse ordering was not observed; with 3 and 5 truths and flat depth curves the test could not establish or exclude it." In A, after the reading, add: "The rule fires in E. coli unless all three truths are perfectly ordered; the order is broken by the sigmulons (5 sets; bootstrap interval of best t 0 to 256). In K562, rho is 0.70 with CLR ranked lowest and 0.20 over attention depths only." |
| M5 | must-fix | B, "Summary, as far as the tables go" | Factual mismatch with the hub-drift table directly above it. V3 is below V0 by pi in K562 (0.328 against 0.336; 0.543 against 0.608), and V2 is below V0 by MI row mean in M3D pairwise (0.130 against 0.169) and in K562 attention (0.331 against 0.413). | Replace "Only V3 (rCLR) moves the hub share down, and only by the MI row mean." with "V3 (rCLR) lowers the hub share by the MI row mean in both datasets and by pi in K562; V2 lowers it by the MI row mean in M3D and in K562 attention." |
| S1 | should-fix | Item 8, Y-E bullet | "The estimate has the opposite sign from the prediction (not positive)" is garbled. | "Inconclusive. The design predicted no gain here (d ≤ 0, the M3D pattern); the estimate is a 4% gain (d = +0.060), but the interval runs from a 6% loss to a 16% gain, so the prediction is neither confirmed nor refuted." |
| S2 | should-fix | Header (line 3) and item 8, last bullet | The header says key results were independently reviewed. Item 8 ends "Not yet independently reviewed", which is ambiguous about whether it covers the whole item. | Header: "… The key results in items 1 to 7 have been checked by an independent statistical review; item 8 (DREAM5, run 11 Nov) and the exploratory follow-ups have had a separate statistical review (12 Nov, `docs/plans/stats_review_2026-11-12.md`)." If it is decided not to cite this review, say "item 8 has not yet been independently reviewed". |
| S3 | should-fix | "In one paragraph"; item 4 "Direct TF→target edges"; Decision 3 | The mechanism "the walk trades edge precision for module recovery" is now contradicted outside the E. coli lineage. The in silico net1 (exact truth) shows edge ratio 1.15 (jackknife interval excludes 0) and module ratio 0.98. Yeast net4 points up on both readouts. The summary paragraph also omits the yeast test. | Item 4: "**Direct TF→target edges:** on E. coli (M3D; DREAM5 network 3, same lineage) attention trails CLR. The trade is not seen in DREAM5 yeast (both readouts point up) or in the in silico network (edge gain 1.15, module 0.98, descriptive)." Add to "In one paragraph": "In the one independent, preregistered test (DREAM5 yeast), the co-membership gain was 5.8% (95% interval 1.8% to 10.0%), below the pre-specified 10% margin, and every method was near chance." Decision 3: scope the claim to E. coli, or to "characterised trade-offs", and say so. |
| S4 | should-fix | Item 5, first bullet: "Attention's small DoRothEA lead over CLR survives hub control" | It survives under the attention-derived stratifier (+0.031 [+0.009, +0.053]). Under follow-up B's MI-row-mean stratifier it is +0.016 [−0.010, +0.042] (`results/robust_z/k562/tables.md`, attention t = 9 minus CLR, MI row mean). | "… survives the attention-derived hub control (it ties GENIE3 and PIDC); under the MI row-mean stratifier of follow-up B (exploratory) the lead is +0.016 [−0.010, +0.042]." |
| S5 | should-fix (positive) | Item 3, hub control | Follow-up B gives a check of the M3D strong-evidence hub claim that does not use the attention walk: attention minus CLR under the MI-row-mean stratifier is +0.182 [+0.124, +0.240] (all evidence +0.031 [+0.002, +0.059]). | Add: "Under a stratifier that does not use the walk (row mean of raw MI; follow-up B, exploratory) attention stays ahead of CLR on M3D strong evidence (+0.18 enrichment [+0.12, +0.24])." This covers attention against CLR only, not GENIE3 or PIDC. |
| S6 | should-fix | "What did not work"; new short ledger | The confirmatory CLR(log P^t) test, NOT ADOPTED on 26 Oct (failure to confirm; P1 p = 0.079 and 0.051), is missing. There is no tally of confirmatory looks. | Add a bullet: "CLR applied to log P^t as a complex fix: preregistered, failed to confirm (26 Oct; EcoCyc/Complex Portal complex superiority p = 0.079 M3D, 0.051 EcoMAC; regulon non-inferiority passed); attention stays the method." Add a 6-row ledger: keff (no variant met the rule), CLR(log P^t) (not adopted), yeast reciprocity (not run, gate failed), scPRINT (superior; weak comparator), DREAM5 (not replicated at the margin), depth log (reading b). |
| S7 | should-fix | B, hub-drift paragraph and summary "The robust variants V1, V2 and V4 drift toward hubs by pi" | For attention the drift is measured at fixed t = 9. The variants mix faster (t* 5 to 7, denser operators). At their own t*, all four variants have a **lower** K562 attention pi share than V0 (0.38 to 0.57 against 0.61). pi is V0's own stationary distribution. The MI row mean is aligned with what the variants change (row centring). No intervals. | Summary: "By pi (from V0's own walk), the pairwise top pairs of V1, V2 and V4 touch V0's hub genes more often than V0's do; for attention this holds at t = 9 but not, in K562, at each variant's own label-free depth (V1 0.571, V2 0.502, V4 0.559 against V0 0.608). Both connectedness rankings are tied to what is being compared (pi to V0's walk, the MI row mean to row centring), so neither is a neutral hub measure, and the shares have no intervals." |
| S8 | should-fix | B, opening paragraph and summary | "No p-values … none should be read into the intervals", yet "intervals excluding 0" is used as a criterion throughout. The pre-read shows 56 of about 735 intervals. | Add to the opening paragraph: "The two `tables.md` files hold about 735 such intervals; this page shows 56. With that many, 'interval excludes 0' is not evidence of a real difference and is used here only to describe." Optionally soften "with intervals excluding 0" to "with intervals on one side of 0". |
| S9 | should-fix | B, summary: "every robust variant is at or above V0 on M3D" | Pooled up, per-regulon down (−0.008 to −0.014 for every robust variant, pairwise). | "… every robust variant is at or above V0 on M3D in pooled enrichment but below V0 in mean per-regulon AP, …" |
| S10 | should-fix | Item 4 table footnote and `sota_table.md` notes | Item 8 now measures what the gene-only intervals omit. On yeast Y-M the experiment-level SE was 77% of the gene-level SE (group level) and 180% (series level). The method-table intervals hold networks fixed, so the "ties" (Bonferroni) are probably too few. | Footnote: "Intervals resample genes only; in the DREAM5 yeast test the experiment-resampling component was of the same size as the gene component or larger, so these intervals understate total uncertainty and some differences marked as wins may be ties." |
| N1 | note | Item 8, first bullet | "10% relative-AUPR margin (d = log2 ratio = 0.1375)" reads as if d were 0.1375. "Design frozen 4 Nov" without the 6 Nov scorer definition. | "… margin delta = log2(1.10) = 0.1375 on d = log2(AUPR ratio). Design frozen 4 Nov; the official-scorer details (trapezoid AUPR, gold-universe box, tie handling) were fixed on 6 Nov from the scoring code, before any score of ours, as the design allowed." |
| N2 | note | Item 8 | The preregistered expected-tie sensitivity gives d = 0.072 (interval 0.017 to 0.128). The primary `ap_fast` counts CLR's tied scores at the end of their block, which slightly favours attention. | Optional sentence in the full report (not the pre-read). |
| N3 | note | Item 8 / report | The exploratory analyses listed in the design (replicate-averaged input, tree module readouts, hub-stratified, t = 3 and 27, per-TF AUROC) are not run. Anything run now is after unblinding. | Label any later run "post-unblinding exploratory". |
| N4 | note | A table, K562 DoRothEA A–B bootstrap "3 to 24.2" | An interpolated percentile on a discrete grid. | Report grid values ("3 to 24") or state the interpolation. |
| N5 | note | A | Coarseness rank is nominal: K562 KD-response modules (median 100 genes) are ranked finer than DoRothEA (median 13). The E. coli bootstrap drops 48 of 1,000 draws (conditional interval). Operators are fixed in the bootstrap, so no experiment resampling. | One line in `docs/depth_log_results.md`. |
| N6 | note | Item 5, scPRINT bullet | The lead is "non-inferior and superior" against a comparator that scores at its own permutation-control level. The caveat is present but comes after the claim. | Lead with the caveat: "scPRINT's zero-shot attention scored at the level of its permutation control on these truths; against it, CLR-NS is non-inferior and superior …". |
| N7 | note | A, "Status of the evidence" | Good, but "Power is small" understates it: reading (a) was nearly unattainable. | Covered by M4. |

---

## Spot-checked numbers (follow-ups A and B, plus item 8)

| # | Pre-read value | Source | Recomputed / file value | Match |
|---|---|---|---|---|
| 1 | E. coli rho 0.500, one-sided P 0.500; K562 rho 0.700, P 0.117 | `results/granularity_log/depth_log_primary.csv` | 0.500 / 0.500; 0.700 / 0.117 | yes |
| 2 | Pooled slope 0.197 [0.078, 0.668]; E. coli 0.186 [0.054, 1.929], 952 draws; K562 0.380 [−0.203, 0.581] | `depth_log_slope.csv` | 0.1973 [0.0782, 0.6681]; 0.1862 [0.0539, 1.9290], 952; 0.3796 [−0.2025, 0.5810] | yes |
| 3 | K562 CORUM: best attention t = 16, AUPR 0.037 below CLR's 0.0894 | `curves_k562.csv` | t = 16 AUPR 0.05225; CLR 0.08943; difference 0.0372 | yes |
| 4 | Sigmulons bootstrap 0 to 256; CollecTRI within-1% 24 to 256; 0 of 8 censored | `depth_log_per_truth.csv` | 0 / 256; 24 / 256; censored all False | yes |
| 5 | Equivalence gate max 5.0e-16, 35 K562 rows | `equivalence_k562.txt` | 4.996e-16, 35 rows, PASS | yes |
| 6 | V1 pairwise CORUM all +0.547 [+0.345, +0.749] | `results/robust_z/k562/contrasts_pooled.csv` | 0.5472 [0.3452, 0.7492] | yes |
| 7 | V1 CORUM all, stratified: pi +0.014 [−0.178, +0.206]; MI row mean +0.483 [+0.234, +0.733] | same | 0.0142 [−0.1779, 0.2062]; 0.4831 [0.2336, 0.7326] | yes |
| 8 | Attention − CLR, CORUM all, V0: −2.88 [−4.32, −1.44] | same | −2.8796 [−4.3155, −1.4437] | yes |
| 9 | Hub drift (all 40 cells), e.g. K562 attention pi V1 0.711 vs V0 0.608; V3 MI 0.034 (M3D pairwise) | `hub_drift.csv` (both) | 0.71087 / 0.60759; 0.03431; all 40 cells agree | yes |
| 10 | Held-out t*: V0 9; V1, V2, V4 6 to 7; V3 5 to 6 | `tstar_summary.csv` | M3D 9/7/6/6/6; K562 9/6/6/5/6 | yes |
| 11 | Row scale vs strong partners: SD 0.86/0.89, IQR 0.72/0.85, MAD 0.71/0.83; top K562 decile SD/IQR 1.35, SD/MAD 1.6 | `row_scale_spearman.csv`, `row_scale_deciles.csv` | 0.863/0.890, 0.717/0.846, 0.706/0.831; 1.351, 1.639 | yes |
| 12 | V0 vs baselines run "within 4e-16" | `check_vs_baselines.csv` | max 3.9e-16 | yes |
| 13 | Y-M d 0.082 [0.026, 0.137], upper 0.1374 vs 0.1375; Y-E d 0.060 [−0.093, 0.213]; series SD 0.040 vs 0.017 | `results/dream5/truth_net4/verdict.csv`, `supplementary_series_jackknife.csv`; replicate SDs recomputed | 0.0817 [0.0260, 0.1374]; 0.0599 [−0.0928, 0.2125]; 0.0402 / 0.0173; SD(jack) 0.022549 = SE_J, SD(refit) 0.017314 = SE_B | yes |
| — | B summary: "Only V3 moves the hub share down, and only by the MI row mean" | `hub_drift.csv` | V3 also down by pi in K562; V2 down by MI row mean in M3D pairwise and K562 attention | **no (M5)** |

All 13 numerical checks match. The one mismatch is an interpretive sentence that contradicts its own table.

---

## What is solid

- **The DREAM5 process is exemplary.** A fixed design with a change log, a sealed truth, a step-6 code review by someone who had seen no score, a scorer gate reproducing all 36 published methods to 1e-3, hash-verified score files and refits, one primary run, and verdicts written by the script before interpretation. The pre-read reports the non-replication with the fixed wording and keeps net3 and net1 descriptive and labelled. Few published benchmark papers meet this standard.
- **The Y-M conclusion "did not replicate at the 10% margin" holds** under every resampling and sensitivity I tried (the point estimate is below delta). Y-E "Inconclusive" holds under all of them.
- **Follow-up A is reproducible and its scope is stated:** the equivalence gate (5e-16), deviations listed with timing, curves byte-identical across reruns, and reading (b) reported literally with the power limits stated.
- **Follow-up B is labelled correctly** (exploratory, not held out, not for the paper). Its code was fixed before results, V0 reproduces the package to 1e-13, and its overall conclusion ("no variant consistently better; not a recommendation") is the right one. It also adds value elsewhere: a stratifier that does not use the walk confirms the M3D strong-evidence attention lead over CLR (S5), and it shows the K562 DoRothEA hub-robust lead is stratifier-dependent (S4).
- **The E. coli headline itself** (+0.022 [+0.019, +0.026], 1,000 of 1,000 nested draws, laboratory splits, replicate-aware bootstrap) is robust to resampling within its data lineage. The pre-read's handling of EcoMAC overlap, of the attention-derived stratifier and of gene-only intervals is candid.
- **Negatives are reported**, including failed gates and failed preregistered tests. The one omission is CLR(log P^t) (S6).

---

## Appendix: computations

- Monte Carlo error of the Y-M upper limit: 20,000 resamples (with replacement) of the 200 truth-jackknife and 100 refit d values, recomputing d + 1.96 × sqrt(SD_J^2 + SD_B^2): mean 0.13716, SD 0.00220, P(≥ delta) 0.438.
- Series-level alternative: d ± 1.96 × sqrt(SE_J^2 + SD_series^2). Y-M [−0.0087, 0.1720]; Y-E [−0.1071, 0.2269].
- Bonferroni across the co-primaries (z = 2.241): Y-M [0.018, 0.145].
- Expected-tie Y-M: d = log2(0.11112 / 0.10569) = 0.0723, interval with the same SE [0.0166, 0.1281].
- Depth slope, least squares across the 7 truths with t_best ≥ 1: 0.197 (SE 0.376), 95% CI −0.768 to 1.163, P = 0.62. With an organism intercept: 0.230 (−0.963 to 1.422). Leave-one-out slopes: −0.217 (operons dropped), 0.324, 0.349, 0.165, 0.141, 0.218, 0.228.
- Hub drift at own t* (`hub_drift.csv`, kind `attT`, pi): M3D V0 0.257, V1 0.305, V2 0.270, V3 0.244, V4 0.283; K562 V0 0.608, V1 0.571, V2 0.502, V3 0.382, V4 0.559.
