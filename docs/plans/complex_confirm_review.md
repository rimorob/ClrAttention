# Independent statistical review of the draft complex test (13 Nov 2026)

Reviewer: independent statistical review (Claude Opus 5.5, reviewer role, run as a worker), no part in drafting. Scope: `docs/plans/complex_confirm_design.md` (DRAFT rev 1, commits 81de68c, 060bb30, def8ca7, 8533d98) before any truth-stage compute. The review now also covers the PI's G0 decision, which the lead relayed during the review (quoted in section 0). Read-only for every project file except this one. **No readout was scored against any truth.** Only counts and projections from stored outputs of spent runs were made. Scratch scripts are in `/tmp/claude-1000/-home-boris-Projects-ClrAttention/4451c478-b0bc-4b2f-84d1-11ecf74c2d08/scratchpad/`:
- `neff.R`: n_eff, the k rule and an 80S variant, from the frozen member file;
- `multi.R`: duplicate-ORF accessions;
- `pp_fix.R`: the power surrogate with zero-positive anchor complexes dropped, written to `ppout/`.

Sources read:
- the design (all of it);
- `docs/preregistration_reciprocity_yeast.md` (disclosure table, SF3, verdict table, change log);
- `docs/plans/reciprocity_heldout_design.md` §7;
- `docs/plans/clr_logpt_interpretation.md` §§3, 4, 5.6, 7;
- `docs/clr_logpt_verdict/*.csv`;
- `docs/plans/dream5_design.md` (with its change log);
- `results/dream5/truth_net4/{verdict,point_estimates,supplementary_series_jackknife}.csv`;
- `docs/plans/stats_review_2026-11-12.md`;
- `docs/explore_stack_results.md` and `docs/STATUS.md`;
- `docs/complex_confirm_selections/*` and `analysis/complex_confirm/*`;
- `results/clr_logpt_test/truth/*/{per_complex_ap,jack_replicates}.csv`;
- `analysis/regulons.R::comembership_pairs` and `analysis/clr_logpt_test.R::jack_dataset`;
- `analysis/explore_stack/stack_lib.R`;
- the raw Complex Portal file (participants and evidence-code columns only).

---

## Executive summary

**Go after fixes.** The draft is careful, honest about its history and rigorous in its protocol: frozen truth with a SHA-256, exact implementation gates, one run, synthetic-only smoke tests, a code review by someone who has seen no score. It correctly refused to treat G0 as resolved.

The PI has now chosen an exploratory estimation study. The draft is still written as a confirmatory test (verdict table, p-value, "eligible for adoption", stopping rule), so it must be rewritten in one pass before FIXED (B1). G0 must be recorded together with its consequences (B2).

Four major items should go into the same pass:
- net4's own experiment-level variance, which is 0.77× (group level) to 1.78× (series level) the truth variance on DREAM5 Y-M. The draft uses the E. coli ratio of 0.41 instead (M1).
- The 3,816 cross-subunit ribosomal pairs scored as negatives (M2).
- An over-calibrated power section (M3).
- The pooled estimate with E. coli (M4).

The label-free gates (Ga, Gb) may run now. Nothing found requires new data, and no finding gives a reason to doubt the frozen truth's integrity.

---

## 0. The PI's G0 decision, and how the lead reads it (addendum requested by the lead)

PI, verbatim (relayed by the lead, 13 Nov): **"yeah, let's do it in the exploratory fashion, use those truths"**.

| Item | Lead's reading | Reviewer's assessment |
|---|---|---|
| 1. Yeast-truth exclusion | Lifted for this test. Complex Portal yeast and net4 may be used. | **Agree.** "Use those truths" is explicit. The 1 Nov exclusion concerned truths, not the net4 matrix. Record three things: (i) the lift covers this study only; (ii) under the 1 Nov per-truth rule ("consumed once any operator has been scored against it"), **Complex Portal 559292 is consumed by this run for every complex readout**; (iii) the other two 1 Nov truths (YEASTRACT, TF edges) are not used and stay as they were. |
| 2. Serial-testing declaration | Not engaged, because the test is exploratory. The design must say plainly that it is not confirmatory and that no later claim may call it that. | **Agree, with a corollary the record must state.** The 1 Nov declaration ("no further … operator change aimed at complexes is tested confirmatorily") **stays in force**. Combined with (ii) above, this means **a3r is not adopted whatever this study shows**: no later confirmatory test of a3r on complexes is possible without a new PI withdrawal and a new untouched truth. The serial-testing ledger must still list this study as one more (exploratory) look at a3r on complexes. |
| 3. Power shortfall | Downgraded to a preregistered estimation study: no CONFIRMED / NOT CONFIRMED labels; estimates with intervals, the per-complex pattern, and which effect sizes are excluded. | **Agree.** Do not replace the verdict with free interpretation. Keep a **mechanical, pre-declared set of compatibility statements**, chosen by the script (section 5). Report **no p-value** for the primary estimate. "Exploratory fashion" plus a fixed, committed plan is stricter than the PI's words require, and that is the right direction. |

**(b) Does calling it exploratory remove the benefit?** No, provided every protocol constraint stays as written. These are:
- the plan committed and FIXED before compute;
- the truth frozen by SHA-256;
- one run;
- one primary estimand at t = 9;
- every pre-declared secondary reported;
- the "does not count" list kept.

The study's real value survives the relabelling. a3r has never been computed on net4, the truth has never been scored, and the plan is fixed, so the estimate of a3r − attention is **selection-free and out of sample**, with known precision. That is worth having. What is lost is only the right to an error-controlled claim and to adoption, and adoption was already off the table (item 2).

The benefit would be lost if "exploratory" were used to relax the protocol. Examples: a re-run at another t; a rule change after the estimate; leading with whichever sensitivity looks best; a post hoc analysis reported without its label. The edits below keep every constraint and change only the inferential layer.

---

## Verdicts on the attack list

### (1) The 1 Nov exclusion and serial testing (G0 items 1–2): **SOUND (as drafted); FIXABLE (disclosure accuracy)**

The draft is honest. It quotes both binding texts, sets out two readings, says only the PI can choose, and holds the truth stage. It did not let the lead resolve G0. Two inaccuracies in the permanent record should be fixed now, even though the PI's decision makes them moot:

- **The paraphrase of reading (a) borrows a clause from another branch.** "It may be used later only under a new preregistration that discloses three things …" belongs to the *complex-power-gate failure* branch of the 1 Nov file (complex membership built, H1 scored, complex truth unscored). Gate F failure is a separate clause: "consumes nothing, because no truth is downloaded before it". Reading (a) is still defensible: "therefore can never use" follows from "single-use", and its premise failed. But it rests on an inference from the premise, not on the disclosure clause.
- **An omitted fact weakens reading (a).** The 1 Nov design memo (`reciprocity_heldout_design.md` §7, item 1) put the choice to the PI explicitly: "the same preregistration must also decide whether it includes … a3r's last test (the interpretation memo's §7.3 asked for yeast) … Whatever is left out cannot use yeast later." The PI then left a3r out. The exclusion was a considered choice, not an accident of wording. The disclosure must cite this.

Options the PI had, and their credibility cost (for the record, now that the decision is made):
- (i) Run as confirmatory under reading (a) and withdraw the 1 Nov declaration. This is the highest cost: a fourth complex-oriented confirmatory test, the second for a3r, run after the PI's own "last" declaration, on the truth the PI had excluded from a3r. A referee would call it testing until it passes.
- (ii) Do not run. No cost, and no information.
- (iii) Run as a preregistered exploratory estimation. **This is what was chosen.** It costs little in credibility if it is labelled honestly and never upgraded. The 1 Nov declaration stays intact, and the selection-free estimate is still obtained.

### (2) Researcher degrees of freedom in the truth: **FIXABLE (one major gap)**

- **Molecule-set expansion: SOUND.** It is the only correct reading of `[ACC1,ACC2](n)`, since both paralogs are genuine subunit members. It was fixed before any score, is documented, and does not touch any spent truth (the E. coli file has no sets, verified in the design and consistent with the builder code). The alternative, dropping the sets, is in effect the pre-declared `noribo` variant, so that alternative is covered.
- **Same-set paralog exclusion: SOUND in principle; minor inconsistency.** It is 48 pairs (0.36% of positives) and cannot move the estimate materially. The run-time path is correct: `comembership_pairs()` treats NA as "no operon", so genes outside a set are never excluded. But the stated rationale (cross-hybridisation of near-identical paralogs) also covers **4 plain participant accessions that encode identical proteins from two ORFs**: P02994 TEF1/TEF2, P10081 TIF1/TIF2, P02309 HHF1/HHF2, P61830 HHT1/HHT2 (`multi.R`; 5 such ORFs in the main truth). These pairs stay positive. Either extend the rule to "ORFs mapped from one UniProt accession form a paralog set" and re-freeze the member file, or state the exception. Both are label-free (m3).
- **Ribosome subunits as separate complexes: not discussed; major (M2).** CPX-1599 (40S) and CPX-1601 (60S) share no gene (Jaccard 0), so the merge never joins them. Their **72 × 53 = 3,816 cross-subunit r-protein pairs are negatives** in the main truth, and they are among the most tightly co-expressed pairs on yeast arrays (one RP regulon). The mitochondrial 54S/37S subunits add 46 × 33 = 1,518 more.
  - This follows the source and matches the E. coli truth (30S CPX-3802 and 50S CPX-3807 are separate there).
  - It is the **opposite of K562 CORUM**, where a3r was selected: there the merged truth contains "Nop56p-associated pre-rRNA complex + Ribosome, cytoplasmic" (115 genes), so cross-subunit pairs were positives.
  - A readout that restores scale across the whole RP clique, which is a3r's mechanism, gains positives within each subunit and false positives across them. **A different, equally defensible rule could move d.**
  - With the cytoplasmic subunits joined (80S), positives rise to 17,255, the largest share to 0.45 and n_eff falls to **4.7**, which is E. coli-like concentration (`neff.R`).
  - The rule should stay as built (source-faithful, precedent-consistent). Before FIXED, disclose the fact and the K562 contrast, and add pre-declared descriptive sensitivities: (a) cross-subunit pairs of the cytoplasmic and mitochondrial ribosomes **excluded** from the universe; (b) subunits **joined** (80S; 54S + 37S).
- **Ribosome list (≥ 30% RP ORFs by the SGD regex): SOUND.** It is the reciprocity file's rule, label-free, and yields exactly the two subunits.
- **Universe = pairs among complex genes: SOUND.** It is the E. coli convention. The arithmetic checks: C(1,602, 2) − 48 = 1,282,353.

### (3) Independence of net4 and decoys: **SOUND; disclosure adequate with two additions**

- I found no a3r code path in `analysis/dream5*`. The committed scoreboard plan (9d79c6a, line 36) excludes net4 and every candidate truth. The untracked scoreboard work belongs to another agent and was not opened. The runner should confirm before job 6010 that no a3r net4 score file exists anywhere in `results/` (m9).
- The attention and CLR networks on net4 are known objects (regulon AUPRs 0.111 and 0.105; `point_estimates.csv`), but their complex AUPRs are not. Knowing that net4 is near chance on edges does not inform a3r − attention on complexes in any exploitable way. Choosing A over B used only label-free facts (data on disk, Kemmeren's gate-F failure).
- Additions to the disclosure:
  - the evidence profile of the truth, which resolves the design's [verify]. The `Evidence Code` column has 409 complexes with ECO:0000353 (physical interaction evidence) and 240 with ECO:0005546/0005547 (curator reconstruction codes; check the exact ECO term labels). None is expression-based.
  - the fact that Complex Portal yeast is consumed after this run (section 0).
- **Keeping the 283 decoys: SOUND.** It is required for the exact Ga part 4 equality with the DREAM5 caches and was fixed in the DREAM5 design. Decoys have mean |r| about 0.03 (DREAM5 design), cannot be truth genes, and enter both readouts through one operator. They do enter a3r's row standardisation as reachable columns, so add the decoy share of P⁹ row mass to Gb's reported (non-gating) diagnostics (N2).

### (4) The power model: **FIXABLE (major: overstated calibration, wrong experiment-variance ratio)**

Arithmetic, recomputed:
- k/√15.08 for the four anchors = 0.00986, 0.00945, 0.01056, 0.00966. With the fractional n_eff of 14.48: 0.0096 to 0.0108.
- Power at +0.005/0.010/0.015/0.025/0.035 with SE 0.0099: 0.073, 0.171, 0.328, 0.714, 0.942.
- MDE = 2.80 × 0.0099 = 0.0277. Combined SE with ratio 0.41 = 0.0107, power 0.647. Surrogate median power 0.408.
- **All of these match the design.** The anchor SEs match the committed verdict files.

Problems:

1. **The k rule rests on about two independent points, not four.** M3D and EcoMAC share the truth (and 37.7% of arrays). K562 and RPE1 share CORUM. The "±5% spread across two organisms, three truths" is mostly agreement within each pair. A two-point power law gives an exponent of 0.445, not 0.5, and SE(15.1) = 0.0103. The projection also extrapolates beyond the anchors' n_eff range (max 10.4). "Best-calibrated" should read "simplest; two independent calibration points".
2. **The k-rule numbers are not produced by any committed script.** `power_projection.R` writes only the surrogate models. The 15.1 figure and the 0.0094–0.0106 range should be written to `power_summary.csv`.
3. **Surrogate defect.** In E. coli, 75 of 129 (M3D) and 71 of 124 (EcoMAC) merged complexes have zero positive pairs, because all their pairs are same-operon. `pc_of()` keeps them. The surrogate gives them weight from gene counts with δ = 0, which distorts the E. coli scale factors, and puts zero draws into the yeast size-class pool. With them dropped (`pp_fix.R`), the net4 main envelope becomes **0.0092–0.0193 (median 0.0128)** instead of 0.0105–0.0226 (median 0.0145). Median power at +0.025 becomes about 0.50 instead of 0.41. The two surrogates also disagree in direction (iid: E. coli anchors give the *lower* SE; empirical: E. coli anchors give the *higher* one), so neither is calibrated. The "median of 8 projections" should not be quoted as a number.
4. **Experiment-level variance on net4 is taken from E. coli when net4's own ratio exists (major, M1).**
   - The design adds SE_B ≈ 0.41 SE_J (the E. coli experiment bootstrap), and argues that truth sampling dominates when it rejects pooling.
   - On this same matrix, the DREAM5 Y-M test gave SE_B(group) = 0.0173 against SE_J = 0.0225 (ratio 0.77), and a series-level SD of 0.0402 (ratio **1.78**) (`results/dream5/truth_net4/verdict.csv`, `supplementary_series_jackknife.csv`; stats review 12 Nov, M1).
   - The ratio is scale-free. In E. coli it was nearly the same for complexes and regulons (0.41 and 0.44), so net4's regulon ratio is the best available guide.
   - Projected total SE: **0.0125 (group level), 0.0203 (series level)**. Expected 95% half-widths: ±0.019 (SE_J alone), ±0.025 (group), **±0.040 (series)**. Power at +0.025 under these: 0.71, 0.52, 0.23.
5. **The expected yeast effect is probably below +0.025.**
   - Scaling the E. coli +0.027 by the ribosome share (0.29/0.60) gives about +0.013, before adding the other large complexes (71% of positives sit in complexes of ≥ 20 genes). The K562 noRP difference was +0.006.
   - The follow-up was triggered by a promising but unconfirmed E. coli estimate, which is a mild winner's-curse setting for the planning effect.
   - The cross-subunit negatives (M2) add uncertainty in both directions.
   - A planning range of +0.010 to +0.020 is more realistic, and there the power is 0.17 to 0.52 even on SE_J alone.

**Can power of 0.41–0.72 support a confirmatory claim?** No. The interpretation memo §7.3 item 2 said: "If the projected SE does not reach the target, the test is not run as confirmatory". On net4's own experiment variance the shortfall is larger than the draft states. **"Estimation study" is the more honest label**, and the PI has now chosen it. Replace the power section with a precision section (expected interval half-widths under SE_J, combined group and combined series).

### (5) The decision rule: **FIXABLE (superseded by the estimation framing)**

As a confirmatory rule it was coherent:
- one primary;
- one-sided α 0.025;
- d ≥ +0.005 as a minimum-effect condition, which never binds at SE 0.0099 because p < 0.025 already implies d ≥ 0.0194;
- the +0.025 cut separating "excluded" from "inconclusive", which is a sensible, pre-declared use of the planning effect.

Its wording table was good. Under estimation it must go. Replace it with **compatibility statements** that are computed by the script and written to `estimate.csv` before anyone reads them. On the primary interval [L, U]:

| Condition (on the primary 95% interval) | Fixed sentence |
|---|---|
| L > 0 | "The interval excludes no difference (preregistered exploratory estimation; not a confirmatory test)." |
| U < 0 | "a3r was lower than attention." |
| L ≤ 0 ≤ U | "The interval includes no difference." |
| U < +0.027 | add: "Effects as large as the E. coli estimate (+0.027) are excluded." |
| L > +0.005 | add: "The estimate exceeds the smallest effect of interest (+0.005)." |
| always | "a3r is not adopted; this study cannot change that (PI decision 12 Nov; G0 13 Nov)." |

The reference value +0.027 is the mean of P1 and P2 (`endpoints_all.csv`: 0.025659 and 0.028503). The draft's "+0.025, about the E. coli estimate" is a rounding.

- **No p-value** for the primary estimate. The interval carries the information, and a p-value invites "exploratory p = 0.01" to be cited as a test.
- **CLR as co-primary: no**, in the confirmatory sense. Under estimation, report a3r − attention (primary), a3r − CLR and attention − CLR **in one table with equal prominence**, with the E. coli S3 trap stated next to them: a3r − CLR has the smaller SE (k_CLR is 0.016–0.027 against 0.037–0.041), so it "excludes 0" more easily, and that is never a basis for "a3r beats CLR" (m7).

### (6) SE choice: **FIXABLE (add the series-level experiment unit; the rest SOUND)**

- **The gene jackknife as the primary truth-sampling SE: SOUND.** It is the paper's standard and the 26 Oct procedure. Genes are the right unit for overlapping complexes: a gene in several complexes moves all its pairs together, and the n_eff analysis shows that this block structure drives the SE. Delete-half with factor 1 is standard for smooth statistics. The U-statistic overstatement (memo §3.2) is conservative and correctly left uncorrected. Halves are drawn from the truth's genes (26 Oct convention), not from all genes (DREAM5 convention), which keeps comparability with E. coli. Keep it.
- **100 → 500 replicates: SOUND.** It cuts the Monte Carlo error of the SD from about 7.1% to 3.2% (1/√(2(B − 1))) and does not change the estimand, as the 12 Nov review (M1) recommended.
- **Normality.** With n_eff of about 15 and two blocks carrying 29% of positives, the replicate distribution of d may be skewed. Report its skewness and a percentile interval beside the normal one, descriptively (m4).
- **Complex-level jackknife as a sensitivity: SOUND as labelled.** It also changes the universe and base rate per replicate, and with two dominant complexes it may be noisy. It must never decide anything.
- **Experiment resampling: major gap (M1).** The draft has only the group-level unit (321 replicate groups). The 12 Nov review showed that chips within a `#Experiment` series are correlated: the series SD was 2.3 × SE_B for Y-M. The series refit set and its hashes exist (`results/dream5/net4/refit_series/`, 59 series, seeds 20264000 + b, DREAM5 change log 2026-11-05 item 4).
  - Add the series-level jackknife (100 replicates, the same halves, hash-checked for CLR and attention).
  - Pre-declare which interval the compatibility statements use. **Recommended:** the total interval d ± 1.96 √(SE_J² + SE_B,series²) is primary, because the estimand generalises over complexes *and* compendium conditions. Report the SE_J-only interval in the same row as the E. coli-comparable line.
  - **Acceptable alternative:** SE_J primary, with every statement also printed under the series-total interval and marked wherever it holds only under SE_J.
  - The draft's own rule ("the two components add; choosing the smaller one would be choosing the answer") argues for the first option on net4.

### (7) Gates, implementation checks, memory: **SOUND with minor additions**

Ga reference values checked against stored files:
- K562 a3r 0.105435083 (7.953×; the file labels this truth `CORUM_raw`, the design says "CORUM_all");
- CORUM a3r − K10 +0.0589 [+0.0334, +0.0844], printing as +0.059 [+0.033, +0.084];
- M3D a3r complexes 0.192091, SC 0.162220, P1 d 0.025659, SE 0.018169;
- `stack_lib.R`, `ideas47_lib.R` and `robust_z_lib.R` have no diff between d8fdb58 and HEAD;
- `SCORE_MANIFEST.tsv` and the net4 caches exist.

Additions:
- **Ga part 4b (m1).** The a3r path (`stack_lib.R`, `rz_operator`) and the attention cache (DREAM5 `build_operator`) are two code paths. Ga part 4 checks the new attention matrix against the cache. Add: (P⁹ + P⁹ᵀ)/2 computed **from the very P⁹ object passed to `sk_a3r`** equals the DREAM5 attention cache to 1e-12. That proves the two readouts share one operator.
- **Fallback B (m2).** Under an estimation framing, falling back to a matrix that failed the reciprocity gate F (mutant reproducibility 0.122 against 0.296 in wild type) gives an estimate of doubtful meaning. Recommended: no fallback, so that a Gb failure on net4 means "not evaluable". Gb is very likely to pass, since the DREAM5 attention operator was already built on net4. If B is kept, the gate-F failure must be printed in the estimate row.
- **Memory (m5).** Declarations follow the 1.4× rule (26, 10, 16, 12, 4 GB). Job 6030 now needs the series refits as well (about double the 3–5 h wall), and a3r adds log P⁹ and Z per worker. The DREAM5 net4 truth stage peaked at 6.4 GB for two methods with refits. Re-estimate 6030 to about 11 GB and declare 16, or measure one refit in 6010 and declare 1.4× that. Keep `gc()` before forks.
- **Overwrite semantics (m5).** Job 6040 rewrites "`verdict.csv` final", while the truth-mode guard requires "no `verdict.csv` exists yet". Have 6020 write `estimate.csv` (primary line, never rewritten) and 6040 write `estimate_final.csv` with the sensitivity lines.
- Nothing else could invalidate a single run, given the guards (FIXED status, member SHA-256, cache SHA-256s, seeds, `--jack 500`, `--t 9`, one run except for a crash or a named-line bug).

### (8) Claims after each outcome, and "last confirmatory test": **FIXABLE**

Under the PI's decision there is no confirmatory verdict. The paper may say, with the compatibility sentence and the interval:

> "In a preregistered exploratory estimation on an independent yeast compendium (DREAM5 network 4) and complex truth (EBI Complex Portal), with the plan fixed before any computation and one run, a3r exceeded [or: did not exceed / was below] attention by d (95% interval L to U). The readout's preregistered E. coli test was not confirmed (26 Oct), and a3r is not adopted."

It belongs in the supplement (Note 15), with at most one main-text clause.

**Never allowed, whatever the outcome:**
- "confirmed", "replicated", "validated", "significant", or a p-value;
- "a3r improves complex recovery";
- any tally with K562, RPE1 or E. coli as converging evidence;
- the pooled yeast + E. coli estimate (M4);
- "a3r beats CLR";
- any regulon claim;
- calling this study confirmatory in any later document.

The interpretation memo's §5.6 list of wording to avoid applies unchanged.

**"Last confirmatory a3r test".** The phrase no longer fits. Replace it with: "This is not a confirmatory test. The 1 Nov declaration (no further confirmatory complex-oriented operator test) stays in force. This is the last evaluation of a3r on a new complex truth in this project." The last sentence is advice: it guards against a sequence of "exploratory" looks. That is wise, because the paper does not need a3r, and the yeast complex truth is consumed after this run.

### (9) Other findings, and numbers checked

- **Pooled descriptive with E. coli (M4).** Secondary item 10 is the likeliest route to an overclaim: "pooled estimate excludes 0" would become the headline. The E. coli estimate is the one that triggered this follow-up, and memo §7.3 item 6 requires the dependence model to be specified in advance. Drop it.
- **Consumption (part of B2).** After this run, no other complex readout can use the yeast complex truth as untouched. Examples: the multi-head operator listed as "[PENDING]" in memo §5.4, and any reciprocity revival. The PI should know this when G0 is recorded.
- **"No zero-variance genes"** is stated while the drafter says no net4 expression value was read. Cite the source (the DREAM5 label-free logs).
- The anchor table's "Merged complexes 129 / 124" for E. coli counts complexes without positives. Say "129 (54 with positives)" and "124 (53)", as in the text.

| # | Design value | Source | Recomputed | Match |
|---|---|---|---|---|
| 1 | net4 main: 262 complexes, 1,602 genes, 1,282,353 pairs, 13,439 positives, n_eff 14.48, largest 0.188, ribosome 0.289 | `yeast_complexes_counts.csv` | same | yes |
| 2 | n_eff anchor-comparable 15.1 | member file, overlap-counted | 15.083 | yes |
| 3 | k rule SE 0.0094–0.0106, central 0.0099 | `power_anchors.csv` k values / √15.08 | 0.00945–0.01056, mean k → 0.00988 | yes |
| 4 | Power 0.07/0.17/0.33/0.72/0.94; MDE 0.028; combined 0.0107 → 0.65 | normal arithmetic | 0.073/0.171/0.328/0.714/0.942; 0.0277; 0.0107 → 0.647 | yes |
| 5 | Surrogate 0.0105–0.0226, median 0.0145, power 0.41 | `power_summary.csv` | same (but see (4)3: 0.0092–0.0193, median 0.0128 after the fix) | yes, defect noted |
| 6 | Anchor SEs 0.0182/0.0174/0.0128/0.0116; k_CLR | `jack_replicates.csv`, `endpoints_all.csv` | 0.018169/0.017445/0.012753/0.011595 | yes |
| 7 | E. coli ribosome about 60% | memo §3.2 (451 of 1,150 noribo) | 60.8% | yes |
| 8 | K562 noRP a3r − attention +0.006 | `metrics_robust.csv` CORUM_noRP | +0.00568 | yes |
| 9 | Ga: K562 0.1054351 / 7.95×; CORUM +0.059 [+0.033, +0.084]; M3D 0.19209, 0.16222, 0.025659, 0.018169 | `t9_jackknife_EXPLORATORY.csv`, `metrics_robust.csv`, `endpoints_all.csv` | same (truth label `CORUM_raw`) | yes |
| 10 | E. coli combined d 0.0271, SE 0.0178 | `endpoints_all.csv` | 0.02708, 0.01781 | yes |
| 11 | Y-E AUPR 0.019–0.022, base rate 0.017 | `truth_net4/point_estimates.csv` | 0.0189–0.0216; 3,940/227,202 = 0.0173 | yes |
| 12 | Member SHA-256 146030f7… | `sha256sum` | same | yes |
| 13 | Library files byte-identical to d8fdb58 | `git diff d8fdb58 HEAD` | empty | yes |
| — | Experiment variance ratio 0.41 used for net4 | net4's own Y-M: 0.77 (group), 1.78 (series) | — | **wrong source (M1)** |

---

## Findings

Severity: **blocker** = must be fixed before FIXED; **major** = a referee would seize on it, or it changes what the numbers mean; **minor** = worth fixing in the same pass.

| ID | Severity | Verdict | Issue | What to change |
|---|---|---|---|---|
| B1 | blocker | FIXABLE | The draft is a confirmatory test (verdict table, p_sup, "one primary test", "eligible for adoption", stopping rule, G0 item 3 rewording), while the PI chose exploratory estimation. | Rewrite the inferential layer per R3–R9. Keep every protocol constraint unchanged. |
| B2 | blocker | FIXABLE | G0 is unrecorded, and its consequences are not stated. | Change-log entry with the verbatim quote, the date and the relay. State: the lift is single-study; Complex Portal 559292 is consumed by this run for all complex readouts; the 1 Nov declaration stays in force; a3r is not adopted whatever the estimate; never to be called confirmatory (R1–R2). |
| M1 | major | FIXABLE | Experiment-level variance is projected from E. coli (0.41), while net4's own ratios are 0.77 (group) and 1.78 (series). Series-level resampling is omitted despite the 12 Nov review. Projected total SE 0.0125–0.0203. | Add the series-level jackknife (59 series, existing refit hashes). Pre-declare which interval the statements use (recommended: SE_J ⊕ SE_B,series). Precision section on net4's ratios (R6, R8). |
| M2 | major | FIXABLE | 40S/60S (and 54S/37S) cross-subunit pairs, 3,816 (+1,518), are negatives. This differs from K562 CORUM (cytoplasmic ribosome merged) and matches E. coli. Undisclosed, no sensitivity. The 80S join gives n_eff 4.7. | Disclose. Add two pre-declared descriptive truths: cross-subunit pairs excluded; subunits joined. Re-freeze the member file (R5). |
| M3 | major | FIXABLE | Power calibration overstated: about two independent anchors; k rule not in a committed script; surrogate keeps zero-positive E. coli complexes (envelope 0.0105–0.0226 → 0.0092–0.0193 when fixed); the median of 8 is meaningless. | Replace with a precision section; move the k rule into `power_projection.R`; fix or caveat the surrogate; drop the "median of 8" figure (R6). |
| M4 | major | FIXABLE | Secondary item 10, the pooled yeast + E. coli estimate, invites "pooled excludes 0". | Delete it; list it under "not computed" (R7). |
| m1 | minor | FIXABLE | The a3r and attention code paths are separate. | Ga part 4b: (P⁹ + P⁹ᵀ)/2 from the P⁹ fed to `sk_a3r` equals the DREAM5 attention cache to 1e-12 (R10). |
| m2 | minor | FIXABLE | Fallback B is a matrix that failed gate F. | Recommended: no fallback (Gb failure → not evaluable). Otherwise print the gate-F failure in the estimate row (R10). |
| m3 | minor | FIXABLE | The paralog rationale also covers 4 identical-protein accessions (TEF1/2, TIF1/2, HHF1/2, HHT1/2) that stay positive. | Extend the rule to "ORFs from one accession" and re-freeze, or state the exception (R5). |
| m4 | minor | FIXABLE | The normal interval may be poor with two dominant blocks. | Report the skewness and percentile interval of replicate d, descriptively (R8). |
| m5 | minor | FIXABLE | Job 6030 memory and runtime do not include series refits or a3r matrices. Job 6040 rewrites `verdict.csv` against the guard. | Re-estimate 6030 (declare ≥ 16 GB or 1.4 × measured). `estimate.csv` (6020, never rewritten) plus `estimate_final.csv` (6040) (R11). |
| m6 | minor | FIXABLE | The reading-(a) paraphrase uses the power-gate clause for the gate-F branch. It omits that the PI explicitly declined a3r on yeast on 1 Nov. | Correct the disclosure; cite `reciprocity_heldout_design.md` §7 item 1 (R2). |
| m7 | minor | FIXABLE | CLR contrasts sit apart from the primary. The S3 trap is not stated. | Put a3r − CLR and attention − CLR in the primary table with the trap sentence (R7). |
| m8 | minor | FIXABLE | Labels and figures: "CORUM_all" (the file says `CORUM_raw`); anchor "129/124 complexes" (54/53 with positives); "+0.025 about the E. coli estimate" (0.027); "no zero-variance genes" without a source; ECO [verify] now resolvable. | Text fixes (R12). |
| m9 | minor | FIXABLE | Scoreboard work is running concurrently. | The 6010 guard confirms that no a3r net4 score exists in `results/` (R10). |
| N1 | note | SOUND | 500 replicates, seeds and draw order. | None. |
| N2 | note | SOUND | Decoys kept. | Add the decoy share of P⁹ row mass to Gb's reported diagnostics. |
| N3 | note | SOUND | The minimum-effect condition d ≥ +0.005 never binds at SE 0.0099. | Moot under estimation. Keep +0.005 only as a reference value. |

---

## REQUIRED edits before FIXED (for one editor, in file order, one pass)

**R1. Title and status.**
- Title: "Preregistration (exploratory estimation): how large is the CLR(log Pᵗ) (a3r) minus CLR attention difference on protein complexes in an untouched yeast compendium?"
- Status line: "DRAFT rev 2 → FIXED on <date>. Exploratory estimation, not a confirmatory test (PI, 13 Nov)."
- Keep the drafter's exposure list. Add: "Reviewer: `docs/plans/complex_confirm_review.md`."

**R2. G0 section.** Replace "must be resolved by the PI" with "Resolved 13 Nov 2026 by the PI" and the verbatim quote: "yeah, let's do it in the exploratory fashion, use those truths" (relayed by the lead). Then state:
- (1) The yeast-truth exclusion of 1 Nov is lifted for this study only. Complex Portal 559292 is consumed by this run for every complex readout (per-truth rule of 1 Nov). YEASTRACT and the TF-edge truth are not used.
- (2) This is not a confirmatory test. The 1 Nov serial-testing declaration stays in force. No later document may call this study confirmatory, a replication or a validation. It is listed in the ledger of looks as an exploratory look at a3r on complexes.
- (3) The power shortfall is accepted by making this a preregistered estimation. No CONFIRMED / NOT CONFIRMED labels and no p-value.
- (4) a3r is not adopted, and this study cannot change that.
- In the item-1 history, correct the reading-(a) paraphrase: the "may be used later … three things" clause belongs to the complex-power-gate branch; under gate F failure the file says "consumes nothing". Add the 1 Nov design memo §7 item 1 quotation ("Whatever is left out cannot use yeast later") and the fact that the PI left a3r out.

**R3. Question and history.**
- Replace "Does a3r beat plain CLR attention…?" with "How large is AUPR(a3r) − AUPR(attention) … ?", estimated with intervals.
- Replace "It is proposed as the last confirmatory test of a3r … stays an explanatory ablation for good" with: "Not confirmatory. The 1 Nov declaration stays in force. This is the last evaluation of a3r on a new complex truth in this project [if the PI agrees; see advice]. a3r stays not adopted."
- Add history item 5: "13 Nov: PI decision, exploratory estimation."

**R4. Candidates table and Choice.** Unchanged, except: under "Pre-declared fallback", either delete fallback B (recommended: a Gb failure on A means "not evaluable") or add "B failed the reciprocity gate F; that failure is printed in the estimate row". Correct "54 complexes" vs the anchor table's 129/124 (state "with positives").

**R5. Truth.**
- Add a paragraph: "The cytoplasmic 40S and 60S (and the mitochondrial 37S and 54S) subunits are separate complexes in the source and share no gene, so their 3,816 (and 1,518) cross-subunit pairs are negatives. This matches the E. coli truth (30S and 50S separate) and differs from K562 CORUM, where a3r was selected (the cytoplasmic ribosome is one merged complex there). Joined, the cytoplasmic ribosome would carry 45% of positives (n_eff 4.7)."
- Add two pre-declared descriptive variants to the builder and the table order: `main_xsub_excl` (cross-subunit pairs of both ribosomes removed from the universe) and `main_ribo_joined` (40S ∪ 60S and 37S ∪ 54S as single complexes).
- Optional, recommended: extend the paralog rule to ORFs mapped from one UniProt accession (P02994, P10081, P02309, P61830).
- Re-run `build_yeast_complex_truth.R`, re-freeze the member SHA-256, and update counts and the Gb part 5 values. This is structure only.
- Replace the truth–data [verify] with the ECO profile: 409 complexes ECO:0000353 (physical interaction), 240 curator reconstruction codes (ECO:0005546/0005547), none expression-based.

**R6. Power → "Precision (structure plus existing outputs only)".**
- Keep the anchor table, adding "(54 / 53 with positives)" for E. coli.
- State that the anchors form about two independent pairs. Give the two-point exponent (0.445; SE 0.0103 at n_eff 15.1).
- Report expected 95% half-widths: SE_J alone ±0.019 (0.0099); with net4's group-level experiment ratio (0.77, DREAM5 Y-M) ±0.025 (0.0125); with the series-level ratio (1.78) ±0.040 (0.0203).
- Give power at +0.025 under these (0.71 / 0.52 / 0.23) only as design sensitivity.
- Delete the E. coli 0.41 combined row as the net4 projection; keep it only as a comparison.
- Surrogates: either fix `pc_of()` to drop zero-positive complexes (corrected envelope 0.0092–0.0193) or caveat the defect. Delete the "median of 8 projections" figure.
- Make `power_projection.R` write the k rule (n_eff overlap-counted 15.08, SE per anchor) to `power_summary.csv`, and re-commit its outputs.
- Delete "This is G0 item 3". Replace the Plain statement with: "a3r − attention will be estimated to about ±0.02 (truth sampling only) to ±0.04 (truth plus experiment-series sampling). The planning effect is unlikely to exceed the E. coli +0.027 and may be about half of it (ribosome share 29% against 60%)."
- Keep the "options considered" list, but soften "the SE is dominated by truth sampling" for pooling A and B: on net4 the experiment component is comparable or larger. The pooling rejection stands on the shared truth and gate F.

**R7. Estimands and reporting rule** (replaces "Estimands and decision rule").
- Primary estimand: d = AUPR(a3r) − AUPR(attention), as defined. No test, and no p-value.
- Report in one table: d with the primary interval (R8), the SE_J-only interval, a3r − CLR, and attention − CLR. Add the sentence: "a3r − CLR has a smaller SE and may exclude 0 where a3r − attention does not; this is never a basis for 'a3r beats CLR' (E. coli S3)."
- Replace the verdict table with the compatibility-statement table of this review's section (5), with reference values +0.005 and +0.027. Delete the "Multiplicity" and "Licenses / Does not license" columns as written. Keep the full list under "Does not license" as "Never allowed".
- The script writes the selected sentences to `estimate.csv` before anyone reads them.
- Rename "What does not count" to "Not the primary estimate; may not be cited in its place". Keep the list.
- Delete secondary item 10 (pooled with E. coli). Add it under "Not computed".

**R8. SE section.**
- Rename "Decisive: SE_J" to "Primary truth-sampling SE: SE_J".
- Add **Sensitivity 3, series-level experiment jackknife**: delete-half over the 59 `#Experiment` series; 100 replicates; seeds 20264000 + b; the same halves as the DREAM5 series refits; CLR and attention refits checked against `results/dream5/net4/refit_series/` hashes; SE_B,series = SD of the paired differences.
- Fix which interval the compatibility statements use:
  - (A, recommended) d ± 1.96 √(SE_J² + SE_B,series²), with the SE_J-only interval in the same row;
  - or (B) SE_J, with every statement also printed under (A) and marked where it holds only under SE_J.
- Report the skewness of the 500 paired replicate differences and their 2.5% and 97.5% percentiles, descriptively.
- The complex-level jackknife and group-level SE_B stay labelled sensitivities.

**R9. "What is NOT claimed" and Limitations.**
- Add: "This is not a confirmatory test and must never be described as one, or as a replication or a validation. No p-value. No adoption. Complex Portal yeast is consumed by this study."
- Update the power limitation to the precision numbers of R6, including net4's experiment-series component.
- Add the cross-subunit negatives (R5).

**R10. Gates.**
- Ga part 2: refer to the K562 truth as `CORUM_raw` (the file's label).
- Add Ga part 4b: "(P⁹ + P⁹ᵀ)/2 computed from the same P⁹ object passed to `sk_a3r` equals `attention__replicate.scores` (SHA-256 or max |diff| ≤ 1e-12)."
- Gb: add the decoy share of P⁹ row mass to the reported diagnostics.
- 6010 guard: no a3r score on net4 exists anywhere under `results/` before the run.
- Fallback: per R4.

**R11. Jobs and outputs.**
- 6030: include the series refits; re-estimate memory (about 11 GB, declare 16, or 1.4 × a measured single refit) and wall time (about 6–10 h).
- 6020 writes `estimate.csv` (primary line, never rewritten); 6040 writes `estimate_final.csv`.
- The truth-mode guard checks for `estimate.csv` instead of `verdict.csv`.
- Protocol step 7: "the estimate row with its fixed sentences" instead of "the verdict row".
- Add: "Analyses not in this file may be run only after `estimate.csv` is committed, and are labelled post hoc."

**R12. Small text fixes.**
- "+0.025 (about the E. coli estimate)" → "+0.027 (mean of P1 and P2)" wherever it is the reference value.
- "No zero-variance genes": cite the source.
- The Independence table's serial-testing row: "exploratory look; not confirmatory; the 1 Nov declaration stays in force".
- Change log: add "2026-11-13: independent review (`complex_confirm_review.md`); rev 2 implements R1–R12" and the G0 entry of R2.

None of R1–R12 requires a score. R5 and R6 require re-running the two design-stage scripts, which only count and project.

---

## Advice to the PI on G0 items 1–3 (the PI decides; the PI has already chosen exploratory estimation)

1. **Item 1 (use of the yeast truths).** Lifting the exclusion for this one exploratory estimation is reasonable. Please confirm that you accept the consequence: once a3r and attention are scored, Complex Portal yeast is consumed for every complex readout in this project, including any later operator (for example the multi-head candidate) that might have wanted it as an untouched truth.
2. **Item 2 (serial testing).** Keep the 1 Nov declaration in force rather than withdrawing it, and state that a3r stays not adopted whatever this study shows. I also advise declaring this the **last evaluation of a3r on a new complex truth** for this paper, exploratory or not. Otherwise a sequence of exploratory looks would recreate the serial-testing problem that the label is meant to avoid. The paper does not need a3r: attention is the method, and a3r's role as the mechanism is secure either way.
3. **Item 3 (power).** The estimation framing is the honest one. The projected precision is worse than the draft says once net4's own experiment-series variance is counted (half-width up to about ±0.04). I advise making the **total interval (truth plus experiment series)** the one the fixed sentences use, with the gene-jackknife-only interval beside it for comparison with E. coli. Report no p-value.

---

## What is solid

- The process design is at the DREAM5 standard:
  - truth frozen by SHA-256 before any score;
  - exact implementation gates against stored outputs (all reference values check);
  - one run, with guards in code;
  - synthetic-only smoke tests;
  - a code review by someone who has seen no score;
  - seeds and draw order fixed.
- G0 was handled correctly: the draft surfaced both binding texts, refused to let the lead resolve them, and recommended waiting a day rather than risking a disowned test.
- The truth build is principled. The molecule-set rule was found while counting structure and disclosed, and it was the only correct reading. All structure counts reproduce.
- The power section is unusually candid for a draft: it says that the memo's SE ≤ 0.009 target is not reached and that NOT CONFIRMED could not mean "no effect". Its arithmetic is correct. The issues are calibration and the source of the experiment-variance ratio, not honesty.
- The "does not count" list, the E. coli history carried into every write-up, and the separation of a3r − CLR as untested are all right, and they carry over unchanged into the estimation framing.
