# Adversarial statistics review: DRAFT yeast reciprocity / H1 preregistration

*1 Nov 2026. Independent reviewer (Claude). Reviewed: `docs/preregistration_reciprocity_yeast.md` (DRAFT, commit 18babdb), `analysis/yeast_prep.R`, `docs/plans/reciprocity_heldout_design.md`. Context: the keff and CLR(log Pᵗ) preregistrations, `docs/plans/clr_logpt_interpretation.md`, `docs/plans/background_width_exploration.md` §8, `analysis/synthetic_width/sim_width.R` / `sim_round2.R` and one round-2 metrics file. No data downloaded, nothing scored, no truth opened, the preregistration not edited. Line numbers refer to the DRAFT preregistration unless a file is named.*

**Overall: NO-GO as written. GO after fixing the three blockers and the SHOULD-FIX items marked (high).** Every fix is pre-truth and cheap. The design is otherwise careful: the disclosure table, the IUT, gatekeeping, git gating, the single-use statement and the "last test" declaration are good.

---

## BLOCKERS

### B1. The regulon truth is very probably not held out from the expression data (P2, H1, S6). Lines 108–113, 21, 33.

- YEASTRACT's "expression evidence" is, by its own description, "mostly that obtained by comparing gene expression in wild-type strains with … TF mutant strains … DNA microarrays" (Teixeira et al. 2018, NAR, [doi:10.1093/nar/gkx842](https://doi.org/10.1093/nar/gkx842), retrieved via PubMed). The Kemmeren 2014 deleteome is the largest such dataset in yeast, and YEASTRACT curates large deletion-profiling studies. The deleteome also absorbs earlier Holstege-lab deletion profiles (van Wageningen et al. 2010 kinases/phosphatases; Lenstra et al. 2011 chromatin), which YEASTRACT may cite separately. Whether Kemmeren-derived associations are actually in the export has not been verified, and the preregistration does not ask.
- If they are, the P2/H1 positive pair (Y, Z) means in part "Y and Z both changed significantly in deletion column X *of this matrix*". Regulator-agnostic co-membership is exactly co-response across columns. So the truth is partly a function of the confirmatory data. The "DNA-binding **and** expression" intersection (line 110) does not fix this, because the expression leg can come from the deleteome alone.
- The responsive-mutant selection (line 64, ≥ 4 limma-significant changes) makes it worse. TF deletions that produced deleteome-derived YEASTRACT associations are exactly the responsive columns that are kept.
- It is not neutral between comparators. It raises absolute AUPRs, and it favours whichever score best captures sparse, few-column co-response. That is uncontrolled for both R vs K10 (P2, S6) and K10 vs CLR (H1). Above all, "held-out truth" (line 20) and "headline replicates in a second organism" (line 244) become false claims. This is the first thing a Nature-family referee will check.
- **Fix (fix it in the text before freezing; the rule must be chosen without seeing the export):**
  1. Primary P2/H1/S6 truth: YEASTRACT associations with DNA-binding evidence **and** expression evidence from at least one reference that is not Kemmeren 2014 (PMID 24766815), not one of its constituent Holstege-lab deletion datasets, and not a reanalysis of them. List the excluded references by PMID in the preregistration. This needs a per-association reference export. YEASTRACT's association search shows references; record the form and parameters.
  2. A fallback, fixed now, if references cannot be obtained: **DNA-binding-only** associations (ChIP/ChIP-chip/ChIP-seq/ChIP-exo/EMSA/footprinting) become primary, decided on availability alone before the export is opened. The current "both" definition is then descriptive only, labelled "possibly circular".
  3. Report the share of primary associations removed by the exclusion. That is truth structure only, so it can go in the truth_structure commit.
- **Complexes (lower risk, but add a rule):** Complex Portal curation rests on physical-interaction evidence, and CYC2008 (2009) predates the deleteome. Still, in `truth_structure`, drop any complex whose only cited evidence is one of the excluded PMIDs, using complextab's cross-reference and evidence columns. Log the count, and record the ECO evidence codes. Fix this rule now.

### B2. Verdict 5 turns "non-inferiority not shown" into "a regulon cost". Lines 241, 226–227.

- Verdict 5 ("P1 passes and P2 fails → small complexes improve **at a regulon cost**; reported as a trade-off") treats a failure to reject the NI null as evidence of harm.
- That is exactly the inference the CLR(log Pᵗ) preregistration avoided. Its table separates "regulon NI not shown → unresolved regulon cost" from "significant loss d ≤ −0.005 → trade-off". The design memo's verdict 5 also carried the significant-loss clause; the DRAFT dropped it.
- The risk is real. P2's SE is projected from *E. coli* keff tweaks (0.0015) with no yeast structure check, and YEASTRACT regulons of up to 500 targets give a few very large co-membership blocks. At SE 0.0025, NI power at d = 0 is 0.52 (memo §1). So verdict 5 could easily be issued under a true d = 0.
- **Fix.**
  - Split the verdict. 5a: P1 passes and P2 NI is not shown, with no significant loss → "not adopted; small-complex gain with an unresolved regulon question (NI not shown)". 5b: P1 passes and there is a significant loss (d ≤ −0.005 and one-sided p < 0.025 for d < 0) → "regulon trade-off".
  - Add to the power section a structure-only regulon n_eff (Σ pairs per regulon, the same formula) to the truth_structure output, reported and not gating.
  - State that a P2 failure under low power is not evidence of a cost.

### B3. The deciding statistic's experiment component is mis-specified for this compendium. Lines 206–213.

Once the file is fixed, the deciding SE is fixed for good. Two concrete defects can bias SE_B, one of them downward:

1. **Exact duplicate columns.** Drawing 700 averaged mutant profiles with replacement puts exact duplicate samples into the B-spline MI estimator: about 37% of columns are repeats, and the effective N is smaller than the nominal one. MI and CLR are nonlinear in N, so the bootstrap distribution of the paired difference is distorted in an unknown direction. The house M3D bootstrap avoided this by drawing *a distinct replicate array per occurrence* (CLR(log Pᵗ) prereg, Statistics). That is impossible here, because the published M is already averaged over the two cultures.
2. **Clustered experiments.** Kemmeren's central result is that deletions of subunits of the same complex or pathway give near-identical profiles (Mediator, SAGA, Rpd3L, …), and the slow-growth signature is shared by many mutants. Resampling mutants as i.i.d. units ignores this clustering and understates experiment-level variance.

- The sum SE_J² + SE_B² is otherwise defensible. In a crossed genes × experiments model, the two one-way components each contain the interaction term, so the sum double-counts it (Owen 2007, the pigeonhole bootstrap) and is conservative. Independence is not needed for that. Keep the combination; fix the experiment component.
- **Fix (pick one now):**
  - (a) A **delete-half experiment jackknife over label-free mutant clusters**. Cluster the 700 profiles by complete linkage on 1 − Pearson r at a fixed cut, say r ≥ 0.5, computed label-free and committed with the selections. Draw 200 random halves of the clusters, refit, and take SE_B = SD of the paired differences. With d = n/2 the delete-d factor is 1. No duplicates arise, and the clustering is respected. It is conservative, because half-N networks are noisier.
  - (b) A cluster bootstrap that, for repeated mutants, draws one of the two culture-level profiles rebuilt from the sample tables. This needs all of the roughly 1,400 mutant sample tables.
  - Either way, also commit the cluster count and the size distribution as a label-free diagnostic.

---

## SHOULD-FIX

### S1 (high). The synthetic evidence and the planning effect are for a different estimand from P1. Lines 27, 102–105, 229.

- The round-2 `small_x` metric (`sim_width.R` `score()`) has **all genes** in the universe and *all non-complex pairs* as negatives, at **t = 9**.
- P1 restricts the universe to pairs of genes in some merged complex and uses each operator's own t\*. "C1's size class, unchanged" (line 102) is true, but the universe is not unchanged.
- Reciprocity's mechanism (unreturned mass from partner-less rows stays on self) acts mostly on pairs that involve non-complex genes, and those pairs are all removed from the P1 universe. The ×1.86–6.24 range, and the +0.010 planning effect that "discounts" it, may not transfer.
- The same synthetic file shows something the disclosure omits. In `std_f10_b1_s1_metrics.csv`, CLR `small_x` is 1877, K10 is 1092 and R is 2089. So in the synthetic data K10 loses heavily to CLR on small complexes, and most of R's "gain" restores CLR-level performance.
- **Fix, synthetic only:**
  - Before freezing, rescore the existing round-2 runs with the P1 construction (universe = complex-gene pairs, larger-complex-only pairs excluded) at t = 9 and at a criterion-A t\* computed on lab halves. Report the ratios in the disclosure table and derive the planning effect from them.
  - Disclose the CLR > K10 synthetic small-complex pattern. A referee will ask "why not just CLR?", and R vs CLR on complexes is only descriptive (line 187).
  - Disclose the round-1 clean-scenario regulon cost of ×0.95 (0/4 seeds; §8.5 item 4 caveat), which line 27 omits.

### S2 (high). Gate F is uncalibrated and does not test what it is said to test. Lines 77–94; `yeast_prep.R` lines 199–223.

- **Trivial pass.** The chance overlap of two top-1% pair sets is 0.01. If WT variation is mostly technical, R_wt sits near 0.01–0.05, and R_mut ≥ 2·R_wt passes almost automatically. The mutant columns were also *selected* for having signal (≥ 4 significant changes against the WT controls), so a pass is nearly built in.
- **Spurious fail.** WT day-to-day variation contains real growth-rate / ESR programs (O'Duibhir et al. 2014), the same programs that dominate many deletion profiles. A high R_wt can then fail F even when the mutant network carries plenty of regulatory signal.
- **Asymmetric variability.** Mutant draws resample columns, WT draws only re-halve a fixed set, so "min R_mut > max R_wt" compares unlike dispersions.
- **Not inferential.** F says nothing about whether R vs K10 is estimable, and it does not "answer the Rosetta concern". The WT controls are not used to calibrate anything in the analysis.
- **Inconsistent fallback.** Line 84 says that if the implementation check fails, F uses rebuilt mutants. The script stops with "not yet implemented" (`yeast_prep.R` line 203).
- **Fix.**
  - Reword F as a sanity check, not a calibration.
  - Use a chance-corrected statistic, R′ = (R − 0.01)/(1 − 0.01).
  - Promote the already-reported "share of mutant top-1% pairs that are also top-1% in the full WT network" to the second gating condition, for example ≤ 0.5 (most of the mutant network is not WT-noise structure). Fix the threshold now and justify it in one sentence.
  - Either implement the rebuilt-mutant branch (download the sample tables for a fixed random 376 mutants) or state that failing the implementation check means verdict 1.
  - A failed F costs nothing (no truth is opened, so yeast is not consumed). Say so, so that F is not seen as an expensive gate.

### S3 (high). S1 decides the default operator, but it sits inside a Holm family with five human and edge tests. Lines 168–177, 238–239.

- Under Holm, whether S1 is rejected depends on S2–S6. S1 may face α/6 ≈ 0.0042, and HepG2/Jurkat outcomes can flip verdict 2 vs 3, the "becomes the method's operator" decision. That is an odd and avoidable dependence.
- An S1 that cannot be evaluated counts as "not shown", which yields verdict 3's "large-complex limitation documented" when nothing was documented.
- **Fix.**
  - Fixed sequence: S1 first at full one-sided α = 0.025, gated on P1 and P2. If S1 is rejected, carry on to Holm over S2–S6 at 0.025. If not, either still run Holm over S2–S6 at 0.025 (the fallback procedure) or stop; pick one and state it.
  - Add verdict wording for "S1 not evaluable".

### S4. The power projection is "conservative" only for the jackknife component. Lines 218–229.

- SE_J ∝ 1/√n_eff is anchored on an *E. coli* a3r contrast, with n_eff ≈ 4 [derived from assumed block shares, not computed]. Plausibly conservative for SE_J.
- The SE_B/SE_J ratio of 0.44 comes from M3D, a 466-condition multi-lab compendium. In a perturbation compendium a small complex's co-expression can hinge on a few deletion columns, so SE_B may match or exceed SE_J. SE_B also does not scale with the truth's n_eff.
- Sensitivity at n_eff = 100:
  - ratio 1: SE ≈ 0.0051, power at +0.010 ≈ 0.50, MDE ≈ 0.014;
  - ratio 2: SE ≈ 0.0080, power ≈ 0.24.
- **Fix.**
  - Add this sensitivity row.
  - Say the 0.72 figure assumes the M3D ratio.
  - Pre-commit to report the realised SE_J and SE_B against the projections. Report the MDE, never observed power.
  - Optionally gate on a label-free SE_B proxy: the bootstrap/jackknife SD of the criterion-A AP difference R vs K10 on the held-out top-1% target, computed at the selection stage.
- No structure projection exists for P2 or H1 at all (see B2).

### S5. No yeast analogue of the operon exclusion. Lines 100, 111.

- The *E. coli* headline truth excludes same-operon pairs because they are trivially co-expressed.
- In yeast, divergent (bidirectional) gene pairs share one intergenic region. A single ChIP peak there is assigned to both genes, which makes them regulon co-members, and they are also co-expressed through shared chromatin. Adjacent tandem pairs are co-expressed too.
- "No adjacency exclusion" (line 100) is reasonable for complexes, not for regulons.
- **Fix.** Exclude from P2/H1/S6 the pairs that share an intergenic region (divergent pairs, from SGD coordinates, label-free). Alternatively keep them in and make the exclusion a descriptive sensitivity. Choose now.

### S6. Deleted-gene masking is incomplete. Lines 71–74; `yeast_prep.R` lines 120–125.

- **Unresolved labels.** These are only logged, so the deleted gene's own value, around −5 log2, stays in its column. That huge outlier co-occurs with the deleted regulator's target response and creates a spurious regulator–target and regulator–responder signal: the selection artefact the brief worries about. Fix: drop the column. Also map through SGD aliases (column V6 of `SGD_features.tab`), since names have changed since 2014.
- **Overlapping and neighbouring ORFs.** Dubious ORFs that overlap the deleted ORF, and the neighbouring-gene effect of the KanMX cassette, leak the deletion into other rows. Fix: also mask, in each mutant's column, ORFs that overlap the deleted ORF or lie within a fixed window of it (SGD coordinates, label-free). Fix the window now, for example 500 bp.
- **Paralogs.** Near-identical paralogs (many RPL/RPS A/B pairs, histones) cross-hybridize, and UniProt sometimes maps one accession to two ORFs (`truth_structure`, `orf_of`). That creates technically co-expressed positives in S1 and the ribosome. Add a descriptive sensitivity that removes pairs mapped from a single accession.
- **Slow growth.** "No slow-growth correction" is a fixed, label-free choice, but give its rationale (consistency with K562 and M3D). Optionally add the authors' growth-corrected matrix as a descriptive sensitivity.

### S7. Protocol enforcement gaps in code and order. Lines 258–276; `yeast_prep.R`.

- `truth_structure` checks that F and the splits are committed, but not that the t\* selections (step 6) or the synthetic R check (step 5) are committed. Add both.
- The test script is reviewed at step 8, *after* the truth member lists sit on disk (step 7). Fix: commit the complete test script, smoke-tested on synthetic or decoy truths, **before** step 7. Record its hash in the truth_structure output. Any change after step 7 is a logged deviation.
- `truth_structure` does not compute the evaluability counts that verdict 1 depends on: P2 and S6 positive pairs, and P1/S1 complexes with positives. Add them.
- The YEASTRACT parser keeps only `m == 1`. If the matrix encodes negative regulation as −1, every repression edge is silently dropped. Use `m != 0`, and log the value table.
- The YEASTRACT export is manual. Require it to be done after the fix commit, and record its timestamp and the form parameters in the committed log.
- Add an explicit line: no R (reciprocity) score is computed on any real truth (K562, *E. coli*, human) before the yeast truth phase. Today this is only a statement of fact (line 16).

### S8. Operator, endpoint and statistics underspecified. Lines 123–131, 102, 113, 172–176, 202–205.

- **Support definition.** "1[Aᵀ > 0]" defines support by numeric positivity, while the memo's words say "in j's cap-50 support". They differ if exp() underflows to 0 inside the cap. Define support as the cap-50 positive-candidate set (an indicator matrix), and unit-test B itself, not only the synthetic AP metrics.
- **Self-loop order.** State that A already carries empty-row self-loops before `recip()` (the outcome is equivalent, but say so).
- **Jackknife "truth's genes".** For P1, say whether this means all merged-complex genes (the P1 universe) or only small-complex genes. Say whether H1 reuses P2's halves.
- **S6 regulator definition.** Define "regulator gene": any YEASTRACT TF, or only regulators with a 5–500 regulon.
- **`make_regulons()` adaptation.** That function expects RegulonDB columns (`bnumber`, `confidence`, `reg_class`). State the yeast mapping.
- **S2/S3.** Do they also need d ≥ +0.005, as P1 does?
- **Human SE.** S2–S5 use jackknife only, which makes them more liberal than the yeast tests. Either inflate by the pre-fixed √(1 + 0.44²) or label them "gene-level uncertainty only".
- **Human t\*.** State the HepG2/Jurkat half-A convention (is t\* selected on half A only?).
- **Ribosome rule.** The "RPL/RPS ORFs" rule (line 182) needs a name → ORF step.

### S9. Verdict edge cases and multiplicity of individual claims. Lines 233–256.

- No verdict covers a failed synthetic R check (step 5). Say: fix it and log it before any real data; no verdict is issued.
- Verdict 4's "significant regulon loss" and verdict 5's claim on P1 alone are individual-endpoint claims outside the IUT. As in the CLR(log Pᵗ) preregistration, say that such claims use Holm over {P1, P2-NI, P2-harm}, or say explicitly that they are per-test-level statements.
- The design memo's condition (a), "run after the Collins readout", is absent from the DRAFT. Either add it as a protocol gate or log that the PI dropped it.

---

## NOTES

- **N1. Seed collision.** The F seed 20261104 equals bootstrap seed b = 1 (20261103 + 1). Harmless, but change one of them.
- **N2. Selection uncertainty.** t\* is fixed inside the bootstrap, so its selection uncertainty is not propagated. This is standard house practice; say so.
- **N3. t\* and N.** t\* is selected on halves of 350 profiles and applied at N = 700 (house practice; the flag rules exist). No truth leakage in depth selection found: criterion A uses only held-out CLR.
- **N4. S6 identifiability.** With the deleted gene masked, the deleteome's causal evidence for X → Y (column X) is invisible to the row-based edge (X, Y). The edge signal then comes only from co-variation of TF X's mRNA across other columns, and many yeast TFs are regulated post-translationally. S is symmetric, so direction is not identifiable. S6 is NI-only, which is acceptable, but state this.
- **N5. Strain artefacts.** Some deletion strains carry aneuploidies or secondary mutations, which give chromosome-wide shifts. These are shared by all operators and mostly hit co-located pairs. Optionally report a same-chromosome-arm sensitivity.
- **N6. H1 depth.** H1 uses t\*_K10, whereas the *E. coli* headline used the transferred t = 9. Pre-state which one is "the headline replication", and report the other descriptively (both are listed at line 145; the label matters).
- **N7. Separate family H.** A separate α for H1 is defensible, since it is a different claim. Say explicitly that the paper does not combine the H and P/S families.
- **N8. Absolute margins.** Absolute AUPR margins of 0.005 at small-complex prevalence near 0.1% mean very different relative effects across truths. Report relative changes and per-complex AP descriptively (already listed).
- **N9. Housekeeping.** The memo's compute total (180–240 core-hours) and the DRAFT's (170–220) differ, and the R check tolerance is 1e-10 in the memo vs 1e-8 in the DRAFT. Trivial.
- **N10. Implementation-check threshold.** The published M is a limma estimate with gene-specific dye correction, so r ≥ 0.99 and median |Δ| ≤ 0.02 against a simple mean may fail for benign reasons. Pre-specify what happens (see S2).

## What a Nature-family referee would attack (in order)

1. Truth–data circularity (B1): "held-out" YEASTRACT regulons built partly from the same deletion arrays. This would also gut H1.
2. Serial testing and post-hoc rule change: the third complex-operator test, C2 relaxed after seeing the tables, the variant added after one seed. Disclosed and declared the last test, which is good. Make sure the synthetic-to-real estimand gap (S1) is not a further, undisclosed degree of freedom.
3. "Why not CLR?" The synthetic data show CLR > K10 on small complexes, and R vs CLR is only descriptive.
4. One lab, one platform, one perturbation design, with the responsive-only selection: generality.
5. Absolute effect sizes of 0.005–0.010 AUPR, and power that depends on an imported SE ratio (S4).
6. The experiment resampling ignores the strong mutant-profile clustering (B3).
7. No DREAM5 yeast or standard TF–target benchmark; no comparison beyond CLR. Probably out of scope; anticipate it in Discussion.

## Ranked list

1. **B1** (BLOCKER): YEASTRACT expression evidence probably derived from the deleteome; P2/H1/S6 not held out.
2. **B3** (BLOCKER): SE_B uses exact duplicate columns and treats clustered mutants as i.i.d.
3. **B2** (BLOCKER): verdict 5 reads failed NI as a regulon cost.
4. **S1** (SHOULD-FIX, high): the synthetic estimand differs from P1; undisclosed CLR > K10 synthetic pattern and round-1 regulon caveat.
5. **S2** (SHOULD-FIX, high): F is uncalibrated, can pass trivially or fail spuriously, and its fallback is unimplemented.
6. **S3** (SHOULD-FIX, high): S1 belongs in a fixed sequence before Holm.
7. **S7** (SHOULD-FIX): gating and order gaps (selections and R check not enforced; test script after truth; evaluability counts; −1 parser; YEASTRACT timing).
8. **S6** (SHOULD-FIX): masking holes (unresolved labels, overlapping/neighbouring ORFs, aliases, paralogs).
9. **S5** (SHOULD-FIX): divergent-pair exclusion for regulons.
10. **S4** (SHOULD-FIX): power sensitivity to the SE_B ratio; P2/H1 projections.
11. **S8 and S9** (SHOULD-FIX): specification and verdict edge cases.
12. **N1–N10** (NOTES).

**Go/no-go:** do not freeze as is. After B1–B3 and S1–S3 are fixed, plus the code items in S7, the file is fit to commit as final. The remaining items can be fixed in the same pass, or logged in the change log before step 2.
