# CLR(log Pᵗ): interpretation of the preregistered verdict

*27 October 2026. Statistical interpretation memo (Claude, Opus, high effort) for the PI. It interprets the verdict computed by `analysis/clr_logpt_test.R` on 26 Oct and committed in `docs/clr_logpt_verdict/` (preregistration hash 9e985b3). No analysis was run and nothing was re-scored for this memo. The numbers come from the committed verdict CSVs and the truth-structure counts in `docs/clr_logpt_selections/ecoli_complexes_counts.csv`. The only new arithmetic is on those numbers: implied bootstrap SEs, power at hypothetical effects and minimum detectable effects. The memo does not edit the preregistration or the paper.*

## 0. Bottom line

1. **The verdict is NOT ADOPTED, and the preregistration fixes it.** One outcome-table row applies: "Complex superiority not shown on M3D or EcoMAC". Attention stays the method. The paper may report the log-domain readout only as (i) the explained mechanism of the K562 loss and (ii) a preregistered extension that was not confirmed.
2. **This is a failure to confirm. It is not evidence of no effect.** The two complex estimates point the predicted way:
   - M3D +0.026 AUPR (95% CI −0.010 to +0.061);
   - EcoMAC +0.029 (−0.006 to +0.063).

   Both are about five times the smallest effect of interest. But the test could reliably detect only effects of about 0.05 AUPR or more, and the intervals include both 0 and the K562 in-sample effect.
3. **The jackknife and the bootstrap disagree, and the disagreement has a simple explanation.** They measure different sources of variance:
   - the gene jackknife measures which genes and complexes are in the truth (SE 0.018);
   - the experiment bootstrap measures which experiments are in the compendium, with the truth fixed (implied SE 0.0075).

   The preregistration made the jackknife decisive and the bootstrap only an extra condition. The bootstrap therefore cannot override the jackknife, and the outcome table already treats the reverse disagreement as "not adopted".
4. **The power analysis had a gap.** The preregistration estimated power only for the regulon arms:
   - its regulon SE assumption (about 0.002) was exactly right (observed 0.0020 and 0.0018);
   - it gave no SE for the complex arms, and those turned out to be the binding constraint (observed 0.018 and 0.017, about nine times larger).

   This must be disclosed. It cannot be used to soften the verdict.
5. **A follow-up is justifiable only as a new preregistered study.** It would need a much larger complex truth that has never been touched, which in practice means another organism. It must report this E. coli test as the first, unconfirmed attempt. A second E. coli compendium would barely help, because the decisive variance comes from the truth, not from the expression data.

## 1. The result in one table

All values are pooled co-membership AUPR, new readout (a3r) minus comparator, at t = 9. The SE comes from the paired delete-half gene jackknife (100 replicates, seed 20261022). Intervals are d ± 1.96 SE. p is one-sided.

| Code | Dataset | Truth | vs | d | SE | 95% CI | p (raw) | Holm p (4 primary) | Rule met? |
|---|---|---|---|---|---|---|---|---|---|
| P1 | M3D | Complex Portal complexes, merged | attention | +0.0257 | 0.0182 | −0.0100, +0.0613 | 0.079 | 0.102 | no (needs p < 0.025) |
| P2 | EcoMAC | same | attention | +0.0285 | 0.0174 | −0.0057, +0.0627 | 0.051 | 0.102 | no |
| P3 | M3D | RegulonDB SC regulons | attention | +0.0043 | 0.0020 | +0.0004, +0.0082 | NI 1.6e-6 | 4.7e-6 | yes (NI, margin 0.005) |
| P4 | EcoMAC | same | attention | +0.0033 | 0.0018 | −0.0001, +0.0068 | NI 1.1e-6 | 4.5e-6 | yes (NI) |

- **Rule 3 (M3D cluster bootstrap, 200 draws):**
  - P1 interval [+0.013, +0.043], lower bound > 0;
  - P3 interval [+0.0033, +0.0068], lower bound > −0.005;
  - both conditions are satisfied.
- **Rule 4:** the K562 reproduction check passed.
- **Verdict (`verdict.csv`):** rule 1 FALSE, rules 2–4 TRUE, adopt FALSE.

## 2. (a) What the preregistration allows us to claim, and what it forbids

### 2.1 The applicable outcome row (verbatim)

> | Complex superiority not shown on M3D or EcoMAC | Not adopted. The K562 rescue is reported as exploratory and not replicated in E. coli. If attention does not lose to CLR on E. coli complexes (descriptive), this is said too (no loss to rescue), but it does not change the verdict. |

The other rows do not apply:

- "Complexes improve, regulon NI not shown": NI was shown.
- "Regulon loss": there was none; both regulon estimates are positive.
- "Gene jackknife passes, cluster bootstrap fails": the reverse happened.
- "Implementation check fails": it passed.

Two further binding sentences:

> If no, the readout is reported as the explained mechanism of the K562 loss. (Question)

> In every case the paper keeps the measured mechanism (softmax budget turns strength into share; the walk repairs the within-row order, not the scale) as an explanation of the K562 loss. (Decision rules)

### 2.2 The "does not count" list (verbatim) and how it applies here

> What does **not** count: a win on K562 alone; a win by a2r or a1 when a3r fails; wins at any depth other than t = 9 (including the committed t\*, which is secondary); wins on a sensitivity truth when the primary truth fails; a win on per-complex AP when pooled fails; a win after swapping the complex source, size rule, merge rule or operon rule; results within ±0.005.

| Item | How it applies |
|---|---|
| a win on K562 alone | K562 is in-sample (the readout was found there among about 20 operators and readouts). Whatever the K562 merged-CORUM number is, it cannot support a complex claim. |
| a2r or a1 when a3r fails | a3r failed rule 1. No a2r result may be presented as "the readout works". |
| depths other than t = 9, including t\* | The held-out t\* results (M3D 6, EcoMAC 5, RPE1 12) are descriptive. A pass at t\* is not a pass. |
| sensitivity truth when the primary fails | The operon-included, size ≥ 2, no-merge and ribosome-excluded truths, the r-protein module coherence and raw CORUM may be shown only as descriptive. None replaces P1 or P2. |
| per-complex AP when pooled fails | Per-complex (macro) AP is descriptive. It cannot replace pooled AUPR. |
| swapping the complex source | **Key for any follow-up:** EcoCyc was the intended source, and Complex Portal was the logged fallback. An EcoCyc result in the same compendia can never stand in for P1 or P2 (section 7). |
| results within ±0.005 | Not relevant to P1 and P2: both are outside ±0.005 but not significant. "Above the minimum effect" does not mean "shown". |

### 2.3 What may be claimed

1. **The readout was not adopted.** Its complex superiority over attention was not shown on either E. coli compendium at the preregistered threshold. Use the preregistration's wording: the K562 rescue is exploratory and "not replicated in E. coli".
2. **Regulon non-inferiority.** The readout is non-inferior to attention on strong-evidence regulons in both compendia (margin 0.005; Holm-adjusted one-sided p < 10⁻⁵).
   - This is a legitimate single-endpoint claim: the preregistration allows single primary-endpoint claims with Holm over the four primary tests.
   - It is a claim about a readout the paper does not adopt. Its value is that it shows the readout has no regulon cost. It must not be turned into a claim that the readout improves regulons:
     - superiority on regulons was not a hypothesis;
     - the P4 interval touches 0.
3. **No detectable loss to rescue on E. coli complexes (required by the outcome row).** Descriptively, attention minus CLR on the merged complex truth is:
   - M3D +0.006 (95% CI −0.021 to +0.032);
   - EcoMAC −0.002 (−0.028 to +0.023).

   The intervals are wide, so this is "no detectable loss", not "attention ties CLR". On K562 (in-sample, raw CORUM), by contrast, CLR beat attention by +0.038 [+0.019, +0.057].
4. **The mechanism.** The paper keeps the measured mechanism as the explanation of the K562 loss: the softmax budget turns strength into share, and the walk repairs the within-row order but not the scale. It must say that this mechanism was measured in exploratory, in-sample K562 analyses.
5. **The estimates themselves.** Point estimates and intervals for P1 and P2 are reported in full. The precision of the test may be reported as a limitation (section 4). It may not be reported as a reason to believe the effect.

### 2.4 What is forbidden

- Any wording that calls the complex gain established, significant, marginal, a trend or "nearly significant". Report P2 as p = 0.051, never "p = 0.05".
- Any claim of a key-secondary result (S1–S5) as tested or significant. The Holm family of S1–S5 is gated on adoption, and the gate did not open. The `p_raw` column in `secondary.csv` must not be quoted (section 6).
- Promoting the bootstrap interval over the jackknife (section 3).
- Combining P1 and P2 after the fact. A pooled z of the two would cross 1.96 only if the two compendia were correlated by less than about 0.2. They share an organism and a truth, and may share arrays, so that assumption is not credible. The pooled test was not preregistered, and it must not be computed for the paper or reported.
- Rerunning the jackknife with more replicates or another seed ("one seed, never rerun with another"), or replacing it with a complex-level resampling after the fact.
- Counting K562, RPE1, M3D and EcoMAC as "four of four datasets point the same way". K562 is in-sample, RPE1 is gated and secondary, and the two E. coli compendia are not independent.

## 3. (b) Why the jackknife and the bootstrap disagree

### 3.1 What each one estimates

| | Paired delete-half gene jackknife (decisive) | Paired cluster bootstrap over experiments (rule 3) |
|---|---|---|
| What is resampled | Truth genes: half of the 561 (M3D) or 544 (EcoMAC) complex genes kept per replicate; pairs with both genes kept are scored | The 466 M3D experiments, drawn with replacement, one replicate array per occurrence; MI, CLR, attention and a3r are refitted in every draw |
| What is held fixed | The fitted networks (expression data) | The truth and the gene universe |
| Source of variance | Which genes, and hence which complexes and pairs, make up the benchmark: variation from sampling the truth | Which experiments make up the compendium: variation from sampling the conditions |
| Question it answers | Would the gain hold on another set of complexes like these? | Would the gain on these complexes hold on another set of experiments like these? |
| Availability | All four primary endpoints | M3D only (P1, P3) |
| SE of the P1 difference | 0.0182 | about 0.0075 (implied: 95% percentile width / 3.92) |
| SE of the P3 difference | 0.0020 | about 0.0009 |

The bootstrap SE is 0.41 times the jackknife SE for P1 and 0.44 times for P3. The ratio is about the same on regulons and on complexes. So the gap is a general property of the two procedures on these data, not something special to the complex endpoint.

- On regulons both procedures agree, because the effect is large relative to either SE.
- On complexes the estimate, +0.026, lies between about 1.4 jackknife SEs and about 3.4 bootstrap SEs. Exactly there, the procedures must disagree.

### 3.2 Why the complex jackknife SE is large

The E. coli complex truth is small and concentrated:

- 54 (M3D) and 53 (EcoMAC) merged complexes contribute positive pairs;
- the ribosome-excluded variant keeps 451 of 1,150 positives (M3D) and 449 of 1,148 (EcoMAC). So about 60% of positive pairs sit in r-protein-dominated merged complexes, while the 54 r-protein genes are about 10% of the truth genes.

A delete-half replicate keeps a varying number of r-protein genes (about 27 ± 3.5), and the number of ribosomal positive pairs grows with the square of that number. So each replicate moves a large, correlated block of positives in or out. That is real uncertainty about generalizing across complexes, and the gene jackknife is designed to capture it. The bootstrap holds this block fixed and cannot see it.

One technical point, for the record only. For a pair-level statistic, a delete-half gene jackknife reproduces the gene-level variance component correctly. Heuristically (treating AUPR as a degree-2 U-statistic), it overstates the purely pair-level component by about a factor of 3, so it may be somewhat conservative. Two limits apply:

- AUPR is not a U-statistic, and the size of this effect here is unknown;
- the gene-level (r-protein block) component is plausibly the dominant one.

This is a design question for a future preregistration (section 7). It is not grounds to adjust the present SE.

### 3.3 Why the preregistration made the jackknife decisive

1. **It is the paper's standard.** Every headline AUPR interval in the paper uses the paired delete-half gene jackknife.
2. **It exists for all four primary endpoints.** EcoMAC has no experiment bootstrap, so a bootstrap-based rule could not have tested P2 or P4.
3. **It targets the right variance for a complex claim.** The claim is about protein complexes in general. Uncertainty about which complexes are in the benchmark is the uncertainty that limits it.
4. **The bootstrap was written in as an additional necessary condition (rule 3), not as an alternative.** Adoption needed both. The outcome table spells out the reverse case (jackknife passes, bootstrap fails, so not adopted). The logic is the same in both directions: two necessary conditions, and failing either one means "not adopted".

### 3.4 What may and may not be done with the disagreement

- **May:** report the bootstrap interval for P1 and P3 in the supplement, labelled as experiment-level robustness with the truth fixed. Say that rule 3 was met, and explain in one sentence why the two intervals differ (different sources of variance).
- **Must say:** the two kinds of uncertainty add up; they are not alternatives. An interval that accounted for both would be wider than the jackknife interval: roughly √(0.018² + 0.0075²) ≈ 0.020 under independence. It would not be narrower.
- **May not:**
  - describe the bootstrap as showing the complex gain "is robust" or "is significant";
  - compute or quote a p-value from the bootstrap;
  - choose the smaller SE;
  - describe the result as "significant under one resampling scheme".
- **May not:** treat the disagreement as a reason to rerun, reweight or replace the jackknife.

## 4. (c) Power: what was assumed and what was observed

| Arm | SE assumed in the preregistration | Observed jackknife SE | Comment |
|---|---|---|---|
| Regulons, M3D (P3) | about 0.002 ("likely", from K562 TF truths at about 60% of the M3D attention–CLR SE of 0.003) | 0.0020 | exactly as assumed |
| Regulons, EcoMAC (P4) | about 0.002 | 0.0018 | as assumed |
| Complexes, M3D (P1) | **none given** | 0.0182 | about 9 times the regulon SE |
| Complexes, EcoMAC (P2) | **none given** | 0.0174 | about 9 times |

The preregistration's power section analysed only the regulon arms, and it named them "the main risk of a false 'not confirmed'". At a true difference of 0, it put NI power at about 0.71 per compendium and the joint power at 0.5–0.71. As it turned out:

- the regulon estimates were positive (+0.003 to +0.004), so NI passed overwhelmingly;
- the binding constraint was the complex arm, for which no SE was projected.

Information to make that projection existed beforehand:

- the K562 exploratory interval for a3r vs attention on CORUM, +0.059 [+0.033, +0.084], implies an SE of about 0.013;
- that SE came from a much larger complex truth, so an E. coli truth of about 54 contributing complexes was bound to give a larger one;
- the E. coli truth counts were logged on 26 Oct, before scoring, and showed this.

This is a design gap, and it should be disclosed in the supplement.

**Sensitivity of the test as run** (one-sided α = 0.025, jackknife SE 0.018; these are design-sensitivity figures, not "observed power"):

| True difference | Power per compendium | Both compendia, if independent | Both, if fully correlated |
|---|---|---|---|
| +0.005 (smallest effect of interest) | 0.05 | 0.002 | 0.05 |
| +0.010 | 0.08 | 0.01 | 0.08 |
| +0.027 (about the observed estimates) | 0.29–0.37 | 0.09–0.14 | 0.29–0.37 |
| +0.035 | 0.49–0.52 | 0.24–0.27 | 0.49–0.52 |
| +0.050 | 0.78–0.82 | 0.62–0.67 | 0.78–0.82 |
| +0.059 (K562 in-sample a3r vs attention, raw CORUM) | 0.90–0.92 | 0.81–0.85 | 0.90–0.92 |

- **Minimum detectable effect** at 80% power per endpoint: about 0.05 AUPR.
- **SE needed** to detect the preregistered minimum effect of 0.005 at 80% power: about 0.0018, ten times smaller than achieved.

How to read this:

1. **The test was insensitive to anything but a large effect.** The design had a reasonable chance to confirm only if the E. coli effect were about as large as the K562 in-sample effect. In-sample effects are optimistic: the readout was chosen among about 20 candidates on the same K562 truths. On top of that, E. coli showed no attention–CLR complex loss to rescue. A smaller E. coli effect was therefore foreseeable.
2. **Low power does not turn a non-significant result into support.** "Underpowered, therefore probably real" is not a valid inference.
   - The data are compatible with no effect and with an effect as large as K562's.
   - Taken alone, the point estimates carry the usual selection-free interpretation: they are unbiased-ish estimates of the E. coli effect. They are not evidence that the effect exceeds 0.
3. **Do not report "observed power" or "post hoc power".** It is a monotone transform of the p-value and adds nothing. Report the minimum detectable effect and the SE instead, as above.

## 5. (d) Recommended wording for the paper

### 5.1 Main text (Results, end of the "Tight cliques" paragraph; 2–3 sentences)

> The loss has a measurable cause: the softmax budget turns each gene's dependence strength into a share, so genes with many interchangeable partners, such as ribosomal proteins, receive diluted weights, and the walk repairs the order within each row but not the scale across rows. A readout that re-standardizes each row of the log walk against its own background (CLR on log Pᵗ), found in exploratory analyses of these K562 data, was then tested under a preregistration on *E. coli* protein complexes and was not confirmed: it exceeded attention by +0.026 AUPR on M3D (95% CI −0.010 to +0.061; one-sided P = 0.079) and +0.029 on EcoMAC (−0.006 to +0.063; P = 0.051), short of the preregistered threshold of P < 0.025 in both, while remaining non-inferior on regulons. On these *E. coli* complexes attention showed no detectable loss to CLR, so attention remains the method and the K562 rescue remains exploratory (Supplementary Note 15).

If space is tight, drop the last sentence's first clause. Keep "not confirmed", both intervals and "remains exploratory".

**Abstract:** no change, and no mention of the readout. The existing sentence on tight co-expression cliques stays true.

**Discussion:** in the "reach has a price" paragraph, one clause:

> ... and it blurs tight cliques such as ribosomal complexes; a log-domain re-normalization that repairs this in K562 was not confirmed in a preregistered test on *E. coli*.

### 5.2 Methods, "Preregistered analyses" (one sentence to add)

> The log-domain readout (`docs/preregistration_clr_logpt.md`, fixed 2026-10-22 after it had been found on K562 among about 20 exploratory readouts) was tested at t = 9 against attention on *E. coli* protein complexes from the EBI Complex Portal (superiority; the preregistered fallback, because EcoCyc files require a licence) and on strong-evidence regulons (non-inferiority, margin 0.005), on M3D and EcoMAC, with adoption only if all four tests passed at one-sided α = 0.025 and an experiment-level cluster bootstrap agreed.

### 5.3 Supplement (new paragraph in Supplementary Note 15, after "Learned attention")

> **Log-domain readout (`docs/preregistration_clr_logpt.md`; fixed 2026-10-22, before the readout was computed on any *E. coli* or RPE1 data and before any *E. coli* complex truth was built).** Exploratory analyses of the K562 CORUM loss (Supplementary Note 11) found that re-standardizing each row of log Pᵗ over reachable pairs and combining the two directions with CLR's Stouffer rule (a3r) recovered complexes without a regulon cost on K562; that readout was chosen among about 20 operators and readouts on the same K562 truths, so K562 is in-sample. The confirmatory test, at the method's transferred depth t = 9, required all of: superiority over attention on merged *E. coli* protein complexes on M3D and EcoMAC (one-sided P < 0.025 and an estimate of at least +0.005 AUPR); non-inferiority on strong-evidence RegulonDB regulons on both (margin 0.005); and, on M3D, a paired cluster bootstrap over experiments (200 draws) with a complex interval above 0 and a regulon interval above −0.005. Inference used the paired delete-half gene jackknife (100 replicates). The complex truth is the EBI Complex Portal *E. coli* K-12 set (taxon 83333; 324 complexes; 2026-01-14 release), the preregistered fallback, used because EcoCyc's flat files require a licence; after flattening, the size filter (3–200 genes) and the overlap merge it gives 54 (M3D) and 53 (EcoMAC) complexes with positive pairs and 1,150 and 1,148 positive pairs (same-operon pairs excluded), about 60% of them in ribosome-dominated complexes. **The readout was not adopted.** Its complex gains over attention were +0.026 on M3D (95% CI −0.010 to +0.061; one-sided P = 0.079) and +0.029 on EcoMAC (−0.006 to +0.063; P = 0.051); both point in the predicted direction but neither reaches the threshold. It was non-inferior on regulons (M3D +0.0043, 95% CI +0.0004 to +0.0082; EcoMAC +0.0033, −0.0001 to +0.0068; Holm-adjusted one-sided P < 10⁻⁵ for both). The experiment bootstrap met its condition (complexes +0.013 to +0.043; regulons +0.003 to +0.007), but it resamples experiments with the truth fixed and so does not capture the variation over complexes that the gene jackknife measures, which is the larger of the two here (SE about 0.0075 against 0.018); it was a necessary, not a sufficient, condition. Descriptively, attention did not lose to CLR on these complexes (M3D +0.006, −0.021 to +0.032; EcoMAC −0.002, −0.028 to +0.023), unlike on K562 CORUM, so the *E. coli* truth offered no attention loss to rescue. The preregistered power analysis covered only the regulon arms, whose standard error it predicted correctly (0.002); the complex arms had a standard error near 0.018, so the test had 80% power only for effects of about 0.05 AUPR or more, and the result is a failure to confirm, not evidence that the readout does not help. Key secondary contrasts (against CLR, and on RPE1) were gated on adoption and were therefore not tested; they are given as intervals only in Supplementary Table [X]. The K562 rescue remains exploratory and is not replicated in *E. coli*.

### 5.4 Supplementary Note 11, "Candidate fix"

Replace the paragraph with:

> **Mechanism and candidate fixes.** Exploratory K562 analyses locate the loss in the softmax budget, which turns dependence strength into a share, so genes with many interchangeable partners (ribosomal proteins in K562 have about 77) receive diluted weights; the walk repairs the within-row order but not the scale across rows. Two preregistered candidates target this: a log-domain readout, which was not confirmed on *E. coli* complexes (Supplementary Note 15), and the multi-head operator [PENDING].

Every number in this text must come from the measured mechanism in `docs/plans/corum_analysis.md`, labelled exploratory. The current sentence ("consistent with the walk trading their sharp pairwise dependences for aggregation") should be replaced by the measured mechanism, as section 4 of that memo recommends. The same applies to the main-text "Tight cliques" paragraph.

### 5.5 Limitations (one or two sentences to add)

> A log-domain readout that explains and repairs the K562 complex loss in exploratory analyses was not confirmed in a preregistered test on *E. coli* complexes; that test had limited precision, because the *E. coli* complex truth has only 53–54 complexes with positive pairs, most of them in the ribosome, and effects below about 0.05 AUPR had less than 80% power, a constraint the preregistered power analysis did not anticipate. Whether the readout improves complex recovery beyond K562 is therefore open.

### 5.6 Wording to avoid

| Avoid | Why | Use instead |
|---|---|---|
| "marginally significant", "trend toward", "nearly significant", "borderline", "p ≈ 0.05" | P2 = 0.051 and P1 = 0.079 did not meet a preregistered α of 0.025 | "not confirmed"; give the CI and the exact P |
| "improves complex recovery in *E. coli*", "rescues complexes" | P1 and P2 not shown | "exceeded attention by +0.026 (CI ...), not confirmed" |
| "robust to resampling experiments" (about complexes) | The bootstrap holds the truth fixed and was a necessary condition only | "met the experiment-bootstrap condition, which does not test variation over complexes" |
| "replicates in human RPE1 cells", "confirmed in a second cell line" | S5 is gated and untested; RPE1 shares the CORUM ribosome content | descriptive interval only (section 6) |
| "beats CLR on complexes" | S3 and S4 are gated and untested | none (or descriptive interval only, in the supplement) |
| "improves regulons", "no regulon cost" without qualification | P3 and P4 were NI tests; P4's CI touches 0 | "non-inferior on regulons (margin 0.005)" |
| "the effect is real but underpowered", "likely a true effect" | A non-significant result is not evidence for the effect | "the test could detect only effects of about 0.05 or more; the question is open" |
| "shown not to work", "no effect", "refuted" | Absence of evidence; both estimates are positive with wide CIs | "not confirmed" / "not replicated in *E. coli*" (the preregistration's own wording) |
| "attention ties CLR on *E. coli* complexes" | Equivalence was not tested; the CIs span ±0.025 | "no detectable loss" |
| "consistent across four datasets" | K562 is in-sample, RPE1 is gated, and the two *E. coli* compendia are not independent | none |
| any statement about t\* results, a2r, a1, per-complex AP or sensitivity truths as support | On the "does not count" list | descriptive, labelled as such, in the supplement only |
| "CLR → attention → CLR(log)" as the method's name or framing | That was conditional on adoption (PI decision 4) | attention stays the method; the readout is "a log-domain re-normalization (not adopted)" |

## 6. (e) How K562 and RPE1 may be mentioned

The key-secondary family (S1–S5) was to be tested "only if adoption holds (gatekeeping, so no α leaks from the primary family)". Adoption did not hold. So S1–S5 were not tested (`secondary.csv`: `interpretable = FALSE`, no Holm p). They may appear only as descriptive intervals. **Never quote their `p_raw` values.** Reporting a p-value for an untested gated contrast would let α leak in practice even though it does not on paper.

| Code | Contrast (descriptive) | d (95% CI) | How to mention |
|---|---|---|---|
| S1 | M3D regulons, a3r vs CLR | +0.026 (+0.020, +0.033) | Not needed. It mostly restates the headline (attention vs CLR) plus P3. If shown, add "descriptive; the headline over CLR stays attributed to attention". |
| S2 | EcoMAC regulons, a3r vs CLR | +0.031 (+0.024, +0.038) | as S1 |
| S3 | M3D complexes, a3r vs CLR | +0.032 (+0.006, +0.057) | Supplementary table only. **Trap:** this CI excludes 0 while P1's does not, because attention ≈ CLR here and the CLR difference has a smaller SE (0.013). It is gated and untested and is no basis for "beats CLR on complexes". |
| S4 | EcoMAC complexes, a3r vs CLR | +0.026 (+0.001, +0.051) | as S3 |
| S5 | RPE1 merged CORUM, a3r vs attention | +0.036 (+0.013, +0.059) | Supplement, one sentence (below). Never "replicated". |

**RPE1, suggested sentence (supplement only):**

> In RPE1 cells (the key secondary set, not tested because the adoption gate was not passed) the difference on overlap-merged CORUM was +0.036 AUPR (95% CI +0.013 to +0.059, descriptive); RPE1 comes from the same Perturb-seq study as K562 and its CORUM content is dominated by the same ribosomal complexes, so it is not an independent test of the complex claim.

Facts behind the caveats:

- the preregistration itself calls RPE1 "a human replication in another cell line, but with the same CORUM content (the ribosome again)";
- the drafter of the preregistration had seen some RPE1 |r| and MI rows (disclosure table).

**K562.** K562 is in-sample: the readout was found on it. The outcome row requires "the K562 rescue is reported as exploratory".

- Use the K562 merged-CORUM number from the confirmatory run (`results/clr_logpt_test/truth/k562ess/metrics.csv`), labelled "in-sample, exploratory". Fact-check it before use; this memo did not read it.
- Do not use the exploratory memo's raw-CORUM figures (7.95×) as if they came from the merged truth.
- The preregistration's rule that a K562 merged-CORUM difference ≤ 0 must be "reported prominently" is formally an after-adoption rule. If that happened, report it anyway.
- The K562 t\* note in the change log applies to any t\* report (t\* chosen on the full half-A gene set, scored on the 5,671-gene subset).

**Not allowed anywhere:**

- a tally of K562 + RPE1 + E. coli as converging evidence;
- RPE1 or K562 in the abstract or the main-text complex claim;
- RPE1 or K562 at t\*, or with a2r/a1, as support.

## 7. (f) Is a follow-up justified, and how would it have to be preregistered?

### 7.1 Is it justified?

On the science, a follow-up is defensible but not urgent.

**For it:**

- both confirmatory estimates point the predicted way and are about five times the smallest effect of interest;
- the test was insensitive below about 0.05;
- the regulon side is clean, so a readout with no regulon cost and a possible complex gain is worth knowing about.

**Against urgency:**

- the paper does not need it: attention stays the method, and the readout's role (the mechanism) is secured either way;
- on E. coli there was no attention–CLR complex loss to rescue, so the readout's motivating premise is so far specific to K562 CORUM, which means the ribosome;
- a follow-up raises a real serial-testing risk: testing until it passes.

**Recommendation:** do not delay the paper. Submit with "not adopted" as written above. If the PI wants the complex question answered, run one new, separately preregistered test, and declare in advance that it is the last confirmatory attempt for this readout.

### 7.2 What cannot serve as a follow-up

| Proposal | Why not |
|---|---|
| Re-test on EcoCyc complexes in M3D or EcoMAC once licensed | Explicitly does not count ("swapping the complex source"; "the unused source may be reported only as a descriptive sensitivity analysis; it never replaces the chosen one"). It also reuses the same expression data and heavily overlapping complexes (the ribosome again), so it is not independent. It may be a descriptive sensitivity analysis only. |
| A second *E. coli* expression compendium with the same complex truth | It barely addresses the problem. The decisive variance is over the truth (which complexes), not over experiments: the bootstrap SE is about 0.4 times the jackknife SE. With the same 53–54 complexes, the SE would stay near 0.017. |
| Re-analysing the current data: complex-level resampling, more jackknife replicates, the RP-excluded or per-complex endpoint, t\*, a2r, pooling P1 and P2 | All on the "does not count" list or forbidden by the inference section; at best data-dependent re-analyses of a test already run |
| Promoting RPE1 (S5) | Gated; same CORUM ribosome; same study as the discovery data |

### 7.3 What a legitimate follow-up would need (a new file, e.g. `docs/preregistration_clr_logpt_2.md`; this file stays unchanged)

1. **New data and a much larger complex truth, never scored with a3r.** In practice this means another organism with many curated complexes and a regulon truth for the no-harm arm. Two examples to check for availability and prior touches:
   - a *S. cerevisiae* expression compendium with Complex Portal or CYC2008 complexes and YEASTRACT regulons;
   - a human compendium other than Replogle et al. with a complex source and an overlap rule that cap the ribosome's share.

   K562 genome-wide would not qualify: it is the same cells, the same screen family and the same CORUM content.
2. **A power analysis for the complex arm, done before scoring.**
   - Project the SE from truth structure alone: the number of complexes contributing positives, and the largest complex's share of positives. Anchor the projection on the E. coli SE of 0.018 at 54 complexes with about 60% of positives in one module.
   - Target effect: no larger than the E. coli estimate (about +0.027). The target must not be the K562 in-sample +0.059.
   - Detecting +0.027 at 80% power needs an SE ≤ 0.009, roughly four times the information of the E. coli truth.
   - If the projected SE does not reach the target, the test is not run as confirmatory.
3. **Fix the inferential unit before scoring, matched to the claim.** The claim is about complexes. Options:
   - keep the gene jackknife for consistency with the paper;
   - or make a complex-level resampling decisive, after checking its calibration on synthetic or permuted truths without labels.

   Either way, name one decisive procedure and state whether any second procedure is a necessary condition (as rule 3 here) or descriptive.
4. **Decide the complex endpoint a priori.** Options: pooled AUPR on the merged truth, a ribosome-capped truth, or per-complex AP. Each needs a written rationale. Disclose that the choice was made after this E. coli failure.
5. **Keep the regulon NI arm** with the same margin of 0.005. The preregistration's SE assumption of 0.002 proved accurate.
6. **Report this E. coli test as the first, unconfirmed attempt in every write-up of the follow-up.**
   - The follow-up cannot change this verdict retroactively.
   - Adoption must rest on the new test passing on its own.
   - Any combined analysis with E. coli must be specified in advance, including how the dependence between the studies is modelled.
7. **Stopping rule.** State that this is the last confirmatory test of a3r. If it fails, the readout stays an explanatory ablation for good.
8. **The asymmetry follow-up (C3, `docs/plans/asymmetry_expert.md`).** Its plan includes "a3r if adopted". a3r was not adopted, so in that preregistration it may be at most a descriptive comparator.

## 8. Checklist for the write-up

- [ ] Main text, Methods, Supplementary Notes 11 and 15 and Limitations use the wording in section 5. The numbers are fact-checked against `docs/clr_logpt_verdict/*.csv` (worker agent).
- [ ] No p-value for S1–S5 appears anywhere.
- [ ] The K562 merged-CORUM number is taken from the confirmatory run's in-sample output and labelled exploratory.
- [ ] Every mention discloses the complex source fallback (Complex Portal, not EcoCyc) and that the readout was found on K562 among about 20 candidates.
- [ ] The "Tight cliques" paragraph and Supplementary Note 11 are rewritten to the measured mechanism (`corum_analysis.md` section 4, item 4).
- [ ] The power-analysis gap (complex arms not covered) is disclosed in the supplement.
