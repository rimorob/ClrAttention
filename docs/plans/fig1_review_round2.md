# Fig. 1 (rebuilt): fresh-eyes review, round 2

Audience: Jim, a senior collaborator who is (a) a big-picture person and (b) a high-impact person; the figure must show the potential of the idea in about 20 seconds at 7.2 in width without claiming more than the data support.

Reviewed: `docs/whitepapers/figures/fig1_method.png/.pdf` (2160 x 2799 px = 7.2 x 9.33 in at 300 dpi; commit f9787a1), `analysis/fig_method.R`, `docs/whitepapers/paper.md` (commit 1bed906: Methods, Results paragraph, legend), `src/core/clr_core.cpp`, `analysis/baselines.R`, `R/calibrate.R`.
Reviewer caveat: I did not build the figure, but I am the same model family as its builder and as the round-1 reviewer, so this is an independent re-derivation, not a human fresh-eyes review. Nothing was viewed in print. Review only: no figure, code or paper edits. Scratch scripts are in the session scratchpad, not in the repo.

## 0. First impression, written before reading the round-1 review (image alone, 7.2 in width)

What I took away in 20 seconds:
- A black badge says "worked example (toy): not results". Title: CLR as attention, compendium in, module map out, no learned parameters.
- Strip: expression matrix goes through four boxes (Dependence, CLR, Softmax, Walk) to a pink gene x gene map with two outlined blocks. The coloured step bars (blue, orange, violet, pink) recur in the four panels below.
- Panels read 1 to 4: a blue matrix, bell curves with z-scores and a formula, an orange matrix turning into a violet one, three pink matrices (t = 1, 3, 9) that get blockier. A bottom-right strip of bars shows 0.69, 0.77, 0.91.
- What I would tell a colleague: "CLR is the first two steps of an attention mechanism; they add softmax and a multi-step walk, with no trained weights, and on a toy the blocks come out cleaner with depth." Input and output are clear. "Added here" is visible only as 7 pt grey italic in the panel corners.

What was confusing or weak:
- Why it might matter is not on the figure. The only hint is three small bars at the bottom right (AP 0.69, 0.77, 0.91), under a 7 pt bold grey caption, after two minutes of reading.
- The top strip does not separate published from added. The four boxes look equal.
- Panel 2 is the densest block (two curves, two z formulas, S formula, matrix, arrows). It is correct and now readable, but it takes about 15 s alone.
- "Output: module map" with outlined blocks looks like the method detected the modules. The outlines are the planted truth.
- Panel 4 bottom is a stack of five 7 pt grey sentences; at 20 s I read none of them.
- Greyscale: the pipeline still reads; the colour-named annotations ("orange tick", "magenta frame") do not.

Then compared with round 1 (`docs/plans/fig1_fresh_eyes_review.md`): all of F1-F13, F15-F18 are addressed (diagonal wording, shared MI axis, mu/sigma, k marked, t = 9 framed, identity colours replaced, assertions extended, toy badge, strip, AP inset). F14 (published vs added bracket) is only half done (see R1 below) and F19 ("potential") is not delivered beyond the swap-in pill.

## 1. Twenty-second test and overclaim check

| Question | Answer |
|---|---|
| Input | Yes (strip, left; the 14 x 60 toy heat strip, "thousands of genes, hundreds+ experiments, many labs"). |
| Output | Yes, but see R2: the outlines are planted truth, not method output. |
| What softmax + walk add | Partly. Panels 3 and 4 say so in small italics ("added here"); the strip does not. |
| Why it might matter (potential) | Weak. One illustrative AP inset (toy) low in panel 4; no pointer to the real-data result; "swap-in slot" pill is the only generality cue. |
| Toy vs result | Good: badge top right, "toy" in the strip, inset and footer. A reader cannot mistake 0.31, 2.8, 0.078 for results. |
| Implies it beats learned attention | No. The legend and footer say "analogy, not equivalence". Good. |
| Implies uniform wins | Risk, see R3: the toy shows the walk helping co-membership AP everywhere, whereas the paper reports losses on tight complexes (ribosome, CORUM) and gains concentrated in loose regulons. |
| Implies a hierarchy of module types | Risk, see R2: "regulons, complexes, operons (generic types)" next to the output; paper Fig. 6 / Supp. Note 12 report no fine-to-coarse ordering, and "no depth helps CORUM". |

## 2. Findings

Severity key: wrong = factually incorrect; misleading = a careful reader draws a false conclusion; unclear = a reader cannot tell; cosmetic. No finding is rated wrong.

### Claims and content

| # | Sev | Location | Finding | Proposed fix |
|---|---|---|---|---|
| R1 | unclear | strip and panel tags (`fig_method.R` L260-285; unused `by_ <- Y(0.7)` at L282) | The "published CLR (2007) / added here" bracket described in the f9787a1 commit message is not drawn; the strip has four equal boxes. Published-vs-added survives only as 7 pt grey italic "CLR (2007)" / "added here" at the right edge of each panel title, where a fast reader does not look. This was round-1 F14 and the single most useful cue for Jim ("what is new"). | Draw the bracket (or a two-tone tint under the four boxes): "published CLR" over boxes 1-2, "added here" over 3-4, in the strip. Delete `by_` if abandoned. |
| R2 | misleading | strip output, `fig_method.R` L287-296 | (a) "outlined: modules, e.g. regulons, complexes, operons (generic types)" with a map whose outlines are the planted blocks (drawn from `chain`, `modB`). Reads as "the method outputs modules and these are the module types it finds". The method outputs a gene x gene map; no module calling is drawn. (b) The paper's own results contradict the generic list: attention loses to CLR on tight complexes (ribosomal, K562), "no depth helps CORUM", depth is not a resolution parameter (Fig. 6, Supp. Note 12), gains concentrate in large loose regulons. | Caption: "gene x gene map; outlines = planted modules in the toy (not detected)". Replace the six-line wrap by "modules, e.g. regulons". Drop "complexes, operons". |
| R3 | misleading | AP inset, panel 4 (L456-475) | Honest as to labels ("one seed, planted modules; illustrative"), and honest in showing K_1 (0.69) below CLR (0.77) and K_3 (0.77) equal to it. But (i) the aggregate gain at t = 9 comes from the tight 4-gene one-factor module, not from the chain where gene k lives: in scratch, per-module AP is chain 0.49 (CLR) -> 0.47 (K_9), tight module 0.38 -> 0.83. So the mechanism told in panels 3-4 (k reached through the walk) is not what produces the bar; and in real data the paper reports tight complexes as the place CLR wins. (ii) The toy is built to favour the walk (AR(1) chain plus factor module, alpha = 0.5), and `stopifnot(AP_K9 > AP_CLR, ...)` (L117) makes the figure conditional on a favourable outcome. (iii) Robustness is fine (scratch: K_9 > CLR in 20 of 20 seeds, mean AP CLR 0.77, K_1 0.65, K_3 0.74, K_9 0.89), so the seed is not cherry-picked, but that is still planted structure. | Keep the inset (it is the only "benefit" cue) but caption: "toy, planted modules, favourable by construction; real-data effect sizes and the cases where it loses: Figs. 2-6". Add a pointer "real data: Fig. 3". Do not print the 0.91 any larger. Optionally state the per-module split in the legend. These per-module and 20-seed numbers are mine (scratch), not in the repo and not paper claims. |
| R4 | misleading | panel tags "CLR (2007)" over panels 1 and 2; legend "the first two published CLR (Faith et al. 2007)" | Step 2 as drawn is not the 2007 rule: S = (z+_ij + z+_ji)/sqrt(2) is the Stouffer combination; `R/calibrate.R` L14 labels sqrt(z_ij^2 + z_ji^2) as "historical 2007", and Methods says "the original Euclidean rule". Step 1 uses a B-spline MI (Daub 2004, common-bin rule), which Methods contrasts with per-gene adaptive bins; I recall (from memory, not checked) that Faith 2007 used a kernel-density MI estimator, so the estimator also differs. The idea (per-row z-score of MI, clipped) is 2007. | Tag "after CLR (2007)" or "CLR (2007), modified: B-spline MI, Stouffer"; legend: "the first two follow published CLR (Faith et al. 2007; here with B-spline MI and Stouffer combination)". |
| R5 | misleading | pill "swap-in slot: any dependence measure (MI, rank correlation)" (L324) | The paper tests MI, |r| and normal-scores correlation (CLR-NS), found the choice matters little on E. coli, and that MI is poor for zero-inflated counts. "Any" is untested; "rank correlation" is shorthand for Blom normal scores, not Spearman. | "swap-in slot: e.g. MI, correlation (tested)". |
| R6 | unclear | panel 4 tag "(cf. attention rollout)"; legend same | No reference in the paper (reference list L232 has none: no Abnar & Zuidema 2020). "Rollout" multiplies layer-specific attention maps; here one fixed A is reused. For Jim the term means nothing. | Add the citation (author to verify) or drop "cf. attention rollout" from figure and legend. Keep "repeated attention with residual". |
| R7 | unclear | "chosen without labels" for t = 9 (panel 4 note, legend) | True of the real procedure (held-out laboratory halves, no gold standard). On the toy, t = 9 is not the output of that procedure; the frame and "the depth used" imply the toy selected it, and the toy AP then peaks at t = 9. | "t = 9: depth of the real method (chosen without labels, held-out laboratories)". |
| R8 | cosmetic | strip pill, title, footer | "no learned parameters" appears three times. It is true for the operator (no trained weights); alpha, k_eff, cap are fixed and t is estimated from the data (without labels). "de novo" and "unsupervised" are jargon for Jim; the paper defines de novo as expression-only. | Say it once. "unsupervised | expression data only". If kept: "no trained weights". |
| R9 | unclear | AP inset | "AP" is never defined in the figure; the bars have no axis (0-1 not shown), and the grey random tick (0.23) is a 1 px mark. The orange tick is labelled only in a separate 7 pt caption. | Caption "average precision of ranking same-module pairs"; label the ticks directly ("CLR 0.77", "random 0.23") or add a 0 / 1 axis. |
| R10 | unclear | colour-named annotations | "orange tick", "magenta frame", "S is the input of step 3" (orange 7 pt bold, ~3.3:1 on white) rely on hue; in greyscale (checked) the tick and the frame are indistinguishable from the bars/ink. | Name by form ("thick tick", "framed"), or direct-label. Darken the orange text a notch. |

### Legibility, layout, size

| # | Sev | Location | Finding | Proposed fix |
|---|---|---|---|---|
| R11 | unclear | whole figure, 7.2 x 9.33 in | Every text element is >= 7 pt nominal and legible at 7.2 in; I found no text-on-text collision. Sub/superscripts (ij, ik, eff, +, T) render at about 70% (~5 pt); z+ (the clipping that matters) has a "+" at ~4.5 pt. At 3.5 in nothing is legible (as before; double-column only). 9.33 in is 237 mm tall: within the usual 247 mm full-page ceiling of Nature-family journals (check the target's guide, not verified here) but leaves no room for the ~260-word legend on the same page, so in practice it is a full-page figure with the legend elsewhere. I would cut to about 7.5-8 in. | Cut order (about 1.4 in in total, none loses content that the legend does not carry): (1) footer four lines (~0.55 in), duplicated by the legend (toy parameters, diagonal convention, shading rule); keep one line "toy; real method: k_eff 10, top 50, t = 9" (the badge and the strip already say toy). (2) Panel 4 notes (L470-475, ~0.4 in): residual, value step, magenta frame go to the legend. (3) Panel 1 "swap-in" pill and dead space under the matrix (~0.3 in). (4) The six-line output wrap in the strip. Do not cut the second curve of panel 2 (it carries the symmetry). |
| R12 | cosmetic | panel 2 curves | sigma_i / sigma_j labels sit on the curve line; D_ij label sits on the curve tail; mu label touches the dashed mu line and the bracket tick. Readable, slightly busy. | Move sigma labels above the bar, D_ij above the dot with a leader. |
| R13 | cosmetic | identity coding | Solid = gene i, dashed = gene j works across D, S, A and the bars. One collision: in each curve panel the dashed vertical black line is mu (not gene j), directly beside the "row j" legend line which is also dashed. Also the dashed j-slot frame in the panel-3 bars leaves stray dashes above and below the bar. | Draw mu as a short solid tick or grey dashed; clip the bar-slot frame to the bar area. |
| R14 | cosmetic | colour key | Blue = dependence, orange = CLR/S, violet = softmax/A, magenta = walk: consistent with `figures/README.md` in strip bars, ramps, arrows, tags and the inset (orange tick = S, magenta bars = K_t). Aqua is correctly absent. No conflict found. The grey random tick is neutral, correct. | None. |
| R15 | unclear | panel 3 k cue | k's S_ik bar is ~8% of the tallest bar (0.22 vs 2.8); the "light = positive but cut" bar is almost invisible, yet the whole panel-3/4 story rests on it. The note text carries it, the picture does not. | Mark the k slot with an arrow/label on the S bar chart ("k: 6th, cut") or add a y-break; or normalise bars to the 6th value for display. |
| R16 | unclear | panel 4 reading order | Formula P = 0.5I + 0.5A and K_t sit under the AP inset, away from the pink arrow that applies them; the small pink arrows between t = 1, 3, 9 are minimal. Fine for a figure that is read with a legend, weak for a 20-second read. | Put "P = 0.5I + 0.5A, t repeats" on or above the 3 -> 4 arrow; K_t definition into the legend. |

### Accuracy

| # | Sev | Finding |
|---|---|---|
| A1 | none | Re-ran a scratch copy of `analysis/fig_method.R` (output redirected to the session scratchpad; committed figures untouched). It printed `D_ij = 0.311; z_ij = 2.22; z_ji = 1.77; S_ij = 2.82; A_ij = 0.64; A_ik = 0; K_t[i,k] = 0.000 / 0.025 / 0.078; mu_i 0.0904, sd_i 0.0994, mu_j 0.1038, sd_j 0.1167; S_ik = 0.219 (6th positive); AP CLR 0.770; K_t 0.693 / 0.773 / 0.906; random 0.231; all assertions passed`. The regenerated PNG is byte-identical to the committed one (`cmp`). Every printed number matches: 0.31, 2.2 / 1.8, 2.8, 0.64, 0.000/0.025/0.078, 0.69/0.77/0.91, CLR 0.77, random 0.23. keff 3 / cap 5: asserted for all 14 rows (exp(entropy) = 3 within 1e-4; cap binds for row i; all rows <= 5 candidates). (2.2 + 1.8)/sqrt(2) = 2.83 -> 2.8 consistent with 2.82. "toy: 14 x 60" matches `G = 14; N = 60`. |
| A2 | cosmetic | Assertions are meaningful: S recomputed in R independently of the C++ (L66-72), row sums, keff for all rows, cap binding, A_ik = A_ki = 0, k is the 6th positive, K_1 = 0 exactly, K_1 < K_3 < K_9 at (i, k), K_t equal to the lines of `baselines.R` (evaluated, not copied) and to repeated squaring, AP unit test. Soft spots: the baselines check reuses the same fitted operator `f$operator` (fine, it tests the walk, not A); and L117 `AP_K9 > AP_CLR` encodes a desired conclusion (see R3). Nothing stops a silent mismatch between label text and numbers beyond the sprintf use, but all labels are sprintf-ed from the computed objects (checked by reading). |
| A3 | none | Methods (paper.md, commit 1bed906) versus code. z+ = max(z, 0) and S = (z+_ij + z+_ji)/sqrt(2): `clr_core.cpp` L270-300 (clip, then Stouffer `(a + b) * inv_sqrt2`; Euclidean is the alternative branch). Sample sd (n - 1): L283. Diagonal: "set to zero and kept among the G values that give mu_i and sigma_i ... S zeroed on the diagonal": L266 zeroes the diagonal before the row mean/var loops over all G entries; L292 zeroes S's diagonal. Match. K_t = 1/2 (P^t + P^t T), P = (1 - alpha) I + alpha A: `baselines.R` L185-187. Legend parameters (k_eff 10, top 50, t = 9; toy 3, 5, 1/3/9) match the figure footer. Not re-verified here: the "rows without candidates keep a self-loop" sentence (see `docs/plans/self_loop_trace.md`), irrelevant to the toy (no empty rows). |
| A4 | cosmetic | Legend versus image: legend says "framed t = 9", "inset gives toy co-membership AP", "normal fit; mu_i, sigma_i marked", k ranks sixth, cap 5: all present in the image. Legend omits the MI < 0 shading, "white = clipped to 0", and the grey slot meaning (all in the image). The legend is about 260 words, long for a figure legend; cuts in R11 shorten it. It still states "(cf. attention rollout)" (R6) and "first two published CLR" (R4). |

## 3. Top 5 findings

1. R1: the published-versus-added bracket is not in the strip (only 7 pt italic corner tags); Jim cannot see what is new in 20 s.
2. R2: output outlines are the planted truth, and "regulons, complexes, operons (generic types)" contradicts the paper's own results on complexes/operons; the map looks like module calls.
3. R3: the AP inset is the only "potential" cue, and it is favourable by construction; its gain comes from the tight module (not the chain/k story), where real data show losses. Caption and pointer needed; plus its low visual weight (R9).
4. R4/R5/R6: step-2 combination and MI estimator are not the 2007 ones; "any dependence measure" untested; "cf. attention rollout" has no reference.
5. R11: 9.33 in is too tall for a figure plus a 260-word legend; cut the footer, the panel-4 notes and the panel-1 pill (about 1.4 in).

## 4. Recommendation

Minor fixes, no restructure. The content is accurate, the toy is honestly badged, and the colour and identity coding are consistent. Do, in order:
1. Add the published/added bracket to the strip (R1); relabel the outline caption and drop "complexes, operons" (R2).
2. Re-caption the AP inset: define AP, direct-label ticks, "toy, favourable by construction; real data: Figs. 2-6" (R3, R9).
3. Reword tags: "after CLR (2007), here with B-spline MI and Stouffer"; "e.g. MI, correlation (tested)"; add or remove the rollout reference (R4-R6); t = 9 as "depth of the real method" (R7).
4. Trim to about 7.9 in (R11) and name annotations by form, not only hue (R10).
5. Optional, for the big-picture reader: in the strip, show plain-CLR S next to K_9 (both already computed) under "published CLR" vs "this work", so the potential is visible in 5 s without a number; keep the toy badge.

## 5. Unverified

- The 2007 MI estimator (kernel density) is from memory; the Euclidean-versus-Stouffer difference is verified from `R/calibrate.R` and Methods.
- Journal size limits (247 mm) are from memory.
- Per-module AP and the 20-seed check are scratch computations (same toy generator, seeds 1-20, 100-400 tie-break repetitions); not in the repo and not paper claims.
- No human or print-size review; the 7.2 in view is the 300 dpi PNG at native size, greyscale and crops were rendered with PIL.
