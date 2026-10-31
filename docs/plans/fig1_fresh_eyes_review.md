# Fig. 1 (method figure): fresh-eyes review

Audience: Jim, a senior collaborator who is (a) a big-picture person and (b) a high-impact person; the figure must show the potential of the idea, not only a worked toy.

Reviewed: `docs/whitepapers/figures/fig1_method.png` (2160 x 2355 px, 7.2 x 7.85 in at 300 dpi), `analysis/fig_method.R`, brief `figures/src/fig1_method.txt`, legend and Methods in `docs/whitepapers/paper.md`, `analysis/baselines.R` and `src/core/clr_core.cpp` (ground truth for the walk and the CLR z-scores), graphical abstract (for comparison), PaperBanana reference `pb_runs/run_20261030_160050_6caa1b/final_output.png`.
Reviewer: a fresh agent context, but the same model family as the figure's builder, so this is an independent re-derivation, not a human fresh-eyes review. Review only: no figure code, paper.md, STATUS.md or tables were touched. No PaperBanana run was made (not needed: the layout reference already exists, and the question is content, not layout). Scratch scripts are in the session scratchpad, not in the repo.

## 1. Verdict in one paragraph

Accuracy is good: every printed number reproduces from the package (a copy of the script run to a scratch output regenerates the committed PNG pixel-identically) and I found nothing numerically wrong. Clarity for a big-picture reader is the weak side. In 20 seconds at 7.2 in the figure answers (ii) partly (the four tags give the transformer correspondence; nothing shows what softmax plus walk buy over plain CLR) and does not answer (i) or (iii) at all: no input is drawn (the figure starts at D), no output is drawn (it ends at K_t bars of one gene), no biological object (module, regulon, complex, gene name) appears, and "toy" status is only in 7 pt footer text. The area is allocated inversely to novelty: panel 1 (plain MI) is the largest element, panels 3-4 (the contribution) are the smallest, and nothing separates published CLR (panels 1-2) from what is added here (panels 3-4), which the graphical abstract does with brackets.

## 2. Legibility and layout (rendered at 7.2 in and shrunk to 3.5 in)

7.2 in (the real double-column target): every text element is readable. Fonts are 7 pt (tags, axis "MI", footer lines, "dark: top 5 kept", "sums to 1", box/triangle note) to 9.5 pt (panel titles); math is 7.5-8 pt, but sub/superscripts (ij, +, eff, t) render at about 70%, i.e. roughly 5 pt, which is below the 7 pt rule for the symbols that carry the content (z⁺_ij, D_ij, K_ik). 7 pt text is at the floor, with nothing to spare.

3.5 in (single column): not legible, by design. All text drops to 3.4-4.6 pt, the footer and every 7 pt label are unreadable, the small panel-3/4 matrices (cells 0.07-0.08 in, about 2.5 pt at 3.5 in) become texture. The figure should be declared double-column only; do not plan a single-column version of this layout.

Greyscale (checked): luminance ramps hold up, the blocks and the traced cell stay readable, but the blue/rose identity of gene i versus gene j is lost (both grey), leaving only row-versus-column position to tell them apart. Acceptable.

Dense areas: panel 2 (two curves, two z formulas, two converging lines, S formula, "z⁺" note, S matrix and its label within 3.5 x 3.7 in) and the bottom of panel 4 (K_ik values, K_t formula, box/triangle note stacked 0.16 in apart above the footer rule). No collisions between text elements found. Panel-3/4 matrices at 1.0-1.1 in are small and the three K_t matrices look alike; the row-i bars below them carry the message, which is right.

Reading path: chips 1 -> 2 (big arrow) work. 2 -> 3 is a 5 in elbow connector (see F11), 3 -> 4 a big pink arrow, then 4 is left without an end. The arrows are not labelled except "softmax" (the pink 3 -> 4 arrow and 1 -> 2 arrow carry no label; "P = 0.5I + 0.5A" is at the far top right of panel 4, away from the arrow that applies it).

## 3. Accuracy: numbers re-derived

Ran the script logic (copy with output redirected). Toy: seed 7, G = 14, N = 60, i = 4, j = 5, k = 9.

| Printed | Re-derived | Status |
|---|---|---|
| D_ij = 0.31 | 0.3110 | ok |
| z_ij = 2.2 | 2.22 (mu_i 0.0904, sd_i 0.0994) | ok |
| z_ji = 1.8 | 1.77 (mu_j 0.1038, sd_j 0.1167) | ok |
| S_ij = 2.8 | 2.82; (2.22 + 1.77)/sqrt(2) = 2.82; printed rounded 2.2 + 1.8 gives 2.83, consistent | ok |
| A_ij = 0.64 | 0.636 | ok |
| K_ik = 0.000 / 0.025 / 0.078 at t = 1, 3, 9 | 0.000 (exact zero) / 0.025 / 0.078 | ok |
| keff = 3 | all 14 rows have exp(entropy) = 3.00 (every row has 4-6 positive candidates, so none is uniform); row i has 5 candidates = cap | ok |
| tau_i (not printed) | 0.77 for row i | n/a |
| matrices | all drawn from the same objects D, S, A, Kt that produced the numbers | ok |

Method consistency with Methods: S = (z⁺_ij + z⁺_ji)/sqrt(2) with each one-sided z clipped before combining (clr_core.cpp lines 285-302: clip, then sum times 1/sqrt2): matches. Sample sd (n-1): matches. Softmax over the top-`cap` positive scores, temperature by bisection so exp(entropy) = keff: matches `softmax_operator`. P = (1-alpha) I + alpha A, alpha = 0.5, K_t = (P^t + P^tT)/2: the script's three lines are copied verbatim from `analysis/baselines.R` lines 185-187 (checked by eye; not tested against it, see F7). The toy differs from the real method in keff 3 (10), cap 5 (50), depths 1/3/9 (9); the footer says so, in 7 pt at the bottom ("real method: 10, 50, 9", the three numbers unnamed).

What is not a number error but should be fixed in the picture or text is listed below.

## 4. Findings

Severity key: wrong = factually incorrect; misleading = a careful reader draws a false conclusion; unclear = a reader cannot tell; cosmetic.

No finding rated "wrong".

### Accuracy and notation

| # | Sev | Location | Finding | Proposed fix |
|---|---|---|---|---|
| F1 | misleading | footer line 3, `fig_method.R` L317 | "grey diagonal = self, excluded" and "normal fits of the row's 14 values (ticks)". In the code (and in `clr_core.cpp`, following the original clr.m) the diagonal is zeroed and then **included** in mu_i and sigma_i, so the self entry counts as a zero among the 14 values (a grey tick at 0 is drawn in each curve panel). "Excluded" is true for the softmax and the matrix colouring, false for the z-score. Paper Methods does not state the convention either. | Footer: "diagonal (self) set to 0 and kept in mu, sigma, as in CLR; drawn grey". Add one sentence to Methods. |
| F2 | misleading | panel 2, `dist_panel` | Each curve panel has its own x scale (0 to 1.08 x that row's max): D_ij = 0.311 sits at 93% of the axis in row i and 84% in row j (row j's max is 0.343). The same value appears at two different positions with no axis numbers, which reads as two different D_ij. | Use one shared MI axis for both curves (0 to 0.37) with 2-3 tick values. |
| F3 | misleading | panel 2 curves | The normal fit is cropped at MI = 0, where its density is still 66-68% of the peak, so the curve starts as a cliff mid-height. The fit puts 18-19% of its mass below 0 (impossible for MI); that is what CLR assumes, but the picture hides it. mu_i and sigma_i, which the formula uses, are not marked anywhere, so the z formula is not "visible" in the drawing. | Mark mu_i as a tick and draw the bracket from mu_i to D_ij labelled "z_ij sigma_i"; extend the axis a little below 0 with a lightly hatched "MI < 0" zone (or label "normal approximation") so the truncation is honest. Keep it small. |
| F4 | unclear | panels 3 and 4, gene k | k = gene 9 is a module mate of i whose CLR score is positive (S_ik = 0.22, the 6th positive entry in row i) but is cut by cap 5, hence A_ik = 0. This is the reason for the whole panel-4 story, yet it is only implied: k is the lightest bar at the right of the S bars (not marked) and is introduced for the first time in panel 4 in a 7 pt note. In the real method (cap 50, keff 10) indirect reach is not the same mechanism; the toy should not suggest it is. | Mark k on the S and A bars in panel 3 (small triangle plus "k: cut by the cap"); introduce k already in panel 1 (a third small tag on the matrix). State in the note: "toy: cap 5 cuts a weak positive edge". |
| F5 | misleading | tags | See table in section 5. Main points: "query-key score" is acceptable only if it says the genes are both query and key (self-attention over genes) and MI is parameter-free and symmetric; "stacked attention layers" overclaims (one fixed A reused t times, no value projection, no learned weights, no nonlinearity); the output is the attention map itself (no A·V value step), and "value" appears nowhere. | See section 5. |
| F6 | unclear | notation across figure/paper/brief | Figure uses K_t for the walk readout; paper.md never names it (it writes "1/2 (P^t + P^tT)" in words), the brief calls it S_t, and S in the paper is the CLR matrix. Figure uses z⁺; paper formula is S_ij = (z_ij + z_ji)/sqrt(2) with the clipping only in words. Paper Results says the z-scoring "plays the role that softmax plays in attention"; the figure labels it "scaling the logits" and gives the softmax its own panel. | Define K_t in paper Results ("readout K_t = 1/2 (P^t + P^tT)") and write z⁺ in the paper formula; pick one analogy wording for CLR (the figure's "per-query normalization of the logits" is the more defensible one) and align the Results sentence. Brief is stale (README already says superseded); stop citing it in the legend. |
| F7 | unclear | assertions, L61-L66 | Meaningful: S_ij and the whole S recomputed independently from D (L61, L63); row sums of A = 1; keff = 3 for row i; cap binds for row i (5 candidates); A_ik = 0; K_1 < K_3 and K_3 > 0 at (i, k). Not covered: keff for the other 13 rows (I checked: all 3.00); K_t is **copied** from `baselines.R`, not computed by calling it, so drift is possible; A_ki = 0 (needed for K_1,ik = 0 exactly, true here, unasserted); K_9 > K_3 at (i, k) (the printed 0.078 versus 0.025 trend); the shown A_ij = 0.64 is the max of the row (the A bar is normalised to it). The K_1 < K_3 assert is nearly vacuous given K_1 = 0. | Add: all-rows keff, A_ki = 0, K_3 < K_9, and one test that Kt equals the package/baselines walk on this toy (source the function if `baselines.R` exposes one, otherwise a `diffuse()` equivalence on an identity "expression"). |
| F8 | misleading | footer, panel 3 note, panel 4 | Toy status is not stated where numbers appear. A reader can mistake D_ij = 0.31, S_ij = 2.8, K_ik = 0.078 for results. "real: 50", "real: 10", "real method: 10, 50, 9" are 7 pt and in three places, the last one with unnamed numbers. | A badge "Worked example (toy: 14 genes, 60 experiments, planted modules)" in the figure's header; footer: "real method: k_eff = 10, top 50, t = 9". |
| F9 | unclear | panel 4, depth | The label-free choice of t (held-out laboratories; median 9) is a headline property and is absent; the three depths 1/3/9 are shown with no statement that 9 is the one used. | Outline the t = 9 matrix and add "t* = 9, chosen on held-out labs, no labels" (7 pt). |
| F10 | cosmetic | various | "MI(i, j)" in panel 1 renders with oversize parentheses (plotmath); "i"/"j" in matrix margins and "gene i" tags are repeated in 4 panels (OK). Bars are normalised to the row maximum, so the bar heights are relative and the A_ij = 0.64 label sits at full height; "sums to 1" is correct but the reader cannot see it. K_t row-i bars share one scale across t (good); matrix colours use per-matrix maxima (0.33, 0.31, 0.29, close enough). | Optional. |

### Reading path and visual design

| # | Sev | Location | Finding | Proposed fix |
|---|---|---|---|---|
| F11 | unclear | the S connector, L234-L242 | The orange elbow leaves the bottom-right of S in panel 2, runs down, then **along the horizontal divider** (orange on grey, so it looks like the divider), crosses the vertical divider, drops, jogs left and lands on the top-right corner of the S matrix in panel 3 (x = 1.62 in on a matrix spanning 0.56-1.68 in), not at its centre, so the arrow tip appears to point at a specific column (gene 13). "same S matrix" is placed in panel 2 at the far end, 4-5 in from the landing. Duplicating S is itself a signal that the layout forces a long return. | Either (a) rearrange so panel 3 sits directly below panel 2's S (snake order 1, 2 / 4, 3 with chips making the order explicit; 3 -> 4 then runs leftward), or (b) keep the grid, drop the duplicate S and let panel 3 show only "row i of S -> row i of A" with a short pointer, or (c) minimum change: land the arrow at the matrix centre-top, route it off the divider line, and move "same S matrix" next to the landing. |
| F12 | misleading | colour use | Identity colours collide with functional colours: gene i = blue (also the dependence/MI ramp; the blue i band sits on the blue D matrix), gene j = rose/pink (also the depth/walk magenta; the pink j band and pink row-j curve look like the walk panel's colour, and the pink 3 -> 4 arrow). README key: blue = dependence, magenta = depth. The code comment says identity colours are "never matrix fills" but the colour confusion is in bands, curves and text, which is where a quick reader looks. | Give gene i / gene j identity colours outside the README key (the key reserves blue, orange, violet, magenta, aqua, grey, yellow, and platforms use red/green/brown/dark grey, so a restyle must pick two dark neutrals or encode row versus column by outline style), or keep the 2007 blue/rose pair only in panel 2 and use black outlines elsewhere. Needs a decision; flag for the restyle. |
| F13 | unclear | area/hierarchy | Panel 1 matrix is 2.3 in (largest element, about 20% of the page) for the least informative step; panels 3-4 (the contribution) have 1.0-1.1 in matrices. | Shrink panel 1 to about 1.6 in, enlarge 3 and 4, and use the freed height for a strip (section 6). |
| F14 | unclear | existing versus new | Nothing marks panels 1-2 as published CLR (2007) and 3-4 as what this work adds, although that is the paper's claim and the graphical abstract does it. | A thin bracket/tint under the chips: "CLR (Faith 2007)" over 1-2, "added here" over 3-4. |
| F15 | unclear | panel 2 linkage | The z values (2.2, 1.8) are not visibly carried into S_ij = 2.8: the converging lines start at the curve axes, the S formula is at the top right, the numbers are bottom left. | Show "(2.2 + 1.8)/sqrt(2) = 2.8" under the formula (it is a computed number; can be asserted). White cells in S are clipped zeros (9 of 14 entries of each row are z < 0): say "white = clipped to 0". |

### Content a big-picture reader needs

| # | Sev | Finding | Proposed fix |
|---|---|---|---|
| F16 | unclear | **No input and no output.** The brief asked for both ("Input: expression matrix... Output: S_t used to recover regulatory modules"); the figure shows neither. The 60 experiments appear only in a 7 pt footer. | Strip (section 6). |
| F17 | unclear | **No biological object.** Planted modules (genes 4-9 and 10-13) are visible as dark blocks in D but unlabelled; no gene names, regulon, complex, operon, TF. i and j are adjacent genes of the chain. | Label the blocks "module (e.g. one regulon / one complex)" at the matrix margin; in the output of the strip show the recovered blocks. |
| F18 | unclear | **What softmax and walk buy over plain CLR is not shown.** The toy supports a small, honest demonstration, computed by me in scratch (not in the repo): row i of plain CLR ranks noise gene g3 (S = 0.86) above the true module mates g8 (0.51) and g9 = k (0.22); at t = 9 the order is reversed (K_9: g8 0.100, g9 0.078, g3 0.044). Co-membership AP in the toy (21 co-members of 91 pairs, base 0.23): CLR S 0.77; K_1 0.70; K_3 0.77; K_9 0.91. One seed, planted structure, 14 genes: illustrative only, must be labelled so and computed inside `fig_method.R` if used, never quoted as evidence. Note the "share of row-i mass on module mates" is not monotone (S 0.89, A 0.95, K_1/3/9 0.92/0.90/0.87), so do not use that statistic. | Add a one-line inset in panel 4, computed in the script: either "noise gene g3 outranks module mates g8 and g9 under CLR, not at t = 9", or "toy co-membership AP: CLR 0.77, t = 9: 0.91 (illustrative)". The rank of k alone is too weak (CLR 6th, t = 9 5th). |
| F19 | unclear | **Potential.** Nothing says what the idea could become: swappable dependence measure (the paper already uses MI and rank correlation, Discussion), per-query width and depth as tunable dials (resolution vs reach, Fig. 6), a template for transferring attention machinery to any gene-gene compatibility, and training-free versus learned attention maps (the comparison is [PENDING], so no outperformance claim can be drawn). | Strip and a plug-in slot, section 6. |

## 5. Transformer correspondence: what is defensible

| Figure label | Verdict | Wording that does not overclaim |
|---|---|---|
| "query-key score" (MI) | acceptable with a qualifier. In a transformer, query and key are learned projections and the score is asymmetric; here each gene is both query and key (self-attention over genes), the score is parameter-free and symmetric. | "pairwise compatibility (stands in for the learned query.key score)" |
| "scaling the logits" (CLR) | defensible: per-row standardization of the logits before the softmax, plus a two-sided combination and clipping that have no transformer counterpart. | "per-query normalization of the logits (z-score vs the gene's own background)". Mention the two-sided Stouffer combination, which makes S symmetric; the asymmetry of A comes from the per-row temperature and top-k. |
| "softmax over the keys" | correct. The per-row temperature set by k_eff (versus a fixed 1/sqrt(d)) is a difference worth one phrase. | keep; add "temperature per query" |
| "stacked attention layers" | overclaims. It is one fixed operator applied t times, no value projection, no MLP, no learned weights; the 0.5 I term is a residual connection. The accurate analogue is attention rollout (Abnar and Zuidema 2020, residual 0.5 I + 0.5 A multiplied across layers; citation to be verified, it is not in the reference list). | "repeated attention with a residual (rollout-like)"; label "0.5 I" as the residual term |
| "value" | absent. The readout is the attention map K_t itself, not A V. | One phrase in the legend: "the readout is the attention map, not a weighted sum of value vectors" |
| "No learned parameters" | correct for the operator. Depth t is the one chosen quantity and is chosen without labels. | keep |

## 6. Big-picture question: what is missing, and where it belongs

Does Jim get (i), (ii), (iii) in 20 seconds? (i) no (no input, no output); (ii) partly (tags yes; benefit over plain CLR no); (iii) no. The figure reads as a competent worked example of a method, not as a statement of why the method matters.

Recommended: a **short strip in Fig. 1**, and **one line of potential in the graphical abstract**, not the reverse.

Why in Fig. 1: input and output are properties of the method figure; a method figure without them cannot be read standalone, and the graphical abstract already has the input (compendium panel) and the results (bars) but cannot carry the 4-step mechanism. Why not "potential" only in the abstract: Jim will look at Fig. 1 first; the strip is where the toy becomes a general operator.

Strip, about 0.9-1.0 in tall across the top of Fig. 1 (recover the height by shrinking panel 1's matrix from 2.3 to about 1.6 in; total height stays about 7.9-8.3 in):

1. Left: a small genes x experiments heat strip, caption "expression compendium: thousands of genes x hundreds to thousands of experiments, many laboratories". Computed from the toy X (so it is real).
2. Middle: the four chips as a compact flow "dependence -> per-gene calibration (CLR, 2007) -> softmax -> multi-step walk (added here)", with a bracket marking published versus added, and the pill "no learned parameters" with an "unsupervised, de novo" word (as in the paper intro).
3. Right: the output, the gene x gene association map K_9 of the toy with its two blocks outlined, captioned "modules: regulons, complexes, operons" (the Methods say regulons and protein complexes are the benchmark; operons are in the brief).
4. A single plug-in note on panel 1's chip: "any dependence measure: MI, rank correlation" (supported by Discussion; do not list learned embeddings as a result).
5. Optional dial: small "t: fine -> coarse" arrow over panel 4 ("fine structure to larger modules"), but the brief's "small t keeps operons, large t regulons" must be checked against Fig. 6c before it is printed; K562 CORUM shows attention losing to CLR on tight complexes, so the dial must not promise complexes.

Potential statement (one line, footer of Fig. 1 or the abstract): "Training-free counterpart of learned attention: the same operations (per-query normalization, softmax, depth) applied to a measured, interpretable gene-gene score." Allowed claims: no training, interpretable at every step, applicable de novo. Not allowed until the scFM comparison is scored [PENDING]: any statement that it matches or beats learned attention. Also keep the honest scope: gains are large on loose E. coli regulons (+16% pooled AUPR over CLR), none on single cell, losses on tight complexes. The strip should show the mechanism and generality, not imply uniform wins.

## 7. Legend versus image

Current legend: "Dependence estimation, per-query context normalization (CLR), softmax with a per-query temperature, and iterated attention, alongside the corresponding transformer components. Drawing brief: figures/src/fig1_method.txt."

Gaps: does not say it is a toy worked example with one traced pair (i, j) and a third gene k; does not define the colours (blue/rose bands = gene i row / gene j column; heat ramps per step); does not give the real parameters; "iterated attention" omits the lazy residual and the symmetrised readout; cites an internal brief that the README calls superseded; names transformer components but not the caveat that this is an analogy (no values, no learned weights). Proposed (about 90 words): "Fig. 1 | CLR as attention, traced through a toy example. A 14-gene toy with planted modules, computed with the package, follows one gene pair (i, j) and a third gene k. (1) Dependence D (mutual information; stands in for the query-key score). (2) CLR: each row of D is z-scored against the gene's own background; the two one-sided scores, clipped at 0, are combined into S (per-query normalization of the logits). (3) Softmax with a per-gene temperature over the top-scoring partners gives the attention matrix A. (4) The lazy walk P = 0.5 I + 0.5 A (residual plus attention) is iterated t times and read out as K_t = (P^t + P^tT)/2; t = 9 is chosen without labels. Toy: k_eff = 3, top 5, t = 1, 3, 9; real method k_eff = 10, top 50, t = 9. No learned parameters; the readout is the attention map itself (no value vectors)."

## 8. PaperBanana reference (inspiration only)

Layout contributions worth keeping from the AI draft: the tag chips beside the panel titles (already adopted), the S -> A matrix pair with bars underneath (adopted), K_t formula placed with the t labels. Content that is invented and must not be reused: a D matrix of four equal-sized blocks on a regular grid with the traced cell outside any block; an S matrix that is near-uniform orange; A as isolated dots; K_t matrices as uniform purple with identical diagonals; bar charts of row i with arbitrary shapes (its row-i bars for K_t happen to carry our 0.000/0.025/0.078 only because those were in the prompt); a "dashed gene k" box on an empty bar chart. The current code figure is already the faithful version; none of the draft should be shipped. For the strip, a PaperBanana run is not needed; the elements are all computable from the toy.

## 9. Recommended next change set (ordered, each small)

1. Add the input -> output strip with the published/added bracket (F14, F16, F17, F19) and the "worked example (toy)" badge (F8); shrink panel 1 (F13).
2. Fix the accuracy-visible items in panel 2: shared MI axis, mu_i/sigma_i marked, honest cut at 0, "(2.2 + 1.8)/sqrt(2) = 2.8" (F2, F3, F15), and the footer wording for the diagonal (F1).
3. Mark k in panel 3 and say why A_ik = 0; outline t = 9 with "chosen without labels"; add a computed "what the walk buys" number or rank swap labelled illustrative (F4, F9, F18).
4. Reword tags ("repeated attention with residual", "pairwise compatibility", legend caveat on values) (F5).
5. Resolve the S connector (F11) and the identity-versus-key colour collision (F12), which need a decision on layout and palette.
6. Extend assertions: all-row keff, A_ki = 0, K_3 < K_9, walk equals the baselines/package walk (F7).
7. Text: define K_t and z⁺ in paper.md, align the "role that softmax plays" sentence, replace the legend (F6, section 7); add the diagonal convention to Methods (F1).
8. Add one line of potential, within the allowed claims, to the graphical abstract or the Fig. 1 footer.

## 10. Unverified or incidental

- `R/clr_attention.R` L348/L568 and Methods both say rows without candidates get a self-loop, but `softmax_operator` (R/packed.R L86-96) documents "rows without candidates stay empty" for the packed path; I did not test which is true for the packed path. The toy has no empty rows, so Fig. 1 is unaffected.
- The depth dial (small t fine structure, large t coarse modules) is from the brief and the abstract's resolution-versus-reach trade; Fig. 6c supports it for regulons versus operons only if it shows what the legend says; I did not check Fig. 6.
- The attention-rollout citation is from my memory, not verified.
- The toy numbers in F18 come from a scratch script run on the toy; they are not in the repo and are not part of any paper claim.
- No human or independent-instance reviewer looked at the figure at print size; the 7.2 in and 3.5 in views are renders of the PNG, not a print proof.
