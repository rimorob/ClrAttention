# Abstract points: completed CLR attention (Paper 1)

*Outline for the PI's rewrite, 2026-10-21. Every result comes from `docs/paper_plan.md`, the record of what survived independent statistical review, with its status in the plan's terms: **worked**, **mixed**, **did not work**; "QC" means checked by the independent reviewer. Where a number comes from somewhere other than the plan (`paper.md`, `sota_table.md`, `results/macro_aupr_v2/macro_aupr.log`), the source is named. The drafts built from these points are in `abstract_draft.md`.*

## (a) The problem and the gap

- **De novo inference is still the hard case.** Recovering regulatory programs from expression data alone, without motifs, binding data or curated priors, remains difficult. On the standard single-cell benchmark most published methods perform close to random (BEELINE; Pratapa et al. 2020).
- **CLR has lasted, within its class.** Context likelihood of relatedness (Faith, Hayete et al. 2007, *PLoS Biol.* 5:e8) is still a standard for de novo, unsupervised inference. It was an off-the-shelf method in DREAM5 and a component of several of its ensemble predictors (Marbach et al. 2012), and it is the strongest of the methods we tested on BEELINE (see c6). Keep this robustness claim inside that class. Prior-informed methods answer a different question. So do regression methods such as GENIE3, which fits each gene on a supplied list of candidate regulators and scored higher than CLR on directed TF→target edges in DREAM5.
- **CLR already computed attention scores.** Attention normalizes each query's compatibilities within its own context. That move is older than neural attention:
  - Nadaraya–Watson kernel regression (Nadaraya 1964; Watson 1964);
  - fast-weight programmers (Schmidhuber 1991), which are formally equivalent to linear attention (Schlag, Irie & Schmidhuber 2021);
  - attention by name (Bahdanau et al. 2014) and the transformer (Vaswani et al. 2017) came later.

  CLR z-scores each gene's dependence on every other gene against that gene's own background. That is per-query normalization. The claim is "the first attention-like operator widely adopted in computational biology" (plan §10, worked), not "the first attention mechanism".
- **The gap.** CLR stopped at the scores. It ranked them as edges and never used a temperature-controlled softmax, aggregation or depth. Over the same two decades, attention research in language modelling turned those components into a mature toolkit. No one has brought that toolkit back to CLR and measured what each piece adds, and what it costs.
- **Why now (context, not our result).** In independent evaluations, learned attention from single-cell foundation models has so far shown little network signal beyond simple baselines. Sources: Kendiukhov 2026; Ahlmann-Eltze et al. 2025; Boiarsky et al. 2024; Kedzierska et al. 2025. Verify these citations first, as `docs/plans/scfm_comparator.md` notes. A label-free, estimated attention operator with a known null is a natural reference point. Kendiukhov 2026 tested perturbation-target prediction, not network recovery; do not conflate the two.

## (b) The idea

- **Complete CLR as attention**, using the components that language modelling made standard, and measure each one:
  - a per-gene softmax temperature over the top 50 CLR scores, with an effective neighbourhood exp(entropy) = 10;
  - iterated attention through the lazy operator P = ½(I + A), read out as ½(Pᵗ + Pᵗᵀ);
  - depth t chosen without labels on held-out laboratories. The criterion is how well one half's attention mass predicts the other half's top-1% CLR pairs.
- **Propagation is the baseline, not the novelty.** Iterated attention is a random walk, related to diffusion kernels and network propagation (Kondor & Lafferty 2002; Zhou & Schölkopf 2004; Cowen et al. 2017). Every comparison includes propagation on uncalibrated scores.
- **Nothing is learned.** There is no corpus, no embeddings and no labels. The CLR scores carry a permutation null.
- **Mechanism (reframed after review, plan §5).** The walk trades pair-level resolution for neighbourhood aggregation. It should therefore help large, loose regulons and hurt tight co-expression cliques (complexes) and direct edges, and that is what we find (c3, c5).

## (c) Supported results, with numbers and review status

**c1. Core gain (M3D *E. coli*).** 4,297 genes; 907 arrays from 466 experiments; RegulonDB 14.5 regulon co-membership; operon pairs excluded; SC = strong/confirmed evidence.
- **Pooled AUPR.** Attention minus CLR is +0.022 [+0.019, +0.026] in a nested bootstrap that re-selects depth in each of 1,000 cluster draws. Attention is ahead in 1,000 of 1,000 draws. Status: **worked, QC**.
  - Relative gain is about 16% (0.022 / 0.136). `paper.md` says 15%, from one split at t* = 8 (0.157 against 0.136). `sota_table.md` gives 1.68 against 1.44 times the base rate at t = 9. Pick one basis and use it everywhere.
  - Over 50 lab splits: +0.022 [+0.017, +0.027]. Status: **worked, QC**.
  - Replicate-aware bootstrap on 907 chips: +0.020 [+0.018, +0.023]. Status: **worked**.
- **Coherence (per-regulon AUROC), SC.** +0.026 [+0.008, +0.045] in the nested bootstrap. Status: **worked, QC**. At all evidence levels it is +0.011 [−0.004, +0.025], not significant. Status: **mixed, QC**.
- **Label-free depth is stable.** Median t* = 9 over 50 splits, and the criterion is within 1% of its maximum for t = 7–11. The softmax candidate cap does not drive the result. Status: **worked, QC**.
- **Not a disguised correction for global programs.** The gain survives removal of 0–20 principal components: M3D +0.017 to +0.022, EcoMAC +0.020 to +0.029. The top components are physiology: growth rate (ribosomal proteins, AUROC 0.98), motility (σ28, FlhDC) and σS. Status: **worked**.

**c2. What the gain needs.**
- **Context normalization is essential.** Softmax propagation on uncalibrated |r| or mutual information recovers about a third of the AUPR gain and none of the coherence gain. Status: **worked, QC**.
- **Shapley decomposition over six components.** AUPR: depth +0.018, CLR normalization +0.013. Coherence: CLR +0.025. Removing the softmax costs the most coherence (−0.043). Status: **worked**.
- **Hard selection does not work as an operator.** The permutation-FDR graph leaves 62–77% of genes isolated. FDR plus the top-10 neighbours gains +0.013 AUPR and no coherence. Status: **did not work**; this supports "the softmax matters".

**c3. Mechanism: resolution against reach (QC).**
- **Helps large, loose regulons.** The edge sits in regulons of ≥ 100 genes (wins 60–81% against CLR). Small regulons are a coin flip. The pooled gain comes from ranking across regulons, not within them. Status: **mixed, QC**.
- **Hurts tight co-expression cliques.**
  - 86% of the K562 CORUM loss comes from 10 ribosome and pre-rRNA complexes (60S AP 0.19 against 0.66 for CLR). Status: **QC**.
  - Knockdown-response modules lose too, and the loss grows with module size (c5).

**c4. Competitors and hub control (*E. coli*).**
- **Raw pooled AUPR: attention is highest.** +0.010 to +0.017 over GENIE3 and +0.008 to +0.028 over CLR; jackknife intervals exclude 0. Status: **mixed, QC**.
  - The gain is concentrated in global regulons: CRP, σS, Fur, IHF and FNR carry 62% of positive pair incidences.
  - The same comparison as fold over the base rate (M3D SC, `sota_table.md`): attention 1.68, GENIE3 1.51, PIDC 1.45, CLR 1.44, |r| 1.23, GRNBoost2 1.18, hub score alone 1.08.
- **Under hub control (within-hub-stratum percentile), attention still leads on SC truths.** Status: **mixed, QC**.
  - AUPR lead on M3D: over CLR 0.015, GENIE3 0.009, PIDC 0.015.
  - AUPR lead on EcoMAC: over CLR 0.020, GENIE3 0.009.
  - At all evidence levels it leads only CLR, and only on EcoMAC (0.009).
  - Hub control removes one third (SC) to three quarters (all evidence) of the raw lead over CLR on M3D.
  - An earlier QC check against CLR on M3D SC gave +0.016 [0.012, 0.020], ahead in all 10 hub deciles in every jackknife replicate. Status: **worked, QC**.
- **Per regulon: parity, not superiority.** The comparison uses mean and median AP differences, Holm-corrected, and requires both the regulon and the gene uncertainty to exclude 0. Status: **mixed, QC**.
  - No difference from CLR, GENIE3 or PIDC on *E. coli* (Holm p ≥ 0.21).
  - Attention beats GRNBoost2 (4 of 4 cells) and |r| (3 of 4) on the mean.
  - Share of regulons improving over CLR under hub control: EcoMAC 66% (p = 0.002); M3D 54% (p = 0.11, not significant). Status: **QC**.
- **Direct edges: attention is not the better GRN inferrer.** Status: **mixed, QC**. The claim is module recovery, not GRN inference.
  - It is last of the top three on direct TF→target pairs (TF-node AUPR ratio on M3D: CLR 6.75, GENIE3 6.54, attention 5.17).
  - GENIE3 wins median per-regulon coherence and per-TF AUROC (EcoMAC 0.690 against 0.596).

**c5. Human (K562 essential-gene Perturb-seq pseudobulk compendium; Replogle et al. 2022).**
- **DoRothEA A–B.** The lead over CLR survives hub control: +0.006 [0.002, 0.010] AUPR. Status: **mixed, QC**. Before hub control the leads were +0.008 (DoRothEA) and +0.009 (CollecTRI). Status: **worked**, but now superseded.
- **CollecTRI.** No signal under hub control; the hub score alone is at least as good as attention. Status: **QC**.
- **TF truths per regulon.** The differences are tiny (|ΔAP| ≤ 0.0025). Status: **QC**.
- **CORUM complexes.** −0.041 AUPR against CLR. Status: **did not work**. Per complex, CLR and GENIE3 beat attention; this is fragile (the delete-50% jackknife interval includes 0).
- **Knockdown-response modules.** Mean per-module AP is higher than attention's for CLR (+0.097), GENIE3 (+0.084) and |r| (+0.091); attention wins only 2% of modules. Status: **mixed, QC**. Two caveats: 62% of the modules are near-duplicates, and they are built from a held-out replicate of the same signal, which favours direct similarity.

**c6. Single cell (BEELINE: 7 datasets × 2 gene sets, 72 dataset–truth combinations).**
- **CLR is the strongest method tested.** Median AUPRC ratio 1.58, against 1.32 for correlation, 1.20 for mutual information and 1.08 for partial correlation. Status: **worked**. GENIE3, GRNBoost2 and PIDC have not yet been run on BEELINE (job 0890, queued).
- **Normal-score CLR beats MI-based CLR in every gene context.** CLR-NS is CLR on |r| of Blom scores. Gain +0.07 to +0.12 in median AUPRC ratio (53–57 of 72). Status: **worked, QC**. QC caveat: the 72 combinations come from 7 datasets, so report per-dataset blocked tests, not a paired p over 72.
- **Attention adds nothing over CLR-NS.** Status: **did not work, QC**. Wider gene contexts give no gain either. Status: **did not work, QC**.
- **Held-out depth transfers.** t* = 5–11 in all 14 cases, bracketing the *E. coli* depth (from `paper.md`; not a separate item in the plan).

**c7. Supplementary-strength results (sound, but not for the abstract).**
- On EcoMAC, platform mixtures leave 8.4% of pairs with far more mutual information than their correlation implies.
- One clean platform (735 arrays) beats all 2,198 arrays pooled.
- Mutual information and correlation are interchangeable on single-platform M3D.
- Stouffer and Euclidean combination are equivalent; B = 100 permutations suffice.
- The packed MI→CLR path uses 5× less peak memory at 8,563 genes.
- Status: **worked**.

## (d) What must not be claimed

- **Priority.** Do not write "CLR was the first attention mechanism": kernel regression (1964) and fast weights (1991) came first.
  - "First widely adopted attention-like operator in computational biology" is the supported claim. Keep "to our knowledge" in the main text.
  - Reviewers may cite topological overlap (Ravasz et al. 2002; WGCNA, Zhang & Horvath 2005), which normalizes shared-neighbour counts by connectivity and predates CLR. Distinguish it in the Introduction.
- **Per-regulon superiority** over CLR, GENIE3 or PIDC. Also not "improves most regulons" on M3D (54%, p = 0.11).
- **Better GRN inference or direct-edge detection.** CLR and GENIE3 are better on direct TF→target edges. GENIE3 wins per-TF AUROC and median coherence. "State of the art" is out.
- **Robustness or superiority outside de novo, unsupervised inference.** That means against prior-informed or supervised methods, or "regression methods are worse".
- **Coherence gains at all evidence levels** (not significant), and **coherence gains against GENIE3**.
- **Hub-controlled leads at all evidence levels** over GENIE3 or PIDC. That lead holds only against CLR, and only on EcoMAC.
- **Human results.**
  - No CollecTRI advantage. No "large" human effects.
  - Do not quote the raw +0.008/+0.009 K562 gains without the hub caveat.
  - Do not write "human cells" as a general claim: only the K562 essential compendium is in. Genome-wide K562 and RPE1 are pending.
  - State the CORUM and knockdown-response losses plainly.
- **Single cell.**
  - Not "attention helps", and not "wider contexts help".
  - Do not call CLR-NS "the best single-cell method" beyond the methods tested.
  - No p-values computed over 72 combinations as if they were independent.
- **EcoMAC as an independent replication.** Its overlap with M3D arrays is unchecked. The held-out-depth and bootstrap analysis has not been run on EcoMAC, and its +0.032 is a point estimate.
- **Depth as a resolution parameter**, discrete scales or plateaus (multiscale pilot: did not work), or the size–depth relation as biology. That relation is a generic diffusion property, and large regulons' extra preference for depth is convergence to hubs.
- **"Calibrated" or "FDR-controlled" attention weights.**
  - The permutation null calibrates CLR scores; iterated attention mass has no null.
  - Hard FDR operators fail as operators, and the FDR bootstrap is still running (job 100).
  - The current `paper.md` Discussion ("attention can be thresholded at a controlled false-discovery rate") overstates this.
- **Estimated attention beats learned attention.** This is pending; see (e1).
- **Perturbation diversity drives performance.** Single-condition K562 control cells scored about as well as the compendium (mixed).
- **Perturbation-target results.** These belong to Paper 2 and are mixed. The current `paper.md` abstract still includes them; remove them there.
- **Mutual information as essential.** The dependence measure matters little on single-platform data.
- **Removal of global components as part of the method.** It was chosen after inspection and is not in the main pipeline.
- **Components that did not work or were not run**, presented as contributions:
  - smarter keff;
  - multi-head attention (not run);
  - value transforms;
  - diffused-profile readout.

## (e) Pending results that would change the abstract

1. **Learned-attention comparator (plan §5, not started; `docs/plans/scfm_comparator.md`).** scPRINT is primary and scGPT secondary; about 2 weeks and under 2 GPU-days; to be preregistered before scoring. Its result fills the [PENDING] sentence:
   - **Foundation-model attention is below CLR-NS and CLR attention:** "Zero-shot attention from single-cell foundation models falls short of both." The calibrated, label-free baseline then becomes a headline point.
   - **The CLR-on-learned-attention hybrid improves the foundation model** (the thesis test): add "and CLR's normalization improves it". This is the strongest Nature Machine Intelligence story: the old operator's key step transfers to learned attention.
   - **Foundation-model attention matches or exceeds CLR-NS:** "matches" or "exceeds both". The closing framing moves from yardstick to complement.
   - **Statistics problem to fix before preregistration.**
     - With n = 7 BEELINE blocks, an exact two-sided Wilcoxon or sign test cannot go below p = 0.0156. That is above the first Holm threshold for the planned four tests (0.05/4 = 0.0125), so as specified no comparison can reach Holm significance.
     - Options: add the K562 single-cell block (n = 8 gives a minimum p of 0.0078), use one-sided tests in a preregistered direction, shrink the family, or make non-inferiority plus consistency across datasets the primary endpoint.
     - Write the [PENDING] sentence to match ("in 7 of 7 datasets"; "within the non-inferiority margin").
2. **Multi-head attention (keff 3/10/30, preregistered).**
   - If it beats keff = 10 on CORUM without losing on regulons (the decision rule), the trade-off sentence becomes: "multi-head attention, another language-model component, recovers tight complexes without losing regulons". This would be the largest possible upgrade to the abstract.
   - If not, the trade-off sentence stands as written.
3. **Other keff variants (O1, O2, O3).** The headline changes only under the preregistered rule.
4. **BEELINE baselines (GENIE3, GRNBoost2, PIDC; job 0890, queued).** These decide whether "normal-score CLR is the strongest method tested on single-cell data" can be said. The PI's summary "CLR-NS best" depends on this.
5. **BEELINE nested cell bootstrap and per-dataset blocked tests.** The bootstrap is job 110, about 29% done; the blocked tests are plan §3, not started. Together they give intervals for every single-cell statement.
6. **Macro rerun with PIDC on EcoMAC (`results/macro_aupr_v2`, 2026-10-21; not yet reviewed).**
   - Per-regulon mean AP favours attention over PIDC on EcoMAC at all evidence levels (Holm p = 0.0066; both intervals exclude 0; attention wins 59%).
   - EcoMAC SC (p = 0.14) and both M3D cells (p = 0.20 and 0.23) do not.
   - If QC confirms this, "on par" in the abstract still holds, but the plan's "Holm p ≥ 0.21" needs updating.
   - The hub-controlled EcoMAC PIDC cell (`sota_table.md`: 1.57 against 1.40) is also not yet in the reviewed plan text.
7. **Hub-control follow-ups.** These are: the macro view under hub control, a method-neutral hub stratifier, a cluster bootstrap for the knockdown-response modules and a no-ribosome sensitivity check. They could soften or sharpen the clique-loss sentence.
8. **EcoMAC core analysis plus the M3D/EcoMAC array-overlap audit.** Together these would license "two compendia".
9. **Genome-wide K562 and RPE1 compendia (jobs 0871, 0872).** These turn "K562" into "human cell lines", or show that the effect is specific to K562.
10. **Remaining human baselines.** The GENIE3 and GRNBoost2 K562 rerun (jobs 0878, 0879) and PIDC on K562 complete the human competitor table.
11. **Remaining analyses.** None of these changes the abstract unless it reverses a result:
    - the FDR bootstrap (job 100), which gives intervals for the hard-threshold control;
    - directed DREAM5-style evaluation and ARACNE;
    - the granularity and size–depth follow-ups;
    - learning curves;
    - review by Jim Collins.

**Venue reality check.**
- **The plan's condition for Nature Machine Intelligence:** the multiscale result and competitive baselines. Multiscale did not work (no plateaus), and the baselines are at parity per regulon.
- **`paper_draft_notes.md` condition:** a single-cell gain and a head-to-head with learned attention. There is no single-cell gain, and the head-to-head is pending.
- **So the Nature Machine Intelligence case now rests on items 1 and 2.** Without them, the 200-word Nature Communications version in `abstract_draft.md` is the realistic submission.
