# Local vs global temperature: real tension or evaluation artefact? (PI memo; exploratory, all truths spent)

**Short answer.** Mostly an artefact of the two evaluation universes, plus one small, real operator property. The two metrics give different labels to the very pairs on which the operators differ. Pooled TF-to-gene AUPR is a head-of-ranking metric, and one TF dominates its head. Per TF and after hub control, the global advantage mostly disappears. **Recommendation:** keep the registered local operator as the method. Treat V2 as a diagnostic that showed what the per-gene temperature throws away, not as a second method. Do not tune any hybrid on these truths. Run three cheap diagnostics (below). Settle the TF-edge question with one pre-specified test on a fresh TF-target truth.

## 1. Diagnosis

**What differs between the universes.**

| | (a) pooled co-membership | (b) TF-to-gene edges |
|---|---|---|
| Positives | pairs that share a regulator (regulator-agnostic regulons: TF, sigma, sRNA, small molecule; the regulator is never a node) or a complex; same-operon pairs removed | the direct TF-target pair itself |
| Pair type of a positive | almost all target-target (co-regulated, two hops through a latent regulator) | regulator-target (one hop) |
| Rows that matter | every annotated gene, both ends | one end is always one of about 200 TFs; non-TF rows enter only as columns |
| Weighting | a regulon or complex of size n contributes about n^2/2 positives, so big modules dominate quadratically | a TF contributes n positives, but AUPR is decided at the head, so the few TFs with the sharpest rows dominate |
| Base rate | much higher (regulon truths 1.4-1.6x is a small lift over a high base) | 0.5-0.6% |

**Mechanism.** The per-gene (local) temperature gives every row exactly 10 effective partners. It erases cross-row magnitude: the Spearman correlation of row-max CLR with row-max weight is 0.20 under local and 0.77 under V2 (recipe ladder, M3D). After the walk, the second CLR row-standardises rows of similar shape. The ranking is row-balanced and propagates through shared neighbours, which is exactly what target-target co-membership rewards.

The global temperature lets high-signal rows collapse to 1-2 partners and lets diffuse rows spread to 26-46. Sharp rows give very large second-CLR z for their top pair, which then rises to the head of a pooled ranking. Those pairs are disproportionately direct edges (lexA-like). Diffuse rows go flat and sink. This helps (b), whose head wants confident one-hop pairs. It hurts (a): a sharp gene's walk stays on its top partner, often an operon mate that (a) excludes, and never reaches its co-regulated genes.

**The decisive existing fact.** Under local tcaCLR the TF-edge head is monopolised by FliZ (b1921): 52 of the top 100 pairs on M3D SC, 3 distinct TFs in the top 100 (V2: 33 TFs, b1921 share 0.32). Of FliZ's 30 highest tcaCLR pairs, 0 are true edges, and 14 (M3D) and 17 (EcoMAC) are two-hop pairs in the flagellar module. These are co-regulated pairs: positives or near-positives in (a), negatives in (b). (Caveat from the project lead: the removal figures below are post-hoc influence diagnostics, picked after seeing the head of the ranking. They are never to be reported as corrected scores. The follow-up is a symmetric leave-one-TF-out over all TFs and a TF-block bootstrap, not a selective removal.) Removing b1921 from the truth lifts local tcaCLR from 5.84 to **7.80** (M3D SC; V2 with b1921 is 7.82) and from 5.56 to **7.03** (EcoMAC SC; V2 6.66). Removing one random TF moves tcaCLR between 4.76 and 6.73 on M3D SC. That range is as wide as the whole V2-local gap (1.98). Other evidence points the same way:
- Median per-TF AUROC favours local: 0.627 vs 0.605 (M3D SC) and 0.599 vs 0.573 (EcoMAC SC).
- On EcoMAC, hub control removes the gain: hub-controlled TF-edge ratio V2 5.45 vs local 5.41 (raw 6.66 vs 5.56). It survives on M3D (6.92 vs 5.76).

**The grid-edge results are not evidence for a global temperature.** At log2 tau = +2 the softmax is nearly flat over each gene's top 50. Effective partners per row are 32.7-49.6 (10th-90th percentile) on M3D and 47.7-49.8 on K562, so every row has about the same breadth again, which is the local regime. The truth-chosen best tau at +2 on the complex rows (K562 CORUM, yeast, Kemmeren, M3D) says that a **flat top-50 neighbourhood** helps complexes. It also beats local keff 48 with cap 192 (Kemmeren 19.6 vs 15.6; K562 CORUM 9.12 vs 8.06), which points to the hard cap of 50, i.e. truncation, not temperature. At V2 itself (heterogeneous breadth), local wins complexes on 5 of 6 datasets (all except HepG2) and wins every regulon co-membership cell.

**Testable predictions.** Each needs only the V2 pair scores rebuilt (a few minutes per dataset with `gt_lib.R`) and rescoring; no new truth, no long jobs.
- **P1. Head artefact.** Score TF-to-gene AUPR for V2 and local with b1921 removed, plus the leave-one-TF-out distribution for both. If the artefact account holds, the V2/local ratio falls from 1.34/1.20 to within 0.9-1.1. A ratio still above 1.15 with b1921 removed means the global effect is broad. This is the critical diagnostic.
- **P2. Sharp rows carry the gain.** Split TF-to-gene AUPR, and per-TF AP/base rate, by tertile of the TF row's effective partners under V2. Prediction: the V2 gain sits in the low-partner tertile, with a loss in the high-partner tertile.
- **P3. Pair-type split.** Split (a) into pairs that involve a TF and pairs that do not, and read the (a) label of local's top-30 b1921 pairs. Prediction: local's advantage on (a) comes from non-TF (target-target) pairs, and FliZ's head pairs are mostly (a)-positives. The label lookup alone can be done from existing files plus the stored tcaCLR scores.

## 2. Can it be resolved?

The two metrics ask different questions ("is this the direct regulator?" vs "are these co-regulated?"). No operator has to win both, and any knob tuned to win both on these truths is a selection on spent data.
- **Hybrid temperature (per-gene tau shrunk toward a global value):** not worth it. It adds a mixing weight that could only be chosen on these truths, and the two endpoints are preferred by different metrics, so the "best" weight would depend on which metric is chosen.
- **Label-free rule to choose local or global:** not possible in principle. The choice depends on the user's question (direct vs co-regulated), not on the data. No unlabelled statistic can know that.
- **Score fusion or rank averaging:** not worth it. The multi-temperature averages (MS1-MS3) are almost never best, and in the width ladder they were 14-26% worse on direct edges (the head gets diluted). A two-operator average will land in between on both metrics.
- **Task-conditional use (local without a TF prior, global with one):** the only candidate worth keeping. Keep it only if P1 shows the gain is broad, and claim it only after the fresh test below. Until then it is a two-operator product built on a selection (V2's tau rule was chosen after seeing M3D/EcoMAC direct edges).
- **Extending the grid above +2:** do not. The limit as tau goes to infinity is a known single operator: uniform weights over the top-50 positive CLR partners (an unweighted kNN-50 graph). At most, score that one limit as a descriptive point. Its question is truncation and cap, not temperature.

## 3. What the paper may say now

**Allowed (exploratory section, descriptive wording):** "Replacing the per-gene temperature by one global temperature, a variant chosen after inspecting the E. coli direct-edge results, ranks known TF-target edges higher in a pooled ranking on both E. coli compendia (1.20-1.34x the registered operator). The registered per-gene temperature ranks co-regulated and co-complex gene pairs higher on 5 of 6 datasets. Per TF and after hub control, the two operators are similar. Much of the pooled difference comes from a single transcription factor (FliZ) whose co-regulated module occupies the head of the local ranking. All results use benchmark truths that had already been used, and no interval or test is reported."

**Not allowed:**
- that a global temperature is better, or that tcaCLR beats CLR on direct edges;
- any recommended temperature, or any reading of the +2-edge optima as support for a global temperature;
- that the method adapts to the task;
- that tcaCLR2 is competitive on direct edges (it is the worst, 4.39 and 4.74);
- any p-value or confidence statement from these files.

In the main text, use only the registered tcaCLR (t = 9) as the method. Keep tcaCLR2 labelled as not preregistered.

**Weak points to state:**
- Every truth is spent, and V2's tau rule was itself selected on M3D and EcoMAC.
- The TF-edge evidence is one truth family (RegulonDB) on two compendia built from overlapping E. coli arrays: one replication, not two.
- 8 of 14 best-tau rows sit at the grid edge, in the near-uniform regime.
- Jackknife or TF-bootstrap intervals are not shown, and leave-one-TF-out variation (about 2 units on M3D SC) is as large as the effect.
- The tcaCLR2 contrast also mixes in depth (t2 = 25 vs t = 9).

**The single decisive experiment.** Freeze everything first: V2 = one tau giving a mean of 10 effective partners, cap 50, t = 9; comparator = registered tcaCLR; scorers as now. Commit before scoring. Data: a fresh expression compendium with an unused, curated TF-target truth. (Decision of the project lead: the B. subtilis data are held in reserve and are not to be touched or scored until released.) My first choice is B. subtilis (the Nicolas et al. condition compendium with SubtiWiki/DBTBS regulations): a fresh prokaryote, closest to the 2007 CLR setting. The yeast net4 matrix with a TF-target truth never scored in this project is the fallback (verify that it is unused). For the complex half, use DepMap co-essentiality with a complex catalogue not yet used in the project (verify; otherwise Complex Portal, stating that the catalogue was used elsewhere).

**Multiplicity.** Fixed-sequence testing, so no extra correction is needed:
- **H1 (primary):** pooled TF-to-gene AUPR, log(V2 / local) > 0, with a 95% TF-block bootstrap interval. Claim it only if the lower bound is above 0 **and** the ratio stays above 1 after removing the TF with the largest share of either score's top 100 (a pre-written robustness condition, not a second test).
- **H2 (tested only if H1 passes, alpha 0.05):** on DepMap complexes, local is not worse than V2, i.e. pooled AUPR(local) at least 0.95 x V2, with a gene-jackknife interval.
- **Descriptive only, no tests:** per-TF median AUROC and AP, hub-controlled ratios, CLR and tcaCLR2 columns, the P2 tertile split, and the flat top-50 limit.

If H1 fails, the paper keeps one operator, and the V2 result is reported as a diagnostic of what row normalisation costs on direct edges.
