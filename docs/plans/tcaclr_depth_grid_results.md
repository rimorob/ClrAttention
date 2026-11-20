# tcaCLR depth/width/clipping grid and adaptive-setting probes: results (exploratory)

Run 19-20 Nov 2026 under `docs/plans/tcaclr_depth_grid_plan.md` (Revision 2). All 11 jobs finished rc=0. Nothing here is confirmatory: every truth in the project is spent, and the registered HepG2 and Kemmeren results (t = 9, keff = 10, cap 50, clipped) are unchanged. Outputs: `results/tcaclr_depth/` (stage1_*, stage1_summary, stage2_*, stage2_summary, adaptive_m3d, adaptive_k562, adaptive_planted_m3d).

## Verdict under the pre-written rule (Revision 2, item 9): no label-free setting is declared

| Requirement | Result |
|---|---|
| Planted-simulation regret <= 5% in >= 80% of configurations | Not met. Detectable datasets only (oracle AUPR >= 5x base, 53-54 of 60): the best rule, the edge rule at eps = 0.01, is within 5% in 65-76% of configurations (tcaCLR 70% at keff 10, 65% at keff 28; aCLR 70% and 76%). Fixed t = 9 is within 5% in 9-20%. |
| Selected scale monotone in module size (Spearman >= 0.7) | Not met, and not attainable: the edge rule does not depend on module size by construction, the likelihood rules give 0.39-0.54, and the oracle itself gives 0.61 (tcaCLR) and 0.48 (aCLR). |
| Permuted negative controls flagged in >= 95% | Edge rule: 95% (19 of 20 per width, the nominal rate of the leave-one-out edge). Graph-signal likelihood, strict gain <= 0: 32% (M3D) and 42% (K562), because permuted gains are positive at about 1e-6; with the 1-SE flag 83% and 92%. Combined flag fails. |
| Split stability within one grid step in >= 7 of 8 splits | Likelihood depth: met at 15 of 16 width settings; K562 keff 20 is 4 of 8. |

## What the label-free readings say (descriptive)

- **Spectral edge.** Many observed eigenvalues of the walk operator sit above the permuted-data edge: 73-101 (M3D) and 10-25 (K562) depending on width. The permuted edge is 0.91 to 0.75 (M3D) and 0.92 to 0.76 (K562) as keff goes 7 to 80, so the depth t* = ln(0.05)/ln(edge) falls from 31 to 10 (M3D) and 38 to 11 (K562); at keff 10 it is 24.5 (M3D) and 29.3 (K562), at keff 28 it is 14.7 and 16.3. Depth and width trade off along a ridge; the registered pair (keff 10, t 9) is on the shallow side of it, by about 1.4 octaves at keff 10.
- **Graph-signal likelihood.** Best depth M3D 24 at keff <= 14 and 18 at keff >= 20; K562 48 at keff 7, 32 at 10-14, 12-13 at keff >= 28. The one-SE rule paired across splits almost never moves the choice; the marginal-SE version is much shallower (M3D 9-16, K562 4-8). N_eff is 12.3 (M3D) and 39.5 (K562), so these are noisy and the marginal-SE column is not clearly better.
- **Stage 1 criteria disagree on width.** Held-out likelihood (C3) picks the top of the grid (keff 80) on both matrices, with t = 4 (M3D) and 9 (K562). Null-adjusted self-reproducibility (C2_adj) picks keff 56, t 16 (M3D) and keff 7, t 48 (K562). Width is not identified label-free; this repeats the earlier width study.
- **Clipping.** Unclipping the second CLR's output (clipO) lowers the null-adjusted reproducibility in most cells (higher in 4% of K562 cells, 14% of M3D cells); unclipping the first score (clipS) raises the likelihood in every cell by about 0.001-0.002 and moves the other criteria by under 0.001. So clipO unclipped is worse on the label-free criteria; clipS is a wash on AUPR and slightly better on likelihood.

## Planted simulation

On the planted modules (permuted M3D background, 3,000 genes, modules of 5-100 genes), the best-cell depth for tcaCLR (best over the two widths) is 32 (median) for modules of 10-100 genes and 9 for 5-gene modules, and fixed t = 9 loses a median 38-48% of the oracle AUPR. The edge rule's chosen depth tracks the oracle (median regret about 0). The diffusion-PMI alternative has much larger regret for every rule (median 0.46-0.92) and was not compared with tcaCLR on AUPR; it stays an exploratory read-out. Caveat for reading this: planted modules are clean, equal-dose, non-overlapping cliques without hubs, which favour deep smoothing; the simulation shows the edge rule can find the planted depth, not that real data want that depth.

## Spent-truth curves (stage 2, descriptive only)

- tcaCLR AUPR at keff 10, clipped, is flat in depth on the real truths. From t = 4 to 32 it moves by under 10% of its value except on K562 CORUM and net4 at t = 4 (down 18% and 30% from t = 9); t = 1 collapses everywhere. Depth 24 relative to depth 9: K562 CORUM -1.3%, net4 -2.8%, Arm A -3.7%, Arm B +0.7%, M3D SC regulons +3.8%. So the edge-rule depth (24-29 at keff 10) is not clearly better or worse than 9 on any truth.
- Truth-chosen best cells (the best of 819, an optimistic upper bound) beat the registered tcaCLR cell by 0.006-0.014 AUPR, 4-10% relative (K562 CORUM +10%, net4 +6%, Arm A +9%, Arm B +4%, M3D SC +7%). The best depths scatter (4, 12, 18, 18, 48, 96) and the best widths scatter (7-40). There is no consistent optimum.
- aCLR is different: its best cells beat its registered cell by 14-69% relative (K562 CORUM +69%, net4 +69%, Arm A +29%, Arm B +14%), with best widths 56-80 and unclipped first score on most truths. The registered aCLR setting (keff 10, t 9, cap 50) is clearly not near its optimum. This is a selection among 819 cells on spent truths, so it is a hypothesis for a fresh test only.

## Reading

1. tcaCLR's t = 9 is not wrong in the sense that matters for the results already reported: on five spent truths the AUPR surface is a broad plateau, and any gain from retuning depth or width is at most a few percent to about 10% as an upper bound.
2. The simulation says depth matters a lot, the real truths say it barely does. That disagreement is the main finding of the probe; the planted structure is too clean to calibrate real depth.
3. No automatic rule passed the pre-written bar. The edge rule is the best candidate (it finds the planted depth and its false alarm rate is nominal) but its prediction on real data is only descriptive.
4. If a retuned setting is ever to be claimed, it needs a registered rule and an unused truth. The cheapest such candidate is the aCLR retune (wide keff, large effect), not tcaCLR.
