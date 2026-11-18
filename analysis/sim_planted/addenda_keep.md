## Addendum 2 (17 Nov 2026): real K562 reachability check (structure only)

Approved by the PI on 17 Nov 2026 (spent data; no method scored against any truth). Script `analysis/sim_planted/real_reach_k562.R`, job 7201 (107 s, 2.8 GB). Matrix: the scoreboard K562 truth-gene matrix (5,671 genes by 707 profiles); operator at the registered settings (aCLR walk, t = 9); complexes from `corum_sets` (CORUM source, size 3 to 200). For member rows of each complex: reachable off-diagonal entries of P^9, and the share of those that lie inside the complex.

Overall: unreachable share of all pairs 0.28%; no row has fewer than 2 reachable entries; median reachable entries per row 5,654 of 5,670.

| Class | Complexes | Median size | Median reachable per member row | Median in-module share | Smallest in-module share | Closed (reachable < 1.5 x size) |
|---|---|---|---|---|---|---|
| ribosomal (30% or more r-protein genes) | 8 | 31.5 | 5,654 | 0.0054 | 0.00035 | 0 |
| 20 or more genes | 22 | 28.5 | 5,654 | 0.0048 | 0.0032 | 0 |
| 7 to 19 | 285 | 10 | 5,654 | 0.0016 | 0.00091 | 0 |
| 3 to 6 | 870 | 3 | 5,654 | 0.00035 | 0.00024 | 0 |

Reading: in real K562 data no complex row is closed. Every member row reaches almost the whole gene universe and the in-module share is about the module's share of the universe, so the reachable-only row statistic has a full background for ribosomal and other large complexes. The reachability-confinement collapse of the simulation (module plus 2 to 3 genes, in-module share 0.96 to 0.97) does not operate on this matrix and operator. Limits: one matrix (K562, 5,671 genes), one operator setting, a structure-only measurement; it says nothing about HepG2 or Kemmeren, where the same descriptor is reported with the registered results, and it does not explain the K562 ribosome gain.
