# Yeast reciprocity preregistration: STOPPED again at protocol step 3 (orientation floor)

*1 Nov 2026, after deviation 1 (743d24e). No truth has been downloaded or opened; yeast is unconsumed.*

## Result (job 0984)

- **Dye check part 1 holds** (1,387 mutant hybridizations). Balanced mutant groups: 685 of 700 (`dye_balance_GSE42527.csv`).
- **Orientation reproduction** (first 10 balanced mutants in C-locale order: abf2, ace2, aco1, ada2, add66, aft1, aim1, apc9, apn2, apq12; `implementation_check.csv`):

  | Rule | Median r |
  |---|---|
  | as deposited | 0.39 |
  | flip_test_ch1 | −0.80 |
  | **flip_ref_ch1 (best)** | **0.793** |

  The best rule's median r is below the preregistered floor of 0.8, so the status is `F_NOT_COMPUTABLE_verdict_1` (`orientation_rule.txt`). **By the fixed preregistration this is verdict 1: gate F is not computable, nothing else runs, and no truth is consumed.**

## Why (label-free metadata diagnosis)

- Every GSE42527 hybridization has channel 1 = Cy5 (1,387/1,387).
- By each sample's own VALUE definition, every deposited value is already log2(mutant/ref):
  - reference in channel 1 → "log2 ratio (Cy3/Cy5)";
  - mutant in channel 1 → "log2 ratio (Cy5/Cy3)".
- Yet the published M agrees with "as deposited" only for some mutants (ada2 r 0.98 and add66 0.93, mostly from the older sub-series). For most `[hs1991]` mutants (ace2, aco1, aft1, apc9, apn2, apq12: r 0.76–0.97) it agrees only after negating the hybridization with the reference in channel 1. abf2 and aim1 fit neither (r ≤ 0.78).
- So the deposited labels or definitions and the published M disagree for part of the series, and no single global rule fits. The preregistered candidate set is global.
- **Wild types (GSE45115):** the 755 VALUE definitions are mutually consistent with log2(wt/ref) under "as deposited" (channel 1 = Cy5 throughout; definitions "(Cy3/Cy5)" with the reference in channel 1, "(Cy5/Cy3)" with wt in channel 1, "test/ref", and "−INV_VALUE … test/ref"). But the mutant discrepancy means the deposited definitions cannot be taken on trust.

## Not done, because of the stop

The wild-type matrix, splits, clusters, gate F and the selections were not run. The **HepG2/Jurkat download (GSE264667, about 22 GB) was not started**, because under verdict 1 from gate F no human arm runs.

## Options for the PI (any change would be deviation 2, before any truth)

- **(a)** Accept verdict 1 for the study. Reciprocity stays a synthetic-only candidate, H1 is not run, and yeast stays unconsumed.
- **(b)** Replace gate F's reliance on rebuilt profiles. Run gate F on the published M (mutants) against wild-type profiles oriented "as deposited" (consistent by their definitions), with the orientation question disclosed. The confirmatory analysis already uses the published M only.
- **(c)** Add a per-hybridization orientation rule derived from sub-series membership, validated on a fixed held-out set of balanced mutants with the same floor.

Option (b) is the smallest change, but it removes the only check that the wild-type orientation is right.
