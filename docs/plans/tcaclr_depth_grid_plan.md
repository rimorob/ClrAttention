# tcaCLR depth, width and clipping grid: exploratory, label-free selection first

Written 19 Nov 2026 before any cell is computed. EXPLORATORY: nothing here is a confirmatory test, and nothing here changes the registered HepG2 and Kemmeren results (frozen at t = 9, keff = 10, cap = 50, clipped). Requested by the PI: tcaCLR's walk depth was set for aCLR and never searched for tcaCLR; the PI asked for a log-scale grid, run both clipped and unclipped.

## Settings searched

| Axis | Values | Note |
|---|---|---|
| depth t | 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256 | the grid of the earlier aCLR depth re-test, about half an octave apart |
| width keff | 7, 10, 14, 20, 28, 40, 56 | the earlier width test (preregistration_keff) chose the top of its grid (40); the axis reaches higher |
| candidate cap | 150 for the whole grid | cap 50 limits keff above about 20; one reference cell keeps cap 50 |
| clip in the pairwise score (clipS) | clipped, unclipped | clipped: S = (max(0, z_ij) + max(0, z_ji))/sqrt(2); unclipped: S = (z_ij + z_ji)/sqrt(2); candidates are still the top-cap partners with S > 0 |
| clip at the tcaCLR output (clipO) | clipped, unclipped | the second CLR's row z (over reachable pairs), combined by Stouffer, clipped at 0 or not |
| lazy weight | 0.5, fixed | trades off against t |

Operators: 7 keff values by 2 clipS values = 14; each carries all 16 depths and both clipO values: 448 tcaCLR cells and 224 aCLR cells (clipO does not apply to aCLR). Reference cell: keff 10, cap 50, t 9, clipped both (reproduces the registered implementation; must match the stored scores to within numerical tolerance, else stop).

## Stage 1: label-free selection (development matrices only: E. coli M3D and K562 essential Perturb-seq; no truth is read)

Splits: M3D by laboratory (the stored laboratory splits, the first 8); K562 by random halves of the knockdowns (8 splits, fixed seed). For every cell, for each direction (half 1 to half 2 and back):
- **C1, agreement:** AUPR of half 1's score for the pairs in the top 1% of half 2's CLR score (the criterion used for aCLR).
- **C2, self-reproducibility:** AUPR of half 1's score for the pairs in the top 1% of half 2's score by the same method (tcaCLR against tcaCLR).
- **C3, held-out log-likelihood:** mean of log P1^t[i, j] over the half-2 top-1% CLR pairs, with P1^t row-normalised over reachable entries; reported with and without a BIC-style penalty 0.5 k log n, k = mean reachable entries per row and n = number of scored pairs.
Each criterion's selected cell is the argmax of its median over splits; the cell with its two nearest neighbours in each axis is shown as the plateau. Ties and plateau width are reported.

**Rule (fixed now).** A label-free setting is declared for a dataset only if C1, C2 and C3 select cells within one grid step of each other on every axis. If they disagree, no label-free setting is declared and the disagreement is the finding. A setting that differs from (keff 10, t 9, clipped) by half an octave or more in depth or width is reported as such. Nothing in stage 2 can override stage 1.

## Stage 2: descriptive curves on truths already used up

Pooled co-membership AUPR (no jackknife; the plan's ap_fast tie convention) for every cell, on: K562 CORUM, the yeast net4 Complex Portal truth, the HepG2 Complex Portal human truth (Arm A) and the Kemmeren CYC2008 truth (Arm B); plus M3D TF regulons. These truths are spent, so the curves are descriptive; the truth-chosen best cell is shown as an upper bound and is never used for a claim, and a depth or width picked here is a selection. Report: AUPR over the grid as heatmaps (depth by keff, one per clip combination), the registered cell, the best cell, the label-free cell from stage 1 (when declared), and the difference between the label-free cell and the registered cell with the gene-jackknife interval of the existing scorer for those two cells only.

## Readings written in advance

- tcaCLR's depth t = 9 is "about right" for a dataset if the declared label-free depth is within half an octave of 9; "not right" otherwise (the direction and size are reported).
- Unclipped is "better" only if it is higher on stage 1 criteria in both development matrices and not lower on stage 2 curves in more than one truth; otherwise the clipped form stays. The two clip switches are read separately.
- Any claim about a changed setting needs a truth that has not been used; none exists in the project now, so this study can only propose a registered setting for a future test.

## Cost and safety

About 14 operators by 16 depths on matrices of 4,297 to 9,622 genes; dense G by G products. Declared per job at the scheduler: at most 40 cores and 60 GB; stage 1 and stage 2 are separate jobs. Seeds fixed in the scripts. Scripts are committed before the first full run; the registered scorer is not modified.
