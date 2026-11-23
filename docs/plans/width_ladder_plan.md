# tcaCLR width ladder and multi-scale scores: exploratory, descriptive

Written 22 Nov 2026 before any real score of the ladder exists. EXPLORATORY and DESCRIPTIVE: every truth used here is spent, nothing here is a test, and no setting of the registered tcaCLR (keff 10, cap 50, t = 9, clipped) or tcaCLR2 changes. Any single width that looks best on these truths is a selection and needs a fresh test; the multi-scale scores below have no selection step. Requested by the PI: the earlier width grid (keff 7 to 80) had its regulon optimum at its lower edge (7), so the ladder covers both ends and all of it is read, with no pick-one-width step.

## Settings

| Axis | Values |
|---|---|
| keff ladder (half-octave steps) | 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128 (13 widths) |
| candidate cap | 4 x keff, at most G - 1 (G = number of genes) |
| lazy weight | 0.5; clipS and clipO TRUE (registered clipping); Stouffer, symmetric read-out |
| depth rule R9 | fixed t = 9 (tcaCLR style) |
| depth rule RE | t from the permutation-edge rule at that width (tcaCLR2 style): edge = mean over 10 per-gene permutations (seeds AD$seed + 7000 + p, as registered) of the largest non-trivial eigenvalue of the symmetrised lazy operator of THAT width and cap; t = max(1, round(ln 0.05 / ln edge)); the edge is recomputed per width on the same ten permuted matrices (one fit per permutation, all widths and the reference operator read from it) |

Width 10 is not on the ladder and the registered cap (50) is not 4 x keff (= 40), so the registered cells are reproduced by a separate reference operator (keff 10, cap 50), used only as the gate and as the comparators tcaCLR (t = 9) and tcaCLR2 (t = t2 of the registered edge file).

Fixed here, before any score (implementation choices that the request did not spell out):
1. If the permuted operator of a width has edge >= 1 - 1e-6 (disconnected, no decay) or ln 0.05 / ln edge > 256, t is set to 256 and the width is flagged `capped` in the t table (the earlier grid's maximum depth is 256). Nothing else is altered.
2. Multi-scale scores use the pair scores over ALL gene pairs (upper triangle). Rank percentile = average-tie rank / number of pairs (the clipped scores have a large tie block at 0). The within-width z (MS2) uses the mean and sd over all pairs.
3. The hub-controlled scope uses ONE stratifier for all cells of a dataset: H = q_i + q_j with q the stationary distribution of the registered operator (keff 10, cap 50), exactly as the scoreboard and the direct-edge diagnostic; it is not recomputed per width. Strata: 10 quantile strata of H (`rz_hub_stratify`, within-stratum percentile).
4. Best width = argmax of the metric over the 13 single widths of a rule (first in case of an exact tie). It is "at an edge" when it is keff 2 or keff 128; "next to an edge" (3 or 96) is also listed. Depth rules are read separately.
5. A missing or failing gate stops the run; the t table, gate tables and ledger are written before any heatmap.

## Multi-scale scores (no tuning, no selection), for each depth rule

| Score | Definition |
|---|---|
| MS1 | mean over the 13 widths of each pair's rank percentile (within the width's score over all pairs) |
| MS2 | mean over the 13 widths of the pair's score standardised within the width (z over all pairs) |
| MS3 | as MS1 over the 8 widths 4, 6, 8, 12, 16, 24, 32, 48 only |

Compared with each single width, CLR (the registered MI-CLR score), tcaCLR (keff 10, cap 50, t 9) and tcaCLR2 (keff 10, cap 50, t2). Differences and ratios to these are shown; no interval, no p-value.

## Datasets, truths, metrics (the stage-2 truth set plus the direct-edge metrics; all spent)

| Dataset | Matrix | Metrics per score |
|---|---|---|
| M3D (4297 genes) and EcoMAC (4189 genes, ComBat) | scoreboard matrices (`sb_load_ecoli`) | pooled AUPR / base rate (raw and hub-controlled) on RegulonDB regulons SC, all and Complex Portal complexes; direct TF-to-target edges (SC, all): pooled AUPR ratio, median per-TF AUROC, hub-controlled pooled ratio (as `analysis/direct_edge_diag`) |
| K562 pseudobulk (5671 genes) | ideas47 set-up | pooled AUPR / base rate (raw and hub-controlled): CORUM, DoRothEA A-B, A-B-C, CollecTRI, KD response |
| yeast net4 (5950), Kemmeren (6169; CYC2008, Arm B), HepG2 (9622; Complex Portal, Arm A) | stage-2 matrices and frozen truth files (hash-guarded) | pooled AUPR / base rate (raw) |

All truth builders, the pooled scorer (`ap_fast`), the per-TF AUROC (`sb_per_tf_median_auroc`), the hub stratifier, the operator builder `gl_make_P`, the lazy walk `gl_walk`, the read-out `gl_zmat` and the edge machinery (`ad_fit`, `ad_permute_genes`, `ad_Psym`, `ad_eigs`, `ed_t2`) are the existing ones; the new code only loops over widths, accumulates the multi-scale scores and collects tables.

## Gate (spent truths, gate only)

Before any ladder cell is accepted, the reference operator through the same code must reproduce, to 1e-9: the registered tcaCLR (keff 10, cap 50, t = 9) cells of `results/tcaclr2_scoreboard/cells_reference_tcaclr_all.csv` (every pooled raw and hub-controlled value, tf_aupr_ratio, tf_auroc_median), the direct-edge numbers of `results/direct_edge_diag/gate_direct_edge.csv` and the hub-controlled direct-edge ratios of `variants_{m3d,ecomac}.csv`; and, with the stored t2, the tcaCLR2 cells of `cells_tcaclr2_all.csv`. The edge of the reference operator recomputed by the ladder edge job must equal the stored edge file (edge to 1e-6, t2 identical). Failure stops the run.

## Outputs (`results/width_ladder/`)

Per dataset: the t table (edge and t per width), the metric of every score and truth (long CSV, one row per setting x truth x scope x metric), the gate tables. Assembled: the response-surface table (width by metric, per rule) with the best width per truth, the edge flags, the ledger below, and descriptive heatmaps and curves (PNG).

## Ledger and multiplicity statement

Every setting scored is a row of the ledger table (dataset, rule, score, n widths, number of metric cells) written by the assemble step. Settings per dataset: 13 widths x 2 rules = 26 single-width settings and 3 multi-scale scores x 2 rules = 6, i.e. 32, plus three comparators (CLR, tcaCLR, tcaCLR2) that are not settings of the ladder. Metric cells: 12 per setting on M3D and EcoMAC (6 pooled + 6 direct-edge), 10 on K562, 1 on each of yeast, Kemmeren, HepG2: 32 x (12 + 12 + 10 + 3) = 1184 setting-by-metric cells. No correction is applied because nothing is tested; the cells are correlated (same pair scores, overlapping truths), the best-width locations are the maximum of 13 values per cell, and a "best" is optimistically biased on spent truths. Nothing here may be quoted as an accuracy result for a width, only as the shape of the response surface.
