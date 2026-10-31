# rCLR on K562 essential Perturb-seq / CORUM: EXPLORATORY pathfinding

Not preregistered, nothing confirmatory, not for paper/ or paper tables. K562 CORUM was already used in the keff test and in analysis/corum_exploratory/, so this is not a held-out dataset. Operators and the five endpoints were fixed before running (script header: `analysis/rclr_k562/rclr_k562.R`); nothing was added after seeing results.

- `truthgenes/` (primary): 5,671 genes, the data behind the earlier CLR 6.39x / attention 3.51x numbers (reproduced exactly: 0.0847 / 0.0465 AUPR, base rate 0.0133).
- `all/` (secondary): 8,559 genes, the keff-test data; cached CLR matches the recomputation exactly (max diff 0) and K10 t = 9 reproduces the keff-test CORUM AUPR 0.0493 (3.72x).
- Each folder: `summary.md` (ratios, enrichments, jackknife spread, label-free check, sanity), `endpoint_metrics.csv`, `label_free_repro.csv`, `sanity_checks.csv`.
- Ribosome = the 8 RP-dominated CORUM complexes (>= 30% RP genes) of analysis/corum_exploratory/corum_complexes.csv; "excl. ribosome" drops them; small = 3-6 genes, large = >= 20 genes (6 of the 28 large complexes are RP-dominated). Endpoint pair universe = gene pairs within the selected complexes' genes, as CORUM_noRP did before.

Ratio of AUPR (= enrichment ratio) to K10 on CLR, primary universe:

| endpoint (base rate) | CLR | rCLR | K10 CLR | K10 rCLR | K25 CLR | K25 rCLR |
|---|---|---|---|---|---|---|
| all CORUM (0.0133) | 1.82 | 1.73 | 1.00 | 1.08 | 1.43 | 1.46 |
| CORUM excl. ribosome (0.0119) | 1.07 | 1.25 | 1.00 | 1.11 | 1.07 | 1.15 |
| ribosome only (0.7502) | 1.02 | 1.02 | 1.00 | 1.01 | 1.01 | 1.02 |
| small 3-6 (0.0033) | 1.07 | 1.11 | 1.00 | 1.10 | 1.00 | 1.09 |
| large >= 20 (0.1263) | 1.19 | 1.24 | 1.00 | 1.08 | 1.08 | 1.19 |
