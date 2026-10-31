# rCLR on K562 essential Perturb-seq, scored on CORUM (EXPLORATORY)

Exploratory pathfinding, not preregistered, not held-out (K562 CORUM was used in the keff test), nothing confirmatory; not for paper/ or paper tables. Script: `analysis/rclr_k562/rclr_k562.R`. Endpoint list and operator list fixed before running.

Data: 5671 genes x 707 knockdown profiles (half A). Ribosome = the 8 RP-dominated CORUM complexes (6 of them have >= 20 genes). Metric: pooled co-membership AUPR; enrichment = AUPR / base rate.

## Ratio of AUPR (= enrichment ratio) to K10 on CLR

| endpoint | CLR | rCLR | K10 on CLR | K10 on rCLR | K25 on CLR | K25 on rCLR |
|---|---|---|---|---|---|---|
| all CORUM (base 0.0133) | 1.82 | 1.73 | 1.00 | 1.08 | 1.43 | 1.46 |
| CORUM excl. ribosome (base 0.0119) | 1.07 | 1.25 | 1.00 | 1.11 | 1.07 | 1.15 |
| ribosome only (base 0.7502) | 1.02 | 1.02 | 1.00 | 1.01 | 1.01 | 1.02 |
| small (3-6 genes) (base 0.0033) | 1.07 | 1.11 | 1.00 | 1.10 | 1.00 | 1.09 |
| large (>=20 genes) (base 0.1263) | 1.19 | 1.24 | 1.00 | 1.08 | 1.08 | 1.19 |

## Enrichment over base rate

| endpoint | CLR | rCLR | K10 on CLR | K10 on rCLR | K25 on CLR | K25 on rCLR |
|---|---|---|---|---|---|---|
| all CORUM | 6.39 | 6.07 | 3.51 | 3.78 | 5.02 | 5.12 |
| CORUM excl. ribosome | 2.18 | 2.57 | 2.05 | 2.27 | 2.20 | 2.36 |
| ribosome only | 1.21 | 1.21 | 1.18 | 1.20 | 1.20 | 1.21 |
| small (3-6 genes) | 2.16 | 2.23 | 2.01 | 2.22 | 2.02 | 2.19 |
| large (>=20 genes) | 2.29 | 2.39 | 1.93 | 2.08 | 2.09 | 2.29 |

## Half-of-genes jackknife of the ratio to K10 on CLR (median [2.5%, 97.5%], 100 reps; descriptive)

| endpoint | CLR | rCLR | K10 on rCLR | K25 on CLR | K25 on rCLR |
|---|---|---|---|---|---|
| all CORUM | 1.78 [1.57, 2.10] | 1.68 [1.51, 1.99] | 1.08 [0.95, 1.25] | 1.42 [1.31, 1.53] | 1.44 [1.29, 1.69] |
| CORUM excl. ribosome | 1.06 [0.97, 1.17] | 1.23 [1.10, 1.42] | 1.10 [1.00, 1.19] | 1.06 [1.00, 1.14] | 1.14 [1.04, 1.22] |
| ribosome only | 1.02 [1.00, 1.06] | 1.03 [1.00, 1.06] | 1.01 [1.00, 1.02] | 1.01 [1.00, 1.02] | 1.02 [1.01, 1.04] |
| small (3-6 genes) | 1.15 [0.75, 1.56] | 1.15 [0.77, 1.54] | 1.10 [0.93, 1.32] | 1.00 [0.85, 1.13] | 1.09 [0.92, 1.26] |
| large (>=20 genes) | 1.19 [1.12, 1.27] | 1.24 [1.17, 1.32] | 1.09 [1.03, 1.14] | 1.08 [1.04, 1.11] | 1.19 [1.14, 1.25] |

## Label-free: reproducible partners (top-50 positive candidates shared by two halves of the profiles)

| operator | frac_zero_shared | frac_ge_k0 | k0 | median_shared | q10 | q90 |
|---|---|---|---|---|---|---|
| CLR | 0.1008640 | 0.6882384 | 3 | 6 | 0 | 23 |
| rCLR | 0.1068595 | 0.6621407 | 3 | 6 | 0 | 23 |

## Sanity checks

| check | value |
|---|---|
| documented CLR CORUM AUPR / enrichment | 0.0847 / 6.39 |
| documented attention (t = 9, K10) CORUM AUPR / enrichment | 0.0465 / 3.51 |
| cor(CLR, rCLR) over upper-triangle pairs | 0.7942 |

