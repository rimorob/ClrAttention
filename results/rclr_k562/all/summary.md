# rCLR on K562 essential Perturb-seq, scored on CORUM (EXPLORATORY)

Exploratory pathfinding, not preregistered, not held-out (K562 CORUM was used in the keff test), nothing confirmatory; not for paper/ or paper tables. Script: `analysis/rclr_k562/rclr_k562.R`. Endpoint list and operator list fixed before running.

Data: 8559 genes x 707 knockdown profiles (half A). Ribosome = the 8 RP-dominated CORUM complexes (6 of them have >= 20 genes). Metric: pooled co-membership AUPR; enrichment = AUPR / base rate.

## Ratio of AUPR (= enrichment ratio) to K10 on CLR

| endpoint | CLR | rCLR | K10 on CLR | K10 on rCLR | K25 on CLR | K25 on rCLR |
|---|---|---|---|---|---|---|
| all CORUM (base 0.0133) | 1.82 | 1.66 | 1.00 | 1.01 | 1.41 | 1.44 |
| CORUM excl. ribosome (base 0.0119) | 1.06 | 1.24 | 1.00 | 1.03 | 1.07 | 1.11 |
| ribosome only (base 0.7502) | 1.02 | 1.02 | 1.00 | 1.01 | 1.01 | 1.02 |
| small (3-6 genes) (base 0.0033) | 1.03 | 1.05 | 1.00 | 1.08 | 1.01 | 1.11 |
| large (>=20 genes) (base 0.1263) | 1.20 | 1.25 | 1.00 | 1.09 | 1.08 | 1.20 |

## Enrichment over base rate

| endpoint | CLR | rCLR | K10 on CLR | K10 on rCLR | K25 on CLR | K25 on rCLR |
|---|---|---|---|---|---|---|
| all CORUM | 6.75 | 6.18 | 3.72 | 3.74 | 5.24 | 5.37 |
| CORUM excl. ribosome | 2.20 | 2.57 | 2.08 | 2.13 | 2.22 | 2.30 |
| ribosome only | 1.21 | 1.21 | 1.19 | 1.20 | 1.20 | 1.21 |
| small (3-6 genes) | 2.19 | 2.23 | 2.13 | 2.31 | 2.15 | 2.36 |
| large (>=20 genes) | 2.32 | 2.40 | 1.92 | 2.10 | 2.09 | 2.31 |

## Half-of-genes jackknife of the ratio to K10 on CLR (median [2.5%, 97.5%], 100 reps; descriptive)

| endpoint | CLR | rCLR | K10 on rCLR | K25 on CLR | K25 on rCLR |
|---|---|---|---|---|---|
| all CORUM | 1.76 [1.54, 2.10] | 1.62 [1.44, 1.94] | 1.01 [0.89, 1.15] | 1.39 [1.30, 1.52] | 1.43 [1.28, 1.66] |
| CORUM excl. ribosome | 1.04 [0.94, 1.16] | 1.21 [1.07, 1.39] | 1.02 [0.88, 1.14] | 1.06 [1.00, 1.13] | 1.10 [0.99, 1.20] |
| ribosome only | 1.02 [1.00, 1.06] | 1.02 [1.00, 1.06] | 1.01 [1.00, 1.02] | 1.01 [1.00, 1.02] | 1.02 [1.01, 1.03] |
| small (3-6 genes) | 1.07 [0.72, 1.50] | 1.11 [0.72, 1.51] | 1.08 [0.85, 1.42] | 1.01 [0.80, 1.17] | 1.10 [0.87, 1.31] |
| large (>=20 genes) | 1.20 [1.13, 1.28] | 1.24 [1.18, 1.33] | 1.09 [1.03, 1.16] | 1.08 [1.05, 1.11] | 1.20 [1.15, 1.27] |

## Label-free: reproducible partners (top-50 positive candidates shared by two halves of the profiles)

| operator | frac_zero_shared | frac_ge_k0 | k0 | median_shared | q10 | q90 |
|---|---|---|---|---|---|---|
| CLR | 0.1544573 | 0.7059236 | 2 | 5 | 0 | 21 |
| rCLR | 0.1657904 | 0.6847763 | 2 | 4 | 0 | 21 |

## Sanity checks

| check | value |
|---|---|
| genes md5 equals keff_test fingerprint | TRUE |
| max |CLR recomputed - cached CLR| | 0 |
| cor(CLR, rCLR) over upper-triangle pairs | 0.796 |

