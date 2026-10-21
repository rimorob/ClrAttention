# Method comparison

*Generated 2026-10-21 06:49 by `analysis/sota_table.R`.* Pooled co-membership AUPR divided by the base rate (1.00 = random), CLR attention at t = 9. **Bold**: best in the row or tied with it (95% interval of the paired jackknife difference includes 0). "-": not run.

## Raw scores

| Dataset | Truth | Attention | CLR | GENIE3 | PIDC | GRNBoost2 | \|r\| | Hub only |
|---|---|---|---|---|---|---|---|---|
| M3D | SC | **1.68** | 1.44 | 1.51 | 1.45 | 1.18 | 1.23 | 1.08 |
| M3D | all | **1.44** | 1.37 | 1.40 | 1.40 | 1.13 | 1.23 | 1.22 |
| EcoMAC | SC | **1.68** | 1.39 | 1.50 | 1.43 | 1.21 | 1.18 | 1.04 |
| EcoMAC | all | **1.47** | 1.33 | 1.38 | 1.39 | 1.14 | 1.19 | 1.20 |
| K562 | DoRothEA AB | **1.09** | 1.06 | **1.08** | **1.07** | 1.02 | 1.05 | 1.00 |
| K562 | DoRothEA ABC | **1.07** | 1.03 | 1.04 | 1.04 | 1.01 | 1.03 | 1.02 |
| K562 | CollecTRI | **1.10** | 1.03 | 1.04 | 1.05 | 1.01 | 1.01 | **1.14** |
| K562 | CORUM | 3.51 | **6.39** | 6.32 | **6.33** | 2.18 | **7.00** | 1.55 |
| K562 | KD response | **2.38** | **2.42** | **2.45** | 2.35 | 1.37 | **2.45** | **2.42** |

## Hub-controlled (within-stratum percentile in 10 strata of pi_i + pi_j)

| Dataset | Truth | Attention | CLR | GENIE3 | PIDC |
|---|---|---|---|---|---|
| M3D | SC | **1.56** | 1.40 | 1.47 | 1.40 |
| M3D | all | **1.35** | 1.33 | **1.37** | 1.34 |
| EcoMAC | SC | **1.57** | 1.37 | 1.47 | 1.40 |
| EcoMAC | all | **1.36** | 1.29 | **1.35** | **1.33** |
| K562 | DoRothEA AB | **1.09** | 1.06 | **1.07** | 1.07 |
| K562 | DoRothEA ABC | **1.06** | 1.02 | 1.04 | 1.03 |
| K562 | CollecTRI | **1.02** | 1.01 | **1.03** | **1.03** |
| K562 | CORUM | 2.08 | **3.40** | **3.47** | **3.82** |
| K562 | KD response | 1.60 | 1.91 | 1.93 | **2.00** |

"?": tie test unavailable (no jackknife replicates saved for this run yet).
