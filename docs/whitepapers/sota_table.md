# Method comparison

*Generated 2026-10-21 10:48 by `analysis/sota_table.R`.* Pooled co-membership AUPR divided by the base rate (1.00 = random); CLR attention at t = 9. **Bold**: best in the row or tied with it (Bonferroni-corrected 95% interval of the paired jackknife difference includes 0); **x**† tied only after the correction. "-": not run.

## Raw scores

| Dataset | Truth | Base rate | Attention | CLR | GENIE3 | PIDC | GRNBoost2 | \|r\| | Hub only |
|---|---|---|---|---|---|---|---|---|---|
| M3D | SC | 0.094 | **1.68** | 1.44 | 1.51 | 1.45 | 1.18 | 1.23 | 1.08 |
| M3D | all | 0.120 | **1.44** | 1.37 | **1.40**† | **1.40**† | 1.13 | 1.23 | 1.22 |
| EcoMAC | SC | 0.094 | **1.68** | 1.39 | 1.50 | 1.43 | 1.21 | 1.18 | 1.04 |
| EcoMAC | all | 0.121 | **1.47** | 1.33 | 1.38 | 1.39 | 1.14 | 1.19 | 1.20 |
| K562 | DoRothEA AB | 0.183 | **1.09** | **1.06**† | **1.08** | **1.07** | 1.02 | **1.05**† | 1.00 |
| K562 | DoRothEA ABC | 0.130 | **1.07** | 1.03 | **1.04**† | **1.04**† | 1.01 | 1.03 | 1.02 |
| K562 | CollecTRI | 0.128 | **1.10** | 1.03 | 1.04 | 1.05 | 1.01 | 1.01 | **1.14** |
| K562 | CORUM | 0.013 | 3.51 | **6.39** | 6.32 | **6.33** | 2.18 | **7.00** | 1.55 |
| K562 | KD response | 0.079 | **2.38** | **2.42** | **2.45** | **2.35**† | 1.37 | **2.45** | **2.42** |

## Hub-controlled (within-stratum percentile in 10 strata of pi_i + pi_j)

| Dataset | Truth | Base rate | Attention | CLR | GENIE3 | PIDC | GRNBoost2 | \|r\| |
|---|---|---|---|---|---|---|---|---|
| M3D | SC | 0.094 | **1.56** | 1.40 | 1.47 | 1.40 | 1.28 | 1.22 |
| M3D | all | 0.120 | **1.35** | 1.33 | **1.37** | 1.34 | 1.18 | 1.21 |
| EcoMAC | SC | 0.094 | **1.57** | 1.37 | 1.47 | 1.40 | 1.19 | 1.16 |
| EcoMAC | all | 0.121 | **1.36** | 1.29 | **1.35** | **1.33** | 1.09 | 1.15 |
| K562 | DoRothEA AB | 0.183 | **1.09** | 1.06 | **1.07** | **1.07**† | 1.02 | **1.05**† |
| K562 | DoRothEA ABC | 0.130 | **1.06** | 1.02 | **1.04**† | **1.03**† | 1.00 | **1.02**† |
| K562 | CollecTRI | 0.128 | **1.02** | 1.01 | **1.03** | **1.03** | 0.94 | 0.99 |
| K562 | CORUM | 0.013 | 2.08 | **3.40** | **3.47** | **3.82** | 1.70 | **3.76** |
| K562 | KD response | 0.079 | 1.60 | 1.91 | 1.93 | **2.00** | 1.08 | 1.88 |

**Notes.**
- Intervals resample genes (delete-half jackknife, 100 replicates) with the fitted networks held fixed; they do not resample experiments or cells, or vary the random seeds of GENIE3/GRNBoost2.
- Depth t = 9 was chosen on M3D with a label-free held-out criterion (lab splits; no truths used), then transferred unchanged to EcoMAC and K562. M3D is therefore not an untouched test set for the depth.
- Hub control: each method's scores are converted to percentiles within 10 strata of pi_i + pi_j, where pi is the stationary distribution of the attention walk (an attention-derived stratifier; a method-neutral one is planned). "Hub only" ranks pairs by pi_i + pi_j alone (a degree baseline).
- K562: the Perturb-seq essential compendium, network built from knockdown half A on the 5,671 genes in at least one literature truth; the knockdown-response truth comes from the other half (B). All K562 TF-truth rows are within 10% of random.
- EcoMAC may share arrays with M3D, so it is not yet an independent replication.
- GRNBoost2 returns a sparse network (most pairs tied at zero), which lowers its pooled AUPR.
