### Equivalence: equivalence_criterion.txt

```
held-out criterion machinery vs results/robust_z/k562/heldout_curves.csv (V0): 96 split-direction-depth rows compared, max |diff| 4.996e-16, tolerance 1e-9 PASS
```

### Equivalence: equivalence_scoring.txt

```
gate 1: t = 9 readout vs results/baselines_k562ess_truthgenes_genie3/metrics.csv (attention): max |diff| pooled AUPR 3.886e-16 over truths DoRothEA_AB,CollecTRI,CORUM_all; tolerance 1e-9 PASS
gate 3: t = 9 enrichment (raw, hub MI, hub pi) vs results/robust_z/k562 V0 att9: max |diff| 4.219e-15 over truths DoRothEA_AB,CORUM_all,CORUM_noRP; tolerance 1e-9 PASS
```

### Knee null (step 0): between-decile statistics of log2 of the decile-median knee, real vs rewired operators (50) and vs label permutations (2,000)

| stratifier | real_sd | null_sd_mean | null_sd_q95 | frac_null_sd_ge_real | real_range | null_range_mean | frac_null_range_ge_real | real_spearman | null_spearman_mean | n_rewirings | real_sd_above_q95 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| mi | 0.175 | 0 | 0 | 0 | 0.415 | 0 | 0 | 0.6963 |  | 50 | True |
| pi | 0.2005 | 0 | 0 | 0 | 0.415 | 0 | 0 | 0.5698 |  | 50 | True |
| mi_labelperm | 0.175 | 0 | 0 | 0 | 0.415 | 0 | 0 | 0.6963 |  | 2000 | True |
| pi_labelperm | 0.2005 | 0 | 0 | 0 | 0.415 | 0 | 0 | 0.5698 |  | 2000 | True |


### Per-decile knees (median of per-gene knee, log2; snapped grid depth) and the rewired-operator values

| stratifier | decile | real_median_log2_knee | real_knee_snapped | rewired_mean_log2_knee | rewired_q025 | rewired_q975 |
|---|---|---|---|---|---|---|
| mi | 1 | 3.585 | 12 | 3.585 | 3.585 | 3.585 |
| mi | 2 | 3.585 | 12 | 3.585 | 3.585 | 3.585 |
| mi | 3 | 3.585 | 12 | 3.585 | 3.585 | 3.585 |
| mi | 4 | 3.585 | 12 | 3.585 | 3.585 | 3.585 |
| mi | 5 | 3.585 | 12 | 3.585 | 3.585 | 3.585 |
| mi | 6 | 3.585 | 12 | 3.585 | 3.585 | 3.585 |
| mi | 7 | 3.585 | 12 | 3.585 | 3.585 | 3.585 |
| mi | 8 | 3.585 | 12 | 3.585 | 3.585 | 3.585 |
| mi | 9 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| mi | 10 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| pi | 1 | 3.585 | 12 | 3.585 | 3.585 | 3.585 |
| pi | 2 | 3.585 | 12 | 3.585 | 3.585 | 3.585 |
| pi | 3 | 3.585 | 12 | 3.585 | 3.585 | 3.585 |
| pi | 4 | 3.585 | 12 | 3.585 | 3.585 | 3.585 |
| pi | 5 | 3.585 | 12 | 3.585 | 3.585 | 3.585 |
| pi | 6 | 3.585 | 12 | 3.585 | 3.585 | 3.585 |
| pi | 7 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| pi | 8 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| pi | 9 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| pi | 10 | 3.585 | 12 | 3.585 | 3.585 | 3.585 |


### Selected depth per decile, stratifier mi (op real; mean criterion over 12 split-directions)

| decile | t_d | t_argmax | t_lo | t_hi | width_octaves | crit_max | se_at_max | splitdir_argmax_median | splitdir_argmax_q25 | splitdir_argmax_q75 |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 12 | 12 | 12 | 12 | 0 | 0.07994 | 0.005055 | 12 | 8 | 12 |
| 2 | 12 | 12 | 12 | 12 | 0 | 0.05613 | 0.002719 | 12 | 12 | 12 |
| 3 | 8 | 12 | 8 | 12 | 0.585 | 0.05626 | 0.002897 | 10 | 8 | 12 |
| 4 | 8 | 8 | 8 | 8 | 0 | 0.0445 | 0.001505 | 8 | 8 | 8 |
| 5 | 8 | 8 | 8 | 8 | 0 | 0.06064 | 0.001823 | 8 | 8 | 8 |
| 6 | 6 | 6 | 6 | 8 | 0.415 | 0.08246 | 0.00219 | 6 | 6 | 6.5 |
| 7 | 3 | 3 | 3 | 4 | 0.415 | 0.1087 | 0.002343 | 3 | 3 | 3 |
| 8 | 3 | 3 | 3 | 3 | 0 | 0.1552 | 0.002902 | 3 | 3 | 3 |
| 9 | 2 | 3 | 2 | 3 | 0.585 | 0.2468 | 0.003043 | 3 | 3 | 3 |
| 10 | 3 | 3 | 3 | 3 | 0 | 0.3986 | 0.003076 | 3 | 3 | 3 |


### Selected depth per decile, stratifier pi (op real; mean criterion over 12 split-directions)

| decile | t_d | t_argmax | t_lo | t_hi | width_octaves | crit_max | se_at_max | splitdir_argmax_median | splitdir_argmax_q25 | splitdir_argmax_q75 |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 6 | 6 | 6 | 8 | 0.415 | 0.04314 | 0.001244 | 6 | 6 | 6 |
| 2 | 4 | 6 | 4 | 6 | 0.585 | 0.03473 | 0.0009983 | 6 | 4 | 6 |
| 3 | 4 | 6 | 4 | 6 | 0.585 | 0.03745 | 0.0009394 | 6 | 4 | 6.5 |
| 4 | 4 | 4 | 3 | 6 | 1 | 0.03618 | 0.0008343 | 5 | 3 | 6 |
| 5 | 4 | 4 | 3 | 6 | 1 | 0.04195 | 0.0008751 | 4 | 3.75 | 4 |
| 6 | 3 | 3 | 3 | 4 | 0.415 | 0.05252 | 0.001247 | 3.5 | 3 | 4 |
| 7 | 3 | 3 | 3 | 4 | 0.415 | 0.1058 | 0.002202 | 3 | 3 | 3 |
| 8 | 3 | 3 | 3 | 3 | 0 | 0.2038 | 0.003703 | 3 | 3 | 3 |
| 9 | 3 | 3 | 3 | 3 | 0 | 0.2909 | 0.005132 | 3 | 3 | 3 |
| 10 | 3 | 3 | 3 | 3 | 0 | 0.4135 | 0.006494 | 3 | 3 | 3 |


### Selected depth per decile on the REWIRED operator (stratifier mi)

| decile | t_d | t_argmax | t_lo | t_hi | width_octaves | crit_max | se_at_max |
|---|---|---|---|---|---|---|---|
| 1 | 24 | 128 | 6 | 128 | 4.415 | 0.00504 | 0.0001386 |
| 2 | 32 | 32 | 8 | 128 | 4 | 0.005042 | 0.0001221 |
| 3 | 32 | 128 | 12 | 128 | 3.415 | 0.005658 | 6.086e-05 |
| 4 | 48 | 96 | 16 | 128 | 3 | 0.005814 | 0.0001101 |
| 5 | 48 | 24 | 16 | 128 | 3 | 0.007744 | 8.855e-05 |
| 6 | 48 | 64 | 16 | 128 | 3 | 0.01009 | 0.0001199 |
| 7 | 48 | 24 | 16 | 128 | 3 | 0.01452 | 0.0001447 |
| 8 | 48 | 64 | 24 | 128 | 2.415 | 0.02369 | 0.0002862 |
| 9 | 48 | 96 | 24 | 128 | 2.415 | 0.04378 | 0.0006932 |
| 10 | 48 | 128 | 24 | 128 | 2.415 | 0.1085 | 0.001682 |


### Selection statistics (SD of log2 t_d over deciles, range, trend, deciles >= 1 octave from t = 9, deciles differing from the global-grid depth)

| op | labelset | sd_log2_t | range_octaves | spearman_decile_vs_log2t | n_deciles_ge1oct_from_9 | n_deciles_differing_from_global_grid | mean_width_octaves | global_grid_t |
|---|---|---|---|---|---|---|---|---|
| real | global |  | 0 |  | 0 | 0 | 0 | 8 |
| real | mi | 0.9339 | 2.585 | -0.9474 | 4 | 7 | 0.2 | 8 |
| real | perm1 | 0.2826 | 0.585 | 0.2659 | 7 | 10 | 0.8585 | 8 |
| real | perm2 | 0.3021 | 0.585 | -0.5685 | 6 | 10 | 0.8415 | 8 |
| real | perm3 | 0.2826 | 0.585 | 0.4179 | 7 | 10 | 0.8415 | 8 |
| real | perm4 | 0.2826 | 0.585 | 0.3419 | 7 | 10 | 0.8585 | 8 |
| real | pi | 0.3302 | 1 | -0.9045 | 9 | 10 | 0.4415 | 8 |
| rewired | mi | 0.3672 | 1 | 0.809 | 10 | 10 | 3.108 | 8 |


### Controls: parameters chosen by the same criterion (global grid depth, PPR alpha, heat-kernel tau)

| family | param | param_argmax | param_lo | param_hi | width_octaves | crit_at_choice | crit_max | n_curves |
|---|---|---|---|---|---|---|---|---|
| ppr | 0.85 | 0.85 | 0.8 | 0.9 | 1 | 0.1527 | 0.1527 | 12 |
| hk | 8 | 8 | 8 | 8 | 0 | 0.1592 | 0.1592 | 12 |
| global_grid | 8 | 8 | 8 | 8 | 0 | 0.1606 | 0.1606 | 12 |


### Pooled AUPR, truth DoRothEA_AB (base rate 0.1830)

| readout | AUPR raw | AUPR hub-strat (MI row mean) | AUPR hub-strat (pi) |
|---|---|---|---|
| t9 | 0.20032 | 0.20346 | 0.19952 |
| global_grid | 0.20060 | 0.20360 | 0.19949 |
| form1_mi | 0.20093 | 0.20189 | 0.19955 |
| form1_pi | 0.19891 | 0.20098 | 0.19724 |
| form2_mi | 0.19903 | 0.20267 | 0.19961 |
| form2_pi | 0.19886 | 0.20249 | 0.19954 |
| form2raw | 0.19732 | 0.20090 | 0.19826 |
| ppr | 0.19949 | 0.20274 | 0.19960 |
| hk | 0.20030 | 0.20340 | 0.19957 |


### AUPR difference to t = 9, truth DoRothEA_AB (paired delete-half gene jackknife, descriptive)

| readout | raw minus t=9 | hub MI minus t=9 | hub pi minus t=9 |
|---|---|---|---|
| global_grid | +0.0003 [-0.0001, +0.0006] | +0.0001 [-0.0002, +0.0005] | -0.0000 [-0.0002, +0.0001] |
| form1_mi | +0.0006 [-0.0014, +0.0027] | -0.0016 [-0.0033, +0.0002] | +0.0000 [-0.0014, +0.0015] |
| form1_pi | -0.0014 [-0.0040, +0.0012] | -0.0025 [-0.0053, +0.0003] | -0.0023 [-0.0038, -0.0008] |
| form2_mi | -0.0013 [-0.0023, -0.0003] | -0.0008 [-0.0019, +0.0003] | +0.0001 [-0.0005, +0.0006] |
| form2_pi | -0.0015 [-0.0026, -0.0004] | -0.0010 [-0.0022, +0.0003] | +0.0000 [-0.0006, +0.0006] |
| form2raw | -0.0030 [-0.0042, -0.0018] | -0.0026 [-0.0039, -0.0012] | -0.0013 [-0.0021, -0.0004] |
| ppr | -0.0008 [-0.0013, -0.0004] | -0.0007 [-0.0012, -0.0002] | +0.0001 [-0.0003, +0.0004] |
| hk | -0.0000 [-0.0002, +0.0001] | -0.0001 [-0.0003, +0.0001] | +0.0000 [-0.0001, +0.0002] |


### Pooled AUPR, truth CORUM_all (base rate 0.0133)

| readout | AUPR raw | AUPR hub-strat (MI row mean) | AUPR hub-strat (pi) |
|---|---|---|---|
| t9 | 0.04654 | 0.03934 | 0.02752 |
| global_grid | 0.04531 | 0.03846 | 0.02744 |
| form1_mi | 0.03746 | 0.03332 | 0.02530 |
| form1_pi | 0.03698 | 0.03255 | 0.02651 |
| form2_mi | 0.04994 | 0.04250 | 0.02808 |
| form2_pi | 0.04929 | 0.04133 | 0.02749 |
| form2raw | 0.04684 | 0.03974 | 0.02694 |
| ppr | 0.04180 | 0.03584 | 0.02788 |
| hk | 0.04479 | 0.03807 | 0.02761 |


### AUPR difference to t = 9, truth CORUM_all (paired delete-half gene jackknife, descriptive)

| readout | raw minus t=9 | hub MI minus t=9 | hub pi minus t=9 |
|---|---|---|---|
| global_grid | -0.0012 [-0.0022, -0.0002] | -0.0009 [-0.0017, -0.0001] | -0.0001 [-0.0004, +0.0002] |
| form1_mi | -0.0091 [-0.0150, -0.0031] | -0.0060 [-0.0102, -0.0018] | -0.0022 [-0.0040, -0.0004] |
| form1_pi | -0.0096 [-0.0159, -0.0032] | -0.0068 [-0.0113, -0.0023] | -0.0010 [-0.0026, +0.0006] |
| form2_mi | +0.0034 [+0.0005, +0.0063] | +0.0032 [-0.0002, +0.0065] | +0.0006 [-0.0006, +0.0017] |
| form2_pi | +0.0028 [+0.0002, +0.0053] | +0.0020 [-0.0005, +0.0045] | -0.0000 [-0.0009, +0.0009] |
| form2raw | +0.0003 [-0.0025, +0.0031] | +0.0004 [-0.0018, +0.0025] | -0.0006 [-0.0018, +0.0006] |
| ppr | -0.0047 [-0.0089, -0.0006] | -0.0035 [-0.0066, -0.0004] | +0.0004 [-0.0006, +0.0013] |
| hk | -0.0017 [-0.0033, -0.0002] | -0.0013 [-0.0025, -0.0001] | +0.0001 [-0.0003, +0.0005] |


### Pooled AUPR, truth CORUM_noRP (base rate 0.0119)

| readout | AUPR raw | AUPR hub-strat (MI row mean) | AUPR hub-strat (pi) |
|---|---|---|---|
| t9 | 0.02432 | 0.02583 | 0.01871 |
| global_grid | 0.02419 | 0.02560 | 0.01877 |
| form1_mi | 0.02240 | 0.02430 | 0.01820 |
| form1_pi | 0.02229 | 0.02312 | 0.01872 |
| form2_mi | 0.02425 | 0.02572 | 0.01855 |
| form2_pi | 0.02392 | 0.02553 | 0.01825 |
| form2raw | 0.02311 | 0.02484 | 0.01797 |
| ppr | 0.02314 | 0.02422 | 0.01902 |
| hk | 0.02406 | 0.02537 | 0.01886 |


### AUPR difference to t = 9, truth CORUM_noRP (paired delete-half gene jackknife, descriptive)

| readout | raw minus t=9 | hub MI minus t=9 | hub pi minus t=9 |
|---|---|---|---|
| global_grid | -0.0001 [-0.0004, +0.0002] | -0.0002 [-0.0006, +0.0001] | +0.0001 [-0.0001, +0.0002] |
| form1_mi | -0.0019 [-0.0038, -0.0000] | -0.0015 [-0.0028, -0.0003] | -0.0005 [-0.0013, +0.0002] |
| form1_pi | -0.0020 [-0.0040, -0.0000] | -0.0027 [-0.0049, -0.0005] | +0.0000 [-0.0007, +0.0007] |
| form2_mi | -0.0001 [-0.0010, +0.0009] | -0.0001 [-0.0008, +0.0006] | -0.0002 [-0.0005, +0.0002] |
| form2_pi | -0.0004 [-0.0013, +0.0005] | -0.0003 [-0.0011, +0.0005] | -0.0005 [-0.0008, -0.0001] |
| form2raw | -0.0012 [-0.0030, +0.0006] | -0.0010 [-0.0027, +0.0007] | -0.0007 [-0.0017, +0.0002] |
| ppr | -0.0012 [-0.0027, +0.0003] | -0.0016 [-0.0032, +0.0000] | +0.0003 [-0.0002, +0.0008] |
| hk | -0.0003 [-0.0008, +0.0003] | -0.0005 [-0.0011, +0.0001] | +0.0001 [-0.0001, +0.0004] |


### Pooled AUPR, truth CollecTRI (base rate 0.1278)

| readout | AUPR raw | AUPR hub-strat (MI row mean) | AUPR hub-strat (pi) |
|---|---|---|---|
| t9 | 0.14021 | 0.14287 | 0.12990 |
| global_grid | 0.13963 | 0.14205 | 0.13006 |
| form1_mi | 0.13924 | 0.14019 | 0.13134 |
| form1_pi | 0.13382 | 0.13508 | 0.12892 |
| form2_mi | 0.14185 | 0.14532 | 0.12934 |
| form2_pi | 0.14210 | 0.14569 | 0.12938 |
| form2raw | 0.14149 | 0.14523 | 0.12904 |
| ppr | 0.14075 | 0.14362 | 0.13000 |
| hk | 0.14008 | 0.14266 | 0.13003 |


### AUPR difference to t = 9, truth CollecTRI (paired delete-half gene jackknife, descriptive)

| readout | raw minus t=9 | hub MI minus t=9 | hub pi minus t=9 |
|---|---|---|---|
| global_grid | -0.0006 [-0.0009, -0.0003] | -0.0008 [-0.0012, -0.0004] | +0.0002 [+0.0000, +0.0003] |
| form1_mi | -0.0010 [-0.0030, +0.0011] | -0.0027 [-0.0045, -0.0009] | +0.0014 [-0.0001, +0.0030] |
| form1_pi | -0.0064 [-0.0090, -0.0038] | -0.0078 [-0.0106, -0.0049] | -0.0010 [-0.0022, +0.0002] |
| form2_mi | +0.0016 [+0.0006, +0.0026] | +0.0025 [+0.0013, +0.0036] | -0.0006 [-0.0010, -0.0001] |
| form2_pi | +0.0019 [+0.0008, +0.0030] | +0.0028 [+0.0015, +0.0042] | -0.0005 [-0.0010, -0.0000] |
| form2raw | +0.0013 [-0.0001, +0.0027] | +0.0024 [+0.0007, +0.0040] | -0.0009 [-0.0015, -0.0002] |
| ppr | +0.0005 [+0.0001, +0.0009] | +0.0007 [+0.0003, +0.0012] | +0.0001 [-0.0002, +0.0004] |
| hk | -0.0001 [-0.0003, +0.0000] | -0.0002 [-0.0004, -0.0000] | +0.0001 [+0.0000, +0.0002] |


### Literal pass rule (hub-stratified diff >= +0.005 and raw diff >= -0.005, in AUPR)

| readout | truth | hub_mi_diff | hub_pi_diff | raw_diff | pass_primary_mi | pass_sensitivity_pi |
|---|---|---|---|---|---|---|
| global_grid | DoRothEA_AB | 0.0001405 | -3.481e-05 | 0.0002741 | False | False |
| global_grid | CollecTRI | -0.000814 | 0.0001584 | -0.0005869 | False | False |
| form1_mi | DoRothEA_AB | -0.001563 | 2.519e-05 | 0.0006035 | False | False |
| form1_mi | CollecTRI | -0.002679 | 0.001443 | -0.0009767 | False | False |
| form1_pi | DoRothEA_AB | -0.002471 | -0.002284 | -0.001413 | False | False |
| form1_pi | CollecTRI | -0.007791 | -0.000985 | -0.006395 | False | False |
| form2_mi | DoRothEA_AB | -0.0007829 | 8.137e-05 | -0.001295 | False | False |
| form2_mi | CollecTRI | 0.002451 | -0.0005575 | 0.001637 | False | False |
| form2_pi | DoRothEA_AB | -0.0009688 | 1.068e-05 | -0.001462 | False | False |
| form2_pi | CollecTRI | 0.002818 | -0.0005166 | 0.001885 | False | False |
| form2raw | DoRothEA_AB | -0.002553 | -0.001269 | -0.003008 | False | False |
| form2raw | CollecTRI | 0.002357 | -0.0008568 | 0.001272 | False | False |
| ppr | DoRothEA_AB | -0.0007132 | 8.013e-05 | -0.0008301 | False | False |
| ppr | CollecTRI | 0.0007469 | 9.609e-05 | 0.0005386 | False | False |
| hk | DoRothEA_AB | -5.968e-05 | 4.06e-05 | -2.646e-05 | False | False |
| hk | CollecTRI | -0.0002109 | 0.0001295 | -0.0001358 | False | False |

