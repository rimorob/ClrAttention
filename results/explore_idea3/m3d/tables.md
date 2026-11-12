### Equivalence: equivalence_criterion.txt

```
held-out criterion machinery vs results/robust_z/m3d/heldout_curves.csv (V0): 128 split-direction-depth rows compared, max |diff| 4.996e-16, tolerance 1e-9 PASS
```

### Equivalence: equivalence_scoring.txt

```
gate 1: t = 9 readout vs results/baselines_m3d/metrics.csv (attention): max |diff| pooled AUPR 2.498e-16 over truths SC,all; tolerance 1e-9 PASS
gate 3: t = 9 enrichment (raw, hub MI, hub pi) vs results/robust_z/m3d V0 att9: max |diff| 3.553e-15 over truths SC,all; tolerance 1e-9 PASS
```

### Knee null (step 0): between-decile statistics of log2 of the decile-median knee, real vs rewired operators (50) and vs label permutations (2,000)

| stratifier | real_sd | null_sd_mean | null_sd_q95 | frac_null_sd_ge_real | real_range | null_range_mean | frac_null_range_ge_real | real_spearman | null_spearman_mean | n_rewirings | real_sd_above_q95 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| mi | 0 | 0 | 0 | 1 | 0 | 0 | 1 |  |  | 50 | False |
| pi | 0.3162 | 0 | 0 | 0 | 1 | 0 | 0 | -0.5222 |  | 50 | True |
| mi_labelperm | 0 | 0 | 0 | 1 | 0 | 0 | 1 |  |  | 2000 | False |
| pi_labelperm | 0.3162 | 0 | 0 | 0 | 1 | 0 | 0 | -0.5222 |  | 2000 | True |


### Per-decile knees (median of per-gene knee, log2; snapped grid depth) and the rewired-operator values

| stratifier | decile | real_median_log2_knee | real_knee_snapped | rewired_mean_log2_knee | rewired_q025 | rewired_q975 |
|---|---|---|---|---|---|---|
| mi | 1 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| mi | 2 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| mi | 3 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| mi | 4 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| mi | 5 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| mi | 6 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| mi | 7 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| mi | 8 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| mi | 9 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| mi | 10 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| pi | 1 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| pi | 2 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| pi | 3 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| pi | 4 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| pi | 5 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| pi | 6 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| pi | 7 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| pi | 8 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| pi | 9 | 4 | 16 | 3.585 | 3.585 | 3.585 |
| pi | 10 | 3 | 8 | 3.585 | 3.585 | 3.585 |


### Selected depth per decile, stratifier mi (op real; mean criterion over 32 split-directions)

| decile | t_d | t_argmax | t_lo | t_hi | width_octaves | crit_max | se_at_max | splitdir_argmax_median | splitdir_argmax_q25 | splitdir_argmax_q75 |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 6 | 8 | 6 | 8 | 0.415 | 0.1709 | 0.002681 | 8 | 6 | 8 |
| 2 | 6 | 6 | 4 | 8 | 1 | 0.1756 | 0.002269 | 6 | 6 | 8 |
| 3 | 4 | 4 | 3 | 6 | 1 | 0.1679 | 0.0017 | 4 | 3 | 6 |
| 4 | 4 | 4 | 3 | 8 | 1.415 | 0.1945 | 0.001622 | 5 | 3 | 6 |
| 5 | 4 | 4 | 3 | 8 | 1.415 | 0.2048 | 0.001496 | 5 | 3 | 6 |
| 6 | 3 | 3 | 2 | 4 | 1 | 0.2018 | 0.00155 | 3 | 3 | 3 |
| 7 | 3 | 3 | 3 | 4 | 0.415 | 0.2093 | 0.001693 | 3 | 3 | 3 |
| 8 | 3 | 3 | 3 | 3 | 0 | 0.1967 | 0.001893 | 3 | 3 | 3 |
| 9 | 2 | 3 | 2 | 3 | 0.585 | 0.2153 | 0.001703 | 3 | 3 | 3 |
| 10 | 3 | 3 | 3 | 3 | 0 | 0.2897 | 0.002608 | 3 | 3 | 3 |


### Selected depth per decile, stratifier pi (op real; mean criterion over 32 split-directions)

| decile | t_d | t_argmax | t_lo | t_hi | width_octaves | crit_max | se_at_max | splitdir_argmax_median | splitdir_argmax_q25 | splitdir_argmax_q75 |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 3 | 3 | 3 | 3 | 0 | 0.04846 | 0.0005394 | 3 | 3 | 3 |
| 2 | 2 | 2 | 2 | 3 | 0.585 | 0.07624 | 0.0007795 | 2 | 2 | 3 |
| 3 | 2 | 2 | 2 | 3 | 0.585 | 0.1028 | 0.001159 | 2 | 2 | 2 |
| 4 | 2 | 2 | 2 | 3 | 0.585 | 0.1286 | 0.00114 | 2.5 | 2 | 3 |
| 5 | 2 | 3 | 2 | 3 | 0.585 | 0.1569 | 0.001687 | 3 | 2 | 3 |
| 6 | 2 | 3 | 2 | 3 | 0.585 | 0.1844 | 0.001761 | 3 | 3 | 3 |
| 7 | 3 | 3 | 3 | 3 | 0 | 0.2251 | 0.001708 | 3 | 3 | 3 |
| 8 | 3 | 3 | 3 | 4 | 0.415 | 0.2816 | 0.002324 | 3 | 3 | 3 |
| 9 | 3 | 3 | 3 | 4 | 0.415 | 0.3689 | 0.003302 | 3 | 3 | 4 |
| 10 | 4 | 4 | 3 | 6 | 1 | 0.5527 | 0.004555 | 3.5 | 3 | 4 |


### Selected depth per decile on the REWIRED operator (stratifier mi)

| decile | t_d | t_argmax | t_lo | t_hi | width_octaves | crit_max | se_at_max |
|---|---|---|---|---|---|---|---|
| 1 | 24 | 32 | 4 | 128 | 5 | 0.005815 | 0.0001205 |
| 2 | 32 | 128 | 8 | 128 | 4 | 0.007397 | 0.000123 |
| 3 | 32 | 64 | 12 | 128 | 3.415 | 0.008177 | 0.0001231 |
| 4 | 48 | 128 | 16 | 128 | 3 | 0.01046 | 0.0001047 |
| 5 | 48 | 128 | 16 | 128 | 3 | 0.01243 | 8.769e-05 |
| 6 | 48 | 128 | 16 | 128 | 3 | 0.01245 | 0.0001051 |
| 7 | 48 | 128 | 16 | 128 | 3 | 0.01468 | 0.0001627 |
| 8 | 48 | 128 | 16 | 128 | 3 | 0.01519 | 0.000223 |
| 9 | 48 | 128 | 24 | 128 | 2.415 | 0.0188 | 0.0002191 |
| 10 | 48 | 128 | 24 | 128 | 2.415 | 0.03866 | 0.0008252 |


### Selection statistics (SD of log2 t_d over deciles, range, trend, deciles >= 1 octave from t = 9, deciles differing from the global-grid depth)

| op | labelset | sd_log2_t | range_octaves | spearman_decile_vs_log2t | n_deciles_ge1oct_from_9 | n_deciles_differing_from_global_grid | mean_width_octaves | global_grid_t |
|---|---|---|---|---|---|---|---|---|
| real | global |  | 0 |  | 0 | 0 | 0.585 | 8 |
| real | mi | 0.4889 | 1.585 | -0.9217 | 8 | 10 | 0.7245 | 8 |
| real | perm1 | 0.2005 | 0.415 | -0.3419 | 10 | 10 | 0.6735 | 8 |
| real | perm2 | 0.1312 | 0.415 | -0.4062 | 10 | 10 | 0.4735 | 8 |
| real | perm3 | 0.2005 | 0.415 | -0.4179 | 10 | 10 | 0.632 | 8 |
| real | perm4 | 0.2143 | 0.415 | -0.2843 | 10 | 10 | 0.649 | 8 |
| real | pi | 0.3732 | 1 | 0.603 | 10 | 10 | 0.4755 | 8 |
| rewired | mi | 0.3672 | 1 | 0.809 | 10 | 10 | 3.225 | 8 |


### Controls: parameters chosen by the same criterion (global grid depth, PPR alpha, heat-kernel tau)

| family | param | param_argmax | param_lo | param_hi | width_octaves | crit_at_choice | crit_max | n_curves |
|---|---|---|---|---|---|---|---|---|
| ppr | 0.9 | 0.95 | 0.9 | 0.95 | 1 | 0.2456 | 0.2478 | 32 |
| hk | 8 | 12 | 8 | 12 | 0.585 | 0.2412 | 0.2413 | 32 |
| global_grid | 8 | 8 | 8 | 12 | 0.585 | 0.239 | 0.239 | 32 |


### Pooled AUPR, truth SC (base rate 0.0943)

| readout | AUPR raw | AUPR hub-strat (MI row mean) | AUPR hub-strat (pi) |
|---|---|---|---|
| t9 | 0.15791 | 0.15855 | 0.14715 |
| global_grid | 0.15656 | 0.15714 | 0.14675 |
| form1_mi | 0.14949 | 0.14893 | 0.14433 |
| form1_pi | 0.14859 | 0.15015 | 0.14502 |
| form2_mi | 0.16427 | 0.16546 | 0.14952 |
| form2_pi | 0.16303 | 0.16425 | 0.14853 |
| form2raw | 0.16249 | 0.16426 | 0.14787 |
| form1_perm1 | 0.14846 | 0.14896 | 0.14430 |
| ppr | 0.16052 | 0.16146 | 0.14882 |
| hk | 0.15712 | 0.15775 | 0.14710 |


### AUPR difference to t = 9, truth SC (paired delete-half gene jackknife, descriptive)

| readout | raw minus t=9 | hub MI minus t=9 | hub pi minus t=9 |
|---|---|---|---|
| global_grid | -0.0014 [-0.0018, -0.0009] | -0.0014 [-0.0019, -0.0010] | -0.0004 [-0.0006, -0.0002] |
| form1_mi | -0.0084 [-0.0113, -0.0055] | -0.0096 [-0.0122, -0.0070] | -0.0028 [-0.0049, -0.0008] |
| form1_pi | -0.0093 [-0.0123, -0.0064] | -0.0084 [-0.0114, -0.0054] | -0.0021 [-0.0046, +0.0004] |
| form2_mi | +0.0064 [+0.0037, +0.0090] | +0.0069 [+0.0042, +0.0096] | +0.0024 [+0.0010, +0.0038] |
| form2_pi | +0.0051 [+0.0030, +0.0072] | +0.0057 [+0.0035, +0.0079] | +0.0014 [+0.0001, +0.0027] |
| form2raw | +0.0046 [+0.0006, +0.0085] | +0.0057 [+0.0016, +0.0098] | +0.0007 [-0.0015, +0.0029] |
| form1_perm1 | -0.0095 [-0.0123, -0.0066] | -0.0096 [-0.0125, -0.0067] | -0.0029 [-0.0046, -0.0011] |
| ppr | +0.0026 [+0.0011, +0.0041] | +0.0029 [+0.0014, +0.0044] | +0.0017 [+0.0009, +0.0024] |
| hk | -0.0008 [-0.0012, -0.0004] | -0.0008 [-0.0012, -0.0004] | -0.0001 [-0.0003, +0.0002] |


### Pooled AUPR, truth all (base rate 0.1205)

| readout | AUPR raw | AUPR hub-strat (MI row mean) | AUPR hub-strat (pi) |
|---|---|---|---|
| t9 | 0.17334 | 0.17359 | 0.16275 |
| global_grid | 0.17224 | 0.17249 | 0.16277 |
| form1_mi | 0.16536 | 0.16537 | 0.16154 |
| form1_pi | 0.16724 | 0.16903 | 0.16371 |
| form2_mi | 0.17884 | 0.17926 | 0.16277 |
| form2_pi | 0.17702 | 0.17745 | 0.16168 |
| form2raw | 0.17622 | 0.17702 | 0.15955 |
| form1_perm1 | 0.16662 | 0.16705 | 0.16334 |
| ppr | 0.17591 | 0.17638 | 0.16266 |
| hk | 0.17262 | 0.17291 | 0.16270 |


### AUPR difference to t = 9, truth all (paired delete-half gene jackknife, descriptive)

| readout | raw minus t=9 | hub MI minus t=9 | hub pi minus t=9 |
|---|---|---|---|
| global_grid | -0.0011 [-0.0014, -0.0008] | -0.0011 [-0.0014, -0.0008] | +0.0000 [-0.0001, +0.0002] |
| form1_mi | -0.0080 [-0.0103, -0.0057] | -0.0082 [-0.0100, -0.0064] | -0.0012 [-0.0032, +0.0008] |
| form1_pi | -0.0061 [-0.0082, -0.0040] | -0.0046 [-0.0067, -0.0024] | +0.0010 [-0.0012, +0.0031] |
| form2_mi | +0.0055 [+0.0040, +0.0070] | +0.0057 [+0.0041, +0.0072] | +0.0000 [-0.0009, +0.0009] |
| form2_pi | +0.0037 [+0.0024, +0.0049] | +0.0039 [+0.0026, +0.0051] | -0.0011 [-0.0019, -0.0002] |
| form2raw | +0.0029 [+0.0006, +0.0052] | +0.0034 [+0.0010, +0.0058] | -0.0032 [-0.0046, -0.0018] |
| form1_perm1 | -0.0067 [-0.0084, -0.0050] | -0.0065 [-0.0083, -0.0048] | +0.0006 [-0.0007, +0.0019] |
| ppr | +0.0026 [+0.0015, +0.0036] | +0.0028 [+0.0017, +0.0039] | -0.0001 [-0.0006, +0.0004] |
| hk | -0.0007 [-0.0009, -0.0005] | -0.0007 [-0.0009, -0.0005] | -0.0001 [-0.0002, +0.0001] |


### Literal pass rule (hub-stratified diff >= +0.005 and raw diff >= -0.005, in AUPR)

| readout | truth | hub_mi_diff | hub_pi_diff | raw_diff | pass_primary_mi | pass_sensitivity_pi |
|---|---|---|---|---|---|---|
| global_grid | SC | -0.001414 | -0.0004075 | -0.001352 | False | False |
| form1_mi | SC | -0.009622 | -0.002824 | -0.008418 | False | False |
| form1_pi | SC | -0.008404 | -0.002134 | -0.009321 | False | False |
| form2_mi | SC | 0.006913 | 0.002369 | 0.006362 | True | False |
| form2_pi | SC | 0.005695 | 0.001377 | 0.005122 | True | False |
| form2raw | SC | 0.005708 | 0.0007193 | 0.004585 | True | False |
| form1_perm1 | SC | -0.009592 | -0.002851 | -0.00945 | False | False |
| ppr | SC | 0.00291 | 0.001664 | 0.002616 | False | False |
| hk | SC | -0.0007974 | -5.05e-05 | -0.0007902 | False | False |

