# Exploratory run of ideation idea 5: self-supervised learned attention on M3D (results)

**Status.** EXPLORATORY. M3D is the project's tuning set and its RegulonDB truths are consumed. Not preregistered, no p-values, no claim for the paper or the pre-read. Plan: `docs/plans/explore_idea5_plan.md` (commit 63aa6ec, fixed before any training or scoring; not edited since). Code: `analysis/idea5/` (`idea5_prep.R`, `idea5_train.py`, `idea5_readout_score.R`, `idea5_summarize.R`, `idea5_ecomac_recon.py`, `idea5_lib.R`). Tables: `results/idea5/{scores,summary}/` and `results/idea5/ecomac_recon.csv`. Intervals are the plan's: difference plus/minus 1.96 SD of the 100 paired delete-half gene-jackknife differences (genes resampled, networks and seeds held fixed).

## 1. Method as run

Matrix: M3D, 4,297 genes x 466 replicate-averaged experiments (the baselines' matrix; not the 907-chip file, plan section 0). One balanced lab split of the 39 labs (`set.seed(20261112)`): 254 and 212 experiments, both directions. One attention layer over genes: logits `L_ij = a_i Z_ij/tau0 + D_ij` (`D = QK'/sqrt(d)`, row-centred; `Z` the CLR score of the training samples; tau0 = 0.683 full, 0.598 and 0.598 for the halves; scratch models have no `Z` term), full softmax over `j != i`, signed-expression values, leave-one-gene-out MSE on standardised expression. Adam, batch 64, 10-step warm-up, epoch cap 600, patience 30 (relative 1e-4), learning rate chosen on tuning seed 0 from {1e-3, 3e-3, 1e-2} by validation MSE (**1e-2 for all four configs, the edge of the grid**). Seeds 1-5, d in {32, 128}, CLR-initialised and scratch: 40 half-trained models, 20 full-data refits (epochs = rounded mean of the two directions' best epochs), 10 untrained scratch controls. Readout exactly as the plan: score matrix `C = pmax(tau0 (L + L')/2, 0)`, `clr:::.softmax_rows(keff = 10, cap = 50)`, lazy walk, `S_t = (P^t + P^t')/2`; label-free depth t\* by the lab-split criterion (AUPR against the held-out half's top-1% CLR pairs, t = 1..32, geometric mean of the two directions). Truths: SC (155 regulons, 1,762,067 pairs, base rate 0.0943) and all evidence (196 regulons, 3,031,397 pairs, 0.1205). Stratifiers (the same for every method): pi of the full-data CLR attention walk (primary) and the MI row mean.

## 2. Gates (all passed, before any truth was scored for the learned models)

| Gate | Result |
|---|---|
| G1 (CLR-initialised step 0 through the learned-logit readout vs the package, full data and each half) | exported D and gamma exactly 0; max abs difference of C 3.6e-15, of the operator 1.7e-15, of the P^9 walk mass 1.2e-16 (limit 1e-8) |
| G3 (ported temperatures vs `clr:::.softmax_rows`) | max difference 0 |
| G2 (scoring code on full-data CLR vs `results/robust_z/m3d` V0, pairwise and t = 9, raw / pi / MI row mean, SC and all) | maximum difference 0 in the 12 AUPRs and in all 1,200 jackknife values (limit 1e-10) |

## 3. Literal pass/fail (plan section 7)

Primary cell: CLR-initialised, d = 32, refit, read out at its own t\*, hub-stratified (pi), SC pooled AUPR, learned minus full-data CLR attention at t = 9 (AUPR 0.1472).

| Seed | learned AUPR | contrast | 95% jackknife interval | class |
|---|---|---|---|---|
| 1 | 0.1292 | -0.0180 | [-0.0232, -0.0128] | loss |
| 2 | 0.1293 | -0.0178 | [-0.0232, -0.0125] | loss |
| 3 | 0.1289 | -0.0183 | [-0.0236, -0.0129] | loss |
| 4 | 0.1281 | -0.0191 | [-0.0246, -0.0135] | loss |
| 5 | 0.1288 | -0.0183 | [-0.0237, -0.0130] | loss |

**Verdict by the rule: learned LOSES** (0 wins, 5 losses, 0 ties, 0 unresolved gains of 5 seeds). Mean contrast -0.0183, SD across seeds 0.0005 (the SD of the jackknife replicates is 0.0027, so the seed-to-seed spread is small against the gene-resampling spread). The MI-row-mean stratifier gives the same verdict (-0.0294, 5 of 5 losses), so the two stratifiers do not disagree.

Descriptive, same rule applied to other cells (SC; mean contrast over seeds, SD across seeds; 5 of 5 seeds are losses in every row):

| Config | read-out | hub (pi) | raw | MI row mean (t\*) |
|---|---|---|---|---|
| CLR-init d = 32 | t\* | -0.0183 (0.0005) | -0.0358 (0.0004) | -0.0294 (0.0004) |
| CLR-init d = 32 | t = 9 | -0.0206 (0.0014) | -0.0304 (0.0015) | |
| CLR-init d = 128 | t\* | -0.0123 (0.0004) | -0.0256 (0.0004) | -0.0235 (0.0005) |
| CLR-init d = 128 | t = 9 | -0.0157 (0.0008) | -0.0257 (0.0008) | |
| scratch d = 32 | t\* | -0.0345 (0.0040) | -0.0411 (0.0046) | -0.0400 (0.0046) |
| scratch d = 32 | t = 9 | -0.0334 (0.0040) | -0.0406 (0.0043) | |
| scratch d = 128 | t\* | -0.0244 (0.0010) | -0.0313 (0.0009) | -0.0293 (0.0009) |
| scratch d = 128 | t = 9 | -0.0234 (0.0008) | -0.0309 (0.0006) | |

All-evidence truth, t\*, hub (pi): CLR-init d = 32 -0.0064 (0.0003), d = 128 -0.0065 (0.0008), scratch d = 32 -0.0236, d = 128 -0.0152 (every one 5 of 5 losses under the plan's thresholds; the interval of the CLR-init all-evidence contrasts is in `results/idea5/summary/contrasts_per_model.csv`). Half-trained models against their matched half-data CLR attention at t = 9 (SC, hub (pi), t\*; 10 models per config, informal count): CLR-init d = 32 -0.0173, d = 128 -0.0124, scratch d = 32 -0.0322, d = 128 -0.0238; all 10 of 10 are below -0.005 in each config (`passfail_half_matched.csv`).

## 4. Levels (pooled AUPR, mean over the 5 seeds [min-max]; CLR_full is the paper's attention)

Refits (full data) and untrained controls. Base rates SC 0.0943, all 0.1205.

**SC, raw**

| model | pairwise | t = 1 | t = 9 | t* |
|---|---|---|---|---|
| CLR_full | 0.1359 [0.1359-0.1359] | 0.1090 [0.1090-0.1090] | 0.1579 [0.1579-0.1579] | - |
| clr_d32_full | 0.1121 [0.1119-0.1122] | 0.1048 [0.1048-0.1049] | 0.1275 [0.1257-0.1297] | 0.1221 [0.1216-0.1226] |
| clr_d128_full | 0.1144 [0.1143-0.1148] | 0.1056 [0.1054-0.1056] | 0.1322 [0.1312-0.1333] | 0.1323 [0.1317-0.1328] |
| scratch_d32_full | 0.1061 [0.1035-0.1070] | 0.0980 [0.0976-0.0984] | 0.1174 [0.1106-0.1215] | 0.1168 [0.1099-0.1216] |
| scratch_d128_full | 0.1100 [0.1095-0.1109] | 0.0990 [0.0988-0.0992] | 0.1270 [0.1264-0.1280] | 0.1267 [0.1257-0.1279] |
| ctrl_scratch_d32 | 0.0945 [0.0940-0.0951] | 0.0943 [0.0942-0.0943] | 0.0965 [0.0926-0.0997] | - |
| ctrl_scratch_d128 | 0.0944 [0.0942-0.0947] | 0.0943 [0.0943-0.0943] | 0.0949 [0.0918-0.0968] | - |

**SC, hub (pi)**

| model | pairwise | t = 1 | t = 9 | t* |
|---|---|---|---|---|
| CLR_full | 0.1324 [0.1324-0.1324] | 0.1073 [0.1073-0.1073] | 0.1472 [0.1472-0.1472] | - |
| clr_d32_full | 0.1134 [0.1132-0.1135] | 0.1180 [0.1179-0.1182] | 0.1265 [0.1246-0.1284] | 0.1289 [0.1281-0.1293] |
| clr_d128_full | 0.1158 [0.1156-0.1160] | 0.1191 [0.1188-0.1193] | 0.1315 [0.1304-0.1324] | 0.1349 [0.1344-0.1354] |
| scratch_d32_full | 0.1066 [0.1041-0.1075] | 0.0935 [0.0924-0.0941] | 0.1137 [0.1074-0.1175] | 0.1127 [0.1065-0.1164] |
| scratch_d128_full | 0.1115 [0.1109-0.1124] | 0.0956 [0.0947-0.0973] | 0.1237 [0.1231-0.1250] | 0.1228 [0.1219-0.1244] |
| ctrl_scratch_d32 | 0.0951 [0.0942-0.0957] | 0.0925 [0.0878-0.0956] | 0.0964 [0.0926-0.0996] | - |
| ctrl_scratch_d128 | 0.0947 [0.0941-0.0954] | 0.0929 [0.0886-0.0978] | 0.0949 [0.0917-0.0967] | - |

**all, raw**

| model | pairwise | t = 1 | t = 9 | t* |
|---|---|---|---|---|
| CLR_full | 0.1648 [0.1648-0.1648] | 0.1330 [0.1330-0.1330] | 0.1733 [0.1733-0.1733] | - |
| clr_d32_full | 0.1421 [0.1418-0.1428] | 0.1299 [0.1297-0.1300] | 0.1482 [0.1470-0.1499] | 0.1457 [0.1452-0.1461] |
| clr_d128_full | 0.1438 [0.1435-0.1441] | 0.1302 [0.1301-0.1304] | 0.1510 [0.1503-0.1523] | 0.1520 [0.1514-0.1532] |
| scratch_d32_full | 0.1337 [0.1313-0.1351] | 0.1239 [0.1234-0.1242] | 0.1452 [0.1395-0.1488] | 0.1449 [0.1389-0.1488] |
| scratch_d128_full | 0.1375 [0.1368-0.1389] | 0.1247 [0.1246-0.1249] | 0.1526 [0.1515-0.1539] | 0.1524 [0.1510-0.1537] |
| ctrl_scratch_d32 | 0.1205 [0.1201-0.1212] | 0.1205 [0.1204-0.1205] | 0.1201 [0.1162-0.1228] | - |
| ctrl_scratch_d128 | 0.1205 [0.1201-0.1207] | 0.1205 [0.1205-0.1205] | 0.1211 [0.1185-0.1227] | - |

**all, hub (pi)**

| model | pairwise | t = 1 | t = 9 | t* |
|---|---|---|---|---|
| CLR_full | 0.1608 [0.1608-0.1608] | 0.1263 [0.1263-0.1263] | 0.1628 [0.1628-0.1628] | - |
| clr_d32_full | 0.1438 [0.1434-0.1444] | 0.1412 [0.1409-0.1413] | 0.1480 [0.1471-0.1493] | 0.1564 [0.1558-0.1567] |
| clr_d128_full | 0.1457 [0.1454-0.1459] | 0.1422 [0.1417-0.1426] | 0.1515 [0.1508-0.1526] | 0.1563 [0.1552-0.1575] |
| scratch_d32_full | 0.1340 [0.1319-0.1351] | 0.1142 [0.1131-0.1151] | 0.1403 [0.1352-0.1432] | 0.1392 [0.1343-0.1415] |
| scratch_d128_full | 0.1393 [0.1388-0.1400] | 0.1162 [0.1154-0.1182] | 0.1486 [0.1475-0.1494] | 0.1475 [0.1461-0.1486] |
| ctrl_scratch_d32 | 0.1213 [0.1203-0.1224] | 0.1191 [0.1132-0.1336] | 0.1200 [0.1159-0.1229] | - |
| ctrl_scratch_d128 | 0.1205 [0.1195-0.1212] | 0.1249 [0.1129-0.1369] | 0.1210 [0.1184-0.1227] | - |

## 5. Training and depth

| Quantity | CLR-init d = 32 | CLR-init d = 128 | scratch d = 32 | scratch d = 128 |
|---|---|---|---|---|
| Mean best epoch, direction h1 / h2 | 410 / 177 | 293 / 173 | 595 / 480 | 582 / 364 |
| Held-out R^2 at step 0 -> best (h1) | 0.615 -> 0.688 | 0.615 -> 0.689 | 0.293 -> 0.548 | 0.293 -> 0.576 |
| Held-out R^2 at step 0 -> best (h2) | 0.610 -> 0.684 | 0.610 -> 0.687 | 0.240 -> 0.522 | 0.240 -> 0.555 |
| Runs whose best epoch is at the cap (>= 595 of 600), of 10 | 0 | 0 | 3 | 3 |
| Refit epochs E (mean), steps | 294, 2350 | 234, 1869 | 538, 4302 | 473, 3787 |
| t\* (range over seeds) | 2 | 3-4 | 10-13 | 10-11 |
| Isolated rows (self-loop only) | 0 | 0 | 0 | 0 |

The half-CLR operators have t\* = 9 and 8 on this split (the project's depth). No t\* is at the scanned maximum of 32. Learning rate 1e-2 was the best of the grid for every config, so the grid edge is binding; 6 of 20 scratch half-trained runs end at the epoch cap (best epoch >= 595 in 3 per d), so scratch models may be under-trained under the plan's cap.

## 6. Post-hoc, label-free EcoMAC check (plan section 10)

Held-out reconstruction R^2 on the 1,370 EcoMAC arrays without an M3D counterpart (subset B, ComBat by platform), 4,056 shared genes, M3D-trained refit models; no truth scored (a trained operator does not depend on arrays, so scoring it against RegulonDB would repeat M3D). CLR attention at step 0 (full-data CLR): 0.311. CLR-initialised refits: d = 32 0.391 (0.390-0.392), d = 128 0.393 (0.392-0.394). Scratch refits: d = 32 0.258, d = 128 0.285. Untrained scratch layers: 0.131. Post hoc; shared genes only; the attention softmax is renormalised over the shared genes for every model.

## 7. Jobs

| Job | Wall time | Peak memory | Note |
|---|---|---|---|
| 5100_i5_prep | 47 s | 1.0 GB (declared 6) | R, 8 threads |
| 5101_i5_train_tune | 27 min | 2.5 GB (declared 8) | GPU 1 (GTX 1060), gate exports + learning-rate selection |
| 5102_i5_train_final | 64 min | 5.4 GB (declared 8) | GPU 1; 70 exported models; per-run 0.2-117 s |
| 5103_i5_score | 19 min, killed | 16.4 GB (declared 14, cap 18.2) | OOM kill during scoring (8 workers); scheduler auto-retry |
| 5103_i5_score_retry | 29 min | 21.1 GB (declared 25) | the run reported here (gates, 13 min of depth curves, 72 scoring tasks) |
| 5104_i5_ecomac_recon | 17 s | negligible | GPU 1 |

GPU memory was about 0.8 GB (nvidia-smi sample on device 1 during tuning).

## 8. Deviations from the plan

1. `5103_i5_score` was declared at 14 GB and was OOM-killed at 16.4 GB; the scheduler retried it at 25 GB. Nothing in the code changed; the retry reused the first attempt's gate table and depth curves (already written to disk) and redid the scoring.
2. `peak_gpu_mb` in the training logs is 0 (it queried device 0 rather than device 1); GPU memory is from nvidia-smi instead.
3. The plan said the seed-0 rerun difference would be reported; it was not rerun (no bit-identity claim is made).
4. `idea5_summarize.R` needed a fix after the run (JSON NaN in refit run records); no result changed. The matched-half pass/fail is an informal count over 10 models per config (the plan applies the 5-seed rule only to the refits).
5. Nothing else: architecture, loss, optimiser, schedule, cap, seeds, stopping rule, readouts, truths, jackknife stream and thresholds are as committed. Transient logits (6.7 GB) were deleted after scoring.

## 9. Reading, limited to the tables

- By the plan's rule the CLR-initialised learned model **loses** to the CLR attention at t = 9 on hub-stratified SC pooled AUPR in 5 of 5 seeds (-0.018, d = 32; -0.012, d = 128), with the seed-to-seed SD (0.0005) well below the jackknife SD (0.0027). The same sign holds for d = 128, scratch initialisation, t = 9, raw and MI-row-mean scoring, and the all-evidence truth.
- Training improves what it optimises and lowers what is scored: held-out reconstruction R^2 rises (0.61 -> 0.69 on the held-out lab half; 0.311 -> 0.391 on EcoMAC non-overlap arrays) while the pairwise score of the logits falls from 0.1324 to 0.113-0.116 (SC, hub (pi)), and the label-free depth drops from 9 to 2-4.
- Scratch models are above their untrained controls (SC pi-stratified, t = 9: 0.114 and 0.124 for d = 32 and 128 vs 0.096 and 0.095) and below the CLR attention (0.147); larger d is closer to CLR than smaller d in every comparison here.
- Caveats from the tables: learning rate at the grid edge, scratch runs near the epoch cap, one lab split, and the depth rule scores operators against the CLR of the held-out half.
