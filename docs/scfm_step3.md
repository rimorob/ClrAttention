# Learned-attention comparator, step 3: production extraction (2026-11-02)

Scope: production extraction of scPRINT (medium-v1.5 zero-shot, arms A and B; readouts R0, R0', R2, R1), the controls, scGPT whole-human, and the secondary scPRINT checkpoints (large-v1, small-v1). **No truth file was opened (BEELINE networks, CollecTRI, DoRothEA, STRING, ChIP, CORUM, Replogle) and nothing was scored.** Preregistration: `docs/preregistration_scfm.md`, change log entry of 2026-11-01 (E1 to E12, commit 17ca611); the step-3 freeze entry (E13) is committed together with this document. Scripts: `analysis/scfm_extract_scprint.py`, `scfm_extract_scgpt.py`, `scfm_extract_r1.py`, `scfm_r1_setup.py`, `scfm_controls.py`, `scfm_shim_crosscheck.py`, `scfm_gate_compare.py`, `scfm_secondary_loadtest.py`, `scfm_qc_dir.py`, `scfm_step3_summary.py`. Small results (QC tables, gate JSONs, the SHA-256 manifest of all 239 network/attention files): `results/scfm_step3/`. Large outputs (3.8 GB, not in the repo): `~/scfm_models/step3/`. Disk: `~/scfm_models` 9.3 GB of the 20 GB budget (it was 4.1 GB before step 3: 3.8 GB outputs, 1.35 GB for the two secondary checkpoints, 0.05 GB lamin instance); nothing pre-existing was deleted, only my own scratch test outputs.

## Result at a glance

| Step | Result |
|---|---|
| 3.0 change-log entry E1-E12 | committed and pushed (17ca611), before any production extraction |
| 3.1 production scripts | done; the unit check reproduces the step-1 hESC_1000 64-cell arrays **bitwise** for both scPRINT (2.2.0 route; max abs diff 0, shape 8x4x1161x1161) and scGPT (12x8x1221x1221) |
| 3.2 gates | all pass: (a) shim vs 2.2.0 on mHSC-E_500 arm A, 1,071 cells: **max abs diff 0.0, bitwise identical**, all 8 layers (tolerance 1e-6); (b) determinism hESC_500 twice: bitwise identical (scPRINT attention, R0, R0', R2, head sets; scGPT attention, both token cosines, hidden-state cosine); (c) large-v1 and small-v1 load strictly under 2.2.0; (d) GNInfer returns a G x G matrix |
| 3.3 scPRINT medium-v1.5, arm A, 14 cases, all cells | done, QC pass, 7 min (job 1003) |
| 3.4 arm B (GEO counts), 8 cases | done, QC pass, 5 min (job 1004) |
| 3.5 R1 (GNInfer), arm A 14 and arm B 8 cases | done, QC pass, 8 min (job 1005); 6 to 26 s per case |
| 3.6 controls | done, QC pass, 63 min (job 1006): permutation seeds 1001-1003; random weights seeds 1-3 in two versions, 14 cases each |
| 3.7 scGPT, 4 human cases, seed 20261021, plus seed-1 rerun of hESC_1000 | done, QC pass, 8 min (job 1007) |
| 3.8 K562 | **not run, blocked** by the no-truth rule (below) |
| 3.9 large-v1, small-v1, arm A, 14 cases, mean (R0) only | done, QC pass, 20 min (job 1008); checkpoints downloaded and hash-verified |
| 3.10 fine-tuning | not done; report below |
| 3.11 this document and E13 | done |

All GPU jobs ran one at a time through the scheduler (`jobs_pc/sched`, `# after:` chain 1001 to 1010), on device 1 (GTX 1060) only, nice 19 inside the per-job scope. Peak host memory 3.5 GB (declared 16 GB), GPU peak allocation 0.2 to 4.5 GB.

## What was run

Environment: scprint 2.2.0 in `~/venvs/scfm_alt` (torch 2.3.1+cu121, borrowed from `~/venvs/scfm`, unchanged); scGPT 0.2.4 in `~/venvs/scfm`. Checkpoints: medium-v1.5 `96ff7c5d...4ff8` (E1), large-v1 `c0c53989...3ef9` and small-v1 `269bc14b...d090` (both downloaded from HF jkobject/scPRINT commit 8e1ad1a1, SHA-256 equal to the HF LFS metadata in E1); scGPT as in E1.

- **Input.** All cells of each case; genes = the case genes with gene-map status ok (scPRINT) or covered by the scGPT vocabulary. Arm A = `2^y - 1` of the BEELINE matrix (step-2 `tpm_beeline`); arm B = GEO counts on genes with counts. Zeros are tokens. The scPRINT encoder applies `log2(1+x)`; no depth token for medium-v1.5 (`depth_atinput` False, 9 prefix tokens).
- **R0.** Per cell, layer and head, softmax of Q K^T / sqrt(d) over gene keys only (identical to stripping the special tokens and renormalising each row), fp32, batch 8, averaged over cells; stored per layer x head (float32, 8 x 4 x G x G) with R0 (mean over all 32 heads), R0' and R2 (last layer, head mean) in `<case>.nets.npz`. Orientation (S = A^T etc.) is applied at scoring.
- **R0'.** Heads with cell-averaged special-token mass < 0.5 (no case needed the 4-lowest fallback; 12 to 15 of 32 heads pass).
- **R1.** `scprint.tasks.GNInfer` (scprint 2.2.0) with the E3 settings: how=given, genes = covered case genes, num_genes = max_len = G, max_cells=0, drop_unexpressed=False, layer=None, head_agg=mean, preprocess=softmax, filtration=none, apc=False, symmetrize=False, known_grn=None, precision "32", dtype float32, forward_mode none, doplot False, comp_attn default, one constant cell-type label; batch_size 8 and num_workers 0 (unlisted in E3). The full kwargs are in every `r1_*/<case>.r1.json` (example: `results/scfm_step3/gates/r1_example_hESC_500.r1.json`). The output gene order equals the vocabulary order and is mapped back to the gene-map order.
- **Controls** (arm A, R0 only): permutation of each gene independently across all cells (seeds 1001-1003); random-weight models built by `scPrint(**checkpoint hyper-parameters)` after `torch.manual_seed(s)` (seeds 1-3), either with the pretrained gene embeddings copied in (`rand_esm`) or with the package's random gene-embedding initialisation (`rand_emb`). ESM2-derived gene-embedding cosine of the checkpoint (data-free prior) for all 14 cases: `~/scfm_models/step3/esm2_cos/`.
- **scGPT**: binning over the case genes, `np.random.seed(20261021)` immediately before the binning of every case; attention of the last layer per head (stored), mean over heads, `<cls>` removed and rows renormalised; token-embedding cosine (raw embedding rows and after the encoder LayerNorm) and final-layer hidden-state cosine averaged over cells.

## QC (all truth-free)

scPRINT medium-v1.5, arm A (`results/scfm_step3/qc_scprint_A.csv`). `split_*` = cells split once with seed 20261021, relative L1 between the two half-means (off-diagonal); `org_acc` = fraction of cells whose organism head predicts the right taxon.

| case | n_cells | G | org_acc_A | org_acc_B | r0p_n_heads | heads_mass_gt_0p9 | mass_median | split_R0_relL1 | split_R0_spearman | split_R0p_relL1 | split_R2_relL1 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| hESC_500 | 758 | 842 | 1.000 | 1.000 | 13 | 18 | 0.943 | 0.108 | 0.977 | 0.166 | 0.272 |
| hESC_1000 | 758 | 1288 | 1.000 | 1.000 | 13 | 16 | 0.908 | 0.111 | 0.974 | 0.180 | 0.285 |
| hHep_500 | 425 | 906 | 1.000 | nan | 12 | 16 | 0.852 | 0.154 | 0.948 | 0.242 | 0.424 |
| hHep_1000 | 425 | 1362 | 1.000 | nan | 13 | 16 | 0.850 | 0.170 | 0.942 | 0.260 | 0.499 |
| mDC_500 | 383 | 798 | 0.995 | nan | 14 | 8 | 0.634 | 0.143 | 0.944 | 0.214 | 0.357 |
| mDC_1000 | 383 | 1247 | 0.995 | nan | 13 | 9 | 0.675 | 0.159 | 0.932 | 0.240 | 0.401 |
| mESC_500 | 421 | 941 | 1.000 | nan | 14 | 13 | 0.668 | 0.100 | 0.981 | 0.140 | 0.238 |
| mESC_1000 | 421 | 1306 | 1.000 | nan | 15 | 12 | 0.618 | 0.109 | 0.977 | 0.139 | 0.247 |
| mHSC-E_500 | 1071 | 660 | 0.668 | 0.685 | 14 | 12 | 0.722 | 0.082 | 0.986 | 0.119 | 0.223 |
| mHSC-E_1000 | 1071 | 1122 | 0.988 | 0.995 | 14 | 12 | 0.681 | 0.093 | 0.983 | 0.136 | 0.253 |
| mHSC-GM_500 | 889 | 584 | 0.984 | 0.994 | 13 | 10 | 0.683 | 0.084 | 0.985 | 0.132 | 0.264 |
| mHSC-GM_1000 | 889 | 1030 | 1.000 | 1.000 | 13 | 10 | 0.602 | 0.096 | 0.982 | 0.149 | 0.268 |
| mHSC-L_500 | 847 | 489 | 0.995 | 0.999 | 14 | 12 | 0.785 | 0.088 | 0.984 | 0.139 | 0.248 |
| mHSC-L_1000 | 847 | 647 | 1.000 | 1.000 | 14 | 12 | 0.720 | 0.089 | 0.984 | 0.133 | 0.223 |

Notes. No NaN or Inf anywhere; row sums of every attention tensor are 1 within 2.4e-7. 14 to 18 of 32 heads put more than 90% of a gene query's mass on the special tokens in the human cases and 8 to 13 in the mouse cases (all values in `sink_mass_A.csv`); layer 1 head 4 and layer 8 head 3 are at about 1.0 in every case. **Flag:** the organism-head accuracy of mHSC-E_500 is 0.67 under arm A and 0.69 under arm B (33% of its cells are classified human); every other case is 0.99 to 1.00 under both arms, including the superset case mHSC-E_1000 (same cells, 1,122 genes). The shim parity for this very case is exact, so this is not an extraction error; it reflects the model's head on a 660-gene input. I did not interpret it further and applied no rule to it. For the secondary checkpoints the organism head is 1.0 in all cases.

Arm A vs arm B (8 cases, common genes; `a_vs_b.csv`). Arm B has one gene fewer in mHSC-E_1000 and mHSC-GM_1000 (a covered gene without deposited counts).

| case | G_A | G_B | G_common | R0_relL1 | R0_spearman | R0p_relL1 | R0p_spearman | R2_relL1 | R2_spearman |
|---|---|---|---|---|---|---|---|---|---|
| hESC_500 | 842 | 842 | 842 | 0.038 | 0.996 | 0.057 | 0.991 | 0.123 | 0.973 |
| hESC_1000 | 1288 | 1288 | 1288 | 0.041 | 0.995 | 0.063 | 0.989 | 0.137 | 0.973 |
| mHSC-E_500 | 660 | 660 | 660 | 0.071 | 0.991 | 0.089 | 0.984 | 0.221 | 0.957 |
| mHSC-E_1000 | 1122 | 1121 | 1121 | 0.074 | 0.991 | 0.097 | 0.984 | 0.235 | 0.960 |
| mHSC-GM_500 | 584 | 583 | 583 | 0.068 | 0.993 | 0.091 | 0.986 | 0.219 | 0.965 |
| mHSC-GM_1000 | 1030 | 1029 | 1029 | 0.076 | 0.991 | 0.103 | 0.983 | 0.233 | 0.959 |
| mHSC-L_500 | 489 | 489 | 489 | 0.079 | 0.990 | 0.115 | 0.978 | 0.223 | 0.947 |
| mHSC-L_1000 | 647 | 647 | 647 | 0.086 | 0.988 | 0.125 | 0.975 | 0.227 | 0.948 |

R1 vs the R0 family (Spearman over off-diagonal entries; `r1_vs_r0.csv`):

| case | G | R1_vs_R0_spearman | R1_vs_R0p_spearman | R1_vs_R2_spearman |
|---|---|---|---|---|
| hESC_500 | 842 | 0.337 | 0.398 | 0.082 |
| hESC_1000 | 1288 | 0.367 | 0.450 | 0.097 |
| hHep_500 | 906 | 0.317 | 0.405 | 0.084 |
| hHep_1000 | 1362 | 0.328 | 0.433 | 0.101 |
| mDC_500 | 798 | 0.439 | 0.555 | 0.349 |
| mDC_1000 | 1247 | 0.435 | 0.540 | 0.363 |
| mESC_500 | 941 | 0.508 | 0.641 | 0.194 |
| mESC_1000 | 1306 | 0.562 | 0.660 | 0.209 |
| mHSC-E_500 | 660 | 0.502 | 0.549 | 0.189 |
| mHSC-E_1000 | 1122 | 0.497 | 0.544 | 0.184 |
| mHSC-GM_500 | 584 | 0.495 | 0.545 | 0.207 |
| mHSC-GM_1000 | 1030 | 0.502 | 0.537 | 0.210 |
| mHSC-L_500 | 489 | 0.410 | 0.440 | 0.101 |
| mHSC-L_1000 | 647 | 0.430 | 0.412 | 0.126 |

Controls (`control_distances.csv`; relative L1 and Spearman against the zero-shot R0 of arm A; the first columns are shown for seeds 1001 / 1 only). Over all cases and seeds: permutation relL1 0.069 to 0.326, Spearman 0.907 to 0.990; random weights relL1 0.726 to 1.090, Spearman -0.121 to 0.143.

| case | G | relL1_vs_uniform | perm1001_relL1 | perm1001_spearman | rand_esm1_relL1 | rand_esm1_spearman | rand_emb1_relL1 | rand_emb1_spearman |
|---|---|---|---|---|---|---|---|---|
| hESC_500 | 842 | 0.868 | 0.130 | 0.977 | 0.868 | -0.012 | 0.868 | -0.001 |
| hESC_1000 | 1288 | 0.903 | 0.134 | 0.977 | 0.903 | -0.015 | 0.903 | -0.004 |
| hHep_500 | 906 | 0.908 | 0.196 | 0.960 | 0.908 | 0.081 | 0.908 | 0.007 |
| hHep_1000 | 1362 | 0.925 | 0.198 | 0.960 | 0.924 | 0.112 | 0.925 | 0.003 |
| mDC_500 | 798 | 0.938 | 0.304 | 0.909 | 0.938 | 0.109 | 0.938 | -0.001 |
| mDC_1000 | 1247 | 1.015 | 0.325 | 0.909 | 1.014 | 0.090 | 1.015 | -0.003 |
| mESC_500 | 941 | 1.051 | 0.168 | 0.964 | 1.051 | 0.082 | 1.051 | -0.008 |
| mESC_1000 | 1306 | 1.089 | 0.179 | 0.963 | 1.089 | 0.084 | 1.089 | -0.003 |
| mHSC-E_500 | 660 | 0.768 | 0.139 | 0.957 | 0.767 | 0.083 | 0.768 | -0.012 |
| mHSC-E_1000 | 1122 | 0.800 | 0.145 | 0.962 | 0.799 | 0.143 | 0.800 | -0.005 |
| mHSC-GM_500 | 584 | 0.732 | 0.085 | 0.988 | 0.731 | 0.108 | 0.732 | -0.008 |
| mHSC-GM_1000 | 1030 | 0.763 | 0.098 | 0.985 | 0.762 | 0.123 | 0.763 | -0.007 |
| mHSC-L_500 | 489 | 0.727 | 0.072 | 0.989 | 0.726 | 0.089 | 0.727 | 0.002 |
| mHSC-L_1000 | 647 | 0.854 | 0.073 | 0.989 | 0.853 | 0.078 | 0.854 | -0.003 |

Reading, descriptive. (i) Permuting each gene independently across cells moves R0 by relL1 0.07 to 0.33, with Spearman 0.91 to 0.99 against the unpermuted R0. The yardstick is the split-half distance (two half-means of the same data, so each half has only half the cells): the permutation distance is 0.8 to 2.1 times the split-half distance (about 1.0 in mHSC-GM, 0.8 in mHSC-L, 1.2 in hESC and hHep, 1.6 to 1.7 in mESC and mHSC-E, 2.0 to 2.1 in mDC). So the expression-dependent component of scPRINT's mean attention is small against the gene-embedding prior, is detectable beyond sampling noise in mDC, mESC and mHSC-E, and is not distinguishable from it in mHSC-GM and mHSC-L; the permutation moves R0 by 7% to 33%, never by more than a third. (ii) Random-weight networks are close to uniform (relL1 against the trained R0 equals the trained network's distance from uniform, 0.73 to 1.09; Spearman about 0 for random gene embeddings and -0.12 to 0.14 for the kept ESM2 embeddings). The random gene embeddings are the package's `nn.Embedding` initialisation (std 1.0) against std 0.069 for the pretrained ones; this scale difference was not tuned.

scGPT (`scgpt_qc.csv`). The `<cls>` mass is at most 0.0012 in the last layer (max over heads) and at most 0.234 in any layer. The seed-1 rerun of hESC_1000 differs from the production seed 20261021 by relative L1 0.0023 for the last-layer head-mean attention (Spearman 1.0000), 0.0004 for the all-layer mean and 0.0001 for the hidden-state cosine.

| case | n_cells | G | seed | infer_s | cls_mass_last_layer_max | cls_mass_max_any_layer | tok_cos_vs_tok_cos_ln_spearman | attn_last_vs_hid_cos_spearman | seed1_vs_prod_attn_last_mean_relL1 | seed1_vs_prod_attn_last_mean_spearman |
|---|---|---|---|---|---|---|---|---|---|---|
| hESC_500 | 758 | 841 | 20261021 | 57.800 | 0.001 | 0.234 | 0.973 | -0.057 | nan | nan |
| hESC_1000 | 758 | 1285 | 20261021 | 126.900 | 0.000 | 0.187 | 0.972 | -0.066 | 0.002 | 1.000 |
| hHep_500 | 425 | 907 | 20261021 | 37.400 | 0.001 | 0.199 | 0.975 | -0.136 | nan | nan |
| hHep_1000 | 425 | 1362 | 20261021 | 78.100 | 0.001 | 0.160 | 0.975 | -0.177 | nan | nan |

## R0' head sets (E13; computed from expression only; layer L, head H, 1-based)

| case | n_heads | rule | heads_layer_head_1based |
|---|---|---|---|
| hESC_500 | 13 | mass<0.5 | L1H1;L1H2;L1H3;L2H3;L3H2;L3H3;L4H1;L4H2;L5H1;L7H1;L7H3;L8H1;L8H2 |
| hESC_1000 | 13 | mass<0.5 | L1H1;L1H2;L1H3;L2H3;L3H2;L3H3;L4H1;L4H2;L5H1;L7H1;L7H3;L8H1;L8H2 |
| hHep_500 | 12 | mass<0.5 | L1H1;L1H2;L1H3;L2H3;L3H2;L4H1;L4H2;L5H1;L7H1;L7H3;L8H1;L8H2 |
| hHep_1000 | 13 | mass<0.5 | L1H1;L1H2;L1H3;L2H3;L3H2;L4H1;L4H2;L5H1;L6H3;L7H1;L7H3;L8H1;L8H2 |
| mDC_500 | 14 | mass<0.5 | L1H1;L1H2;L1H3;L2H3;L3H2;L3H4;L4H1;L4H2;L5H1;L6H3;L7H1;L7H3;L8H1;L8H2 |
| mDC_1000 | 13 | mass<0.5 | L1H1;L1H2;L1H3;L2H3;L3H2;L4H1;L4H2;L5H1;L6H3;L7H1;L7H3;L8H1;L8H2 |
| mESC_500 | 14 | mass<0.5 | L1H1;L1H2;L1H3;L2H3;L2H4;L3H2;L4H1;L4H2;L5H1;L5H2;L7H1;L7H3;L8H1;L8H2 |
| mESC_1000 | 15 | mass<0.5 | L1H1;L1H2;L1H3;L2H1;L2H3;L2H4;L3H2;L4H1;L4H2;L5H1;L5H2;L7H1;L7H3;L8H1;L8H2 |
| mHSC-E_500 | 14 | mass<0.5 | L1H1;L1H2;L1H3;L2H3;L2H4;L3H2;L3H3;L4H1;L4H2;L5H1;L7H1;L7H3;L8H1;L8H2 |
| mHSC-E_1000 | 14 | mass<0.5 | L1H1;L1H2;L1H3;L2H3;L2H4;L3H2;L3H3;L4H1;L4H2;L5H1;L7H1;L7H3;L8H1;L8H2 |
| mHSC-GM_500 | 13 | mass<0.5 | L1H1;L1H2;L1H3;L2H3;L2H4;L3H2;L4H1;L4H2;L5H1;L7H1;L7H3;L8H1;L8H2 |
| mHSC-GM_1000 | 13 | mass<0.5 | L1H1;L1H2;L1H3;L2H3;L2H4;L3H2;L4H1;L4H2;L5H1;L7H1;L7H3;L8H1;L8H2 |
| mHSC-L_500 | 14 | mass<0.5 | L1H1;L1H2;L1H3;L2H3;L3H2;L3H3;L4H1;L4H2;L5H1;L7H1;L7H3;L8H1;L8H2;L8H4 |
| mHSC-L_1000 | 14 | mass<0.5 | L1H1;L1H2;L1H3;L2H3;L3H2;L3H3;L4H1;L4H2;L5H1;L7H1;L7H3;L8H1;L8H2;L8H4 |

The same table is `results/scfm_step3/r0p_head_sets.csv`; the per-layer x head special-token mass and entropy of every case are in `sink_mass_A.csv`. The sets are stable: the 12 heads L1H1, L1H2, L1H3, L2H3, L3H2, L4H1, L4H2, L5H1, L7H1, L7H3, L8H1 and L8H2 are in every one of the 14 sets; the others vary by case (L2H1, L2H4, L3H3, L3H4, L5H2, L6H3, L8H4).

## Unexpected findings and deviations (all logged here and in E13)

1. **`populate_my_ontology()` raised a TypeError** in scdataloader 1.0.1 (`records[-4][0]` in its mouse developmental-stage step). The cell-type, organism, sex, ethnicity, assay, tissue and developmental-stage steps before it had completed in the separate lamin instance `~/scfm_models/lamin_r1` (isolated with `LAMIN_SETTINGS_DIR`, so `~/.lamin` is untouched). The mouse developmental-stage step was skipped (GNInfer reads only the organism), and the Disease and Gene steps that follow it were run by `analysis/scfm_r1_setup.py`; the Gene registry had to be filled in chunks of 400 because sqlite fails with "too many SQL variables", and only with the 44,756 Ensembl ids of the scPRINT vocabulary (all 44,756 are in bionty's public Ensembl table, none missing). GNInfer's gene order was checked against the vocabulary order (equal in all cases).
2. **GNInfer sends batches to the bare device string `"cuda"`** (the current device), so a first test with an explicit `cuda:1` model crashed on a cuda:0/cuda:1 mismatch. That single failed test touched device 0 for a moment (a few MiB, process ended with an error; device 0 shows its normal ~1.0 GB desktop usage). The R1 runner now requires `CUDA_VISIBLE_DEVICES=1` (physical device 0 invisible) and asserts that the only visible card is the GTX 1060; the 22 production R1 runs were made that way.
3. **R1 row sums are 0.25 to 0.62**, not 1: the package drops the special-token rows and columns after the softmax without renormalising (the R0 definition renormalises). This is the package behaviour fixed in E3, not changed.
4. **large-v1 and small-v1 need a depth token** (`depth_atinput` True; 8 prefix tokens; 16 x 4 heads, d=512 and 4 x 2 heads, d=128; `attn_bias` none, no cell transformer; vocabulary identical to medium-v1.5). The load-test forward in step 3.2 first ran without the depth token, was corrected, and was rerun; production feeds `req_depth` = per-cell sum of the model input over the case genes (what scdataloader's Collator computes for the adata it is given). This choice is not in the preregistration; the secondary models are descriptive.
5. **npz files contain string object arrays** (`genes`, `ensembl`); load with `allow_pickle=True`. The network arrays themselves are float32 and need no pickling. The files were not rewritten, to keep the hashes below.
6. Arm B uses one gene fewer than arm A in mHSC-E_1000 (1,121 vs 1,122) and mHSC-GM_1000 (1,029 vs 1,030): a covered gene is absent from the counts deposit. Scoring uses each contrast's G-intersection anyway.
7. The control permutation and random-weight models store only R0 (the preregistered control readout) to stay inside the disk budget.
8. mHSC-E_500 organism-head accuracy (see QC): flagged, not a gate.

## Not done

**3.8 K562 is blocked by the no-truth rule.** The preregistered gene set is "CollecTRI TFs with >= 10 targets detected in >= 5% of control cells, plus the 1,000 most variable detected genes". The first part needs the CollecTRI regulon sizes (`data/human_truth/collectri.tsv`, a truth file under this task's rules). `results/baselines_k562ess/expression.tsv` cannot be reused: its columns are perturbation pseudobulks (`10_AARS2_P1P2_...`), not control cells. The raw Replogle single-cell file is also a truth source (knockdown-response modules) and was not opened. A superset gene list cannot replace the preregistered one, because attention depends on the token set. A PI decision is needed, for example: allow reading CollecTRI only to count targets per TF for the gene-set definition (no scoring), or provide the gene list. The extractors take any gene list, so the run itself is about 20 GPU-minutes once the list exists.

**3.10 Fine-tuning (report only, not started).** What it needs: scPRINT 2.2.0 trains through `scprint fit` (LightningCLI) on a scdataloader `DataModule` that reads a lamindb `Collection`; the sqlite lamin instance fails on bulk gene registration (see deviation 1), so it needs either a postgres-backed instance or a custom `LightningDataModule` that feeds our cells with the same Collator. "Package default losses" for scPRINT are the ZINB expression loss, classification and, with `mvc_decoder` set, the MVC loss (the medium checkpoint has `inner product`), with `lr` 1e-4; the 10-epoch, seed 20261021, no-truth setting of E10 is expressible, and the 7 datasets (4,794 cells per gene set) would take minutes per epoch on the 1060. For scGPT there is no packaged default recipe in the repository beyond its tutorials (GEP and GEPC losses), so "package default losses" would have to be defined by us. Estimate: about 1 day of engineering, 2 to 4 GPU-hours, as in the draft; the main risk is the data module. It does not block scoring of the primary and may finish after it.

## SHA-256 of the network and attention files

All 239 files (`.npy`, `.npz` of every output directory, including the controls, the secondary checkpoints and the ESM2 cosines) are listed with size and SHA-256 in `results/scfm_step3/sha256_manifest.csv`; the SHA-256 of that manifest is `e27a8422f5cf868459c694407a0cd0e900a18ce8b3aaf0c33ee179260c2e3f15`. The primary and descriptive network files (relative to `~/scfm_models/step3/`):

| path | sha256 |
|---|---|
| ctrl_perm_s1001/hESC_1000.nets.npz | 6189a707259cd0938dd93c87288f416243f024cc7de62fc79cbac45459e5db39 |
| ctrl_perm_s1001/hESC_500.nets.npz | c4f969698e06aa0e2a60b5fbd1a46f78d3adf527bb3881fa2fe7f6ce43b6cb6c |
| ctrl_perm_s1001/hHep_1000.nets.npz | 9dc7667128ec99f1401c801769d1f344a4eda05e5ac5693b12da72b192d03583 |
| ctrl_perm_s1001/hHep_500.nets.npz | a206ea3bcc350660199983c613ee2eaaaad9b7570b3902f4d9414abb2a317995 |
| ctrl_perm_s1001/mDC_1000.nets.npz | b42a826ce3810d05af7fe7a34c0194be9c5dcde5d71f1be564a067e1d6bff7ee |
| ctrl_perm_s1001/mDC_500.nets.npz | ec50614bec37730d409003469324f7772e690fa9ef7c9c8d5acd10c3a2d216a3 |
| ctrl_perm_s1001/mESC_1000.nets.npz | e1e0c3f20936087a61da750110dff0732ba99b1e14b5663c4dc8d37653387021 |
| ctrl_perm_s1001/mESC_500.nets.npz | 0ea9113ab519a673fc0a3884db31a22d9f6f94bc72deddd7095fa6efa1192a60 |
| ctrl_perm_s1001/mHSC-E_1000.nets.npz | fcb42d0663084b59205212fa6c3a89c9e864e96d52471a36973b143b7e028246 |
| ctrl_perm_s1001/mHSC-E_500.nets.npz | acb622e4865c8c78c1ec5b7287fd45b1d6265909b8fd084fc5da8bf078ce5878 |
| ctrl_perm_s1001/mHSC-GM_1000.nets.npz | 21bcac31bedccd5fdc64f901c0822d16deaf2b678f91409ea141abfb8a45f427 |
| ctrl_perm_s1001/mHSC-GM_500.nets.npz | a48157e8aa2b63dca7f9f5790153efd9854c4e264f7abc8224de19cc368a2a0a |
| ctrl_perm_s1001/mHSC-L_1000.nets.npz | 9b78f8041fee1a5325873e06c0b1aa5b985c58e38187ea5cecfff72826bab36f |
| ctrl_perm_s1001/mHSC-L_500.nets.npz | 5e318465bca65603c07c0251f36468f256eb16b32abe1039b8c62c40446f1d92 |
| ctrl_perm_s1002/hESC_1000.nets.npz | 9b18787f9c60aabcf7a7c618dd0469a833258b5eaca0e54edba174bc8d9c8b04 |
| ctrl_perm_s1002/hESC_500.nets.npz | 55659d293fbfbba87e7702b5aeec1424832f611d273f57b6b9390ad204cc8fad |
| ctrl_perm_s1002/hHep_1000.nets.npz | a2b6cdf2f35d074fef32ef3ffb8ebf25b8d08ce806a2cb2f99c4e8efe80c772f |
| ctrl_perm_s1002/hHep_500.nets.npz | 90476c4b4bfa76f3cf46bdbcf990cc842ba3cc0babf5bf3ffcc188ab371b956c |
| ctrl_perm_s1002/mDC_1000.nets.npz | 879c8f73bd7d98a80a836c9c6642170c6397baa8b434a2efae32213d6224c9b6 |
| ctrl_perm_s1002/mDC_500.nets.npz | ff0327a62ed4a21e00c252f38e8cd522fe9ff2e0e334f946be1eab2cbae61f75 |
| ctrl_perm_s1002/mESC_1000.nets.npz | 6a06b8fc6c5c254fc993079fdcdd618871f696d7e6c19fa070e57f16caa1bfe6 |
| ctrl_perm_s1002/mESC_500.nets.npz | 15dc1ff61ac51c629cda84cfdd8e7ce8f13676dfe2944e8579646e67f64949c2 |
| ctrl_perm_s1002/mHSC-E_1000.nets.npz | ec49c59767e27171d06e564ef9017e0b3aa7ba2d43101c1c8cd26d6c70683f89 |
| ctrl_perm_s1002/mHSC-E_500.nets.npz | 0267e6085c5bf563830bf819006a70ea838bc05bca7dfaa2b9811c620ae8f7de |
| ctrl_perm_s1002/mHSC-GM_1000.nets.npz | b08656de417cb1f3a9cd36cdfd115afd96be62cbb265de87474cd69666357831 |
| ctrl_perm_s1002/mHSC-GM_500.nets.npz | e6cdc82a44ef78a2bf14ee7105b14fc5fe880c2653670db91c654f7949ceca39 |
| ctrl_perm_s1002/mHSC-L_1000.nets.npz | d4e014abd81bfa1dde71babbaec7d93b11e931143e232e11642bb24d2f82f49a |
| ctrl_perm_s1002/mHSC-L_500.nets.npz | 64f27a171ceb7d59b50b2f5563bb0dd5886d62c52efc0caba92c27205e7e7dbc |
| ctrl_perm_s1003/hESC_1000.nets.npz | ba75488f53deb6fcbd885843da7fb12e12bb95eb7ec4b669eb98cede32d060c8 |
| ctrl_perm_s1003/hESC_500.nets.npz | 87029227b50bc3e042d59fa06069878c416e588acf60e4a97e9e252ff3ce2314 |
| ctrl_perm_s1003/hHep_1000.nets.npz | c55af4a36de240fa7d80587bfbc09e22752ae693ad9d4ce6efe46c4ee521fdb2 |
| ctrl_perm_s1003/hHep_500.nets.npz | bc5cd4d932dea9c10932c32ad35d7e25a987df91887a0feb4c5732a0110f60be |
| ctrl_perm_s1003/mDC_1000.nets.npz | 9ac841f4338a90435f00695dcbfefc6e26da9269c49668a0dcb9032a8f8d75c8 |
| ctrl_perm_s1003/mDC_500.nets.npz | 0c879ceff3f70440675a2d7e31059867f82dd74c51fa8428097fb6340147c2f0 |
| ctrl_perm_s1003/mESC_1000.nets.npz | fb00392339bbfea4db5746dee9d2994087306df98bed6e163514903a45a9d4cf |
| ctrl_perm_s1003/mESC_500.nets.npz | 0c6cee6d47141ff01a1bf41aab2795ccecd9c933af0cc8e3e838c204f749042f |
| ctrl_perm_s1003/mHSC-E_1000.nets.npz | 93ba58f9ff2b0885a0ceb87f6ee2ae887a84c300c5dff6689c161f75c6e39151 |
| ctrl_perm_s1003/mHSC-E_500.nets.npz | 700dfc14639d13f012dd314e7613f373c66cc1c4ea3cc146c700a3213f748172 |
| ctrl_perm_s1003/mHSC-GM_1000.nets.npz | 49efcc30f69d5e5655420cef9fdb7a079be87bd5d9d672f03aba6e2d517b4eff |
| ctrl_perm_s1003/mHSC-GM_500.nets.npz | 78d2d3c5c9184829ffac9cd8b1d0feff4bc0e469993ecb365d3c8965b0e2c92c |
| ctrl_perm_s1003/mHSC-L_1000.nets.npz | 5be6daa3a635e4e00bd791af664f18e7c70c44f95356111fa4284a2ab019cbe8 |
| ctrl_perm_s1003/mHSC-L_500.nets.npz | 0c16b28af4fe3c1ac797c77711ccbac1845ce85c1663eac7033e1c1014dd0859 |
| ctrl_rand_emb_s1/hESC_1000.nets.npz | 3e95936ab72c1fff9d85eb71c5c55050aba6767ec2e0a140e9c85201dd1c6511 |
| ctrl_rand_emb_s1/hESC_500.nets.npz | 3c327837089c2cf37bbdaebd62e650d3bc693654ae42b3417d6e203302ed380d |
| ctrl_rand_emb_s1/hHep_1000.nets.npz | 519b48d181d532fb2ac02e99c540547cd25d1fe80b6fa3014e39c9b1d1ffb306 |
| ctrl_rand_emb_s1/hHep_500.nets.npz | 564b35e785f76cf2eb86db7fa13ab389bf828c40f70f1c04e3c9ae3a64990a3e |
| ctrl_rand_emb_s1/mDC_1000.nets.npz | 6dc25ef2dcd59120ad9587bd44f650fddb2c6793d51ebc083e3f4fc1fb3fdd87 |
| ctrl_rand_emb_s1/mDC_500.nets.npz | 3d57fd59a132bdbe1489b58ca62713c2519cdc9908227d8c29f4affee204b113 |
| ctrl_rand_emb_s1/mESC_1000.nets.npz | d2d23fc461d32b9e0a23259e3d1f45f64151cb1e99f0e2cb14c11a7c8d2e7006 |
| ctrl_rand_emb_s1/mESC_500.nets.npz | 77ece5809b868cef7d1bbf1b8d7c39530f220af0465e2cfe8032429c7912f688 |
| ctrl_rand_emb_s1/mHSC-E_1000.nets.npz | b546a4dcc84f9dee812a201571ac278af784b2d6a05bf4a5a87cdb5edcfe7b12 |
| ctrl_rand_emb_s1/mHSC-E_500.nets.npz | 95624a6c81987158dc0d8e9ff6aed7f7b1384788b2b053f85e3b4fdb6b3bb1f3 |
| ctrl_rand_emb_s1/mHSC-GM_1000.nets.npz | 689f2c0467f6172c540c56711882e8fa1bb37ce64e69fab178e75add5d15c665 |
| ctrl_rand_emb_s1/mHSC-GM_500.nets.npz | ae583c0a2eb3517bb09bfd61817134ee90c61e043cb8a45f775b60ea9b03aa17 |
| ctrl_rand_emb_s1/mHSC-L_1000.nets.npz | 6ad864db517328414f2fcdf94bab541acc7d9e5de836617e1f4f8491c658b9c3 |
| ctrl_rand_emb_s1/mHSC-L_500.nets.npz | 535f29f759ef37679c8a5d12eae9847910d5195908e2c8258470575010292fe3 |
| ctrl_rand_emb_s2/hESC_1000.nets.npz | fed201530cd48079298e0a26e14d830d7b9305fd9ca6a6ea0f6fefb5057d51ad |
| ctrl_rand_emb_s2/hESC_500.nets.npz | 580baff4f13a1a60117c889237b94358383c8f96eaa34a2df43f1d8711343abf |
| ctrl_rand_emb_s2/hHep_1000.nets.npz | 3327baff68601fb07e59af6e491359a785f4f9a18adb70dd338399879df25109 |
| ctrl_rand_emb_s2/hHep_500.nets.npz | b88d1ebfc9db33e6d734e546383cf456973ed6c7e5a6d317dd2ba87d72d02e85 |
| ctrl_rand_emb_s2/mDC_1000.nets.npz | 01fd2f1c7f5577e15bfae1c570e27f1cd764b38eec56f4d5697214d2779b974e |
| ctrl_rand_emb_s2/mDC_500.nets.npz | 6f2625963141a1d96718388bc9dd9870f57a5aa9675393a9b5ae118cd351b11d |
| ctrl_rand_emb_s2/mESC_1000.nets.npz | a392ce7b7f3094f1bbab17ca65661a0689e56ea835f9fff23418f00e78db9358 |
| ctrl_rand_emb_s2/mESC_500.nets.npz | 182d00ac7039b2f04d4bb738c0e9a4545b3f04b0f5e7430a1ba94d76ac89c6e8 |
| ctrl_rand_emb_s2/mHSC-E_1000.nets.npz | 3850a8dc728b095aca72d8078e26b6751f08a929b56c06693deeb62201298e92 |
| ctrl_rand_emb_s2/mHSC-E_500.nets.npz | 80ca7cc0d54162e2a49ba07482f93975beddff23888257552e52fc760b560b1c |
| ctrl_rand_emb_s2/mHSC-GM_1000.nets.npz | ec465dca2b80b1ba9b8cb1e5dd47acae859f8eba290cf8b58e72cb6fe71a0352 |
| ctrl_rand_emb_s2/mHSC-GM_500.nets.npz | 12bc3eaa7c48339d2de682b5b9fb201a74f20d3718d5632bdfc5ea213aef0ded |
| ctrl_rand_emb_s2/mHSC-L_1000.nets.npz | 8d9b4a38ca08c93e09001bf7f001a0aee0e09acdb588a39f6d4b8db6e7f81f67 |
| ctrl_rand_emb_s2/mHSC-L_500.nets.npz | c1bd1ff66e4a4e1e35dd6af787760f2e903ab7e4547e3590392ecd1a9b8b8cd6 |
| ctrl_rand_emb_s3/hESC_1000.nets.npz | 526c2fa5c9f0f61852c7c60e7f2cdf386c007905f5e224ca70e329c2cee11341 |
| ctrl_rand_emb_s3/hESC_500.nets.npz | 3e2739a2a7c76d789a7d70f550cc5aa926429ea9a3a3bf0b184f16165be90946 |
| ctrl_rand_emb_s3/hHep_1000.nets.npz | 90d9413923816558f646b957d334a7f167b6a01d6ec3d11b6442f954bb45f8d0 |
| ctrl_rand_emb_s3/hHep_500.nets.npz | 34ac9dd99472f258edecfdc8094d962113d156a07f0af92a5816911e70c3232d |
| ctrl_rand_emb_s3/mDC_1000.nets.npz | fbbef15613a71667ba4c1ac0194d05737ca89905e743a4be1d8a0c343a9255c9 |
| ctrl_rand_emb_s3/mDC_500.nets.npz | a383121b369ebd060aed53d39304c3fb5ce65aec92a79e7b586c194aa268f158 |
| ctrl_rand_emb_s3/mESC_1000.nets.npz | 9f5b4791bf054d8c31b524e0adb0ba2bdd074ef4065ce674f24948b1934d08f7 |
| ctrl_rand_emb_s3/mESC_500.nets.npz | f5c99a68c4980095cfc2f307adcc13ee0475573a8d856667330095136b11ddd3 |
| ctrl_rand_emb_s3/mHSC-E_1000.nets.npz | 69ec5638217e79fdabdd296e3d658c2bb1c0d5238be9f1b10dedc0ae1fa8c6e2 |
| ctrl_rand_emb_s3/mHSC-E_500.nets.npz | 25de8532fd6a19203236a2e5b962bc549f8a355e67caa45f542a048c267e66c3 |
| ctrl_rand_emb_s3/mHSC-GM_1000.nets.npz | df6562ab811600ac62e8a735c67b035f67432b5cc0ad0b42bfa9b04cdbee7e38 |
| ctrl_rand_emb_s3/mHSC-GM_500.nets.npz | 18a34bf597bf1bac71261b1ed41e33cce05840ab4c84adf41b7f1451932f9e81 |
| ctrl_rand_emb_s3/mHSC-L_1000.nets.npz | a53eb4b11ad9fe71088e81f88cdd075e28b78f8fb513dd092bdd33a156578ce1 |
| ctrl_rand_emb_s3/mHSC-L_500.nets.npz | 72a392a2aa6196f6c6692e48bd9929a4acdeba7b8aafcd7170d439380eea8018 |
| ctrl_rand_esm_s1/hESC_1000.nets.npz | 423dbc51b9d3bf916c8263a1380b1ce51dd1104e2a6bf3d963a155e1e0ac0b84 |
| ctrl_rand_esm_s1/hESC_500.nets.npz | abfd62f27f3012a735a2f1a9e8be8d5ebb027e22f8bcee1df401f2ad1146ceb2 |
| ctrl_rand_esm_s1/hHep_1000.nets.npz | c1584ff466ff6b0697f28f35fd1b814c686d660c4c4300976ab6de4bcea4ddba |
| ctrl_rand_esm_s1/hHep_500.nets.npz | 55621401e533a36f400971d8f12abf0693220a83294839f77f25563d1dd6efd2 |
| ctrl_rand_esm_s1/mDC_1000.nets.npz | 727073797e4c61a6bd8f77f1df917ac949720721196bc9e309ec4820db9595b6 |
| ctrl_rand_esm_s1/mDC_500.nets.npz | e0176d84efd49c5835b8ebb3be5418d3219b190f4dc1a54ccffdb5fddddf42ed |
| ctrl_rand_esm_s1/mESC_1000.nets.npz | bf175877da4e5c8f25de567a152e6c983b4b5b85f9b9c70db28706f7824072d7 |
| ctrl_rand_esm_s1/mESC_500.nets.npz | 53317b2a1ae6af81733fa64fb4d402a5f98efdde4e75a519480fcb446b582551 |
| ctrl_rand_esm_s1/mHSC-E_1000.nets.npz | 886ac16ff85476ddcdc6532d0701eb798600699dc342d709a6ca108cfdef84b1 |
| ctrl_rand_esm_s1/mHSC-E_500.nets.npz | f862ae8c95a5b749c08acebd2abd58f7eb20aae7e610df99747a7a61abc3e6c9 |
| ctrl_rand_esm_s1/mHSC-GM_1000.nets.npz | cfc2d05070dc1d8418358619f67633362f3a84de037350a74020f02bcb384c01 |
| ctrl_rand_esm_s1/mHSC-GM_500.nets.npz | 660393509e9aeef87f6da47d2746edff9e77565540a0af129eab838bb3379dba |
| ctrl_rand_esm_s1/mHSC-L_1000.nets.npz | ad11a3290706973ba1e0d6ce743c12d98bfffc1b5a725a656be23de86cb3660c |
| ctrl_rand_esm_s1/mHSC-L_500.nets.npz | 12d5f3e1dceb5947c3c0e305465265eaaa315d0702f74feaff6fea548561014a |
| ctrl_rand_esm_s2/hESC_1000.nets.npz | d48190df8e1458d1374345dcfc90d1a3d4caf26ee98606582085146804fd7240 |
| ctrl_rand_esm_s2/hESC_500.nets.npz | 902a464680bc2e04071dec23cc1a10f7d2fb703be5d02551fe75568c65af72b2 |
| ctrl_rand_esm_s2/hHep_1000.nets.npz | 791d52f6298023a8f0eaf785d2f33b6c8bdfc5ea480453b8978578bf6fcd7b89 |
| ctrl_rand_esm_s2/hHep_500.nets.npz | 2d240d3c2e20c68281aa0e4c7591deb1f53d8a8ba37b5805a1df6e516fc40e99 |
| ctrl_rand_esm_s2/mDC_1000.nets.npz | c094147c21aa9b4081456197cc4601e5610f656d9d957ff963ad09c3fd62b901 |
| ctrl_rand_esm_s2/mDC_500.nets.npz | 1d0ce117187b2fba582b0acb8e8a4943b42658d1e7baa66d573cdaebb30addd8 |
| ctrl_rand_esm_s2/mESC_1000.nets.npz | d1c584103e8394f550083ff196845628b65c25a0aff45e90bb49f9dc4cf78f08 |
| ctrl_rand_esm_s2/mESC_500.nets.npz | d4cedfe45e4e1f196dbfa6c42fafebac4f587c820191560b7d96be49d6027473 |
| ctrl_rand_esm_s2/mHSC-E_1000.nets.npz | 54eebc93874a1982d3ba3d04dcfa36d1161f94ab84fec018b735fb83117b717e |
| ctrl_rand_esm_s2/mHSC-E_500.nets.npz | 00c94412d0564028566845d5a25d3d93e58a42674620a69518b8f135bc2c6631 |
| ctrl_rand_esm_s2/mHSC-GM_1000.nets.npz | 9c69cb95056385f9831a2c6b0c474af97a0421615030076c11676174726f630a |
| ctrl_rand_esm_s2/mHSC-GM_500.nets.npz | 6d440e880d94adf1ceed7511e8111c7181115c99955be1db2a887a42e6cb39fa |
| ctrl_rand_esm_s2/mHSC-L_1000.nets.npz | d6eff7cd86a2cfe60b3a36044dc9b84291fe38abe4031f61ee49ca8bd5f60bd1 |
| ctrl_rand_esm_s2/mHSC-L_500.nets.npz | e776efebe5abc09008bb4029079d0b06cdd820142aba94c44aa58541dd007814 |
| ctrl_rand_esm_s3/hESC_1000.nets.npz | 206ff4b46e49b63c4eb880dbcf01ebc97dcedbbc12601f18ee30d4de89276e43 |
| ctrl_rand_esm_s3/hESC_500.nets.npz | 74579e8aefef9da2429e747c34fda9d3f1bdd04d0b9a9ad07dba4bcd0d6e3e26 |
| ctrl_rand_esm_s3/hHep_1000.nets.npz | 70f6ff284ab34562f7a8907f3e0842e8eb2906c448eae7d6ec9f546ab2e06b09 |
| ctrl_rand_esm_s3/hHep_500.nets.npz | 104619b8059e2c5825cdc831dc9adc7e876aa78b97bc0536b1d79fd6d5505413 |
| ctrl_rand_esm_s3/mDC_1000.nets.npz | 7b37acb1bf3057c692973df6e92fd9653c372ab8cb55cc39f3be5a5e99c8725c |
| ctrl_rand_esm_s3/mDC_500.nets.npz | 3a94d406cc1fe233940990ea00d488109fa1a0d204fcb02eea2230c23ed5f98e |
| ctrl_rand_esm_s3/mESC_1000.nets.npz | 4f7aed1fb18b922a68a4c9f6fcfc2a880033fd02934134e1ea3e08de5d759068 |
| ctrl_rand_esm_s3/mESC_500.nets.npz | a331204f0bf1fed04f07cc8e2b8611f0fb7987e7083bc809eb090b4ff4f7d179 |
| ctrl_rand_esm_s3/mHSC-E_1000.nets.npz | 5d0a36f364bd40f421f054c531adf0384fab1ae97c62181f78232470b1d2af2c |
| ctrl_rand_esm_s3/mHSC-E_500.nets.npz | eaee2f90fa315d8922a2cb0b9c9edee96c28f195bf1196bc56936d1a75110f9d |
| ctrl_rand_esm_s3/mHSC-GM_1000.nets.npz | b671a963e0b311a14a30e5470fe879b71d7f8d9d1d712045d5e90b098ef19313 |
| ctrl_rand_esm_s3/mHSC-GM_500.nets.npz | b1c1abed4e916016bc7956c25e84a6284d66af0ec33659e0d114129cfbb08b16 |
| ctrl_rand_esm_s3/mHSC-L_1000.nets.npz | 2c41ae820ceed5332f36a8f9d9d7a58a4e2f9c73589717354ca70a286173f3f7 |
| ctrl_rand_esm_s3/mHSC-L_500.nets.npz | 7e005a1f05c48e1d412b4d03350511ac3e1d211f7987a0c82fe7846d58a3a0bf |
| r1_A/hESC_1000.r1.npz | 77258e10711f826de63381d205fe35982cdc23f341a65e7cc1f89124b289636d |
| r1_A/hESC_500.r1.npz | 47f6d13f38d693d86c08d4e26d4e994a4d11331ae6073b18ed6381de78591137 |
| r1_A/hHep_1000.r1.npz | 8df7e92c8ac85c306c57eeeef73257fa56774e2725cfdaed63cd73b6106281ac |
| r1_A/hHep_500.r1.npz | e31d06802c6a36206d3989184fd9bc77ee4329d8704ef39791be50aadc833682 |
| r1_A/mDC_1000.r1.npz | 5c5eec0f81c45a347e8ffb36c94c50be9a8c4dbbf010bba2d2402a26389b2d48 |
| r1_A/mDC_500.r1.npz | f9630a755ee96280e22196fda86dbabca40c690846516f2024815eeb90b0e2e2 |
| r1_A/mESC_1000.r1.npz | d1e5dda7468de0936d434807e628ec247ad819bcf165bc3ad695298c0704fb24 |
| r1_A/mESC_500.r1.npz | 84a6ffc8b0ea55c6a68d79c31ceffad086b9061053a3cbb9c304c9ffcbb94477 |
| r1_A/mHSC-E_1000.r1.npz | 4d8818714d6b00b12e2369c5a2a353ef4ab4be318d4d805f90f10902e5dc7b00 |
| r1_A/mHSC-E_500.r1.npz | 8f5f387c48171ea16941a2576d8ebefca9b438417277acddd9b75d0d8a0eb259 |
| r1_A/mHSC-GM_1000.r1.npz | b8df3f1dfb04359b89d784ee764a20a7ce99a58ad3005263d53c8b151d06ddfe |
| r1_A/mHSC-GM_500.r1.npz | d5b0322fc5787376626f4f055b8db76e543e6c99f106c120ae465e49b27f9a42 |
| r1_A/mHSC-L_1000.r1.npz | d2fcf91c92b946ca2a48069deabb74c01bba9e28026be8b5216c1066da8fe37c |
| r1_A/mHSC-L_500.r1.npz | dfa065640d6f8f2dfb318e586bc522ece9e3ebfa4b40a4b41cd3c121f9004675 |
| r1_B/hESC_1000.r1.npz | 3cedf96e30bd010dfb0a474595b60cd689503251f370a8a148836fba62f1c00b |
| r1_B/hESC_500.r1.npz | 8be71d1120d7319133be80fbb954e616089227764a31da199cdb75c016bd61c5 |
| r1_B/mHSC-E_1000.r1.npz | 491eee3ec1906bd8271c493605151d87f1fcca3ce8ac9d3535a945ad9ff0ae62 |
| r1_B/mHSC-E_500.r1.npz | b3d447d61b7c5bcd8b693971bf5b6214866fc8720df409da5cbf220934aa492b |
| r1_B/mHSC-GM_1000.r1.npz | e631a767ebf4da7a72d85912071386d3133b978dc6905ad7656391061ab3c7e9 |
| r1_B/mHSC-GM_500.r1.npz | 04a995775b6bec69f1ffd89129493f333038e0d9855f547f06d30d06c50b9fde |
| r1_B/mHSC-L_1000.r1.npz | 5684a4a2756c8725ed1887c1fc89ee651fa1654b158d07f829bf309a70a8752f |
| r1_B/mHSC-L_500.r1.npz | 8a65b9b7e6b6ffca954f28df4ed168c620a7dc4ec01dfb7b46c730dcd8ed568a |
| scgpt_A/hESC_1000.nets.npz | 45970b8d99cb683ee654c4ba2cc75809a31a5c98769c507bc098a994332e1559 |
| scgpt_A/hESC_1000_seed1.nets.npz | 2d58f4d6808b7a045fd7580be3bb8590a7b481a2fa1baa2a72f9ade2e22f2b7b |
| scgpt_A/hESC_500.nets.npz | ec91914be710e46c3d0481efaf09bd5c4a068366710e6ecdba6866614c609768 |
| scgpt_A/hHep_1000.nets.npz | eef6b4d17487b526abe6fb561c77c440a20a9d84cd1f21888e4f012a2ef68531 |
| scgpt_A/hHep_500.nets.npz | 1d666b9acaaa3fe598c183ce33b94ab8d4588d715a8ffd53ae9f2f625a29720f |
| scprint_large_A/hESC_1000.nets.npz | 37b51e2e77857abf6f2978994aeb2ea05469fe497b80d5c1a7af3d5854e14b67 |
| scprint_large_A/hESC_500.nets.npz | efd2650c4004e3f9d94aa407256aa91b835eec8b1ef77e364a22bb2b70fc80ad |
| scprint_large_A/hHep_1000.nets.npz | 9845cb03bbf0d6c72bff436bbba8f7b1c33aa330ee77d89e707d22af3b4bb553 |
| scprint_large_A/hHep_500.nets.npz | c233ad975af54c7367139d8b22b67c963c207fa6d0c1637a84e2e4ca6fb76148 |
| scprint_large_A/mDC_1000.nets.npz | e8d735534eaf89977ce4ad42baa90778b9a83d4a3c4853c82338e20bb598ab98 |
| scprint_large_A/mDC_500.nets.npz | e4dc7dc7953149c922a0b1b48f80f46838be02d57242149b5a48419f08b37d2b |
| scprint_large_A/mESC_1000.nets.npz | 1d0b8029465afa5ce8df6d791b9a79cbac043f2c5e3e8f46584eadcd64b40b12 |
| scprint_large_A/mESC_500.nets.npz | 0129ee14f8488beca1c226881c458e3376eaa7e174babc89f8893c67fb9fdb58 |
| scprint_large_A/mHSC-E_1000.nets.npz | 68d5d4620978ae605440c35f0811a5da9176850eeab5ec5b454f497e05090259 |
| scprint_large_A/mHSC-E_500.nets.npz | 8c8933059e7e1cc6a185771990fbbcbc467c9b2255af6dfdc2e26e9b53b59a9c |
| scprint_large_A/mHSC-GM_1000.nets.npz | 3c05558eaee5319b9106c14f2a3af7aab9d0b3abec1440de0761a76dd8d81ca0 |
| scprint_large_A/mHSC-GM_500.nets.npz | 39b03c0eec2c28f9c80c772eacd01c81a0d91c5ffe02d7c66ef8721dfd69a42c |
| scprint_large_A/mHSC-L_1000.nets.npz | a91b8ecb5588d5c4a6c53954ab0798c90f28a8e6d54b0013b86f3eca8a7a0b41 |
| scprint_large_A/mHSC-L_500.nets.npz | bce8b03f2223c85ad329473ff96738e3b7495b59f9788c4d745ed801c0ca33f8 |
| scprint_medium_A/hESC_1000.attn.npy | 65f5ddfbb24aaf004f5eaaa0288a8ada974e322e062e511a35170997c8364acc |
| scprint_medium_A/hESC_1000.nets.npz | c2c6f2b2ad61feea913cbe88bbdcc0b1bc9d3a315114e5866b5cbc4c8edc25e0 |
| scprint_medium_A/hESC_500.attn.npy | 4ea8d13011774f227d37c16c1844a5eb5dbadd8684ec5d3b0c45f14b5de552b5 |
| scprint_medium_A/hESC_500.nets.npz | 6f90ef4d3b2a9611f2a54ed5a1b02b2afad3ac26c152ea07a9f5b4909cfd5f88 |
| scprint_medium_A/hHep_1000.attn.npy | 377cf12784035f795b82952e4935a07914b51ee22c2c52e8b0f350776890d16f |
| scprint_medium_A/hHep_1000.nets.npz | 42f4d23b7bca0dd5f88c584075bd379ca42137270c35d7620502ce1cc7191b60 |
| scprint_medium_A/hHep_500.attn.npy | 56d67fe115ee6d14010783c5f6f30eecf89268615009354b973e800551f80b8e |
| scprint_medium_A/hHep_500.nets.npz | 5e6e4bf9f729d23edc098774585b3d4164c205fe1ded1229479599aeec30e893 |
| scprint_medium_A/mDC_1000.attn.npy | 7ad6a3e23011e4efddf02d582616e5a2053b7c85ae110fa3443504cea0d1df03 |
| scprint_medium_A/mDC_1000.nets.npz | 3b9575d31b068c2199767e51bc33c1089857b87c97ad962cebad6ec306bfcbd3 |
| scprint_medium_A/mDC_500.attn.npy | 15cb44d6ba59eba794a73dd2a07f5b171ad6f1b515d66faa6b25b6c4bad39844 |
| scprint_medium_A/mDC_500.nets.npz | 3daa10cc2d38688d3dae6e868710cd24206eeaa4a80197f7d907c370ceeb5283 |
| scprint_medium_A/mESC_1000.attn.npy | 37f3992d46209cd7bad0560f4d48d3c2425c7db29db36af78760d90e26bcc742 |
| scprint_medium_A/mESC_1000.nets.npz | 91f77f1a3f9ca38938be2f45d539ba54ad8d38325d7e9491d79495a298d3c932 |
| scprint_medium_A/mESC_500.attn.npy | 05104095872be00bf4e32e824985ae480ac2f3bf7394c3c02d2392f94e390678 |
| scprint_medium_A/mESC_500.nets.npz | eb58d00a20b83c3c4a0fd04a2039e5d18578bab36281d1363f2027df3380aebe |
| scprint_medium_A/mHSC-E_1000.attn.npy | 403dd389dca4be6ffbbb4efda01b9f567876c32bb74bfa228e993e14fab97e61 |
| scprint_medium_A/mHSC-E_1000.nets.npz | 6226f021e2237e7f805d66c92457a3d5e2d49826b1dc5eef9732de6d52b0617b |
| scprint_medium_A/mHSC-E_500.attn.npy | 05c5876e7c060de0ae3b60bc9d9c7fe6de5f0414b30e94ffbd2708f63f526bd0 |
| scprint_medium_A/mHSC-E_500.nets.npz | 3c265b5f140b900246dfa911d4615786d99ed118ad751cb28b9eb77097059c4a |
| scprint_medium_A/mHSC-GM_1000.attn.npy | 49121034646a453be9e4faa19f1c36ed3085a73281e0f92f66239813a402e2df |
| scprint_medium_A/mHSC-GM_1000.nets.npz | 83966e1975eece4311ebcdfd25cd3d52418c2864dfe8ce4b9ae38340a5572f39 |
| scprint_medium_A/mHSC-GM_500.attn.npy | 84567e180855cc2048e1ebede8e8345cbbf8ed8fc69c84d28076e2015892ec8d |
| scprint_medium_A/mHSC-GM_500.nets.npz | 0d1cce74e4c82314b4dbb31637a8f9d6c26f08c33e75507d4489c478cefd8a3e |
| scprint_medium_A/mHSC-L_1000.attn.npy | cccf8e86d7ef8853c94dbcd1095ec8c3d6c3b442a348e25a68ae632d26dc536b |
| scprint_medium_A/mHSC-L_1000.nets.npz | 95ac20ffa3e80bcf448ee15b7c8d2ffa8800945c60481fcd688291c2f9cebbd5 |
| scprint_medium_A/mHSC-L_500.attn.npy | 7f2810c0c6cc61ed9dc45ec712d5b91a74f47148db00684a06eaa407c0f8ccee |
| scprint_medium_A/mHSC-L_500.nets.npz | 5d2f57b6349ecb07bbb7691f49c1476b7fc4c313af32833929560c0683751921 |
| scprint_medium_B/hESC_1000.attn.npy | 904f36d4f73676bd5c7824b4fba49c135d0252d5032c6e475f03f4406c44978a |
| scprint_medium_B/hESC_1000.nets.npz | efa775d1751b1a0c8fdabbac940f0ac4900e41d041a50eaaa42190d265c9a8c2 |
| scprint_medium_B/hESC_500.attn.npy | 11c9c1c53b4de08bd50a49aec17d867cd8d9fe9473134f35071fae5fa4e0c30d |
| scprint_medium_B/hESC_500.nets.npz | 532c838193c29a1e3ec9b671dcda73bb2843221bed8dd0e323c99fef85754735 |
| scprint_medium_B/mHSC-E_1000.attn.npy | 067a5f65692d7a48374a574768ae19cf990728077fcd32b5dc0871bdffb95d40 |
| scprint_medium_B/mHSC-E_1000.nets.npz | 45a37065f4d9e54dc8e6cab574be5e17fc0c94a31e15fa4a730f0aa1f9985b3b |
| scprint_medium_B/mHSC-E_500.attn.npy | 0838452bf8af3bfcf9b55dae81dedf65c5e73ac21749099cbf50b2afe8eadc82 |
| scprint_medium_B/mHSC-E_500.nets.npz | f09ab53b8d05e252944d7d96966aa31db63c1eadf60201a6f1f549758e753f5c |
| scprint_medium_B/mHSC-GM_1000.attn.npy | 7ae5734f886a153f54af8e3ee2b38468671bd2733fe39cf14eeaafdc706ae54b |
| scprint_medium_B/mHSC-GM_1000.nets.npz | 43ab2b6d7a1413f7ff8b7db3de6c334172bb3c9af566feb7bd648e5149206d70 |
| scprint_medium_B/mHSC-GM_500.attn.npy | c1423b3e7b36f3ffd66398fddbe63fb245910db1a12382b07a68fc52e89d19e8 |
| scprint_medium_B/mHSC-GM_500.nets.npz | 617f910af59984a44da384f6357d818b26db1af0ac2c63984f9921589ff62456 |
| scprint_medium_B/mHSC-L_1000.attn.npy | 6651f2df91eec903183d2f6683e79533d9bdf9a663c94b5374956d6255820760 |
| scprint_medium_B/mHSC-L_1000.nets.npz | e75eeeed503d906d32abf5822e89123dc3f1f1765197168a769c3f5492d64d80 |
| scprint_medium_B/mHSC-L_500.attn.npy | 8c913613dbd714b924bcdc16476763333b83cb5c1983c0be19991bf8abe1d6fe |
| scprint_medium_B/mHSC-L_500.nets.npz | e0c4f1d01c89030d05c9e4427fd3651642c41ce4e22e5d6f4ceb6b9cdede48b7 |
| scprint_small_A/hESC_1000.nets.npz | f6d3a2b4206e3dece3c8bc7f5e8bc48403ad86d7a4f8c111c8a8b8c6464448ef |
| scprint_small_A/hESC_500.nets.npz | 0fdafafb897d7b6aa3018ec7a5300a355d92417b3aaa040d0f2994e4aee42485 |
| scprint_small_A/hHep_1000.nets.npz | 33348d5168bcba30c405ba8df60ca24dd65628c2a543faf749bb26bcfea2c6c9 |
| scprint_small_A/hHep_500.nets.npz | 0bb02b37db43b8bdd7582cc1b55ee79e7fca7f775aad9c29138c61d4497fbf98 |
| scprint_small_A/mDC_1000.nets.npz | c2b398567aa6fa8282ce50dc1274e987036bda5c8cb4d335aed3d281f6d722da |
| scprint_small_A/mDC_500.nets.npz | 97530875e78c8176d3f73c8319c79b3b0b66152ad8c023bccf8edaec16bebba9 |
| scprint_small_A/mESC_1000.nets.npz | d09f793427ba40eb1e6ec993e957d1327f6cb9b158495bd79f3b360053dac84f |
| scprint_small_A/mESC_500.nets.npz | ee69013588352917b60e2097e26066f5ac6ecb99662163339dadd8cee6af4e82 |
| scprint_small_A/mHSC-E_1000.nets.npz | b85be5441c6bf1d247041590ad23406e0b264bbfc0f5057af0d401fa87cc662d |
| scprint_small_A/mHSC-E_500.nets.npz | f1a94ed388679d4a2e850ee8e220bf903633c7c66f1bf96d1dbd3819fef4f296 |
| scprint_small_A/mHSC-GM_1000.nets.npz | 18f53d665e8412835372d4627f478f36d4f3087608c7805d4ac39572d3e645ab |
| scprint_small_A/mHSC-GM_500.nets.npz | 46edd57b1f335b946b1150d9e14e120e941e5f361450838e0394054f50280ffc |
| scprint_small_A/mHSC-L_1000.nets.npz | 8c0462ed58cf940af21860ad2896f1db0d450ca08496be467563f62418efb68a |
| scprint_small_A/mHSC-L_500.nets.npz | c5b9b4fbff2814dfcc646cde823d6e4eea0f82b607ed657c456d2455a4819720 |

## Reproduce

```
source ~/venvs/scfm_alt/bin/activate; REPO=/home/boris/Projects/ClrAttention
cd ~/scfm_models/lamin_instance               # smoke lamin instance auto-connects (scPRINT imports need it)
python -I $REPO/analysis/scfm_extract_scprint.py --cases all    --arm A --out ~/scfm_models/step3/scprint_medium_A --device 1     # 3.3
python -I $REPO/analysis/scfm_extract_scprint.py --cases counts --arm B --out ~/scfm_models/step3/scprint_medium_B --device 1     # 3.4
# 3.5: cd ~/scfm_models/lamin_r1/inst; CUDA_VISIBLE_DEVICES=1 LAMIN_SETTINGS_DIR=~/scfm_models/lamin_r1/home python -I $REPO/analysis/scfm_extract_r1.py --cases all --arm A --out ~/scfm_models/step3/r1_A   (and counts/B)
# 3.6: --control perm --seed 1001..1003 ; --control rand_esm|rand_emb --seed 1..3 ; analysis/scfm_controls.py esm2cos|qc
# 3.7: ~/venvs/scfm/bin/python -I $REPO/analysis/scfm_extract_scgpt.py --out ~/scfm_models/step3/scgpt_A --device 1 ; same with --cases hESC_1000 --binning-seed 1 --tag seed1
# 3.9: --ckpt-name small-v1|large-v1 --store-attn no
python -I $REPO/analysis/scfm_step3_summary.py --out $REPO/results/scfm_step3
```
The exact job files are `jobs_pc/done/1001 ... 1010_scfm_*.sh` (job logs next to them).
