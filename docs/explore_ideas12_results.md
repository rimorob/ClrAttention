# Exploratory runs of ideation-memo ideas 1 and 2 (K562 in-sample; M3D lab-split stage)

*12 November 2026. EXPLORATORY, IN-SAMPLE, NOT PREREGISTERED. K562 truths are spent for confirmation; no p-values are computed; no variant is selected for the paper; nothing here goes into a paper or pre-read. The plan (`docs/plans/explore_ideas12_plan.md`, commit 9d905ec) was committed before any experiment compute; deviations are listed in section 7. Source memo: `docs/plans/fable_ideation.md` section 2, ideas 1 and 2. Code: `analysis/explore_ideas12/`. Tables: below and `docs/explore_ideas12_tables/` (CSV). Enrichment = pooled AUPR / base rate. Intervals are d +/- 1.96 SD of the paired delete-half gene jackknife (100 replicates, the `robust_z` stream). Raw numbers; no winners marked.*

## 1. Method as run

- **Data, truths, scoring:** as the plan section 0 (K562 essential half A, 5,671 genes x 707 profiles; CORUM all/noRP/RP/small/large, DoRothEA A-B, KD response, CollecTRI; raw and two hub-stratified scopes, stratifiers pi of the V0 K10 walk and MI row mean; contrasts versus CLR and versus K10 attention at t = 9).
- **Idea 1:** attention graph W = (A + A')/2 on each of 12 half-datasets (6 `robust_z` depth-stage splits of the 707 profiles; 353/354 profiles each; `halves_info.csv`) and on the full data; Leiden (modularity with resolution, `leidenalg` 0.12.0, `python-igraph` 1.0.0 from the existing `~/venvs/scfm`, nothing installed) on gamma = 2^-6..2^6 (13), seeds 1-5, `n_iterations = -1`: 780 observed + 780 rewired-null + 65 full-data = 1,625 partitions. Rewired null: degree-preserving edge swaps (10 x edges), weights permuted. Weights: R (excess of nVI-null over nVI-observed, normalised), U (1/13), G (proportional to gamma). Consensus frequency over 60 half-runs per gamma (5 for FULL).
- **Idea 2:** detector rho from the top-20 positive-CLR partner overlap between K562 essential half A and a second compendium (hypergeometric tail p, rho = min(1, ln p / ln 0.05)); gated operator B_ij = A_ij rho_i rho_j, B_ii = 1 - sum, lazy walk, t = 9. Second compendia: gwps half A (1,553 profiles; universe 5,484 genes), RPE1 half A (1,023 profiles; universe 4,976). Controls: hard gate, half-split detector, five permuted-rho operators.
- **M3D stage** (run because CONS-R passed the idea-1 rule): same construction on M3D (4,297 genes, 466 experiments), 6 lab splits, RegulonDB SC and all-evidence regulons.

## 2. Equivalence gate

| dataset | reference | n_checks | max_abs_diff_AUPR | max_abs_diff_jackknife_replicates | all_pass_1e-9 |
|---|---|---|---|---|---|
| K562 | robust_z V0 cache + baselines metrics.csv | 50 | 3.9e-16 | 0 | True |
| M3D | robust_z m3d V0 cache | 12 | 0 | 0 | True |

K562: 50 checks (7 `robust_z` truths x 3 scopes x CLR/K10 against the stored V0 caches, point estimate and all 100 replicates; plus CLR/attention AUPR for DoRothEA A-B, CORUM, KD response and CollecTRI against `results/baselines_k562ess_truthgenes_genie3/metrics.csv`). Maximum absolute AUPR difference 3.9e-16 (tolerance 1e-9); jackknife replicates identical. (In the final gate run the CLR and K10 scorings were read from the cache written by the first scoring pass of the same code; MI, CLR, operator and walk were recomputed and all comparisons were executed. See section 7.) The M3D gate has 12 checks (2 truths x 3 scopes x CLR/K10 against `results/robust_z/m3d/cache`).

## 3. Jobs and peak memory

| job | what | declared GB | measured peak GB | ended (12 Nov) |
|---|---|---|---|---|
| 5000_ei12_gate | K562 equivalence gate, caches, CLR/K10 scoring | 20 | 3.7 (final run); 8.9 and 14.3 (two earlier runs that stopped on code errors) | 10:14 |
| 5001_ei12_i2diag | idea-2 detectors and diagnostic | 16 | 3.1 (2.3 in the earlier run that stopped on a code error) | 10:21 |
| 5002_ei12_halves | 12 half operators, graphs, top-20 lists | 14 | 2.5 | 10:21 |
| 5003_ei12_leiden | 1,625 Leiden partitions, nVI, null, counts | 6 | 2.0 | 11:31 (70 min) |
| 5004_ei12_score1 | idea-1 scoring (7 readouts + 13 single-gamma) | 20 | 13.6 | 11:45 |
| 5005_ei12_score2 | idea-2 scoring (9 gated operators) | 20 | 15.0 | 12:06 |
| 5010_ei12_m3d_prep | M3D gate, references, lab splits, half graphs | 14 | 5.6 | 14:57 |
| 5011_ei12_m3d_leiden | M3D Leiden sweep | 5 | 1.5 | 15:34 (36 min) |
| 5012_ei12_m3d_score | M3D scoring | 14 | 4.4 | 15:37 |

Peaks from `jobs_pc/runner.log` END lines. No OOM kill. Disk: `results/explore_ideas12/` 0.52 GB kept (caches 0.36 GB for K562, graphs, memberships, tables; the 1.9 GB of uint8 count matrices were deleted after scoring and are regenerable from `memberships.npz` by `leiden_sweep.py`); small CSVs and tables copied to `docs/explore_ideas12_tables/` (0.3 MB). A row for `docs/disk_budget.md` was drafted but that file has another worker's uncommitted edits, so it is not committed by this work (row text in section 7).

## 4. Idea 1 (K562): literal pass/fail

Rule (plan, copied from the memo): "Pass if CORUM raw >= 6.4x (CLR) with DoRothEA within -0.005 of K10; fail if CORUM < K10 (3.51x) or DoRothEA -0.01." CORUM = CORUM_all pooled raw enrichment; dDoRothEA = DoRothEA A-B pooled raw AUPR minus K10's (0.20032).

| readout | CORUM_all_enrichment | K10_CORUM_all_enrichment | dDoRothEA_AB_aupr_vs_K10 | verdict |
|---|---|---|---|---|
| CONS-R | 6.6665 | 3.5106 | -0.0011 | PASS |
| CONS-U | 5.9915 | 3.5106 | 0.0035 | between |
| CONS-G | 4.4912 | 3.5106 | 0.0022 | between |
| FULL-R | 2.6012 | 3.5106 | -0.0097 | FAIL |
| FULL-U | 2.5589 | 3.5106 | -0.0054 | FAIL |
| NULL-U | 1.05 | 3.5106 | -0.0184 | FAIL |
| NULL-R | 1.0024 | 3.5106 | -0.0173 | FAIL |

- **CONS-R (primary): PASS** (CORUM_all 6.67x; dDoRothEA -0.0011 [-0.0069, +0.0046]).
- CONS-U (uniform weights): 5.99x, dDoRothEA +0.0035: under the rule it is neither PASS nor FAIL ("between"). The plan flags "disguise risk realised" if the uniform form *fails* the rule; it does not fail. CORUM_all of the uniform form is 5.99x against K10 3.51x.
- FULL-R / FULL-U (full-data graph, 5 runs per gamma) and the rewired-operator consensus (NULL-U, NULL-R): FAIL by the rule.
- Memo falsifier (i), CONS-R minus CLR on CORUM_all (raw AUPR): +0.0037 [-0.0065, +0.0138]; on CORUM_RP: -0.0145 [-0.0309, +0.0019] (both intervals include 0).

### Resolution grid, reproducibility and rewired null (K562)

| gamma | nvi_obs | nvi_null | excess | w_R | w_U | w_G | n_comm_obs | largest_obs | n_comm_null | largest_null | n_comm_full | largest_full | singleton_frac_obs |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.01562 | 0 | 0 | 0 | 0 | 0.07692 | 0.0001221 | 1 | 5671 | 1 | 5671 | 1 | 5671 | 0 |
| 0.03125 | 0 | 0 | 0 | 0 | 0.07692 | 0.0002442 | 1 | 5671 | 1 | 5671 | 1 | 5671 | 0 |
| 0.0625 | 0 | 0 | 0 | 0 | 0.07692 | 0.0004883 | 1 | 5671 | 1 | 5671 | 1 | 5671 | 0 |
| 0.125 | 0.0003754 | 0 | 0 | 0 | 0.07692 | 0.0009767 | 1.1 | 5670 | 1 | 5671 | 1 | 5671 | 0 |
| 0.25 | 0.1102 | 0 | 0 | 0 | 0.07692 | 0.001953 | 2.45 | 4097 | 1 | 5671 | 2.4 | 3328 | 0 |
| 0.5 | 0.3161 | 0.03191 | 0 | 0 | 0.07692 | 0.003907 | 8.05 | 1574 | 3.567 | 5514 | 7.8 | 1538 | 0 |
| 1 | 0.5054 | 0.8442 | 0.3387 | 0.4408 | 0.07692 | 0.007813 | 21.67 | 734.1 | 52.22 | 192.6 | 16.6 | 796.4 | 0 |
| 2 | 0.628 | 0.8417 | 0.2138 | 0.2782 | 0.07692 | 0.01563 | 52.35 | 346.8 | 115.9 | 85.27 | 43 | 364 | 0 |
| 4 | 0.6533 | 0.7602 | 0.1069 | 0.1391 | 0.07692 | 0.03125 | 111 | 166.7 | 197.7 | 48.95 | 93.2 | 181.2 | 0 |
| 8 | 0.6186 | 0.6714 | 0.05279 | 0.06869 | 0.07692 | 0.06251 | 208.3 | 102.6 | 308.1 | 31.33 | 179.4 | 109 | 0 |
| 16 | 0.5572 | 0.5858 | 0.02853 | 0.03712 | 0.07692 | 0.125 | 360.3 | 62.87 | 456.9 | 21.08 | 324 | 63.4 | 0 |
| 32 | 0.4843 | 0.5017 | 0.01747 | 0.02273 | 0.07692 | 0.25 | 583 | 37.42 | 667.3 | 15.2 | 536.8 | 39.4 | 0 |
| 64 | 0.4103 | 0.4206 | 0.01034 | 0.01346 | 0.07692 | 0.5001 | 883.5 | 24.18 | 956.2 | 11.08 | 842.6 | 24 | 1.469e-05 |

nVI = variation of information between the half-1 and half-2 partitions (observed) or between the two independently rewired halves (null), divided by ln G; excess = max(0, nVI_null - nVI_obs). w_R is nonzero for gamma = 1 to 64 only.

### Single-gamma consensus, pooled raw enrichment (K562, point estimates only, observed halves)

| gamma | CORUM_RP | CORUM_all | CORUM_large | CORUM_noRP | CORUM_small | CollecTRI | DoRothEA_AB | KD_response |
|---|---|---|---|---|---|---|---|---|
| 0.01562 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |
| 0.03125 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |
| 0.0625 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |
| 0.125 | 1.0 | 1.0 | 1.01 | 1.0 | 1.0 | 1.0 | 1.0 | 0.97 |
| 0.25 | 1.18 | 1.59 | 1.71 | 1.16 | 1.13 | 1.07 | 1.09 | 1.18 |
| 0.5 | 1.2 | 4.65 | 2.27 | 1.73 | 1.55 | 1.0 | 1.07 | 1.55 |
| 1 | 1.18 | 6.62 | 2.48 | 2.16 | 1.92 | 1.03 | 1.08 | 1.65 |
| 2 | 1.18 | 6.69 | 2.45 | 2.27 | 2.1 | 1.04 | 1.08 | 1.67 |
| 4 | 1.17 | 6.56 | 2.36 | 2.17 | 2.17 | 1.03 | 1.06 | 1.59 |
| 8 | 1.15 | 5.32 | 2.17 | 2.2 | 2.08 | 1.03 | 1.04 | 1.45 |
| 16 | 1.12 | 4.09 | 1.87 | 2.13 | 2.05 | 1.02 | 1.02 | 1.31 |
| 32 | 1.07 | 2.97 | 1.57 | 1.71 | 1.84 | 1.01 | 1.01 | 1.21 |
| 64 | 1.04 | 2.1 | 1.35 | 1.5 | 1.72 | 1.01 | 1.01 | 1.13 |

## 5. Idea 2 (K562): diagnostic and literal pass/fail

Diagnostic (label-free, run before any truth). Plan rule: fraction of genes with rho < 0.5 on gwps <= 0.20 proceed; >= 0.30 stop; between: proceed, flagged. (Memo reference value: 10.1% by the half-split detector, a top-50 Poisson rule; the memo's "40%" is the stop case.)

| detector | n_profiles_Y | n_universe | n_outside | frac_rho_lt_0.5 | frac_rho_lt_0.5_all_genes | frac_overlap0 | frac_p_le_0.05 | median_overlap | mean_overlap | chance_overlap |
|---|---|---|---|---|---|---|---|---|---|---|
| gwps | 1553.0 | 5484 | 187 | 0.2259 | 0.2185 | 0.2259 | 0.6431 | 3 | 3.9705 | 0.073 |
| rpe1 | 1023.0 | 4976 | 695 | 0.5629 | 0.4939 | 0.5629 | 0.2484 | 0 | 1.3014 | 0.0804 |
| half_split | nan | 5671 | 0 | 0.3204 | 0.3509 | 0.3204 | 0.5209 | 2 | 3.1257 | 0.0705 |

gwps fraction rho < 0.5 = 0.2259: **proceed, flagged "not near 10.1%"**. (The rho < 0.5 set coincides with overlap 0 in all four detectors by construction.) The RPE1 fraction (0.563) is above the stop level, but the plan applies the rule to gwps only; RPE1 and the half-split detector are scored as secondary readouts.

Rho summaries of the scored gates:

| readout | frac_rho_lt_0.5 | mean_rho | n_rho_eq_1 | n_rho_eq_0 |
|---|---|---|---|---|
| G-SOFT | 0.2185 | 0.7669 | 3714 | 1239 |
| G-HARD | 0.3451 | 0.6549 | 3714 | 1957 |
| HALF-SOFT | 0.3509 | 0.663 | 2096 | 749 |
| R-SOFT | 0.4939 | 0.4818 | 1931 | 2801 |
| PERM-1 | 0.2185 | 0.7669 | 3714 | 1239 |
| PERM-2 | 0.2185 | 0.7669 | 3714 | 1239 |
| PERM-3 | 0.2185 | 0.7669 | 3714 | 1239 |
| PERM-4 | 0.2185 | 0.7669 | 3714 | 1239 |
| PERM-5 | 0.2185 | 0.7669 | 3714 | 1239 |

Pass rule (copied from the memo): "Pass: small complexes (3-6) >= +0.005 and DoRothEA >= -0.005", operationalised on G-SOFT as pooled raw AUPR minus K10's.

| readout | dCORUM_small_aupr_vs_K10 | dDoRothEA_AB_aupr_vs_K10 | verdict | diagnostic_flag |
|---|---|---|---|---|
| G-SOFT | -0.00018 | -0.00185 | FAIL | proceed_flag |

- **G-SOFT: FAIL** (CORUM_small -0.00018, needs >= +0.005; DoRothEA A-B -0.0019, satisfies >= -0.005). The memo's +0.005 is about x1.75 of K10's CORUM_small AUPR (0.00663); the largest CORUM_small difference among all nine gated readouts is +0.0000 (G-HARD).
- Full tables (all readouts, three scopes, versus K10 and versus CLR) are in section 9.

## 6. Conditional M3D lab-split stage (idea 1; no pass rule in the memo; descriptive)

Run because CONS-R passed. Reference AUPR: CLR 0.1359 / 0.1648, K10 0.1579 / 0.1733 on SC / all (raw). Weights and single-gamma tables:

| gamma | nvi_obs | nvi_null | excess | w_R | w_U | w_G | n_comm_obs | largest_obs | n_comm_null | largest_null | n_comm_full | largest_full | singleton_frac_obs |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.01562 | 0 | 0 | 0 | 0 | 0.07692 | 0.0001221 | 1 | 4297 | 1 | 4297 | 1.4 | 4269 | 0 |
| 0.03125 | 0.006108 | 0 | 0 | 0 | 0.07692 | 0.0002442 | 1.383 | 4271 | 1 | 4297 | 1.6 | 4255 | 0 |
| 0.0625 | 0.005763 | 0 | 0 | 0 | 0.07692 | 0.0004883 | 1.867 | 4238 | 1 | 4297 | 2.2 | 4223 | 0 |
| 0.125 | 0.03598 | 0 | 0 | 0 | 0.07692 | 0.0009767 | 2.583 | 4022 | 1 | 4297 | 2.6 | 3880 | 0 |
| 0.25 | 0.2232 | 0 | 0 | 0 | 0.07692 | 0.001953 | 6.85 | 1622 | 1 | 4297 | 7.6 | 1547 | 0 |
| 0.5 | 0.337 | 0.294 | 0 | 0 | 0.07692 | 0.003907 | 13.18 | 935.8 | 4.7 | 1951 | 14.2 | 891.8 | 0 |
| 1 | 0.4414 | 0.851 | 0.4096 | 0.2587 | 0.07692 | 0.007813 | 25.98 | 451.8 | 52.03 | 145.1 | 26.4 | 475.8 | 0 |
| 2 | 0.4959 | 0.8299 | 0.334 | 0.211 | 0.07692 | 0.01563 | 50.43 | 242.5 | 107.1 | 68.72 | 50.6 | 235 | 0 |
| 4 | 0.5051 | 0.7456 | 0.2405 | 0.1519 | 0.07692 | 0.03125 | 98.28 | 150.5 | 178.9 | 41.38 | 98.4 | 149.8 | 0 |
| 8 | 0.4702 | 0.6576 | 0.1874 | 0.1184 | 0.07692 | 0.06251 | 180.7 | 86.93 | 272.9 | 27.32 | 188.6 | 85.8 | 0 |
| 16 | 0.4121 | 0.5714 | 0.1593 | 0.1006 | 0.07692 | 0.125 | 311.6 | 43.6 | 401.6 | 18.77 | 325.2 | 45 | 0 |
| 32 | 0.3483 | 0.4856 | 0.1373 | 0.08675 | 0.07692 | 0.25 | 493.5 | 32.17 | 583 | 13.33 | 503.6 | 33.8 | 0 |
| 64 | 0.2877 | 0.4027 | 0.115 | 0.07262 | 0.07692 | 0.5001 | 726.6 | 22.32 | 832.2 | 10.03 | 735.2 | 22 | 0 |

Single-gamma pooled raw enrichment (M3D):

| gamma | SC | all |
|---|---|---|
| 0.01562 | 1.0 | 1.0 |
| 0.03125 | 1.04 | 1.02 |
| 0.0625 | 1.03 | 0.99 |
| 0.125 | 1.14 | 1.06 |
| 0.25 | 1.24 | 1.14 |
| 0.5 | 1.44 | 1.29 |
| 1 | 1.53 | 1.36 |
| 2 | 1.52 | 1.38 |
| 4 | 1.44 | 1.35 |
| 8 | 1.35 | 1.27 |
| 16 | 1.23 | 1.18 |
| 32 | 1.16 | 1.12 |
| 64 | 1.11 | 1.07 |

Full tables (raw and both hub-stratified scopes, versus K10 and versus CLR) are in section 10.

## 7. Deviations and notes

1. **Code errors, reruns (no result affected).** The first gate run stopped (rc 1) at scoring because `H_att` had been removed before the lazy argument was forced; the second gate run completed both scorings and stopped (rc 1) in the comparison code (`names()` of a vector element); the third run reused the cached scorings and passed. The first diagnostic run stopped (rc 1) on the same `names()` mistake and was rerun. Fixes were to code only (`ideas12_lib.R`, `stage0_gate.R`, `stage_i2diag.R`); no readout definition changed. The final-run gate therefore compares cached CLR/K10 scorings (same code, same data) rather than a fresh scoring; MI, CLR, operator and walk were recomputed.
2. Idea-2 universe on gwps is 5,484 genes (plan text: 5,494 from raw gene-name overlap; 5,484 after the missing-value and zero-range filters of `perturbseq_bulk.R`).
3. The five PERM readouts are five random permutations of the same G-SOFT rho vector, so their rho summaries are identical by construction.
4. NULL-R uses the observed R weights on the null counts, as planned.
5. The plan's disguise flag is worded "CONS-U failing the same rule"; CONS-U is "between" (not PASS, not FAIL); both readings are stated in section 4.
6. The half-split detector's rho < 0.5 fraction here (0.32, top-20 hypergeometric rule on 353/354-profile halves) is a different definition from the memo's 10.1% (top-50, Poisson, repro_n < k0).
7. No macro (per-regulon) scoring and no hub-drift diagnostic were run (not in the plan).
8. Jobs 5000 and 5001 were submitted twice under the same name after code fixes.
9. Draft row for `docs/disk_budget.md`: `results/explore_ideas12/` | 0.52 GB kept (rds caches 0.36 GB, graphs 0.14 GB, memberships); 1.9 GB of count matrices deleted after scoring | exploratory ideas 1-2 | Yes (regenerable from the committed code; `memberships.npz` regenerate counts).

## 8. Reading (limited to the tables)

- K562: CONS-R pooled raw enrichment on CORUM_all is 6.67 against CLR 6.39 and K10 3.51; the difference to CLR is +0.0037 AUPR [-0.0065, +0.0138], to K10 +0.0418 [+0.0204, +0.0632]. CONS-U 5.99, CONS-G 4.49, FULL-R 2.60, FULL-U 2.56, NULL-U 1.05, NULL-R 1.00.
- CONS-R DoRothEA A-B AUPR 0.1992 (K10 0.2003, CLR 0.1940); CollecTRI -0.0081 [-0.0131, -0.0030] versus K10 (raw), +0.0006 [-0.0015, +0.0028] versus CLR; KD response -0.0519 versus K10 and -0.0558 versus CLR (intervals exclude 0), in every consensus readout (enrichment 1.59-1.72 against K10 2.38 and CLR 2.42).
- Hub-stratified CORUM_all, CONS-R: enrichment 3.67 (pi stratifier; K10 2.08, CLR 3.40) and 5.75 (MI stratifier; K10 2.97, CLR 5.03).
- CORUM_small: CONS-R +0.0014 [-0.0021, +0.0050] over K10 (intervals include 0 in every consensus readout); the rewired-operator consensus is at or below base rate on every truth.
- The reproducibility weights are zero for gamma <= 0.5 and largest at gamma = 1 (0.44) and 2 (0.28); single-gamma CORUM_all enrichment is 6.62 at gamma = 1, 6.69 at gamma = 2, 6.56 at gamma = 4 and falls to 2.1 at gamma = 64.
- Idea 2: G-SOFT CORUM_all 3.66x (K10 3.51, CLR 6.39), CORUM_small 1.96 (K10 2.01), DoRothEA A-B 1.08 (K10 1.09). The five permuted-rho operators are below K10 on CORUM_all (2.00-2.36x). G-HARD, HALF-SOFT and R-SOFT: CORUM_all 3.80, 3.78, 4.05; KD response 2.66, 2.67, 2.10.
- M3D (raw, versus K10): CONS-R -0.0051 [-0.0099, -0.0003] on SC and -0.0014 [-0.0052, +0.0023] on all; versus CLR +0.0170 [+0.0112, +0.0227] and +0.0071 [+0.0028, +0.0114]. Hub-stratified (pi stratifier), CONS-R versus K10 +0.0078 [+0.0035, +0.0121] (SC), +0.0114 [+0.0077, +0.0152] (all). Hub-stratified (MI stratifier) CONS-R versus K10 -0.0014 [-0.0063, +0.0035] (SC) and +0.0010 [-0.0027, +0.0047] (all). The rewired-operator consensus is at about base rate (SC 1.01x, all 1.01x).
- The memo's expected effects were CORUM raw +0.02 to +0.05 over K10 (observed +0.042 for CONS-R), DoRothEA -0.003 to 0 (observed -0.0011), large regulons -0.005 to -0.01 (CollecTRI -0.008); the table gives no entry for "large loose regulons" other than CollecTRI and DoRothEA.

## 9. Full K562 tables

### Idea 1 (CONS = halves, FULL = full-data graph, NULL = rewired operators)

#### Pooled enrichment (AUPR / base rate), raw

| truth | CLR | K10 | CONS-R | CONS-U | CONS-G | FULL-R | FULL-U | NULL-U | NULL-R |
|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB (base 0.1830) | 1.06 | 1.09 | 1.09 | 1.11 | 1.11 | 1.04 | 1.07 | 0.99 | 1.00 |
| CORUM_all (base 0.0133) | 6.39 | 3.51 | 6.67 | 5.99 | 4.49 | 2.60 | 2.56 | 1.05 | 1.00 |
| KD_response (base 0.0787) | 2.42 | 2.38 | 1.72 | 1.71 | 1.59 | 1.32 | 1.41 | 1.09 | 1.00 |
| CORUM_noRP (base 0.0119) | 2.18 | 2.05 | 2.47 | 2.44 | 2.38 | 1.59 | 1.65 | 1.02 | 1.00 |
| CORUM_RP (base 0.7502) | 1.21 | 1.18 | 1.19 | 1.21 | 1.20 | 1.14 | 1.16 | 0.98 | 1.00 |
| CORUM_small (base 0.0033) | 2.16 | 2.01 | 2.45 | 2.48 | 2.52 | 1.53 | 1.61 | 1.00 | 0.98 |
| CORUM_large (base 0.1263) | 2.29 | 1.93 | 2.55 | 2.49 | 2.39 | 1.92 | 1.92 | 1.00 | 1.00 |
| CollecTRI (base 0.1278) | 1.03 | 1.10 | 1.03 | 1.05 | 1.04 | 1.02 | 1.02 | 0.99 | 1.00 |

#### Pooled AUPR, raw

| truth | CLR | K10 | CONS-R | CONS-U | CONS-G | FULL-R | FULL-U | NULL-U | NULL-R |
|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB | 0.19395 | 0.20032 | 0.19919 | 0.20385 | 0.20255 | 0.19060 | 0.19496 | 0.18197 | 0.18302 |
| CORUM_all | 0.08471 | 0.04654 | 0.08838 | 0.07943 | 0.05954 | 0.03448 | 0.03392 | 0.01392 | 0.01329 |
| KD_response | 0.19084 | 0.18698 | 0.13508 | 0.13425 | 0.12478 | 0.10404 | 0.11108 | 0.08563 | 0.07883 |
| CORUM_noRP | 0.02590 | 0.02432 | 0.02930 | 0.02899 | 0.02819 | 0.01885 | 0.01951 | 0.01212 | 0.01187 |
| CORUM_RP | 0.90638 | 0.88797 | 0.89188 | 0.90631 | 0.90030 | 0.85259 | 0.87050 | 0.73610 | 0.74950 |
| CORUM_small | 0.00712 | 0.00663 | 0.00806 | 0.00818 | 0.00830 | 0.00505 | 0.00529 | 0.00329 | 0.00322 |
| CORUM_large | 0.28925 | 0.24321 | 0.32145 | 0.31448 | 0.30239 | 0.24187 | 0.24294 | 0.12653 | 0.12681 |
| CollecTRI | 0.13150 | 0.14021 | 0.13214 | 0.13440 | 0.13285 | 0.13029 | 0.13037 | 0.12662 | 0.12781 |

#### Pooled enrichment (AUPR / base rate), hub-stratified (pi of V0 attention)

| truth | CLR | K10 | CONS-R | CONS-U | CONS-G | FULL-R | FULL-U | NULL-U | NULL-R |
|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB (base 0.1830) | 1.06 | 1.09 | 1.09 | 1.11 | 1.11 | 1.05 | 1.07 | 0.99 | 1.00 |
| CORUM_all (base 0.0133) | 3.40 | 2.08 | 3.67 | 3.46 | 3.25 | 2.45 | 2.37 | 1.03 | 1.03 |
| KD_response (base 0.0787) | 1.91 | 1.60 | 1.63 | 1.59 | 1.60 | 1.79 | 1.49 | 1.05 | 1.07 |
| CORUM_noRP (base 0.0119) | 1.71 | 1.58 | 1.93 | 1.92 | 2.03 | 1.64 | 1.64 | 1.01 | 1.01 |
| CORUM_RP (base 0.7502) | 1.21 | 1.18 | 1.19 | 1.20 | 1.20 | 1.16 | 1.16 | 0.98 | 1.01 |
| CORUM_small (base 0.0033) | 1.78 | 1.61 | 2.02 | 2.01 | 2.07 | 1.54 | 1.58 | 1.00 | 0.97 |
| CORUM_large (base 0.1263) | 2.18 | 1.94 | 2.29 | 2.24 | 2.24 | 1.87 | 1.88 | 1.01 | 1.00 |
| CollecTRI (base 0.1278) | 1.01 | 1.02 | 1.05 | 1.06 | 1.06 | 1.09 | 1.03 | 0.98 | 1.01 |

#### Pooled AUPR, hub-stratified (pi of V0 attention)

| truth | CLR | K10 | CONS-R | CONS-U | CONS-G | FULL-R | FULL-U | NULL-U | NULL-R |
|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB | 0.19385 | 0.19952 | 0.19933 | 0.20279 | 0.20246 | 0.19291 | 0.19649 | 0.18168 | 0.18319 |
| CORUM_all | 0.04508 | 0.02752 | 0.04871 | 0.04583 | 0.04303 | 0.03250 | 0.03139 | 0.01370 | 0.01363 |
| KD_response | 0.15037 | 0.12590 | 0.12810 | 0.12551 | 0.12576 | 0.14107 | 0.11705 | 0.08253 | 0.08443 |
| CORUM_noRP | 0.02024 | 0.01871 | 0.02293 | 0.02274 | 0.02405 | 0.01949 | 0.01940 | 0.01202 | 0.01200 |
| CORUM_RP | 0.90395 | 0.88445 | 0.89350 | 0.89933 | 0.89931 | 0.87389 | 0.87385 | 0.73388 | 0.75989 |
| CORUM_small | 0.00585 | 0.00529 | 0.00666 | 0.00662 | 0.00681 | 0.00508 | 0.00522 | 0.00329 | 0.00321 |
| CORUM_large | 0.27577 | 0.24444 | 0.28961 | 0.28311 | 0.28306 | 0.23556 | 0.23762 | 0.12692 | 0.12676 |
| CollecTRI | 0.12845 | 0.12990 | 0.13355 | 0.13584 | 0.13518 | 0.13879 | 0.13177 | 0.12528 | 0.12917 |

#### Pooled enrichment (AUPR / base rate), hub-stratified (MI row mean)

| truth | CLR | K10 | CONS-R | CONS-U | CONS-G | FULL-R | FULL-U | NULL-U | NULL-R |
|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB (base 0.1830) | 1.10 | 1.11 | 1.11 | 1.12 | 1.11 | 1.07 | 1.07 | 1.00 | 1.00 |
| CORUM_all (base 0.0133) | 5.03 | 2.97 | 5.75 | 5.10 | 3.97 | 2.48 | 2.40 | 1.04 | 1.02 |
| KD_response (base 0.0787) | 2.19 | 2.18 | 1.69 | 1.69 | 1.56 | 1.29 | 1.41 | 1.06 | 1.03 |
| CORUM_noRP (base 0.0119) | 2.30 | 2.18 | 2.58 | 2.51 | 2.39 | 1.69 | 1.70 | 1.02 | 1.00 |
| CORUM_RP (base 0.7502) | 1.21 | 1.18 | 1.20 | 1.20 | 1.20 | 1.15 | 1.16 | 0.98 | 1.01 |
| CORUM_small (base 0.0033) | 2.32 | 2.08 | 2.58 | 2.65 | 2.66 | 1.64 | 1.67 | 1.00 | 0.98 |
| CORUM_large (base 0.1263) | 2.32 | 1.87 | 2.53 | 2.43 | 2.34 | 2.01 | 1.92 | 1.00 | 1.01 |
| CollecTRI (base 0.1278) | 1.06 | 1.12 | 1.04 | 1.05 | 1.04 | 1.02 | 1.02 | 1.00 | 1.00 |

#### Pooled AUPR, hub-stratified (MI row mean)

| truth | CLR | K10 | CONS-R | CONS-U | CONS-G | FULL-R | FULL-U | NULL-U | NULL-R |
|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB | 0.20052 | 0.20346 | 0.20228 | 0.20501 | 0.20292 | 0.19499 | 0.19653 | 0.18274 | 0.18234 |
| CORUM_all | 0.06673 | 0.03934 | 0.07617 | 0.06767 | 0.05262 | 0.03287 | 0.03182 | 0.01375 | 0.01349 |
| KD_response | 0.17224 | 0.17174 | 0.13314 | 0.13324 | 0.12314 | 0.10131 | 0.11102 | 0.08374 | 0.08094 |
| CORUM_noRP | 0.02724 | 0.02583 | 0.03058 | 0.02973 | 0.02837 | 0.02002 | 0.02012 | 0.01211 | 0.01188 |
| CORUM_RP | 0.90554 | 0.88247 | 0.89678 | 0.90322 | 0.89825 | 0.86262 | 0.87019 | 0.73352 | 0.75712 |
| CORUM_small | 0.00765 | 0.00686 | 0.00851 | 0.00873 | 0.00875 | 0.00541 | 0.00551 | 0.00329 | 0.00321 |
| CORUM_large | 0.29338 | 0.23637 | 0.32008 | 0.30700 | 0.29491 | 0.25440 | 0.24288 | 0.12674 | 0.12699 |
| CollecTRI | 0.13543 | 0.14287 | 0.13297 | 0.13425 | 0.13242 | 0.13072 | 0.13071 | 0.12717 | 0.12724 |

#### Difference in AUPR vs K10 attention (t = 9) [95% paired delete-half gene jackknife], raw

| truth | CONS-R | CONS-U | CONS-G | FULL-R | FULL-U | NULL-U | NULL-R |
|---|---|---|---|---|---|---|---|
| DoRothEA_AB | -0.0011 [-0.0069, +0.0046] | +0.0035 [-0.0020, +0.0091] | +0.0022 [-0.0033, +0.0078] | -0.0097 [-0.0155, -0.0039] | -0.0054 [-0.0111, +0.0003] | -0.0184 [-0.0256, -0.0112] | -0.0173 [-0.0238, -0.0108] |
| CORUM_all | +0.0418 [+0.0204, +0.0632] | +0.0329 [+0.0151, +0.0507] | +0.0130 [+0.0036, +0.0224] | -0.0121 [-0.0209, -0.0033] | -0.0126 [-0.0217, -0.0035] | -0.0326 [-0.0472, -0.0180] | -0.0333 [-0.0480, -0.0185] |
| KD_response | -0.0519 [-0.0605, -0.0433] | -0.0527 [-0.0614, -0.0441] | -0.0622 [-0.0705, -0.0539] | -0.0829 [-0.0915, -0.0744] | -0.0759 [-0.0844, -0.0673] | -0.1013 [-0.1105, -0.0922] | -0.1081 [-0.1179, -0.0984] |
| CORUM_noRP | +0.0050 [+0.0002, +0.0097] | +0.0047 [-0.0000, +0.0094] | +0.0039 [-0.0003, +0.0081] | -0.0055 [-0.0091, -0.0019] | -0.0048 [-0.0084, -0.0012] | -0.0122 [-0.0170, -0.0075] | -0.0125 [-0.0172, -0.0077] |
| CORUM_RP | +0.0039 [-0.0161, +0.0239] | +0.0183 [-0.0056, +0.0423] | +0.0123 [-0.0152, +0.0398] | -0.0354 [-0.0643, -0.0065] | -0.0175 [-0.0498, +0.0149] | -0.1519 [-0.2050, -0.0987] | -0.1385 [-0.1882, -0.0888] |
| CORUM_small | +0.0014 [-0.0021, +0.0050] | +0.0015 [-0.0022, +0.0053] | +0.0017 [-0.0021, +0.0055] | -0.0016 [-0.0042, +0.0010] | -0.0013 [-0.0039, +0.0013] | -0.0033 [-0.0063, -0.0003] | -0.0034 [-0.0064, -0.0004] |
| CORUM_large | +0.0782 [+0.0539, +0.1026] | +0.0713 [+0.0478, +0.0947] | +0.0592 [+0.0371, +0.0813] | -0.0013 [-0.0206, +0.0179] | -0.0003 [-0.0198, +0.0192] | -0.1167 [-0.1522, -0.0812] | -0.1164 [-0.1536, -0.0792] |
| CollecTRI | -0.0081 [-0.0131, -0.0030] | -0.0058 [-0.0105, -0.0012] | -0.0074 [-0.0122, -0.0025] | -0.0099 [-0.0149, -0.0050] | -0.0098 [-0.0147, -0.0050] | -0.0136 [-0.0189, -0.0083] | -0.0124 [-0.0175, -0.0073] |

#### Difference in enrichment vs K10 attention (t = 9) [95% interval], raw

| truth | CONS-R | CONS-U | CONS-G | FULL-R | FULL-U | NULL-U | NULL-R |
|---|---|---|---|---|---|---|---|
| DoRothEA_AB | -0.01 [-0.04, +0.03] | +0.02 [-0.01, +0.05] | +0.01 [-0.02, +0.04] | -0.05 [-0.08, -0.02] | -0.03 [-0.06, +0.00] | -0.10 [-0.14, -0.06] | -0.09 [-0.13, -0.06] |
| CORUM_all | +3.16 [+1.54, +4.77] | +2.48 [+1.14, +3.82] | +0.98 [+0.27, +1.69] | -0.91 [-1.57, -0.25] | -0.95 [-1.64, -0.27] | -2.46 [-3.56, -1.36] | -2.51 [-3.62, -1.40] |
| KD_response | -0.66 [-0.77, -0.55] | -0.67 [-0.78, -0.56] | -0.79 [-0.90, -0.68] | -1.05 [-1.16, -0.95] | -0.96 [-1.07, -0.86] | -1.29 [-1.40, -1.17] | -1.37 [-1.50, -1.25] |
| CORUM_noRP | +0.42 [+0.02, +0.82] | +0.39 [-0.00, +0.79] | +0.33 [-0.03, +0.68] | -0.46 [-0.76, -0.16] | -0.41 [-0.71, -0.10] | -1.03 [-1.43, -0.63] | -1.05 [-1.45, -0.65] |
| CORUM_RP | +0.01 [-0.02, +0.03] | +0.02 [-0.01, +0.06] | +0.02 [-0.02, +0.05] | -0.05 [-0.09, -0.01] | -0.02 [-0.07, +0.02] | -0.20 [-0.27, -0.13] | -0.18 [-0.25, -0.12] |
| CORUM_small | +0.43 [-0.65, +1.52] | +0.47 [-0.67, +1.61] | +0.51 [-0.65, +1.66] | -0.48 [-1.27, +0.30] | -0.41 [-1.19, +0.38] | -1.01 [-1.92, -0.10] | -1.04 [-1.95, -0.13] |
| CORUM_large | +0.62 [+0.43, +0.81] | +0.56 [+0.38, +0.75] | +0.47 [+0.29, +0.64] | -0.01 [-0.16, +0.14] | -0.00 [-0.16, +0.15] | -0.92 [-1.21, -0.64] | -0.92 [-1.22, -0.63] |
| CollecTRI | -0.06 [-0.10, -0.02] | -0.05 [-0.08, -0.01] | -0.06 [-0.10, -0.02] | -0.08 [-0.12, -0.04] | -0.08 [-0.12, -0.04] | -0.11 [-0.15, -0.06] | -0.10 [-0.14, -0.06] |

#### Difference in AUPR vs CLR [95% paired delete-half gene jackknife], raw

| truth | CONS-R | CONS-U | CONS-G | FULL-R | FULL-U | NULL-U | NULL-R |
|---|---|---|---|---|---|---|---|
| DoRothEA_AB | +0.0052 [+0.0027, +0.0078] | +0.0099 [+0.0059, +0.0139] | +0.0086 [+0.0051, +0.0121] | -0.0033 [-0.0060, -0.0007] | +0.0010 [-0.0023, +0.0043] | -0.0120 [-0.0163, -0.0077] | -0.0109 [-0.0142, -0.0077] |
| CORUM_all | +0.0037 [-0.0065, +0.0138] | -0.0053 [-0.0150, +0.0045] | -0.0252 [-0.0403, -0.0101] | -0.0502 [-0.0729, -0.0276] | -0.0508 [-0.0743, -0.0273] | -0.0708 [-0.1017, -0.0399] | -0.0714 [-0.1024, -0.0404] |
| KD_response | -0.0558 [-0.0642, -0.0473] | -0.0566 [-0.0653, -0.0478] | -0.0661 [-0.0751, -0.0570] | -0.0868 [-0.0959, -0.0777] | -0.0798 [-0.0887, -0.0708] | -0.1052 [-0.1152, -0.0952] | -0.1120 [-0.1229, -0.1011] |
| CORUM_noRP | +0.0034 [-0.0004, +0.0072] | +0.0031 [-0.0006, +0.0068] | +0.0023 [-0.0013, +0.0058] | -0.0071 [-0.0107, -0.0034] | -0.0064 [-0.0101, -0.0027] | -0.0138 [-0.0189, -0.0087] | -0.0140 [-0.0191, -0.0090] |
| CORUM_RP | -0.0145 [-0.0309, +0.0019] | -0.0001 [-0.0224, +0.0223] | -0.0061 [-0.0368, +0.0247] | -0.0538 [-0.0894, -0.0181] | -0.0359 [-0.0719, +0.0002] | -0.1703 [-0.2313, -0.1092] | -0.1569 [-0.2135, -0.1003] |
| CORUM_small | +0.0009 [-0.0015, +0.0034] | +0.0011 [-0.0015, +0.0036] | +0.0012 [-0.0011, +0.0035] | -0.0021 [-0.0050, +0.0008] | -0.0018 [-0.0047, +0.0010] | -0.0038 [-0.0071, -0.0006] | -0.0039 [-0.0072, -0.0007] |
| CORUM_large | +0.0322 [+0.0164, +0.0480] | +0.0252 [+0.0104, +0.0400] | +0.0131 [-0.0050, +0.0312] | -0.0474 [-0.0681, -0.0267] | -0.0463 [-0.0678, -0.0248] | -0.1627 [-0.2093, -0.1162] | -0.1624 [-0.2108, -0.1141] |
| CollecTRI | +0.0006 [-0.0015, +0.0028] | +0.0029 [-0.0002, +0.0060] | +0.0014 [-0.0014, +0.0041] | -0.0012 [-0.0033, +0.0009] | -0.0011 [-0.0035, +0.0013] | -0.0049 [-0.0081, -0.0016] | -0.0037 [-0.0059, -0.0014] |

#### Difference in enrichment vs CLR [95% interval], raw

| truth | CONS-R | CONS-U | CONS-G | FULL-R | FULL-U | NULL-U | NULL-R |
|---|---|---|---|---|---|---|---|
| DoRothEA_AB | +0.03 [+0.01, +0.04] | +0.05 [+0.03, +0.08] | +0.05 [+0.03, +0.07] | -0.02 [-0.03, -0.00] | +0.01 [-0.01, +0.02] | -0.07 [-0.09, -0.04] | -0.06 [-0.08, -0.04] |
| CORUM_all | +0.28 [-0.49, +1.04] | -0.40 [-1.13, +0.34] | -1.90 [-3.04, -0.76] | -3.79 [-5.50, -2.08] | -3.83 [-5.61, -2.06] | -5.34 [-7.67, -3.01] | -5.39 [-7.73, -3.05] |
| KD_response | -0.71 [-0.82, -0.60] | -0.72 [-0.83, -0.61] | -0.84 [-0.95, -0.72] | -1.10 [-1.22, -0.99] | -1.01 [-1.13, -0.90] | -1.34 [-1.46, -1.21] | -1.42 [-1.56, -1.28] |
| CORUM_noRP | +0.29 [-0.03, +0.61] | +0.26 [-0.05, +0.57] | +0.19 [-0.11, +0.49] | -0.60 [-0.91, -0.28] | -0.54 [-0.85, -0.23] | -1.16 [-1.59, -0.74] | -1.18 [-1.61, -0.76] |
| CORUM_RP | -0.02 [-0.04, +0.00] | -0.00 [-0.03, +0.03] | -0.01 [-0.05, +0.03] | -0.07 [-0.12, -0.02] | -0.05 [-0.10, +0.00] | -0.23 [-0.31, -0.15] | -0.21 [-0.28, -0.13] |
| CORUM_small | +0.28 [-0.47, +1.04] | +0.32 [-0.45, +1.09] | +0.36 [-0.34, +1.06] | -0.63 [-1.50, +0.24] | -0.56 [-1.42, +0.31] | -1.16 [-2.15, -0.18] | -1.19 [-2.17, -0.20] |
| CORUM_large | +0.26 [+0.13, +0.38] | +0.20 [+0.08, +0.32] | +0.10 [-0.04, +0.25] | -0.38 [-0.54, -0.21] | -0.37 [-0.54, -0.20] | -1.29 [-1.66, -0.92] | -1.29 [-1.67, -0.90] |
| CollecTRI | +0.01 [-0.01, +0.02] | +0.02 [-0.00, +0.05] | +0.01 [-0.01, +0.03] | -0.01 [-0.03, +0.01] | -0.01 [-0.03, +0.01] | -0.04 [-0.06, -0.01] | -0.03 [-0.05, -0.01] |

#### Difference in AUPR vs K10 attention (t = 9) [95% paired delete-half gene jackknife], hub-stratified (pi of V0 attention)

| truth | CONS-R | CONS-U | CONS-G | FULL-R | FULL-U | NULL-U | NULL-R |
|---|---|---|---|---|---|---|---|
| DoRothEA_AB | -0.0002 [-0.0043, +0.0039] | +0.0033 [-0.0006, +0.0071] | +0.0029 [-0.0005, +0.0064] | -0.0066 [-0.0113, -0.0019] | -0.0030 [-0.0069, +0.0008] | -0.0178 [-0.0234, -0.0123] | -0.0163 [-0.0207, -0.0120] |
| CORUM_all | +0.0212 [+0.0088, +0.0335] | +0.0183 [+0.0078, +0.0288] | +0.0155 [+0.0082, +0.0228] | +0.0050 [+0.0004, +0.0095] | +0.0039 [-0.0004, +0.0082] | -0.0138 [-0.0194, -0.0082] | -0.0139 [-0.0195, -0.0083] |
| KD_response | +0.0022 [-0.0035, +0.0079] | -0.0004 [-0.0062, +0.0054] | -0.0001 [-0.0053, +0.0050] | +0.0152 [+0.0087, +0.0216] | -0.0088 [-0.0145, -0.0032] | -0.0434 [-0.0497, -0.0370] | -0.0415 [-0.0481, -0.0349] |
| CORUM_noRP | +0.0042 [+0.0015, +0.0069] | +0.0040 [+0.0014, +0.0067] | +0.0053 [+0.0027, +0.0080] | +0.0008 [-0.0009, +0.0025] | +0.0007 [-0.0011, +0.0025] | -0.0067 [-0.0090, -0.0044] | -0.0067 [-0.0090, -0.0045] |
| CORUM_RP | +0.0090 [-0.0136, +0.0317] | +0.0149 [-0.0144, +0.0441] | +0.0149 [-0.0141, +0.0438] | -0.0106 [-0.0373, +0.0161] | -0.0106 [-0.0424, +0.0212] | -0.1506 [-0.2045, -0.0966] | -0.1246 [-0.1730, -0.0761] |
| CORUM_small | +0.0014 [-0.0009, +0.0036] | +0.0013 [-0.0009, +0.0036] | +0.0015 [-0.0006, +0.0036] | -0.0002 [-0.0011, +0.0007] | -0.0001 [-0.0010, +0.0008] | -0.0020 [-0.0032, -0.0008] | -0.0021 [-0.0033, -0.0009] |
| CORUM_large | +0.0452 [+0.0197, +0.0706] | +0.0387 [+0.0153, +0.0620] | +0.0386 [+0.0179, +0.0594] | -0.0089 [-0.0269, +0.0091] | -0.0068 [-0.0247, +0.0111] | -0.1175 [-0.1505, -0.0845] | -0.1177 [-0.1516, -0.0838] |
| CollecTRI | +0.0036 [+0.0007, +0.0066] | +0.0059 [+0.0030, +0.0089] | +0.0053 [+0.0025, +0.0080] | +0.0089 [+0.0048, +0.0130] | +0.0019 [-0.0010, +0.0047] | -0.0046 [-0.0084, -0.0008] | -0.0007 [-0.0040, +0.0025] |

#### Difference in enrichment vs K10 attention (t = 9) [95% interval], hub-stratified (pi of V0 attention)

| truth | CONS-R | CONS-U | CONS-G | FULL-R | FULL-U | NULL-U | NULL-R |
|---|---|---|---|---|---|---|---|
| DoRothEA_AB | -0.00 [-0.02, +0.02] | +0.02 [-0.00, +0.04] | +0.02 [-0.00, +0.03] | -0.04 [-0.06, -0.01] | -0.02 [-0.04, +0.00] | -0.10 [-0.13, -0.07] | -0.09 [-0.11, -0.07] |
| CORUM_all | +1.60 [+0.67, +2.53] | +1.38 [+0.59, +2.17] | +1.17 [+0.62, +1.72] | +0.38 [+0.03, +0.72] | +0.29 [-0.03, +0.62] | -1.04 [-1.47, -0.62] | -1.05 [-1.47, -0.63] |
| KD_response | +0.03 [-0.04, +0.10] | -0.00 [-0.08, +0.07] | -0.00 [-0.07, +0.06] | +0.19 [+0.11, +0.27] | -0.11 [-0.18, -0.04] | -0.55 [-0.63, -0.47] | -0.53 [-0.61, -0.44] |
| CORUM_noRP | +0.36 [+0.13, +0.58] | +0.34 [+0.11, +0.56] | +0.45 [+0.23, +0.67] | +0.07 [-0.08, +0.21] | +0.06 [-0.09, +0.21] | -0.56 [-0.76, -0.37] | -0.57 [-0.76, -0.38] |
| CORUM_RP | +0.01 [-0.02, +0.04] | +0.02 [-0.02, +0.06] | +0.02 [-0.02, +0.06] | -0.01 [-0.05, +0.02] | -0.01 [-0.06, +0.03] | -0.20 [-0.27, -0.13] | -0.17 [-0.23, -0.10] |
| CORUM_small | +0.42 [-0.27, +1.10] | +0.41 [-0.28, +1.09] | +0.46 [-0.18, +1.10] | -0.06 [-0.32, +0.20] | -0.02 [-0.30, +0.25] | -0.61 [-0.98, -0.23] | -0.63 [-1.00, -0.27] |
| CORUM_large | +0.36 [+0.16, +0.56] | +0.31 [+0.12, +0.49] | +0.31 [+0.14, +0.47] | -0.07 [-0.21, +0.07] | -0.05 [-0.20, +0.09] | -0.93 [-1.19, -0.67] | -0.93 [-1.20, -0.66] |
| CollecTRI | +0.03 [+0.01, +0.05] | +0.05 [+0.02, +0.07] | +0.04 [+0.02, +0.06] | +0.07 [+0.04, +0.10] | +0.01 [-0.01, +0.04] | -0.04 [-0.07, -0.01] | -0.01 [-0.03, +0.02] |

#### Difference in AUPR vs CLR [95% paired delete-half gene jackknife], hub-stratified (pi of V0 attention)

| truth | CONS-R | CONS-U | CONS-G | FULL-R | FULL-U | NULL-U | NULL-R |
|---|---|---|---|---|---|---|---|
| DoRothEA_AB | +0.0055 [+0.0027, +0.0083] | +0.0089 [+0.0048, +0.0131] | +0.0086 [+0.0050, +0.0122] | -0.0009 [-0.0052, +0.0034] | +0.0026 [-0.0009, +0.0062] | -0.0122 [-0.0165, -0.0078] | -0.0107 [-0.0138, -0.0075] |
| CORUM_all | +0.0036 [-0.0027, +0.0100] | +0.0007 [-0.0050, +0.0064] | -0.0020 [-0.0086, +0.0045] | -0.0126 [-0.0218, -0.0034] | -0.0137 [-0.0234, -0.0040] | -0.0314 [-0.0467, -0.0160] | -0.0315 [-0.0468, -0.0161] |
| KD_response | -0.0223 [-0.0282, -0.0163] | -0.0249 [-0.0315, -0.0182] | -0.0246 [-0.0313, -0.0180] | -0.0093 [-0.0153, -0.0033] | -0.0333 [-0.0400, -0.0266] | -0.0678 [-0.0757, -0.0599] | -0.0659 [-0.0742, -0.0577] |
| CORUM_noRP | +0.0027 [+0.0005, +0.0049] | +0.0025 [+0.0004, +0.0046] | +0.0038 [+0.0015, +0.0061] | -0.0007 [-0.0026, +0.0011] | -0.0008 [-0.0026, +0.0009] | -0.0082 [-0.0109, -0.0056] | -0.0082 [-0.0108, -0.0056] |
| CORUM_RP | -0.0105 [-0.0351, +0.0142] | -0.0046 [-0.0347, +0.0255] | -0.0046 [-0.0376, +0.0284] | -0.0301 [-0.0660, +0.0059] | -0.0301 [-0.0668, +0.0066] | -0.1701 [-0.2315, -0.1086] | -0.1441 [-0.1992, -0.0890] |
| CORUM_small | +0.0008 [-0.0009, +0.0025] | +0.0008 [-0.0010, +0.0025] | +0.0010 [-0.0006, +0.0025] | -0.0008 [-0.0026, +0.0010] | -0.0006 [-0.0024, +0.0012] | -0.0026 [-0.0046, -0.0005] | -0.0026 [-0.0047, -0.0006] |
| CORUM_large | +0.0138 [-0.0039, +0.0315] | +0.0073 [-0.0090, +0.0236] | +0.0073 [-0.0108, +0.0254] | -0.0402 [-0.0602, -0.0202] | -0.0381 [-0.0594, -0.0169] | -0.1489 [-0.1947, -0.1030] | -0.1490 [-0.1958, -0.1022] |
| CollecTRI | +0.0051 [+0.0026, +0.0076] | +0.0074 [+0.0038, +0.0110] | +0.0067 [+0.0036, +0.0099] | +0.0103 [+0.0061, +0.0146] | +0.0033 [+0.0006, +0.0061] | -0.0032 [-0.0066, +0.0003] | +0.0007 [-0.0020, +0.0034] |

#### Difference in enrichment vs CLR [95% interval], hub-stratified (pi of V0 attention)

| truth | CONS-R | CONS-U | CONS-G | FULL-R | FULL-U | NULL-U | NULL-R |
|---|---|---|---|---|---|---|---|
| DoRothEA_AB | +0.03 [+0.01, +0.05] | +0.05 [+0.03, +0.07] | +0.05 [+0.03, +0.07] | -0.01 [-0.03, +0.02] | +0.01 [-0.00, +0.03] | -0.07 [-0.09, -0.04] | -0.06 [-0.08, -0.04] |
| CORUM_all | +0.27 [-0.20, +0.75] | +0.06 [-0.37, +0.49] | -0.15 [-0.65, +0.34] | -0.95 [-1.64, -0.25] | -1.03 [-1.77, -0.30] | -2.37 [-3.52, -1.21] | -2.37 [-3.53, -1.22] |
| KD_response | -0.28 [-0.36, -0.21] | -0.32 [-0.40, -0.23] | -0.31 [-0.40, -0.23] | -0.12 [-0.19, -0.04] | -0.42 [-0.51, -0.34] | -0.86 [-0.96, -0.76] | -0.84 [-0.94, -0.73] |
| CORUM_noRP | +0.23 [+0.04, +0.41] | +0.21 [+0.03, +0.39] | +0.32 [+0.13, +0.52] | -0.06 [-0.22, +0.10] | -0.07 [-0.22, +0.08] | -0.69 [-0.92, -0.47] | -0.69 [-0.91, -0.48] |
| CORUM_RP | -0.01 [-0.05, +0.02] | -0.01 [-0.05, +0.03] | -0.01 [-0.05, +0.04] | -0.04 [-0.09, +0.01] | -0.04 [-0.09, +0.01] | -0.23 [-0.31, -0.14] | -0.19 [-0.27, -0.12] |
| CORUM_small | +0.24 [-0.26, +0.75] | +0.23 [-0.30, +0.77] | +0.29 [-0.19, +0.77] | -0.23 [-0.78, +0.32] | -0.19 [-0.74, +0.35] | -0.78 [-1.39, -0.16] | -0.80 [-1.42, -0.19] |
| CORUM_large | +0.11 [-0.03, +0.25] | +0.06 [-0.07, +0.19] | +0.06 [-0.09, +0.20] | -0.32 [-0.48, -0.16] | -0.30 [-0.47, -0.13] | -1.18 [-1.54, -0.82] | -1.18 [-1.55, -0.81] |
| CollecTRI | +0.04 [+0.02, +0.06] | +0.06 [+0.03, +0.09] | +0.05 [+0.03, +0.08] | +0.08 [+0.05, +0.11] | +0.03 [+0.00, +0.05] | -0.02 [-0.05, +0.00] | +0.01 [-0.02, +0.03] |

#### Difference in AUPR vs K10 attention (t = 9) [95% paired delete-half gene jackknife], hub-stratified (MI row mean)

| truth | CONS-R | CONS-U | CONS-G | FULL-R | FULL-U | NULL-U | NULL-R |
|---|---|---|---|---|---|---|---|
| DoRothEA_AB | -0.0012 [-0.0070, +0.0046] | +0.0015 [-0.0037, +0.0068] | -0.0005 [-0.0059, +0.0048] | -0.0085 [-0.0143, -0.0026] | -0.0069 [-0.0126, -0.0013] | -0.0207 [-0.0282, -0.0132] | -0.0211 [-0.0278, -0.0144] |
| CORUM_all | +0.0368 [+0.0164, +0.0573] | +0.0283 [+0.0119, +0.0448] | +0.0133 [+0.0042, +0.0223] | -0.0065 [-0.0135, +0.0005] | -0.0075 [-0.0144, -0.0007] | -0.0256 [-0.0358, -0.0154] | -0.0259 [-0.0360, -0.0157] |
| KD_response | -0.0386 [-0.0457, -0.0315] | -0.0385 [-0.0455, -0.0315] | -0.0486 [-0.0555, -0.0417] | -0.0704 [-0.0781, -0.0628] | -0.0607 [-0.0681, -0.0533] | -0.0880 [-0.0972, -0.0788] | -0.0908 [-0.1003, -0.0813] |
| CORUM_noRP | +0.0048 [+0.0003, +0.0092] | +0.0039 [-0.0005, +0.0083] | +0.0025 [-0.0016, +0.0067] | -0.0058 [-0.0097, -0.0020] | -0.0057 [-0.0097, -0.0018] | -0.0137 [-0.0191, -0.0084] | -0.0139 [-0.0193, -0.0086] |
| CORUM_RP | +0.0143 [-0.0072, +0.0358] | +0.0207 [-0.0067, +0.0482] | +0.0158 [-0.0142, +0.0458] | -0.0199 [-0.0506, +0.0108] | -0.0123 [-0.0465, +0.0220] | -0.1490 [-0.2014, -0.0965] | -0.1254 [-0.1725, -0.0782] |
| CORUM_small | +0.0017 [-0.0019, +0.0052] | +0.0019 [-0.0018, +0.0056] | +0.0019 [-0.0019, +0.0057] | -0.0015 [-0.0039, +0.0010] | -0.0013 [-0.0038, +0.0011] | -0.0036 [-0.0065, -0.0007] | -0.0036 [-0.0065, -0.0008] |
| CORUM_large | +0.0837 [+0.0548, +0.1126] | +0.0706 [+0.0444, +0.0968] | +0.0585 [+0.0344, +0.0827] | +0.0180 [-0.0035, +0.0395] | +0.0065 [-0.0136, +0.0266] | -0.1096 [-0.1409, -0.0784] | -0.1094 [-0.1412, -0.0775] |
| CollecTRI | -0.0099 [-0.0153, -0.0045] | -0.0086 [-0.0135, -0.0037] | -0.0104 [-0.0156, -0.0052] | -0.0122 [-0.0178, -0.0065] | -0.0122 [-0.0174, -0.0069] | -0.0157 [-0.0215, -0.0100] | -0.0156 [-0.0213, -0.0100] |

#### Difference in enrichment vs K10 attention (t = 9) [95% interval], hub-stratified (MI row mean)

| truth | CONS-R | CONS-U | CONS-G | FULL-R | FULL-U | NULL-U | NULL-R |
|---|---|---|---|---|---|---|---|
| DoRothEA_AB | -0.01 [-0.04, +0.03] | +0.01 [-0.02, +0.04] | -0.00 [-0.03, +0.03] | -0.05 [-0.08, -0.01] | -0.04 [-0.07, -0.01] | -0.11 [-0.15, -0.07] | -0.12 [-0.15, -0.08] |
| CORUM_all | +2.78 [+1.24, +4.32] | +2.14 [+0.89, +3.38] | +1.00 [+0.32, +1.68] | -0.49 [-1.01, +0.04] | -0.57 [-1.09, -0.05] | -1.93 [-2.70, -1.16] | -1.95 [-2.72, -1.18] |
| KD_response | -0.49 [-0.58, -0.40] | -0.49 [-0.58, -0.40] | -0.62 [-0.71, -0.53] | -0.89 [-0.99, -0.80] | -0.77 [-0.87, -0.68] | -1.12 [-1.23, -1.00] | -1.15 [-1.27, -1.03] |
| CORUM_noRP | +0.40 [+0.03, +0.78] | +0.33 [-0.04, +0.70] | +0.21 [-0.13, +0.56] | -0.49 [-0.82, -0.16] | -0.48 [-0.81, -0.15] | -1.16 [-1.61, -0.71] | -1.18 [-1.62, -0.73] |
| CORUM_RP | +0.02 [-0.01, +0.05] | +0.03 [-0.01, +0.06] | +0.02 [-0.02, +0.06] | -0.03 [-0.07, +0.01] | -0.02 [-0.06, +0.03] | -0.20 [-0.27, -0.13] | -0.17 [-0.23, -0.10] |
| CORUM_small | +0.50 [-0.56, +1.57] | +0.57 [-0.55, +1.69] | +0.57 [-0.58, +1.73] | -0.44 [-1.17, +0.29] | -0.41 [-1.15, +0.33] | -1.08 [-1.96, -0.21] | -1.11 [-1.98, -0.23] |
| CORUM_large | +0.66 [+0.43, +0.89] | +0.56 [+0.35, +0.77] | +0.46 [+0.27, +0.65] | +0.14 [-0.03, +0.31] | +0.05 [-0.11, +0.21] | -0.87 [-1.12, -0.62] | -0.87 [-1.12, -0.61] |
| CollecTRI | -0.08 [-0.12, -0.03] | -0.07 [-0.11, -0.03] | -0.08 [-0.12, -0.04] | -0.10 [-0.14, -0.05] | -0.10 [-0.14, -0.05] | -0.12 [-0.17, -0.08] | -0.12 [-0.17, -0.08] |

#### Difference in AUPR vs CLR [95% paired delete-half gene jackknife], hub-stratified (MI row mean)

| truth | CONS-R | CONS-U | CONS-G | FULL-R | FULL-U | NULL-U | NULL-R |
|---|---|---|---|---|---|---|---|
| DoRothEA_AB | +0.0018 [-0.0006, +0.0042] | +0.0045 [+0.0017, +0.0073] | +0.0024 [+0.0001, +0.0047] | -0.0055 [-0.0080, -0.0030] | -0.0040 [-0.0066, -0.0014] | -0.0178 [-0.0226, -0.0130] | -0.0182 [-0.0215, -0.0148] |
| CORUM_all | +0.0094 [-0.0008, +0.0196] | +0.0009 [-0.0080, +0.0098] | -0.0141 [-0.0258, -0.0024] | -0.0339 [-0.0512, -0.0165] | -0.0349 [-0.0530, -0.0168] | -0.0530 [-0.0768, -0.0292] | -0.0532 [-0.0770, -0.0295] |
| KD_response | -0.0391 [-0.0441, -0.0341] | -0.0390 [-0.0440, -0.0341] | -0.0491 [-0.0541, -0.0441] | -0.0709 [-0.0772, -0.0647] | -0.0612 [-0.0669, -0.0556] | -0.0885 [-0.0970, -0.0800] | -0.0913 [-0.1000, -0.0826] |
| CORUM_noRP | +0.0033 [-0.0002, +0.0068] | +0.0025 [-0.0009, +0.0058] | +0.0011 [-0.0020, +0.0043] | -0.0072 [-0.0110, -0.0034] | -0.0071 [-0.0109, -0.0033] | -0.0151 [-0.0206, -0.0097] | -0.0154 [-0.0208, -0.0099] |
| CORUM_RP | -0.0088 [-0.0264, +0.0089] | -0.0023 [-0.0273, +0.0226] | -0.0073 [-0.0392, +0.0246] | -0.0429 [-0.0803, -0.0055] | -0.0354 [-0.0732, +0.0025] | -0.1720 [-0.2343, -0.1097] | -0.1484 [-0.2035, -0.0933] |
| CORUM_small | +0.0009 [-0.0018, +0.0035] | +0.0011 [-0.0016, +0.0038] | +0.0011 [-0.0014, +0.0036] | -0.0022 [-0.0050, +0.0005] | -0.0021 [-0.0049, +0.0006] | -0.0044 [-0.0077, -0.0010] | -0.0044 [-0.0078, -0.0011] |
| CORUM_large | +0.0267 [+0.0109, +0.0425] | +0.0136 [+0.0002, +0.0270] | +0.0015 [-0.0137, +0.0168] | -0.0390 [-0.0580, -0.0200] | -0.0505 [-0.0695, -0.0315] | -0.1666 [-0.2109, -0.1224] | -0.1664 [-0.2110, -0.1218] |
| CollecTRI | -0.0025 [-0.0046, -0.0003] | -0.0012 [-0.0036, +0.0012] | -0.0030 [-0.0051, -0.0009] | -0.0047 [-0.0073, -0.0022] | -0.0047 [-0.0066, -0.0028] | -0.0083 [-0.0115, -0.0050] | -0.0082 [-0.0102, -0.0061] |

#### Difference in enrichment vs CLR [95% interval], hub-stratified (MI row mean)

| truth | CONS-R | CONS-U | CONS-G | FULL-R | FULL-U | NULL-U | NULL-R |
|---|---|---|---|---|---|---|---|
| DoRothEA_AB | +0.01 [-0.00, +0.02] | +0.02 [+0.01, +0.04] | +0.01 [+0.00, +0.03] | -0.03 [-0.04, -0.02] | -0.02 [-0.04, -0.01] | -0.10 [-0.12, -0.07] | -0.10 [-0.12, -0.08] |
| CORUM_all | +0.71 [-0.06, +1.48] | +0.07 [-0.60, +0.74] | -1.06 [-1.95, -0.18] | -2.55 [-3.86, -1.25] | -2.63 [-4.00, -1.27] | -4.00 [-5.79, -2.20] | -4.02 [-5.81, -2.22] |
| KD_response | -0.50 [-0.56, -0.43] | -0.50 [-0.56, -0.43] | -0.62 [-0.69, -0.56] | -0.90 [-0.98, -0.82] | -0.78 [-0.85, -0.71] | -1.12 [-1.23, -1.02] | -1.16 [-1.27, -1.05] |
| CORUM_noRP | +0.28 [-0.01, +0.58] | +0.21 [-0.07, +0.49] | +0.09 [-0.17, +0.36] | -0.61 [-0.93, -0.29] | -0.60 [-0.92, -0.28] | -1.28 [-1.74, -0.81] | -1.29 [-1.75, -0.84] |
| CORUM_RP | -0.01 [-0.04, +0.01] | -0.00 [-0.04, +0.03] | -0.01 [-0.05, +0.03] | -0.06 [-0.11, -0.01] | -0.05 [-0.10, +0.00] | -0.23 [-0.31, -0.15] | -0.20 [-0.27, -0.12] |
| CORUM_small | +0.26 [-0.55, +1.07] | +0.33 [-0.50, +1.16] | +0.33 [-0.44, +1.11] | -0.68 [-1.53, +0.16] | -0.65 [-1.48, +0.19] | -1.33 [-2.34, -0.31] | -1.35 [-2.36, -0.33] |
| CORUM_large | +0.21 [+0.09, +0.34] | +0.11 [+0.00, +0.21] | +0.01 [-0.11, +0.13] | -0.31 [-0.46, -0.16] | -0.40 [-0.55, -0.25] | -1.32 [-1.67, -0.97] | -1.32 [-1.67, -0.96] |
| CollecTRI | -0.02 [-0.04, -0.00] | -0.01 [-0.03, +0.01] | -0.02 [-0.04, -0.01] | -0.04 [-0.06, -0.02] | -0.04 [-0.05, -0.02] | -0.06 [-0.09, -0.04] | -0.06 [-0.08, -0.05] |



### Idea 2

#### Pooled enrichment (AUPR / base rate), raw

| truth | CLR | K10 | G-SOFT | G-HARD | HALF-SOFT | R-SOFT | PERM-1 | PERM-2 | PERM-3 | PERM-4 | PERM-5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB (base 0.1830) | 1.06 | 1.09 | 1.08 | 1.08 | 1.09 | 1.04 | 1.06 | 1.06 | 1.04 | 1.04 | 1.05 |
| CORUM_all (base 0.0133) | 6.39 | 3.51 | 3.66 | 3.80 | 3.78 | 4.05 | 2.22 | 2.12 | 2.36 | 2.29 | 2.00 |
| KD_response (base 0.0787) | 2.42 | 2.38 | 2.59 | 2.66 | 2.67 | 2.10 | 1.76 | 1.77 | 1.76 | 1.67 | 1.68 |
| CORUM_noRP (base 0.0119) | 2.18 | 2.05 | 2.13 | 2.15 | 2.21 | 2.08 | 1.56 | 1.50 | 1.58 | 1.68 | 1.43 |
| CORUM_RP (base 0.7502) | 1.21 | 1.18 | 1.19 | 1.20 | 1.20 | 1.22 | 1.09 | 1.10 | 1.10 | 1.09 | 1.10 |
| CORUM_small (base 0.0033) | 2.16 | 2.01 | 1.96 | 2.01 | 1.96 | 1.93 | 1.68 | 1.60 | 1.68 | 1.47 | 1.57 |
| CORUM_large (base 0.1263) | 2.29 | 1.93 | 1.93 | 1.91 | 1.96 | 2.00 | 1.56 | 1.48 | 1.56 | 1.54 | 1.46 |
| CollecTRI (base 0.1278) | 1.03 | 1.10 | 1.11 | 1.11 | 1.12 | 1.06 | 1.04 | 1.08 | 1.06 | 1.02 | 1.05 |

#### Pooled AUPR, raw

| truth | CLR | K10 | G-SOFT | G-HARD | HALF-SOFT | R-SOFT | PERM-1 | PERM-2 | PERM-3 | PERM-4 | PERM-5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB | 0.19395 | 0.20032 | 0.19847 | 0.19793 | 0.20022 | 0.19080 | 0.19354 | 0.19476 | 0.18993 | 0.18947 | 0.19203 |
| CORUM_all | 0.08471 | 0.04654 | 0.04851 | 0.05039 | 0.05013 | 0.05364 | 0.02941 | 0.02807 | 0.03128 | 0.03036 | 0.02648 |
| KD_response | 0.19084 | 0.18698 | 0.20368 | 0.20958 | 0.20994 | 0.16537 | 0.13833 | 0.13930 | 0.13863 | 0.13179 | 0.13209 |
| CORUM_noRP | 0.02590 | 0.02432 | 0.02522 | 0.02548 | 0.02616 | 0.02466 | 0.01848 | 0.01782 | 0.01871 | 0.01989 | 0.01699 |
| CORUM_RP | 0.90638 | 0.88797 | 0.89585 | 0.89832 | 0.89752 | 0.91325 | 0.81824 | 0.82518 | 0.82410 | 0.81437 | 0.82150 |
| CORUM_small | 0.00712 | 0.00663 | 0.00645 | 0.00663 | 0.00646 | 0.00636 | 0.00552 | 0.00526 | 0.00554 | 0.00484 | 0.00518 |
| CORUM_large | 0.28925 | 0.24321 | 0.24318 | 0.24090 | 0.24718 | 0.25204 | 0.19663 | 0.18706 | 0.19665 | 0.19462 | 0.18375 |
| CollecTRI | 0.13150 | 0.14021 | 0.14144 | 0.14136 | 0.14304 | 0.13567 | 0.13236 | 0.13798 | 0.13590 | 0.13049 | 0.13399 |

#### Pooled enrichment (AUPR / base rate), hub-stratified (pi of V0 attention)

| truth | CLR | K10 | G-SOFT | G-HARD | HALF-SOFT | R-SOFT | PERM-1 | PERM-2 | PERM-3 | PERM-4 | PERM-5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB (base 0.1830) | 1.06 | 1.09 | 1.08 | 1.08 | 1.10 | 1.02 | 1.05 | 1.06 | 1.03 | 1.03 | 1.04 |
| CORUM_all (base 0.0133) | 3.40 | 2.08 | 2.25 | 2.31 | 2.33 | 2.56 | 1.56 | 1.55 | 1.69 | 1.66 | 1.51 |
| KD_response (base 0.0787) | 1.91 | 1.60 | 1.64 | 1.69 | 1.68 | 1.45 | 1.36 | 1.36 | 1.40 | 1.36 | 1.37 |
| CORUM_noRP (base 0.0119) | 1.71 | 1.58 | 1.72 | 1.74 | 1.81 | 1.75 | 1.33 | 1.32 | 1.37 | 1.42 | 1.28 |
| CORUM_RP (base 0.7502) | 1.21 | 1.18 | 1.19 | 1.19 | 1.19 | 1.21 | 1.09 | 1.10 | 1.10 | 1.09 | 1.11 |
| CORUM_small (base 0.0033) | 1.78 | 1.61 | 1.69 | 1.67 | 1.73 | 1.69 | 1.45 | 1.39 | 1.40 | 1.33 | 1.38 |
| CORUM_large (base 0.1263) | 2.18 | 1.94 | 1.92 | 1.86 | 1.95 | 1.95 | 1.53 | 1.47 | 1.56 | 1.53 | 1.46 |
| CollecTRI (base 0.1278) | 1.01 | 1.02 | 1.03 | 1.04 | 1.03 | 1.00 | 0.99 | 1.02 | 1.02 | 0.99 | 1.01 |

#### Pooled AUPR, hub-stratified (pi of V0 attention)

| truth | CLR | K10 | G-SOFT | G-HARD | HALF-SOFT | R-SOFT | PERM-1 | PERM-2 | PERM-3 | PERM-4 | PERM-5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB | 0.19385 | 0.19952 | 0.19822 | 0.19680 | 0.20102 | 0.18610 | 0.19242 | 0.19333 | 0.18890 | 0.18896 | 0.18985 |
| CORUM_all | 0.04508 | 0.02752 | 0.02984 | 0.03059 | 0.03095 | 0.03397 | 0.02063 | 0.02060 | 0.02238 | 0.02203 | 0.02000 |
| KD_response | 0.15037 | 0.12590 | 0.12887 | 0.13294 | 0.13186 | 0.11413 | 0.10694 | 0.10738 | 0.10999 | 0.10673 | 0.10811 |
| CORUM_noRP | 0.02024 | 0.01871 | 0.02041 | 0.02069 | 0.02141 | 0.02081 | 0.01572 | 0.01569 | 0.01630 | 0.01683 | 0.01524 |
| CORUM_RP | 0.90395 | 0.88445 | 0.89032 | 0.89305 | 0.89208 | 0.91088 | 0.81691 | 0.82277 | 0.82800 | 0.81774 | 0.83214 |
| CORUM_small | 0.00585 | 0.00529 | 0.00555 | 0.00549 | 0.00569 | 0.00556 | 0.00478 | 0.00459 | 0.00463 | 0.00439 | 0.00453 |
| CORUM_large | 0.27577 | 0.24444 | 0.24257 | 0.23536 | 0.24683 | 0.24609 | 0.19380 | 0.18505 | 0.19731 | 0.19311 | 0.18448 |
| CollecTRI | 0.12845 | 0.12990 | 0.13144 | 0.13273 | 0.13211 | 0.12727 | 0.12608 | 0.13089 | 0.13061 | 0.12651 | 0.12899 |

#### Pooled enrichment (AUPR / base rate), hub-stratified (MI row mean)

| truth | CLR | K10 | G-SOFT | G-HARD | HALF-SOFT | R-SOFT | PERM-1 | PERM-2 | PERM-3 | PERM-4 | PERM-5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB (base 0.1830) | 1.10 | 1.11 | 1.11 | 1.11 | 1.13 | 1.06 | 1.07 | 1.07 | 1.04 | 1.04 | 1.06 |
| CORUM_all (base 0.0133) | 5.03 | 2.97 | 3.14 | 3.26 | 3.25 | 3.36 | 1.91 | 1.88 | 2.06 | 2.10 | 1.76 |
| KD_response (base 0.0787) | 2.19 | 2.18 | 2.33 | 2.41 | 2.40 | 1.95 | 1.65 | 1.69 | 1.66 | 1.59 | 1.59 |
| CORUM_noRP (base 0.0119) | 2.30 | 2.18 | 2.30 | 2.37 | 2.38 | 2.23 | 1.57 | 1.54 | 1.61 | 1.78 | 1.48 |
| CORUM_RP (base 0.7502) | 1.21 | 1.18 | 1.19 | 1.19 | 1.19 | 1.21 | 1.09 | 1.10 | 1.10 | 1.08 | 1.09 |
| CORUM_small (base 0.0033) | 2.32 | 2.08 | 2.07 | 2.15 | 2.07 | 2.04 | 1.74 | 1.62 | 1.75 | 1.54 | 1.64 |
| CORUM_large (base 0.1263) | 2.32 | 1.87 | 1.86 | 1.82 | 1.89 | 1.90 | 1.54 | 1.44 | 1.53 | 1.52 | 1.41 |
| CollecTRI (base 0.1278) | 1.06 | 1.12 | 1.15 | 1.16 | 1.17 | 1.09 | 1.04 | 1.09 | 1.08 | 1.02 | 1.06 |

#### Pooled AUPR, hub-stratified (MI row mean)

| truth | CLR | K10 | G-SOFT | G-HARD | HALF-SOFT | R-SOFT | PERM-1 | PERM-2 | PERM-3 | PERM-4 | PERM-5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB | 0.20052 | 0.20346 | 0.20381 | 0.20392 | 0.20631 | 0.19344 | 0.19538 | 0.19605 | 0.19072 | 0.19021 | 0.19401 |
| CORUM_all | 0.06673 | 0.03934 | 0.04158 | 0.04317 | 0.04315 | 0.04453 | 0.02530 | 0.02494 | 0.02730 | 0.02784 | 0.02339 |
| KD_response | 0.17224 | 0.17174 | 0.18375 | 0.18977 | 0.18912 | 0.15329 | 0.12986 | 0.13310 | 0.13045 | 0.12529 | 0.12516 |
| CORUM_noRP | 0.02724 | 0.02583 | 0.02730 | 0.02806 | 0.02824 | 0.02646 | 0.01866 | 0.01832 | 0.01907 | 0.02108 | 0.01753 |
| CORUM_RP | 0.90554 | 0.88247 | 0.88950 | 0.89080 | 0.89017 | 0.91122 | 0.81630 | 0.82782 | 0.82463 | 0.81266 | 0.81647 |
| CORUM_small | 0.00765 | 0.00686 | 0.00683 | 0.00709 | 0.00683 | 0.00672 | 0.00573 | 0.00534 | 0.00577 | 0.00508 | 0.00539 |
| CORUM_large | 0.29338 | 0.23637 | 0.23495 | 0.22981 | 0.23883 | 0.24032 | 0.19455 | 0.18238 | 0.19288 | 0.19224 | 0.17805 |
| CollecTRI | 0.13543 | 0.14287 | 0.14726 | 0.14882 | 0.14910 | 0.13940 | 0.13318 | 0.13903 | 0.13740 | 0.13094 | 0.13507 |

#### Difference in AUPR vs K10 attention (t = 9) [95% paired delete-half gene jackknife], raw

| truth | G-SOFT | G-HARD | HALF-SOFT | R-SOFT | PERM-1 | PERM-2 | PERM-3 | PERM-4 | PERM-5 |
|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB | -0.0019 [-0.0061, +0.0024] | -0.0024 [-0.0068, +0.0021] | -0.0001 [-0.0040, +0.0038] | -0.0095 [-0.0143, -0.0047] | -0.0068 [-0.0120, -0.0016] | -0.0056 [-0.0102, -0.0009] | -0.0104 [-0.0152, -0.0056] | -0.0109 [-0.0163, -0.0054] | -0.0083 [-0.0136, -0.0030] |
| CORUM_all | +0.0020 [+0.0005, +0.0035] | +0.0039 [+0.0013, +0.0064] | +0.0036 [+0.0020, +0.0052] | +0.0071 [+0.0018, +0.0124] | -0.0171 [-0.0258, -0.0084] | -0.0185 [-0.0291, -0.0079] | -0.0153 [-0.0242, -0.0063] | -0.0162 [-0.0256, -0.0068] | -0.0201 [-0.0303, -0.0098] |
| KD_response | +0.0167 [+0.0137, +0.0198] | +0.0226 [+0.0177, +0.0275] | +0.0230 [+0.0204, +0.0255] | -0.0216 [-0.0302, -0.0130] | -0.0486 [-0.0569, -0.0403] | -0.0477 [-0.0549, -0.0405] | -0.0483 [-0.0569, -0.0398] | -0.0552 [-0.0636, -0.0467] | -0.0549 [-0.0626, -0.0471] |
| CORUM_noRP | +0.0009 [-0.0002, +0.0020] | +0.0012 [-0.0004, +0.0027] | +0.0018 [+0.0008, +0.0029] | +0.0003 [-0.0019, +0.0026] | -0.0058 [-0.0093, -0.0024] | -0.0065 [-0.0105, -0.0025] | -0.0056 [-0.0094, -0.0018] | -0.0044 [-0.0074, -0.0015] | -0.0073 [-0.0109, -0.0038] |
| CORUM_RP | +0.0079 [-0.0029, +0.0187] | +0.0104 [-0.0014, +0.0221] | +0.0096 [-0.0005, +0.0196] | +0.0253 [-0.0023, +0.0529] | -0.0697 [-0.1096, -0.0299] | -0.0628 [-0.0992, -0.0264] | -0.0639 [-0.0963, -0.0314] | -0.0736 [-0.1053, -0.0419] | -0.0665 [-0.1052, -0.0277] |
| CORUM_small | -0.0002 [-0.0009, +0.0005] | +0.0000 [-0.0006, +0.0006] | -0.0002 [-0.0009, +0.0006] | -0.0003 [-0.0015, +0.0009] | -0.0011 [-0.0031, +0.0008] | -0.0014 [-0.0030, +0.0003] | -0.0011 [-0.0025, +0.0003] | -0.0018 [-0.0049, +0.0013] | -0.0015 [-0.0031, +0.0002] |
| CORUM_large | -0.0000 [-0.0055, +0.0055] | -0.0023 [-0.0096, +0.0050] | +0.0040 [-0.0009, +0.0088] | +0.0088 [-0.0048, +0.0224] | -0.0466 [-0.0716, -0.0216] | -0.0562 [-0.0816, -0.0307] | -0.0466 [-0.0682, -0.0249] | -0.0486 [-0.0705, -0.0267] | -0.0595 [-0.0853, -0.0336] |
| CollecTRI | +0.0012 [-0.0024, +0.0048] | +0.0011 [-0.0024, +0.0047] | +0.0028 [-0.0007, +0.0064] | -0.0045 [-0.0083, -0.0007] | -0.0079 [-0.0126, -0.0031] | -0.0022 [-0.0071, +0.0026] | -0.0043 [-0.0089, +0.0002] | -0.0097 [-0.0151, -0.0044] | -0.0062 [-0.0110, -0.0015] |

#### Difference in enrichment vs K10 attention (t = 9) [95% interval], raw

| truth | G-SOFT | G-HARD | HALF-SOFT | R-SOFT | PERM-1 | PERM-2 | PERM-3 | PERM-4 | PERM-5 |
|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB | -0.01 [-0.03, +0.01] | -0.01 [-0.04, +0.01] | -0.00 [-0.02, +0.02] | -0.05 [-0.08, -0.03] | -0.04 [-0.07, -0.01] | -0.03 [-0.06, -0.00] | -0.06 [-0.08, -0.03] | -0.06 [-0.09, -0.03] | -0.05 [-0.07, -0.02] |
| CORUM_all | +0.15 [+0.03, +0.26] | +0.29 [+0.10, +0.48] | +0.27 [+0.15, +0.39] | +0.54 [+0.14, +0.93] | -1.29 [-1.95, -0.64] | -1.39 [-2.19, -0.59] | -1.15 [-1.82, -0.48] | -1.22 [-1.93, -0.51] | -1.51 [-2.29, -0.74] |
| KD_response | +0.21 [+0.17, +0.25] | +0.29 [+0.22, +0.35] | +0.29 [+0.26, +0.32] | -0.27 [-0.38, -0.16] | -0.62 [-0.72, -0.51] | -0.61 [-0.70, -0.51] | -0.61 [-0.72, -0.51] | -0.70 [-0.81, -0.59] | -0.70 [-0.80, -0.60] |
| CORUM_noRP | +0.08 [-0.02, +0.17] | +0.10 [-0.03, +0.23] | +0.16 [+0.07, +0.24] | +0.03 [-0.16, +0.22] | -0.49 [-0.79, -0.20] | -0.55 [-0.89, -0.21] | -0.47 [-0.80, -0.15] | -0.37 [-0.62, -0.13] | -0.62 [-0.92, -0.32] |
| CORUM_RP | +0.01 [-0.00, +0.02] | +0.01 [-0.00, +0.03] | +0.01 [-0.00, +0.03] | +0.03 [-0.00, +0.07] | -0.09 [-0.15, -0.04] | -0.08 [-0.13, -0.04] | -0.09 [-0.13, -0.04] | -0.10 [-0.14, -0.06] | -0.09 [-0.14, -0.04] |
| CORUM_small | -0.05 [-0.26, +0.15] | +0.00 [-0.19, +0.19] | -0.05 [-0.29, +0.19] | -0.08 [-0.45, +0.28] | -0.34 [-0.93, +0.25] | -0.42 [-0.91, +0.08] | -0.33 [-0.75, +0.09] | -0.54 [-1.48, +0.39] | -0.44 [-0.94, +0.05] |
| CORUM_large | -0.00 [-0.04, +0.04] | -0.02 [-0.08, +0.04] | +0.03 [-0.01, +0.07] | +0.07 [-0.04, +0.18] | -0.37 [-0.57, -0.17] | -0.44 [-0.65, -0.24] | -0.37 [-0.54, -0.20] | -0.38 [-0.56, -0.21] | -0.47 [-0.68, -0.27] |
| CollecTRI | +0.01 [-0.02, +0.04] | +0.01 [-0.02, +0.04] | +0.02 [-0.01, +0.05] | -0.04 [-0.07, -0.01] | -0.06 [-0.10, -0.02] | -0.02 [-0.06, +0.02] | -0.03 [-0.07, +0.00] | -0.08 [-0.12, -0.03] | -0.05 [-0.09, -0.01] |

#### Difference in AUPR vs CLR [95% paired delete-half gene jackknife], raw

| truth | G-SOFT | G-HARD | HALF-SOFT | R-SOFT | PERM-1 | PERM-2 | PERM-3 | PERM-4 | PERM-5 |
|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB | +0.0045 [+0.0002, +0.0089] | +0.0040 [-0.0002, +0.0082] | +0.0063 [+0.0016, +0.0109] | -0.0032 [-0.0071, +0.0008] | -0.0004 [-0.0061, +0.0052] | +0.0008 [-0.0042, +0.0059] | -0.0040 [-0.0096, +0.0016] | -0.0045 [-0.0094, +0.0005] | -0.0019 [-0.0073, +0.0035] |
| CORUM_all | -0.0362 [-0.0549, -0.0175] | -0.0343 [-0.0522, -0.0164] | -0.0346 [-0.0529, -0.0162] | -0.0311 [-0.0469, -0.0153] | -0.0553 [-0.0808, -0.0298] | -0.0566 [-0.0830, -0.0303] | -0.0534 [-0.0786, -0.0282] | -0.0544 [-0.0800, -0.0287] | -0.0582 [-0.0841, -0.0323] |
| KD_response | +0.0128 [+0.0046, +0.0211] | +0.0187 [+0.0101, +0.0273] | +0.0191 [+0.0112, +0.0270] | -0.0255 [-0.0351, -0.0159] | -0.0525 [-0.0634, -0.0417] | -0.0515 [-0.0625, -0.0406] | -0.0522 [-0.0622, -0.0422] | -0.0590 [-0.0692, -0.0489] | -0.0587 [-0.0689, -0.0486] |
| CORUM_noRP | -0.0007 [-0.0033, +0.0020] | -0.0004 [-0.0031, +0.0022] | +0.0003 [-0.0023, +0.0028] | -0.0012 [-0.0045, +0.0020] | -0.0074 [-0.0116, -0.0033] | -0.0081 [-0.0129, -0.0032] | -0.0072 [-0.0115, -0.0029] | -0.0060 [-0.0101, -0.0020] | -0.0089 [-0.0129, -0.0049] |
| CORUM_RP | -0.0105 [-0.0322, +0.0111] | -0.0081 [-0.0298, +0.0137] | -0.0089 [-0.0279, +0.0102] | +0.0069 [-0.0161, +0.0299] | -0.0881 [-0.1343, -0.0419] | -0.0812 [-0.1271, -0.0353] | -0.0823 [-0.1239, -0.0407] | -0.0920 [-0.1341, -0.0499] | -0.0849 [-0.1292, -0.0406] |
| CORUM_small | -0.0007 [-0.0036, +0.0023] | -0.0005 [-0.0034, +0.0025] | -0.0007 [-0.0034, +0.0020] | -0.0008 [-0.0035, +0.0020] | -0.0016 [-0.0044, +0.0011] | -0.0019 [-0.0052, +0.0015] | -0.0016 [-0.0044, +0.0012] | -0.0023 [-0.0052, +0.0007] | -0.0019 [-0.0046, +0.0007] |
| CORUM_large | -0.0461 [-0.0655, -0.0267] | -0.0483 [-0.0665, -0.0301] | -0.0421 [-0.0615, -0.0226] | -0.0372 [-0.0564, -0.0180] | -0.0926 [-0.1296, -0.0557] | -0.1022 [-0.1380, -0.0664] | -0.0926 [-0.1240, -0.0612] | -0.0946 [-0.1276, -0.0617] | -0.1055 [-0.1403, -0.0707] |
| CollecTRI | +0.0099 [+0.0054, +0.0144] | +0.0099 [+0.0058, +0.0140] | +0.0115 [+0.0065, +0.0166] | +0.0042 [+0.0001, +0.0082] | +0.0009 [-0.0038, +0.0055] | +0.0065 [+0.0015, +0.0115] | +0.0044 [-0.0010, +0.0098] | -0.0010 [-0.0057, +0.0037] | +0.0025 [-0.0027, +0.0076] |

#### Difference in enrichment vs CLR [95% interval], raw

| truth | G-SOFT | G-HARD | HALF-SOFT | R-SOFT | PERM-1 | PERM-2 | PERM-3 | PERM-4 | PERM-5 |
|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB | +0.02 [+0.00, +0.05] | +0.02 [-0.00, +0.04] | +0.03 [+0.01, +0.06] | -0.02 [-0.04, +0.00] | -0.00 [-0.03, +0.03] | +0.00 [-0.02, +0.03] | -0.02 [-0.05, +0.01] | -0.02 [-0.05, +0.00] | -0.01 [-0.04, +0.02] |
| CORUM_all | -2.73 [-4.14, -1.32] | -2.59 [-3.94, -1.24] | -2.61 [-3.99, -1.23] | -2.34 [-3.53, -1.15] | -4.17 [-6.09, -2.25] | -4.27 [-6.26, -2.28] | -4.03 [-5.93, -2.13] | -4.10 [-6.04, -2.16] | -4.39 [-6.35, -2.44] |
| KD_response | +0.16 [+0.06, +0.27] | +0.24 [+0.13, +0.35] | +0.24 [+0.14, +0.34] | -0.32 [-0.45, -0.20] | -0.67 [-0.80, -0.53] | -0.65 [-0.79, -0.52] | -0.66 [-0.79, -0.54] | -0.75 [-0.88, -0.62] | -0.75 [-0.88, -0.62] |
| CORUM_noRP | -0.06 [-0.28, +0.17] | -0.04 [-0.26, +0.19] | +0.02 [-0.20, +0.24] | -0.11 [-0.38, +0.17] | -0.63 [-0.98, -0.28] | -0.68 [-1.09, -0.27] | -0.61 [-0.97, -0.25] | -0.51 [-0.85, -0.17] | -0.75 [-1.09, -0.41] |
| CORUM_RP | -0.01 [-0.04, +0.01] | -0.01 [-0.04, +0.02] | -0.01 [-0.04, +0.01] | +0.01 [-0.02, +0.04] | -0.12 [-0.18, -0.06] | -0.11 [-0.17, -0.05] | -0.11 [-0.17, -0.05] | -0.12 [-0.18, -0.07] | -0.11 [-0.17, -0.05] |
| CORUM_small | -0.20 [-1.10, +0.70] | -0.15 [-1.04, +0.75] | -0.20 [-1.02, +0.62] | -0.23 [-1.07, +0.61] | -0.49 [-1.32, +0.35] | -0.57 [-1.58, +0.45] | -0.48 [-1.33, +0.37] | -0.69 [-1.58, +0.20] | -0.59 [-1.39, +0.21] |
| CORUM_large | -0.36 [-0.52, -0.21] | -0.38 [-0.53, -0.24] | -0.33 [-0.49, -0.18] | -0.29 [-0.45, -0.14] | -0.73 [-1.03, -0.44] | -0.81 [-1.09, -0.53] | -0.73 [-0.98, -0.48] | -0.75 [-1.01, -0.49] | -0.84 [-1.11, -0.56] |
| CollecTRI | +0.08 [+0.04, +0.11] | +0.08 [+0.05, +0.11] | +0.09 [+0.05, +0.13] | +0.03 [+0.00, +0.06] | +0.01 [-0.03, +0.04] | +0.05 [+0.01, +0.09] | +0.03 [-0.01, +0.08] | -0.01 [-0.04, +0.03] | +0.02 [-0.02, +0.06] |

#### Difference in AUPR vs K10 attention (t = 9) [95% paired delete-half gene jackknife], hub-stratified (pi of V0 attention)

| truth | G-SOFT | G-HARD | HALF-SOFT | R-SOFT | PERM-1 | PERM-2 | PERM-3 | PERM-4 | PERM-5 |
|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB | -0.0013 [-0.0078, +0.0052] | -0.0027 [-0.0100, +0.0045] | +0.0015 [-0.0045, +0.0075] | -0.0134 [-0.0202, -0.0066] | -0.0071 [-0.0124, -0.0018] | -0.0062 [-0.0111, -0.0013] | -0.0106 [-0.0159, -0.0054] | -0.0106 [-0.0158, -0.0054] | -0.0097 [-0.0153, -0.0040] |
| CORUM_all | +0.0023 [+0.0004, +0.0042] | +0.0031 [+0.0006, +0.0056] | +0.0034 [+0.0016, +0.0052] | +0.0065 [+0.0029, +0.0100] | -0.0069 [-0.0104, -0.0034] | -0.0069 [-0.0111, -0.0027] | -0.0051 [-0.0086, -0.0017] | -0.0055 [-0.0088, -0.0022] | -0.0075 [-0.0116, -0.0034] |
| KD_response | +0.0030 [-0.0003, +0.0063] | +0.0070 [+0.0022, +0.0118] | +0.0060 [+0.0031, +0.0088] | -0.0118 [-0.0183, -0.0053] | -0.0190 [-0.0247, -0.0133] | -0.0185 [-0.0236, -0.0135] | -0.0159 [-0.0213, -0.0105] | -0.0192 [-0.0247, -0.0136] | -0.0178 [-0.0231, -0.0125] |
| CORUM_noRP | +0.0017 [+0.0003, +0.0031] | +0.0020 [+0.0001, +0.0039] | +0.0027 [+0.0014, +0.0040] | +0.0021 [+0.0001, +0.0041] | -0.0030 [-0.0048, -0.0012] | -0.0030 [-0.0052, -0.0008] | -0.0024 [-0.0043, -0.0005] | -0.0019 [-0.0035, -0.0003] | -0.0035 [-0.0052, -0.0018] |
| CORUM_RP | +0.0059 [-0.0062, +0.0179] | +0.0086 [-0.0047, +0.0219] | +0.0076 [-0.0035, +0.0187] | +0.0264 [-0.0044, +0.0572] | -0.0675 [-0.1113, -0.0238] | -0.0617 [-0.1008, -0.0225] | -0.0565 [-0.0877, -0.0252] | -0.0667 [-0.0984, -0.0351] | -0.0523 [-0.0933, -0.0114] |
| CORUM_small | +0.0003 [-0.0003, +0.0009] | +0.0002 [-0.0005, +0.0009] | +0.0004 [-0.0001, +0.0009] | +0.0003 [-0.0006, +0.0012] | -0.0005 [-0.0014, +0.0003] | -0.0007 [-0.0015, +0.0001] | -0.0007 [-0.0013, +0.0000] | -0.0009 [-0.0021, +0.0003] | -0.0008 [-0.0013, -0.0002] |
| CORUM_large | -0.0019 [-0.0122, +0.0085] | -0.0091 [-0.0215, +0.0033] | +0.0024 [-0.0064, +0.0111] | +0.0017 [-0.0161, +0.0194] | -0.0506 [-0.0752, -0.0261] | -0.0594 [-0.0845, -0.0342] | -0.0471 [-0.0681, -0.0262] | -0.0513 [-0.0718, -0.0309] | -0.0600 [-0.0840, -0.0360] |
| CollecTRI | +0.0015 [-0.0033, +0.0064] | +0.0028 [-0.0023, +0.0080] | +0.0022 [-0.0024, +0.0068] | -0.0026 [-0.0076, +0.0023] | -0.0038 [-0.0086, +0.0009] | +0.0010 [-0.0038, +0.0058] | +0.0007 [-0.0041, +0.0055] | -0.0034 [-0.0082, +0.0014] | -0.0009 [-0.0056, +0.0038] |

#### Difference in enrichment vs K10 attention (t = 9) [95% interval], hub-stratified (pi of V0 attention)

| truth | G-SOFT | G-HARD | HALF-SOFT | R-SOFT | PERM-1 | PERM-2 | PERM-3 | PERM-4 | PERM-5 |
|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB | -0.01 [-0.04, +0.03] | -0.01 [-0.05, +0.02] | +0.01 [-0.02, +0.04] | -0.07 [-0.11, -0.04] | -0.04 [-0.07, -0.01] | -0.03 [-0.06, -0.01] | -0.06 [-0.09, -0.03] | -0.06 [-0.09, -0.03] | -0.05 [-0.08, -0.02] |
| CORUM_all | +0.18 [+0.03, +0.32] | +0.23 [+0.04, +0.42] | +0.26 [+0.12, +0.39] | +0.49 [+0.22, +0.75] | -0.52 [-0.78, -0.25] | -0.52 [-0.84, -0.21] | -0.39 [-0.65, -0.13] | -0.41 [-0.66, -0.16] | -0.57 [-0.88, -0.26] |
| KD_response | +0.04 [-0.00, +0.08] | +0.09 [+0.03, +0.15] | +0.08 [+0.04, +0.11] | -0.15 [-0.23, -0.07] | -0.24 [-0.31, -0.17] | -0.24 [-0.30, -0.17] | -0.20 [-0.27, -0.13] | -0.24 [-0.31, -0.17] | -0.23 [-0.29, -0.16] |
| CORUM_noRP | +0.14 [+0.03, +0.26] | +0.17 [+0.01, +0.33] | +0.23 [+0.12, +0.34] | +0.18 [+0.01, +0.34] | -0.25 [-0.40, -0.10] | -0.25 [-0.44, -0.07] | -0.20 [-0.36, -0.04] | -0.16 [-0.29, -0.03] | -0.29 [-0.44, -0.15] |
| CORUM_RP | +0.01 [-0.01, +0.02] | +0.01 [-0.01, +0.03] | +0.01 [-0.00, +0.02] | +0.04 [-0.01, +0.08] | -0.09 [-0.15, -0.03] | -0.08 [-0.13, -0.03] | -0.08 [-0.12, -0.03] | -0.09 [-0.13, -0.05] | -0.07 [-0.12, -0.02] |
| CORUM_small | +0.08 [-0.10, +0.26] | +0.06 [-0.14, +0.26] | +0.12 [-0.03, +0.27] | +0.08 [-0.18, +0.35] | -0.16 [-0.41, +0.10] | -0.21 [-0.45, +0.03] | -0.20 [-0.40, +0.00] | -0.27 [-0.63, +0.08] | -0.23 [-0.41, -0.05] |
| CORUM_large | -0.01 [-0.10, +0.07] | -0.07 [-0.17, +0.03] | +0.02 [-0.05, +0.09] | +0.01 [-0.13, +0.15] | -0.40 [-0.60, -0.21] | -0.47 [-0.67, -0.27] | -0.37 [-0.54, -0.21] | -0.41 [-0.57, -0.24] | -0.47 [-0.66, -0.28] |
| CollecTRI | +0.01 [-0.03, +0.05] | +0.02 [-0.02, +0.06] | +0.02 [-0.02, +0.05] | -0.02 [-0.06, +0.02] | -0.03 [-0.07, +0.01] | +0.01 [-0.03, +0.05] | +0.01 [-0.03, +0.04] | -0.03 [-0.06, +0.01] | -0.01 [-0.04, +0.03] |

#### Difference in AUPR vs CLR [95% paired delete-half gene jackknife], hub-stratified (pi of V0 attention)

| truth | G-SOFT | G-HARD | HALF-SOFT | R-SOFT | PERM-1 | PERM-2 | PERM-3 | PERM-4 | PERM-5 |
|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB | +0.0044 [-0.0003, +0.0090] | +0.0030 [-0.0022, +0.0081] | +0.0072 [+0.0030, +0.0113] | -0.0078 [-0.0136, -0.0019] | -0.0014 [-0.0075, +0.0046] | -0.0005 [-0.0058, +0.0048] | -0.0050 [-0.0111, +0.0012] | -0.0049 [-0.0105, +0.0007] | -0.0040 [-0.0097, +0.0017] |
| CORUM_all | -0.0152 [-0.0263, -0.0042] | -0.0145 [-0.0252, -0.0038] | -0.0141 [-0.0249, -0.0033] | -0.0111 [-0.0206, -0.0016] | -0.0244 [-0.0375, -0.0114] | -0.0245 [-0.0381, -0.0108] | -0.0227 [-0.0356, -0.0098] | -0.0231 [-0.0361, -0.0100] | -0.0251 [-0.0384, -0.0118] |
| KD_response | -0.0215 [-0.0277, -0.0153] | -0.0174 [-0.0237, -0.0112] | -0.0185 [-0.0245, -0.0125] | -0.0362 [-0.0441, -0.0284] | -0.0434 [-0.0519, -0.0350] | -0.0430 [-0.0515, -0.0345] | -0.0404 [-0.0481, -0.0326] | -0.0436 [-0.0513, -0.0360] | -0.0423 [-0.0498, -0.0347] |
| CORUM_noRP | +0.0002 [-0.0016, +0.0019] | +0.0005 [-0.0015, +0.0024] | +0.0012 [-0.0005, +0.0029] | +0.0006 [-0.0018, +0.0029] | -0.0045 [-0.0066, -0.0024] | -0.0045 [-0.0072, -0.0018] | -0.0039 [-0.0063, -0.0016] | -0.0034 [-0.0058, -0.0010] | -0.0050 [-0.0072, -0.0028] |
| CORUM_RP | -0.0136 [-0.0362, +0.0089] | -0.0109 [-0.0333, +0.0115] | -0.0119 [-0.0314, +0.0077] | +0.0069 [-0.0158, +0.0296] | -0.0870 [-0.1374, -0.0367] | -0.0812 [-0.1291, -0.0333] | -0.0760 [-0.1169, -0.0350] | -0.0862 [-0.1300, -0.0424] | -0.0718 [-0.1178, -0.0259] |
| CORUM_small | -0.0003 [-0.0019, +0.0013] | -0.0004 [-0.0020, +0.0013] | -0.0002 [-0.0018, +0.0014] | -0.0003 [-0.0021, +0.0015] | -0.0011 [-0.0029, +0.0008] | -0.0013 [-0.0033, +0.0008] | -0.0012 [-0.0031, +0.0006] | -0.0015 [-0.0034, +0.0005] | -0.0013 [-0.0029, +0.0003] |
| CORUM_large | -0.0332 [-0.0554, -0.0111] | -0.0404 [-0.0618, -0.0190] | -0.0289 [-0.0519, -0.0060] | -0.0297 [-0.0528, -0.0066] | -0.0820 [-0.1193, -0.0447] | -0.0907 [-0.1280, -0.0535] | -0.0785 [-0.1113, -0.0456] | -0.0827 [-0.1165, -0.0489] | -0.0913 [-0.1260, -0.0566] |
| CollecTRI | +0.0030 [-0.0009, +0.0068] | +0.0043 [+0.0001, +0.0085] | +0.0037 [-0.0002, +0.0076] | -0.0012 [-0.0056, +0.0033] | -0.0024 [-0.0072, +0.0024] | +0.0024 [-0.0026, +0.0075] | +0.0022 [-0.0033, +0.0076] | -0.0019 [-0.0071, +0.0032] | +0.0005 [-0.0049, +0.0060] |

#### Difference in enrichment vs CLR [95% interval], hub-stratified (pi of V0 attention)

| truth | G-SOFT | G-HARD | HALF-SOFT | R-SOFT | PERM-1 | PERM-2 | PERM-3 | PERM-4 | PERM-5 |
|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB | +0.02 [-0.00, +0.05] | +0.02 [-0.01, +0.04] | +0.04 [+0.02, +0.06] | -0.04 [-0.07, -0.01] | -0.01 [-0.04, +0.03] | -0.00 [-0.03, +0.03] | -0.03 [-0.06, +0.01] | -0.03 [-0.06, +0.00] | -0.02 [-0.05, +0.01] |
| CORUM_all | -1.15 [-1.98, -0.32] | -1.09 [-1.90, -0.29] | -1.07 [-1.88, -0.25] | -0.84 [-1.55, -0.12] | -1.84 [-2.83, -0.86] | -1.85 [-2.88, -0.82] | -1.71 [-2.69, -0.74] | -1.74 [-2.72, -0.76] | -1.89 [-2.89, -0.89] |
| KD_response | -0.27 [-0.35, -0.19] | -0.22 [-0.30, -0.14] | -0.24 [-0.31, -0.16] | -0.46 [-0.56, -0.36] | -0.55 [-0.66, -0.44] | -0.55 [-0.65, -0.44] | -0.51 [-0.61, -0.41] | -0.55 [-0.65, -0.46] | -0.54 [-0.63, -0.44] |
| CORUM_noRP | +0.01 [-0.14, +0.16] | +0.04 [-0.13, +0.20] | +0.10 [-0.05, +0.24] | +0.05 [-0.15, +0.25] | -0.38 [-0.56, -0.21] | -0.38 [-0.61, -0.16] | -0.33 [-0.53, -0.13] | -0.29 [-0.49, -0.09] | -0.42 [-0.61, -0.24] |
| CORUM_RP | -0.02 [-0.05, +0.01] | -0.01 [-0.04, +0.02] | -0.02 [-0.04, +0.01] | +0.01 [-0.02, +0.04] | -0.12 [-0.18, -0.05] | -0.11 [-0.17, -0.04] | -0.10 [-0.16, -0.05] | -0.11 [-0.17, -0.06] | -0.10 [-0.16, -0.03] |
| CORUM_small | -0.09 [-0.59, +0.41] | -0.11 [-0.60, +0.39] | -0.05 [-0.53, +0.43] | -0.09 [-0.63, +0.46] | -0.33 [-0.89, +0.23] | -0.38 [-0.99, +0.23] | -0.37 [-0.93, +0.18] | -0.44 [-1.03, +0.14] | -0.40 [-0.89, +0.09] |
| CORUM_large | -0.26 [-0.44, -0.09] | -0.32 [-0.49, -0.15] | -0.23 [-0.41, -0.05] | -0.24 [-0.42, -0.05] | -0.65 [-0.94, -0.35] | -0.72 [-1.01, -0.42] | -0.62 [-0.88, -0.36] | -0.65 [-0.92, -0.39] | -0.72 [-1.00, -0.45] |
| CollecTRI | +0.02 [-0.01, +0.05] | +0.03 [+0.00, +0.07] | +0.03 [-0.00, +0.06] | -0.01 [-0.04, +0.03] | -0.02 [-0.06, +0.02] | +0.02 [-0.02, +0.06] | +0.02 [-0.03, +0.06] | -0.02 [-0.06, +0.02] | +0.00 [-0.04, +0.05] |

#### Difference in AUPR vs K10 attention (t = 9) [95% paired delete-half gene jackknife], hub-stratified (MI row mean)

| truth | G-SOFT | G-HARD | HALF-SOFT | R-SOFT | PERM-1 | PERM-2 | PERM-3 | PERM-4 | PERM-5 |
|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB | +0.0004 [-0.0053, +0.0060] | +0.0005 [-0.0058, +0.0067] | +0.0029 [-0.0018, +0.0076] | -0.0100 [-0.0161, -0.0039] | -0.0081 [-0.0141, -0.0020] | -0.0074 [-0.0127, -0.0021] | -0.0127 [-0.0186, -0.0069] | -0.0132 [-0.0195, -0.0070] | -0.0094 [-0.0156, -0.0033] |
| CORUM_all | +0.0022 [+0.0003, +0.0042] | +0.0038 [+0.0009, +0.0068] | +0.0038 [+0.0021, +0.0055] | +0.0052 [+0.0011, +0.0093] | -0.0140 [-0.0206, -0.0075] | -0.0144 [-0.0225, -0.0064] | -0.0120 [-0.0185, -0.0055] | -0.0115 [-0.0178, -0.0053] | -0.0160 [-0.0236, -0.0083] |
| KD_response | +0.0120 [+0.0079, +0.0161] | +0.0180 [+0.0122, +0.0238] | +0.0174 [+0.0141, +0.0207] | -0.0185 [-0.0269, -0.0100] | -0.0419 [-0.0499, -0.0339] | -0.0386 [-0.0452, -0.0321] | -0.0413 [-0.0493, -0.0332] | -0.0465 [-0.0546, -0.0383] | -0.0466 [-0.0540, -0.0391] |
| CORUM_noRP | +0.0015 [-0.0001, +0.0030] | +0.0022 [-0.0000, +0.0045] | +0.0024 [+0.0011, +0.0038] | +0.0006 [-0.0020, +0.0033] | -0.0072 [-0.0111, -0.0033] | -0.0075 [-0.0121, -0.0029] | -0.0068 [-0.0110, -0.0025] | -0.0048 [-0.0081, -0.0015] | -0.0083 [-0.0125, -0.0041] |
| CORUM_RP | +0.0070 [-0.0055, +0.0196] | +0.0083 [-0.0051, +0.0218] | +0.0077 [-0.0038, +0.0192] | +0.0287 [-0.0018, +0.0592] | -0.0662 [-0.1076, -0.0247] | -0.0547 [-0.0925, -0.0168] | -0.0578 [-0.0898, -0.0259] | -0.0698 [-0.1010, -0.0387] | -0.0660 [-0.1078, -0.0242] |
| CORUM_small | -0.0000 [-0.0007, +0.0006] | +0.0002 [-0.0004, +0.0009] | -0.0000 [-0.0007, +0.0007] | -0.0001 [-0.0011, +0.0009] | -0.0011 [-0.0029, +0.0007] | -0.0015 [-0.0032, +0.0002] | -0.0011 [-0.0026, +0.0004] | -0.0018 [-0.0049, +0.0013] | -0.0015 [-0.0029, -0.0001] |
| CORUM_large | -0.0014 [-0.0093, +0.0065] | -0.0066 [-0.0164, +0.0032] | +0.0025 [-0.0043, +0.0092] | +0.0039 [-0.0114, +0.0193] | -0.0418 [-0.0644, -0.0192] | -0.0540 [-0.0777, -0.0303] | -0.0435 [-0.0626, -0.0244] | -0.0441 [-0.0635, -0.0247] | -0.0583 [-0.0823, -0.0344] |
| CollecTRI | +0.0044 [-0.0005, +0.0093] | +0.0059 [+0.0006, +0.0113] | +0.0062 [+0.0017, +0.0108] | -0.0035 [-0.0086, +0.0017] | -0.0097 [-0.0152, -0.0042] | -0.0038 [-0.0095, +0.0018] | -0.0055 [-0.0108, -0.0001] | -0.0119 [-0.0184, -0.0055] | -0.0078 [-0.0135, -0.0021] |

#### Difference in enrichment vs K10 attention (t = 9) [95% interval], hub-stratified (MI row mean)

| truth | G-SOFT | G-HARD | HALF-SOFT | R-SOFT | PERM-1 | PERM-2 | PERM-3 | PERM-4 | PERM-5 |
|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB | +0.00 [-0.03, +0.03] | +0.00 [-0.03, +0.04] | +0.02 [-0.01, +0.04] | -0.05 [-0.09, -0.02] | -0.04 [-0.08, -0.01] | -0.04 [-0.07, -0.01] | -0.07 [-0.10, -0.04] | -0.07 [-0.11, -0.04] | -0.05 [-0.09, -0.02] |
| CORUM_all | +0.17 [+0.02, +0.31] | +0.29 [+0.07, +0.51] | +0.29 [+0.16, +0.41] | +0.39 [+0.08, +0.70] | -1.06 [-1.56, -0.56] | -1.09 [-1.69, -0.48] | -0.91 [-1.40, -0.42] | -0.87 [-1.34, -0.40] | -1.20 [-1.78, -0.63] |
| KD_response | +0.15 [+0.10, +0.20] | +0.23 [+0.16, +0.30] | +0.22 [+0.18, +0.26] | -0.23 [-0.34, -0.13] | -0.53 [-0.63, -0.43] | -0.49 [-0.57, -0.41] | -0.52 [-0.63, -0.42] | -0.59 [-0.69, -0.49] | -0.59 [-0.69, -0.50] |
| CORUM_noRP | +0.12 [-0.01, +0.26] | +0.19 [-0.00, +0.38] | +0.20 [+0.09, +0.32] | +0.05 [-0.17, +0.27] | -0.60 [-0.93, -0.27] | -0.63 [-1.02, -0.25] | -0.57 [-0.93, -0.21] | -0.40 [-0.68, -0.12] | -0.70 [-1.05, -0.35] |
| CORUM_RP | +0.01 [-0.01, +0.03] | +0.01 [-0.01, +0.03] | +0.01 [-0.01, +0.03] | +0.04 [-0.00, +0.08] | -0.09 [-0.14, -0.03] | -0.07 [-0.12, -0.02] | -0.08 [-0.12, -0.03] | -0.09 [-0.13, -0.05] | -0.09 [-0.14, -0.03] |
| CORUM_small | -0.01 [-0.20, +0.18] | +0.07 [-0.13, +0.27] | -0.01 [-0.22, +0.20] | -0.04 [-0.35, +0.27] | -0.34 [-0.89, +0.20] | -0.46 [-0.98, +0.06] | -0.33 [-0.78, +0.12] | -0.54 [-1.47, +0.39] | -0.45 [-0.87, -0.02] |
| CORUM_large | -0.01 [-0.07, +0.05] | -0.05 [-0.13, +0.03] | +0.02 [-0.03, +0.07] | +0.03 [-0.09, +0.15] | -0.33 [-0.51, -0.15] | -0.43 [-0.62, -0.24] | -0.34 [-0.50, -0.19] | -0.35 [-0.50, -0.20] | -0.46 [-0.65, -0.27] |
| CollecTRI | +0.03 [-0.00, +0.07] | +0.05 [+0.00, +0.09] | +0.05 [+0.01, +0.08] | -0.03 [-0.07, +0.01] | -0.08 [-0.12, -0.03] | -0.03 [-0.07, +0.01] | -0.04 [-0.08, -0.00] | -0.09 [-0.14, -0.04] | -0.06 [-0.11, -0.02] |

#### Difference in AUPR vs CLR [95% paired delete-half gene jackknife], hub-stratified (MI row mean)

| truth | G-SOFT | G-HARD | HALF-SOFT | R-SOFT | PERM-1 | PERM-2 | PERM-3 | PERM-4 | PERM-5 |
|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB | +0.0033 [-0.0023, +0.0089] | +0.0034 [-0.0027, +0.0095] | +0.0058 [+0.0004, +0.0112] | -0.0071 [-0.0129, -0.0012] | -0.0051 [-0.0115, +0.0012] | -0.0045 [-0.0101, +0.0011] | -0.0098 [-0.0159, -0.0036] | -0.0103 [-0.0155, -0.0051] | -0.0065 [-0.0124, -0.0007] |
| CORUM_all | -0.0252 [-0.0412, -0.0091] | -0.0236 [-0.0393, -0.0078] | -0.0236 [-0.0394, -0.0078] | -0.0222 [-0.0366, -0.0078] | -0.0414 [-0.0619, -0.0209] | -0.0418 [-0.0632, -0.0204] | -0.0394 [-0.0596, -0.0193] | -0.0389 [-0.0591, -0.0187] | -0.0433 [-0.0642, -0.0224] |
| KD_response | +0.0115 [+0.0041, +0.0190] | +0.0175 [+0.0092, +0.0258] | +0.0169 [+0.0098, +0.0240] | -0.0190 [-0.0276, -0.0104] | -0.0424 [-0.0511, -0.0337] | -0.0391 [-0.0475, -0.0308] | -0.0418 [-0.0499, -0.0336] | -0.0470 [-0.0556, -0.0383] | -0.0471 [-0.0552, -0.0390] |
| CORUM_noRP | +0.0001 [-0.0028, +0.0029] | +0.0008 [-0.0023, +0.0039] | +0.0010 [-0.0018, +0.0038] | -0.0008 [-0.0040, +0.0024] | -0.0086 [-0.0131, -0.0040] | -0.0089 [-0.0142, -0.0036] | -0.0082 [-0.0127, -0.0036] | -0.0062 [-0.0103, -0.0020] | -0.0097 [-0.0141, -0.0054] |
| CORUM_RP | -0.0160 [-0.0404, +0.0083] | -0.0147 [-0.0395, +0.0100] | -0.0154 [-0.0366, +0.0059] | +0.0057 [-0.0216, +0.0329] | -0.0892 [-0.1392, -0.0393] | -0.0777 [-0.1270, -0.0284] | -0.0809 [-0.1239, -0.0379] | -0.0929 [-0.1365, -0.0492] | -0.0891 [-0.1377, -0.0405] |
| CORUM_small | -0.0008 [-0.0036, +0.0020] | -0.0006 [-0.0034, +0.0023] | -0.0008 [-0.0036, +0.0020] | -0.0009 [-0.0037, +0.0018] | -0.0019 [-0.0047, +0.0009] | -0.0023 [-0.0058, +0.0011] | -0.0019 [-0.0047, +0.0010] | -0.0026 [-0.0058, +0.0007] | -0.0023 [-0.0049, +0.0004] |
| CORUM_large | -0.0584 [-0.0828, -0.0341] | -0.0636 [-0.0873, -0.0398] | -0.0546 [-0.0791, -0.0300] | -0.0531 [-0.0755, -0.0306] | -0.0988 [-0.1363, -0.0613] | -0.1110 [-0.1468, -0.0752] | -0.1005 [-0.1322, -0.0688] | -0.1011 [-0.1336, -0.0686] | -0.1153 [-0.1503, -0.0804] |
| CollecTRI | +0.0118 [+0.0061, +0.0176] | +0.0134 [+0.0075, +0.0193] | +0.0137 [+0.0075, +0.0199] | +0.0040 [-0.0017, +0.0097] | -0.0023 [-0.0075, +0.0030] | +0.0036 [-0.0019, +0.0091] | +0.0020 [-0.0036, +0.0076] | -0.0045 [-0.0101, +0.0011] | -0.0004 [-0.0061, +0.0054] |

#### Difference in enrichment vs CLR [95% interval], hub-stratified (MI row mean)

| truth | G-SOFT | G-HARD | HALF-SOFT | R-SOFT | PERM-1 | PERM-2 | PERM-3 | PERM-4 | PERM-5 |
|---|---|---|---|---|---|---|---|---|---|
| DoRothEA_AB | +0.02 [-0.01, +0.05] | +0.02 [-0.01, +0.05] | +0.03 [+0.00, +0.06] | -0.04 [-0.07, -0.01] | -0.03 [-0.06, +0.01] | -0.02 [-0.06, +0.01] | -0.05 [-0.09, -0.02] | -0.06 [-0.08, -0.03] | -0.04 [-0.07, -0.00] |
| CORUM_all | -1.90 [-3.11, -0.69] | -1.78 [-2.96, -0.59] | -1.78 [-2.97, -0.59] | -1.67 [-2.76, -0.59] | -3.13 [-4.67, -1.58] | -3.15 [-4.76, -1.54] | -2.97 [-4.50, -1.45] | -2.93 [-4.46, -1.41] | -3.27 [-4.85, -1.69] |
| KD_response | +0.15 [+0.05, +0.24] | +0.22 [+0.12, +0.33] | +0.21 [+0.12, +0.30] | -0.24 [-0.35, -0.13] | -0.54 [-0.65, -0.43] | -0.50 [-0.60, -0.39] | -0.53 [-0.63, -0.43] | -0.60 [-0.71, -0.49] | -0.60 [-0.70, -0.50] |
| CORUM_noRP | +0.00 [-0.24, +0.25] | +0.07 [-0.19, +0.33] | +0.08 [-0.15, +0.32] | -0.07 [-0.33, +0.20] | -0.72 [-1.11, -0.34] | -0.75 [-1.20, -0.31] | -0.69 [-1.07, -0.30] | -0.52 [-0.87, -0.17] | -0.82 [-1.19, -0.45] |
| CORUM_RP | -0.02 [-0.05, +0.01] | -0.02 [-0.05, +0.01] | -0.02 [-0.05, +0.01] | +0.01 [-0.03, +0.04] | -0.12 [-0.19, -0.05] | -0.10 [-0.17, -0.04] | -0.11 [-0.17, -0.05] | -0.12 [-0.18, -0.07] | -0.12 [-0.18, -0.05] |
| CORUM_small | -0.25 [-1.10, +0.60] | -0.17 [-1.05, +0.71] | -0.25 [-1.10, +0.60] | -0.28 [-1.12, +0.56] | -0.58 [-1.44, +0.27] | -0.70 [-1.75, +0.35] | -0.57 [-1.44, +0.30] | -0.78 [-1.77, +0.21] | -0.69 [-1.49, +0.12] |
| CORUM_large | -0.46 [-0.66, -0.27] | -0.50 [-0.69, -0.32] | -0.43 [-0.63, -0.24] | -0.42 [-0.60, -0.24] | -0.78 [-1.08, -0.49] | -0.88 [-1.16, -0.60] | -0.80 [-1.05, -0.55] | -0.80 [-1.06, -0.54] | -0.91 [-1.19, -0.64] |
| CollecTRI | +0.09 [+0.05, +0.14] | +0.10 [+0.06, +0.15] | +0.11 [+0.06, +0.16] | +0.03 [-0.01, +0.08] | -0.02 [-0.06, +0.02] | +0.03 [-0.01, +0.07] | +0.02 [-0.03, +0.06] | -0.04 [-0.08, +0.01] | -0.00 [-0.05, +0.04] |



## 10. Full M3D tables (idea 1, lab splits)

#### Pooled enrichment (AUPR / base rate), raw

| truth | CLR | K10 | CONS-R | CONS-U | CONS-G | NULL-U |
|---|---|---|---|---|---|---|
| SC (base 0.0943) | 1.44 | 1.68 | 1.62 | 1.61 | 1.61 | 1.01 |
| all (base 0.1205) | 1.37 | 1.44 | 1.43 | 1.41 | 1.39 | 1.01 |

#### Pooled AUPR, raw

| truth | CLR | K10 | CONS-R | CONS-U | CONS-G | NULL-U |
|---|---|---|---|---|---|---|
| SC | 0.13587 | 0.15791 | 0.15282 | 0.15153 | 0.15154 | 0.09513 |
| all | 0.16479 | 0.17334 | 0.17190 | 0.17036 | 0.16800 | 0.12154 |

#### Pooled enrichment (AUPR / base rate), hub-stratified (pi of V0 attention)

| truth | CLR | K10 | CONS-R | CONS-U | CONS-G | NULL-U |
|---|---|---|---|---|---|---|
| SC (base 0.0943) | 1.40 | 1.56 | 1.64 | 1.59 | 1.63 | 1.01 |
| all (base 0.1205) | 1.33 | 1.35 | 1.45 | 1.41 | 1.43 | 1.01 |

#### Pooled AUPR, hub-stratified (pi of V0 attention)

| truth | CLR | K10 | CONS-R | CONS-U | CONS-G | NULL-U |
|---|---|---|---|---|---|---|
| SC | 0.13238 | 0.14715 | 0.15494 | 0.15018 | 0.15330 | 0.09513 |
| all | 0.16079 | 0.16275 | 0.17419 | 0.16982 | 0.17180 | 0.12130 |

#### Pooled enrichment (AUPR / base rate), hub-stratified (MI row mean)

| truth | CLR | K10 | CONS-R | CONS-U | CONS-G | NULL-U |
|---|---|---|---|---|---|---|
| SC (base 0.0943) | 1.50 | 1.68 | 1.67 | 1.63 | 1.61 | 1.01 |
| all (base 0.1205) | 1.41 | 1.44 | 1.45 | 1.42 | 1.40 | 1.01 |

#### Pooled AUPR, hub-stratified (MI row mean)

| truth | CLR | K10 | CONS-R | CONS-U | CONS-G | NULL-U |
|---|---|---|---|---|---|---|
| SC | 0.14140 | 0.15855 | 0.15717 | 0.15350 | 0.15143 | 0.09539 |
| all | 0.16989 | 0.17359 | 0.17460 | 0.17088 | 0.16812 | 0.12169 |

#### Difference in AUPR vs K10 attention (t = 9) [95% paired delete-half gene jackknife], raw

| truth | CONS-R | CONS-U | CONS-G | NULL-U |
|---|---|---|---|---|
| SC | -0.0051 [-0.0099, -0.0003] | -0.0064 [-0.0112, -0.0016] | -0.0064 [-0.0115, -0.0012] | -0.0628 [-0.0731, -0.0524] |
| all | -0.0014 [-0.0052, +0.0023] | -0.0030 [-0.0066, +0.0007] | -0.0053 [-0.0088, -0.0019] | -0.0518 [-0.0596, -0.0440] |

#### Difference in enrichment vs K10 attention (t = 9) [95% interval], raw

| truth | CONS-R | CONS-U | CONS-G | NULL-U |
|---|---|---|---|---|
| SC | -0.05 [-0.11, -0.00] | -0.07 [-0.12, -0.02] | -0.07 [-0.12, -0.01] | -0.67 [-0.78, -0.56] |
| all | -0.01 [-0.04, +0.02] | -0.02 [-0.06, +0.01] | -0.04 [-0.07, -0.02] | -0.43 [-0.49, -0.37] |

#### Difference in AUPR vs CLR [95% paired delete-half gene jackknife], raw

| truth | CONS-R | CONS-U | CONS-G | NULL-U |
|---|---|---|---|---|
| SC | +0.0170 [+0.0112, +0.0227] | +0.0157 [+0.0099, +0.0214] | +0.0157 [+0.0099, +0.0214] | -0.0407 [-0.0479, -0.0335] |
| all | +0.0071 [+0.0028, +0.0114] | +0.0056 [+0.0012, +0.0099] | +0.0032 [-0.0011, +0.0075] | -0.0433 [-0.0500, -0.0365] |

#### Difference in enrichment vs CLR [95% interval], raw

| truth | CONS-R | CONS-U | CONS-G | NULL-U |
|---|---|---|---|---|
| SC | +0.18 [+0.12, +0.24] | +0.17 [+0.11, +0.23] | +0.17 [+0.11, +0.23] | -0.43 [-0.51, -0.36] |
| all | +0.06 [+0.02, +0.09] | +0.05 [+0.01, +0.08] | +0.03 [-0.01, +0.06] | -0.36 [-0.41, -0.30] |

#### Difference in AUPR vs K10 attention (t = 9) [95% paired delete-half gene jackknife], hub-stratified (pi of V0 attention)

| truth | CONS-R | CONS-U | CONS-G | NULL-U |
|---|---|---|---|---|
| SC | +0.0078 [+0.0035, +0.0121] | +0.0030 [-0.0011, +0.0072] | +0.0061 [+0.0024, +0.0099] | -0.0520 [-0.0606, -0.0435] |
| all | +0.0114 [+0.0077, +0.0152] | +0.0071 [+0.0032, +0.0109] | +0.0091 [+0.0058, +0.0123] | -0.0415 [-0.0481, -0.0348] |

#### Difference in enrichment vs K10 attention (t = 9) [95% interval], hub-stratified (pi of V0 attention)

| truth | CONS-R | CONS-U | CONS-G | NULL-U |
|---|---|---|---|---|
| SC | +0.08 [+0.04, +0.13] | +0.03 [-0.01, +0.08] | +0.07 [+0.03, +0.11] | -0.55 [-0.64, -0.46] |
| all | +0.09 [+0.06, +0.13] | +0.06 [+0.03, +0.09] | +0.08 [+0.05, +0.10] | -0.34 [-0.40, -0.29] |

#### Difference in AUPR vs CLR [95% paired delete-half gene jackknife], hub-stratified (pi of V0 attention)

| truth | CONS-R | CONS-U | CONS-G | NULL-U |
|---|---|---|---|---|
| SC | +0.0226 [+0.0164, +0.0287] | +0.0178 [+0.0120, +0.0236] | +0.0209 [+0.0153, +0.0265] | -0.0373 [-0.0437, -0.0309] |
| all | +0.0134 [+0.0089, +0.0179] | +0.0090 [+0.0046, +0.0134] | +0.0110 [+0.0068, +0.0152] | -0.0395 [-0.0460, -0.0329] |

#### Difference in enrichment vs CLR [95% interval], hub-stratified (pi of V0 attention)

| truth | CONS-R | CONS-U | CONS-G | NULL-U |
|---|---|---|---|---|
| SC | +0.24 [+0.17, +0.30] | +0.19 [+0.13, +0.25] | +0.22 [+0.16, +0.28] | -0.40 [-0.46, -0.33] |
| all | +0.11 [+0.07, +0.15] | +0.07 [+0.04, +0.11] | +0.09 [+0.06, +0.13] | -0.33 [-0.38, -0.27] |

#### Difference in AUPR vs K10 attention (t = 9) [95% paired delete-half gene jackknife], hub-stratified (MI row mean)

| truth | CONS-R | CONS-U | CONS-G | NULL-U |
|---|---|---|---|---|
| SC | -0.0014 [-0.0063, +0.0035] | -0.0050 [-0.0098, -0.0003] | -0.0071 [-0.0124, -0.0018] | -0.0632 [-0.0734, -0.0530] |
| all | +0.0010 [-0.0027, +0.0047] | -0.0027 [-0.0065, +0.0011] | -0.0055 [-0.0091, -0.0019] | -0.0519 [-0.0598, -0.0440] |

#### Difference in enrichment vs K10 attention (t = 9) [95% interval], hub-stratified (MI row mean)

| truth | CONS-R | CONS-U | CONS-G | NULL-U |
|---|---|---|---|---|
| SC | -0.01 [-0.07, +0.04] | -0.05 [-0.10, -0.00] | -0.08 [-0.13, -0.02] | -0.67 [-0.78, -0.56] |
| all | +0.01 [-0.02, +0.04] | -0.02 [-0.05, +0.01] | -0.05 [-0.08, -0.02] | -0.43 [-0.50, -0.37] |

#### Difference in AUPR vs CLR [95% paired delete-half gene jackknife], hub-stratified (MI row mean)

| truth | CONS-R | CONS-U | CONS-G | NULL-U |
|---|---|---|---|---|
| SC | +0.0158 [+0.0102, +0.0214] | +0.0121 [+0.0067, +0.0175] | +0.0100 [+0.0049, +0.0152] | -0.0460 [-0.0532, -0.0389] |
| all | +0.0047 [+0.0009, +0.0085] | +0.0010 [-0.0029, +0.0049] | -0.0018 [-0.0056, +0.0020] | -0.0482 [-0.0554, -0.0410] |

#### Difference in enrichment vs CLR [95% interval], hub-stratified (MI row mean)

| truth | CONS-R | CONS-U | CONS-G | NULL-U |
|---|---|---|---|---|
| SC | +0.17 [+0.11, +0.23] | +0.13 [+0.07, +0.19] | +0.11 [+0.05, +0.16] | -0.49 [-0.56, -0.41] |
| all | +0.04 [+0.01, +0.07] | +0.01 [-0.02, +0.04] | -0.01 [-0.05, +0.02] | -0.40 [-0.46, -0.34] |


