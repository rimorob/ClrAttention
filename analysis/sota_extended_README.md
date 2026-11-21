# Extended SOTA comparison (tcaCLR and tcaCLR2 columns)

`analysis/sota_extended.R` extends `analysis/sota_table.R` / `docs/whitepapers/sota_table.md` with tcaCLR (registered, t = 9) and tcaCLR2 (depth t2 of the permuted-data edge rule, `docs/plans/tcaclr2_definition.md`) on the same 18 rows (9 raw, 9 hub-controlled): M3D SC and all, EcoMAC SC and all, K562 DoRothEA AB, ABC, CollecTRI, CORUM, KD response.

**Status and limits.** tcaCLR2 is exploratory (outside the preregistration; the depth-24 and depth-32 grid curves had been seen on the spent truths). All intervals resample genes only (delete-half jackknife, 100 replicates, fitted networks held fixed); they do not resample experiments or cells or vary the seeds of GENIE3 / GRNBoost2. No interpretation is made in the outputs.

## What is computed
Per matrix, nine columns: aCLR (= "Attention" of the old table, symmetrised P^9), tcaCLR (second CLR on log P^9, `sk_a3r`), tcaCLR2 (same at t2: M3D 25, EcoMAC 18, K562 29), CLR, GENIE3, PIDC, GRNBoost2, |r|, Hub only (raw rows only; a degree baseline, never competes for a mark). Pooled co-membership AUPR / base rate, raw and hub-controlled (within-stratum percentile in 10 strata of pi_i + pi_j, pi of the attention walk, as before). The jackknife halves are those of `analysis/baselines.R` (seed 20261016, drawn per truth in its order, 100 draws; the first `--jack` are used), so the six old methods are paired exactly with their stored replicates. Nothing is re-implemented: operator `rz_operator`/`rz_lazy_P` (checked equal to the package operator), the registered sequential chain (`tsb_walk`; K562 continues the chain of `ideas47_lib.R` as `score_scoreboard.R`), `sk_a3r` / `a3r`, `rz_stationary`, `rz_hub_stratify`, truths of `sb_load_ecoli` and `ideas47_lib.R`, `ap_fast`. GENIE3, PIDC and GRNBoost2 are the stored `results/baselines_*/{genie3,pidc,grnboost2}_scores.rds`. The mark logic is the one of `sota_table.R` (copied into the assemble step because `sota_table.R` runs on `source`); the assemble gate shows it reproduces the old marks.

Marks (per row and scope): best method in the set, Bonferroni over (set size - 1) comparisons, SE = sd of the paired replicate differences. Three sets: `eight` (aCLR, tcaCLR, tcaCLR2, CLR, GENIE3, PIDC, GRNBoost2, |r|), `six` (without the two new), `legacy` (the six, plus Hub only in raw rows = the convention of the old table; used for the reproduction gate), and `eight_plus_hub` (raw only).

Paired differences (`paired_differences.csv`): tcaCLR - CLR, tcaCLR - aCLR, tcaCLR2 - CLR, tcaCLR2 - aCLR, tcaCLR2 - tcaCLR, for all 18 rows; point difference of pooled AUPR (and divided by the base rate), SE = sd of the replicate differences, interval = difference +/- 1.96 SE (normal, delete-d factor 1, as `baselines.R`), z, `interval_excludes_0`. Uncorrected intervals (no multiplicity adjustment across rows or comparisons).

## Run (repo root)
`Rscript analysis/sota_extended.R --matrix m3d|ecomac|k562` (one job each), then `--matrix assemble`. Options: `--out DIR` (default `results/sota_extended`), `--jack N` (default 100, 2 to 100), `--tag`, `--threads` (8), `--matrices` (assemble). Job files `analysis/tcaclr2/jobs/7550_sota_ext_{m3d,ecomac,k562,assemble}.sh` (the directory is gitignored: force-add; assemble has a comma-separated `# after:` line).

## Inputs
`results/baselines_{m3d,ecomac,k562ess_truthgenes_genie3}/` (`*_scores.rds`, `metrics.csv`, `jack_replicates.csv`, `expression.tsv` for K562), the raw expression and truth data read by `sb_load_ecoli` / `ideas47_lib.R`, `results/tcaclr2/edge_m3d_tf_regulons.tsv`, `results/tcaclr2/edge_k562_corum.tsv`, `results/tcaclr2_scoreboard/edge_ecomac.tsv` (t2), `results/tcaclr2_scoreboard/cells_reference_tcaclr_all.csv`, `cells_tcaclr2_all.csv` (gate), `docs/whitepapers/sota_table.md` (assemble gate).

## Outputs (`results/sota_extended/`; an existing file is never overwritten)
Per matrix: `rows_<m>.csv` (dataset, scope, truth, method, AUPR, base rate, enrichment, t2), `jack_<m>.csv` (replicate AUPRs), `gate_<m>.csv` (43 checks), `run_<m>.csv` (wall, peak RSS, per-method seconds), `sota_extended_<m>.log`; K562 also leaves `work_k562/run.log`. On a failed gate only `gate_<m>_FAILED.csv` is written. Assemble: `table_extended.md` / `.csv` (values, marks for every set), `marks_by_row.csv`, `paired_differences.csv`, `gate_assemble.csv` (or `gate_assemble_FAILED.csv`).

## Gate (a job stops before writing results on any failure)
1. tcaCLR (t = 9) and tcaCLR2 pooled AUPR/base, raw and hub_att, against `cells_reference_tcaclr_all.csv` / `cells_tcaclr2_all.csv`, 1e-9.
2. Old six methods and Hub only: AUPR and base rate against `metrics.csv`, every replicate against `jack_replicates.csv`, raw and hub-stratified, 1e-9.
3. Assemble: printed value (two decimals) and bold / dagger mark of every old-table cell against `docs/whitepapers/sota_table.md` (legacy convention); any difference stops the step. A deliberately altered copy of the table was rejected in a test.

## Declared resources, measured and smoke times (`nice 19`, 8 threads/cores, one job at a time)
| job | declared | measured |
|---|---|---|
| 7550_sota_ext_m3d | mem 6 GB, cpus 8 | 229 s at 100 replicates (178 s at 10), peak RSS 3.84 GB (3.20 GB at 10) |
| 7550_sota_ext_ecomac | mem 6 GB, cpus 8 | 234 s at 10 replicates, peak RSS 3.5 GB; about 290 s and 4.0 GB expected at 100 (M3D: +51 s, +0.6 GB from 10 to 100) |
| 7550_sota_ext_k562 | mem 13 GB, cpus 8 | 742 s at 100 replicates, peak RSS 9.3 GB (9.7 GB by `time -v`) |
| 7550_sota_ext_assemble | mem 2 GB, cpus 1 | seconds |

Wall time if run one after the other: about 21 min (about 13 min for K562, the longest). Smoke runs (M3D 10 and 100 replicates, EcoMAC 10, K562 100) are in `results/sota_extended/smoke/` (tags `_j100`, `_full`, `_j10`); every gate passed (43 checks each, max abs difference 4.9e-15 M3D, 4.4e-15 EcoMAC, 5.1e-15 K562; `gate_assemble` 26 cells for M3D alone, 91 cells for M3D + K562). The M3D and K562 smoke runs are full-size, 100-replicate runs and could be copied to the production names if desired; they were kept apart so as not to occupy the production file names.
