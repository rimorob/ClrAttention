# Scoreboard figure (graphical abstract page 3): data matrix and compute plan

Written 13 Nov 2026 and committed BEFORE any compute. Purpose: a descriptive scoreboard of every metric we have, wins and losses, for Attention (CLR attention, t = 9, adopted), CLR, a3r, a3r2, a3r3 (exploratory), GENIE3, PIDC, GRNBoost2, |r|. Nothing here is a test: no p-values, no winner declared, ties are not marked (the method table marks Bonferroni-tied cells). K562 and M3D are tuning sets; a3r was picked among about 20 read-outs by its K562 scores; the a3r family is exploratory.

## 1. Data matrix (`docs/whitepapers/figures/scoreboard_data.csv`)

One row per dataset x truth x metric x scope x method. Columns: dataset, truth, metric, scope (raw, hub_att = the pi_i + pi_j stratifier of `docs/whitepapers/sota_table.md`, hub_mi where available), method, enrichment_over_random (pooled AUPR / base rate, or the stated metric), base_rate, source_file, status in {from_existing, computed_now, not_run}, note. Built by `analysis/scoreboard_matrix.R`, which reads (never recomputes) existing outputs and the four new result folders `results/scoreboard/{m3d,ecomac,k562,beeline}/`.

Existing sources (status from_existing):

| block | source |
|---|---|
| Attention, CLR, GENIE3, PIDC, GRNBoost2, \|r\| pooled AUPR, raw and hub_att (M3D SC/all, EcoMAC SC/all, K562 DoRothEA AB/ABC, CollecTRI, CORUM, KD response) | `results/baselines_{m3d,ecomac,k562ess_truthgenes_genie3}/metrics.csv` (the generator input of `analysis/sota_table.R`; `<method>` and `<method>_hubstrat`), cross-checked against the printed `sota_table.md` |
| K10 (= attention), CLR, a3r, a3r2, a3r3: K562 DoRothEA AB, CollecTRI, CORUM_all/noRP/RP/small/large/merged/merged_noRP, KD response; M3D regulons SC/all and complexes; scopes raw, hub_att, hub_mi | `results/explore_stack/{k562,m3d}/pooled_levels.csv` |
| M3D and EcoMAC complexes: attention, CLR, a3r, \|r\| | `results/clr_logpt_test/truth/{m3d,ecomac}/metrics.csv` and `docs/clr_logpt_verdict/endpoints_all.csv` |
| M3D and EcoMAC direct TF -> target edges: tf_aupr_ratio, tf_auroc_median (six old methods) | `results/baselines_{m3d,ecomac}/metrics.csv` columns `tf_aupr_ratio`, `tf_auroc_median` (RegulonDB TF-node benchmark, 199 TFs, SC and all universes) |
| per-regulon mean AP (six old methods), M3D SC/all | `results/baselines_m3d/per_regulon.csv` (regulons with >= 10 positive pairs, `ap` column of `analysis/baselines.R::regulon_ap_table`, ties in expectation); cross-checked with `results/macro_aupr_v3/macro_levels.csv` |
| BEELINE median AUPRC ratio (|r|, CLR, CLR-NS, GENIE3, PIDC, GRNBoost2) | `results/beeline_baselines/metrics.csv` and `summary.csv` (median over the complete primary case x truth cells: 14 cases x BEELINE's primary truths = 44 cells; the cell-type ChIP-only median over the 14 cases is carried as a second row) |

Per-regulon mean AP scale: enrichment = mean over regulons of AP divided by the mean of the regulons' base rates (1.00 = random), regulons with >= 10 positive pairs; `mean_ap_norm` ((AP - b)/(1 - b), the macro-table convention) and the plain mean AP are carried as extra metric rows in the csv. TF-edge AUROC is carried with chance 0.50.

Cross-check demanded by the lead: attention and CLR in `pooled_levels.csv` (K10 = attention) must equal `sota_table.md` to the printed precision; every mismatch is reported in the results message and logged in the csv note (none is hidden).

## 2. Cells computed now (status computed_now), after a gate

Same code paths as `analysis/explore_stack/stack_lib.R` (`sk_a3r`, `sk_layer`, layer-k operators), `analysis/explore_ideas47/ideas47_lib.R` (K562 set-up, truths, `rz_hub_stratify`), `analysis/baselines.R` (truths, TF-node universe, per-regulon table, `evaluate`), `analysis/beeline.R` / `beeline_baselines.R` (case set-up, CLR-NS, `score_truth`). No new scoring definitions. No jackknife (descriptive point estimates only; no intervals are shown on the figure).

| # | cells | job |
|---|---|---|
| 1 | M3D direct TF -> target edges: `tf_aupr_ratio`, `tf_auroc_median`, SC and all universes, for a3r, a3r2, a3r3 (also recomputed for attention and CLR as the gate) | `5600_sb_m3d` |
| 2 | M3D per-regulon mean AP, SC and all, for a3r, a3r2, a3r3 (and attention, CLR, \|r\| recomputed as the gate) | `5600_sb_m3d` |
| 3 | EcoMAC regulons SC and all, raw and hub_att, for a3r, a3r2, a3r3 (and attention, CLR as the gate); direct TF-edge metrics on EcoMAC as a by-product (same universes as `baselines.R`) | `5610_sb_ecomac` |
| 4 | K562 DoRothEA ABC raw and hub_att (and hub_mi) for a3r, a3r2, a3r3, attention, CLR | `5620_sb_k562` |
| 5 | BEELINE: for the a3r family on the CLR-NS operator (14 cases), median AUPRC ratio over the primary cells and over the cell-type ChIP cases | `5630_sb_beeline` |

Not computed (status not_run, shown as a dash): anything on the DREAM5 network-4 truth (held-out test, run once); any candidate truth of the upcoming confirmatory complex test (yeast Complex Portal, CYC2008, Kemmeren or any new complex truth: not opened, not scored); RPE1 (no script reuse planned); hub-controlled versions of direct edges, per-regulon AP and BEELINE (no existing stratified run); GENIE3/PIDC/GRNBoost2 on complexes; a3r-family and competitor cells that no existing file covers.

## 3. Equivalence gates (each job aborts before writing anything if one fails)

Tolerance 1e-9 absolute on AUPR (the same numerical path), 1e-12 on matrices.

- M3D (job 1): attention (K10 = S9) and CLR pooled AUPR on regulons SC/all, raw and hub_att, equal `results/baselines_m3d/metrics.csv` (aupr, attention_hubstrat, clr_hubstrat) and `pooled_levels.csv`; a3r, a3r2, a3r3 pooled AUPR (raw, hub_att) equal `results/explore_stack/m3d/pooled_levels.csv`; `tf_aupr_ratio`, `tf_auroc_median` of attention, CLR and |r| equal `results/baselines_m3d/metrics.csv`; per-regulon table of attention, CLR, |r| equal `per_regulon.csv` (AP, n_pos, base_rate).
- EcoMAC (job 2): attention, CLR pooled AUPR raw and hub_att on SC/all equal `results/baselines_ecomac/metrics.csv`; a3r on SC equals `docs/clr_logpt_verdict/endpoints_all.csv` (0.16140379470775...); TF-edge metrics of attention, CLR equal `metrics.csv`.
- K562 (job 3): CLR, K10, a3r, a3r2, a3r3 pooled AUPR on all ten truths of `pooled_levels.csv` x raw, hub_att, hub_mi equal the stored values; DoRothEA ABC attention and CLR raw equal `results/baselines_k562ess_truthgenes_genie3/metrics.csv`; DoRothEA ABC attention_hubstrat and clr_hubstrat compared with the pi stratifier of this job (any difference reported, not hidden: the baselines stratifier and `rz_hub_stratify` should coincide).
- BEELINE (job 4): per case, |r|, CLR and CLR-NS AUPRC ratio and EPR (primary truths) equal `results/beeline_baselines/metrics.csv`; attention on CLR-NS at the held-out t* of `results/beeline_v2/heldout_depths.csv` equals `att_clr_absr_blomavg` of `results/beeline_blom_test/metrics.csv` to 1e-9 on `ap_fast` (the scoring of that file); generic layer 1 (package operator vs `rz_operator`) 1e-12.

## 4. Design choices for BEELINE (not previously made by the project; most literal choice, not tuned)

- Operator: CLR-NS (Blom normal scores with average ranks, |r|, CLR Stouffer; `ns_clr` of `beeline.R`), the same scores the BEELINE attention of `beeline_blom_test.R` (`att_clr_absr_blomavg`) uses.
- Attention (the BEELINE method) keeps its protocol: held-out depth t* per case from `results/beeline_v2/heldout_depths.csv` (the depth was chosen label-free on cell halves). It is the Attention cell of the BEELINE row. Attention at t = 9 is also computed and carried in the csv.
- a3r family: t = 9 at every layer (as everywhere in the stacked-a3r plan), same softmax (keff 10, cap 50), lazy 0.5, readout `sk_a3r` unchanged (row statistics over reachable pairs), a3r2 and a3r3 stacked on the previous layer's z-scale matrix exactly as `explore_stack_plan.md`. No depth tuning. Therefore attention uses t* and the a3r family t = 9 in the BEELINE row; the figure footnote says so.
- Scoring: `score_truth` of `beeline_baselines.R` (pr_summary AUPRC ratio, tie-aware EPR), the cases and truths of `beeline.R` (7 datasets x TFs+500 and TFs+1000). Median over the complete primary case x truth cells (as `summary.csv`), and over the 14 cell-type ChIP cases.

## 5. Code paths and jobs (scheduler `jobs_pc/sched`; nice 19, capped scopes at 1.3 x declared; `gc()` before forking)

New code: `analysis/scoreboard/sb_m3d.R`, `sb_ecomac.R`, `sb_k562.R`, `sb_beeline.R`, `sb_lib.R` (shared gate helpers); matrix builder `analysis/scoreboard_matrix.R`; figure `analysis/fig_graphical_abstract_p3.R`. Outputs in `results/scoreboard/`.

| job | what | est. peak | declared (1.4 x) | cpus |
|---|---|---|---|---|
| `5600_sb_m3d` | MI/CLR, a3r layers 1-3, 4,297 genes; pooled AUPR raw/hub_att for SC, all; TF edges; per-regulon AP; gates | 6 GB | 9 GB | 8 |
| `5610_sb_ecomac` | ComBat, MI/CLR, layers 1-3; SC, all raw/hub_att; TF edges; gates | 7 GB | 10 GB | 8 |
| `5620_sb_k562` | sources `ideas47_lib.R` (10 truths, 5,671 genes), layers 1-3, pooled AUPR (no jackknife) 11 truths x 3 scopes | 9 GB | 13 GB | 8 |
| `5630_sb_beeline` | 14 cases, G <= 1,528 genes, small dense matrices | 3 GB | 5 GB | 4 |

Memory basis: G x G double matrix 148 MB (M3D), about 150 MB (EcoMAC), 257 MB (K562); about 12 live at once; the scoring uses upper-triangle vectors only. The earlier K562 jobs peaked at 20 GB only because of 8 forked jackknife workers, which are not used here. Total declared 37 GB (utilisation about 40% of the 90 GB budget, jobs run concurrently). An OOM is reported; the scheduler retries once at 1.5 x.

Disk: no matrices are kept; result CSVs and logs only (< 5 MB); row in `docs/disk_budget.md`.

## 6. Figure (`analysis/fig_graphical_abstract_p3.R`)

R grid, conventions and palette of `fig_graphical_abstract.R` / `_p2.R` (plot_theme.R, Liberation Sans, 2400 x 1350 px png plus 16 x 9 pdf; `docs/whitepapers/figures/graphical_abstract_p3.{png,pdf}`). Title "Every metric, wins and losses". Rows = metrics grouped by dataset, columns = nine methods (a3r family dashed, labelled "exploratory"; attention labelled "adopted"; no column is labelled "this work"). Each cell: the number and a bar from random (1.00 tick), length proportional to value - 1 scaled within the row; the highest value of the row is bold with a dot (ties not marked); missing values: grey dash (footnote "not run"). Raw block and hub-controlled block (split into two figures if the type would fall below 9 pt).
