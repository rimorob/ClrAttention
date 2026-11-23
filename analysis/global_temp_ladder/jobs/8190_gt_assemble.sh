. jobs_pc/common.sh
# mem: 4
# cpus: 2
# after: 8110_gt_score_m3d,8111_gt_score_ecomac,8112_gt_score_k562,8113_gt_score_yeast,8114_gt_score_kemmeren,8115_gt_score_hepg2
# Global-temperature ladder, assemble: response-surface tables, best temperature per truth with edge flags, effective-partner table, t table, ledger, gate summary, heatmaps and curves (PNG) in results/global_temp_ladder/. Datasets whose cells are missing (optional HepG2 / Kemmeren jobs not run) are skipped. About 1 min, under 2 GB.
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/global_temp_ladder/gt_assemble.R
