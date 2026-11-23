. jobs_pc/common.sh
# mem: 4
# cpus: 2
# after: 8010_wl_score_m3d,8011_wl_score_ecomac,8012_wl_score_k562,8013_wl_score_yeast,8014_wl_score_kemmeren,8015_wl_score_hepg2
# Width ladder, assemble: response-surface tables, best width per truth with edge flags, t table, ledger, heatmaps and curves (PNG) in results/width_ladder/. Datasets whose cells are missing (optional HepG2 / Kemmeren jobs not run) are skipped. About 1 min, under 2 GB.
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/width_ladder/wl_assemble.R
