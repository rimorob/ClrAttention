. jobs_pc/common.sh
# mem: 10
# cpus: 8
# after: 8000_wl_edge_m3d
# Width ladder (docs/plans/width_ladder_plan.md), score stage: single widths 2 to 128 under R9 (t = 9) and RE (t from the edge table), multi-scale scores MS1-3 for both rules, comparators CLR / tcaCLR / tcaCLR2, every metric of the plan on the spent truths (the reference gate passed earlier: results/width_ladder/<ds>/gate.csv is required). Exploratory, descriptive. No forked workers (BLAS threads only). Estimated about 15 min wall, main-process peak about 4 GB (smoke: 3.6 GB at 3 widths).
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/width_ladder/wl_run.R --mode score --ds m3d --threads 8
