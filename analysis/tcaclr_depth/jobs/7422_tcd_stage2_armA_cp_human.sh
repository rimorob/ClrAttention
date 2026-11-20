. jobs_pc/common.sh
# mem: 46
# cpus: 8
# after: 7410_tcd_stage1_m3d,7411_tcd_stage1_k562
# tcaCLR depth/width/clip grid, stage 2: pooled co-membership AUPR (ap_fast) of every cell (aCLR and tcaCLR, both clip switches, 17 depths x 8 widths) on the spent truth armA_cp_human. Descriptive; the truth-chosen best cell is an upper bound, never a claim. Estimated about 1 h wall, 3 workers.
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
Rscript analysis/tcaclr_depth/stage2_truth_curves.R --truth armA_cp_human --workers 3 --threads 8 --out results/tcaclr_depth/stage2_armA_cp_human
