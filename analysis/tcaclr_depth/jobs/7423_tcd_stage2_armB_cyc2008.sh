. jobs_pc/common.sh
# mem: 24
# cpus: 8
# after: 7410_tcd_stage1_m3d,7411_tcd_stage1_k562
# tcaCLR depth/width/clip grid, stage 2: pooled co-membership AUPR (ap_fast) of every cell (aCLR and tcaCLR, both clip switches, 17 depths x 8 widths) on the spent truth armB_cyc2008. Descriptive; the truth-chosen best cell is an upper bound, never a claim. Estimated about 15 min wall, 4 workers.
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
Rscript analysis/tcaclr_depth/stage2_truth_curves.R --truth armB_cyc2008 --workers 4 --threads 8 --out results/tcaclr_depth/stage2_armB_cyc2008
