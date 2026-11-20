. jobs_pc/common.sh
# mem: 40
# cpus: 16
# tcaCLR depth/width/clip grid, stage 1 (label-free C1 diagnostic, C2 null-adjusted, C3 blocked held-out likelihood) on M3D (4,297 genes x 466): the first 8 stored laboratory splits, 16 operators + the reference operator, 17 depths. No truth is read. Exploratory (plan Revision 2). Estimated 2.5 h wall, 4 workers x 9 GB.
set -e
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=1
Rscript analysis/tcaclr_depth/stage1_labelfree.R --matrix results/baselines_m3d/expression.tsv --kind lab --splits 8 --workers 4 --threads 4 --out results/tcaclr_depth/stage1_m3d
