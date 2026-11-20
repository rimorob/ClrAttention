. jobs_pc/common.sh
# mem: 12
# cpus: 16
# tcaCLR adaptive-setting study on E. coli M3D (4,297 genes x 466): spectral-edge depth rule (10 permutations, 8 widths), blocked held-out graph-signal likelihood over the 17 depths x 8 widths on the 8 stored laboratory splits, and the per-gene-permuted negative control (20 controls, widths 10/28/80). No truth is read. Exploratory (plan Revision 2 items 3, 8, 11). Estimated 1.5 h wall, 4 workers x 4 threads, measured peak 1.7 GB per worker.
set -e
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=1
Rscript analysis/tcaclr_depth/adaptive/run_adaptive.R --dataset m3d --splits 8 --workers 4 --threads 4 --out results/tcaclr_depth/adaptive_m3d
