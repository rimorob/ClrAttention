. jobs_pc/common.sh
# mem: 20
# cpus: 16
# tcaCLR adaptive-setting study on K562 essential Perturb-seq (5,671 genes x 707 knockdowns): spectral-edge depth rule, blocked held-out graph-signal likelihood on 8 splits by k-means knockdown cluster (k = round(N_eff)), and the per-gene-permuted negative control. No truth is read. Exploratory (plan Revision 2 items 3, 8, 11). Estimated 2 h wall, 4 workers x 4 threads, measured peak 2.9 GB per worker.
set -e
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=1
Rscript analysis/tcaclr_depth/adaptive/run_adaptive.R --dataset k562 --splits 8 --workers 4 --threads 4 --out results/tcaclr_depth/adaptive_k562
