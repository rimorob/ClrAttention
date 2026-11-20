. jobs_pc/common.sh
# mem: 12
# cpus: 8
# tcaCLR planted-module simulation: 3,000 random genes of the M3D matrix, each gene's values permuted across experiments (real null), modules of 5/10/20/50/100 genes planted at correlation 0.1-0.6 with a linear and a nonlinear link (60 datasets); CLR, aCLR, tcaCLR and the exploratory diffusion-PMI read-out over keff 10/28 x 7 depths; oracle best cell, regret of the edge rule and of the graph-signal likelihood rule. Synthetic data only. Exploratory (plan Revision 2 items 9, 11, 12). Estimated 1 h wall, 4 workers x 2 threads, measured peak 1.9 GB per worker.
set -e
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
Rscript analysis/tcaclr_depth/adaptive/planted_sim.R --matrix results/baselines_m3d/expression.tsv --dataset M3D --genes 3000 --workers 4 --threads 2 --out results/tcaclr_depth/adaptive_planted_m3d
