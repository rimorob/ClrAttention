. jobs_pc/common.sh
# mem: 62
# cpus: 16
# after: 7410_tcd_stage1_m3d
# tcaCLR depth/width/clip grid, stage 1 on K562 essential Perturb-seq (5,671 genes x 707 knockdowns): 8 splits by knockdown cluster (Ward on 1 - r, number of clusters = round N_eff), seed 20261120. No truth is read. Exploratory (plan Revision 2). Estimated 5 h wall, 4 workers x 15 GB (use --workers 3 and mem 48 for a lower footprint, about 7 h).
set -e
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=1
Rscript analysis/tcaclr_depth/stage1_labelfree.R --matrix results/baselines_k562ess_truthgenes_genie3/expression.tsv --kind random --splits 8 --workers 4 --threads 4 --out results/tcaclr_depth/stage1_k562
