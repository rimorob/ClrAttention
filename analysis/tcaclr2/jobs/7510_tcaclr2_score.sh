. jobs_pc/common.sh
# mem: 10
# cpus: 8
# after: 7500_tcaclr2_edge_depth
# tcaCLR2 step 3: pooled co-membership AUPR (ap_fast, the stage-2 truths and code of the depth grid) of the registered tcaCLR (t = 9, keff 10, cap 50) and tcaCLR2 (t2 from 7500) side by side on the five spent truths, the ratio tcaCLR2 / tcaCLR, and the gene delete-half jackknife for those two cells only; results/tcaclr2/summary.csv. Exploratory; reading rule in the definition note. Estimated about 15 min wall, one dataset at a time.
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
analysis/tcaclr2/run_tcaclr2.sh score k562_corum net4_cp_yeast armA_cp_human armB_cyc2008 m3d_tf_regulons
