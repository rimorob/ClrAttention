. jobs_pc/common.sh
# mem: 8
# cpus: 8
# after: 7301_tc_prep_armB
# Arm B: CLR, aCLR, tcaCLR, |r| and the hub variable written as packed score files (frozen plan section 3). No truth is read.
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
Rscript analysis/tcaclr_confirm/make_scores.R --mode fit --expr data/yeast/deleteome_all_1484.tsv --sha docs/tcaclr_confirm_selections/armB_expression.sha256 --out results/tcaclr_confirm/armB/scores --threads 8
