. jobs_pc/common.sh
# mem: 20
# cpus: 4
# after: 7341_tc_grnboost2_armB
# Arm B: pack PIDC, GENIE3 and GRNBoost2 into the score directory (each symmetrised by max).
set -e
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript analysis/tcaclr_confirm/make_scores.R --mode pack --expr data/yeast/deleteome_all_1484.tsv --sha docs/tcaclr_confirm_selections/armB_expression.sha256 --out results/tcaclr_confirm/armB/scores --work results/tcaclr_confirm/armB
