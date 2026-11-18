. jobs_pc/common.sh
# mem: 40
# cpus: 4
# after: 7340_tc_grnboost2_armA
# Arm A: pack PIDC, GENIE3 and GRNBoost2 into the score directory (each symmetrised by max).
set -e
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript analysis/tcaclr_confirm/make_scores.R --mode pack --expr results/tcaclr_confirm/armA/expression.tsv --sha docs/tcaclr_confirm_selections/armA_expression.sha256 --out results/tcaclr_confirm/armA/scores --work results/tcaclr_confirm/armA
