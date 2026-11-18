. jobs_pc/common.sh
# mem: 16
# cpus: 8
# Arm A: CLR, aCLR, tcaCLR, |r| and the hub variable written as packed score files (frozen plan section 3). No truth is read.
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
Rscript analysis/tcaclr_confirm/make_scores.R --mode fit --expr results/tcaclr_confirm/armA/expression.tsv --sha docs/tcaclr_confirm_selections/armA_expression.sha256 --out results/tcaclr_confirm/armA/scores --threads 8
