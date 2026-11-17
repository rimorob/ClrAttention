. jobs_pc/common.sh
# mem: 30
# cpus: 16
# after: 7131_tc_scale_score_A
# tcaCLR refit scale test, SYNTHETIC data only: one half-refit (CLR, aCLR, tcaCLR, |r|) at Arm A size (9,000 genes, 750 of 1,500 profiles) and Arm B size (6,169 genes, 642 of 1,283 profiles).
set -e
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
Rscript analysis/tcaclr_confirm/scale_test.R --mode refit --G 6169 --n_prof 1283 --workers 16 --out results/tcaclr_confirm/scale_B
Rscript analysis/tcaclr_confirm/scale_test.R --mode refit --G 9000 --n_prof 1500 --workers 16 --out results/tcaclr_confirm/scale_A
