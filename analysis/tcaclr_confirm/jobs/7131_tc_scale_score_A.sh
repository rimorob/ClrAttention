. jobs_pc/common.sh
# mem: 28
# cpus: 16
# after: 7111_tc_val_net4_refit_series,7121_tc_val_net4_refit_chip
# tcaCLR scoring code scale test, SYNTHETIC data only: Arm-A-sized score directory (9,000 genes, truth of 3,500 genes, 7 methods + 4 read-outs), gates + score + report stage timings and peak memory.
set -e
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript analysis/tcaclr_confirm/scale_test.R --mode score --G 9000 --n_truth 3500 --workers 16 --out results/tcaclr_confirm/scale_A
rm -rf results/tcaclr_confirm/scale_A/scores
