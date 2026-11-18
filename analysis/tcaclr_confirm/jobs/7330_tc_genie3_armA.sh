. jobs_pc/common.sh
# mem: 48
# cpus: 32
# after: 7320_tc_pidc_armA
# Arm A: GENIE3, all genes regulators, 1,000 trees, seed 20261016, 32 workers (estimate 5.0 to 5.8 h, declare 48 GB).
set -e
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
Rscript analysis/tcaclr_confirm/run_genie3.R results/tcaclr_confirm/armA 32
