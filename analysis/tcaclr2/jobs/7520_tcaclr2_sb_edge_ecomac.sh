. jobs_pc/common.sh
# mem: 14
# cpus: 8
# tcaCLR2 on the scoreboard rows, step 1a (docs/plans/tcaclr2_definition.md, rule unchanged): permuted-data edge and depth t2 of the EcoMAC matrix (ComBat-adjusted, 2,198 arrays), 10 per-gene permutations, no truth read. Exploratory. Output results/tcaclr2_scoreboard/edge_ecomac.tsv. Measured: 2 permutations with 2 workers x 2 threads took 119 s, worker peak 2.3 GB, main process 2.5 GB (ComBat); estimate for 4 workers x 2 threads about 6 min and about 12 GB.
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/tcaclr2/edge_scoreboard.R --id ecomac --workers 4 --threads 2
