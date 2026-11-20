. jobs_pc/common.sh
# mem: 14
# cpus: 8
# tcaCLR2 on the scoreboard rows, step 1b: permuted-data edge and depth t2 of the RPE1 matrix (8,747 genes, half A of the profiles, as clr_logpt_test.R), 10 permutations, no truth read. Exploratory. Output results/tcaclr2_scoreboard/edge_rpe1.tsv. Measured: 2 permutations with 2 workers x 4 threads took 148 s, worker peak 4.0 GB, main 0.65 GB; estimate for 10 permutations on 2 workers about 12 min and about 9 GB.
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/tcaclr2/edge_scoreboard.R --id rpe1 --workers 2 --threads 4
