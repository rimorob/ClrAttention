. jobs_pc/common.sh
# mem: 3
# cpus: 2
# tcaCLR2 on the scoreboard rows, step 1c: edge and t2 of the 14 BEELINE cases (CLR-NS operator, the one extension of the rule: the operator is built from CLR-NS because that is the score the scoreboard's BEELINE cells use; flagged in the report), 10 permutations per case. Exploratory. Output results/tcaclr2_scoreboard/edge_beeline_<dataset>_<N>.tsv. Measured: all 14 cases about 4 min, peak under 1 GB.
set -e
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/tcaclr2/edge_scoreboard.R --id beeline --workers 1 --threads 1
