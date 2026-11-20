. jobs_pc/common.sh
# mem: 10
# cpus: 8
# after: 7521_tcaclr2_sb_edge_rpe1
# tcaCLR2 on the scoreboard rows, step 2d: RPE1 CORUM-merged cell (raw + hub-controlled) for tcaCLR2 (t2 from 7521) and the tcaCLR t = 9 gate. Exploratory. Output results/tcaclr2_scoreboard/cells_rpe1.csv, gate_rpe1.csv. Measured: 437 s, peak RSS 7.1 GB, 2 cells gated, max difference 2e-14.
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/tcaclr2/score_scoreboard.R --family rpe1 --threads 8
