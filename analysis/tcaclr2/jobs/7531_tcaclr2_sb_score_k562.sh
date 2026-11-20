. jobs_pc/common.sh
# mem: 8
# cpus: 8
# tcaCLR2 on the scoreboard rows, step 2b: K562 cells (11 truths, raw + hub-controlled) for tcaCLR2 (t2 from results/tcaclr2/edge_k562_corum.tsv) and the tcaCLR t = 9 gate. Exploratory. Output results/tcaclr2_scoreboard/cells_k562.csv, gate_k562.csv. Measured: 234 s, peak RSS 5.3 GB, 22 cells gated, max difference 5e-15.
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/tcaclr2/score_scoreboard.R --family k562 --threads 8
