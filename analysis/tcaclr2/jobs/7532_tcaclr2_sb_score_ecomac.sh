. jobs_pc/common.sh
# mem: 5
# cpus: 8
# after: 7520_tcaclr2_sb_edge_ecomac
# tcaCLR2 on the scoreboard rows, step 2c: EcoMAC cells (regulons SC / all and complexes raw + hub-controlled, direct TF edges, per-regulon AP) for tcaCLR2 (t2 from 7520) and the tcaCLR t = 9 gate. Exploratory. Output results/tcaclr2_scoreboard/cells_ecomac.csv, gate_ecomac.csv. Measured: 158 s, peak RSS 3.15 GB, 16 cells gated, max difference 4e-14.
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/tcaclr2/score_scoreboard.R --family ecomac --threads 8
