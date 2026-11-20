. jobs_pc/common.sh
# mem: 5
# cpus: 8
# tcaCLR2 on the scoreboard rows, step 2a: M3D cells (regulons SC / all and complexes raw + hub-controlled, direct TF edges, per-regulon AP) for tcaCLR2 (t2 from results/tcaclr2/edge_m3d_tf_regulons.tsv) and for tcaCLR at t = 9 as the gate against the stored scoreboard values (stops on a difference above 1e-9). Exploratory. Output results/tcaclr2_scoreboard/cells_m3d.csv, gate_m3d.csv. Measured: 130 s, peak RSS 2.95 GB, 16 cells gated, max difference 4e-14.
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/tcaclr2/score_scoreboard.R --family m3d --threads 8
